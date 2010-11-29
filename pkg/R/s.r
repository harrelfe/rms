survfit.coxph <- function (formula, newdata, se.fit = TRUE, conf.int = 0.95, individual = FALSE, 
    type, vartype, conf.type = c("log", "log-log", "plain", "none"), 
    censor = TRUE, ...) 
{
    Call <- match.call()
    Call[[1]] <- as.name("survfit")
    object <- formula
    if (missing(type)) {
        temp1 <- c("exact", "breslow", "efron")
        survtype <- match(object$method, temp1)
    }
    else {
        temp1 <- c("kalbfleisch-prentice", "aalen", "efron", 
            "kaplan-meier", "breslow", "fleming-harrington", 
            "greenwood", "tsiatis", "exact")
        survtype <- match(match.arg(type, temp1), temp1)
        survtype <- c(1, 2, 3, 1, 2, 3, 1, 2, 3)[survtype]
    }
    if (missing(vartype)) {
        vartype <- survtype
    }
    else {
        temp2 <- c("greenwood", "aalen", "efron", "tsiatis")
        vartype <- match(match.arg(vartype, temp2), temp2)
        if (vartype == 4) 
            vartype <- 2
    }
    if (!se.fit) 
        conf.type <- "none"
    else conf.type <- match.arg(conf.type)
    if (is.null(object$y) || is.null(object[["x"]]) || !is.null(object$call$weights) || 
        !is.null(attr(object$terms, "specials")$strata) || !is.null(attr(object$terms, 
        "offset"))) {
        mf <- model.frame(object)
    }
    else mf <- NULL
    if (is.null(object[["y"]])) 
        y <- model.response(mf)
    else y <- object[["y"]]
    if (is.null(object[["x"]])) 
        x <- model.matrix.coxph(object, mf = mf)
    else x <- object[["x"]]
    n <- nrow(y)
    if (n != object$n[1] || nrow(x) != n) 
        stop("Failed to reconstruct the original data set")
    if (is.null(mf)) 
        wt <- rep(1, n)
    else {
        wt <- model.weights(mf)
        if (is.null(wt)) 
            wt <- rep(1, n)
    }
    type <- attr(y, "type")
    if (type != "right" && type != "counting") 
        stop("Cannot handle \"", type, "\" type survival data")
    if (is.null(mf)) 
        offset <- 0
    else {
        offset <- model.offset(mf)
        if (is.null(offset)) 
            offset <- 0
    }
    Terms <- object$terms
    temp <- untangle.specials(Terms, "strata")
    if (length(temp$terms) == 0) 
        strata <- rep(0L, n)
    else strata <- mf[[temp$vars]]
    if (is.null(x) || ncol(x) == 0) {
        x <- matrix(0, nrow = n)
        coef <- 0
        varmat <- matrix(0, 1, 1)
        risk <- rep(exp(offset - mean(offset)), length = n)
    }
    else {
        varmat <- object$var
        coef <- ifelse(is.na(object$coefficients), 0, object$coefficients)
        xcenter <- object$means
        if (is.null(object$frail)) {
            x <- scale(x, center = xcenter, scale = FALSE)
            risk <- c(exp(x %*% coef + offset - mean(offset)))
        }
        else {
            x <- x[, !is.na(match(dimnames(x)[[2]], names(coef))), 
                drop = F]
            risk <- exp(object$linear.predictor)
            x <- scale(x, center = xcenter, scale = FALSE)
        }
    }
    if (individual) {
        if (missing(newdata)) 
            stop("The newdata argument must be present")
        if (!is.null(object$frail)) 
            stop("The newdata argument is not supported for sparse frailty terms")
        if (!is.data.frame(newdata)) 
            stop("Newdata must be a data frame")
        temp <- untangle.specials(Terms, "cluster")
        if (length(temp$vars)) 
            Terms2 <- Terms[-temp$terms]
        else Terms2 <- Terms
        mf2 <- model.frame(Terms2, newdata, xlev = object$xlevels)
        temp <- untangle.specials(Terms2, "strata")
        if (length(temp$vars) > 0) {
            strata2 <- strata(mf2[temp$vars], shortlabel = TRUE)
            strata2 <- factor(strata2, levels = levels(strata))
            if (any(is.na(strata2))) 
                stop("New data set has strata levels not found in the original")
            Terms2 <- Terms2[-temp$terms]
        }
        else strata2 <- rep(0, nrow(mf2))
        x2 <- model.matrix(Terms2, mf2)[, -1, drop = FALSE]
        if (length(x2) == 0) 
            x2 <- matrix(0, nrow = nrow(mf2), ncol = 1)
        else x2 <- scale(x2, center = xcenter, scale = FALSE)
        offset2 <- model.offset(mf2)
        if (length(offset2) > 0) 
            offset2 <- offset2 - mean(offset)
        else offset2 <- 0
        y2 <- model.extract(mf2, "response")
        if (attr(y2, "type") != type) 
            stop("Survival type of newdata does not match the fitted model")
    }
    else {
        if (missing(newdata)) {
            x2 <- matrix(0, nrow = 1, ncol = ncol(x))
            offset2 <- 0
        }
        else {
            if (!is.null(object$frail)) 
                stop("The newdata argument is not supported for sparse frailty terms")
            if (!is.data.frame(newdata)) {
                if (is.list(newdata)) 
                  newdata <- data.frame(newdata)
                else if (is.numeric(newdata) && length(newdata) == 
                  length(object$coefficients)) {
                  if (is.null(names(newdata))) 
                    names(newdata) <- names(object$coefficients)
                  newdata <- data.frame(as.list(newdata))
                }
                else stop("Invalid newdata object")
            }
            Terms2 <- delete.response(Terms)
            temp <- untangle.specials(Terms2, "cluster")
            if (length(temp$vars)) 
                Terms2 <- Terms2[-temp$terms]
            temp <- untangle.specials(Terms2, "strata")
            if (length(temp$vars)) 
                Terms2 <- Terms2[-temp$terms]
            mf2 <- model.frame(Terms2, newdata, xlev = object$xlevels)
            x2 <- model.matrix(Terms2, mf2)[, -1, drop = FALSE]
            x2 <- scale(x2, center = xcenter, scale = FALSE)
            offset2 <- model.offset(mf2)
            if (length(offset2) > 0) 
                offset2 <- offset2 - mean(offset)
            else offset2 <- 0
        }
    }
    newrisk <- exp(c(x2 %*% coef) + offset2)
    if (is.factor(strata)) 
        ustrata <- levels(strata)
    else ustrata <- sort(unique(strata))
    nstrata <- length(ustrata)
    survlist <- vector("list", nstrata)
    if (se.fit) 
        varhaz <- vector("list", nstrata)
    for (i in 1:nstrata) {
        indx <- which(strata == ustrata[i])
        survlist[[i]] <- survival:::agsurv(y[indx, , drop = F], x[indx, 
            , drop = F], wt[indx], risk[indx], survtype, vartype)
    }
    if (!individual) {
        cumhaz <- unlist(lapply(survlist, function(x) x$cumhaz))
        varhaz <- unlist(lapply(survlist, function(x) cumsum(x$varhaz)))
        nevent <- unlist(lapply(survlist, function(x) x$n.event))
        ndeath <- unlist(lapply(survlist, function(x) x$ndeath))
        xbar <- t(matrix(unlist(lapply(survlist, function(x) t(x$xbar))), 
            nrow = ncol(x)))
        hazard <- unlist(lapply(survlist, function(x) x$hazard))
        if (survtype == 1) 
            surv <- unlist(lapply(survlist, function(x) cumprod(x$surv)))
        else surv <- exp(-cumhaz)
        if (is.matrix(x2) && nrow(x2) > 1) {
            surv <- outer(surv, newrisk, "^")
            varh <- matrix(0, nrow = length(varhaz), ncol = nrow(x2))
            for (i in 1:nrow(x2)) {
                dt <- outer(cumhaz, x2[i, ], "*") - xbar
                varh[, i] <- (varhaz + rowSums((dt %*% varmat) * 
                  dt)) * newrisk[i]^2
            }
        }
        else {
            surv <- surv^newrisk
            dt <- outer(cumhaz, c(x2)) - xbar
            varh <- (varhaz + rowSums((dt %*% varmat) * dt)) * 
                newrisk^2
        }
        result <- list(n = as.vector(table(strata)), time = unlist(lapply(survlist, 
            function(x) x$time)), n.risk = unlist(lapply(survlist, 
            function(x) x$n.risk)), n.event = nevent, n.censor = unlist(lapply(survlist, 
            function(x) x$n.censor)), surv = surv, type = type)
        if (nstrata > 1) {
            result$strata <- unlist(lapply(survlist, function(x) length(x$n.risk)))
            names(result$strata) <- ustrata
        }
    }
    else {
        ntarget <- nrow(x2)
        surv <- vector("list", ntarget)
        n.event <- n.risk <- n.censor <- varh1 <- varh2 <- time <- surv
        stemp <- match(strata2, ustrata)
        timeforward <- 0
        for (i in 1:ntarget) {
            slist <- survlist[[stemp[i]]]
            indx <- which(slist$time > y2[i, 1] & slist$time <= 
                y2[i, 2])
            time[[i]] <- diff(c(y2[i, 1], slist$time[indx]))
            time[[i]][1] <- time[[i]][1] + timeforward
            timeforward <- y2[i, 2] - max(slist$time[indx])
            if (survtype == 1) 
                surv[[i]] <- slist$surv[indx]^newrisk[i]
            else surv[[i]] <- slist$hazard[indx] * newrisk[i]
            n.event[[i]] <- slist$n.event[indx]
            n.risk[[i]] <- slist$n.risk[indx]
            n.censor[[i]] <- slist$n.censor[indx]
            dt <- outer(slist$cumhaz[indx], x2[i, ]) - slist$xbar[indx, 
                , drop = F]
            varh1[[i]] <- slist$varhaz[indx] * newrisk[i]^2
            varh2[[i]] <- rowSums((dt %*% varmat) * dt) * newrisk[i]^2
        }
        varh <- cumsum(unlist(varh1)) + unlist(varh2)
        if (survtype == 1) 
            surv <- cumprod(unlist(surv))
        else surv <- exp(-cumsum(unlist(surv)))
        result <- list(n = as.vector(table(strata)[stemp[1]]), 
            time = cumsum(unlist(time)), n.risk = unlist(n.risk), 
            n.event = unlist(n.event), n.censor = unlist(n.censor), 
            surv = surv, type = type)
    }
    if (!censor) {
        kfun <- function(x, keep) {
            if (is.matrix(x)) 
                x[keep, , drop = F]
            else if (length(x) == length(keep)) 
                x[keep]
            else x
        }
        keep <- (result$n.event > 0)
        if (nstrata > 1) {
            temp <- rep(ustrata, result$strata)  #each = result$strata)
            result$strata <- c(table(temp[keep]))
        }
        result <- lapply(result, kfun, keep)
        if (se.fit) 
            varh <- kfun(varh, keep)
    }
    if (se.fit) {
        result$std.err <- sqrt(varh)
        zval <- qnorm(1 - (1 - conf.int)/2, 0, 1)
        if (conf.type == "plain") {
            temp1 <- result$surv + zval * result$std.err * result$surv
            temp2 <- result$surv - zval * result$std.err * result$surv
            result <- c(result, list(upper = pmin(temp1, 1), 
                lower = pmax(temp2, 0), conf.type = "plain", 
                conf.int = conf.int))
        }
        if (conf.type == "log") {
            xx <- ifelse(result$surv == 0, 1, result$surv)
            temp1 <- ifelse(result$surv == 0, 0 * result$std.err, 
                exp(log(xx) + zval * result$std.err))
            temp2 <- ifelse(result$surv == 0, 0 * result$std.err, 
                exp(log(xx) - zval * result$std.err))
            result <- c(result, list(upper = pmin(temp1, 1), 
                lower = temp2, conf.type = "log", conf.int = conf.int))
        }
        if (conf.type == "log-log") {
            who <- (result$surv == 0 | result$surv == 1)
            xx <- ifelse(who, 0.1, result$surv)
            temp1 <- exp(-exp(log(-log(xx)) + zval * result$std.err/log(xx)))
            temp1 <- ifelse(who, result$surv + 0 * result$std.err, 
                temp1)
            temp2 <- exp(-exp(log(-log(xx)) - zval * result$std.err/log(xx)))
            temp2 <- ifelse(who, result$surv + 0 * result$std.err, 
                temp2)
            result <- c(result, list(upper = temp1, lower = temp2, 
                conf.type = "log-log", conf.int = conf.int))
        }
    }
    result$call <- Call
    if (is.R()) 
        class(result) <- c("survfit.cox", "survfit")
    else oldClass(result) <- "survfit.cox"
    result
}

