# Design  FEH 1Aug90, re-written 21Oct91
#
# Augments S formula language to include:
# 
#	name	- name[i] = name of ith original variable in x
#	label	- label[i] = label of ith original variable (=name if none)
#	assume	- assume(original x)
#	assume.code - coded version of assume (1-9, 9=added interaction)
#	parms	- parms(original x)
#		  for interaction effects parms[[i]] is a matrix with dim
#		  3 x (1+# interaction terms).  First element in pair
#		  is 1 if first factor is represented as an expanded
#		  non-linear term, 0 otherwise (this applies to polynomial,
#		  lspline, rcspline, scored).  Second element applies to
#		  second factor in interaction effect.  Third element
#		  applies to third factor (0 if second order interaction)
#		  First column contains factor numbers involved in interaction.
#	limits  - limits(original x)
#	values	- For continuous variables with <=10 unique values, is
#		  vector of values.  NULL otherwise.
#	interactions - 3 x k matrix of factor numbers
#
# Cannot have interactions between two stratification factors. 
#
#


Design <- function(mf, allow.offset=TRUE, intercept=1) {

  Terms <- attr(mf, "terms")
  Term.labels <- attr(Terms,'term.labels')
  iscluster <- if(length(Term.labels))
    substring(Term.labels,1,8)=='cluster('  else FALSE
  ## Handle cluster() for cph
  if(any(iscluster)) {
    clustername <- Term.labels[iscluster]
    cluster     <- mf[[clustername]]
    mf[[clustername]] <- NULL
    Terms       <- Terms[!iscluster]
    Term.labels <- Term.labels[!iscluster]
  }
  else cluster <- NULL

  ## Problem: offsets not included in term.labels in R

  ##For some reason, model frame sometimes has a blank name if using %ia%

  namx <- names(mf)

  if(any(namx=="")) {
    namx <- names(mf) <- c(namx[1],Term.labels)
    dimnames(mf)[[2]] <- namx
    dimnames(attr(Terms,"factors"))[[1]] <- namx
  }

  wts <- if(any(namx=='(weights)'))(1:length(namx))[namx=='(weights)']
  else 0

  if(length(Terms)==0) inner.name <- NULL
  else {
    ##Handle case where a function has two arguments that are names,
    ##e.g. rcs(x,knots) -> want x only
    inner.name <- unique(var.inner(Terms))
    ##Note: these exclude interaction terms and %ia% terms
  }
  
  response.pres <- attr(Terms, 'response') > 0
  
  offs <- attr(Terms, "offset")
  if(!length(offs)) offs <- 0
  if(offs>0 & !allow.offset)
	stop("offset variable not allowed in formula")


  factors <- attr(Terms, "factors")
  if(length(factors) && response.pres) factors <- factors[-1,,drop=FALSE]

  attr(Terms, "intercept") <- intercept
  fname <- flabel <- name <- strt <- asm <- len <- 
    fname.incl.dup <- ia <- funits <- NULL
  parm <- nonlinear <- limits <- values <- list()

  scol<-1
  colnam <- list()

  XDATADIST <- .Options$datadist
  if(length(XDATADIST)) {
    if(!exists(XDATADIST))
      stop(paste("dataset",XDATADIST,
                 "not found for options(datadist=)"))
    datadist <- eval(as.name(XDATADIST))
    Limits <- datadist$limits
    Limnames <- dimnames(Limits)[[2]]
  }
  
  nc <- 0

  options(Design.attr=NULL, TEMPORARY=FALSE)
  ##Used internally by asis, rcs, ...

  anyfactors <- ncol(mf) > 1*response.pres
  i1.noia <- 0
  if(anyfactors)for(i in (response.pres+1):ncol(mf)) {
    if(i != offs && i !=wts) {
      i1 <- i - response.pres
      xi <- mf[[i]]
      z <- attributes(xi)
      assu <- z$assume.code
      if(!length(assu) || assu!=9) i1.noia <- i1.noia+1
      if(!length(assu)) {
        ## Not processed w/asis,et
        nam <- inner.name[i1.noia]
        lab <- attr(xi, "label")
        ord <- is.ordered(xi) && all.is.numeric(levels(xi))
        if(!length(lab) || lab=="") lab <- nam
        if(ord) {
          xi <- scored(xi, name=nam, label=lab)
          attr(mf[,i],"contrasts") <- attr(xi,"contrasts")
        }
        else if(is.character(xi) | is.factor(xi)) {
          if(is.ordered(xi) &&
             .Options$contrasts[2]!='contr.treatment')
            warning(paste('Variable',nam,'is an ordered factor.\n',
                          'You should set options(contrasts=c("contr.treatment","contr.treatment"))\nor Design will not work properly.'))
          xi <- catg(xi, name=nam, label=lab)
        }
        else if(is.matrix(xi)) xi <- matrx(xi, name=nam, label=lab)
        else xi <- asis(xi, name=nam, label=lab)
        z <- c(z,attributes(xi))
      }

      za <- z$assume.code
      zname <- z$name

      fname.incl.dup <- c(fname.incl.dup, zname)
      if(!length(fname) || !any(fname==zname)) { # unique factor
        nc <- nc+1
        fname <- c(fname,zname)
        flabel <- c(flabel,z$label)
        asm <- c(asm,za)
        colnam[[i1]] <- z$colnames
        if(za != 8 && length(colnam)) name <- c(name, colnam[[i1]])
        if(za != 9) {
          funits <- c(funits, if(length(z$units))z$units else '')
          if(length(z$parms)) parm[[zname]] <- z$parms
          if(length(XDATADIST)) {
            limits[[zname]] <- if(any(Limnames==zname)) {
              j <- match(zname, Limnames, 0) #require EXACT match
              Limits[,j[j>0]]
            }
            else rep(NA,7)
            j <- match(zname, names(datadist$values), 0)
            if(j>0) {
              values[[zname]] <- datadist$values[[j]]
              l1 <- levels(xi); l2 <- datadist$values[[j]]
              if(length(l1) && ((length(l1) != length(l2)) ||
                                any(sort(l1) != sort(l2))))
                warning(paste('Variable',zname,'has levels',paste(l1,collapse=' '),
                              'which do not match levels given to datadist (',
                              paste(l2,collapse=' '),'). datadist values ignored.'))
              values[[zname]] <- l1
            }
          }
        }
        
        if(length(nonl <- z$nonlinear)) nonlinear[[zname]] <- nonl

        if(za == 9) {
          iia <- match(z$ia, fname)
          if(any(is.na(iia)))stop(paste(paste(z$ia,collapse=" "),
                                        "cannot be used in %ia% since not listed as main effect"))
          ia <- cbind(ia, c(iia,0))
          parms <- rbind(z$parms,0)
          parms[,1] <- c(iia,0)
          if(length(parms)) parm[[zname]] <- parms
        }
      }
      nrows <- if(is.matrix(xi))nrow(xi) else length(xi)
    }
  }
  
  ##Save list of which factors where %ia% interactions
  ## (before adding automatic ias)

  which.ia <- (1:length(asm))[asm==9]

  ##Add automatically created interaction terms
  
  if(anyfactors) {
    nrf <- if(!length(factors)) 0 else nrow(factors)

    ## S-Plus, if only offset in model, has factors as 2 rows 0 cols
    if(nrf || length(fname.incl.dup))
      if((nrf-(offs > 0)) != length(fname.incl.dup))
        stop("program logic error 1")
    if(length(factors)) for(i in 1:ncol(factors)) {
      f <- factors[,i]
      j <- (1:length(f))[f>0]
      nia <- length(j)
      if(nia>1) {
        fn <- fname.incl.dup[j]
        jf <- match(fn,fname.incl.dup)
        if(any(is.na(jf))) stop("program logic error 2")
        nc <- nc + 1
        asm <- c(asm,9)
        if(nia==2)ialab <- paste(fn[1],"*",fn[2])
        else if(nia==3)ialab <- paste(fn[1],"*",fn[2],"*",fn[3])
        else stop("interaction term not second or third order")
        fname <- c(fname, ialab)
        flabel <- c(flabel, ialab)
        if(sum(asm[jf]==8)>1)
          stop("cannot have interaction between two strata factors")
        nn <- list()
        for(k in 1:nia) {
          if(asm[jf[k]]==5 | asm[jf[k]]==8)
            nn[[k]] <- paste(fn[k],"=",parm[[fname[jf[k]]]][-1],sep="")
	      else if(asm[jf[k]]==7) {
            nn[[k]] <- c(fn[k],
                         paste(fn[k],"=",
                               parm[[fname[jf[k]]]][c(-1,-2)],sep=""))
          }
	      else nn[[k]] <- colnam[[jf[k]]]
        }
        if(nia==2) nn[[3]] <- ""
        parms <- jf
        if(length(jf)==2) parms <- c(parms, 0)
        nonlin <- NULL
        nl1 <- nonlinear[[fname[jf[1]]]]
        nl2 <- nonlinear[[fname[jf[2]]]]
        ##Strata factors don't have nonlinear duplicated for # levels - 1
        if(asm[jf[1]]==8)
          nl1 <- rep(FALSE, length(parm[[fname[jf[1]]]])-1)
        if(asm[jf[2]]==8)
          nl2 <- rep(FALSE, length(parm[[fname[jf[2]]]])-1)
        if(nia==2) nl3 <- FALSE
        else if(asm[jf[3]]==8)
          nl3 <- rep(FALSE, length(parm[[fname[jf[3]]]])-1)
        else nl3 <- nonlinear[[fname[jf[3]]]]
        n1 <- nn[[1]]
        n2 <- nn[[2]]
        n3 <- nn[[3]]
        ##model.matrix makes auto-products move first variable fastest, etc.
        for(j3 in 1:length(n3)) {
          for(j2 in 1:length(n2)) {
            for(j1 in 1:length(n1)) {
              parms <- cbind(parms,c(nl1[j1],nl2[j2],nl3[j3]))
              nonlin <- c(nonlin, nl1[j1] | nl2[j2] | nl3[j3])
              if(nia==2)name <- c(name, paste(n1[j1],"*",n2[j2]))
              else name <- c(name, paste(n1[j1],"*",n2[j2],"*",n3[j3]))
            }
          }
        }
        
        ##If was 2-way interaction and one of the factors was restricted %ia%,
        ##adjust indicators
        k <- match(jf, which.ia, 0)
        if(any(k > 0)) {
          if(nia==3)
            stop("cannot have 2-way interaction with an %ia% interaction")
          k <- jf[k>0]
          wparm <- parms[,1]==k; wparm[3] <- TRUE
          parms[wparm,] <- parm[[fname[k]]][1:2,,drop=FALSE]
          jf <- parms[,1]
          nonlin <- apply(parms, 2, any)[-1]
        }
        
        if(length(jf)==2) {jf <- c(jf, 0); parms[3,] <- 0}
        ia <- cbind(ia, jf)
        if(length(parms)) parm[[ialab]] <- parms
        if(length(nonlin)) nonlinear[[ialab]] <- nonlin
        
      }
    }
  }
  
  if(anyfactors) {
    if(length(XDATADIST))
      limits <- structure(limits,
                          row.names=c("Low:effect","Adjust to",
                            "High:effect",
                            "Low:prediction","High:prediction",
                            "Low","High"),
                          class="data.frame")
    ##data.frame converts variables always NA to factor!
    
    if(length(funits) != sum(asm!=9)) warning('program logic warning 1')
    else names(funits) <- fname[asm!=9]
    
    atr <- list(name=fname, label=flabel, units=funits, colnames=name,
                assume=c("asis","polynomial","lspline","rcspline","category",
                  "","scored","strata","interaction","matrix")[asm],
                assume.code=as.integer(asm), parms=parm, limits=limits,
                values=values,nonlinear=nonlinear,
                interactions=structure(ia,dimnames=NULL))
    
    nact <- attr(mf, 'na.action')
    if(length(nact) && length(nmiss <- nact$nmiss)) {
      jia <- grep('%ia%',names(nmiss), fixed=TRUE)
      if(length(jia)) nmiss <- nmiss[-jia]
      jz <- which(names(nmiss) != '(weights)' & ! grepl('offset\\(', names(nmiss)))
      if(response.pres) jz <- jz[jz > 1]
      names(nmiss)[jz] <- fname[asm != 9]
      attr(mf, 'na.action')$nmiss <- nmiss
    }
  }
  
  else atr <- list(name=NULL, assume=NULL, assume.code=NULL, parms=NULL)
  
  attr(mf, 'Design') <- atr
  attr(mf, 'terms')  <- Terms
  if(length(cluster)) attr(mf, 'cluster') <- cluster
  mf
}
