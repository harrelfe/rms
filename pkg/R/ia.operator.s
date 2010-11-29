#ia.operator.s - restricted interaction operators for use with rms
#F. Harrell  8 Nov 91

#Set up proper attributes for a restricted interaction for a model
#such as y ~ rcs(x1) + rcs(x2) + x1 %ia% x2 or x1 %ia% rcs(x2)
#or rcs(x1) %ia% x2

"%ia%" <- function(x1, x2)
{
  a1 <- attributes(x1)
  a2 <- attributes(x2)
  nam <- as.character(sys.call())[-1]
  
  redo <- function(x,nam)
    {
      if(is.null(attr(x,"assume.code")))
        {
          if(!is.null(oldClass(x)) && oldClass(x)[1]=="ordered")
            x <- scored(x, name=nam)
          else if(is.character(x) | is.category(x)) x <- catg(x, name=nam)
          else if(is.matrix(x)) x <- matrx(x, name=nam)
          else x <- asis(x, name=nam)
        }
      ass <- attr(x,"assume.code")
      nam <- attr(x,"name")
      
      if(ass==5)
        {
          colnames <- attr(x,"colnames")
          len <- length(attr(x,"parms"))-1	}
      else
        if(ass==8)
          {
            prm <- attr(x,"parms")
            colnames <- paste(nam,"=",prm[-1],sep="")
            len <- length(prm)-1
          }
        else if(ass==7)
          {
            prm <- attr(x,"parms")
            colnames <- c(nam,paste(nam,"=",prm[-(1:2)],sep=""))
            len <- length(prm)-1
          }
        else
          {
            if(is.null(ncol(x)))
              {
                len <- 1
                colnames <- nam		}
            else
              {
                colnames <- dimnames(x)[[2]]
                len <- ncol(x)
              }
          }
      
      attr(x,"colnames") <- colnames
      attr(x,"len") <- len
      if(ass==8) attr(x,"nonlinear") <- rep(FALSE, len)
      x
    }

  x1 <- redo(x1,nam[1])
  x2 <- redo(x2,nam[2])
  a1 <- attributes(x1)
  a2 <- attributes(x2)
  n1 <- a1$colnames
  n2 <- a2$colnames
  nl1 <- a1$nonlinear
  nl2 <- a2$nonlinear
  as1 <- a1$assume.code
  as2 <- a2$assume.code

  l1 <- a1$len
  l2 <- a2$len
  if(any(nl1) & any(nl2))	nc <- l1+l2-1   
  else nc <- l1*l2
  if(is.matrix(x1)) nr <- nrow(x1) 
  else nr <- length(x1)
  x <- matrix(single(1),nrow=nr,ncol=nc)
  name <- character(nc)
  parms <- matrix(integer(1),nrow=2,ncol=nc+1)
  nonlinear <- logical(nc)
  
  k <- 0
  if(!is.factor(x1)) x1 <- as.matrix(x1)
  if(!is.factor(x2)) x2 <- as.matrix(x2)
  for(i in 1:l1)
    {
      if(as1==5 | as1==8) x1i <- oldUnclass(x1)==(i+1)
      else x1i <- x1[,i]
      for(j in 1:l2)
        {
          ##Remove doubly nonlinear terms
          if(nl1[i] & nl2[j]) break
          k <- k + 1
          if(as2==5 | as2==8) x2j <- oldUnclass(x2)==(j+1)
          else x2j <- x2[,j]
          x[,k] <- x1i * x2j
          name[k] <- paste(n1[i],"*",n2[j])
          parms[,k+1] <- c(nl1[i],nl2[j])
          nonlinear[k] <- nl1[i] | nl2[j]
        }
    }

  dimnames(x) <- list(NULL, name)
  attr(x,"ia") <- c(a1$name, a2$name)
  attr(x,"parms") <- parms
  attr(x,"nonlinear") <- nonlinear
  attr(x,"assume.code") <- 9
  attr(x,"name") <- paste(a1$name,"*",a2$name)
  attr(x,"label") <- attr(x,"name")
  attr(x,"colnames") <- name
  attr(x,"class") <- "rms"
  x
}
