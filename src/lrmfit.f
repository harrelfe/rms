      SUBROUTINE LRMFIT(BETA,nxin,IDX,X,R,offset,U,V,LOGLIK,
     &  NOBS,NX,numy,NVI,VSQ,c,wv,deltab,pivot,opts,ftable,penalt,wt)
C
C     SUBROUTINE TO COMPUTE MAXIMUM LIKELIHOOD ESTIMATES OF THE PARAMETERS
C     OF AN ORDINAL LOGISTIC MODEL, ALONG WITH THE -2 LOG MAXIMUM LIKELIHOOD
C     (LOGLIK) AND THE -2 LOG LIKELIHOOD WITH SLOPES=0. SEE LLOGIT ROUTINE
C     FOR DESCRIPTION OF ARGUMENTS.
C       Fits coefficients for x(i,j), j=idx(1),...,idx(nxin).
C       nx is total # col for x (some col may be ignored in idx).
C       VSQ is nvi*nvi uncompressed version of V (output)
C       EPS (input) is singularity criterion, e.g. 1D-7
C       C(NVI), wv(2*nvi), and DELTAB(NVI) are scratch vectors
C       pivot(NVI) is integer logical vector
C       wt is a DOUBLE PRECISION vector of case weights
C       If V is singular, pivot(nvi) has offending column #
C       Uses all observations in analysis.
C      numy is rounded sum of weights for each y level if .NOT. normwt
C
C
C     REQUIRES SUBROUTINES LLOGIT, ainvb (this uses qr-decomposition
C       for solving system of equations using Fortran routines called
C       by S function solve).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION BETA(nvi),U(nvi),V(nvi*(nvi+1)/2),LOGLIK,
     &     C(nvi), wv(2*nvi),
     &     deltab(nvi),vsq(nvi,nvi),penalt(nvi,nvi),wt(NOBS)
      DOUBLE PRECISION opts(12),prc(4)
      DOUBLE PRECISION X(NOBS,NX),offset(NOBS)
C     Dimensions of X, penalt, wt not needed until inside LLOGIT
      INTEGER IDX(nxin), R(NOBS), ftable(501,nvi-nxin+1),
     &     numy(nvi-nxin+1), pivot(nvi)
      LOGICAL DVRG,ofpres,piter,normwt
      eps   = opts(1)
      dlike = opts(2)
      maxit = opts(3)
      piter = opts(4) .EQ. 1d0
      ofpres= opts(5) .EQ. 1d0
      normwt= opts(12) .EQ. 1d0
      opts(6)=0d0
      opts(7)=nvi
      DVRG=.FALSE.
      kint=nvi - nxin
      KINT1=KINT + 1
      nxm=nx
C     Don't mess up dimension statements in next routine
      IF(nxm .EQ. 0) nxm=1
C     
C     NEWTON-RAPHSON ITERATIONS TO SOLVE FOR MLEs
C     
      OLDL=1D30
      maxit1=maxit - 1
      if(maxit .LT. 2)GO TO 730
      DO i=1,nvi
         deltab(i)=0D0
      ENDDO
      curstp=1D0
      DO 700 ITER=1,maxit1
         CALL LLOGIT(BETA,IDX,X,R,offset,U,V,C,LOGLIK,NOBS,NOBS,
     &        nxm,ofpres,NVI,KINT,DVRG,ftable,.FALSE.,penalt,wt,normwt)
         dmax = 0D0
         DO i=1,nvi
            dmax = DMAX1(dmax, DABS(U(i)))
         ENDDO
         IF(piter) THEN
            prc(1)=LOGLIK
            prc(2)=CURSTP
            prc(3)=OLDL - LOGLIK
            prc(4)=dmax
            CALL dblepr('-2LL,Step,delta ll,max u', 24, prc, 4)
            CALL dblepr('u', 1, U, nvi)
         ENDIF
C     Added 2Jul10: was trying to step halve if initial estimates were MLEs
C     Allow for a tiny worsening of LL without step-halving if
C     Max absolute first derivative is small
         IF(dmax .LT. 1D-9 .AND. DABS(loglik - oldl) .LT. 
     &      0.1D0*dlike) GO TO 730
C         IF(ITER .EQ. 1 .AND. dmax .LT. 1D-12) GO TO 730
         IF(loglik .GT. oldl.or.dvrg) THEN
            IF(iter .EQ. 1) THEN
               opts(6)=1D0
               RETURN
            ENDIF
C     Try step-halving
            curstp=curstp/2D0
            DO j=1,nvi
               beta(j)=beta(j) - curstp*deltab(j)
            ENDDO
            GO TO 700
         ENDIF
         curstp=1D0
C     Compute V inverse * U by solving system of equations
         CALL ainvb(v, u, deltab, nvi, eps, nrank, pivot, vsq, c, wv)
         IF(nrank .LT. nvi) THEN
            opts(7)=nrank
            RETURN
         ENDIF
C     UPDATE BETA ESTIMATES
         DO 630 I=1,NVI
 630        BETA(I)=BETA(I) + deltab(I)
C     call dblepr('v',1,v,nvi*(nvi+1)/2)
C     call dblepr('deltab',6,deltab,nvi)
C     call dblepr('beta',4,beta,nvi)
C     SEE IF CONVERGENCE OBTAINED
            IF (DABS(OLDL - LOGLIK) .LE. dlike) GO TO 730
            OLDL=LOGLIK
 700     CONTINUE
         IF(maxit .GT. 2) THEN
            opts(6)=1d0
            RETURN
         ENDIF
C     COMPUTE LOGLIK AND DERIVATIVES AT LAST ESTIMATE
C     NOTE: V IS NOT INVERTED.
 730     CONTINUE
         CALL LLOGIT(BETA,IDX,X,R,offset,U,V,C,LOGLIK,NOBS,NOBS,
     &        nxm,ofpres,NVI,KINT,DVRG,ftable,.TRUE.,penalt,wt,normwt)
         CALL gcorr(ftable,kint,numy,nx,cindex,somer,gamma,taua)
         opts(8)=cindex
         opts(9)=somer
         opts(10)=gamma
         opts(11)=taua
         DO i=1,nvi
            DO j=1,nvi
               vsq(i,j)=v(isub(i,j))
            ENDDO
         ENDDO
         IF(dvrg)opts(6)=1d0
         RETURN
         END
      SUBROUTINE LLOGIT(BETA, IDX, X, R, offset, U, V, C,
     & LOGLIK, NOBS, NMAX, nxm, ofpres, NVI, KINT, DVRG,
     & ftable, calcc, penalt, wt, normwt)
C
C     ROUTINE TO CALCULATE LOGISTIC LOG LIKELIHOOD AND ITS DERIVATIVES AT BETA
C     FOR VARIABLES IDX(1),IDX(2),...,IDX(NVI-KINT).
C     ASSUMES THAT BETA IS SET
C     UP TO CONTAIN THE CURRENT PARAMETER ESTIMATES CONTIGUOUSLY,
C     WITH THE FIRST KINT PARAMETERS BEING INTERCEPTS AND THE REMAINING
C     NVI-KINT BEING SLOPES FOR THE XS.
C     THE DESIGN MATRIX IS X(NMAX,NX). THE RESPONSE VECTOR IS R(NMAX).
C     THE FIRST DERIVATIVES ARE IN U(NVI) AND THE SECOND DERIVATIVES ARE
C     IN V(NVI*(NVI+1)/2) STORED AS THE LOWER TRIANGLE.
C     C(NVI-KINT) IS A WORKING ARRAY.
C     LOGLIK IS -2 LOG LIKELIHOOD.
C     DVRG IS SET TO .TRUE. IF DIVERGENCE OCCURS.
C     Computes rank correlation measures in calcc=.TRUE.
C     Penalt is the penalty matrix (center of quadratic form)
C     Penalties for intercept terms assumed to be zero
C     Log likelihood is penalized by 0.5 times beta'*penalt*beta
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION BETA(nvi),U(nvi),V(nvi*(nvi+1)/2),
     & C(nvi), LOGLIK, penalt(nvi,nvi), wt(nmax)
      DOUBLE PRECISION X(NMAX,nxm),offset(nmax)
      INTEGER IDX(nxm),R(NMAX),ftable(501,KINT+1)
      LOGICAL DVRG,ofpres,calcc,normwt
      PROB(BB)=1D0/(1D0 + DEXP(-DMIN1(DMAX1(BB,-30D0),30D0)))
      NV=NVI
      BX=0
      KINT1=KINT + 1
      NVIX=NVI - KINT
      DVRG=.FALSE.
      LOGLIK=0D0
      OFF=0D0
      incobs=1
      DO I=1,NV
         IF(i .LE. kint) THEN
            ii=i
         ELSE
            ii=idx(i - kint) + kint
         ENDIF
         U(I)=0D0
         DO j=1,NV
            IF(j .LE. kint) THEN
               jj=j
            ELSE
               jj=idx(j - kint) + kint
            ENDIF
            U(I)=U(I) - penalt(ii,jj)*beta(j)
            v(isub(i,j))=penalt(ii,jj)
         ENDDO
         loglik=loglik - 0.5D0*penalt(ii,ii)*(beta(i)**2)
         IF(nv .GT. 1) THEN
            DO j=(i + 1),nv
               IF(j .LE. kint) THEN
                  jj=j
               ELSE
                  jj=idx(j - kint) + kint
               ENDIF
               loglik=loglik - penalt(ii,jj)*beta(i)*beta(j)
            ENDDO
         ENDIF
      ENDDO
      NE=NV*(NV + 1)/2
C
C       Get middle category for predicted prob.
C
      IF(calcc) THEN
         mid=kint1 / 2
         DO i=1, 501
            DO j=1, kint1
               ftable(i,j)=0
            ENDDO
         ENDDO
      ENDIF
C
C     For efficiency handle special case KINT=1
C
      IF(kint .EQ. 1) THEN
         DO n=1,nobs
            iy=r(n)
            IF(ofpres)off=offset(n)
            w=wt(n)
            bx=0D0
            DO i=1,nv
               IF(i .EQ. 1) THEN
                  c(i)=1D0
               ELSE
                  c(i)=x(n,idx(i - 1))
               ENDIF
               bx=bx + beta(i)*c(i)
            ENDDO
            cpiy=PROB(bx + off)
            IF(calcc) THEN
               ipp=500d0*cpiy + 1d0
               if(.NOT. normwt) incobs=w + .5D0
               ftable(ipp,iy + 1)=ftable(ipp,iy + 1) + incobs
            ENDIF
            viy=cpiy*(1D0 - cpiy)
            IF(iy .EQ. 0) THEN
               IF(cpiy.ge.1d0) GO TO 74
               loglik=loglik + DLOG(1D0 - cpiy)*w
               dy=0D0
            ELSE
               IF(cpiy.le.0d0) GO TO 74
               loglik=loglik + DLOG(cpiy)*w
               dy=1D0
            ENDIF
C     Add to first and second derivatives
            k=0
            DO i=1,nv
               xi=c(i)
               u(i)=u(i) + (dy - cpiy)*xi*w
               DO j=1,i
                  k=k + 1
                  v(k)=v(k) + viy*xi*c(j)*w
               ENDDO
            ENDDO
         ENDDO
         GO TO 210
      ENDIF     
      DO 200 N=1,NOBS
         iy=r(n)
         w=wt(n)
         IF(ofpres)OFF=offset(N)
         IF(nvix .GT. 0) THEN
            bx=0D0
            DO I=1,NVIX
               C(I)=X(N,IDX(I))
               bx=bx + beta(kint + i)*c(i)
            ENDDO
         ENDIF
C     COMPUTE EXCEEDENCE PROBABILITIES AND CHECK FOR DIVERGENCE
         IF(calcc) THEN
            ipp=500d0*prob(bx + beta(mid) + off) + 1d0
            if(.NOT. normwt)incobs=w + .5D0
            ftable(ipp,iy + 1)=ftable(ipp,iy + 1) + incobs
         ENDIF
         IF(iy.NE.0) THEN
            CPIY=PROB(bx + beta(iy) + off)
            VIY=CPIY*(1D0 - CPIY)
            IF(IY .EQ. KINT)GO TO 78
         ENDIF
         CPIY1=PROB(BX + BETA(IY + 1) + off)
         VIY1=CPIY1*(1D0 - CPIY1)
         IF(CPIY1 .LT. 1D0)GO TO 80
         GO TO 74
 78      IF(CPIY .EQ. 0D0)GO TO 74
C     ADD TO LOG-LIKELIHOOD
 80      IF(IY .EQ. 0)LOGLIK=LOGLIK + DLOG(1D0 - CPIY1)*w
         IF(IY .EQ. KINT)LOGLIK=LOGLIK + DLOG(CPIY)*w
         IF(IY .GT. 0 .AND. IY .LT. KINT)LOGLIK=LOGLIK +
     &      DLOG(CPIY - CPIY1)*w
C     ADD TO FIRST AND SECOND DERIVATIVES
         K=0
         DO 180 I=1,NV
            XI=1D0
            IF(I .GT. KINT)XI=C(I - KINT)
C     COMPUTE DERIVATIVES OF ALPHA(IY) + XBETA WRT XI AND
C     DERIV OF ALPHA(IY + 1) + XBETA WRT XI
            DIIY=0D0
            IF(I .EQ. IY .OR. I .GT. KINT)DIIY=XI
            DIIY1=0D0
            IF(I .EQ. IY + 1 .OR. I .GT. KINT)DIIY1=XI
            IF(IY .GT. 0)GO TO 110
            U(I)=U(I) - CPIY1*DIIY1*w
            GO TO 130
 110        IF(IY .LT. KINT)GO TO 120
            U(I)=U(I) + (1D0 - CPIY)*DIIY*w
            GO TO 130
C     0<IY<KINT
 120        PIY=CPIY - CPIY1
            IF(VIY .GT. 0D0 .OR. VIY1 .GT. 0D0)
     &        U(I)=U(I) + (VIY*DIIY - VIY1*DIIY1)/PIY * w
 130        DO 160 J=1,I
               K=K + 1
               XJ=1D0
               IF(J .GT. KINT)XJ=C(J - KINT)
C     COMPUTE DERIVATIVE OF ALPHA(IY) + XBETA WRT XJ AND
C     DERIV OF ALPHA(IY + 1) + XBETA WRT XJ
               DJIY=0D0
               IF(J .EQ. IY .OR. J .GT. KINT)DJIY=XJ
               DJIY1=0D0
               IF(J .EQ. IY + 1 .OR. J .GT. KINT)DJIY1=XJ
               IF(IY .GT. 0)GO TO 140
               V(K)=V(K) + VIY1*DIIY1*DJIY1 * w
               GO TO 160
 140           IF(IY .LT. KINT)GO TO 150
               V(K)=V(K) + VIY*DIIY*DJIY * w
               GO TO 160
C     0<IY<KINT
 150  IF(VIY .GT. 0D0 .OR. VIY1 .GT. 0D0)V(K)=V(K) - (
     & PIY*(VIY*(1D0 - 2D0*CPIY)*DIIY*DJIY
     &     -VIY1*(1D0 - 2D0*CPIY1)*DIIY1*DJIY1
     & ) - (VIY*DIIY - VIY1*DIIY1)*(VIY*DJIY - VIY1*DJIY1))/PIY/PIY * w
 160        CONTINUE
 180     CONTINUE
 200  CONTINUE
 210  LOGLIK=-2D0*LOGLIK
      RETURN
 74   dvrg=.TRUE.
      RETURN
      END
      SUBROUTINE gcorr(ftable, kint, numy, nx, c, somer,
     &     gamma, taua)
      INTEGER ftable(501, kint + 1), numy(kint + 1)
      DOUBLE PRECISION fn, c, somer, gamma, taua, con, dis, tie,
     &   k, ic, nn
      kint1=kint + 1
      con=0d0
      dis=0d0
      tie=0d0
      c=.5d0
      somer=0d0
      gamma=0d0
      taua=0d0
      if(nx .EQ. 0)RETURN
      nn=0
      DO i1=1,kint1
         nn=nn + numy(i1)
      ENDDO
      DO i1=1,kint
         i1p1=i1 + 1
         DO j1=1,501
            k=ftable(j1,i1)
            IF(k .GT. 0d0) THEN
               j1p1=j1 + 1
               DO i2=i1p1,kint1
                  ic=0d0
                  if(j1.lt.501) then
                     DO j2=j1p1,501
                        ic=ic + ftable(j2,i2)
                     ENDDO
                  ENDIF
                  con=con + k*ic
                  dis=dis + k*(numy(i2) - ic - ftable(j1,i2))
                  tie=tie + k*ftable(j1,i2)
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      fn=con + dis + tie
      c=(con + tie*.5d0)/fn
      somer=(con - dis)/fn
      gamma=0d0
      if(con + dis .GT. 0) gamma=(con - dis)/(con + dis)
      taua=(con - dis)/(nn*(nn - 1d0)/2D0)
      RETURN
      END
