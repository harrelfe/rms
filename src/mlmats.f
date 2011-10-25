      FUNCTION isub(i,j)
C-----------------------------------------------------------------------------
C     Computes subscript in lower triangular matrix corresponding to (i,j)
C-----------------------------------------------------------------------------
      INTEGER i,j,isub
      IF(i-j)10,10,20
10    isub=i+j*(j-1)/2
      RETURN
20    isub=j+i*(i-1)/2
      RETURN
      END
        SUBROUTINE sqtria(vsq,vtri,n,k)
C----------------------------------------------------------------------------
C       k=1 : converts n x n square symmetric matrix vsq to lower triangular
C               form and stores result in vtri
C       k=2 : converts lower triangular matrix vtri to n x n uncompressed
C               square matrix
C       F. Harrell 6Sep90
C----------------------------------------------------------------------------
        DOUBLE PRECISION vsq(n,n),vtri(*)
        IF(k.EQ.1) THEN
                l=0
                DO i=1,n
                DO j=1,i
                l=l+1
                vtri(l)=vsq(i,j)
                END DO
                END DO
          ELSE
                DO i=1,n
                DO j=1,n
                vsq(i,j)=vtri(isub(i,j))
                END DO
                END DO
                ENDIF
        RETURN
        END
      subroutine inner(b,x,n,z)
C-----------------------------------------------------------------------------
C     Computes dot product of b and x, each of length n, returns result in z
C-----------------------------------------------------------------------------
      DOUBLE PRECISION b(1),x(1),z
      z=0D0
        DO i=1,n
        z=z+b(i)*x(i)
        end do
      return
      end
      SUBROUTINE SPROD(M,V,P,N)
C-----------------------------------------------------------------------------
C     MULTIPLIES N*N SYMMETRIC MATRIX M STORED IN COMPRESSED FORMAT BY
C     THE N*1 VECTOR V AND RETURNS THE N*1 VECTOR PRODUCT P
C-----------------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION M(1),V(1),P(1)
      DO 20 I=1,N
      PI=0D0
      II=I*(I-1)/2
      DO 10 J=1,N
      IF(I-J)2,4,4
    2 IR=I+J*(J-1)/2
      GO TO 10
    4 IR=J+II
   10 PI=PI+M(IR)*V(J)
      P(I)=PI
   20 CONTINUE
      RETURN
      END
      SUBROUTINE AVA(A,V,P,N)
C-----------------------------------------------------------------------------
C     V IS AN N X N SYMMETRIC MATRIX AND A IS AN N X 1 VECTOR.
C     THIS ROUTINE RETURNS P=A'VA
C-----------------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION  A(1),V(1)
      P=0D0
      K=0
      DO 10 I=1,N
         AI=A(I)
         DO 20 J=1,I
            K=K+1
            IF (I.EQ.J) THEN
               P=P+AI*AI*V(K)
            ELSE
               P=P+2D0*AI*A(J)*V(K)
            ENDIF
20       CONTINUE
10    CONTINUE
      RETURN
      END
      SUBROUTINE avia(a,v,p,n,idx,nidx,nrank,eps,vsub,wv1,wv2,wv3,
     &                  wv4,pivot)
C----------------------------------------------------------------------------
C       V is an n x n symmetric matrix and a is an n x 1 vector.
C       Returns P=a' v**-1 a and nrank=rank(v), where 
C       a=a(idx(i),i=1,...,nidx), v=(v(idx(i),idx(i),i=1,...,nidx).
C       vsub is nidx x nidx scratch matrix and wv1-wv4 are scratch
C       vectors of length nidx (except for wv3 which is 2*nidx).  
C   pivot is scratch integer vector
C       of length nidx.  eps is singularity criterion, e.g. 1d-7.
C       Uses Fortran routines dqr (see S function qr) and dqrsl1 
C   (see S function solve).  In R these are dqrdc2 and dqrsl (args
C      differ too).
C
C       F. Harrell 20 Nov 90
C
C----------------------------------------------------------------------------
        DOUBLE PRECISION a(n),wv1(nidx),wv2(nidx),wv3(*),wv4(nidx),
     &          v(n,n),eps,vsub(nidx,nidx),p
        INTEGER idx(nidx),pivot(nidx),dim(2)
        k=nidx
C       CALL intpr("k",1,k,1)
        dim(1)=k
        dim(2)=k
                DO i=1,k
                wv4(i)=a(idx(i))
                pivot(i)=i
                        DO j=1,k
                        vsub(i,j)=v(idx(i),idx(j))
                        ENDDO
                ENDDO
C       CALL dblepr('wv4',3,wv4,k)
C       CALL dblepr('vsub',4,vsub,k*k)
        nrank=k
C        CALL dqr(vsub,dim,pivot,wv2,eps,wv3,nrank)
        CALL dqrdc2(vsub,dim,dim,dim,eps,nrank,wv2,pivot,wv3)
C       CALL intpr('nrank',5,nrank,1)
        IF(nrank.LT.k)RETURN
                DO i=1,k
                wv3(i)=wv4(i)
                ENDDO
        j=1
        i=100
C        CALL dqrsl1(vsub,dim,wv2,nrank,wv4,1,wv3,wv1,i,j)
        CALL dqrsl(vsub,dim,dim,nrank,wv2,wv4,wv3,wv1,wv1,
     &             wv3,wv3,i,j)
        p=0d0
                DO i=1,k
                p=p+wv4(i)*wv1(i)
                ENDDO
C       CALL intpr('dim',3,dim,2)
C       CALL dblepr('vsub',4,vsub,k*k)
C       CALL dblepr('wv1',3,wv1,k)
C       CALL dblepr('wv4',3,wv4,k)
C       CALL dblepr('p',1,p,1)
        RETURN
        END
      SUBROUTINE AVIA2(A,V,P,N,idx,nidx,nrank,eps,vsub,s,swept)
C----------------------------------------------------------------------------
C     V IS AN N X N SYMMETRIC square MATRIX AND A IS AN 
C     N X 1 VECTOR.
C     THIS ROUTINE RETURNS P=a' vinverse a and nrank=rank(v) where
C     a=A(idx(i),i=1,...,nidx), v=V(idx(i),idx(i),i=1,...,nidx).
C     S(nidx) is DOUBLE PRECISION scratch vector, SWEPT(nidx) is LOGICAL scratch 
C     vector, VSUB(nidx*(nidx+1)/2) is DOUBLE PRECISION scratch vector
C     eps is singularity criterion, e.g. 1D-6
C
C     F. Harrell 6 Sep90
C----------------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION  A(n),V(n,n),s(nidx),vsub(*)
      INTEGER idx(nidx)
      LOGICAL swept(nidx)
        l=0
                DO i=1,nidx
                swept(i)=.FALSE.
                idxi=idx(i)
C       Initialize s vector to diagonal elements
                s(i)=v(idxi,idxi)
                DO j=1,i
                l=l+1
                vsub(l)=v(idxi,idx(j))
                END DO
                END DO
        nrank=0
                DO i=1,nidx
                CALL GSWEEP(s,vsub,i,lsing,nidx,eps,swept,ifault)
                IF(lsing.EQ.0)nrank=nrank+1
                ENDDO
      P=0D0
      K=0
      DO 10 I=1,Nidx
C       Singularities are like parameter never appeared
        IF(swept(i)) THEN 
                AI=A(idx(i))
            ELSE
                AI=0D0
                ENDIF
         DO 20 J=1,I
            K=K+1
            IF (I.EQ.J) THEN
               P=P+AI*AI*Vsub(K)
            ELSE
               P=P+2D0*AI*A(idx(J))*Vsub(K)
            ENDIF
20       CONTINUE
10    CONTINUE
C       gsweep returns negative of inverse
        P=-P
      RETURN
      END
        SUBROUTINE ainvb(a, b, aib, k, tol, irank, pivot, 
     &                   wv1, wv2, wv3)
C-----------------------------------------------------------------------
C       Uses same Fortran subroutines as S function solve to accurately
C       compute aib=a inverse * b, for k x k symmetric matrix a stored in
C       lower triangular form and k x 1 vector b.  wv1(k,k), wv2(k), wv3(2*k)
C       are DOUBLE PRECISION scratch arrays and pivot(k) is INTEGER scratch vector. 
C       tol is tolerance, e.g. 1d-7.
C       IF irank (output) < k, result is not computed.  Index of singular
C       column will be stored in pivot(k) if irank<k.
C-----------------------------------------------------------------------
        DOUBLE PRECISION a(1),b(k),aib(k),tol,wv1(k,k),wv2(k),wv3(*)
        INTEGER pivot(k),dim(2)
        CALL sqtria(wv1, a, k, 2)
        dim(1)=k
        dim(2)=k
        DO i=1,k
                pivot(i)=i
                ENDDO
        irank=k
C        CALL dqr(wv1, dim, pivot, wv2, tol, wv3, irank)
        CALL dqrdc2(wv1,dim,dim,dim,tol,irank,wv2,pivot,wv3)
C       CALL dblepr('wv1',3,wv1,10)
C       CALL intpr('pivot',5,pivot,k)
C       CALL dblepr('wv2',3,wv2,k)
C       CALL dblepr('wv3',3,wv3,k)
C       CALL intpr('irank',5,irank,1)
        IF(irank.LT.k)RETURN
        DO i=1,k
                wv3(i)=b(i)
                ENDDO
        j=1
        i=100
C        CALL dqrsl1(wv1, dim, wv2, irank, b, 1, wv3, aib, i, j)
        CALL dqrsl(wv1,dim,dim,irank,wv2,b,wv3,aib,aib,
     &             wv3,wv3,i,j)
C       CALL intpr('wv1',3,wv1,10)
C       CALL intpr('dim',3,dim,2)
C       CALL dblepr('wv2',3,wv2,k)
C       CALL intpr('irank',5,irank,1)
C       CALL dblepr('b',1,b,k)
C       CALL dblepr('wv3',3,wv3,k)
C       CALL dblepr('aib',3,aib,k)
C       CALL intpr('i',1,i,1)
C       CALL intpr('j',1,j,1)
        RETURN
        END
        SUBROUTINE matinv(x, n, ne, idx, swept, lswept, t, s, nrank, 
     &          eps,negate)
C-----------------------------------------------------------------------
C       Uses subroutine GINV to invert n*n symmetric matrix x stored in
C       full form, and returns the result in x and rank of the matrix
C       in nrank.  If collinearities are detected, nrank will be <n
C       and the appropriate row and column of x are set to 0D0.
C       t(n(n+1)/2) (DOUBLE PRECISION), s(n) (DOUBLE PRECISION) are scratch areas.
C       Inversion is done with respect to elements idx(1)...idx(ne).
C       Swept should be initialized to all .FALSE. before the first call
C       to matinv for x.  Swept and lswept are n-element LOGICAL
C       Eps is singularity criterion, e.g. .0001
C       Negate is LOGICAL - .TRUE. to negate inverted part so that it
C       is immediately usable.  Use .FALSE. if further inversions will
C       be done.
C-----------------------------------------------------------------------
        IMPLICIT DOUBLE PRECISION (a-h, o-z)
        DOUBLE PRECISION x(n,n), t(1), s(1), eps
        INTEGER idx(ne)
        LOGICAL swept(1),lswept(1),negate
        LOGICAL logeq
        k=0
         DO i=1,n
         lswept(i)=swept(i)
          DO j=1,i
          k=k+1
          t(k)=x(i,j)
          END DO
         END DO
        CALL ginv(t,s,l,ne,idx,n,eps,negate,swept,nrank,ifault)
         DO i=1,n
          DO j=1,n
          x(i,j)=t(isub(i,j))
          ENDDO
         ENDDO
         DO i=1,ne
          DO j=1,ne
          ie=idx(i)
          je=idx(j)
          IF(logeq(swept(ie),lswept(ie)).or.logeq(swept(je),
     &       lswept(je)))
     &          x(ie,je)=0d0
          END DO
         END DO
        RETURN
        END
        LOGICAL FUNCTION logeq(a,b)
        LOGICAL a,b,aa,bb
        INTEGER ia,ib
        EQUIVALENCE (aa,ia),(bb,ib)
        aa=a
        bb=b
        logeq=ia.eq.ib
        return
        end
      SUBROUTINE GINV(T, S, L, ne, idx, N, E, NEG, SWEPT, NRANK, IFAULT)
C
C     USES SUBROUTINE GSWEEP TO CALCULATE THE INVERSE OF AN NXN SYMMETRIC
C     MATRIX T STORED LOWER TRIANGULAR ROW BY ROW, SHORTEST ROW FIRST, IN
C     LOCATIONS 1 TO N*(N+1)/2.  Inverse is taken with respect to diagonal
C     pivots idx(1),...,idx(ne).  For total inverse set idx(1)=1,...
C     idx(n)=n.  Set idx(1)=0 (or use scalar idx) to pretend that
C     idx(1)=1,...,idx(ne)=ne were specified.
C
C     S IS A SCRATCH ARRAY OF LENGTH N.
C     L (OUTPUT) IS 0 IF THE MATRIX WAS INVERTED SUCCESSFULLY. IF IT WAS
C     NOT, L IS SET TO THE INDEX OF THE FIRST ROW/COLUMN THAT WAS
C     COLLINEAR.
C     E (INPUT) IS A SMALL POSITIVE CONSTANT SPECIFYING THE LEVEL AT
C     WHICH A MULTIPLE CORRELATION COEFFICIENT IS ASSUMED TO BE CLOSE
C     ENOUGH TO UNITY TO INDICATE COLLINEARITY.
C     SWEPT (INPUT AND OUTPUT) is a LOGICAL array with N elements.
C     Before inverting, it should be set to .FALSE., and should not
C     be otherwise changed by the caller.
C     NRANK is the rank of the portion of the matrix inverted
C     IFAULT (OUTPUT) IS SET TO 1 IF AN ARGUMENT IS IMPROPER, 0 OTHERWISE.
C     When a singularity is found, SWEPT is left with its original value.
C     When you are finished inverting all necessary parts of the matrix,
C     you must negate all the pertinent elements, unless NEG is .true.
C     This causes all elements for the idx portion to be negated to be
C     in usual form (this assumes you are through calling ginv for the
C     matrix).
C
      DOUBLE PRECISION T(*),S(N),E
      INTEGER idx(ne)
      LOGICAL SWEPT(N)
      LOGICAL NEG,trick
      DATA ZERO/0D0/
      trick=idx(1).EQ.0
      L=0
      IFAULT=1
      IF (N.LT.1.OR.E.LT.ZERO) RETURN
      IFAULT=0
      J=0
      DO 10 I=1,N
      J=J+I
   10 S(I)=T(J)
C
      NRANK=0
      DO j=1,ne
        IF(trick) THEN
                je=j
             ELSE
                je=idx(j)
          ENDIF
        CALL GSWEEP(S, T, je, LL, N, E, SWEPT, IFAULT)
C       IF (IFAULT.NE.0) RETURN
        IF(LL.EQ.0)NRANK=NRANK+1
        IF (LL.GT.0.AND.L.EQ.0)L=LL
        ENDDO
      IF(NEG)THEN
        DO i=1,ne
        IF(trick) THEN
                ie=i
             ELSE
                ie=idx(i)
           ENDIF
         DO j=i,ne
         IF(trick) THEN
                je=j
             ELSE
                je=idx(j)
           ENDIF
         t(isub(ie,je))=-t(isub(ie,je))
         end do
        end do
       end if
      RETURN
      END
      SUBROUTINE GSWEEP(S, T, K, L, N, E, SWEPT, IFAULT)
C
C     Clark: ALGORITHM AS 178  APPLIED STATISTICS (1982) VOL. 31, NO. 2 
C     Improvements by Ridout and Cobby, Applied Statistics 1989 AS R78
C
C     PERFORMS GAUSS-JORDAN PIVOT FOR ROW/COL K IN NXN WORKING ARRAY
C     STORED LOWER TRIANGLE ONLY ROW BY ROW, SHORTEST ROW FIRST, IN
C     LOCATIONS 1 TO N(N+1)/2 OF T.
C     S should be set to diagonal elements of T on input.
C     Modified F. Harrell 25Sep90 to allow E=0 to effectively turn
C     off singularity checking.
C
      DOUBLE PRECISION A,B,E,S(N),T(*),ZERO,ONE
      LOGICAL SWEPT(N)
      DATA ZERO,ONE /0.0D0,1.0E0/
C
      IFAULT=1
      IF (N.LT.1.OR.K.LT.1.OR.K.GT.N.OR.E.LT.ZERO)RETURN
      IFAULT=0
C
C     PARAMETERS IN RANGE SO TEST FOR COLLINEARITY
C
      L=K
      KK=K*(K+1)/2
      IF(SWEPT(K).AND.T(KK).LT.ZERO) GOTO 20
      IF(SWEPT(K).AND.T(KK).GT.ZERO)GOTO 95
      IF(T(KK).LT.ZERO)GOTO 95
C       Following was .LE. - changed FEH 25Sep90
      IF(T(KK).LT.E*S(K))RETURN
      II=0
      IK=KK-K
      DO 10 L=1,N
      II=II+L
      IK=IK+1
      IF (L.GT.K) IK=IK+L-2
      IF(.NOT.SWEPT(L).AND.T(II).GE.ZERO)GOTO 10
      IF(.NOT.SWEPT(L))GOTO 95
      IF(T(II).GT.ZERO)GOTO 95
      IF (ONE/ (T(IK) * T(IK) / T(KK) - T(II)) .LT. E*S(L)) RETURN
   10 CONTINUE
C
C     NO COLLINEARITY SO UPDATE TRIANGLE
C
   20 L=0
      T(KK)=-ONE/T(KK)
      A=DABS(T(KK))
      IK=KK-K
      IJ=0
      DO 90 I=1,N
      IK=IK+1
      IF (I-K) 50,30,40
   30 IJ=IJ+K
      GO TO 90
   40 IK=IK+I-2
   50 B=T(IK)
      IF (T(KK).LT.ZERO) B=-B
      T(IK)=A*T(IK)
      JK=KK-K
      DO 80 J=1,I
      IJ=IJ+1
      JK=JK+1
      IF (J-K) 70,80,60
   60 JK=JK+J-2
   70 T(IJ)=T(IJ)+B*T(JK)
   80 CONTINUE
   90 CONTINUE
      SWEPT(K)=.NOT.SWEPT(K)
      RETURN
95    IFAULT=2
      RETURN
      END

