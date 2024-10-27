! Converted from Ratfor 2024-10-26 using ChatGPT and the following request:
! Convert the following Ratfor code to Fortran 2018 using Fortran ISO environment without using module

SUBROUTINE ormuv(n, p, kint, nx, x, y, pr, fpa, fpb, fppa, fppb, u, v, ja, ia, l, lia, kk)
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPLICIT NONE
  
    INTEGER, INTENT(IN) :: n, p, kint, nx, l, lia
    INTEGER, INTENT(IN) :: y(n)
    INTEGER, INTENT(OUT) :: ja(l), ia(lia), kk(p)
    REAL(DP), INTENT(IN) :: x(n, nx), pr(n), fpa(n), fpb(n), fppa(n), fppb(n)
    REAL(DP), INTENT(OUT) :: u(p), v(l)
  
    INTEGER :: i, j, k, m, z, iv, nkk, mi, ki, ik, mm
    REAL(DP) :: uk, vmk, a, pa, pb, ppa, ppb, w
  
    ! Initialize the score vector u to zero
    u = 0.0_DP
  
    ! Calculate score vector
    DO k = 1, kint
      uk = 0.0_DP
      DO j = 1, n
        z = y(j)
        IF (z - 1 == k) THEN
          uk = uk + fpa(j) / pr(j)
        ELSE IF (z == k) THEN
          uk = uk - fpb(j) / pr(j)
        END IF
      END DO
      u(k) = uk
    END DO
  
    IF (nx > 0) THEN
      DO k = (kint + 1), p
        u(k) = SUM((fpa(:) - fpb(:)) * x(:, k - kint) / pr(:))
      END DO
    END IF
  
    ! Initialize the information matrix index and fill with zeros
    iv = 0
    v(:) = 0.0_DP
  
    DO m = 1, p
      ! Define kk array based on conditions
      IF (kint > 1) THEN
        IF (m == 1) THEN
          nkk = 2
          kk(1) = 1
          kk(2) = 2
        ELSE IF (m > 1 .AND. m < kint) THEN
          nkk = 3
          kk(1) = m - 1
          kk(2) = m
          kk(3) = m + 1
        ELSE IF (m == kint) THEN
          nkk = 2
          kk(1) = m - 1
          kk(2) = m
        ELSE
          nkk = kint
          kk(1:kint) = [(i, i = 1, kint)]
        END IF
        DO mm = (kint + 1), p
          nkk = nkk + 1
          kk(nkk) = mm
        END DO
      ELSE
        nkk = p
        kk(1:p) = [(i, i = 1, p)]
      END IF
  
      ! Calculate elements of the information matrix
      DO ik = 1, nkk
        k = kk(ik)
        vmk = 0.0_DP
  
        DO j = 1, n
          z = y(j)
          pa = fpa(j)
          pb = fpb(j)
          ppa = fppa(j)
          ppb = fppb(j)
          w = 1.0_DP / (pr(j) * pr(j))
  
          IF (m <= kint .AND. k <= kint) THEN
            a = -w * (pa * ld(z - 1 == m) - pb * ld(z == m)) * &
                (pa * ld(z - 1 == k) - pb * ld(z == k)) + &
                (ppa * ld(z - 1 == m) * ld(m == k) - ppb * ld(z == m) * ld(m == k)) / pr(j)
  
          ELSE IF (m > kint .AND. k > kint) THEN
            a = x(j, m - kint) * x(j, k - kint) / pr(j) * &
                (-1.0_DP / pr(j) * (pa - pb) * (pa - pb) + ppa - ppb)
  
          ELSE
            mi = MAX(m, k)
            ki = MIN(m, k)
            a = x(j, mi - kint) / pr(j) * &
                (-1.0_DP / pr(j) * (pa - pb) * (pa * ld(z - 1 == ki) - pb * ld(z == ki)) + &
                ppa * ld(z - 1 == ki) - ppb * ld(z == ki))
          END IF
  
          vmk = vmk + a
        END DO
  
        iv = iv + 1
        v(iv) = -vmk
        IF (kint > 1) THEN
          ja(iv) = k
          IF (ik == 1) ia(m) = iv
        END IF
      END DO
    END DO
  
    IF (kint > 1) ia(p + 1) = iv + 1

    ! Internal function for logical decision
    CONTAINS
      REAL(DP) FUNCTION ld(a)
        LOGICAL, INTENT(IN) :: a
        ld = MERGE(1.0_DP, 0.0_DP, a)
      END FUNCTION ld
  

END SUBROUTINE ormuv
  