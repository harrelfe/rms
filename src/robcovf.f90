! Converted from Ratfor 2024-10-26 using ChatGPT and the following request:
! Convert the following Ratfor code to Fortran 2018 using Fortran ISO environment without using module

SUBROUTINE robcovf(n, p, nc, start, len, u, s, w)
    ! Use the ISO Fortran environment
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPLICIT NONE
  
    ! Arguments
    INTEGER, INTENT(IN) :: n, p, nc
    INTEGER, INTENT(IN) :: start(nc), len(nc)
    REAL(DP), INTENT(IN) :: u(n, p)
    REAL(DP), INTENT(OUT) :: s(p), w(p, p)
  
    ! Local variables
    INTEGER :: i, j, k
  
    ! Initialize w to zero
    w = 0.0_DP
  
    ! Loop over clusters
    DO k = 1, nc
      ! Initialize s and v to zero for each cluster
    !  s = 0.0_DP
    !  v = 0.0_DP
  
      ! Accumulate within-cluster sum of u into s
      s = SUM(u(start(k):start(k) + len(k) - 1, :), DIM=1)  ! Vectorized accumulation
  
      do i=1,p
        do j=1,p
          w(i,j) = w(i,j) + s(i) * s(j)
        end do
      end do
      
    END DO
  
    RETURN
  END SUBROUTINE robcovf
  