subroutine lrmll(n, k, p, x, y, offset, wt, penmat, alpha, beta, logL, u, &
                 ha, hb, hab, what, debug, penhess, salloc)
! n      : # observations
! k      : # intercepts
! p      : # x columns
! x      : covariate matrix
! y      : integer outcome vector with values 0 : k
! offset : n-vector of offsets
! wt     : n-vector of case weights
! penmat : p x p penalty matrix (ignored for Hessian in case x was QR-transformed)
! alpha  : k-vector of intercepts
! beta   : p-vector of regression coefficients
! logL   : returned -2 LL if what=1 or 3
! u      : returned k + p-vector of 1st derivatives if what=2
! ha     : returned k x 2 matrix for alpha portion of Hessian for sparse
!          representation with the Matrix R package; first column is the 
!          diagonal, second is the superdiagonal (last element ignored)
! hb     : returned p x p matrix for beta portion of Hessian
! hab    : returned k x p matrix for alpha x beta portion of Hessian
!!! hess   : returned k + p x k + p matrix of 2nd derivatives if what=3
!          should be zero length if what is not 3
! what   : 1 for -2LL only, 2 for gradient, 3 for hessian
! debug  : 1 to print input parameters/sizes and quit, 2 to also print Hessian info
!          and not quit, 0 otherwise
! penhess: 1 to penalize Hessian with penmat, 0 to ignore penmat for what=3
     
  use, intrinsic :: ISO_FORTRAN_ENV, only: dp => real64, int32
  implicit none
  integer(int32), intent(in)  :: n, y(n), k, p, what, debug, penhess
  real(dp),       intent(in)  :: x(n, p), offset(n), wt(n), penmat(p, p), alpha(k), beta(p)
  real(dp),       intent(out) :: logL, u(k + p), &
                                 ha(k, 2), hb(p, p), hab(k, p)
  integer(int32), intent(out) :: salloc
  integer(int32)  :: i, j, l, c, n0, nk, nb, ii, j1
  real(dp), allocatable :: lp(:), ww(:), d(:), &
                           p1(:), p2(:), v1(:), v2(:), w1(:), w2(:)
  ! which obs have y=0, y=k, 0<y<k
  ! Suffix b = "between"
  integer(int32), allocatable :: i0(:), ik(:), ib(:)

  if(debug > 0) then
    call intpr('n', 1, n, 1)
    call intpr('k', 1, k, 1)
    call intpr('p', 1, p, 1)
    call intpr('x', 1, size(x), 1)
    call intpr('y', 1, size(y), 1)
    call intpr('offset', 6, size(offset), 1)
    call intpr('wt', 2, size(wt), 1)
    call intpr('penmat', 6, size(penmat), 1)
    call intpr('alpha', 5, size(alpha), 1)
    call intpr('beta', 4, size(beta), 1)
    call intpr('ha', 2, size(ha), 1)
    call intpr('what', 4, what, 1)
    call intpr('debug', 5, debug, 1)
    call dblepr('alpha', 5, alpha, k)
    call dblepr('beta', 4, beta, p)
    call dblepr('penmat', 6, penmat, p * p)
  end if
 
  n0 = count(y == 0)
  nk = count(y == k)
  nb = count(y > 0 .and. y < k)

  allocate(lp(n), ww(p), i0(n0), ik(nk), ib(nb), &
           p1(n), p2(n), d(n), v1(n), v2(n), w1(n), w2(n), stat=salloc)

  i0 = pack([(i, i=1,n)], y == 0)             ! row numbers for which y==0
  ik = pack([(i, i=1,n)], y == k)
  ib = pack([(i, i=1,n)], y > 0 .and. y < k)

  lp = offset
  if(p > 0) lp = lp + matmul(x, beta)

  ! Model:
  ! Pr(Y = 0) = 1 - expit(alpha(1) + lp) = expit(-alpha(1) - lp)
  ! Pr(Y = k) =     expit(alpha(k) + lp)
  ! Pr(Y = j) = expit(alpha(j) + lp) - expit(alpha(j+1) + lp), 0 < j < k
  
  p1     = 0_dp
  p2     = 0_dp
  p1(i0) = expit(- alpha(1)       - lp(i0))   ! 1 - expit(alpha(1) + lp)
  p1(ib) = expit(alpha(y(ib))     + lp(ib))
  p2(ib) = expit(alpha(y(ib) + 1) + lp(ib))
  p1(ik) = expit(alpha(k)         + lp(ik))
  d      = p1 - p2

  logL = -2_dp * sum(wt * log(d)) +  dot_product(beta, matmul(penmat, beta))

  ! Derivative of log expit(x) is expit(-x)
  ! First derivatives of log(d) = log(p1) for Y=0:
  ! d = expit(-alpha(1) - lp)
  ! alpha(1): - expit(alpha(1)+lp) = -(1 - d)
  ! beta    : - expit(alpha(1)+lp) x = -(1 - d) x
  ! For Y=k
  ! d = expit(alpha(k) + lp)
  ! alpha(k): expit(-alpha(k)-lp) = 1 - d
  ! beta    : exp(-alpha(k)-lp) x = (1 - d) x
  ! For 0 < Y < k
  !  D log(d) = (1/d) (p1' D() - p2' D())
  !  D p1 or p2 = D expit(x) = p' = p(1-p) D()
  !  => [p1(1-p1) D() - p2(1-p2) D()] / d
  ! () = argument to expit()

  if(what /= 1_int32) then
    v1  = p1 * (1_dp - p1)
    v2  = p2 * (1_dp - p2)
  end if

  ! Gradient (per-parameter score vector)
  if(what == 2_int32) then
    u = 0_dp

    ! All obs with y=0
    ! The derivative of log expit(x) wrt x is expit(-x)
    ! Prob element is expit(-alpha(1) - lp)
    u(1) = - sum(wt(i0) * (1_dp - d(i0)))
    if(p > 0) then
      do l = 1, p
        u(k + l) = - sum(wt(i0) * x(i0, l) * (1_dp - d(i0)))
      end do
    end if
    ! All obs with y=k
    ! Prob element is expit(alpha(k) + lp)
    u(k) = u(k) + sum(wt(ik) * (1_dp - d(ik)))
    if(p > 0) then
      do l = 1, p
        u(k + l) = u(k + l) + sum(wt(ik) * x(ik, l) * (1_dp - d(ik)))
      end do
    end if
    ! All obs with 0 < y < k
    if(nb > 0) then
      do ii = 1, nb
        i = ib(ii)
        j = y(i)
        ! For p1, D() = 1 for alpha(j), 0 for alpha(j+1)
        ! For p2, D() = 0 for alpha(j), 1 for alpha(j+1)
        u(j)     = u(j)     + wt(i) * v1(i) / d(i)
        u(j + 1) = u(j + 1) - wt(i) * v2(i) / d(i)
        if(p > 0) then
          do l = 1, p
            u(k + l) = u(k + l) + wt(i) * x(i, l) * (v1(i) - v2(i)) / d(i)
          end do
        end if
      end do
    end if

    ! Add derivative of penalty function -0.5 b'Pb = -Pb
    if(p > 0) u((k + 1) : (k + p)) = u((k + 1) : (k + p)) - matmul(penmat, beta)
  end if


  ! Hessian
  ! For large k, hess is not storage-efficient because it stores all the zeros.
  ! It is computationally efficient because no terms are computed that are zero.
  ! I.e., computations respect the tri-band diagonal form of hess.

  if(what == 3_int32) then
    ha  = 0.0_dp
    hb  = 0.0_dp
    hab = 0.0_dp
    w1  = v1 * (1_dp - 2_dp * p1)
    w2  = v2 * (1_dp - 2_dp * p2)

    ! Second derivative of log d is (f''(d) x d - f'(d) x f'(d)) / d^2
    ! f(x) = Pr(Y = j); f'(d) is given above

    do i = 1, n
      j  = y(i)
      j1 = max(j, 1)
      if(j == 0 .or. j == k) then
        ha(j1, 1) = ha(j1, 1) - wt(i) * v1(i)
        if(p > 0) then
          do l = 1, p
            hab(j1, l) = hab(j1, l) - wt(i) * v1(i) * x(i, l)
            do c = l, p
            hb(l, c) = hb(l, c) - wt(i) * x(i, l) * x(i, c) * v1(i)
            end do
          end do
        end if
      else   ! 0 < y(i) < k
        ha(j1,     1) = ha(j1,     1) + wt(i) * (w1(i) * d(i) - v1(i) ** 2) / d(i) **2
        ha(j1 + 1, 1) = ha(j1 + 1, 1) + wt(i) * (-w2(i) * d(i) - v2(i) ** 2) / d(i) ** 2
        ha(j1,     2) = ha(j1,     2) + wt(i) * v1(i) * v2(i) / d(i) ** 2
        if(p > 0) then
          do l = 1, p
            hab(j1,     l) = hab(j1,     l) + wt(i) * x(i, l) * (  w1(i) * d(i) - &
                                                          v1(i) * (v1(i) - v2(i))) / d(i) ** 2
            hab(j1 + 1, l) = hab(j1 + 1, l) + wt(i) * x(i, l) * &
                                                          (- w2(i) * d(i) + v2(i) * (v1(i) - v2(i))) / d(i) ** 2
            do c = l, p
              hb(l, c) = hb(l, c) + &
                wt(i) * x(i, l) * x(i, c) * ((w1(i) - w2(i)) * d(i) - (v1(i) - v2(i)) ** 2) / d(i) ** 2
            end do
          end do
        end if
      end if
    end do

    ! Finish symmetric matrix
    do l=1, p - 1
      do c = l + 1, p
        hb(c, l) = hb(l, c)
      end do
    end do
    if(debug > 0) call intpr('hess A', 6, 0, 1)
    ! To add derivative of penalty function -0.5 b'Pb = -Pb :
    if(p > 0 .and. penhess > 0) hb = hb - penmat
  end if

if(debug > 0) call intpr('hab B', 5, size(hab), 1)

deallocate(lp, ww, i0, ik, ib, d, &
           p1, p2, v1, v2, w1, w2)

return

contains

  real(dp) function expit(x) result(r)
  real(dp), intent(in) :: x(:)
  allocatable :: r(:)
  ! allocate(r(size(x)))
  r = 1.0_dp / (1.0_dp + exp(-x))
  end function expit

  ! Indicator function: true/false becomes 1d0, 0d0
  real(dp) function ld(x)
    logical, intent(in) :: x
    if(x) then
      ld = 1.0_dp
    else
      ld = 0.0_dp
    end if
    end function ld

    
end subroutine lrmll
