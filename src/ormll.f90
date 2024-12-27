subroutine ormll(n, k, p, x, y, offset, wt, penmat, link, alpha, beta, logL, u, &
                 ha, hb, hab, &
                 what, debug, penhess, salloc)

! n      : # observations
! k      : # intercepts
! p      : # x columns
! x      : covariate matrix
! y      : integer outcome vector with values 0 : k
! offset : n-vector of offsets
! wt     : n-vector of case weights
! penmat : p x p penalty matrix (ignored for Hessian in case x was QR-transformed)
! link   : link function/distribution famiily, 1-5 (see definition of cdf below)
! alpha  : k-vector of intercepts
! beta   : p-vector of regression coefficients
! logL   : -2 LL
! u      : returned k + p-vector of 1st derivatives
! ha     : returned k x 2 matrix for alpha portion of Hessian for sparse
!          representation with the Matrix R package; first column is the 
!          diagonal, second is the superdiagonal (last element ignored)
! hb     : returned p x p matrix for beta portion of Hessian
! hab    : returned k x p matrix for alpha x beta portion of Hessian
! what   : 1 for -2 LL, 2 for -2 LL and gradient, 3 for those plus Hessian
! debug  : 1 to print input parameters/sizes and quit, 2 to also print Hessian info
!          and not quit, 0 otherwise
! penhess: 1 to penalize Hessian with penmat, 0 to ignore penmat
! salloc : 0 if dynamic array allocation succeeded, > 0 if not
     
  use, intrinsic :: ISO_FORTRAN_ENV, only: dp => real64, int32
  implicit none
  integer(int32), intent(in)  :: n, y(n), k, p, what, debug, penhess, link
  real(dp),       intent(in)  :: x(n, p), offset(n), wt(n), penmat(p, p), & 
                                 alpha(k), beta(p)
  real(dp),       intent(out) :: logL, u(k + p), &
                                 ha(k, 2), hb(p, p), hab(k, p)
  integer(int32), intent(out) :: salloc
    
  integer(int32)  :: i, j, l, c, n0, nk, nb, ii, j1
  real(dp)                    :: w, z
  real(dp), allocatable :: lp(:), d(:), &
                           p1(:), p2(:), pdf1(:), pdf2(:), dpdf1(:), dpdf2(:)
  ! which obs have y=0, y=k, 0<y<k
  ! Suffix b = "between"
  integer(int32), allocatable :: i0(:), ik(:), ib(:), y1(:)
  
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
    call intpr('ha', 4, size(ha), 1)
    call intpr('debug', 5, debug, 1)
    call dblepr('alpha', 5, alpha, k)
    call dblepr('beta', 4, beta, p)
    call dblepr('penmat', 6, penmat, p * p)
  end if

  n0 = count(y == 0)
  nk = count(y == k)
  nb = count(y > 0 .and. y < k)

  allocate(lp(n), i0(n0), ik(nk), ib(nb), y1(n), &
           p1(n), p2(n), d(n), pdf1(n), pdf2(n), &
           dpdf1(n), dpdf2(n), stat=salloc)

  if(salloc /= 0) return

  i0 = pack([(i, i=1,n)], y == 0)             ! row numbers for which y==0
  ik = pack([(i, i=1,n)], y == k)
  ib = pack([(i, i=1,n)], y > 0 .and. y < k)

  lp = offset
  if(p > 0) lp = lp + matmul(x, beta)

  ! Model:
  ! Pr(Y = 0) = 1 - F(alpha(1) + lp) = F(-alpha(1) - lp) (logit, probit only)
  ! Pr(Y = k) =     F(alpha(k) + lp)
  ! Pr(Y = j) = F(alpha(j) + lp) - F(alpha(j+1) + lp), 0 < j < k

  ! p1(i0) = 1.0_dp - cdf(alpha(1)+ lp(i0), link)
  ! p1(ib) = cdf(alpha(y(ib))     + lp(ib), link)
  ! p1(ik) = cdf(alpha(k)         + lp(ik), link)
  ! Defer 1 - p1 for i0 until derivatives are computed
  y1     = max(y, 1_int32)
  p1     = cdf(alpha(y1)        + lp,     link)
  p2     = 0_dp
  p2(ib) = cdf(alpha(y(ib) + 1) + lp(ib), link)

  ! Compute all derivatives of cdf corresponding to p1 and p2
  ! Variables ending with 2 only apply to 0 < Y < k so are left at 0
  ! for Y=0 or k  

  ! pdf1(i0) =   pdf(alpha(1)       + lp(i0), 1.0_dp - p1(i0), link)   ! first term -
  ! pdf1(ib) =   pdf(alpha(y(ib))   + lp(ib), p1(ib), link)
  ! pdf1(ik) =   pdf(alpha(k)       + lp(ik), p1(ik), link)
  pdf1     = pdf(alpha(y1)       + lp,     p1,      link)  ! for vectors, max = R pmax
  pdf2     = 0_dp
  pdf2(ib) = pdf(alpha(y(ib) + 1) + lp(ib), p2(ib), link)

  ! Compute all second derivatives of cdf corresponding to p1 and p2
  dpdf2     = 0_dp
  
  ! dpdf1(i0) = dpdf(alpha(1)         + lp(i0), p1(i0), pdf1(i0), link)
  ! dpdf1(ib) = dpdf(alpha(y(ib))     + lp(ib), p1(ib), pdf1(ib), link)
  ! dpdf1(ik) = dpdf(alpha(k)         + lp(ik), p1(ik), pdf1(ik), link)
  dpdf1     = dpdf(alpha(y1)        + lp,     p1,     pdf1,     link)
  dpdf2(ib) = dpdf(alpha(y(ib) + 1) + lp(ib), p2(ib), pdf2(ib), link)

  p1(i0) = 1.0_dp - p1(i0)   ! Handle Y=0 case now that derivatives are computed
  d      = p1 - p2

  if(debug > 0) then
    call dblepr('p1', 2, p1, size(p1))
    call dblepr('p2', 2, p2, size(p2))
    call dblepr('pdf1', 4, pdf1, size(pdf1))
    call dblepr('pdf2', 4, pdf2, size(pdf2))
    call dblepr('dpdf1', 5, dpdf1, size(dpdf1))
    call dblepr('dpdf2', 5, dpdf2, size(dpdf2))
  end if

  logL = -2_dp * sum(wt * log(d)) +  dot_product(beta, matmul(penmat, beta))

  u   = 0_dp
  ha  = 0.0_dp
  hb  = 0.0_dp
  hab = 0.0_dp

  if(what == 1) then
    deallocate(lp, i0, ik, ib, y1, d, &
               p1, p2, pdf1, pdf2, dpdf1, dpdf2)
    return
  end if

  ! Y=0
  ! Probability: 1 - cdf(alpha(1) + lp) = p1
  !  D log p /D theta = - 1 / p1 * pdf
  ! Y=k
  ! Probability: cdf(alpha(k) + lp) = p1
  !  D log p / D theta = 1 / p1 * pdf
  ! 0 < Y < k
  ! D log(d) = (1/d) (p1' D() - p2' D())
  !  D p1 or p2 = pdf() D()
  !  => [pdf(alpha(y) + lp) D( ) + pdf(alpha(y+1) + lp) D( )] / d
  ! () = argument to pdf()

  ! Gradient (score vector)
  
  ! All obs with y=0
  u(1) = - sum(wt(i0) * pdf1(i0) / p1(i0))
  if(p > 0) then
    do l = 1, p
      u(k + l) = - sum(wt(i0) * pdf1(i0) * x(i0, l) / p1(i0))
    end do
  end if
  ! All obs with y=k
  u(k) = u(k) + sum(wt(ik) * pdf1(ik) / p1(ik))
  if(p > 0) then
    do l = 1, p
      u(k + l) = u(k + l) + sum(wt(ik) * pdf1(ik) * x(ik, l) / p1(ik))
    end do
  end if
  ! All obs with 0 < y < k
  if(nb > 0) then
    do ii = 1, nb
      i = ib(ii)   ! original row # of ii'th observation with 0 < y < k
      j = y(i)
      ! For p1, D() = 1 for alpha(j), 0 for alpha(j+1)
      ! For p2, D() = 0 for alpha(j), 1 for alpha(j+1)
      u(j)     = u(j)     + wt(i) * pdf1(i) / d(i)
      u(j + 1) = u(j + 1) - wt(i) * pdf2(i) / d(i)

      if(p > 0) then
        do l = 1, p
          u(k + l) = u(k + l) + wt(i) * x(i, l) * (pdf1(i) - pdf2(i)) / d(i)
        end do
      end if
    end do
  end if

  ! Add derivative of penalty function -0.5 b'Pb = -Pb
  if(p > 0) u((k + 1) : (k + p)) = u((k + 1) : (k + p)) - matmul(penmat, beta)

  if(what == 2) then
    deallocate(lp, i0, ik, ib, y1, d, &
               p1, p2, pdf1, pdf2, dpdf1, dpdf2)
    return
end if

  
  ! Hessian
  ! Create 3 matrices: ha (k x 2), hb (p x p), hab (k x p)
  ! ha represents a sparse matrix to be served up to the Matrix package in R
  ! It is tri-band diagonal, and since it is symmetric we only need to compute
  ! two bands (diagonal ha[, 1] and superdiagonal ha[, 2] with ha[k, 2] irrelevant.)
  ! For derivation of the Hessian components see
  ! https://fharrell.com/post/mle#sec-hessformula
  
  
  do i = 1, n
    j  = y(i)
    j1 = max(j, 1)
    w  = wt(i) * 1.0_dp / d(i) ** 2
    if(j == 0 .or. j == k) then
      if(j == 0) then
        z = - (d(i) * dpdf1(i) + pdf1(i) ** 2)
      else
        z =    d(i) * dpdf1(i) - pdf1(i) ** 2
      end if
      z = w * z
      ha(j1, 1) = ha(j1, 1) + z
      if(p > 0) then
        do l = 1, p
          hab(j1, l) = hab(j1, l) + x(i, l) * z
          do c = l, p
            hb(l, c) = hb(l, c) + x(i, l) * x(i, c) * z
          end do
        end do 
      end if
    else    ! 0 < Y < k
      ! D alpha(j)^2:
      ha(j, 1) = ha(j, 1) + w * (d(i) * dpdf1(i) - pdf1(i) ** 2)
      ! D alpha(j + 1)^2:
      ha(j + 1, 1) = ha(j + 1, 1) - w * &
          (d(i) * dpdf2(i) + pdf2(i) ** 2)
      ! alpha(j1), alpha(j1 + 1):
      ha(j, 2) = ha(j, 2) + w * pdf1(i) * pdf2(i)
      if(p > 0) then
        do l = 1, p
          ! D alpha(j)^2:
          hab(j, l) = hab(j, l) + w * x(i, l) * &
             (d(i) * dpdf1(i) - pdf1(i) * (pdf1(i) - pdf2(i)))
          ! D alpha(j+1)^2:
          hab(j + 1, l) = hab(j + 1, l) - w * x(i, l) * &
              (d(i) * dpdf2(i) - pdf2(i) * (pdf1(i) - pdf2(i)))   ! was + pdf2(i) by ChatGPT
          do c = l, p
            hb(l, c) = hb(l, c) + w * x(i, l) * x(i, c) * &
              (d(i) * (dpdf1(i) - dpdf2(i)) - (pdf1(i) - pdf2(i)) ** 2)
          end do
        end do
      end if
    end if
  end do

  ! Finish symmetric matrix
  if(p > 0) then
    do l = 1, p - 1
    do c = l + 1, p
      hb(c, l) = hb(l, c)
    end do
  end do
  end if

  if(debug > 0) call intpr('hess A', 6, 0, 1)
  ! To add derivative of penalty function -0.5 b'Pb = -Pb :
  if(p > 0 .and. penhess > 0) hb = hb - penmat

if(debug > 0) then
  call dblepr('ha',  2, ha, size(ha))
  call dblepr('hb',  2, hb, size(hb))
  call dblepr('hab', 3, hab, size(hab))
end if

deallocate(lp, i0, ik, ib, y1, d, &
           p1, p2, pdf1, pdf2, dpdf1, dpdf2)

return

contains

  ! Compute CDF per link function given x
  real(dp) function cdf(x, link) result(p)
  real(dp),       intent(in) :: x(:)
  integer(int32), intent(in) :: link
  allocatable :: p(:)

  select case(link)
  case(1)              ! logistic
    p = 1.0_dp / (1.0_dp + exp(- x))
  case(2)              ! probit
    p = 0.5_dp * (1.0_dp + erf(x / sqrt(2.0_dp)))
  case(3)
    p = exp(-exp(-x))  ! loglog
  case(4)              ! complementary loglog
    p = 1 - exp(-exp(x))
  case(5)              ! Cauchy
    p = (1.0_dp / 3.14159265358979323846_dp) * atan(x) + 0.5_dp
  end select

  end function cdf

  ! Compute probability density function (derivative of cdf) given x and cdf f
  ! cdf is used as extra input to save time
  ! Note: f is the value returned from cdf() in pure form.  Likewise for deriv
  ! in dpdf.  For example if you computed 1 - cdf( ), f=cdf( ) not 1 - cdf( ).
  real(dp) function pdf(x, f, link) result(p)
    real(dp),       intent(in) :: x(:), f(:)
    integer(int32), intent(in) :: link
    allocatable                :: p(:)

    select case(link)
    case(1)
      p = f * (1_dp - f)
    case(2)
      p = (1_dp / sqrt(2_dp * 3.14159265358979323846_dp)) * exp(- x * x / 2.0_dp)
    case(3)
      p = exp(-x - exp(-x))
    case(4)
      p = exp(x - exp(x))
    case(5)
      p = (1.0_dp / 3.14159265358979323846_dp) / (1_dp + x * x)
    end select

    end function pdf

  ! Compute 2nd derivative of cdf (derivative of pdf) given x, cdf, pdf
  real(dp) function dpdf(x, f, deriv, link) result(p)
  real(dp),       intent(in) :: x(:), f(:), deriv(:)
  integer(int32), intent(in) :: link
  allocatable                :: p(:)

  select case(link)
  case(1)
    p = f * (1_dp - 3_dp * f + 2 * f * f)
  case(2)
    p = - deriv * x
  case(3)
    p = deriv * (-1_dp + exp(-x))
  case(4)
    p = deriv * (1_dp - exp(x))
  case(5)
    p = -2_dp * x * ((1 + x * x) ** (-2_dp)) / 3.14159265358979323846_dp
  end select

end function dpdf

end subroutine ormll
