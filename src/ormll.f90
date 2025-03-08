subroutine ormll(n, k, p,  &
                 x, y, y2, offset, wt, penmat, link, alpha, beta, logL, u, &
                 d, ha, hb, hab,                 &
                 intcens, row, col, ai, nai, ne, &
                 urow, ucol, um, nu, nuu,        &
                 what, debug, penhess, salloc)

! n      : # observations
! k      : # intercepts
! p      : # x columns
! ne     : # elements set aside for sparse intercept
!          hessian when interval censoring is present
! nu     : # elements set aside for sparse score matrix
!           Set to zero to not compute the score matrix
!           To be safe, set nu = n * (2 + p)
! intcens: 1 if any interval censoring, 0 otherwose
! x      : covariate matrix
! y      : integer outcome vector with values 0 : k
! y2     : for censoring, y=y2 => uncensored; y2 is 0 : k
!           y  =  -1 for left censored observation
!           y2 = k+1 for right censored observation
!           0 <= y < y1 <= k for interval-censored data
! offset : n-vector of offsets
! wt     : n-vector of case weights
! penmat : p x p penalty matrix (ignored for Hessian in case x was QR-transformed)
! link   : link function/distribution famiily, 1-5 (see definition of cdf below)
! alpha  : k-vector of intercepts
! beta   : p-vector of regression coefficients
! logL   : -2 LL
! u      : returned k + p-vector of 1st derivatives
! d      : returned per-observation probability element used in log likelihood calculation
! ha     : returned k x 2 matrix for alpha portion of Hessian for sparse
!          representation with the Matrix R package; first column is the 
!          diagonal, second is the superdiagonal (last element ignored)
! hb     : returned p x p matrix for beta portion of Hessian
! hab    : returned k x p matrix for alpha x beta portion of Hessian
! nai    : number of elements set aside for row, col, ai when intcens=1
! row    : returned vector of row numbers for intercept hessian with interval censoring
! col    : returned column numbers
! ai     : returned intercept hessian entry when intcens=1
! ne     : returned number of row, col, ai actually used
! urow   : returned integer vector of row numbers for score matrix
! ucol   : returned integer vector of column numbers for score matrix
! um     : returned real vector of score matrix elements going with urow, ucol
! nuu    : number of score matrix elements actually computed
! what   : 1 for -2 LL, 2 for -2 LL and gradient, 3 for those plus Hessian
! debug  : 1 to print input parameters/sizes, 2 to also print Hessian info
!          0 otherwise
! penhess: 1 to penalize Hessian with penmat, 0 to ignore penmat
! salloc : 0 if dynamic array allocation succeeded, > 0 if not,
!          999 if negative or zero Y=j probability encountered
!          998 if censored data configuration encountered that is not implemented
!          997 if hessian needs more than 1000000 elements due to the variety
!              of interval-censored values
!          996 if nu > 0 but is not large enough to hold all score vector elements
     
  use, intrinsic :: ISO_FORTRAN_ENV, only: dp => real64, int32
  implicit none
  integer(int32), intent(in)  :: n, y(n), y2(n), k, p, what, debug, penhess, link, &
                                 nai, intcens, nu
  real(dp),       intent(in)  :: x(n, p), offset(n), wt(n), penmat(p, p), & 
                                 alpha(k), beta(p)
  real(dp),       intent(out) :: logL, d(n), u(k + p), &
                                 ha(k * (1 - intcens), 2), hb(p, p), hab(k, p), &
                                 ai(nai), um(nu)
  integer(int32), intent(out) :: row(nai), col(nai), ne, urow(nu), ucol(nu), nuu, salloc
    
  integer(int32)  :: i, j, l, c, nb, j2, a, b, nbad, il(1)
  real(dp)                    :: w, z, g1, g2
  real(dp), allocatable :: lp(:), sgn(:), &
                           p1(:), p2(:), pdf1(:), pdf2(:), dpdf1(:), dpdf2(:)
  ! ib: which obs have y=0, y=k, 0<y<k or are interval censored
  ! Suffix b = "between", suffix a = "alpha"
  integer(int32), allocatable :: ib(:), ia(:), ia2(:), ibad(:)
  
  if(debug > 0) then
    call intpr('n,k,p,ic,nai,x,y,o,w,pen,a,b,ha,d,r,c,ai', 40, &
      [n, k, p, intcens, nai, size(x), size(y), size(offset), size(wt), size(penmat), &
       size(alpha), size(beta), size(ha), size(d), size(row), size(col), size(ai)],   17)
    call dblepr('alpha',  5, alpha,  k)
    call dblepr('beta',   4, beta,   p)
    call dblepr('penmat', 6, penmat, p * p)
  end if

  nb = count((y > 0 .and. y < k .and. y == y2) .or. & ! uncensored Y, 0 < Y < k
             (y >= 0 .and. y2 <= k .and. y /= y2))    ! or interval censored

  allocate(lp(n), ib(nb), ia(n), ia2(n), sgn(n), &
           p1(n), p2(n), pdf1(n), pdf2(n), &
           dpdf1(n), dpdf2(n), stat=salloc)

  if(salloc /= 0) return

  ! Compute observation numbers of uncensored data with 0 <= y <= k or
  ! interval censored observations involving two alphas
  ! Interval censored [a, k] just involves one alpha
  ib = pack([(i, i=1,n)], (y > 0 .and. y < k .and. y == y2) .or. &
                          (y >= 0 .and. y2 < k .and. y /= y2)      )

  lp = offset
  if(p > 0) lp = lp + matmul(x, beta)

  ! Model:
  ! Pr(Y = 0) = 1 - F(alpha(1) + lp)
  ! Pr(Y = k) = F(alpha(k) + lp)
  ! Pr(Y = j) = F(alpha(j) + lp) - F(alpha(j+1) + lp), 0 < j < k, uncensored
  ! Pr(Y < j) = 1 - F(alpha(j) + lp)     ! left censored at j
  ! Pr(Y > j) = F(alpha(j+1) + lp)       ! right censored at j, j < k
  ! Pr(Y > k) = F(alpha(k) + lp)         ! right censored at k
  !                                      ! => reinterpret F(alpha(k) + lp) as P(Y >= k) ??
  ! The first F() corresponds to p1.  For 0 < y < k p2 corresponds to second F()
  ! General formula:
  !  ia  = index of involved alpha for p1 term
  !  s   = -1 for Y=0 or left censoring, +1 otherwise
  !  [s = -1] + s * F(alpha(ia) + lp) - [0 < Y < k, uncensored] * F(alpha(ia + 1) + lp)
  !  [s = -1] + s * p1 - [0 < Y < k, uncensored] * p2
  !
  ! For interval censored observations [j, l], 0 <= j,l <= k
  !   Pr(j <= Y <= l) = F(alpha(j) + lp) - F(alpha(l + 1) + lp), j > 0, l < k
  !   Pr(0 <= Y <= l) = 1 - Pr(Y > l) = 1 - F(alpha(l + 1) + lp), l < k (s = -1)
  !   Pr(j <= Y <= k) = Pr(Y >= j) = F(alpha(j) + lp), j > 0, j < k

  ! Compute ia : which alpha is involved (first alpha if 0 < y < k & uncensored)
  !         ia2: second alpha involved, will always have negated F(); ia2=0 if no second alpha
  !         ia goes with p1, ia2 goes with p2, take p1 - p2
  !         s  : 1.0 when prob is F(), -1.0 when prob is 1 - F() 

  sgn = 1_dp
  ia2 = 0_int32
  do i = 1, n
    a  = y (i)
    b  = y2(i)
    if(a == b) then       ! uncensored
      if(a == 0) then
        ia(i)  = 1        ! alpha(1)
        sgn(i) = -1_dp
      else if(a == k) then
        ia(i) = k
      else                ! 0 < a < k
        ia(i)  = a        ! alpha number of p1; for p2 is ia + 1
        ia2(i) = a + 1
      end if
    else                       ! a not= b: censored
      if(a == -1_int32) then   ! left censored
        ia(i)  = max(b, 1)
        sgn(i) = -1_dp
      else if(b > k) then ! right censored
        ! It is possible that the highest right-censored a-value is at a=k
        ! In that case the intercept involved is a=k and the interpretation
        ! of the fitted model for the highest value of y (y=k) is
        ! P(Y >= k | X) instead of P(Y == k | X)
        ! If right censoring occurs when b=k the observation is treated the
        ! same as an uncensored point in the likelihood calculations
        ia(i) = min(a + 1, k)
      else if(a == 0 .and. b < k) then      ! interval censored [0, b]
        ia(i)  = b + 1
        sgn(i) = -1_dp
      else if(a > 0 .and. b == k) then  ! interval censored [a, k]
        ia(i) = a
      else if(a > 0 .and. b < k .and. a < b) then
        ia(i)  = a
        ia2(i) = b + 1
      else
        salloc = 998
        return
      end if
    end if
  end do

  ! Compute first probability component without applying sgn, as this will become an
  ! argument for derivative functions to reduce execution time
  p1  = cdf(alpha(ia) + lp, link)
  ! Compute second probability component for in-between y observations
  p2  = 0_dp
  if(nb > 0) p2(ib) = cdf(alpha(ia2(ib)) + lp(ib), link)
  ! Compute probability element for likelihood
  d = merge(p1 - p2, 1_dp - p1, sgn == 1_dp)

  if(debug > 0) then
    call intpr('ia',           2, ia,          size(ia))
    call intpr('ia2',          3, ia2,         size(ia2))
    call dblepr('alpha(ia)',   9, alpha(ia),   size(ia))
    call dblepr('alpha(ia2)', 10, alpha(ia2),  size(ia2))
    call dblepr('sgn',         3, sgn,         size(sgn))
    call dblepr('lp',          2, lp,          size(lp))
  end if

  nbad = count(d <= 0_dp)
  if(nbad > 0_int32) then
    if(debug > 0) then
      allocate(ibad(nbad))
      ibad = pack([(i, i=1,n)], d <= 0_dp)
      call intpr('Zero or negative probability for observations ', 45, ibad, nbad)
      call intpr('Intercept involved', 18, ia(ibad), nbad)
      if(any(ia2(ibad) > 0)) call intpr('2nd Intercept involved', 22, ia2(ibad), nbad)
      call intpr('y',    1, y  (ibad),   nbad)
      call intpr('y2',   2, y2 (ibad),   nbad)
      call dblepr('d',   1, d  (ibad),   nbad)
      call dblepr('p1',  2, p1 (ibad),   nbad)
      call dblepr('p2',  2, p2 (ibad),   nbad)
      call dblepr('sgn', 3, sgn(ibad),   nbad)
      deallocate(ibad)
    end if
    salloc = 999_int32
    deallocate(lp, ib, ia, ia2, sgn, p1, p2, pdf1, pdf2, dpdf1, dpdf2)
    return
  end if

  if(debug > 0) then
    call dblepr('alpha', 5, alpha, k)
    call dblepr('beta',  4, beta,  p)
    call dblepr('sgn',   3, sgn,   n)
    call dblepr('p1',    2, p1,    size(p1))
    call dblepr('p2',    2, p2,    size(p2))
    call dblepr('d',     1, d,     size(d))
  end if

  logL = -2_dp * sum(wt * log(d)) +  dot_product(beta, matmul(penmat, beta))

  u   = 0_dp
  ha  = 0_dp
  hb  = 0_dp
  hab = 0_dp

  if(what == 1) then
    deallocate(lp, ib, ia, ia2, sgn, &
               p1, p2, pdf1, pdf2, dpdf1, dpdf2)
    return
  end if

  ! Probability: [s = 1] + s * cdf(alpha(ia) + lp) - p2 = d
  !            = [s = 1] + s * p1 - p2 (p2 = 0 for L/R censored y or y=0,k)
  ! D log d / D theta = s * pdf(alpha(ia) + lp) D() / d - pdf(alpha(ia2) + lp) D() / d
  ! For L/R censored y or y=0, k the second term is ignored
  ! D() = D(argument to pdf) / D theta
  ! D log d / D theta = [s * pdf1 - pdf2] / d

  ! Gradient (score vector)

  pdf1     = pdf(alpha(ia)      + lp,     p1,     link)
  pdf2     = 0_dp
  pdf2(ib) = pdf(alpha(ia2(ib)) + lp(ib), p2(ib), link)
  nuu      = 0_int32

  if(debug > 0) then
    call dblepr('pdf1',  4, pdf1,  size(pdf1))
    call dblepr('pdf2',  4, pdf2,  size(pdf2))
  end if

  if(nu > 1) then
    urow = 0_int32
    ucol = 0_int32
    um   = 0_dp
  end if

  do i = 1, n
    j    = ia(i)    ! subscript of applicable alpha
    j2   = ia2(i)   ! subscript of second alpha, 0 if not there
    w    = wt(i) / d(i)
    g1   = pdf1(i) * sgn(i)
    g2   = pdf2(i)
    u(j) = u(j) + w * g1
    if(j2 > 0) u(j2) = u(j2) - w * g2
    if(p > 0) then
      do l = 1, p
        u(k + l) = u(k + l) + w * (g1 - g2) * x(i, l)
      end do
    end if

    if(nu > 0) then  ! compute sparse score matrix elements
      if(nuu + 2_int32 + p > nu) then
        salloc = 996_int32
        deallocate(lp, ib, ia, ia2, sgn, p1, p2, pdf1, pdf2, dpdf1, dpdf2)
        return
      end if
      nuu = nuu + 1
      urow(nuu) = i
      ucol(nuu) = j
      um(nuu)   = w * g1
      if(j2 > 0) then
        nuu = nuu + 1
        urow(nuu) = i
        ucol(nuu) = j2
        um(nuu)   = - w * g2
      end if
      if(p > 0) then
        do l = 1, p
          nuu = nuu + 1
          urow(nuu) = i
          ucol(nuu) = k + l
          um(nuu)   = w * (g1 - g2) * x(i, l)
        end do
      end if
    end if

  end do

  ! Add derivative of penalty function -0.5 b'Pb = -Pb
  ! Ignored at present for mscore
  if(p > 0) u((k + 1) : (k + p)) = u((k + 1) : (k + p)) - matmul(penmat, beta)

  if(debug > 0) call dblepr('u', 1, u, k + p)

  if(what == 2) then
    deallocate(lp, ib, ia, sgn, &
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
  ! If there are any interval-censored observatons, ha is instead computed
  ! as ai using general sparse form
  
  ! Compute all second derivatives of cdf corresponding to p1 and p2
  
  ! dpdf1(i0) = dpdf(alpha(1)         + lp(i0), p1(i0), pdf1(i0), link)
  ! dpdf1(ib) = dpdf(alpha(y(ib))     + lp(ib), p1(ib), pdf1(ib), link)
  ! dpdf1(ik) = dpdf(alpha(k)         + lp(ik), p1(ik), pdf1(ik), link)

  ! For the sometimes-interval-censoring form of ha,
  ! initialize row and col for intercept hessian elements that are always used
  ! E.g. if k=4 we use (row,col) (1,1), (2,2), (3,3), (4,4), (l,2), (2,3), (3,4)

  ne = 0
  if(intcens == 1) then
    row(1 : k) = [(i, i=1, k)]    ! diagonal elements
    col(1 : k) = row(1 : k)
    row((k + 1) : (2 * k - 1)) = [(i, i=1, k - 1)]   ! minor diagonal above major one
    col((k + 1) : (2 * k - 1)) = [(i, i=2, k)]
    ne = 2_int32 * k - 1_int32
    ai(1 : ne)  = 0_dp
    end if

  dpdf1     = dpdf(alpha(ia)        + lp,     p1,     pdf1,     link)
  dpdf2     = 0_dp
  dpdf2(ib) = dpdf(alpha(ia2(ib)) + lp(ib), p2(ib), pdf2(ib), link)

  if(debug > 0) then
    call dblepr('dpdf1', 5, dpdf1, size(dpdf1))
    call dblepr('dpdf2', 5, dpdf2, size(dpdf2))
  end if

  
  do i = 1, n
    a  = y(i)
    j  = ia(i)     ! intercept involved in p1
    j2 = ia2(i)    ! intercept involved in p2 (0 if not there)
    w  = wt(i) * 1_dp / d(i) ** 2
    z  = w * (sgn(i) * d(i) * dpdf1(i) - pdf1(i) ** 2)
    ! Intercept-only part of hessian, starting with the no-interval-censored case
    if(intcens == 0) then
      ha(j, 1) = ha(j, 1) + z
      if(j2 /= 0) then
        ha(j2, 1) = ha(j2, 1) - w * (d(i) * dpdf2(i) + pdf2(i) ** 2)  ! diagonal
        ha(j,  2) = ha(j, 2)  + w * pdf1(i) * pdf2(i)                 ! super diagonal
      end if
    else                 ! some interval censored observations
      ai(j) = ai(j) + z  ! diagonal element; place in ai already allocated
      if(j2 > 0) then    ! second intercept involved
        ai(j2) = ai(j2) - w * (d(i) * dpdf2(i) + pdf2(i) ** 2)    ! diagonal element
        if(j2 == (j + 1)) then   ! adjacent intercepts involved and place already allocated
          ai(k + j2 - 1) = ai(k + j2 - 1) + w * pdf1(i) * pdf2(i) ! super diagonal
        else
          ! involves 2 non-adjacent intercepts; may have to allocate new position in ai
          if(any((row(1 : ne) == j) .and. (col(1 : ne) == j2))) then
            il = findloc(row(1 : ne), j, mask = (col(1 : ne) == j2))
            l  = il(1)
            ai(l) = ai(l) + w * pdf1(i) * pdf2(i)
          else   ! add a new entry
            ne = ne + 1
            if(ne > nai) then
              salloc = 997
              return
            end if
            row(ne) = j
            col(ne) = j2
            ai(ne)  = w * pdf1(i) * pdf2(i)
          end if
        end if
      end if
    end if

    if(p == 0) cycle

    if(j2 == 0) then      ! only one intercept involved
      do l = 1, p
        hab(j, l) = hab(j, l) + x(i, l) * z
        do c = l, p
          hb(l, c) = hb(l, c) + x(i, l) * x(i, c) * z
        end do
      end do 
    else                    ! two intercepts involved
      do l = 1, p
        ! D alpha(j)^2:
        hab(j, l) = hab(j, l) + w * x(i, l) * &
               (d(i) * dpdf1(i) - pdf1(i) * (pdf1(i) - pdf2(i)))
        ! D alpha(j+1)^2:
        hab(j2, l) = hab(j2, l) - w * x(i, l) * &
              (d(i) * dpdf2(i) - pdf2(i) * (pdf1(i) - pdf2(i)))
        do c = l, p
          hb(l, c) = hb(l, c) + w * x(i, l) * x(i, c) * &
              (d(i) * (dpdf1(i) - dpdf2(i)) - (pdf1(i) - pdf2(i)) ** 2)
        end do
      end do
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

  if(debug > 1) call intpr1('hess A', 6, 0)

  ! To add derivative of penalty function -0.5 b'Pb = -Pb :
  if(p > 0 .and. penhess > 0) hb = hb - penmat

  if(debug > 1) then
    call dblepr('ha',  2, ha,  size(ha))
    call dblepr('hb',  2, hb,  size(hb))
    call dblepr('hab', 3, hab, size(hab))
  end if

deallocate(lp, ib, ia, ia2, sgn, &
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
    p = 0.5_dp * (1.0_dp + erf(x / 1.414213562373095_dp))
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
