subroutine ormll(n, k, p,  &
                 x, ia, ia2, sgn, ib, nb, offset, wt, penmat, link, alpha, beta, logL, u, &
                 d, ha, hb, hab,                 &
                 intcens, row, col, ai, nai, ne, &
                 urow, ucol, um, nu, nuu,        &
                 what, debug, penhess, salloc)

! Claude Sonnet 5 2026-08-30
! Modified to accept ia, ia2, sgn, ib, nb as precomputed inputs (produced
! once per fit by ormidx) instead of recomputing them internally on every
! call -- this bookkeeping is a pure function of y, y2, k and does not
! depend on alpha, beta, or offset, so recomputing it every
! Newton-Raphson iteration was wasted work. cdf/pdf/dpdf are now
! supplied by the orm_links module instead of an internal contains
! block, so the same link-function code can be shared with ormeta.
! All downstream math (from the p1/p2 construction onward) is
! byte-for-byte unchanged from the original.
!
! Claude Sonnet 5 2026-09-19
! All link-specific arithmetic (d, pdf, dpdf and the score/Hessian
! coefficients built from them) now comes from orm_links::ormatoms.
! This routine only accumulates.  ib and nb are no longer used (ormatoms
! selects two-intercept observations from ia2 > 0) but stay in the
! argument list so the .Fortran() call and init.c are unchanged.

! n      : # observations
! k      : # intercepts
! p      : # x columns
! ia     : intercept index for the "p1" probability component (from ormidx)
! ia2    : intercept index for "p2" component, 0 if none (from ormidx)
! sgn    : +1 if p1 enters as F(), -1 if as 1-F() (from ormidx)
! ib     : observation numbers having a p2 component; only ib(1:nb) used (from ormidx)
! nb     : count of such observations (from ormidx)
! ne     : # elements set aside for sparse intercept
!          hessian when interval censoring is present
! nu     : # elements set aside for sparse score matrix
!           Set to zero to not compute the score matrix
!           To be safe, set nu = n * (2 + p)
! intcens: 1 if any interval censoring, 0 otherwose
! x      : covariate matrix
! (y, y2 removed -- see note below at the ia/ia2/sgn declaration)
! offset : n-vector of offsets
! wt     : n-vector of case weights
! penmat : p x p penalty matrix (ignored for Hessian in case x was QR-transformed)
! link   : link function/distribution famiily, 1-5 (see orm_links module)
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
  use orm_links, only: ormatoms, wlogsum
  implicit none
  ! Claude Sonnet 5 2026-08-30          4 lines
  ! y, y2 removed entirely: every substantive use in this subroutine
  ! (the likelihood/gradient/Hessian computation) was already replaced
  ! by ia, ia2, sgn, ib, nb (from ormidx) in the prior revision -- the
  ! only remaining uses were debug-only diagnostics, now adapted below
  ! to print ia/ia2/sgn instead (a deterministic, information-equivalent
  ! transform of y/y2). A future y-dependent-effects (YDE / partial PO)
  ! extension should route new information through ia/ia2 (already
  ! cutpoint-indexed), not y/y2.
  integer(int32), intent(in)  :: n, k, p, what, debug, penhess, link, &
                                 nai, intcens, nu, nb, ia(n), ia2(n), ib(nb)
  real(dp),       intent(in)  :: x(n, p), sgn(n), offset(n), wt(n), penmat(p, p), &
                                 alpha(k), beta(p)
  real(dp),       intent(out) :: logL, d(n), u(k + p), &
                                 ha(k * (1 - intcens), 2), hb(p, p), hab(k, p), &
                                 ai(nai), um(nu)
  integer(int32), intent(out) :: row(nai), col(nai), ne, urow(nu), ucol(nu), nuu, salloc

  integer(int32)  :: i, j, l, c, j2, nbad, il(1)
  real(dp), allocatable :: lp(:), s1(:), s2(:), hee(:), h11(:), h22(:), h12(:), &
                           h1e(:), h2e(:)
  integer(int32), allocatable :: ibad(:)

  if(debug > 0) then
    call intpr('n,k,p,ic,nai,x,o,w,pen,a,b,ha,d,r,c,ai', 38, &
      [n, k, p, intcens, nai, size(x), size(offset), size(wt), size(penmat), &
       size(alpha), size(beta), size(ha), size(d), size(row), size(col), size(ai)],   16)
    call dblepr('alpha',  5, alpha,  k)
    call dblepr('beta',   4, beta,   p)
    call dblepr('penmat', 6, penmat, p * p)
  end if

  allocate(lp(n), s1(n), s2(n), hee(n), h11(n), h22(n), h12(n), h1e(n), h2e(n), &
           stat=salloc)

  if(salloc /= 0) return

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
  !
  ! ia, ia2, sgn, ib, nb are now supplied by ormidx (called once per fit)
  ! rather than computed here on every call.
  !
  ! Probability element d, with x1 = alpha(ia) + lp and x2 = alpha(ia2) + lp:
  !
  !   d = F(x1) - F(x2)       ia2 > 0  (0 < Y < k uncensored, or interval censored
  !                                     with two intercepts; sgn = +1)
  !   d = F(x1)               ia2 = 0, sgn = +1
  !   d = 1 - F(x1)           ia2 = 0, sgn = -1
  !
  ! (The earlier code formed d = [s = -1] + s * p1 - p2 with p1 = F(x1), p2 = F(x2),
  ! s = sgn.)  ormatoms never forms 1 - F or F(x1) - F(x2) by subtracting rounded
  ! numbers.  1 - F comes from flipping x (G(x) = F(-x) for the logistic, probit
  ! and Cauchy; loglog <-> cloglog), and the two-intercept d uses the exact width
  ! alpha(ia) - alpha(ia2) (orm_links::celldiff).  So d is accurate however narrow
  ! the cell or far into the tail x1 and x2 are.
  !
  ! Log likelihood:  ell_i = wt_i * log(d_i),
  !                  logL  = -2 * sum(ell_i) + beta' penmat beta

  ! d, and (per what) the score and Hessian coefficients, for all observations.
  ! lp is the common linear predictor.  The y-dependent parts e1 and e2 are omitted
  ! (0) until y-dependent effects exist.
  call ormatoms(link, n, k, alpha, ia, ia2, sgn, lp, wt, what, &
                d, s1, s2, hee, h11, h22, h12, h1e, h2e, salloc)

  if(debug > 0) then
    call intpr('ia',           2, ia,          size(ia))
    call intpr('ia2',          3, ia2,         size(ia2))
    call dblepr('alpha(ia)',   9, alpha(ia),   size(ia))
    call dblepr('sgn',         3, sgn,         size(sgn))
    call dblepr('lp',          2, lp,          size(lp))
  end if

  if(salloc == 999_int32) then     ! zero, negative or NaN probability element
    if(debug > 0) then
      nbad = count(.not. (d > 0_dp))
      allocate(ibad(nbad))
      ibad = pack([(i, i=1,n)], .not. (d > 0_dp))
      call intpr('Zero or negative probability for observations ', 45, ibad, nbad)
      call intpr('Intercept involved', 18, ia(ibad), nbad)
      if(any(ia2(ibad) > 0)) call intpr('2nd Intercept involved', 22, ia2(ibad), nbad)
      call dblepr('d',   1, d  (ibad),   nbad)
      call dblepr('sgn', 3, sgn(ibad),   nbad)
      deallocate(ibad)
    end if
    return
  end if
  if(salloc /= 0_int32) return

  if(debug > 0) then
    call dblepr('alpha', 5, alpha, k)
    call dblepr('beta',  4, beta,  p)
    call dblepr('sgn',   3, sgn,   n)
    call dblepr('d',     1, d,     size(d))
  end if

  logL = -2_dp * wlogsum(n, wt, d) +  dot_product(beta, matmul(penmat, beta))

  ! Outputs that the requested level does not compute are left as passed in (R passes
  ! zeros) rather than being zeroed here: what = 1 returns logL and d only, what = 2
  ! also u, what = 3 also ha/ai, hb, hab.
  if(what == 1) return
  u = 0_dp

  ! Score = gradient of sum(ell), where logL = -2 * sum(ell) + penalty.
  !
  ! Notation: pdf1 = F'(x1), pdf2 = F'(x2), dpdf1 = F''(x1), dpdf2 = F''(x2),
  !           r1 = pdf1 / d, r2 = pdf2 / d.  D() = derivative w.r.t. a parameter.
  !
  !   D log d / D alpha(ia)  =  sgn * pdf1 / d
  !   D log d / D alpha(ia2) = -pdf2 / d                   (only if ia2 > 0)
  !   D log d / D beta_l     =  x_l * (sgn * pdf1 - pdf2) / d    (pdf2 = 0 if ia2 = 0)
  !
  ! (each is pdf(argument) * D(argument) / d, the argument being alpha or alpha + x beta.)
  ! ormatoms returns  s1 = wt * sgn * r1  and  s2 = -wt * r2  (0 if ia2 = 0), so
  !   u(ia) += s1,   u(ia2) += s2,   u(k + l) += x_l * (s1 + s2).
  ! s1 + s2 = D ell / D eta, eta = a common shift of the linear predictor.
  ! The penalty gradient, -penmat beta, is subtracted after the loop.

  nuu = 0_int32

  if(nu > 1) then
    urow = 0_int32
    ucol = 0_int32
    um   = 0_dp
  end if

  do i = 1, n
    j    = ia(i)    ! subscript of applicable alpha
    j2   = ia2(i)   ! subscript of second alpha, 0 if not there
    u(j) = u(j) + s1(i)
    if(j2 > 0) u(j2) = u(j2) + s2(i)
    if(p > 0) then
      do l = 1, p
        u(k + l) = u(k + l) + (s1(i) + s2(i)) * x(i, l)
      end do
    end if

    if(nu > 0) then  ! compute sparse score matrix elements
      if(nuu + 2_int32 + p > nu) then
        salloc = 996_int32
        return
      end if
      nuu = nuu + 1
      urow(nuu) = i
      ucol(nuu) = j
      um(nuu)   = s1(i)
      if(j2 > 0) then
        nuu = nuu + 1
        urow(nuu) = i
        ucol(nuu) = j2
        um(nuu)   = s2(i)
      end if
      if(p > 0) then
        do l = 1, p
          nuu = nuu + 1
          urow(nuu) = i
          ucol(nuu) = k + l
          um(nuu)   = (s1(i) + s2(i)) * x(i, l)
        end do
      end if
    end if

  end do

  ! Add derivative of penalty function -0.5 b'Pb = -Pb
  ! Ignored at present for mscore
  if(p > 0) u((k + 1) : (k + p)) = u((k + 1) : (k + p)) - matmul(penmat, beta)

  if(debug > 0) call dblepr('u', 1, u, k + p)

  if(what == 2) return

  ha  = 0_dp
  hb  = 0_dp
  hab = 0_dp

  ! Hessian
  ! Create 3 matrices: ha (k x 2), hb (p x p), hab (k x p)
  ! ha represents a sparse matrix to be served up to the Matrix package in R
  ! It is tri-band diagonal, and since it is symmetric we only need to compute
  ! two bands (diagonal ha[, 1] and superdiagonal ha[, 2] with ha[k, 2] irrelevant.)
  ! For derivation of the Hessian components see
  ! https://fharrell.com/post/mle#sec-hessformula
  ! If there are any interval-censored observatons, ha is instead computed
  ! as ai using general sparse form

  ! Second derivatives of ell = wt * log(d).  Let a1 = alpha(ia), a2 = alpha(ia2),
  ! x1 = a1 + eta, x2 = a2 + eta (eta = lp), r = pdf / d, dpdf = F''.
  !
  ! One intercept (ia2 = 0):
  !   h11 = h1e = hee = wt * ( sgn * dpdf1 / d - r1^2 )        (h22 = h12 = h2e = 0)
  ! Two intercepts (ia2 > 0, sgn = +1):
  !   h11 =  wt * ( dpdf1 / d - r1^2 )                         d2 ell / d a1^2
  !   h22 = -wt * ( dpdf2 / d + r2^2 )                         d2 ell / d a2^2
  !   h12 =  wt * r1 * r2                                      d2 ell / d a1 d a2
  !   h1e =  wt * ( dpdf1 / d - r1 * (r1 - r2) )               d2 ell / d a1 d eta
  !   h2e = -wt * ( dpdf2 / d - r2 * (r1 - r2) )               d2 ell / d a2 d eta
  !   hee =  wt * ( (dpdf1 - dpdf2) / d - (r1 - r2)^2 )        d2 ell / d eta^2
  ! Since D/D eta = D/D a1 + D/D a2, every link satisfies
  !   h1e = h11 + h12,   h2e = h12 + h22,   hee = h11 + 2 h12 + h22.
  ! (This assumes both terms share one linear predictor.  With y-dependent effects,
  ! x1 = a1 + lp + e1 and x2 = a2 + lp + e2 with e1 = z*gamma(ia), e2 = z*gamma(ia2);
  ! the accumulation then uses h11, h12, h22 with the two design vectors instead of
  ! hee, h1e, h2e, and ormatoms should be given e1 and e2 separately (not pre-added
  ! into lp) so it can form the cell width from them.  See the ormatoms header in
  ! orm_links.f90.)
  !
  ! Logistic: dpdf = pdf * (1 - 2F), and everything collapses to
  !   one intercept :  h11 = h1e = hee = -wt * pdf1
  !   two intercepts:  h1e = -wt * pdf1,  h2e = -wt * pdf2,  hee = -wt * (pdf1 + pdf2),
  !                    h12 = wt * pdf1 * pdf2 / d^2,  h11 = h1e - h12,  h22 = h2e - h12
  ! These have no cancellation.  The general formulas subtract numbers of size
  ! 1/width^2 when the cell is narrow (r1, r2 ~ 1/width), so for links 2-5 h1e and hee
  ! carry absolute errors per unit weight of about 1e-3 at width 1e-5 and 1e-7 at width
  ! 1e-3.  That affects Newton's path and standard errors, not logL or the score.

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

  if(debug > 1) then
    call dblepr('hee', 3, hee, size(hee))
  end if

  ! h11, h22, h12, h1e, h2e, hee are the second derivatives defined above, returned
  ! by ormatoms one observation at a time.  For one-intercept observations
  ! h1e = hee = h11 and the rest are 0.  They are accumulated into the Hessian of
  ! sum(ell), i.e. of -logL/2 without the penalty (orm.fit negates it to form the
  ! information matrix):
  !   ha(j, 1) or ai(j)      += h11         diagonal, alpha part
  !   ha(j2, 1) or ai(j2)    += h22
  !   ha(j, 2) or ai(k+j2-1) += h12         super diagonal (j2 = j + 1)
  !   new (j, j2) element of ai += h12      non-adjacent intercepts, interval censoring only
  !   hab(j, l)  += x(i, l) * h1e           hab(j2, l) += x(i, l) * h2e
  !   hb(l, c)   += x(i, l) * x(i, c) * hee
  ! (the penalty, -penmat, is added to hb at the end when penhess > 0.)

  do i = 1, n
    j  = ia(i)     ! intercept involved in p1
    j2 = ia2(i)    ! intercept involved in p2 (0 if not there)
    ! Intercept-only part of hessian, starting with the no-interval-censored case
    if(intcens == 0) then
      ha(j, 1) = ha(j, 1) + h11(i)
      if(j2 /= 0) then
        ha(j2, 1) = ha(j2, 1) + h22(i)   ! diagonal
        ha(j,  2) = ha(j, 2)  + h12(i)   ! super diagonal
      end if
    else                 ! some interval censored observations
      ai(j) = ai(j) + h11(i)  ! diagonal element; place in ai already allocated
      if(j2 > 0) then    ! second intercept involved
        ai(j2) = ai(j2) + h22(i)                  ! diagonal element
        if(j2 == (j + 1)) then   ! adjacent intercepts involved and place already allocated
          ai(k + j2 - 1) = ai(k + j2 - 1) + h12(i) ! super diagonal
        else
          ! involves 2 non-adjacent intercepts; may have to allocate new position in ai
          il = findloc(row(1 : ne), j, mask = (col(1 : ne) == j2))
          l  = il(1)
          if(l > 0) then
            ai(l) = ai(l) + h12(i)
          else   ! add a new entry
            ne = ne + 1
            if(ne > nai) then
              salloc = 997
              return
            end if
            row(ne) = j
            col(ne) = j2
            ai(ne)  = h12(i)
          end if
        end if
      end if
    end if

    if(p == 0) cycle

    do l = 1, p
      hab(j, l) = hab(j, l) + x(i, l) * h1e(i)
      if(j2 > 0) hab(j2, l) = hab(j2, l) + x(i, l) * h2e(i)
      do c = l, p
        hb(l, c) = hb(l, c) + x(i, l) * x(i, c) * hee(i)
      end do
    end do
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

return

end subroutine ormll
