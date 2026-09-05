subroutine ormeta(n, k, link, alpha, lp1, lp2, wt, ia, ia2, sgn, ib, nb, &
                   logd, g, h, s1, s2, debug, salloc)
! Claude Sonnet 5 2026-08-30
! Extended to also return s1, s2: the per-observation score
! contributions to alpha(ia(i)) and alpha(ia2(i)) SEPARATELY, rather
! than only the combined g = s1 + s2. g/h are unchanged and still the
! right (combined) quantities for random-intercept mode-finding and
! the outer Newton step's gradient direction (where both p1 and p2
! move together, since a random intercept shifts eta identically on
! both sides). s1/s2 exist for a different purpose: computing the
! sparse "missing information" correction to the final information
! matrix (Louis 1982), which needs each cluster's own score
! contribution broken out by WHICH alpha index it credits, since a
! cluster's active alpha-index set (at most ~2x its observation
! count, however large k is) is exactly what keeps that correction
! sparse -- independent of k -- for realistic cluster configurations.
!
! s1, s2 are cheap to add: both are algebraic rearrangements of
! pdf1/pdf2/d/sgn/wt, already computed as intermediates for g/h below.
! No new cdf/pdf/dpdf evaluations, no new O(n) passes.
!
! Lean per-observation companion to ormll, for use by orm.rfit's
! per-cluster random-intercept mode-finding (scalar Newton-Raphson per
! cluster) and adaptive Gauss-Hermite node evaluation. Deliberately
! omits everything ormll needs only for the O(k) intercept Hessian and
! O(p^2) beta Hessian -- no x, beta, penmat, ha, hb, hab -- since those
! are held fixed while this routine is called (once per AGQ node, once
! per mode-finding Newton step). Shares cdf/pdf/dpdf with ormll via the
! orm_links module, and shares ia/ia2/sgn/ib/nb (from ormidx) with
! ormll -- both computed once per fit, not recomputed here.
!
! lp1, lp2 are taken as separate arguments (rather than a single lp)
! so this routine needs no change when y-dependent effects (YDE/
! partial PO) are eventually added: the caller builds
!   lp1 = offset + x*beta + [yde term at ia]  (+ trial cluster shift)
!   lp2 = offset + x*beta + [yde term at ia2] (+ trial cluster shift)
! and for a random-intercept-only fit (today), simply passes lp1=lp2.
!
! n      : # observations
! k      : # intercepts
! link   : link function/distribution family, 1-5 (see orm_links module)
! alpha  : k-vector of (fixed, current) intercepts
! lp1    : n-vector, linear predictor for the "p1" probability component
! lp2    : n-vector, linear predictor for the "p2" component (only
!          lp2(ib) is read; irrelevant elsewhere)
! wt     : n-vector of case weights
! ia     : intercept index for the "p1" probability component (from ormidx)
! ia2    : intercept index for "p2" component, 0 if none (from ormidx)
! sgn    : +1 if p1 enters as F(), -1 if as 1-F() (from ormidx)
! ib     : observation numbers having a p2 component; only ib(1:nb) used (from ormidx)
! nb     : count of such observations (from ormidx)
! logd   : returned n-vector, wt(i) * log(d(i)) -- grouped-sum by
!          cluster gives each cluster's log-likelihood at this trial u
! g      : returned n-vector, d(logL)/d(eta_i) -- grouped-sum by cluster
!          gives each cluster's score with respect to its own u_j
!          (since d(lp1_i)/d(u_j) = d(lp2_i)/d(u_j) = 1 for i in cluster j)
! h      : returned n-vector, d2(logL)/d(eta_i^2) -- grouped-sum by
!          cluster gives each cluster's curvature with respect to u_j
! s1     : returned n-vector, per-observation score contribution
!          credited to alpha(ia(i)) alone (g(i) = s1(i) + s2(i))
! s2     : returned n-vector, per-observation score contribution
!          credited to alpha(ia2(i)) alone; 0 where ia2(i) = 0
! debug  : 1 to print bad-probability diagnostics, 0 otherwise
! salloc : 0 on success, matching ormll's convention:
!          > 0 if dynamic array allocation failed,
!          999 if a negative or zero Y=j probability was encountered
!          (signals the caller to step-halve, same as ormll's what=1 path)

  use, intrinsic :: ISO_FORTRAN_ENV, only: dp => real64, int32
  use orm_links, only: cdf, pdf, dpdf
  implicit none
  integer(int32), intent(in)  :: n, k, link, ia(n), ia2(n), ib(nb), nb, debug
  real(dp),       intent(in)  :: alpha(k), lp1(n), lp2(n), wt(n), sgn(n)
  real(dp),       intent(out) :: logd(n), g(n), h(n), s1(n), s2(n)
  integer(int32), intent(out) :: salloc

  integer(int32) :: i, nbad
  real(dp), allocatable       :: p1(:), p2(:), pdf1(:), pdf2(:), dpdf1(:), dpdf2(:), d(:)
  integer(int32), allocatable :: ibad(:)

  salloc = 0
  allocate(p1(n), p2(n), pdf1(n), pdf2(n), dpdf1(n), dpdf2(n), d(n), stat=salloc)
  if(salloc /= 0) return

  ! Probability element for likelihood -- same formula as ormll, but
  ! p1 and p2 are now evaluated at possibly-different linear predictors
  p1  = cdf(alpha(ia) + lp1, link)
  p2  = 0_dp
  if(nb > 0) p2(ib) = cdf(alpha(ia2(ib)) + lp2(ib), link)
  d   = merge(p1 - p2, 1_dp - p1, sgn == 1_dp)

  nbad = count(d <= 0_dp)
  if(nbad > 0_int32) then
    if(debug > 0) then
      allocate(ibad(nbad))
      ibad = pack([(i, i=1,n)], d <= 0_dp)
      call intpr('ormeta: zero or negative probability, observations', 51, ibad, nbad)
      call intpr('Intercept involved',                                  18, ia(ibad),  nbad)
      if(any(ia2(ibad) > 0)) &
        call intpr('2nd Intercept involved',                            22, ia2(ibad), nbad)
      deallocate(ibad)
    end if
    salloc = 999_int32
    deallocate(p1, p2, pdf1, pdf2, dpdf1, dpdf2, d)
    return
  end if

  logd = wt * log(d)

  ! Score coefficient (analogous to ormll's g1 - g2, evaluated at eta = 1
  ! rather than eta = x(i,l), i.e. this IS d(logL)/d(eta_i))
  pdf1     = pdf(alpha(ia)      + lp1,     p1,     link)
  pdf2     = 0_dp
  pdf2(ib) = pdf(alpha(ia2(ib)) + lp2(ib), p2(ib), link)

  s1 = (wt / d) * sgn * pdf1
  s2 = 0_dp
  s2(ib) = -(wt(ib) / d(ib)) * pdf2(ib)
  g  = s1 + s2

  ! Hessian coefficient (analogous to ormll's z / two-intercept hb
  ! multiplier, again evaluated at eta = 1). Branches exactly as ormll
  ! does: when ia2(i) = 0, sgn(i) can be +-1 and must appear; when
  ! ia2(i) > 0, sgn(i) is always +1 by construction in ormidx (only the
  ! single-alpha branches ever set sgn = -1), so the two-intercept
  ! formula below correctly omits it rather than assuming it away.
  dpdf1     = dpdf(alpha(ia)      + lp1,     p1,     pdf1,     link)
  dpdf2     = 0_dp
  dpdf2(ib) = dpdf(alpha(ia2(ib)) + lp2(ib), p2(ib), pdf2(ib), link)

  h = merge( &
        (wt / d**2) * (d * (dpdf1 - dpdf2) - (pdf1 - pdf2)**2), &  ! ia2 > 0
        (wt / d**2) * (sgn * d * dpdf1 - pdf1**2),               &  ! ia2 = 0
        ia2 > 0)

  deallocate(p1, p2, pdf1, pdf2, dpdf1, dpdf2, d)

end subroutine ormeta
