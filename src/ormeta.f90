subroutine ormeta(n, k, link, alpha, lp1, lp2, wt, ia, ia2, sgn, ib, nb, &
                   logd, g, h, s1, s2, what, debug, salloc)
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
! Internally orm_links::ormatoms takes the common part lp and the y-dependent
! parts e1, e2 separately, so that the width of a cell alpha(ia) - alpha(ia2) +
! (e1 - e2) can be formed from small terms.  Here lp = lp1 and e2 = lp2 - lp1
! (all zero today).  When YDE is added, the R caller should pass the common
! part and the yde parts instead of pre-summed lp1, lp2, and this routine
! should pass them on to ormatoms as lp, e1, e2.
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
! what   : which outputs to compute, a sum of
!            1 = logd
!            2 = g, s1, s2
!            4 = h
!          so 7 = everything (the behavior before this argument existed).  Outputs
!          not requested are not touched: they keep the values passed in (R passes
!          zeros).  Exception: s1 and s2 are by-products of h and are filled whenever
!          2 or 4 is set.  The call sites need
!            log-likelihood at trial nodes (agqStep, surrogate) : 1
!            random-effect mode finding (clusterModeFind)        : 6  (g and h)
!            missing-information correction                       : 3  (logd, g, s1, s2)
!          Every level computes d for all observations, checks it, and returns
!          salloc = 999 exactly as before; asking for less skips the derivative
!          arithmetic (and, without bit 1, the n log() calls).
! debug  : 1 to print bad-probability diagnostics, 0 otherwise
! salloc : 0 on success, matching ormll's convention:
!          > 0 if dynamic array allocation failed, 995 if what is not 1..7,
!          999 if a negative or zero Y=j probability was encountered
!          (signals the caller to step-halve, same as ormll's what=1 path)

  use, intrinsic :: ISO_FORTRAN_ENV, only: dp => real64, int32
  use orm_links, only: ormatoms
  implicit none
  integer(int32), intent(in)  :: n, k, link, nb, what, debug, ia(n), ia2(n), ib(nb)
  real(dp),       intent(in)  :: alpha(k), lp1(n), lp2(n), wt(n), sgn(n)
  real(dp),       intent(inout) :: logd(n), g(n), h(n), s1(n), s2(n)
  integer(int32), intent(out) :: salloc

  integer(int32) :: i, nbad, lvl
  real(dp), allocatable       :: d(:), e2(:)
  integer(int32), allocatable :: ibad(:)

  ! Claude Sonnet 5 2026-09-19
  ! d, s1, s2 and the eta^2 curvature h now come straight from
  ! orm_links::ormatoms (all link arithmetic lives there).  ib and nb are no
  ! longer used but remain in the argument list so callers are unchanged.
  salloc = 0
  if(what < 1_int32 .or. what > 7_int32) then
    salloc = 995_int32
    return
  end if
  allocate(d(n), stat=salloc)
  if(salloc /= 0) return

  ! ormatoms levels are nested (1 = d, 2 = + s1 s2, 3 = + h); map the requested bits
  lvl = 1_int32
  if(iand(what, 2_int32) /= 0_int32) lvl = 2_int32
  if(iand(what, 4_int32) /= 0_int32) lvl = 3_int32

  ! d, and the derivatives of ell = wt * log(d), from ormatoms.  d is F(x1) - F(x2),
  ! F(x1) or 1 - F(x1) exactly as in ormll (see the comments there), except
  ! x1 = alpha(ia) + lp1 and x2 = alpha(ia2) + lp2 may use different linear predictors.
  !   logd = wt * log(d)          summed by cluster = cluster log likelihood at this trial shift
  !   g    = s1 + s2 = D ell / D eta:   s1 = wt * sgn * pdf1 / d,  s2 = -wt * pdf2 / d
  !          (ormll's u contribution, with eta = 1 instead of x(i,l))
  !   h    = D2 ell / D eta2 = hee, ormll's hb multiplier
  !          one intercept: wt * (sgn * dpdf1 / d - (pdf1/d)^2)
  !          two intercepts: wt * ((dpdf1 - dpdf2)/d - ((pdf1 - pdf2)/d)^2)
  !          logistic: -wt * pdf1  or  -wt * (pdf1 + pdf2)
  ! ia2 > 0 always has sgn = +1 (ormidx sets sgn = -1 only in one-intercept cases), so
  ! the two-intercept formulas omit sgn.
  if(all(lp1 == lp2)) then                  ! always, until YDE exists
    call ormatoms(link, n, k, alpha, ia, ia2, sgn, lp1, wt, lvl, &
                  d, s1, s2, hee=h, salloc=salloc)
  else
    allocate(e2(n), stat=salloc)
    if(salloc /= 0) return
    e2 = lp2 - lp1
    call ormatoms(link, n, k, alpha, ia, ia2, sgn, lp1, wt, lvl, &
                  d, s1, s2, hee=h, salloc=salloc, e2=e2)
  end if

  if(salloc == 999_int32) then
    if(debug > 0) then
      nbad = count(.not. (d > 0_dp))
      allocate(ibad(nbad))
      ibad = pack([(i, i=1,n)], .not. (d > 0_dp))
      call intpr('ormeta: zero or negative probability, observations', 51, ibad, nbad)
      call intpr('Intercept involved',                                  18, ia(ibad),  nbad)
      if(any(ia2(ibad) > 0)) &
        call intpr('2nd Intercept involved',                            22, ia2(ibad), nbad)
      deallocate(ibad)
    end if
    return
  end if
  if(salloc /= 0) return

  if(iand(what, 1_int32) /= 0_int32) logd = wt * log(d)
  if(iand(what, 2_int32) /= 0_int32) g    = s1 + s2

end subroutine ormeta
