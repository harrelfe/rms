module orm_links
! Link-function primitives shared by ormll and ormeta, and ormatoms, which
! turns them into the per-observation likelihood, score and Hessian pieces.
! All link-specific arithmetic in the package lives in this file.
!
! Links: 1 logistic, 2 probit, 3 loglog  F = exp(-exp(-x)),
!        4 cloglog F = 1 - exp(-exp(x)),  5 Cauchy
!
! Two ideas make everything accurate for any x, however far into a tail:
!
!  1. Never form a SMALL 1 - F (or F) by subtraction.  linkpieces returns F,
!     G = 1 - F and the pdf, each with full relative accuracy.  The small side
!     comes from FLIPPING x:
!       logistic, probit, Cauchy :  G(x) = F(-x)
!       loglog  <-> cloglog      :  G_loglog(x) = F_cloglog(-x) and vice versa
!     evaluated with erfc / expm1-type forms.  For probit and Cauchy only the
!     small side is computed (one erfc or atan); the large side is 1 minus it,
!     exact to a rounding because it is at least 0.5.  A one-intercept cell
!     probability is F or G directly.
!
!  2. Never form F(x1) - F(x2) from two rounded arguments when the cell is
!     narrow.  celldiff takes the exact width delta = alpha(ia) - alpha(ia2)
!     (+ e1 - e2, the y-dependent parts of the two linear predictors) and uses
!     an exact difference formula where one exists
!     (logistic, loglog, cloglog, Cauchy).  For probit, which has none, it
!     differences the small-valued tails (G when x2 > 0, else F), which removes
!     the collapse in the far tail but keeps a cancellation error of about
!     eps/delta.
!
! Claude Sonnet 5 2026-09-19 : rewritten (was cdf, pdf, dpdf array primitives)

  use, intrinsic :: ISO_FORTRAN_ENV, only: dp => real64, int32
  implicit none
  private
  public :: linkpieces, linkdpdf, celldiff, ormatoms, wlogsum

  real(dp), parameter :: pi     = 3.14159265358979323846_dp
  real(dp), parameter :: sqrt2  = 1.41421356237309504880_dp
  real(dp), parameter :: isq2pi = 0.39894228040143267794_dp    ! 1 / sqrt(2 pi)
  integer(int32), parameter :: NEEDF = 1_int32, NEEDG = 2_int32, NEEDP = 4_int32, &
                               NEEDALL = 7_int32

contains

  ! expm1(x) = exp(x) - 1 without cancellation.  Fortran has no expm1 intrinsic;
  ! this keeps the package free of any C-library dependency.
  !   |x| < 0.02 : Taylor series through x^9/10!, error below 1e-21 relative, and
  !                 no exp or log call.  This is the usual case for a narrow cell
  !                 (x = width) and for tiny tail probabilities.
  !   otherwise   : Kahan's (u - 1) x / log(u) with u = exp(x), which is exact to a
  !                 couple of ulps; expm1u takes u when the caller already has it.
  pure function expm1s(x) result(r)
    real(dp), intent(in) :: x
    real(dp) :: r
    r = x * (1.0_dp + x / 2.0_dp * (1.0_dp + x / 3.0_dp * (1.0_dp + x / 4.0_dp * &
        (1.0_dp + x / 5.0_dp * (1.0_dp + x / 6.0_dp * (1.0_dp + x / 7.0_dp * &
        (1.0_dp + x / 8.0_dp * (1.0_dp + x / 9.0_dp))))))))
  end function expm1s

  pure function expm1u(x, u) result(r)
    real(dp), intent(in) :: x, u
    real(dp) :: r, um1
    if(abs(x) < 0.02_dp) then
      r = expm1s(x)
    else if(u > huge(u)) then
      r = u
    else
      um1 = u - 1.0_dp
      if(um1 == -1.0_dp) then
        r = -1.0_dp
      else
        r = um1 * x / log(u)
      end if
    end if
  end function expm1u

  pure function expm1(x) result(r)
    real(dp), intent(in) :: x
    real(dp) :: r
    if(abs(x) < 0.02_dp) then
      r = expm1s(x)
    else
      r = expm1u(x, exp(x))
    end if
  end function expm1

  ! F = cdf, g = 1 - cdf (small side never by subtraction), p = pdf, and for the
  ! extreme-value links t = exp(-x) (loglog) or exp(x) (cloglog), which linkdpdf
  ! and celldiff reuse.  need is a sum of NEEDF, NEEDG, NEEDP; pieces that are not
  ! requested may be returned as 0 (some come for free and are always returned).
  pure subroutine linkpieces(x, link, need, f, g, p, t)
    real(dp),       intent(in)  :: x
    integer(int32), intent(in)  :: link, need
    real(dp),       intent(out) :: f, g, p, t
    real(dp) :: e, den

    f = 0.0_dp ; g = 0.0_dp ; p = 0.0_dp ; t = 0.0_dp
    select case(link)
    case(1)                                    ! logistic
      e   = exp(-abs(x))
      den = 1.0_dp / (1.0_dp + e)
      if(x >= 0.0_dp) then
        f = den ;      g = e * den
      else
        f = e * den ;  g = den
      end if
      p = f * g
    case(2)                                    ! probit: one erfc, the other side is 1 - it
      if(x >= 0.0_dp) then
        g = 0.5_dp * erfc(x / sqrt2)
        f = 1.0_dp - g
      else
        f = 0.5_dp * erfc(-x / sqrt2)
        g = 1.0_dp - f
      end if
      if(iand(need, NEEDP) /= 0) p = isq2pi * exp(-0.5_dp * x * x)
    case(3)                                    ! loglog: F = exp(-exp(-x))
      if(x < -708.0_dp) then
        f = 0.0_dp ; g = 1.0_dp ; t = huge(t)
      else
        t = exp(-x)
        f = exp(-t)
        if(iand(need, NEEDG) /= 0) g = -expm1u(-t, f)
        if(iand(need, NEEDP) /= 0) p = t * f
      end if
    case(4)                                    ! cloglog: F = 1 - exp(-exp(x))
      if(x > 708.0_dp) then
        f = 1.0_dp ; g = 0.0_dp ; t = huge(t)
      else
        t = exp(x)
        g = exp(-t)
        if(iand(need, NEEDF) /= 0) f = -expm1u(-t, g)
        if(iand(need, NEEDP) /= 0) p = t * g
      end if
    case(5)                                    ! Cauchy: one atan, the other side is 1 - it
      if(x > 0.0_dp) then
        g = atan(1.0_dp / x) / pi
        f = 1.0_dp - g
      else if(x < 0.0_dp) then
        f = atan(-1.0_dp / x) / pi
        g = 1.0_dp - f
      else
        f = 0.5_dp ; g = 0.5_dp
      end if
      if(iand(need, NEEDP) /= 0) p = 1.0_dp / (pi * (1.0_dp + x * x))
    end select
  end subroutine linkpieces

  ! Derivative of the pdf, links 2-5, given x, the pdf p and t from linkpieces.
  ! (The logistic never needs it: its derivatives have closed forms in ormatoms.)
  pure function linkdpdf(x, p, t, link) result(q)
    real(dp),       intent(in) :: x, p, t
    integer(int32), intent(in) :: link
    real(dp) :: q
    q = 0.0_dp
    if(p == 0.0_dp) return
    select case(link)
    case(2)
      q = - x * p
    case(3)
      q = p * expm1u(-x, t)                    ! pdf * (exp(-x) - 1),  t = exp(-x)
    case(4)
      q = - p * expm1u(x, t)                   ! pdf * (1 - exp(x)),   t = exp(x)
    case(5)
      q = -2.0_dp * pi * x * p * p             ! 1 / (1 + x^2) = pi p
    end select
  end function linkdpdf

  ! sinh(z) for the half-width z = dl/2 of a cell.  |z| < 0.02: Taylor series through
  ! z^9/9!, relative error below 1e-21 and no library call (a narrow cell is the usual
  ! case for continuous Y).  Otherwise the intrinsic.
  pure function sinhh(z) result(r)
    real(dp), intent(in) :: z
    real(dp) :: r, z2
    if(abs(z) < 0.02_dp) then
      z2 = z * z
      r  = z * (1.0_dp + z2 / 6.0_dp * (1.0_dp + z2 / 20.0_dp * (1.0_dp + z2 / 42.0_dp * &
           (1.0_dp + z2 / 72.0_dp))))
    else
      r = sinh(z)
    end if
  end function sinhh

  ! d = F(x1) - F(x2), x1 = x2 + dl, dl the EXACT width.  Inputs are the
  ! linkpieces of x1 and x2.  Exact-difference forms:
  !   logistic : 2 sinh(dl/2) sqrt(pdf1 pdf2)   (series for narrow cells, see sinhh)
  !   loglog   : F1 (-expm1(-t1 expm1(dl))),  t1 = exp(-x1)
  !   cloglog  : G2 (-expm1(-t2 expm1(dl))),  t2 = exp(x2)
  !   Cauchy   : atan2(dl, 1 + x1 x2) / pi
  ! Otherwise (probit; and for dl > 30 where sinh/expm1 could overflow) the
  ! small-valued tails are differenced: G2 - G1 when x2 > 0, else F1 - F2.
  pure function celldiff(link, dl, x1, x2, f1, g1, f2, g2, p1, p2, t1, t2) result(dd)
    integer(int32), intent(in) :: link
    real(dp),       intent(in) :: dl, x1, x2, f1, g1, f2, g2, p1, p2, t1, t2
    real(dp) :: dd

    if(link == 5_int32) then                   ! Cauchy: atan(x1) - atan(x2) = atan2(dl, 1 + x1 x2)
      dd = atan2(dl, 1.0_dp + x1 * x2) / pi
      return
    end if
    if(dl <= 30.0_dp) then
      select case(link)
      case(1)
        dd = 2.0_dp * sinhh(0.5_dp * dl) * sqrt(p1) * sqrt(p2) ;  return
      case(3)
        dd = f1 * (-expm1(-t1 * expm1(dl))) ;                    return
      case(4)
        dd = g2 * (-expm1(-t2 * expm1(dl))) ;                    return
      end select
    end if
    if(x2 > 0.0_dp) then
      dd = g2 - g1
    else
      dd = f1 - f2
    end if
  end function celldiff

  ! sum(wt * log(d)) accumulated in blocks of 512: a plain running sum over n = 250,000
  ! terms of size ~10 (total ~3e6) loses about 1.8e-8 absolutely; summing blocks
  ! first reduces that to about 1e-9, i.e. the resolution of a double at that size.
  pure function wlogsum(n, wt, d) result(s)
    integer(int32), intent(in) :: n
    real(dp),       intent(in) :: wt(n), d(n)
    real(dp) :: s, t
    integer(int32) :: i, i0, i1
    s = 0.0_dp
    do i0 = 1_int32, n, 512_int32
      i1 = min(i0 + 511_int32, n)
      t  = 0.0_dp
      do i = i0, i1
        t = t + wt(i) * log(d(i))
      end do
      s = s + t
    end do
  end function wlogsum

  ! Likelihood-only evaluation (what = 1): the same d as celldiff / linkpieces, but
  ! each link computes only the pieces it uses.  Loglog needs only point 1, cloglog
  ! only point 2, Cauchy neither (atan2 form); probit and logistic need both, at
  ! one transcendental call each.  Because linkpieces returns identical values
  ! whatever is requested, d here is bit-identical to d from what >= 2.
  pure function oned(link, sgn, x) result(dd)
    integer(int32), intent(in) :: link
    real(dp),       intent(in) :: sgn, x
    real(dp) :: dd, f, g, p, t
    if(sgn == 1.0_dp) then
      call linkpieces(x, link, NEEDF, f, g, p, t)
      dd = f
    else
      call linkpieces(x, link, NEEDG, f, g, p, t)
      dd = g
    end if
  end function oned

  pure function celld(link, dl, x1, x2) result(dd)
    integer(int32), intent(in) :: link
    real(dp),       intent(in) :: dl, x1, x2
    real(dp) :: dd, f1, g1, p1, t1, f2, g2, p2, t2

    f1 = 0.0_dp ; g1 = 0.0_dp ; p1 = 0.0_dp ; t1 = 0.0_dp
    f2 = 0.0_dp ; g2 = 0.0_dp ; p2 = 0.0_dp ; t2 = 0.0_dp
    if(link /= 5_int32 .and. dl > 30.0_dp) then      ! wide cell: tails of both points
      call linkpieces(x1, link, NEEDF + NEEDG, f1, g1, p1, t1)
      call linkpieces(x2, link, NEEDF + NEEDG, f2, g2, p2, t2)
    else
      select case(link)
      case(1, 2)
        call linkpieces(x1, link, NEEDP, f1, g1, p1, t1)
        call linkpieces(x2, link, NEEDP, f2, g2, p2, t2)
      case(3)
        call linkpieces(x1, link, NEEDF, f1, g1, p1, t1)
      case(4)
        call linkpieces(x2, link, NEEDG, f2, g2, p2, t2)
      end select
    end if
    dd = celldiff(link, dl, x1, x2, f1, g1, f2, g2, p1, p2, t1, t2)
  end function celld

  ! ---------------------------------------------------------------------
  ! ormatoms: probability element d and the derivatives of ell = wt * log(d).
  !
  ! Observation i has  d = F(x1) - F(x2)   (ia2 > 0, sgn = +1), or
  !                    d = F(x1) or 1 - F(x1)   (ia2 = 0, sgn = +1 or -1),
  !   x1 = alpha(ia)  + lp + e1,   x2 = alpha(ia2) + lp + e2.
  ! lp is the linear predictor common to every cumulative probability (offset +
  ! X beta + any random effect); e1 and e2 are the parts specific to the two
  ! intercepts, i.e. the y-dependent effects z * gamma(ia) and z * gamma(ia2).
  ! Derivatives are with respect to a1 = alpha(ia), a2 = alpha(ia2), and a
  ! common shift eta of both linear predictors.
  !
  ! link   : 1-5            n, k : # observations, # intercepts
  ! alpha  : k intercepts (descending)
  ! ia, ia2, sgn : from ormidx (ia2 = 0 if none)
  ! lp     : n-vector, the common linear predictor
  ! wt     : case weights
  ! what   : 1 = d only; 2 = also s1, s2; 3 = also hee and the optional h*
  ! d      : returned probability element (unweighted); always assigned
  ! s1, s2, hee and the optional h* are assigned only for the observations and levels
  ! of what that call for them; otherwise they are left exactly as passed in.
  ! s1, s2 : d ell / d a1, d ell / d a2 (s2 = 0 if ia2 = 0); s1 + s2 = d ell / d eta
  ! hee    : d2 ell / d eta2
  ! h11, h22, h12 : d2 ell / d a1^2, d a2^2, d a1 d a2
  ! h1e, h2e      : d2 ell / d a1 d eta, d a2 d eta      (h11..h2e optional)
  ! salloc : 0, or 999 if some d is not strictly positive or is NaN (d is still
  !          filled in for every observation so callers can report offenders)
  ! e1, e2 : optional n-vectors, the y-dependent parts of the linear predictors
  !          for the first and second intercept (default 0; e2 is read only when
  !          ia2 > 0).  The width of the cell is delta = alpha(ia) - alpha(ia2) +
  !          (e1 - e2), formed HERE from the small terms.  Callers must pass the
  !          parts separately and must not pre-add them into lp: two numbers
  !          lp + e1 and lp + e2 are each rounded to about ulp(lp) (2e-15 for
  !          |lp| near 10), so their difference has lost the information that
  !          matters, roughly 2e-10 relative for a cell of width 1e-5.  Only the
  !          exact-width forms use the width (logistic, loglog, cloglog, Cauchy);
  !          probit differences the tails and cannot benefit.
  ! With e1 /= e2, s1 and s2 are still D ell / D x1 and D ell / D x2 and h11, h12,
  ! h22 the second derivatives in x1 and x2, which is what yde parameters need
  ! (D x1 / D gamma = z1, D x2 / D gamma = z2).  hee, h1e and h2e describe a
  ! COMMON shift of x1 and x2 (beta, random intercepts) and are valid only for that.
  !
  ! For the logistic the second derivatives have closed forms free of any
  ! cancellation (verified against 450-digit arithmetic):
  !   one intercept :  hee = h11 = h1e = -wt pdf1
  !   two intercepts:  hee = -wt (pdf1 + pdf2),  h1e = -wt pdf1,  h2e = -wt pdf2,
  !                    h12 = wt pdf1 pdf2 / d^2,  h11 = h1e - h12,  h22 = h2e - h12
  ! For links 2-5 they are the direct formulas in pdf and dpdf.
  ! ---------------------------------------------------------------------
  subroutine ormatoms(link, n, k, alpha, ia, ia2, sgn, lp, wt, what, &
                      d, s1, s2, hee, h11, h22, h12, h1e, h2e, salloc, e1, e2)
    integer(int32), intent(in)  :: link, n, k, ia(n), ia2(n), what
    real(dp),       intent(in)  :: alpha(k), sgn(n), lp(n), wt(n)
    real(dp),       intent(out)   :: d(n)
    real(dp),       intent(inout) :: s1(n), s2(n), hee(n)
    real(dp),       intent(inout), optional :: h11(n), h22(n), h12(n), h1e(n), h2e(n)
    integer(int32), intent(out) :: salloc
    real(dp),       intent(in),  optional :: e1(n), e2(n)

    integer(int32) :: i
    logical        :: bad, logistic, has1, has2
    real(dp)       :: x1, x2, f1, g1, p1, t1, f2, g2, p2, t2, dd, dl, w, sg, r1, r2, &
                      q1, q2, z, a12, e1i, e2i

    salloc   = 0_int32
    bad      = .false.
    logistic = link == 1_int32
    has1     = present(e1)
    has2     = present(e2)

    do i = 1, n
      w  = wt(i)
      x1  = alpha(ia(i)) + lp(i)
      e1i = 0.0_dp
      if(has1) then
        e1i = e1(i) ; x1 = x1 + e1i
      end if

      if(what == 1_int32) then                       ! likelihood only: d, nothing else
        if(ia2(i) == 0_int32) then
          dd = oned(link, sgn(i), x1)
        else
          x2  = alpha(ia2(i)) + lp(i)
          e2i = 0.0_dp
          if(has2) then
            e2i = e2(i) ; x2 = x2 + e2i
          end if
          dd = celld(link, (alpha(ia(i)) - alpha(ia2(i))) + (e1i - e2i), x1, x2)
        end if
        d(i) = dd
        if(.not. (dd > 0.0_dp)) bad = .true.
        cycle
      end if
      call linkpieces(x1, link, NEEDALL, f1, g1, p1, t1)

      if(ia2(i) == 0_int32) then                     ! ---- one intercept
        sg = sgn(i)
        dd = merge(f1, g1, sg == 1.0_dp)
        d(i) = dd
        if(.not. (dd > 0.0_dp)) then
          bad = .true. ; cycle
        end if
        if(what < 2) cycle
        s1(i) = w * sg * p1 / dd
        s2(i) = 0.0_dp
        if(what < 3) cycle
        if(logistic) then
          z = - w * p1
        else
          q1 = linkdpdf(x1, p1, t1, link)
          r1 = p1 / dd
          z  = w * (sg * q1 / dd - r1 ** 2)
        end if
        hee(i) = z
        if(present(h11)) then
          h11(i) = z ;      h1e(i) = z
          h22(i) = 0.0_dp ; h12(i) = 0.0_dp ; h2e(i) = 0.0_dp
        end if

      else                                           ! ---- two intercepts
        x2  = alpha(ia2(i)) + lp(i)
        e2i = 0.0_dp
        if(has2) then
          e2i = e2(i) ; x2 = x2 + e2i
        end if
        call linkpieces(x2, link, NEEDALL, f2, g2, p2, t2)
        dl = (alpha(ia(i)) - alpha(ia2(i))) + (e1i - e2i)
        dd = celldiff(link, dl, x1, x2, f1, g1, f2, g2, p1, p2, t1, t2)
        d(i) = dd
        if(.not. (dd > 0.0_dp)) then
          bad = .true. ; cycle
        end if
        if(what < 2) cycle
        r1    = p1 / dd
        r2    = p2 / dd
        s1(i) =   w * r1
        s2(i) = - w * r2
        if(what < 3) cycle
        if(logistic) then
          a12    = w * r1 * r2
          hee(i) = - w * (p1 + p2)
          if(present(h11)) then
            h12(i) = a12
            h1e(i) = - w * p1
            h2e(i) = - w * p2
            h11(i) = h1e(i) - a12
            h22(i) = h2e(i) - a12
          end if
        else
          ! written with the ratios r = pdf/d so nothing forms d**2 (which underflows)
          q1 = linkdpdf(x1, p1, t1, link) / dd
          q2 = linkdpdf(x2, p2, t2, link) / dd
          hee(i) = w * ((q1 - q2) - (r1 - r2) ** 2)
          if(present(h11)) then
            h11(i) =   w * (q1 - r1 ** 2)
            h22(i) = - w * (q2 + r2 ** 2)
            h12(i) =   w * r1 * r2
            h1e(i) =   w * (q1 - r1 * (r1 - r2))
            h2e(i) = - w * (q2 - r2 * (r1 - r2))
          end if
        end if
      end if
    end do

    if(bad) salloc = 999_int32
  end subroutine ormatoms

end module orm_links
