subroutine ormidx(n, k, y, y2, ia, ia2, sgn, ib, nb, salloc)
! Precomputes, once per fit, the per-observation intercept-index
! bookkeeping used by ormll and (in a later addition) ormeta. This is a
! pure function of y, y2, k -- independent of alpha, beta, offset -- so
! it need be called once per Newton-Raphson procedure, not once per
! iteration. Extracted verbatim from ormll.f90 (indexing block).
!
! Claude Sonnet 5 2026-08-30          entire file
!
! n      : # observations
! k      : # intercepts
! y      : integer outcome vector, values 0:k, or -1 for left censored
! y2     : y2 = k+1 for right censored; y=y2 => uncensored;
!          0 <= y < y2 <= k for interval-censored
! ia     : returned intercept index for the "p1" probability component
! ia2    : returned intercept index for "p2" (0 if none)
! sgn    : returned +1 if p1 enters as F(), -1 if as 1-F()
! ib     : returned observation numbers having a p2 component;
!          only ib(1:nb) is meaningful -- caller must dimension ib >= n
! nb     : returned count of such observations
! salloc : 0 on success, 998 if a censoring configuration is unhandled
!
! NOTE: the two masks below intentionally match ormll.f90's original
! code exactly, including an existing asymmetry between the nb-count
! mask (y2 <= k) and the ib-pack mask (y2 < k). Preserved as-is pending
! confirmation this was intentional in the original implementation.

  use, intrinsic :: ISO_FORTRAN_ENV, only: dp => real64, int32
  implicit none
  integer(int32), intent(in)  :: n, k, y(n), y2(n)
  integer(int32), intent(out) :: ia(n), ia2(n), ib(n), nb, salloc
  real(dp),       intent(out) :: sgn(n)

  integer(int32) :: i, a, b

  salloc = 0

  ! Claude Sonnet 5 2026-08-30
  ! Second, deeper fix on top of the y2<=k -> y2<k correction above:
  ! both masks' second clause used "y >= 0", which incorrectly also
  ! matches the (a==0, b<k) interval-censoring branch below (lower
  ! bound at the very bottom category) -- that branch sets ONLY ia
  ! (with sgn=-1), structurally identical to left-censoring, and never
  ! sets ia2. The only branches that ever produce ia2>0 are the
  ! uncensored-interior case and the (a>0, b<k, a<b) interval-censored
  ! case -- both require a>0 (y>0), not merely a>=0. Confirmed via a
  ! concrete failure: with the old "y>=0" mask, ib included
  ! observations whose actual ia2 was 0, and code consuming ib assumed
  ! every ib-indexed observation has ia2>0 (as documented), silently
  ! producing wrong results (a length-mismatch recycling warning
  ! downstream, in this case) rather than an obvious crash.
  nb = count((y > 0 .and. y < k .and. y == y2) .or. & ! uncensored Y, 0 < Y < k
             (y > 0 .and. y2 <  k .and. y /= y2))    ! or interval censored

  ib(1:nb) = pack([(i, i=1,n)], (y > 0 .and. y < k .and. y == y2) .or. &
                                (y > 0 .and. y2 < k .and. y /= y2)      )

  ! Compute ia : which alpha is involved (first alpha if 0 < y < k & uncensored)
  !         ia2: second alpha involved, will always have negated F(); ia2=0 if no second alpha
  !         ia goes with p1, ia2 goes with p2, take p1 - p2
  !         sgn: 1.0 when prob is F(), -1.0 when prob is 1 - F()

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

end subroutine ormidx
