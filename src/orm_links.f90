module orm_links
! Link-function primitives (cdf, pdf, dpdf) shared by ormll and (in a
! later addition) ormeta. Extracted verbatim from ormll.f90's internal
! contains block -- no numerical changes, just made independently
! usable by other Fortran translation units via `use orm_links`.
!
! Claude Sonnet 5 2026-08-30          entire file

  use, intrinsic :: ISO_FORTRAN_ENV, only: dp => real64, int32
  implicit none
  private
  public :: cdf, pdf, dpdf

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

end module orm_links
