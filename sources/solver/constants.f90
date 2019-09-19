module constants

  implicit none

  save

  ! ------------------------------------
  ! This module holds global constants.
  ! ------------------------------------

  ! PRECISION
  integer, parameter :: dp = KIND(0.0d0)

  ! GLOBAL CONSTANTS
  real(kind=dp), parameter :: bignum    = huge( 1.0d0 )
  real(kind=dp), parameter :: macheps   = epsilon( 1.0d0 )
  real(kind=dp), parameter :: macheps12 = sqrt( macheps )
  real(kind=dp), parameter :: macheps13 = macheps ** ( 1.0d0 / 3.0d0 )
  real(kind=dp), parameter :: macheps23 = macheps ** ( 2.0d0 / 3.0d0 )
  real(kind=dp), parameter :: macheps34 = macheps ** ( 3.0d0 / 4.0d0 )

end module constants
