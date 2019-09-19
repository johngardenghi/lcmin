module variables

  use constants, only: dp

  implicit none

  save

  ! ------------------------------------
  ! This module holds global variables.
  ! ------------------------------------

  ! LCMIN global parameters
  character(len=80) :: ini_file  ! Output file for initial point
  character(len=80) :: out_file  ! LCMIN output file
  character(len=80) :: sol_file  ! Output file for solution found by LCMIN

  integer :: hnnzmax      ! Maximum size of nonzeros of the Hessian
  integer :: hslunit      ! Stream for HSL output
  integer :: maxhslfix    ! Maximum iterations trying to run an HSL routine
  integer :: maxinnerit   ! Maximum inner iterations
  integer :: maxouterit   ! Maximum outer iterations
  integer :: iniunit      ! Stream for initial point output file
  integer :: outunit      ! Stream for LCMIN output (besides screen)
  integer :: solunit      ! Stream for solution output file

  logical :: scalea     ! Decides whether scaling the constraints
  logical :: scalef     ! Decides whether scaling the objective function
  logical :: printinit  ! Print inner iteration information
  logical :: printtl    ! Print a table line with the results from LCMIN run

  real(kind=dp) :: chiini    ! Initial perturbation on the PD system
  real(kind=dp) :: chimin    ! Minimum value to perturbation on PD system
  real(kind=dp) :: epstol    ! Optimality tolerance
  real(kind=dp) :: epsfeas   ! Feasibility tolerance
  real(kind=dp) :: gamma     ! Armijo constant
  real(kind=dp) :: gmax      ! Constant for constraints and obj function scaling
  real(kind=dp) :: imu       ! The initial value for the barrier parameter
  real(kind=dp) :: infty     ! Numerical representation of infinity
  real(kind=dp) :: kappad    ! Value used in the computation of the barrier function
  real(kind=dp) :: kappaeps  ! Tolerance to optimality error on barrier subproblem
  real(kind=dp) :: kappaip1  ! Constant for perturbation of box constraints
  real(kind=dp) :: kappaip2  ! Constant for perturbation of box constraints
  real(kind=dp) :: kappamu   ! Constant to update the barrier parameter
  real(kind=dp) :: kappaz    ! Safeguard for the approximation of primal Hessian
  real(kind=dp) :: kchibar   ! Constant to compute the perturbation on the PD system
  real(kind=dp) :: kchimin   ! Constant to compute the perturbation on the PD system
  real(kind=dp) :: kchimax   ! Constant to compute the perturbation on the PD system
  real(kind=dp) :: smax      ! Constant to scale the optimality error
  real(kind=dp) :: thetamu   ! Constant to update the barrier parameter
  real(kind=dp) :: tmin      ! Minimum value of fraction to boundary
  real(kind=dp) :: vartol    ! Variation tolerance
  
  ! Problem data scalars
  integer :: nprob      ! Original number of variables
  integer :: nl, nu     ! Number of lower and upper bounded variables

  ! Problem data arrays
  integer, allocatable, dimension(:) :: cind    ! Indexes of the lines that remains in A
  integer, allocatable, dimension(:) :: find    ! New index of each x_i after fixed variables elimination
  integer, allocatable, dimension(:) :: lind    ! Indexes of lower bounded variables
  integer, allocatable, dimension(:) :: lnu     ! Indexes of lower, but not upper, bounded variables
  integer, allocatable, dimension(:) :: nfind   ! Holds the original index of non-fixed vars
  integer, allocatable, dimension(:) :: uind    ! Indexes of upper bounded variables
  integer, allocatable, dimension(:) :: unl     ! Indexes of upper, but not lower, bounded variables

  ! Counters
  integer :: bfevalcnt     ! Barrier function evaluations
  integer :: bgevalcnt     ! Barrier function gradient evaluations
  integer :: fevalcnt      ! Objective function evaluations
  integer :: gevalcnt      ! Objective function gradient evaluations
  integer :: hevalcnt      ! Objective function Hessian evaluations
  integer :: inneraccitcnt ! Total inner iterations
  integer :: inneritcnt    ! Inner iterations in one outer iteration
  integer :: outeritcnt    ! Outer iterations
  integer :: pertnw        ! Perturbations on the NW of the primal-dual system
  integer :: pertse        ! Perturbations on the SE of the primal-dual system
  integer :: smallsd       ! Small search directions

  real :: pretime   ! Time for preprocessing
  real :: opttime   ! Time for optimization

  ! Scaling
  real(kind=dp)                            :: sf  ! Scale factor for the objective function
  real(kind=dp), allocatable, dimension(:) :: sa  ! Scale factors for the matrix of the constraints

end module variables
