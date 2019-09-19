! ------------------------------------------------------------------
! File : toyprob.f90
! ------------------------------------------------------------------

! ------------------------------------------------------------------
! Problem routines definition
! ------------------------------------------------------------------

subroutine inidim(n, m, annz, hnnzmax)

  implicit none
  
  ! SCALAR ARGUMENTS
  integer, intent(out) :: n, m, annz, hnnzmax

  ! =================================================================
  ! This subroutine defines the main scalars arguments to the solver.
  ! Please modify this according to your problem.
  ! =================================================================
  !
  ! Parameters:
  !
  ! On entry: none.
  !
  ! On return:
  !
  ! n         integer,
  !           number of variables.
  !
  ! m         integer,
  !           number of constraints (except the lower and upper
  !           bounds).
  !
  ! annz      integer,
  !           number of nonzero entries on the jacobian of the
  !           constraints.
  !
  ! hnnzmax   integer,
  !           number of nonzero entries on the hessian matrix (or a
  !           upper bound for it).
  ! =================================================================

  n       = 4
  m       = 1
  annz    = 2
  hnnzmax = 3

end subroutine inidim

! ------------------------------------------------------------------

subroutine inip(n, x, l, u, m, b, lambda, annz, alin, acol, aval, flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(out) :: flag
  integer, intent(in)  :: m, n, annz

  ! ARRAY ARGUMENTS
  real(kind=8), dimension(n),    intent(out) :: x, l, u
  real(kind=8), dimension(m),    intent(out) :: b, lambda
  real(kind=8), dimension(annz), intent(out) :: aval
  integer,      dimension(annz), intent(out) :: alin, acol

  ! =================================================================
  ! This subroutine defines some problem data to the solver.
  ! Please modify it according to your problem.
  ! =================================================================
  !
  ! Parameters:
  !
  ! On entry:
  !
  ! n        integer,
  !          number of variables.
  !
  ! m        integer,
  !          number of constraints (except the lower and upper
  !          bounds).
  !
  ! annz     integer,
  !          number of nonzero entries on the jacobian of the constraints.
  !
  ! On return:
  !
  ! x        real, kind of double precision, array of size n,
  !          initial point.
  !
  ! l        real, kind of double precision, array of size n,
  !          lower bounds.
  !
  ! u        real, kind of double precision, array of size n,
  !          upper bounds.
  !
  ! b        real, kind of double precision, array of size m,
  !          right hand side of the constraints.
  !
  ! lambda   real, kind of double precision, array of size m,
  !          initial estimation of Lagrange multipliers.
  !
  ! alin     
  ! acol     integer arrays of size annz,
  !          the kth element of the jacobian of the constraints a_ij
  !          is represented by alin(k) = i and acol(k) = j
  !
  ! aval     real, kind of double precision, array of size annz,
  !          aval(k) contains the value of the jacobian of the
  !          constraints whose indexes were held in alin(k) and
  !          acol(k).
  !
  ! flag     integer,
  !          return 0 if all the settings were successful, or another
  !          value otherwise.
  !
  ! =================================================================

  ! LOCAL SCALAR
  integer :: i

  ! Initial point
  x = (/ - 2.0d0, - 4.0d0, 1.40d+01, 1.0d+01 /)
  ! x = (/ -12.0d0, - 20.0d0, 10.0d0 /)

  ! Lower and upper bounds
  l = (/ -1.0d+01, - 1.0d+01, 0.0d0, 0.0d0 /)
  u = (/  1.0d+01,   1.0d+01, 0.0d0, 0.0d0 /)

  ! Initial estimation of lagrange multipliers
  lambda = (/ (0.0d0, i=1, m) /)

  ! Jacobian of the constraints
  alin = (/     1,     1 /)
  acol = (/     2,     1 /)
  aval = (/ 1.0d0, 2.0d0 /)

  ! alin = (/     1,       1,     2,     2,     3,     3 /)
  ! acol = (/     1,       2,     2,     3,     2,     3 /)
  ! aval = (/ 1.0d0, - 1.0d0, 2.0d0, 1.0d0, 4.0d0, 2.0d0 /)

  b = (/ 2.0d0, 6.0d0 /)

  ! b = (/ 2.0d0, 6.0d0, 11.0d0 /)

  flag = 0

end subroutine inip

! ------------------------------------------------------------------

subroutine evalf(n, x, f, flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)  :: n
  integer,      intent(out) :: flag
  real(kind=8), intent(out) :: f

  ! ARRAY ARGUMENTS
  real(kind=8), dimension(n), intent(in) :: x

  f = (2 * x(1)**2) + (x(2)**3 * x(3))

  flag = 0

end subroutine evalf

! ------------------------------------------------------------------

subroutine evalg(n, x, g, flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)   :: n
  integer, intent(out) :: flag

  ! ARRAY ARGUMENTS
  real(kind=8), dimension(n), intent(in)  :: x
  real(kind=8), dimension(n), intent(out) :: g

  g = (/ 4 * x(1), 3.0d0 * (x(2) ** 2.0d0) * x(3), x(2) ** 3 /)

  flag = 0

end subroutine evalg

! ------------------------------------------------------------------

subroutine evalh(n, x, hnnzmax, hlin, hcol, hval, hnnz, flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)  :: n, hnnzmax
  integer, intent(out) :: hnnz, flag

  ! ARRAY ARGUMENTS
  real(kind=8), dimension(n),       intent(in)  :: x
  real(kind=8), dimension(hnnzmax), intent(out) :: hval
  integer,      dimension(hnnzmax), intent(out) :: hlin, hcol

  hlin = (/     1,               2,           3 /)
  hcol = (/     1,               2,           2 /)
  hval = (/ 4.0d0, 6 * x(2) * x(3), 3 * x(2)**2 /)

  hnnz = 3

  if(hnnz > hnnzmax) then
     flag = - 1
  else
     flag = 0
  end if

end subroutine evalh
