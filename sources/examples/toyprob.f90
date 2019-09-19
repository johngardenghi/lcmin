! --------------------------------------------------------------------
! File : toyprob.f90
! --------------------------------------------------------------------

! --------------------------------------------------------------------
! Simple codification example
! --------------------------------------------------------------------

program toyprob

  implicit none

  ! LOCAL SCALARS
  integer :: n, m, annz, hnnzmax
  integer :: allocstat, flag
  
  ! LOCAL ARRAYS
  integer,      allocatable, dimension(:) :: arow, acol
  real(kind=8), allocatable, dimension(:) :: aval, b, l, lambda, u, x

  external :: evalf, evalg, evalh

  annz = 4
  m    = 2
  n    = 3

  allocate( acol(annz), arow(annz), aval(annz), b(m), l(n), lambda(m), &
       u(n), x(n), stat = allocstat )
  if( allocstat .ne. 0 ) then
     write (*,*) "Memory allocation error, unable to continue."
     stop
  end if

  ! Initial point
  x = (/ - 12.0d0, 20.0d0, 10.0d0 /)

  ! Lower and upper bounds
  l = (/ - 10.0d0, - 10.0d0, -1.0d+20 /)
  u = (/   10.0d0,   10.0d0,  1.0d+20 /)

  ! Jacobian of the constraints
  arow = (/     1,       1,     2,     2 /)
  acol = (/     1,       2,     2,     3 /)
  aval = (/ 1.0d0, - 1.0d0, 2.0d0, 1.0d0 /)

  b = (/ 2.0d0, 6.0d0 /)

  hnnzmax = 3

  call lcmin(n, x, l, u, m, b, lambda, arow, acol, aval, annz, &
       hnnzmax, 1.0d-08, 1.0d-09, evalf, evalg, evalh, flag)

  deallocate( acol, arow, aval, b, l, lambda, u, x )

end program toyprob

! --------------------------------------------------------------------

subroutine evalf(n, x, f, flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)  :: n
  integer,      intent(out) :: flag
  real(kind=8), intent(out) :: f

  ! ARRAY ARGUMENTS
  real(kind=8), dimension(n), intent(in) :: x

  ! ==================================================================
  ! This subroutine evaluates the objective function value at a given
  ! point x
  ! ==================================================================
  !
  ! PARAMETERS
  !
  ! On entry:
  ! 
  ! n        integer, scalar
  !          number of variables.
  !
  ! x        double precision, array,
  !          the point where the objective function should be 
  !          evaluated.
  !
  ! On exit:
  !
  ! f        double precision, scalar,
  !          the value of the objective function at x.
  !
  ! flag     integer, scalar,
  !          the return of the routine (0 if sucessful)
  !
  ! ==================================================================

  f = (2 * x(1)**2) + (x(2)**3 * x(3))

  flag = 0

end subroutine evalf

! --------------------------------------------------------------------

subroutine evalg(n, x, g, flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)   :: n
  integer, intent(out) :: flag

  ! ARRAY ARGUMENTS
  real(kind=8), dimension(n), intent(in)  :: x
  real(kind=8), dimension(n), intent(out) :: g

  ! ==================================================================
  ! This subroutine evaluates the objective function gradient at a
  ! given point x
  ! ==================================================================
  !
  ! PARAMETERS
  !
  ! On entry:
  ! 
  ! n        integer, scalar
  !          number of variables.
  !
  ! x        double precision, array,
  !          the point where the objective function gradient should
  !          be evaluated.
  !
  ! On exit:
  !
  ! g        double precision, array,
  !          the gradient of the objective function at x.
  !
  ! flag     integer, scalar,
  !          the return of the routine (0 if sucessful)
  !
  ! ==================================================================

  g = (/ 4 * x(1), 3.0d0 * (x(2) ** 2.0d0) * x(3), x(2) ** 3 /)

  flag = 0

end subroutine evalg

! --------------------------------------------------------------------

subroutine evalh(n, x, hnnzmax, hrow, hcol, hval, hnnz, flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)  :: n, hnnzmax
  integer, intent(out) :: hnnz, flag

  ! ARRAY ARGUMENTS
  real(kind=8), dimension(n),       intent(in)  :: x
  real(kind=8), dimension(hnnzmax), intent(out) :: hval
  integer,      dimension(hnnzmax), intent(out) :: hrow, hcol

  ! ==================================================================
  ! This subroutine evaluates the objective function Hessian at a
  ! given point x
  ! ==================================================================
  !
  ! PARAMETERS
  !
  ! On entry:
  ! 
  ! n        integer, scalar
  !          number of variables.
  !
  ! x        double precision, array,
  !          the point where the objective function gradient should
  !          be evaluated.
  !
  ! hnnzmax  integer, scalar,
  !          the actual size of hrow, hcol and hval. The limit of
  !          nonzeros elements of the Hessian.
  !
  ! On exit:
  !
  ! hrow     integer, array,
  ! hcol     integer, array,
  ! hval     double precision, array,
  !          the Hessian of the objective function evaluated at x,
  !          in sparse coordinate format.
  !
  ! hnnz     double precision, array,
  !          the actual number of nonzeros elements in the Hessian.
  !
  ! flag     integer, scalar,
  !          the return of the routine (0 if sucessful)
  !
  ! ==================================================================

  hrow = (/     1,               2,           3 /)
  hcol = (/     1,               2,           2 /)
  hval = (/ 4.0d0, 6 * x(2) * x(3), 3 * x(2)**2 /)

  hnnz = 3

  if(hnnz > hnnzmax) then
     flag = - 1
  else
     flag = 0
  end if

end subroutine evalh
