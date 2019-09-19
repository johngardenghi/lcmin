module cutest_subs

  implicit none

  save

  integer :: n_

contains

  subroutine cutest_cevalf(n, x, f, flag)

    ! SCALAR ARGUMENTS
    integer,      intent(in)  :: n
    integer,      intent(out) :: flag
    real(kind=8), intent(out) :: f

    ! ARRAY ARGUMENTS
    real(kind=8), dimension(n), intent(in) :: x

    ! LOCAL SCALARS
    integer :: dum, stat

    call cutest_cofg( stat, n_, x(1:n_), f, dum, .false. )

    flag = - stat

  end subroutine cutest_cevalf

  ! ------------------------------------------------------------------

  subroutine cutest_cevalg(n, x, g, flag)

    ! SCALAR ARGUMENTS
    integer, intent(in)  :: n
    integer, intent(out) :: flag

    ! ARRAY ARGUMENTS
    real(kind=8), dimension(n), intent(in)  :: x
    real(kind=8), dimension(n), intent(out) :: g

    ! LOCAL SCALARS
    integer      :: stat
    real(kind=8) :: f

    call cutest_cofg( stat, n_, x(1:n_), f, g(1:n_), .true. )

    g(n_+1:n) = 0.0d0

    flag = - stat

  end subroutine cutest_cevalg

  ! ------------------------------------------------------------------

  subroutine cutest_cevalh(n, x, hnnzmax, hrow, hcol, hval, hnnz, flag)

    ! SCALAR ARGUMENTS
    integer, intent(in)  :: n, hnnzmax
    integer, intent(out) :: hnnz, flag

    ! ARRAY ARGUMENTS
    real(kind=8), dimension(n),       intent(in)  :: x
    real(kind=8), dimension(hnnzmax), intent(out) :: hval
    integer,      dimension(hnnzmax), intent(out) :: hrow, hcol

    ! LOCAL SCALARS
    integer :: stat

    call cutest_cish( stat, n_, x(1:n_), 0, hnnz, hnnzmax, hval, hrow, hcol )

    if(hnnz .gt. hnnzmax) then
       flag = - 1
    else
       flag = - stat
    end if

  end subroutine cutest_cevalh

  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------

  subroutine cutest_uevalf(n, x, f, flag)

    ! SCALAR ARGUMENTS
    integer,      intent(in)  :: n
    integer,      intent(out) :: flag
    real(kind=8), intent(out) :: f

    ! ARRAY ARGUMENTS
    real(kind=8), dimension(n), intent(in) :: x

    ! LOCAL SCALARS
    integer :: dum, stat

    call cutest_ufn( stat, n, x(1:n), f )

    flag = - stat

  end subroutine cutest_uevalf

  ! ------------------------------------------------------------------

  subroutine cutest_uevalg(n, x, g, flag)

    ! SCALAR ARGUMENTS
    integer, intent(in)  :: n
    integer, intent(out) :: flag

    ! ARRAY ARGUMENTS
    real(kind=8), dimension(n), intent(in)  :: x
    real(kind=8), dimension(n), intent(out) :: g

    ! LOCAL SCALARS
    integer      :: stat
    real(kind=8) :: f

    call cutest_ugr( stat, n, x(1:n), g(1:n) )

    flag = - stat

  end subroutine cutest_uevalg

  ! ------------------------------------------------------------------

  subroutine cutest_uevalh(n, x, hnnzmax, hrow, hcol, hval, hnnz, flag)

    ! SCALAR ARGUMENTS
    integer, intent(in)  :: n, hnnzmax
    integer, intent(out) :: hnnz, flag

    ! ARRAY ARGUMENTS
    real(kind=8), dimension(n),       intent(in)  :: x
    real(kind=8), dimension(hnnzmax), intent(out) :: hval
    integer,      dimension(hnnzmax), intent(out) :: hrow, hcol

    ! LOCAL SCALARS
    integer :: stat

    call cutest_ush( stat, n, x(1:n), hnnz, hnnzmax, hval, hrow, hcol )

    if(hnnz .gt. hnnzmax) then
       flag = - 1
    else
       flag = - stat
    end if

  end subroutine cutest_uevalh

end module cutest_subs

! ====================================================================
! ====================================================================

program lcminma_cutest

  use cutest_subs
  use lcmin

  implicit none

  ! LOCAL PARAMETERS
  integer, parameter :: dp = kind(0.0d0)

  ! LOCAL TYPE
  type(lcmin_constraints) :: constraints
  type(lcmin_params)      :: param
  type(lcmin_variables)   :: variables

  ! LOCAL SCALARS
  character(len=10) :: pname
  integer           :: allocstat, annz, hnnzmax, i, inequatn, info, m, &
                       n, stat
  
  ! LOCAL ARRAYS
  character(len=10), allocatable, dimension(:) :: cnames, vnames
  logical,           allocatable, dimension(:) :: equatn, linear
  real(kind=dp),     allocatable, dimension(:) :: cl, cu,   &
                                                  tmp, &
                                                  zeros

  ! Initialize CUTEst processing
  open(10, file="OUTSDIF.d", form="formatted", status="old")

  rewind 10

  ! Retrieves the number of variables and constraints
  call cutest_cdimen( stat, 10, n, m )

  ! If the problem is constrained
  if ( m .gt. 0 ) then

     allocate( constraints%b(m), constraints%l(n), constraints%u(n), &
               cl(m), cnames(m), cu(m), equatn(m), linear(m),        &
               variables%lambda(m), variables%x(n), vnames(n),       &
               zeros(n), stat=allocstat )
     if ( allocstat .ne. 0 ) then
        write (*, *) "It was not possible to allocate memory."
        stop
     end if

     ! Retrieves the initial point, the lower and upper bounds, the
     ! initial estimation for the lagrange multipliers and the
     ! constraints information
     call cutest_csetup( stat, 10, 20, 11, n, m, variables%x,     &
          constraints%l, constraints%u, variables%lambda, cl, cu, &
          equatn, linear, 1, 0, 0 )

     close(10)

     ! Retrieves the number of nonzero entries on the Lagrangian Hessian
     ! (= Hessian of the objective function, since we are considering
     ! only linear constraints)
     call cutest_cdimsh( stat, hnnzmax )

     ! Retrieves the number of nonzero entries on the gradient of the
     ! objective function and the gradient of the constraints
     call cutest_cdimsj( stat, annz )

     ! Retrieves the problem name
     call cutest_cnames( stat, n, m, pname, vnames, cnames )

     ! Checks if is there any nonlinear constraints. If there is, then
     ! LCMIN is not applicable.
     if( .not. all( linear ) ) then
        write (*, *) "The CUTEst problem " // pname // &
             " is not a linearly constrained problem."
        stop
     else
        write (*, 10) pname
     end if

     ! Allocate memory for the jacobian of the constraints
     allocate( constraints%arow(annz+m), constraints%acol(annz+m), &
          constraints%aval(annz+m), stat=allocstat )
     if( allocstat .ne. 0 ) then
        write (*, *) "It was not possible to allocate memory."
        stop
     end if

     ! Retrieves the jacobian of the constraints
     zeros(1:n) = 0.0d0

     call cutest_ccfsg( stat, n, m, zeros, constraints%b, annz, annz, &
          constraints%aval, constraints%acol, constraints%arow, .true. )

     ! Stores the right hand side of the constraints
     constraints%b(1:m) = - constraints%b(1:m)

     ! Stores original number of variables
     n_ = n

     ! Counts the number of inequalities constraints
     inequatn = count( equatn(1:m) .eqv. .false. )

     ! Verifies if slack variables are needed
     if( inequatn .gt. 0 ) then
        ! Reallocates memory for l, u and x
        call move_alloc( constraints%l, tmp )

        allocate( constraints%l(n+inequatn), stat=allocstat )
        if( allocstat .ne. 0 ) then
           write (*, *) "It was not possible to allocate memory."
           stop
        end if

        constraints%l(1:n) = tmp
        deallocate( tmp )

        call move_alloc( constraints%u, tmp )

        allocate( constraints%u(n+inequatn), stat = allocstat )
        if( allocstat .ne. 0 ) then
           write (*, *) "It was not possible to allocate memory."
           stop
        end if

        constraints%u(1:n) = tmp
        deallocate( tmp )

        call move_alloc( variables%x, tmp )

        allocate( variables%x(n+inequatn), stat = allocstat )
        if( allocstat .ne. 0 ) then
           write (*, *) "It was not possible to allocate memory."
           stop
        end if

        variables%x(1:n) = tmp
        deallocate( tmp )

        do i = 1, m
           if( equatn(i) .eqv. .false. ) then
              n    = n    + 1
              annz = annz + 1

              ! Adds the slack variable to the ith constraint
              ! Replace 
              !          cl <= a_i(x) <= cu
              ! by
              !          a_i(x) - s = b and cl <= s <= cu
              constraints%arow(annz) = i
              constraints%acol(annz) = n
              constraints%aval(annz) = - 1.0d0

              constraints%l(n) = cl(i)
              constraints%u(n) = cu(i)
           end if
        end do

        write (*, 20) n_, n - n_
     end if

     constraints%annz = annz

     allocate(variables%zl(n+inequatn), variables%zu(n+inequatn), &
          stat=allocstat)
     if( allocstat .ne. 0 ) then
        write (*, *) "It was not possible to allocate memory."
        stop
     end if

     deallocate( cl, cnames, cu, equatn, linear, vnames, zeros )

     write (*, 30)

     ! Initialize default LCMIN parameters
     call lcmin_load( param )

     ! Set custom LCMIN parameters values
     param%hnnzmax   = hnnzmax
     param%printinit = .true.
     param%printtl   = .true.

     ! Call solver
     call lcmin_optimize( n, m, param, constraints, variables, &
          cutest_cevalf, cutest_cevalg, cutest_cevalh, info )

  ! The problem is unconstrained or box constrained
  else

     constraints%annz = 0

     allocate( constraints%arow(0), constraints%acol(0), &
               constraints%aval(0), constraints%b(0),    &
               constraints%l(n), constraints%u(n),       &
               variables%lambda(0), variables%x(n),      &
               variables%zl(n), variables%zu(n),         &
               stat=allocstat )

     if ( allocstat .ne. 0 ) then
        write (*, *) "It was not possible to allocate memory."
        stop
     end if

     ! Retrieves the initial point, the lower and upper bounds, the
     ! initial estimation for the lagrange multipliers and the
     ! constraints information
     call cutest_usetup( stat, 10, 20, 11, n, variables%x, &
          constraints%l, constraints%u )

     close(10)

     ! Stores original number of variables
     n_ = n

     ! Retrieves the number of nonzero entries on the Hessian of the
     ! objective function
     call cutest_udimsh( stat, hnnzmax )

     ! Retrieves the problem name
     call cutest_probname( stat, pname )

     write (*, 10) pname

     write (*, 30)

     ! Initialize default LCMIN parameters
     call lcmin_load( param )

     ! Set custom LCMIN parameters values
     param%hnnzmax    = hnnzmax
     param%printinit  = .true.
     param%printtl    = .true.
     param%maxinnerit = 50000

     ! Call solver
     call lcmin_optimize( n, m, param, constraints, variables, &
          cutest_uevalf, cutest_uevalg, cutest_uevalh, info )

  end if

  deallocate( constraints%acol, constraints%arow, constraints%aval,   &
       constraints%b, constraints%l, constraints%u, variables%lambda, &
       variables%x, variables%zl,variables%zu, stat=allocstat )

  ! NONEXECUTABLE STATEMENTS
10 format(/, 1X, "-----------------------------------------", &
          /, 1X, "Solving the CUTEst problem ", A10)

20 format(/, 1X, "This problem has inequality constraints.", &
          /, 1X, "Slacks variables will be used.",           &
         //, 1X, "Original number of variables : ", I7,      &
          /, 1X, "Number of slacks to be used  : ", I7)

30 format(   1X, "-----------------------------------------")

end program lcminma_cutest
