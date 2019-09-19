module lcmin

  use constants, only: dp

  implicit none

  public :: lcmin_constraints, lcmin_load, lcmin_optimize, &
            lcmin_params, lcmin_variables

  private

  type lcmin_constraints
     ! SCALAR VARIABLES
     integer :: annz

     ! ARRAY VARIABLES
     integer,       allocatable, dimension(:) :: acol, arow
     real(kind=dp), allocatable, dimension(:) :: aval, b, l, u
  end type lcmin_constraints

  type lcmin_variables
     ! ARRAY VARIABLES
     real(kind=dp), allocatable, dimension(:) :: lambda, x, zl, zu
  end type lcmin_variables

  type lcmin_params
     ! Character parameters
     character(len=80) :: initial_file   ! Output file for the initial point
     character(len=80) :: output_file    ! LCMIN output file
     character(len=80) :: solution_file  ! Output file for the solution found by LCMIN

     ! Integer parameters
     integer :: hnnzmax      ! Maximum size of nonzeros of the Hessian
     integer :: hslunit      ! Stream for HSL output
     integer :: maxhslfix    ! Maximum iterations trying to run an HSL routine
     integer :: maxinnerit   ! Maximum inner iterations
     integer :: maxouterit   ! Maximum outer iterations
     integer :: iniunit      ! Stream for initial point output file
     integer :: outunit      ! Stream for LCMIN output (besides screen)
     integer :: solunit      ! Stream for solution output file

     ! Logical parameters
     logical :: scalea     ! Decides whether scaling the constraints
     logical :: scalef     ! Decides whether scaling the objective function
     logical :: printinit  ! Print inner iteration information
     logical :: printtl    ! Print a table line with the results from LCMIN run

     ! Double precision parameters
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
  end type lcmin_params

contains

  ! -----------------------------------------------------------------------
  ! ROUTINE OUTPUT
  !
  ! flag:
  !     0   KKT optimality conditions sufficiently satisfied
  !     1   LCMIN converged to an acceptable feasible point
  !     2   LCMIN converged to an infeasible point
  !   - 1   maximum outer iterations exceeded
  !   - 2   lack of progress by variation between iterands
  !   - 3   lack of progress by maximum inner iterations
  !   - 4   the iterates are diverging - problem may be unbounded
  !   - 5   it was not possible to find an strict feasible initial point
  !   - 6   the feasible region is void
  !   - 7   inertia correction perturbation is too large
  !   - 8   there was a fatal error on a user routine (evalf, evalg or evalh)
  !   - 9   HSL routine error
  !   - 10  memory allocation error
  ! -----------------------------------------------------------------------
  
  subroutine lcmin_load(param)

    ! TYPE ARGUMENTS
    type(lcmin_params), intent(out) :: param
    
    ! ----------------------------------------------------------------
    ! This routine defines the default controlling parameters for the
    ! solver.
    ! ----------------------------------------------------------------

    ! Character parameters
    param%initial_file  = "initial_point.txt"
    param%output_file   = "lcmin.out"
    param%solution_file = "solution.txt"

    ! Integer parameters
    param%hnnzmax    = 100000
    param%hslunit    = 50
    param%maxhslfix  = 10
    param%maxinnerit = 200
    param%maxouterit = 50
    param%iniunit    = 40
    param%outunit    = 20
    param%solunit    = 30

    ! Logical parameters
    param%scalea    = .true.
    param%scalef    = .true.
    param%printinit = .false.
    param%printtl   = .false.

    ! Double precision parameters
    param%chiini   = 1.0d-04
    param%chimin   = 1.0d-20
    param%epstol   = 1.0d-08
    param%epsfeas  = 1.0d-08
    param%gamma    = 1.0d-04
    param%gmax     = 1.0d+02
    param%imu      = 1.0d-01
    param%infty    = 1.0d+20
    param%kappad   = 1.0d-04
    param%kappaeps = 1.0d+01
    param%kappaip1 = 1.0d-02
    param%kappaip2 = 1.0d-02
    param%kappamu  = 2.0d-01
    param%kappaz   = 1.0d+10
    param%kchibar  = 1.0d+02
    param%kchimin  = 1.0d0 / 3.0d0
    param%kchimax  = 8.0d0
    param%smax     = 1.0d+02
    param%thetamu  = 1.5d+00
    param%tmin     = 9.9d-01
    param%vartol   = 1.0d-12

    ! NON EXECUTABLE STATEMENTS
100 format(/, 1X, "ERROR (lcmin_load): Unable to allocate memory.")

  end subroutine lcmin_load

  ! ------------------------------------------------------------------

  subroutine lcmin_optimize(n, m, param, cons, vars, uevalf, uevalg, &
       uevalh, flag)

    use constants
    use constraints
    use initial_point
    use lcmin_solver
    use seval
    use testing
    use user_subs
    use util
    use variables

    ! SCALAR ARGUMENTS
    integer, intent(in)  :: m, n
    integer, intent(out) :: flag

    ! PROCEDURE ARGUMENTS
    procedure(evalf) :: uevalf
    procedure(evalg) :: uevalg
    procedure(evalh) :: uevalh

    ! TYPE ARGUMENTS
    type(lcmin_constraints), intent(in)    :: cons
    type(lcmin_params),      intent(in)    :: param
    type(lcmin_variables),   intent(inout) :: vars

    ! LOCAL SCALARS
    integer :: allocstat, annz, annzrank, i, mred, nred
    real    :: start, finish

    ! LOCAL ARRAYS
    integer,                    dimension(cons%annz) :: acol, arow
    real(kind=dp),              dimension(cons%annz) :: aval
    real(kind=dp),              dimension(m)         :: b
    real(kind=dp),              dimension(n)         :: l, laggrad, u
    real(kind=dp), allocatable, dimension(:)         :: lambda, xred, zl, zu

    ! TEMP
    real(kind=dp) :: compl, compu
    real(kind=dp), dimension(n) :: tmpgrad, tmpzl, tmpzu, tmpatlambda, tmp
    real(kind=dp), dimension(m) :: tmplambda

    flag = 0

    ! ------------------------------------------------------------------
    !                           PREPROCESSING                          
    ! ------------------------------------------------------------------
    
    ! Start timing for preprocessing
    call cpu_time( start )

    ! Set up LCMIN global parameters
    chiini     = param%chiini
    chimin     = param%chimin
    epsfeas    = param%epsfeas
    epstol     = param%epstol
    gamma      = param%gamma
    gmax       = param%gmax
    hnnzmax    = param%hnnzmax
    hslunit    = param%hslunit
    imu        = param%imu
    infty      = param%infty
    ini_file   = param%initial_file
    kappad     = param%kappad
    kappaeps   = param%kappaeps
    kappaip1   = param%kappaip1
    kappaip2   = param%kappaip2
    kappamu    = param%kappamu
    kappaz     = param%kappaz
    kchibar    = param%kchibar
    kchimin    = param%kchimin
    kchimax    = param%kchimax
    maxhslfix  = param%maxhslfix
    maxinnerit = param%maxinnerit
    maxouterit = param%maxouterit
    out_file   = param%output_file
    outunit    = param%outunit
    printinit  = param%printinit
    printtl    = param%printtl
    scalea     = param%scalea
    scalef     = param%scalef
    smax       = param%smax
    sol_file   = param%solution_file
    solunit    = param%solunit
    thetamu    = param%thetamu
    tmin       = param%tmin
    vartol     = param%vartol

    ! Open output files
    open(iniunit, file=ini_file)
    open(solunit, file=sol_file)
    open(outunit, file=out_file)

    ! Writes initial information
    write (      *, 7000) n, m
    write (outunit, 7000) n, m

    ! Allocates an array to store a copy of x
    allocate ( x_user(n), stat=allocstat )
    if ( allocstat .ne. 0 ) then
       write (      *, 1000)
       write (outunit, 1000)

       flag = - 10
       go to 900
    end if

    x_user(1:n) = vars%x(1:n)

    ! Writes user initial point to the file
    write (iniunit, 2000)
    write (iniunit, 9010) (i, vars%x(i), i = 1, n)

    ! Copies the original entries of A
    annz = cons%annz

    arow(1:annz) = cons%arow(1:annz)
    acol(1:annz) = cons%acol(1:annz)
    aval(1:annz) = cons%aval(1:annz)

    b(1:m) = cons%b(1:m)

    call set_constraints( m, cons%b, annz, cons%arow, cons%acol, cons%aval )

    ! Copies the original bounds
    l(1:n) = cons%l(1:n)
    u(1:n) = cons%u(1:n)

    ! Sets the user subroutines
    subf => uevalf
    subg => uevalg
    subh => uevalh

    ! Saves the original number of variables
    nprob = n

    ! Computes the scale factor for A and scales it
    call scale_constraints( m, b, annz, arow, acol, aval, flag )
    if ( flag .ne. 0 ) go to 900

    ! call print_octave( m, n, la, ar, ac, av, 123, 'A', flag)
    ! if ( flag .ne. 0 ) go to 900

    ! call print_octave( m, bc, 123, 'b', flag)
    ! if ( flag .ne. 0 ) go to 900

    ! Check constraints for implicit equalities
    if ( m .gt. 0 ) then
       call check_implicit_equalities( n, l, u, m, b, annz, arow, acol, &
            aval, flag )
       if ( flag .ne. 0 ) go to 900
    end if

    ! Check bound constraints for fixed variables and reorder the
    ! jacobian of the constraints to eliminate these if there is any.
    call check_constraints_fxv( m, b, n, x_user, l, u, annz, arow, &
         acol, aval, nred, flag )
    if ( flag .ne. 0 ) go to 900

    ! Allocates LCMIN variables
    allocate( xred(nred), stat=allocstat )
    if ( allocstat .ne. 0 ) then
       write (      *, 1000)
       write (outunit, 1000)

       flag = - 10
       go to 900
    end if

    ! Check bound constraints for free variables. lind and uind hold the
    ! indexes of non-free variables from below and non-free variables
    ! from above, respectively.
    call check_constraints_frv( nred, l(nfind), u(nfind), flag )
    if ( flag .ne. 0 ) go to 900

    ! Checks the jacobian of the constraints for linear dependent lines
    call check_constraints_rank( nred, m, b, annz, arow, acol, &
         aval, annzrank, mred, flag )
    if ( flag .ne. 0 ) go to 900

    ! Allocates LCMIN variables
    allocate( lambda(mred), zl(nl), zu(nu), stat=allocstat )
    if ( allocstat .ne. 0 ) then
       write (      *, 1000)
       write (outunit, 1000)

       flag = - 10
       go to 900
    end if

    ! Computes the initial point
    if ( m .gt. 0 .or.  nl .gt. 0 .or. nu .gt. 0 ) then

       call ipoint( nred, x_user(nfind), l(nfind), u(nfind), mred, b, &
            annz, arow, acol, aval, xred, flag )
       if ( flag .ne. 0 ) go to 900

       x_user(nfind(1:nred)) = xred

    else
       xred(1:nred) = x_user(nfind(1:nred))
    end if

    ! Writes LCMIN initial point to the file
    write (iniunit, 2010)
    write (iniunit, 9010) (i, x_user(i), i = 1, n)

900 continue

    ! End timing for preprocessing
    call cpu_time( finish )

    pretime = finish - start

    if ( flag .ne. 0 ) then
    
       ! Write null tabline
       if ( printtl ) then
          open ( 350, file="lcmin-tabline.out" )

          ! Writes LCMIN flag
          write (350, fmt='(I2)') flag

          write (350, 650) n, m, bignum, bignum, bignum, bignum, &
                           bignum, 0, 0, 0, 0, 0, 0, 0, pretime, 0.0

          close (350)
       end if

    else

       ! ------------------------------------------------------------------
       !                           OPTIMIZATION                           
       ! ------------------------------------------------------------------

       ! Writes initial information
       write(      *, 7010) nred, mred
       write(outunit, 7010) nred, mred

       ! Call solver
       call lcmin_solve( nred, xred, l(nfind), u(nfind), mred, b, &
            annz, arow, acol, aval, lambda, zl, zu, flag )

       ! COMPUTES THE USER PRIMAL AND DUAL VARIABLES
       vars%x(1:n) = x_user(1:n)

       ! if ( allocated( cind ) ) then
       !    vars%lambda(1:m) = 0.0d0

       !    vars%lambda(cind(1:mred)) = lambda(1:mred) / sf
       !    if ( scalea ) then
       !       vars%lambda(cind(1:mred)) = sa(cind(1:mred)) * vars%lambda(cind(1:mred))
       !       ! vars%lambda(1:m) = sa(1:m) * vars%lambda(1:m)
       !    end if
       ! else
       !    vars%lambda(1:m) = lambda(1:m) / sf
       !    if ( scalea ) then 
       !      vars%lambda(1:m) = sa(1:m) * vars%lambda(1:m)
       !    end if
       ! end if

       ! vars%zl(1:n) = 0.0d0
       ! vars%zl(nfind(lind(1:nl))) = zl(1:nl) / sf

       ! vars%zu(1:n) = 0.0d0
       ! vars%zu(nfind(uind(1:nu))) = zu(1:nu) / sf

       ! Computes the gradient of the Lagrangian to adjust zl and zu
       ! call subg( n, vars%x, laggrad, flag )
       ! gevalcnt = gevalcnt + 1

       ! tmp = smvp( n, cons%annz, cons%acol, cons%arow, cons%aval, vars%lambda )
       ! tmpgrad = laggrad
       ! tmpatlambda = 0.0d0

       ! ! open (555, file='test.txt')
       ! ! read (555, fmt='(1P, D24.16)') (tmpgrad(nfind(i)), i = 1, nred)
       ! ! read (555, fmt='(1P, D24.16)') (tmpatlambda(nfind(i)), i = 1, nred)
       ! ! close (555)

       ! ! print '(/)'
       ! ! print *, maxval( abs( tmpgrad - laggrad ) )
       ! ! print *, maxval( abs( tmplambda(1:mred) - vars%lambda(cind(1:mred)) ) )
       ! ! print *, maxval( abs( tmplambda(1:mred) - vars%lambda(cind(1:mred)) ) )
       ! ! print *, maxval( abs( tmpatlambda(nfind) - tmp(nfind) ) )
       ! ! print '(/)'

       ! laggrad(1:n) = laggrad(1:n) &
       !      + smvp( n, cons%annz, cons%acol, cons%arow, cons%aval, vars%lambda )

       ! ! do i = 1, n
       ! !    if ( find(i) .eq. 0 .and. abs( cons%l(i) - cons%u(i) ) .gt. macheps ) print *, i, abs( laggrad(i) )
       ! ! end do

       ! do i = 1, n
       !    if ( find(i) .eq. 0 .and. abs( cons%l(i) - cons%u(i) ) .le. macheps ) then
       !    ! if ( find(i) .eq. 0 ) then
       !    ! if ( find(i) .eq. 0 .and. vars%x(i) .gt. cons%l(i) .and. vars%x(i) .lt. cons%u(i) ) then
       !       if ( laggrad(i) .gt. 0.0d0 ) then
       !          vars%zl(i) = laggrad(i)
       !       else
       !          vars%zu(i) = - laggrad(i)
       !       end if
       !    end if
       ! end do

       ! -------------> TESTING
       ! laggrad(1:n) = laggrad(1:n) - vars%zl(1:n) + vars%zu(1:n)

       ! print '(/, A, 1P, D24.16)', "Lagrangian gradient norm ......: ", maxval( abs( laggrad(1:n) ) )
       ! print '(A, 1P, D24.16)', "Feasibility ...................: ", &
       !      maxval( abs( smvp( m, cons%annz, cons%arow, cons%acol, cons%aval, vars%x ) - cons%b ) )

       ! compl = 0.0d0
       ! compu = 0.0d0
       ! do i = 1, n
       !    if ( cons%l(i) .gt. - infty ) then
       !       compl = max( compl, ( vars%x(i) - cons%l(i) ) * vars%zl(i) )
       !    end if

       !    if ( cons%u(i) .lt.   infty ) then
       !       compu = max( compu, ( cons%u(i) - vars%x(i) ) * vars%zu(i) )
       !    end if
       ! end do

       ! print '(A, 1P, D24.16)', "Lower bound complementarity ...: ", compl
       ! print '(A, 1P, D24.16)', "Upper bound complementarity ...: ", compu
       ! --------------------

    end if

    ! ------------------------------------------------------------------
    !                     WRITES THE FINAL RESULTS                     
    ! ------------------------------------------------------------------

    ! Write the LCMIN exit flag
    write (      *, 7020)
    write (outunit, 7020)

    select case ( flag )

    case (2)
       write (      *, 780)
       write (outunit, 780)

    case (1)
       write (      *, 790)
       write (outunit, 790)

    case (0)
       write (      *, 800)
       write (outunit, 800)

    case (-1)
       write (      *, 810)
       write (outunit, 810)

    case (-2)
       write (      *, 820)
       write (outunit, 820)

    case (-3)
       write (      *, 830)
       write (outunit, 830)

    case (-4)
       write (      *, 840)
       write (outunit, 840)

    case (-5)
       write (      *, 850)
       write (outunit, 850)

    case (-6)
       write (      *, 860)
       write (outunit, 860)

    case (-7)
       write (      *, 870)
       write (outunit, 870)

    case (-8)
       write (      *, 880)
       write (outunit, 880)

    case (-9)
       write (      *, 890)
       write (outunit, 890)

    case (-10)
       write (      *, 8100)
       write (outunit, 8100)

    end select

    ! Writes the statistics
    if( inneraccitcnt .gt. 0 ) then
       write (      *, 9000) pretime, opttime, inneraccitcnt,        &
                             outeritcnt, fevalcnt, gevalcnt,         &
                             hevalcnt, bfevalcnt, bgevalcnt,         &
                             real(inneraccitcnt) / real(outeritcnt), &
                             real(fevalcnt) / real(inneraccitcnt),   &
                             real(fevalcnt) / real(outeritcnt),      &
                             real(gevalcnt) / real(inneraccitcnt),   &
                             real(gevalcnt) / real(outeritcnt),      &
                             real(hevalcnt) / real(inneraccitcnt),   &
                             real(hevalcnt) / real(outeritcnt),      &
                             real(bfevalcnt) / real(inneraccitcnt),  &
                             real(bfevalcnt) / real(outeritcnt),     &
                             real(bgevalcnt) / real(inneraccitcnt),  &
                             real(bgevalcnt) / real(outeritcnt)

       write (outunit, 9000) pretime, opttime, inneraccitcnt,        &
                             outeritcnt, fevalcnt, gevalcnt,         &
                             hevalcnt, bfevalcnt, bgevalcnt,         &
                             real(inneraccitcnt) / real(outeritcnt), &
                             real(fevalcnt) / real(inneraccitcnt),   &
                             real(fevalcnt) / real(outeritcnt),      &
                             real(gevalcnt) / real(inneraccitcnt),   &
                             real(gevalcnt) / real(outeritcnt),      &
                             real(hevalcnt) / real(inneraccitcnt),   &
                             real(hevalcnt) / real(outeritcnt),      &
                             real(bfevalcnt) / real(inneraccitcnt),  &
                             real(bfevalcnt) / real(outeritcnt),     &
                             real(bgevalcnt) / real(inneraccitcnt),  &
                             real(bgevalcnt) / real(outeritcnt)

    else
       if ( outeritcnt .gt. 0 ) then
          write (      *, 9000) pretime, opttime, inneraccitcnt, outeritcnt, &
                                fevalcnt, gevalcnt, hevalcnt, bfevalcnt,     &
                                bgevalcnt,                                   &
                                real(inneraccitcnt) / real(outeritcnt),      &
                                0.0_dp, real(fevalcnt) / real(outeritcnt),   &
                                0.0_dp, real(gevalcnt) / real(outeritcnt),   &
                                0.0_dp, real(hevalcnt) / real(outeritcnt),   &
                                0.0_dp, real(bfevalcnt) / real(outeritcnt),  &
                                0.0_dp, real(bgevalcnt) / real(outeritcnt)

          write (outunit, 9000) pretime, opttime, inneraccitcnt, outeritcnt, &
                                fevalcnt, gevalcnt, hevalcnt, bfevalcnt,     &
                                bgevalcnt,                                   &
                                real(inneraccitcnt) / real(outeritcnt),      &
                                0.0_dp, real(fevalcnt) / real(outeritcnt),   &
                                0.0_dp, real(gevalcnt) / real(outeritcnt),   &
                                0.0_dp, real(hevalcnt) / real(outeritcnt),   &
                                0.0_dp, real(bfevalcnt) / real(outeritcnt),  &
                                0.0_dp, real(bgevalcnt) / real(outeritcnt)
       else
          write (      *, 9000) pretime, opttime, inneraccitcnt, outeritcnt, &
                                fevalcnt, gevalcnt, hevalcnt, bfevalcnt,     &
                                bgevalcnt, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,   &
                                0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,      &
                                0.0_dp, 0.0_dp

          write (outunit, 9000) pretime, opttime, inneraccitcnt, outeritcnt, &
                                fevalcnt, gevalcnt, hevalcnt, bfevalcnt,     &
                                bgevalcnt, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,   &
                                0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,      &
                                0.0_dp, 0.0_dp
       end if
    end if
    
    ! Writes the solution to the file
    ! TODO: fix lambda to deal with rank-deficient A
    write (solunit, 9020)
    write (solunit, 9010) (i, vars%x(i), i = 1, n)
    write (solunit, 9030)
    write (solunit, 9010) (i, vars%lambda(i), i = 1, m)

    close( outunit )
    close( solunit )

    deallocate( find, lambda, lind, lnu, nfind, sa, uind, unl, xred, &
         x_user, zl, zu, stat=allocstat )

    if ( allocated( cind ) ) deallocate( cind, stat=allocstat )

    nullify( subf, subg, subh )

    ! NONEXECUTABLE STATEMENTS
650 format(2(I6, 3X), 1P, 2(D24.16, 3X), 3(D11.4, 3X), 7(I7, 3X), 0P, 2(F9.2, 3X), "&")

780 format(/, 5X, "Convergence to an infeasible point.")

790 format(/, 5X, "Acceptable feasible optimal solution found.")

800 format(/, 5X, "Optimal solution found.")

810 format(/, 5X, "Maximum outer iterations exceeded.")

820 format(/, 5X, "Lack of progress.")

830 format(/, 5X, "Maximum inner iterations exceeded.")

840 format(/, 5X, "Iterates diverging: the problem may be unbounded.")

850 format(/, 5X, "Unable to find an interior initial point.", &
           /, 5X, "The interior of the feasible region may be void.")

860 format(/, 5X, "The problem may be infeasible.")

870 format(/, 5X, "Inertia correction failed.")

880 format(/, 5X, "There was a fatal error on a evaluation subroutine.")

890 format(/, 5X, "There was a fatal error on a HSL subroutine.")

8100 format(/, 5X, "Memory allocation error.")

1000 format(/, 1X, "ERROR (lcmin): Unable to allocate memory.")

1010 format(/, 1X, "ERROR (lcmin): Unable to deallocate memory.")

2000 format(/, 1X, "USER INITIAL POINT:", //, 3X, "INDEX", 4X, "X(INDEX)")

2010 format(/, 1X, "LCMIN INITIAL POINT:", //, 3X, "INDEX", 4X, "X(INDEX)")

7000 format(/, 1X, "Number of variables            : ", I7, &
            /, 1X, "Number of equality constraints : ", I7)

7010 format(/, 1X, "Starting LCMIN.", &
            /, 1X, "Number of variables   : ", I7, &
            /, 1X, "Number of constraints : ", I7)

7020 format(/, 1X, "Flag of LCMIN:")

9000 format(/, 1X, "====================================================", &
            /, 1X, "                     STATISTICS                     ", &
            /, 1X, "====================================================", &
            /, 1X, "GENERAL:                                            ", &
           //, 1X, "          Preprocessing time (in seconds): ", F9.2,    &
            /, 1X, "           Optimization time (in seconds): ", F9.2,    &
           //, 1X, "                     Inner Iterations (IIT): ", I7,    &
            /, 1X, "                     Outer Iterations (OIT): ", I7,    &
            /, 1X, "                 Function Evaluations (FEV): ", I7,    &
            /, 1X, "        Function Gradient Evaluations (GEV): ", I7,    &
            /, 1X, "         Function Hessian Evaluations (HEV): ", I7,    &
            /, 1X, "         Barrier Function Evaluations (BFE): ", I7,    &
            /, 1X, "Barrier Function Gradient Evaluations (BGE): ", I7,    &
            /, 1X, "====================================================", &
            /, 1X, "AVERAGE:                                            ", &
           //, 1X, "Inner Iterations per Outer Iterations  : ", F11.4,     &
            /, 1X, "FEV per IIT: ", F11.4, 4X, "FEV per OIT: ", F11.4,     &
            /, 1X, "GEV per IIT: ", F11.4, 4X, "GEV per OIT: ", F11.4,     &
            /, 1X, "HEV per IIT: ", F11.4, 4X, "HEV per OIT: ", F11.4,     &
            /, 1X, "BFE per IIT: ", F11.4, 4X, "BFE per OIT: ", F11.4,     &
            /, 1X, "BGE per IIT: ", F11.4, 4X, "BGE per OIT: ", F11.4,     &
            /, 1X, "====================================================", /)

9010 format(1X, I7, 2X, 1P, D24.16)

9020 format(/, 1X, "FINAL POINT:", //, 3X, "INDEX", 4X, "X(INDEX)")

9030 format(/, 1X, "LAGRANGE MULTIPLIERS:", //, 3X, "INDEX", 4X, "LAMBDA(INDEX)")

  end subroutine lcmin_optimize

end module lcmin
