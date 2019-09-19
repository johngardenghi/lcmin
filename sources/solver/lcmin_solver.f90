module lcmin_solver

  implicit none

  private

  public :: lcmin_solve

contains

  subroutine lcmin_solve(n, x, l, u, m, b, annz, arow, acol, aval, &
       lambda, zl, zu, flag)

    use constants
    use constraints
    use initial_point
    use seval
    use util
    use variables

    ! SCALAR ARGUMENTS
    integer, intent(in)  :: annz, m, n
    integer, intent(out) :: flag

    ! ARRAY ARGUMENTS
    integer,       dimension(annz), intent(in)    :: acol, arow
    real(kind=dp), dimension(annz), intent(in)    :: aval
    real(kind=dp), dimension(m),    intent(in)    :: b
    real(kind=dp), dimension(n),    intent(in)    :: l, u
    real(kind=dp), dimension(n),    intent(inout) :: x
    real(kind=dp), dimension(m),    intent(out)   :: lambda
    real(kind=dp), dimension(nl),   intent(out)   :: zl
    real(kind=dp), dimension(nu),   intent(out)   :: zu

    ! LOCAL SCALARS
    integer       :: cnt
    real          :: curtime

    real(kind=dp) :: complnorm, compunorm, func, lambdanorm,          &
                     laggradnorm, minsa, mu, outerr, scafeas, sd, sl, &
                     su, unscafeas, zlnorm, zunorm

    ! LOCAL ARRAYS
    real(kind=dp), dimension(n) :: grad, laggrad

    ! Start timing for optimization
    call cpu_time( opttime )

    ! COUNTERS INITIALIZATION
    inneraccitcnt = 0
    inneritcnt    = 0
    outeritcnt    = 0
    fevalcnt      = 0
    gevalcnt      = 0
    hevalcnt      = 0
    bfevalcnt     = 0
    bgevalcnt     = 0
    smallsd       = 0

    ! PARAMETERS INITIALIZATION
    if ( nl .gt. 0 .or. nu .gt. 0 ) then
       mu = imu
    else
       mu = epstol / kappaeps
    end if

    sd = 1.0d0
    sl = 1.0d0
    su = 1.0d0

    zl(1:nl) = 1.0d0
    zu(1:nu) = 1.0d0

    ! Computes the scale factor for f
    call sevalg( .false., 1.0d0, grad, flag )
    if(flag .lt. 0) return

    if ( scalef ) then
       sf = max( macheps12, min( 1.0d0, gmax / infnorm(grad(1:n)) ) )

       if ( sf .ne. 1.0d0 ) grad(1:n) = sf * grad(1:n)
    else
       sf = 1.0d0
    end if

    ! Computes the function value at the current point
    call sevalf( scalef, sf, func, flag )
    if ( flag .ne. 0 ) return

    ! Estimates initial lagrange multipliers
    if( m .gt. 0 ) then
       call imult( m, n, grad, zl, zu, annz, arow, acol, aval, lambda )
    end if

    ! Writes the initial lagrange multipliers to the file
    write(iniunit, 200)
    write(iniunit, 100) ( cnt, lambda(cnt), cnt = 1, m )

    close(iniunit)

    ! Writes scaling information
    if ( m .gt. 0 .and. scalea ) then
       minsa = minval(sa, 1)
    else
       minsa = 1.0d0
    end if

    write (      *, 300) sf, minsa
    write (outunit, 300) sf, minsa

    ! Norm of the lagrangian gradient at the initial point
    laggrad(1:n) = grad(1:n) + smvp( n, annz, acol, arow, aval, lambda(1:m) )
    laggrad(lind(1:nl)) = laggrad(lind(1:nl)) - zl(1:nl)
    laggrad(uind(1:nu)) = laggrad(uind(1:nu)) + zu(1:nu)
    laggradnorm = infnorm(laggrad(1:n))

    ! Required values to compute the optimality error at the initial point
    lambdanorm = norm1( lambda(1:m) )
    zlnorm     = norm1( zl(1:nl) )
    zunorm     = norm1( zu(1:nu) )

    if( (m + nl + nu) .gt. 0 ) then
       sd = max(smax, (lambdanorm + zlnorm + zunorm) / (m + nl + nu)) / smax
    end if

    if( nl .gt. 0 ) then
       sl = max(smax, zlnorm / nl) / smax
    end if

    if( nu .gt. 0 ) then
       su = max(smax, zunorm / nu) / smax
    end if

    ! This flag means first iteration
    flag = 50

    ! ------------------------------------------------------------------
    !                              MAIN LOOP
    ! ------------------------------------------------------------------
    do

       ! =============================================================
       !                     VERIFY STOP CRITERIA
       ! =============================================================
       
       ! Complementarity norm
       complnorm = infnorm( ( x(lind(1:nl)) - l(lind(1:nl)) ) * zl(1:nl) )
       compunorm = infnorm( ( u(uind(1:nu)) - x(uind(1:nu)) ) * zu(1:nu) )

       ! Optimality error for the main problem (mu = 0)
       outerr = max( laggradnorm / sd, complnorm / sl, compunorm / su )

       ! Optimality
       if ( outerr .le. epstol ) then
          call print_iteration("C")
          flag = 0
          exit
       end if

       ! Number of iterations
       if ( outeritcnt .gt. maxouterit ) then
          call print_iteration("M")
          flag = - 1
          exit
       end if

       ! =============================================================

       ! Print inner iteration information
       select case ( flag )

       case (0)
          call print_iteration("C")

       case (-3)
          call print_iteration("M")
          if( abs( mu - 1.0d-01 * epstol ) .le. macheps ) exit

       case (-4)
          call print_iteration("D")
          exit

       case (50)
          call print_iteration(" ")

       case default
          exit

       end select

       ! =============================================================
       !                         INNER SOLVER
       ! =============================================================

       call lcmin_barrier( mu, n, x, l, u, zl, zu, m, b, lambda, annz, &
            arow, acol, aval, func, grad, laggradnorm, sd, sl, su, flag )

       ! =============================================================

       ! Updates barrier parameter
       mu = max( 1.0d-01 * epstol, min( kappamu * mu, mu ** thetamu ) )

       ! Increments iteration counter
       outeritcnt = outeritcnt + 1

    end do

    ! Evaluates the scaled constraints feasibility of the final point
    scafeas = eval_feas( x_user, .true. )

    ! Evaluates the unscaled constraints feasibility of the final point
    unscafeas = eval_feas( x_user, .false. )

    ! End timing for optimization
    call cpu_time( curtime )
    opttime = curtime - opttime

    ! Writes the final KKT results
    write (      *, 500) func, func / sf, laggradnorm, laggradnorm / sf, &
                         scafeas, unscafeas, complnorm, complnorm / sf,  &
                         compunorm, compunorm / sf

    write (outunit, 500) func, func / sf, laggradnorm, laggradnorm / sf, &
                         scafeas, unscafeas, complnorm, complnorm / sf,  &
                         compunorm, compunorm / sf

    ! Checks final feasibility and reports
    if ( flag .eq. 0 .and. unscafeas .gt. epsfeas ) then
       if ( unscafeas .le. sqrt(epsfeas) ) then
          ! LCMIN converged to an acceptable feasible point
          flag = 1
       else
          flag = 2
       end if
    end if

    ! Print table line
    if ( printtl ) then
       open ( 350, file="lcmin-tabline.out" )

       ! Writes LCMIN flag
       write (350, fmt='(I2)') flag
       write (350, 600) n, m, func, func/sf, scafeas, unscafeas, &
                        outerr, fevalcnt, gevalcnt, hevalcnt,    &
                        bfevalcnt, bgevalcnt, inneraccitcnt,     &
                        outeritcnt, pretime, opttime
       close ( 350 )
    end if

    ! NON EXECUTABLE STATEMENTS
100 format(1X, I7, 2X, 1P, D24.16)

200 format(/, 1X, "ESTIMATE OF LAGRANGE MULTIPLIERS:", //, 3X, "INDEX", &
               4X, "LAMBDA(INDEX)")

300 format(/, 1X, "Scale factor for the objective function   : ", 1P, D7.1, &
           /, 1X, "Smallest scale factor for the constraints : ", 1P, D7.1)

500 format(   1X, "==========================================================", &
          //, 1X, "FINAL RESULTS:",                                           &
           /, 1X, "                  ", 9X, "scaled", 9X, 2X, 8X, "unscaled", &
           /, 1X, " Objective Funct: ", 1P, D24.16, 2X, 1P, D24.16,           &
           /, 1X, " Lagrangian Grad: ", 1P, D24.16, 2X, 1P, D24.16,           &
           /, 1X, "Final Point Feas: ", 1P, D24.16, 2X, 1P, D24.16,           &
           /, 1X, "Lower Bound Comp: ", 1P, D24.16, 2X, 1P, D24.16,           &
           /, 1X, "Upper Bound Comp: ", 1P, D24.16, 2X, 1P, D24.16)

600 format(2(I6, 3X), 1P, 2(D24.16, 3X), 3(D11.4, 3X), 7(I7, 3X), 0P, &
           2(F9.2, 3X))

  contains

    subroutine print_iteration(reason)

      ! SCALAR ARGUMENT
      character(len=1), intent(in) :: reason

      if ( outeritcnt .gt. 0 .and. printinit ) then
         write (      *, '(/)')
         write (outunit, '(/)')
      end if

      ! Print iteration information
      if ( mod(outeritcnt, 10) .eq. 0 ) then
         write (      *, 800)
         write (outunit, 800)
      end if

      write (      *, 810) outeritcnt, mu, func / sf, func, &
                           outerr, inneraccitcnt, reason
      write (outunit, 810) outeritcnt, mu, func / sf, func, &
                           outerr, inneraccitcnt, reason

      ! NONEXECUTABLE STATEMENTS
800   format(/, 1X, "==========================================================", &
             /, 1X, "outer", 1X, "barrier", 1X, "  objective  ", 1X, "    scaled   ", &
                1X, "  outer  ", 1X, "inner",                                         &
             /, 1X, " ite ", 1X, " param ", 1X, "   function  ", 1X, "  obj funct  ", &
                1X, " opt err ", 1X, " ite ")

810   format(1X, I5, 1X, 1P, D7.1, 2(1X, 1P, D13.6), 1X, 1P, D9.3, 1X, I5, A1)

    end subroutine print_iteration

  end subroutine lcmin_solve

  ! ****************************************************************
  ! ****************************************************************

  subroutine lcmin_barrier(mu, n, x, l, u, zl, zu, m, b, lambda, annz, &
       arow, acol, aval, func, grad, laggradnorm, sd, sl, su, flag)

    use constants
    use kktsys_solver
    use seval
    use util
    use variables

    ! SCALAR ARGUMENTS
    integer,       intent(in)    :: annz, m, n
    integer,       intent(out)   :: flag
    real(kind=dp), intent(in)    :: mu
    real(kind=dp), intent(inout) :: func
    real(kind=dp), intent(out)   :: laggradnorm, sd, sl, su

    ! ARRAY ARGUMENTS
    integer,       dimension(annz), intent(in)    :: arow, acol
    real(kind=dp), dimension(annz), intent(in)    :: aval
    real(kind=dp), dimension(m),    intent(in)    :: b
    real(kind=dp), dimension(n),    intent(in)    :: l, u
    real(kind=dp), dimension(n),    intent(inout) :: grad
    real(kind=dp), dimension(m),    intent(inout) :: lambda
    real(kind=dp), dimension(n),    intent(inout) :: x
    real(kind=dp), dimension(nl),   intent(inout) :: zl
    real(kind=dp), dimension(nu),   intent(inout) :: zu

    ! LOCAL SCALARS
    integer       :: btcnt, i, hnnz
    real          :: curtime, itfinish, itstart
    real(kind=dp) :: alpha, alphamax, alphazl, alphazu, barfunc,  &
                     barfuncn, complnorm, compunorm, gtd, inerr,  &
                     lambdanorm, snormsup, tau, xnormsup, zlnorm, &
                     zunorm

    ! LOCAL ARRAYS
    integer,       dimension(hnnzmax) :: hcol, hrow
    real(kind=dp), dimension(hnnzmax) :: hval
    real(kind=dp), dimension(m)       :: dlambda
    real(kind=dp), dimension(n)       :: bgrad, dx, laggrad, xnew
    real(kind=dp), dimension(nl)      :: dzl, lgap
    real(kind=dp), dimension(nu)      :: dzu, ugap

    ! Initialization
    inneritcnt = 0

    xnormsup = 0.0d0

    ! ------------------------------------------------------------------
    !                              MAIN LOOP
    ! ------------------------------------------------------------------
    do
       ! Updates lgap and ugap
       lgap(1:nl) = x(lind(1:nl)) - l(lind(1:nl))
       ugap(1:nu) = u(uind(1:nu)) - x(uind(1:nu))

       ! =============================================================
       !                     VERIFY STOP CRITERIA
       ! =============================================================

       ! Required values to compute the optimality error
       lambdanorm = norm1( lambda(1:m) )
       zlnorm     = norm1( zl(1:nl) )
       zunorm     = norm1( zu(1:nu) )

       if( (m + nl + nu) .gt. 0 ) then
          sd = max(smax, (lambdanorm + zlnorm + zunorm) / (m + nl + nu)) / smax
       end if

       if( nl .gt. 0 ) then
          sl = max(smax, zlnorm / nl) / smax
       end if

       if( nu .gt. 0 ) then
          su = max(smax, zunorm / nu) / smax
       end if

       ! Norm of the lagrangian gradient
       laggrad(1:n) = grad(1:n) + smvp( n, annz, acol, arow, aval, lambda(1:m) )
       laggrad(lind(1:nl)) = laggrad(lind(1:nl)) - zl(1:nl)
       laggrad(uind(1:nu)) = laggrad(uind(1:nu)) + zu(1:nu)
       laggradnorm = infnorm(laggrad(1:n))

       complnorm = infnorm( lgap(1:nl) * zl(1:nl) - mu )
       compunorm = infnorm( ugap(1:nu) * zu(1:nu) - mu )

       ! Optimality error for the barrier subproblem
       inerr = max( laggradnorm / sd, complnorm / sl, compunorm / su )

       ! Print current iteration information
       call print_iteration()

       ! Optimality condition
       if ( inerr .le. kappaeps * mu ) then
          flag = 0
          return

       ! Number of inner iterations
       else if ( inneritcnt .gt. maxinnerit ) then
          flag = - 3
          return

       ! Norm sup of the iterates
       else if ( xnormsup .ge. infty ) then
          flag = - 4
          return

       end if

       ! =============================================================

       ! Evalutes the Hessian at the current point
       call sevalh( n, scalef, sf, hnnz, hrow, hcol, hval, flag )

       ! Evaluates the barrier gradient (needed by kktsys and
       ! backtracking)
       call sevalbg( mu, n, x, grad, l, u, bgrad, flag )
       if ( flag .lt. 0 ) return

       ! Solve the KKT system to find the directions
       call kktsys( mu, n, x, grad, lgap, ugap, m, lambda, zl, zu, &
            annz, arow, acol, aval, hnnz, hrow, hcol, hval, dx,    &
            dlambda, dzl, dzu, flag )
       if(flag .lt. 0) return

       ! Fraction to boundary
       tau = max(tmin, 1 - mu)

       ! Computes step sizes
       alphamax = 1.0d0
       alphazl  = 1.0d0
       alphazu  = 1.0d0

       do i = 1, nl
          if ( dx(lind(i)) .lt. 0.0d0 ) then
             alphamax = min( alphamax, - tau * lgap(i) / dx(lind(i)) )
          end if

          if ( dzl(i) .lt. 0.0d0 ) then
             alphazl = min( alphazl, - tau * zl(i) / dzl(i) )
          end if
       end do
    
       do i = 1, nu
          if ( dx(uind(i)) .gt. 0.0d0 ) then
             alphamax = min( alphamax, tau * ugap(i) / dx(uind(i)) )
          end if

          if ( dzu(i) .lt. 0.0d0 ) then
             alphazu = min( alphazu, - tau * zu(i) / dzu(i) )
          end if
       end do

       ! Backtracking
       btcnt = 0
       alpha = alphamax

       if ( maxval( abs(dx(1:n)) / ( 1.0d0 + abs(x(1:n)) ) ) .lt. 1.0d+01 * macheps ) then
          if ( ( smallsd .ge. 2 ) .and. ( ( mu - 1.0d-01 * epstol ) .le. macheps ) ) then
             flag = - 2
             return
          end if

          smallsd = smallsd + 1

          ! Updates the trial point with the maximum step size
          xnew(1:n) = x(1:n) + alpha * dx(1:n)
          call seval_setpoint( n, xnew )

       else
          smallsd = 0

          call sevalbf( mu, n, x, l, u, .false., scalef, sf, func, &
               barfunc, flag, lgap, ugap )
          if ( flag .ne. 0 ) return

          gtd = dot_product( bgrad(1:n), dx(1:n) )

          do
             ! Updates the trial point
             xnew(1:n) = x(1:n) + alpha * dx(1:n)

             call seval_setpoint( n, xnew )

             call sevalbf( mu, n, xnew, l, u, .true., scalef, sf, func, &
                  barfuncn, flag )
             if ( flag .ne. 0 ) return

             ! Armijo condition
             if ( barfuncn .le. barfunc + gamma * alpha * gtd ) exit

             alpha = 5.0d-01 * alpha
             btcnt = btcnt + 1
          end do
       end if

       ! Computes x sup norm
       xnormsup = infnorm( x(1:n) )

       ! Updates the point
       x(1:n)      = xnew(1:n)
       lambda(1:m) = lambda(1:m) + dlambda(1:m)
       zl(1:nl)    =    zl(1:nl) + alphazl * dzl(1:nl)
       zu(1:nu)    =    zu(1:nu) + alphazu * dzu(1:nu)

       ! Safeguard to the values of zl and zu
       zl(1:nl) = max( min( zl(1:nl), kappaz * (mu / lgap(1:nl)) ), (1.0d0 / kappaz) * (mu / lgap(1:nl)) )
       zu(1:nu) = max( min( zu(1:nu), kappaz * (mu / ugap(1:nu)) ), (1.0d0 / kappaz) * (mu / ugap(1:nu)) )

       ! Evaluates the objective function value and gradient of the
       ! objective funcion at the new point
       call sevalg( scalef, sf, grad, flag )
       if(flag .lt. 0) return

       ! Increments the inner iteration counter
       inneritcnt    = inneritcnt + 1
       inneraccitcnt = inneraccitcnt + 1

    end do

  contains

    subroutine print_iteration()

      ! LOCAL SCALARS
      real          :: ittime
      real(kind=dp) :: minxl, minux, scafeas, unscafeas

      ! Stop timing for inner iteration
      call cpu_time( itfinish )

      if ( inneritcnt .ne. 0 ) then

         ! Prints inner iteration information
         if ( printinit ) then
            if ( mod(inneritcnt-1, 10) .eq. 0 ) then
               write (      *, 400)
               write (outunit, 400)
            end if

            scafeas = infnorm( smvp(m, annz, arow, acol, aval, x) - b )

            minxl = 0.0d0
            minux = 0.0d0

            if ( nl .gt. 0 ) minxl = minval(lgap,1)
            if ( nu .gt. 0 ) minux = minval(ugap,1)

            write (      *, 450) inneritcnt, mu, func, scafeas, minxl, &
                                 minux, pertnw, pertse, norm1(dx),     &
                                 alphamax, alpha, inerr,               &
                                 itfinish - itstart

            write (outunit, 450) inneritcnt, mu, func, scafeas, minxl, &
                                 minux, pertnw, pertse, norm1(dx),     &
                                 alphamax, alpha, inerr,               &
                                 itfinish - itstart

         end if

         ! Print table line
         if ( printtl ) then
            open ( 350, file="lcmin-tabline.out" )

            scafeas = infnorm( smvp(m, annz, arow, acol, aval, x) - b )

            if ( scalea ) then
               unscafeas = infnorm( ( smvp(m, annz, arow, acol, aval, x) - b ) / sa(1:m) )
            else
               unscafeas = scafeas
            end if

            ! Computes current optimization time
            call cpu_time( curtime )

            ! Writes LCMIN flag
            write (350, fmt='(I2)') 50
            write (350, 600) n, m, func, func/sf, scafeas, unscafeas, &
                             inerr, fevalcnt, gevalcnt, hevalcnt,     &
                             bfevalcnt, bgevalcnt, inneraccitcnt,     &
                             outeritcnt, pretime, curtime-opttime
            close ( 350 )
         end if
      end if

      ! Start timing for inner iteration
      call cpu_time( itstart )

      ! NONEXECUTABLE STATEMENTS
400 format(/, 4X, "inner", 1X, "barrier", 1X, "  scaled  ", 1X, "  scaled  ",       &
              1X, "  min  ", 1X, "  min  ", 1X, " NW ", 1X, " SE ", 1X, "   dx   ", &
              1X, "  alpha ", 1X, "  alpha ", 1X, "  inner  ", 1X, "time (s)",      &

           /, 4X, " ite ", 1X, " param ", 1X, " obj func ", 1X, "   feas   ",       &
              1X, "  x-l  ", 1X, "  u-x  ", 1X, "pert", 1X, "pert", 1X, "  norm  ", &
              1X, "   max  ", 1X, "        ", 1X, " opt err ")
                                                                                                                                                            
450 format(4X, I5, 1X, 1P, D7.1, 2(1X, 1P, D10.3), 2(1X, 1P, D7.1), 2(1X, I4), &
           3(1X, 1P, D8.2), 1X, 1P, D9.3, 1X, 0P, F8.2)

600 format(2(I6, 3X), 1P, 2(D24.16, 3X), 3(D11.4, 3X), 7(I7, 3X), 0P, &
           2(F9.2, 3X), "&")

    end subroutine print_iteration

  end subroutine lcmin_barrier

end module lcmin_solver
