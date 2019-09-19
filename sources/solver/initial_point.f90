module initial_point
  
  implicit none

contains

  subroutine ipoint(n, xst, l, u, m, b, annz, arow, acol, aval, sol, &
       flag)

    use algsolver, only: algsolve
    use constants, only: dp, macheps34
    use util,      only: infnorm, smvp

    use variables, only: uepsfeas => epsfeas, cind, kappaip1,         &
                         kappaip2, lind, nl, nu, outunit, sa, scalea, &
                         uind
    
    ! SCALAR ARGUMENTS
    integer, intent(in)  :: annz, m, n
    integer, intent(out) :: flag

    ! ARRAY ARGUMENTS
    integer,       dimension(annz), intent(in)  :: arow, acol
    real(kind=dp), dimension(annz), intent(in)  :: aval
    real(kind=dp), dimension(m),    intent(in)  :: b
    real(kind=dp), dimension(n),    intent(in)  :: l, u, xst
    real(kind=dp), dimension(n),    intent(out) :: sol

    ! ----------------------------------------------------------------
    ! This subroutine tries to compute a strict feasible initial point
    ! to start LCMIN by solving the feasibility problem
    !
    !                    min  | Ax - b |
    !                    s.t. l+leps <= x <= u-ueps
    !
    ! using ALGENCAN. If we are not able to find a feasible point, we
    ! declare the problem is infeasible and stop. If we are not able
    ! to find a strict feasible initial point, we declare the interior
    ! region is void and also stop the process.
    ! ----------------------------------------------------------------

    ! LOCAL SCALARS
    logical       :: origfeas
    real(kind=dp) :: algepsfeas, epsfeas, epsfeas12, factor, feas,  &
                     feasp, minleps, minleps12, minueps, minueps12, &
                     minxl, minux, scafeas
    
    ! LOCAL ARRAYS
    real(kind=dp), dimension(annz) :: unsaval
    real(kind=dp), dimension(m)    :: unsb
    real(kind=dp), dimension(n)    :: difful, lp, up, solp
    real(kind=dp), dimension(nl)   :: leps
    real(kind=dp), dimension(nu)   :: ueps

    write (      *, 300)
    write (outunit, 300)

    ! Real feasibility required
    epsfeas   = max( 1.0d-01 * uepsfeas, macheps34 )
    epsfeas12 = sqrt( epsfeas )

    ! Feasibility to be asked to Algencan
    algepsfeas = epsfeas

    ! if ( scalea ) then
    !    algepsfeas = max( min( 1.0d-01, minval( sa(1:m), 1 ) ) * uepsfeas, macheps34 )
    ! else
    !    algepsfeas = epsfeas
    ! end if

    minleps12 = 0.0d0
    minueps12 = 0.0d0

    ! Unscales the constraints
    if ( scalea ) then
       unsaval(1:annz) = aval(1:annz) / sa(cind(arow(1:annz)))
       unsb(1:m) = b(1:m) / sa(cind(1:m))
    else
       unsaval(1:annz) = aval(1:annz)
       unsb(1:m) = b(1:m)
    end if

    ! Computes the difference between l and u
    difful(1:n) = u(1:n) - l(1:n)

    ! This factor controls the perturbations of the box constraints
    factor = 1.0d0

    origfeas = .false.

900 continue

    ! Computes perturbation on bounds to define the interior region
    leps = min(factor * kappaip1 * max(1.0d0, abs(l(lind))), &
         factor * kappaip2 * difful(lind))
    ueps = min(factor * kappaip1 * max(1.0d0, abs(u(uind))), &
         factor * kappaip2 * difful(uind))

    minleps = minval(leps, 1)
    minueps = minval(ueps, 1)

    if ( nl .gt. 0 ) minleps12 = max( minleps ** 1.5d0, macheps34 )
    if ( nu .gt. 0 ) minueps12 = max( minueps ** 1.5d0, macheps34 )

    minxl = minleps
    minux = minueps

    lp = l
    up = u
    lp(lind) = l(lind) + leps
    up(uind) = u(uind) - ueps

    ! Computes the initial point using Algencan
    sol = xst

    write (*, 310)
    write (outunit, 310)

    ! COMPUTES THE INITIAL POINT USING TRUNCATED NEWTON UNTIL REACH,
    ! AT LEAST, THE SQUARE ROOT OF THE FEASIBILITY REQUIRED
    call algsolve( n, sol, lp, up, m, unsb, annz, arow, acol, unsaval, &
         algepsfeas ** 5.0d-01, 'tn', flag )

    ! Evaluates the feasibility of the initial point
    feas = infnorm( smvp(m, annz, arow, acol, unsaval, sol) - unsb )

    ! If the feasibility os not the desired, uses Newton method in Algencan
    if ( feas .gt. epsfeas ) then
       write (      *, 320)
       write (outunit, 320)

       ! COMPUTES THE INITIAL POINT USING NEWTON METHOD UNTIL REACH,
       ! AT LEAST, THE FEASIBILITY REQUIRED
       call algsolve( n, sol, lp, up, m, unsb, annz, arow, acol, &
            unsaval, algepsfeas, 'nw', flag )

       ! Evaluates the feasibility of the initial point
       feas = infnorm( smvp( m, annz, arow, acol, unsaval, sol ) - unsb )
    end if
 
    ! Checks for feasibility of the bounds for the initial point found
    if(nl .gt. 0) minxl = minval(sol(lind) - l(lind), 1)
    if(nu .gt. 0) minux = minval(u(uind) - sol(uind), 1)

    if ( feas .gt. epsfeas .and. .not. origfeas ) then

       write (      *, 350) feas, epsfeas
       write (outunit, 350) feas, epsfeas

       if(nl .gt. 0) then
          write (      *, 360) minxl, minleps12
          write (outunit, 360) minxl, minleps12
       end if

       if(nu .gt. 0) then
          write (      *, 370) minux, minueps12
          write (outunit, 370) minux, minueps12
       end if

       write (      *, 200)
       write (outunit, 200)

       ! Save prior solutions
       feasp = feas
       solp(1:n) = sol(1:n)

       ! TRY TO SOLVE THE FEASIBILITY PROBLEM WITH THE ORIGINAL BOUNDS
       call algsolve( n, sol, l, u, m, unsb, annz, arow, acol, unsaval, &
            algepsfeas, 'tn', flag )

       ! Evaluates the feasibility of the initial point
       feas = infnorm( smvp(m, annz, arow, acol, unsaval, sol) - unsb )

       ! Checks for feasibility of the bounds for the initial point found
       if(nl .gt. 0) minxl = minval(sol(lind) - l(lind), 1)
       if(nu .gt. 0) minux = minval(u(uind) - sol(uind), 1)

       ! IF THE ORIGINAL PROBLEM IS INFEASIBLE
       if ( feas .gt. epsfeas ) then

          ! If the original problem is not acceptably feasible, stop
          if ( feas .gt. epsfeas12 ) then
             write (      *, 350) feas, epsfeas
             write (outunit, 350) feas, epsfeas

             flag = - 6
             return

          else if ( minxl .lt. minleps12 .or. minux .lt. minueps12 ) then

             write (      *, 350) feas, epsfeas
             write (outunit, 350) feas, epsfeas

             if(nl .gt. 0) then
                write (      *, 360) minxl, minleps12
                write (outunit, 360) minxl, minleps12
             end if

             if(nu .gt. 0) then
                write (      *, 370) minux, minueps12
                write (outunit, 370) minux, minueps12
             end if

             ! The original problem is acceptable feasible but the
             ! point found is not interior.
             origfeas = .true.

             factor = min( kappaip1, kappaip2 ) * factor

             write (      *, 450)
             write (      *, 500)

             write (outunit, 450)
             write (outunit, 500)

             go to 900

          else
             write (      *, 210)
             write (outunit, 210)

          end if

       else if ( minxl .lt. minleps12 .or. minux .lt. minueps12 ) then

          write (      *, 350) feas, epsfeas
          write (outunit, 350) feas, epsfeas

          if(nl .gt. 0) then
             write (      *, 360) minxl, minleps12
             write (outunit, 360) minxl, minleps12
          end if

          if(nu .gt. 0) then
             write (      *, 370) minux, minueps12
             write (outunit, 370) minux, minueps12
          end if

          ! At this point, the feasibility of the original problem has
          ! just been checked, and it is feasible.
          origfeas = .true.

          factor = min( kappaip1, kappaip2 ) * factor

          write (      *, 450)
          write (      *, 500)

          write (outunit, 450)
          write (outunit, 500)

          go to 900

       end if

    else if ( feas .gt. epsfeas .or. minxl .lt. minleps12 &
                                .or. minux .lt. minueps12 ) then

       ! At this point, the feasibility of the original point has
       ! already been checked, but with the perturbations on the box
       ! constraints, something gone wrong.

       if ( feas .le. epsfeas12 .and. minxl .ge. minleps12 &
                                .and. minux .ge. minueps12 ) then

          ! In this case, the point found is acceptably interior and
          ! feasible.
          write (      *, 210)
          write (outunit, 210)

       else
          factor = min( kappaip1, kappaip2 ) * factor

          write (      *, 350) feas, epsfeas
          write (outunit, 350) feas, epsfeas

          if(nl .gt. 0) then
             write (      *, 360) minxl, minleps12
             write (outunit, 360) minxl, minleps12
          end if

          if(nu .gt. 0) then
             write (      *, 370) minux, minueps12
             write (outunit, 370) minux, minueps12
          end if

          if( factor * min( kappaip1, kappaip2 ) .le. macheps34 ) then
             flag = - 5
             return
          else
             write (      *, 500)
             write (outunit, 500)
            
             go to 900
          end if
       end if
    end if

    ! Writes success message
    write (      *, 400)
    write (outunit, 400)

    write (      *, 350) feas, epsfeas
    write (outunit, 350) feas, epsfeas

    if(nl .gt. 0) then
       write (      *, 360) minxl, minleps12
       write (outunit, 360) minxl, minleps12
    end if

    if(nu .gt. 0) then
       write (      *, 370) minux, minueps12
       write (outunit, 370) minux, minueps12
    end if

    ! NONEXECUTABLE STATEMENTS
200 format(/, 1X, "The initial point found is not feasible and/or is not interior.", &
           /, 1X, "We will check if there is any feasible point for the original constraints.")

210 format(/, 1X, "The initial point found is acceptably feasible. LCMIN will continue.")

300 format(/, 1X, "Computing a suitable initial point to start LCMIN.")

310 format(/, 1X, "Adjusting feasibility using Algencan with Truncated Newton Method.")
320 format(/, 1X, "Adjusting feasibility using Algencan with Newton Method.")

350 format(/, 36X, "        achieved        ", 1X, "         desired        ", &
           /, 1X, "Unscaled feasibility ..............: ", 1P, D24.16, 1X, D24.16)
360 format(   1X, "Smallest distance to lower bound ..: ", 1P, D24.16, 1X, D24.16)
370 format(   1X, "Smallest distance to upper bound ..: ", 1P, D24.16, 1X, D24.16)

400 format(/, 1X, "Initial point computed successfully.")

450 format(/, 1X, "The original problem is feasible.")

500 format(/, 1X, "We will try to compute a new initial point with smaller perturbations.")

  end subroutine ipoint

  ! ------------------------------------------------------------------------

  subroutine imult(m, n, grad, zl, zu, annz, arow, acol, aval, lambda)

    use constants, only: dp
    use lsys_solver
    use variables, only: lind, nl, nu, uind

    ! SCALAR ARGUMENTS
    integer, intent(in) :: annz, m, n

    ! ARRAY ARGUMENTS
    integer,       dimension(annz), intent(in)  :: arow, acol
    real(kind=dp), dimension(annz), intent(in)  :: aval
    real(kind=dp), dimension(n),    intent(in)  :: grad
    real(kind=dp), dimension(nl),   intent(in)  :: zl
    real(kind=dp), dimension(nu),   intent(in)  :: zu
    real(kind=dp), dimension(m),    intent(out) :: lambda

    ! ----------------------------------------------------------------
    ! This subroutine computes an initial estimation to the lagrange
    ! multipliers.
    ! ----------------------------------------------------------------

    ! LOCAL SCALARS
    integer       :: cnt, dim, dum, flag
    real(kind=dp) :: residnorm

    ! LOCAL ARRAYS
    integer,       dimension(annz+n) :: msysrow, msyscol
    real(kind=dp), dimension(annz+n) :: msysval
    real(kind=dp), dimension(n+m)    :: rhs
    
    dim = n+m

    msysrow(1:n) = (/ (    cnt, cnt = 1, n) /)
    msyscol(1:n) = (/ (    cnt, cnt = 1, n) /)
    msysval(1:n) = (/ (- 1.0d0, cnt = 1, n) /)

    msysrow(n+1:n+annz) = arow + n
    msyscol(n+1:n+annz) = acol
    msysval(n+1:n+annz) = aval

    rhs(1:n)     = - grad
    rhs(lind)    = rhs(lind) + zl
    rhs(uind)    = rhs(uind) - zu
    rhs(n+1:n+m) = (/ (0.0d0, cnt = 1, m) /)

    ! Initialize
    call lsys_initialize()

    ! Analyse
    call lsys_analyse( dim, annz+n, msysrow, msyscol, msysval, flag )
    if ( flag .lt. 0 ) return

    ! Factorize
    call lsys_factorize( dim, dum, dum, flag )
    if ( flag .lt. 0 ) return

    ! Solve
    call lsys_solve( dim, rhs, flag )
    if ( flag .lt. 0 ) return

    lambda = rhs(n+1:n+m)

  end subroutine imult

end module initial_point
