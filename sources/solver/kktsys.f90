module kktsys_solver

  use constants, only: dp

  implicit none

  save

  ! INERTIA CORRECTION PARAMETERS
  integer       :: conscor = 0
  real(kind=dp) :: lchi    = 0.0d0

contains

  subroutine kktsys(mu, n, x, g, lgap, ugap, m, lambda, zl, zu, annz, &
       arow, acol, aval, hnnz, hrow, hcol, hval, dx, dlambda, dzl,    &
       dzu, flag)

    use constants, only: macheps, macheps12
    use lsys_solver
    use util,      only: smvp
    use variables, only: chiini, chimin, kchibar, kchimax, kchimin, &
                         lind, nl, nu, pertnw, pertse, uind

    ! SCALAR ARGUMENTS
    integer,       intent(in)  :: annz, hnnz, m, n
    integer,       intent(out) :: flag
    real(kind=dp), intent(in)  :: mu

    ! ARRAY ARGUMENTS
    integer,       dimension(annz), intent(in)  :: arow, acol
    integer,       dimension(hnnz), intent(in)  :: hrow, hcol
    real(kind=dp), dimension(annz), intent(in)  :: aval
    real(kind=dp), dimension(hnnz), intent(in)  :: hval
    real(kind=dp), dimension(n),    intent(in)  :: g, x
    real(kind=dp), dimension(m),    intent(in)  :: lambda
    real(kind=dp), dimension(nl),   intent(in)  :: lgap, zl
    real(kind=dp), dimension(nu),   intent(in)  :: ugap, zu
    real(kind=dp), dimension(n),    intent(out) :: dx
    real(kind=dp), dimension(m),    intent(out) :: dlambda
    real(kind=dp), dimension(nl),   intent(out) :: dzl
    real(kind=dp), dimension(nu),   intent(out) :: dzu

    ! ----------------------------------------------------------------
    ! This subroutine solves the Primal-Dual system. It builds the
    ! coefficient matrix and calls the MA57 solver.
    ! ----------------------------------------------------------------

    ! LOCAL SCALARS
    integer       :: allocstat, cnt, msysnnz, nneigv, rank, shift
    real(kind=dp) :: chi, chiold

    ! LOCAL ARRAYS
    integer,       dimension(annz+hnnz+m+n) :: msysrow, msyscol
    real(kind=dp), dimension(annz+hnnz+m+n) :: msysval
    real(kind=dp), dimension(m+n)           :: rhs
    real(kind=dp), dimension(nl)            :: sigmal
    real(kind=dp), dimension(nu)            :: sigmau
    real(kind=dp), dimension(n)             :: xlinv, xuinv

    flag = 0

    ! The system matrix has to have the null diagonal entries (that
    ! hessian probably does not have) in order to sum the necessary
    ! terms in the diagonal
    msysnnz = hnnz + annz + n

    xlinv(1:n) = 0.0d0
    xuinv(1:n) = 0.0d0

    xlinv(lind) = 1.0d0 / lgap(1:nl)
    xuinv(uind) = 1.0d0 / ugap(1:nu)

    sigmal(1:nl) = xlinv(lind) * zl(1:nl)
    sigmau(1:nu) = xuinv(uind) * zu(1:nu)

    ! Builds the system reduced matrix, of dimension n+m

    ! The first n positions of msysrow, msyscol and msysval will hold
    ! the diagonal entries of the hessian, including the null
    ! ones. The following m will hold the diagonal of southeast.
    msysrow(1:n) = (/ (cnt, cnt=1, n) /)
    msyscol(1:n) = (/ (cnt, cnt=1, n) /)
    msysval(1:n) = 0.0d0

    shift = 1
    do cnt = 1, hnnz
       if ( hrow(cnt) .eq. hcol(cnt) ) then
          msysval(hrow(cnt)) = hval(cnt)
          msysnnz = msysnnz - 1
       else
          msysrow(n+shift) = hrow(cnt)
          msyscol(n+shift) = hcol(cnt)
          msysval(n+shift) = hval(cnt)
          shift = shift + 1
       end if
    end do

    ! Increase the northwest with the primal-dual hessian
    msysval(lind) = msysval(lind) + sigmal(1:nl)
    msysval(uind) = msysval(uind) + sigmau(1:nu)

    msysrow(msysnnz-annz+1:msysnnz) = acol(1:annz)
    msyscol(msysnnz-annz+1:msysnnz) = arow(1:annz) + n
    msysval(msysnnz-annz+1:msysnnz) = aval(1:annz)

    rhs(1:n) = - ( g(1:n) - mu * ( xlinv(1:n) - xuinv(1:n) ) )
    rhs(n+1:n+m) = 0.0d0

    ! Counts the perturbations on the primal-dual system
    pertnw = 0
    pertse = 0

    ! Initialize
    call lsys_initialize()

    ! Analyse
    call lsys_analyse( n+m, msysnnz, msysrow, msyscol, msysval, flag )
    if ( flag .lt. 0 ) return

    ! Factorize
    if ( conscor .lt. 3 ) then
       call lsys_factorize( n+m, nneigv, rank, flag )
       if ( flag .lt. 0 ) return
    else
       rank   = n+m-1
       nneigv = m
    end if

    ! INERTIA CORRECTION
    if ( rank .lt. n+m .or. nneigv .ne. m ) then

       ! If there are less that m negative eigenvalues, perturb the
       ! southeast of the matrix
       if( nneigv .lt. m ) then
          msysrow(msysnnz+1:msysnnz+m) = (/ (cnt+n, cnt = 1, m) /)
          msyscol(msysnnz+1:msysnnz+m) = (/ (cnt+n, cnt = 1, m) /)
          msysval(msysnnz+1:msysnnz+m) = - macheps12 * (mu ** (1.0d0/4.0d0))
          msysnnz = msysnnz + m

          pertse = pertse + 1

          ! Analyse
          call lsys_analyse( n+m, msysnnz, msysrow, msyscol, msysval, flag )
          if ( flag .lt. 0 ) return
       end if

       if( abs( lchi ) .le. macheps ) then
          chi = chiini
       else
          chi = max( chimin, kchimin * lchi )
       end if

       pertnw  = pertnw + 1
       conscor = conscor + 1

       ! Factorize
       call lsys_factorize( n+m, nneigv, rank, flag, n=n, diag=chi )
       if ( flag .lt. 0 ) return

       do while ( rank .lt. n+m .or. nneigv .ne. m )

          chiold = chi
          
          if( abs(lchi) .le. macheps ) then
             chi = kchibar * chi
          else
             chi = kchimax * chi
          end if

          pertnw = pertnw + 1

          if ( isnan(chi-chiold) ) stop

          ! Factorize
          call lsys_factorize( n+m, nneigv, rank, flag, n=n, diag=(chi-chiold) )
          if ( flag .lt. 0 ) return

       end do
       
       lchi = chi
    else
       conscor = 0
    end if

    ! Solve the system
    call lsys_solve( n+m, rhs, flag )
    if ( flag .lt. 0 ) return

    ! Set the directions
    dx(1:n)      = rhs(1:n)
    dlambda(1:m) = rhs(n+1:n+m) - lambda(1:m)
    dzl(1:nl)    = - (zl(1:nl) - mu * xlinv(lind) + sigmal(1:nl) * dx(lind))
    dzu(1:nu)    = -  zu(1:nu) + mu * xuinv(uind) + sigmau(1:nu) * dx(uind)

    ! if ( maxval( abs( sigmal(1:nl) * dx(lind) ) ) .lt. 10.0d0 * macheps .or. &
    !      maxval( abs( sigmau(1:nu) * dx(uind) ) ) .lt. 10.0d0 * macheps ) then
    !    print *, "ALGUEM PEQUENO AQUI"
    ! end if

    ! NONEXECUTABLE STATEMENTS
100 format(/, 1X, "ERROR (kktsys): Unable to allocate memory.")

200 format(/, 1X, "WARNING (kktsys): Adding perturbations to the southeast", &
           /, 1X, "                  of the primal-dual system.")

  end subroutine kktsys

end module kktsys_solver
