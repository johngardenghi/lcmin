module seval

  use constants, only: bignum, dp

  implicit none

  save

  ! Global pointer to user variable x
  real(kind=dp), allocatable, dimension(:) :: x_user  

contains

  subroutine seval_setpoint(n, x)

    use variables, only: nfind

    ! SCALAR ARGUMENT
    integer, intent(in) :: n

    ! ARRAY ARGUMENT
    real(kind=dp), dimension(n), intent(in) :: x

    x_user(nfind(1:n)) = x(1:n)

  end subroutine seval_setpoint

  ! -------------------------------------------------------------------

  subroutine sevalf(scale, sf, f, flag)

    use user_subs, only: subf
    use variables, only: fevalcnt, nprob, outunit

    ! SCALAR ARGUMENTS
    logical,       intent(in)  :: scale
    integer,       intent(out) :: flag
    real(kind=dp), intent(in)  :: sf
    real(kind=dp), intent(out) :: f

    ! LOCAL SCALAR
    integer       :: i
    real(kind=dp) :: ftmp

    call subf( nprob, x_user, ftmp, flag )

    if( flag .ne. 0 ) then
       write (      *, 100) flag
       write (outunit, 100) flag

       flag = - 8
       return
    end if

    if( isnan(ftmp) ) then
       write (      *, 110)
       write (outunit, 110)

       flag = - 8
       return
    end if

    if( abs(ftmp) .ge. bignum ) then
       write (      *, 120)
       write (outunit, 120)
    end if

    f = ftmp

    if ( scale ) then
       f = sf * f
    end if

    fevalcnt = fevalcnt + 1
 
    ! NONEXECUTABLE STATEMENTS
100 format(/, 1X, "ERROR (sevalf): There was a problem while evaluating the value of f.", &
           /, 1X, "                Your subroutine reports ", I3, ".")

110 format(/, 1X, "ERROR (sevalf): The objective function value is NaN.")

120 format(/, 1X, "WARNING (sevalf): The objective function value is -Inf or Inf.")

  end subroutine sevalf

  ! -------------------------------------------------------------------

  subroutine sevalg(scale, sf, g, flag)

    use user_subs, only: subg
    use variables, only: gevalcnt, nfind, nprob, outunit

    implicit none

    ! SCALAR ARGUMENT
    integer,       intent(out) :: flag
    logical,       intent(in)  :: scale
    real(kind=dp), intent(in)  :: sf

    ! ARRAY ARGUMENTS
    real(kind=dp), dimension(:), intent(out) :: g

    ! LOCAL ARRAY
    real(kind=dp), dimension(nprob) :: grad

    call subg( nprob, x_user, grad, flag )

    if( flag .ne. 0 ) then
       write (      *, 100) flag
       write (outunit, 100) flag

       flag = - 8
       return
    end if

    if( any( isnan( grad(nfind) ) ) .or. any( abs( grad(nfind) ) .ge. bignum ) ) then
       write (      *, 110)
       write (outunit, 110)

       flag = - 8
       return
    end if

    if( scale ) then
       g = sf * grad(nfind)
    else
       g = grad(nfind)
    end if

    gevalcnt = gevalcnt + 1

    ! NONEXECUTABLE STATEMENTS
100 format(/, 1X, "ERROR (sevalg): There was a problem while evaluating the gradient of f.", &
           /, 1X, "                Your subroutine reports ", I3, ".")

110 format(/, 1X, "ERROR (sevalg): The objective function gradient is -Inf, Inf or NaN.")

  end subroutine sevalg

  ! -------------------------------------------------------------------

  subroutine sevalh(n, scale, sf, hnnz, hrow, hcol, hval, flag)

    use user_subs, only: subh
    use util,      only: swap
    use variables, only: find, hevalcnt, hnnzmax, nprob, outunit

    ! SCALAR ARGUMENTS
    integer,       intent(in)  :: n
    integer,       intent(out) :: flag, hnnz
    logical,       intent(in)  :: scale
    real(kind=dp), intent(in)  :: sf

    ! ARRAY ARGUMENTS
    integer,       dimension(hnnzmax), intent(out) :: hcol, hrow
    real(kind=dp), dimension(hnnzmax), intent(out) :: hval

    ! LOCAL SCALAR
    integer :: cnt

    call subh( nprob, x_user, hnnzmax, hrow, hcol, hval, hnnz, flag )

    if(flag .ne. 0) then
       write (      *, 100) flag
       write (outunit, 100) flag

       flag = - 8
       return
    end if

    ! If there are fixed variables, reorder hessian to ignore
    ! derivatives to respect to fixed variables
    if( nprob .gt. n ) then
       cnt = 1
       do while( cnt .le. hnnz )
          ! If any variable is a fixed one, ignored the value of the
          ! hessian, swapping with the last elements and decrementing
          ! hnnz
          if ( ( find(hrow(cnt)) .eq. 0 ) .or. &
               ( find(hcol(cnt)) .eq. 0 ) ) then
             call swap( hrow(cnt), hrow(hnnz) )
             call swap( hcol(cnt), hcol(hnnz) )
             call swap( hval(cnt), hval(hnnz) )
             hnnz = hnnz - 1
          else
             hrow(cnt) = find(hrow(cnt))
             hcol(cnt) = find(hcol(cnt))
             cnt = cnt + 1
          end if
       end do
    end if

    if ( any( isnan( hval ) ) .or. any( abs( hval ) .ge. bignum ) ) then
       write (      *, 110)
       write (outunit, 110)

       flag = - 8
       return
    end if

    ! Scales the hessian
    if ( scale ) then
       hval(1:hnnz) = sf * hval(1:hnnz)
    end if

    hevalcnt = hevalcnt + 1

    ! NONEXECUTABLE STATEMENTS
100 format(/, 1X, "ERROR (sevalh): There was a problem while evaluating the hessian of f.", &
           /, 1X, "                Your subroutine reports ", I3, ".")

110 format(/, 1X, "ERROR (sevalh): The objective function hessian is -Inf, Inf or NaN.")

  end subroutine sevalh

  ! -------------------------------------------------------------------

  subroutine sevalbf(mu, n, x, l, u, evalf, scale, s, f, bf, flag, &
       lgap, ugap)
    
    use constants, only: dp
    use variables, only: bfevalcnt, kappad, lind, lnu, nl, nprob, nu, &
                         outunit, uind, unl

    ! SCALAR ARGUMENTS
    integer,       intent(in)    :: n
    integer,       intent(out)   :: flag
    logical,       intent(in)    :: evalf, scale
    real(kind=dp), intent(in)    :: mu, s
    real(kind=dp), intent(inout) :: f
    real(kind=dp), intent(out)   :: bf
    
    ! ARRAY ARGUMENTS
    real(kind=dp), dimension(n),  intent(in)           :: l, u, x
    real(kind=dp), dimension(nl), intent(in), optional :: lgap
    real(kind=dp), dimension(nu), intent(in), optional :: ugap
  
    ! -----------------------------------------------------
    ! This subroutine evaluates the barrier function at x
    ! -----------------------------------------------------

    if ( evalf ) then
       call sevalf( scale, s, f, flag )
       if ( flag .ne. 0 ) return
    end if

    if ( present( lgap ) .and. present( ugap ) ) then
       bf = f                                            &
            -          mu * ( sum( log( lgap(1:nl) ) )   &
                            + sum( log( ugap(1:nu) ) ) ) &
            + kappad * mu * ( sum( x(lnu) - l(lnu) )     &
                            + sum( u(unl) - x(unl) ) )
    else
       bf = f                                                   &
            -          mu * ( sum( log( x(lind) - l(lind) ) )   &
                            + sum( log( u(uind) - x(uind) ) ) ) &
            + kappad * mu * ( sum( x(lnu) - l(lnu) )            &
                            + sum( u(unl) - x(unl) ) )
    end if

    if(isnan(bf)) then
       write (      *, 100)
       write (outunit, 100)

       flag = - 8
       return
    end if

    bfevalcnt = bfevalcnt + 1

    ! NONEXECUTABLE STATEMENTS
100 format(/, 1X, "ERROR (evalbf): The barrier function value is NaN.", /)

  end subroutine sevalbf

  ! -------------------------------------------------------------------

  subroutine sevalbg(mu, n, x, g, l, u, bg, flag)

    use constants, only: dp
    use variables, only: bgevalcnt, kappad, lind, lnu, nl, nu, &
                         outunit, uind, unl

    ! SCALAR ARGUMENTS
    integer,       intent(in)  :: n
    integer,       intent(out) :: flag
    real(kind=dp), intent(in)  :: mu

    ! ARRAY ARGUMENTS
    real(kind=dp), dimension(n),  intent(in)  :: g, l, u, x
    real(kind=dp), dimension(n),  intent(out) :: bg

    ! ----------------------------------------------------------------
    ! This subroutine evaluates the barrier function gradient at x
    ! given the gradient of the objective function g at x.
    ! ----------------------------------------------------------------

    bg(1:n)  = g(1:n)
    bg(lind) = bg(lind) - mu * ( 1.0d0 / ( x(lind) - l(lind) ) )
    bg(uind) = bg(uind) + mu * ( 1.0d0 / ( u(uind) - x(uind) ) )
    bg(lnu)  = bg(lnu)  + kappad * mu
    bg(unl)  = bg(unl)  - kappad * mu

    if( any( isnan( bg(1:n) ) ) ) then
       write (      *, 100)
       write (outunit, 100)

       flag = - 8
       return
    end if

    bgevalcnt = bgevalcnt + 1

    ! NONEXECUTABLE STATEMENTS
100 format(/, 1X, "ERROR (evalbg): The barrier function gradient is NaN.")

  end subroutine sevalbg

end module seval
