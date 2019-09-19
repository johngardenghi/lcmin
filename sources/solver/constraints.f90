module constraints

  use constants, only: dp

  implicit none

  save

  ! GLOBAL SCALAR
  integer :: ne

  ! GLOBAL POINTERS
  integer,       pointer, dimension(:), private :: ac, ar
  real(kind=dp), pointer, dimension(:), private :: av, bc

contains

  subroutine set_constraints(m, b, annz, arow, acol, aval)

    ! SCALAR ARGUMENT
    integer, intent(in) :: annz, m

    ! ARRAY ARGUMENTS
    integer,       dimension(annz), intent(in), target :: acol, arow
    real(kind=dp), dimension(annz), intent(in), target :: aval
    real(kind=dp), dimension(m),    intent(in), target :: b

    ne = annz

    ac => acol
    ar => arow
    av => aval
    bc => b

  end subroutine set_constraints

  ! -------------------------------------------------------------------

  function eval_feas(x, scaled) result( feas )

    use util,      only: infnorm, smvp
    use variables, only: sa, scalea

    ! FUNCTION RESULT
    real(kind=dp) :: feas

    ! SCALAR ARGUMENT
    logical, intent(in) :: scaled

    ! ARRAY ARGUMENT
    real(kind=dp), dimension(:), intent(in) :: x

    if ( .not. scaled .or. .not. scalea ) then
       feas = infnorm( smvp( size(sa), ne, ar, ac, av, x ) - bc )
    else
       feas = infnorm( sa * ( smvp( size(sa), ne, ar, ac, av, x ) - bc ) )
    end if

  end function eval_feas

  ! -------------------------------------------------------------------

  subroutine check_implicit_equalities(n, l, u, m, b, annz, arow, &
       acol, aval, flag)

    use constants, only: macheps, macheps12, macheps34
    use util,      only: swap
    use variables, only: infty, hslunit, outunit

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)  :: annz, m, n
    integer, intent(out) :: flag

    ! ARRAY ARGUMENTS
    integer,       dimension(annz), intent(in)    :: arow, acol
    real(kind=dp), dimension(annz), intent(in)    :: aval
    real(kind=dp), dimension(m),    intent(in)    :: b
    real(kind=dp), dimension(n),    intent(inout) :: l, u

    ! ------------------------------------------------------------------
    ! This subroutine detects implicit equalities in the box
    ! constraints by intervalar analysis in the linear equality
    ! constraints. Then, l and u are altered so that the implicit
    ! equalities rises with l(i) = u(i).
    !
    ! CALLS : MC59AD
    ! ------------------------------------------------------------------

    ! LOCAL SCALARS
    integer       :: allocstat, i, j, k, la, liw, nfx, nvar
    logical       :: changed
    real(kind=dp) :: cl, cu, rhs

    ! LOCAL ARRAYS
    integer, allocatable, dimension(:) :: iw

    integer,       dimension(10)   :: icntl, info
    integer,       dimension(m+1)  :: ip
    integer,       dimension(annz) :: irn, jcn
    real(kind=dp), dimension(annz) :: a
    real(kind=dp), dimension(n)    :: lorig, lprev, ltmp, newl, newu, uorig, &
                                      uprev, utmp, val

    ! EXTERNAL SUBROUTINES
    external :: mc59ad

    flag = 0
    liw = max(m,n) + 1

    allocate(iw(liw), stat=allocstat)
    if(allocstat .ne. 0) then
       write (*, 100)
       write (outunit, 100)

       flag = - 10
       return
    end if

    la = 0 
    do i = 1, annz
       if( abs( aval(i) ) .gt. macheps ) then
          la = la + 1
          irn(la) = arow(i)
          jcn(la) = acol(i)
          a(la)   = aval(i)
       end if
    end do

    open(hslunit, file='mc59.out')

    ! Sort the entries of A by rows
    icntl(1) = 0
    icntl(2) = 0
    icntl(3) = 0
    icntl(4) = hslunit
    icntl(5) = hslunit
    icntl(6) = 0
    call mc59ad(icntl, m, n, la, jcn, la, irn, la, a, m+1, ip, liw, iw, info)

    close(hslunit)

    if(info(1) .lt. 0) then
       write (*, 200)
       write (outunit, 200)

       flag = - 9
       return
    end if

    ! Computes the number of fixed variables
    nfx = count( u(1:n) - l(1:n) .le. macheps )

    lorig(1:n) = l(1:n)
    uorig(1:n) = u(1:n)

    do
       lprev(1:n) = l(1:n)
       uprev(1:n) = u(1:n)

       ! Check each constraint
       do i = 1, m
          rhs = b(i)

          nvar = ip(i+1) - ip(i)
          ! val(1:nvar) = a(ip(i):ip(i+1)-1)
          ! ltmp(1:nvar) = l(jcn(ip(i):ip(i+1)-1))
          ! utmp(1:nvar) = u(jcn(ip(i):ip(i+1)-1))

          do k = 1, nvar
             if( a(ip(i)+k-1) .lt. 0.0d0 ) then
                val(k)  = - a(ip(i)+k-1)
                ltmp(k) = - l(jcn(ip(i)+k-1))
                utmp(k) = - u(jcn(ip(i)+k-1))
                call swap( ltmp(k), utmp(k) )
             else
                val(k)  = a(ip(i)+k-1)
                ltmp(k) = l(jcn(ip(i)+k-1))
                utmp(k) = u(jcn(ip(i)+k-1))
             end if
          end do

          ! Computes new bounds newl(1:nvar) and newu(1:nvar) for the
          ! nvar variables of the i-th constraint
          call newBounds( nvar, ltmp, utmp, val, rhs, newl, newu )

          where( a(ip(i):ip(i+1)-1) .lt. 0 )
             l(jcn(ip(i):ip(i+1)-1)) = max( l(jcn(ip(i):ip(i+1)-1)), - newu(1:nvar) )
             u(jcn(ip(i):ip(i+1)-1)) = min( u(jcn(ip(i):ip(i+1)-1)), - newl(1:nvar) )
          elsewhere
             l(jcn(ip(i):ip(i+1)-1)) = max( l(jcn(ip(i):ip(i+1)-1)), newl(1:nvar) )
             u(jcn(ip(i):ip(i+1)-1)) = min( u(jcn(ip(i):ip(i+1)-1)), newu(1:nvar) )
          end where

          if( any( ( u(jcn(ip(i):ip(i+1)-1)) - l(jcn(ip(i):ip(i+1)-1)) ) .lt. - macheps12 ) ) then
             flag = - 6
             return
          end if

       end do

       if ( max( norm2( lprev(1:n) - l(1:n) ), norm2( uprev(1:n) - u(1:n) ) ) .le. macheps12 ) then
          where ( abs( u(1:n) - l(1:n) ) .le. macheps34 )
             u(1:n) = l(1:n)
          end where

          exit
       end if
    end do
          
    ! Computes the number of fixed variables found by the algorithm
    nfx = count( abs( u(1:n) - l(1:n) ) .le. macheps ) - nfx

    ! Restores original values of the bounds that were not identified
    ! as implicit equalities
    where( abs( u(1:n) - l(1:n) ) .gt. macheps )
       l(1:n) = lorig(1:n)
       u(1:n) = uorig(1:n)
    end where

    if( nfx .eq. 0 ) then
       write (*, 300)
       write (outunit, 300)
    else
       write (*, 350) nfx
       write (outunit, 350) nfx
    end if

999 deallocate(iw)

    ! NON EXECUTABLE STATEMENTS
100 format(/, 1X, "ERROR (implicit_equalities): Unable to allocate memory")
    
200 format(/, 1X, "ERROR (implicit_equalities): MC59AD executed with errors. Please", &
           /, 1X, "                             check mc59.out file for details.")

300 format(/, 1X, "No implicit equalities found.")

350 format(/, 1X, "Number of implicit equalities found: ", I6)

  contains
    
    subroutine newBounds(nvar, l, u, val, rhs, newl, newu)

      ! SCALAR ARGUMENTS
      integer,       intent(in)  :: nvar
      real(kind=dp), intent(in)  :: rhs

      ! ARRAY ARGUMENTS
      real(kind=dp), dimension(nvar), intent(in)  :: l, u, val
      real(kind=dp), dimension(nvar), intent(out) :: newl, newu

      ! LOCAL SCALARS
      integer       :: infvar, k, minusinf, plusinf
      real(kind=dp) :: cl, cu

      minusinf = count( l .le. - infty )
      plusinf  = count( u .ge.   infty )

      ! Computes the lower bound
      if ( minusinf .eq. 0 ) then
         cl = - rhs + sum( val(1:nvar) * l(1:nvar) )
         newu(1:nvar) = ( l(1:nvar) - cl / val(1:nvar) )

      else if( minusinf .eq. 1 ) then
         cl = - rhs

         do k = 1, nvar
            if( l(k) .gt. - infty ) then
               cl = cl + val(k) * l(k)
               newu(k) = u(k)
            else
               infvar = k
            end if
         end do

         newu(infvar) = - cl / val(infvar)

      else
         newu(1:nvar) = u(1:nvar)
      end if

      ! Computes the upper bound
      if ( plusinf .eq. 0 ) then
         cu = - rhs + sum( val(1:nvar) * u(1:nvar) )
         newl(1:nvar) = ( u(1:nvar) - cu / val(1:nvar) )

      else if( plusinf .eq. 1 ) then
         cu = - rhs

         do k = 1, nvar
            if( u(k) .lt. infty ) then
               cu = cu + val(k) * u(k)
               newl(k) = l(k)
            else
               infvar = k
            end if
         end do

         newl(infvar) = - cu / val(infvar)

      else
         newl(1:nvar) = l(1:nvar)
      end if

    end subroutine newBounds

  end subroutine check_implicit_equalities

  ! -------------------------------------------------------------------

  subroutine scale_constraints(m, b, annz, arow, acol, aval, flag)

    use constants, only: macheps34
    use variables, only: epsfeas, gmax, outunit, sa, scalea

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)  :: annz, m
    integer, intent(out) :: flag

    ! ARRAY ARGUMENTS
    integer,       dimension(annz), intent(in)    :: arow, acol
    real(kind=dp), dimension(annz), intent(inout) :: aval
    real(kind=dp), dimension(m),    intent(inout) :: b

    ! ----------------------------------------------------------------
    ! This subroutine scales the jacobian of the constraints (matrix
    ! A)
    ! 
    ! ALLOCATES : sa
    ! ----------------------------------------------------------------

    ! LOCAL SCALAR
    integer       :: allocstat, i
    real(kind=dp) :: minscale, norm

    flag = 0

    ! The minimum value for scaling the constraints
    minscale = macheps34 / epsfeas

    ! Allocates the sa array
    allocate( sa(m), stat=allocstat )
    if( allocstat .ne. 0 ) then
       write (      *, 100)
       write (outunit, 100)

       flag = - 10
       return
    end if

    if ( scalea ) then

       ! Computes the inf norm of each line of A
       sa(1:m) = 0.0d0
       do i = 1, annz
          sa(arow(i)) = max( abs(aval(i)), sa(arow(i)) )
       end do

       ! Computes the scale factors
       do i = 1, m
          if(sa(i) .ne. 0) then
             sa(i) = max( min( 1.0d0, gmax / sa(i) ), minscale )
          else
             sa(i) = 1.0d0
          end if
       end do

       ! Scale the constraints
       aval(1:annz) = sa(arow(1:annz)) * aval(1:annz)
       b(1:m)       = sa(1:m) * b(1:m)

    else
       sa(1:m) = 1.0d0
    end if

    ! NONEXECUTABLE STATEMENTS
100 format(/, 1X, "ERROR (scale_constraints): Unable to allocate memory.")

  end subroutine scale_constraints

  ! ------------------------------------------------------------------------

  subroutine check_constraints_rank(n, m, b, annz, arow, acol, aval, &
       annzrank, mred, flag)

    use constants, only: macheps, macheps34
    use hsl_ma48_double
    use hsl_zd11_double
    use testing
    use util,      only: swap
    use variables, only: cind, epsfeas, hslunit, maxhslfix, outunit
    
    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: m, n
    integer, intent(inout) :: annz
    integer, intent(out)   :: annzrank, flag, mred

    ! ARRAY ARGUMENTS
    integer,       dimension(annz), intent(inout) :: arow, acol
    real(kind=dp), dimension(annz), intent(inout) :: aval
    real(kind=dp), dimension(m),    intent(inout) :: b

    ! ----------------------------------------------------------------
    ! This subroutines verifies A for linear dependent lines and
    ! reorder its elements if necessary.
    !
    ! ALLOCATES : cind
    !
    ! CALLS : HSL_MA48, MC58ID, MC58AD, MC59AD
    ! ----------------------------------------------------------------

    ! LOCAL SCALARS

    integer       :: allocstat, cnt, i, la, lirn, liw, ljcn, nnz, &
                     rank, trans
    real(kind=dp) :: ptolA, ptolAb

    ! LOCAL ARRAYS
    integer,       allocatable, dimension(:)   :: ip, irn, iw, jcn
    integer,                    dimension(20)  :: icntl, info
    integer,                    dimension(n+1) :: cols
    integer,                    dimension(m)   :: newrow, rows
    real(kind=dp), allocatable, dimension(:)   :: a
    real(kind=dp),              dimension(m)   :: coef
    real(kind=dp),              dimension(n)   :: rhs
    real(kind=dp),              dimension(10)  :: cntl, rinfo

    ! LOCAL HSL_MA48 TYPES
    type(zd11_type)    :: msys
    type(ma48_control) :: control
    type(ma48_ainfo)   :: ainfo
    type(ma48_finfo)   :: finfo
    type(ma48_sinfo)   :: sinfo
    type(ma48_factors) :: factors

    ! EXTERNAL SUBROUTINES
    external :: mc58id, mc58ad, mc59ad

    flag = 0

    if( annz .eq. 0 ) then
       mred = 0

       allocate( cind(mred), stat=allocstat )
       if ( allocstat .ne. 0 ) then
          write (      *, 100)
          write (outunit, 100)

          flag = - 10
          return
       end if

       return
    end if

    ! Open output stream for HSL messages
    open(hslunit, file='mc58.out')

    la   = 4 * (annz + m)
    lirn = 4 * (annz + m)
    ljcn = 3 * (annz + m)
    liw  = 6 * (m+n+1) + max( m, n+1 )

    allocate( irn(lirn), jcn(ljcn), a(la), iw(liw), &
         stat=allocstat )

    if(allocstat .ne. 0) then
       write (      *, 100)
       write (outunit, 100)

       flag = - 10
       return
    end if

    nnz = annz
    irn(1:nnz) = arow(1:annz)
    jcn(1:nnz) = acol(1:annz)
    a(1:nnz)   = aval(1:annz)

    ! Sets controlling parameters
    call mc58id(cntl, icntl)
    icntl(1) = hslunit
    icntl(2) = hslunit
    icntl(3) = hslunit
    icntl(4) = 2

    ! cntl(3) = epsfeas

    call smc58ad( m, n, rows )
    if(flag .lt. 0) go to 999

    ptolA = cntl(3) * rinfo(2)

    annzrank = 0
    mred = rank

    if ( rank .lt. m ) then
       write (      *, 110)
       write (outunit, 110)

       ! Now double check to see if Ax = b is compatible, by computing
       ! the rank of the matrix (A b)
       nnz = annz
       irn(1:nnz) = arow(1:annz)
       jcn(1:nnz) = acol(1:annz)
       a(1:nnz)   = aval(1:annz)

       do i = 1, m
          if ( abs( b(i) ) .gt. macheps ) then
             nnz = nnz + 1
             irn(nnz) = i
             jcn(nnz) = n+1
             a(nnz)   = b(i)
          end if
       end do

       ! Check rank of the matrix ( A | b )
       call smc58ad( m, n+1, newrow )
       if ( flag .lt. 0 ) go to 999

       ptolAb = cntl(3) * rinfo(2)

       ! If the rank of the augmented matrix is different of the rank
       ! of A, then the system Ax = b is incompatible
       if ( rank .ne. mred ) then
          flag = - 6
       end if

       ! Copies b and lambda in a
       a(annz+1:annz+m) = b(1:m)

       ! Allocates cind
       allocate( cind(mred), stat=allocstat )
       if ( allocstat .ne. 0 ) then
          write (      *, 100)
          write (outunit, 100)

          flag = - 10
          return
       end if

       ! It will hold the new values for lines of A 
       ! Also uses the loop to reorder elements from b
       do i = 1, mred
          newrow(rows(i)) = i
          cind(i) = rows(i)
          b(i) = a(annz+rows(i))
       end do

       do i = mred+1, m
          newrow(rows(i)) = mred-i
          b(i) = a(annz+rows(i))
       end do

       ! Reorder elements of A
       cnt = 1
       annzrank = 0
       do while ( cnt .le. annz )
          if ( newrow(arow(cnt)) .lt. 0 ) then
             arow(cnt) = - newrow(arow(cnt))
             call swap(arow(cnt), arow(annz))
             call swap(acol(cnt), acol(annz))
             call swap(aval(cnt), aval(annz))
             annz = annz - 1
             annzrank = annzrank + 1
          else
             arow(cnt) = newrow(arow(cnt))
             cnt = cnt + 1
         end if
       end do

       ! If the rank estimator detected different ranks between A and
       ! (A b), try another strategy to check Ax=b compatibility
       if ( flag .eq. - 6 .and. annzrank .gt. 0 ) then

          ! TODO: place MA48 in a module in order to make the presence
          ! of MA48 optional

          allocate( ip(m-mred+1), stat=allocstat )
          if( allocstat .ne. 0 ) then
             write (      *, 100)
             write (outunit, 100)

             flag = - 10
             return
          end if

          ! Sort the last elements of A by rows
          if( hslunit .gt. 0 ) open( hslunit, file='mc59.out' )

          irn(1:annzrank) = arow(annz+1:annz+annzrank)
          jcn(1:annzrank) = acol(annz+1:annz+annzrank)
          a(1:annzrank)   = aval(annz+1:annz+annzrank)

          icntl(1) = 0
          icntl(2) = 0
          icntl(3) = 0
          icntl(4) = hslunit
          icntl(5) = hslunit
          icntl(6) = 0

          call mc59ad( icntl(1:10), m-mred, n, annzrank, jcn,     &
               annzrank, irn, annzrank, a, m-mred+1, ip, liw, iw, &
               info(1:10) )

          if ( hslunit .gt. 0 ) close( hslunit )

          if ( info(1) .lt. 0 ) then
             write (      *, 200) info(1)
             write (outunit, 200) info(1)

             flag = - 20
             return
          end if

          allocate( msys%col(annz), msys%row(annz), msys%val(annz), &
               stat=allocstat )
          if ( allocstat .ne. 0 ) then
             write (      *, 100)
             write (outunit, 100)

             flag = - 10
             return
          end if

          msys%m  = mred
          msys%n  = n
          msys%ne = annz

          msys%row(1:annz) = arow(1:annz)
          msys%col(1:annz) = acol(1:annz)
          msys%val(1:annz) = aval(1:annz)

          ! Initialize
          call ma48_initialize( factors, control )

          ! Analyse and factorsorize
          call ma48_analyse( msys, factors, control, ainfo, &
               finfo )

          if ( ainfo%flag .lt. 0 ) then
             write (      *, 300) ainfo%flag
             write (outunit, 300) ainfo%flag

             flag = - 9
             return
          else if ( finfo%flag .lt. 0 ) then
             write (      *, 300) finfo%flag
             write (outunit, 300) finfo%flag

             flag = - 9
             return
          end if

          ! Solve the annzrank systems
          do i = 1, m-mred
             ! The rhs will be the i-th eliminated line from A
             rhs(1:n) = 0.0d0
             rhs(jcn(ip(i):ip(i+1)-1)) = a(ip(i):ip(i+1)-1)

             ! Solve A^T coef = rhs
             call ma48_solve( msys, factors, rhs, coef(1:mred), &
                  control, sinfo, trans )

             if ( abs( b(mred+i) - dot_product( b(1:mred), coef(1:mred) ) ) .gt. macheps34 ) then
                flag = - 6
                return
             end if
          end do

          flag = 0

          write (      *, 115) mred, ptolA, rank, ptolAb
          write (outunit, 115) mred, ptolA, rank, ptolAb

          call ma48_finalize( factors, control, allocstat )

          deallocate( ip, msys%col, msys%row, msys%val )
       end if

       write (      *, 120) m-mred
       write (outunit, 120) m-mred

    else

       ! Allocates cind
       allocate( cind(m), stat=allocstat )
       if ( allocstat .ne. 0 ) then
          write (      *, 100)
          write (outunit, 100)

          flag = - 10
          return
       end if

       cind(1:m) = (/ ( i, i = 1, m ) /)

    end if

999 deallocate( irn, jcn, a, iw )

    close(hslunit)

    ! NONEXECUTABLE STATEMENTS
100 format(/, 1X, "ERROR(check_constraints): Unable to allocate memory.")

110 format(/, 1X, "WARNING: The matrix A does not have full rank." &
           /, 1X, "         The preprocessor will fix it.")

115 format(/, 1X, "WARNING: The rank estimator reports that Ax = b is incompatible.",   &
           /, 1X, "         MC59 reports:",                                             &
           /, 1X, "              RANK( A )   = ", I7, ", PIVOT TOLERANCE = ", 1P, D7.1, &
           /, 1X, "              RANK( A b ) = ", I7, ", PIVOT TOLERANCE = ", 1P, D7.1, &
           /, 1X, "         Nevertheless, it seems the rank is overestimated."          &
           /, 1X, "         LCMIN will continue.")

120 format(/, 1X, "The linear dependent lines of A were removed.", &
           /, 1X, "Number of redundant constraints : ", I7)

200 format(/, 1X, "ERROR (check_constraints_rank): MC59AD executed with errors. It reports ", I3, ".")

300 format(/, 1X, "ERROR (check_constraints_rank): MA48 executed with errors. It reports ", I3, ".")

  contains

    subroutine smc58ad( m, n, rows )

      ! SCALAR ARGUMENTS
      integer, intent(in) :: m, n

      ! ARRAY ARGUMENT
      integer, dimension(m), intent(out) :: rows

      ! TMP
      integer :: i
      
      ! Verifies the linear independent rows and columns (we are
      ! interested only in rows)
      call mc58ad( m, n, nnz, la, a, lirn, irn, ljcn, jcn, cntl, icntl, &
           liw, iw, info, rinfo, rank, rows, cols )

      ! TESTING STATEMENT
      ! Print the indexes of the lines of a linear independent submatrix
      if ( rank .ne. m ) then
         open ( 123, file='li_rows.txt' )
         write ( 123, * ) m, rank
         do i = 1, rank
            write ( 123, * ) rows(i)
         end do
         close ( 123 )
      end if

      ! Error handling for mc58ad
      cnt = 0
      do
         if(info(1) >= 0) then
            exit
         end if

         if(cnt > maxhslfix) then
            write(*, 101) info(1)
            write(outunit, 101) info(1)

            flag = - 9
            return
         end if

         select case (info(1))

         case (-4) ! LA is too small
            ! Reallocates the vector
            deallocate(a, irn, jcn)
            la   = max( info(4), la ) + annz
            lirn = max( info(4), la ) + annz
            ljcn = max( info(4), la ) + annz
            allocate(a(la), irn(lirn), jcn(ljcn), stat=allocstat)
            if(allocstat /= 0) then
               write (*, 100) la
               flag = - 10
               return
            end if
            
         case (-5) ! LIRN is too small
            ! Reallocates the vector
            deallocate(irn)
            lirn = info(5) + annz
            allocate(irn(lirn), stat=allocstat)
            if(allocstat /= 0) then
               write (*, 100) lirn
               flag = - 10
               return
            end if

         case (-6) ! LJCN is too small
            ! Reallocates the vector
            deallocate(jcn)
            ljcn = info(6) + annz
            allocate(jcn(ljcn), stat=allocstat)
            if(allocstat /= 0) then
               write (*, 100) ljcn
               flag = - 10
               return
            end if

         case (-7) ! LIW is too small
            ! Reallocates the vector
            deallocate(iw)
            liw = info(7)
            allocate(iw(liw), stat=allocstat)
            if(allocstat /= 0) then
               write (*, 100) liw
               flag = - 10
               return
            end if
            
         end select

         ! Restores default values of matrix A
         irn(1:annz) = arow
         jcn(1:annz) = acol
         a(1:annz) = aval

         ! Call MC58AD again
         call mc58ad(m, n, nnz, la, a, lirn, irn, ljcn, jcn, cntl, icntl, &
              liw, iw, info, rinfo, rank, rows, cols)

         cnt = cnt + 1
      end do

      ! NONEXECUTABLE STATEMENTS
101   format(/, 1X, "ERROR (smc58ad): MC58AD executed with errors. It reports ", I3, ".")

100   format(/, 1X, "ERROR (smc58ad): Unable to allocate memory.", &
             /, 1X, "                 We are trying to allocate a vector with length ", I10, ".", &
             /, 1X, "                 There may be problems in the execution of MC58AD routine.", & 
             /, 1X, "                 Check mc58.out for details.")

    end subroutine smc58ad

  end subroutine check_constraints_rank

  ! ------------------------------------------------------------------------

  subroutine check_constraints_fxv(m, b, n, x, l, u, annz, arow, acol, &
       aval, nred, flag)

    use constants, only: macheps
    use util,      only: swap
    use variables, only: find, nfind, outunit

    ! SCALAR ARGUMENT
    integer, intent(in)    :: m, n
    integer, intent(inout) :: annz
    integer, intent(out)   :: flag, nred
    
    ! ARRAY ARGUMENTS
    integer,       dimension(annz), intent(inout) :: arow, acol
    real(kind=dp), dimension(annz), intent(inout) :: aval
    real(kind=dp), dimension(m),    intent(inout) :: b
    real(kind=dp), dimension(n),    intent(inout) :: x, l, u

    ! ----------------------------------------------------------------
    ! This subroutine checks the bound constraints for fixed
    ! variables. Fixed variables are eliminated from the problem, so a
    ! rearrangement of the jacobian of the equality constraints is
    ! also required. Arrays containing the index of all fixed and
    ! non-fixed variables are stored.
    !
    ! ALLOCATES : find, nfind
    ! ----------------------------------------------------------------

    ! LOCAL SCALARS
    integer :: allocstat, annzfv, i

    ! LOCAL ARRAY
    real(kind=dp), dimension(n) :: work

    flag = 0

    ! Counts the number of non-fixed variables
    nred = count( (u(1:n) - l(1:n)) .gt. macheps )

    allocate( find(n), nfind(nred), stat=allocstat )
    if( allocstat .ne. 0 ) then
       write(      *, 300)
       write(outunit, 300)

       flag = - 10
       return
    end if

    ! nfind(j) holds the original index of the non-fixed variable x_j
    !  find(i) holds the new index of the x_i, if x_i is a non-fixed
    !          variable, or 0 otherwise
    nred = 0
    do i = 1, n
       if( u(i) - l(i) .gt. macheps ) then
          nred = nred + 1
          nfind(nred) = i
          find(i) = nred
       else
          x(i) = l(i)
          find(i) = 0
       end if
    end do

    annzfv = 0
    if( n .gt. nred ) then

       ! Reorder the jacobian of the constraints due to fixed variables
       ! elimination
       i = 1
       do while( i .le. annz )
          if( find(acol(i)) .eq. 0 ) then
             b(arow(i)) = b(arow(i)) - aval(i) * l(acol(i))
             call swap(arow(i), arow(annz))
             call swap(acol(i), acol(annz))
             call swap(aval(i), aval(annz))
             annz   = annz   - 1
             annzfv = annzfv + 1
          else
             acol(i) = find(acol(i))
             i = i + 1
          end if
       end do

       write (*, 100) n-nred
       write (outunit, 100) n-nred
    else
       write (*, 150)
       write (outunit, 150)
    end if

    ! NONEXECUTABLE STATEMENTS
100 format(/, 1X, "The fixed variables were removed.", &
           /, 1X, "Number of fixed variables     : ", I7)

110 format(   1X, "Number of removed constraints : ", I7)

150 format(/, 1X, "There are no fixed variables.")

300 format(/, 1X, "ERROR (check_constraints_fxv): Unable to allocate memory.")

  end subroutine check_constraints_fxv
  
  ! ------------------------------------------------------------------------

  subroutine check_constraints_frv(n, l, u, flag)

    use variables, only: infty, lind, lnu, nl, nu, outunit, uind, unl

    ! SCALAR ARGUMENTS
    integer, intent(in)  :: n
    integer, intent(out) :: flag

    ! ARRAY ARGUMENTS
    real(kind=dp), dimension(n), intent(in) :: l, u

    ! ----------------------------------------------------------------
    ! This subroutine checks bound constraints for free
    ! variables. Free variables, just require a rearrangement of the
    ! lower and upper bound arrays. The index of non-free variables
    ! from below and above are stored.
    !
    ! ALLOCATES : lind, lnu, uind, unl
    ! ----------------------------------------------------------------

    ! LOCAL SCALARS
    integer :: allocstat, i, nlnu, nunl

    flag = 0

    ! Checks lower and upper bounds for free variables.
    ! lind(i) = j if x_j is a non-free variable from below.
    ! uind(i) = j if x_j is a non-free variable from above.
    nl   = 0
    nu   = 0
    nlnu = 0
    nunl = 0
    do i = 1, n
       if( l(i) .gt. -infty ) then
          nl = nl + 1
          if( u(i) .ge. infty ) nlnu = nlnu + 1
       end if

       if( u(i) .lt. infty ) then
          nu = nu + 1
          if( l(i) .le. -infty ) nunl = nunl + 1
       end if
    end do

    allocate(lind(nl), uind(nu), lnu(nlnu), unl(nunl), &
         stat=allocstat)
    if(allocstat .ne. 0) then
       flag = - 10
       return
    end if

    nl   = 0
    nu   = 0
    nlnu = 0
    nunl = 0
    do i = 1, n
       if( l(i) .gt. -infty ) then
          nl = nl + 1
          lind(nl) = i

          if( u(i) .ge. infty ) then
             nlnu = nlnu + 1
             lnu(nlnu) = i
          end if
       end if

       if( u(i) .lt. infty ) then
          nu = nu + 1
          uind(nu) = i

          if( l(i) .le. -infty ) then
             nunl = nunl + 1
             unl(nunl) = i
          end if
       end if
    end do

    write (*, 100) nl, nu
    write (outunit, 100) nl, nu

    ! NONEXECUTABLE STATEMENTS
100 format(/, 1X, "Number of lower bounded variables: ", I7, &
           /, 1X, "Number of upper bounded variables: ", I7)

110 format(/, 1X, "ERROR (check_constraints_frv): Unable to allocate memory.")

  end subroutine check_constraints_frv

end module constraints
