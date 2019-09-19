module testing

  implicit none

  private

  public :: print_octave

  interface print_octave
     module procedure print_spmatrix, print_array
  end interface print_octave

contains

  subroutine print_array(m, vec, unit, file, flag)

    use constants, only: dp
    use variables, only: outunit

    ! SCALAR ARGUMENTS
    character(len=*), intent(in)  :: file
    integer,          intent(in)  :: m, unit
    integer,          intent(out) :: flag

    ! ARRAY ARGUMENT
    real(kind=dp), dimension(m), intent(in) :: vec

    ! ----------------------------------------------------------------
    ! This subroutine prints a sparse array vec in octave sparse
    ! format.
    ! ----------------------------------------------------------------

    ! LOCAL SCALAR
    integer :: i

    open ( unit, file=file )

    write (unit, *) "# name: ", file
    write (unit, *) "# type: sparse matrix"
    write (unit, *) "# nnz: ", m
    write (unit, *) "# rows: ", m
    write (unit, *) "# columns: ", 1

    write (unit, fmt='(I5, 1X, I5, 1X, 0P, F24.8)') ( i, 1, vec(i), i = 1, m )

    close ( unit )

  end subroutine print_array

  ! ------------------------------------------------------------------

  subroutine print_spmatrix(m, n, mnnz, mrow, mcol, mval, unit, file, &
       flag)

    use constants, only: dp, macheps
    use variables, only: outunit

    ! SCALAR ARGUMENTS
    character(len=*), intent(in)  :: file
    integer,          intent(in)  :: m, mnnz, n, unit
    integer,          intent(out) :: flag

    ! ARRAY ARGUMENTS
    integer,       dimension(mnnz), intent(in) :: mcol, mrow
    real(kind=dp), dimension(mnnz), intent(in) :: mval

    ! ----------------------------------------------------------------
    ! This subroutine prints a sparse matrix M in octave sparse
    ! format.
    !
    ! CALLS : MC59AD
    ! ----------------------------------------------------------------

    ! LOCAL SCALARS
    integer :: allocstat, i, j, liw

    ! LOCAL ARRAYS
    integer,                    dimension(10)   :: icntl, info
    integer,                    dimension(n+1)  :: ip
    integer,       allocatable, dimension(:)    :: iw
    integer,                    dimension(mnnz) :: col, row
    real(kind=dp),              dimension(mnnz) :: val

    ! EXTERNAL SUBROUTINES
    external :: mc59ad

    flag = 0

    liw = max(m,n) + 1
    allocate( iw(liw), stat=allocstat )
    if ( allocstat .ne. 0 ) then
       write (      *, 100)
       write (outunit, 100)

       flag = - 10
       return
    end if

    row(1:mnnz) = mrow(1:mnnz)
    col(1:mnnz) = mcol(1:mnnz)
    val(1:mnnz) = mval(1:mnnz)

    open ( unit, file=file )

    write (unit, *) "# name: ", file
    write (unit, *) "# type: sparse matrix"
    write (unit, *) "# nnz: ", mnnz
    write (unit, *) "# rows: ", m
    write (unit, *) "# columns: ", n

    if ( n .gt. 1 ) then
       icntl(1) = 0
       icntl(2) = 0
       icntl(3) = 0
       icntl(4) = 6
       icntl(5) = 6
       icntl(6) = 0

       call mc59ad( icntl, n, m, mnnz, row, mnnz, col, mnnz, val, n+1, &
            ip, liw, iw, info )

       if ( info(1) .lt. 0 ) then
          flag = - 9
          go to 999
       end if

       do j = 1, n
          do i = ip(j), ip(j+1)-1
             if ( abs( val(i) ) .gt. macheps ) then
                write (unit, fmt='(I5, 1X, I5, 1X, 0P, F24.16)') row(i), j, val(i)
             end if
          end do
       end do
    end if

    close (unit)

999 deallocate(iw)

    ! NONEXECUTABLE STATEMENTS
100 format(/, 1X, "ERROR (print_spmatrix): Unable to allocate memory.")

  end subroutine print_spmatrix

end module testing
