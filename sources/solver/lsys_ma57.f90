module lsys_solver

  use constants, only: dp
  use hsl_ma57_double
  use variables, only: hslunit, outunit

  implicit none

  save

  private

  public :: lsys_analyse, lsys_factorize, lsys_initialize, lsys_solve

  ! LOCAL STRUCTURES
  type(zd11_type)    :: matrix
  type(ma57_control) :: cntl
  type(ma57_factors) :: fact
  type(ma57_ainfo)   :: ainfo
  type(ma57_finfo)   :: finfo
  type(ma57_sinfo)   :: sinfo

contains

  subroutine lsys_initialize()

    ! Open output stream for HSL messages
    open(hslunit, file='ma57.out')

    ! Initialize default controlling parameters
    call ma57_initialize(fact, cntl)
    cntl%lp = hslunit
    cntl%mp = hslunit
    cntl%wp = hslunit

  end subroutine lsys_initialize

  ! --------------------------------------------------------------------

  subroutine lsys_analyse(nsys, mnnz, mrow, mcol, mval, flag)

    ! SCALAR ARGUMENTS
    integer, intent(in)  :: mnnz, nsys
    integer, intent(out) :: flag

    ! ARRAY ARGUMENTS
    integer,       dimension(mnnz), intent(in) :: mrow, mcol
    real(kind=dp), dimension(mnnz), intent(in) :: mval

    ! LOCAL SCALARS
    integer :: allocstat

    flag = 0

    ! Deallocate arrays if allocated
    if( allocated( matrix%row ) ) then
       deallocate( matrix%col, matrix%row, matrix%val, stat=allocstat )
    end if

    ! Allocate arrays
    allocate(matrix%row(mnnz), matrix%col(mnnz), matrix%val(mnnz), &
         stat=allocstat)
    if(allocstat .ne. 0) then
       write (*, 100)
       write (outunit, 100)

       flag = - 10
       return
    end if

    ! Initialize the system matrix data
    matrix%n   = nsys
    matrix%ne  = mnnz

    matrix%row(1:mnnz) = mrow(1:mnnz)
    matrix%col(1:mnnz) = mcol(1:mnnz)
    matrix%val(1:mnnz) = mval(1:mnnz)

    ! Analyse
    call ma57_analyse(matrix, fact, cntl, ainfo)
    if(ainfo%flag .lt. 0) then
       write (      *, 120)
       write (outunit, 120)

       flag = - 9
       return
    end if

    ! NON EXECUTABLE STATEMENTS
100 format(/, 1X, "ERROR (lsys_analyse): Unable to allocate memory.")

120 format(/, 1X, "ERROR (lsys_analyse): MA57 error. Check ma57.out for details.")

  end subroutine lsys_analyse

  ! --------------------------------------------------------------------

  subroutine lsys_factorize(nsys, nneigv, rank, flag, n, diag)

    ! SCALAR ARGUMENTS
    integer,       intent(in)            :: nsys
    integer,       intent(in),  optional :: n
    integer,       intent(out)           :: flag, nneigv, rank
    real(kind=dp), intent(in),  optional :: diag

    flag = 0

    ! If necessary, updates the diagonal
    if( present( diag ) .and. present( n ) ) then
       matrix%val(1:n) = matrix%val(1:n) + diag
    end if

    ! Factorize
    call ma57_factorize(matrix, fact, cntl, finfo)
    if(finfo%flag .lt. 0) then
       write (      *, 100)
       write (outunit, 100)

       flag = - 9
       return
    end if

    rank   = finfo%rank
    nneigv = finfo%neig

    ! NON EXECUTABLE STATEMENTS
100 format(/, 1X, "ERROR (lsys_factorize): MA57 error. Check ma57.out for details.")

  end subroutine lsys_factorize

  ! --------------------------------------------------------------------

  subroutine lsys_solve(nsys, rhs, flag)

    ! SCALAR ARGUMENTS
    integer, intent(in)  :: nsys
    integer, intent(out) :: flag

    ! ARRAY ARGUMENTS
    real(kind=dp), dimension(nsys), intent(inout) :: rhs

    ! ----------------------------------------------------------------
    ! This subroutine solves a symmetric linear system using MA57. It
    ! also contains a inertia corrector, used exclusively by the KKT
    ! system, setting 'inertia' argument as a positive integer equals
    ! m (the correct value for info(24)).
    ! ----------------------------------------------------------------

    ! LOCAL SCALARS
    integer       :: allocstat, i, m, n
    real(kind=dp) :: chi, chiold

    ! LOCAL ARRAYS
    real(kind=dp), dimension(nsys) :: sol

    ! Solve the linear system
    sol(1:nsys) = rhs(1:nsys)
    call ma57_solve(matrix, fact, sol, cntl, sinfo)
    if(sinfo%flag .lt. 0) then
       write (      *, 120)
       write (outunit, 120)

       flag = - 9
       go to 999
    end if

    ! Perform iterative refinement
    call ma57_solve(matrix, fact, sol, cntl, sinfo, rhs=rhs, iter=1)
    if(sinfo%flag .lt. 0) then
       write (      *, 120)
       write (outunit, 120)

       flag = - 9
       go to 999
    end if
       
    ! Return the solution on rhs
    rhs(1:nsys) = sol(1:nsys)

999 continue

    ! Finalize MA57 structures
    call ma57_finalize( fact, cntl, allocstat )

    close ( hslunit )

    deallocate( matrix%col, matrix%row, matrix%val )

    ! NONEXECUTABLE STATEMENTS
100 format(/, 1X, "ERROR (lsys_solve): Unable to allocate memory.")

120 format(/, 1X, "ERROR (lsys_solve): MA57 error. Check ma57.out for details.")

  end subroutine lsys_solve

end module lsys_solver
