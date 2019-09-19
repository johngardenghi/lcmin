module util

  use constants, only: dp

  implicit none

  public :: infnorm, norm1, smvp, swap

  private

  interface infnorm
     module procedure infnorm_matrix, infnorm_vector
  end interface infnorm

  interface swap
     module procedure swap_dp, swap_int
  end interface swap

contains

  function smvp(lines, mnnz, mrow, mcol, mval, vec, symm)

    ! FUNCTION TYPE
    real(kind=dp), dimension(lines) :: smvp

    ! SCALAR ARGUMENT
    integer, intent(in)           :: mnnz, lines
    logical, intent(in), optional :: symm

    ! ARRAY ARGUMENTS
    integer,       dimension(mnnz), intent(in) :: mrow, mcol
    real(kind=dp), dimension(mnnz), intent(in) :: mval
    real(kind=dp), dimension(:),    intent(in) :: vec

    ! -------------------------------------------------------------
    ! This function returns the result of the product of a sparse
    ! matrix by a vector of dimension lines
    ! -------------------------------------------------------------

    ! LOCAL SCALAR
    integer :: i

    smvp(1:lines) = 0.0d0

    do i = 1, mnnz
       smvp(mrow(i)) = smvp(mrow(i)) + ( mval(i) * vec(mcol(i)) )

       if( present( symm ) ) then
          if( mcol(i) .ne. mrow(i) ) then
             smvp(mcol(i)) = smvp(mcol(i)) + ( mval(i) * vec(mrow(i)) )
          end if
       end if
    end do

  end function smvp

  ! --------------------------------------------------------------------

  subroutine swap_int(v1, v2)
    
    ! SCALAR ARGUMENTS
    integer, intent(inout) :: v1, v2

    ! ------------------------------------------------------------------
    ! This subroutine swaps the integer values from variables v1 and
    ! v2.
    ! ------------------------------------------------------------------

    ! LOCAL SCALAR
    integer :: temp

    temp = v1
    v1   = v2
    v2   = temp

  end subroutine swap_int

  ! --------------------------------------------------------------------

  subroutine swap_dp(v1, v2)
    
    ! SCALAR ARGUMENTS
    real(kind=dp), intent(inout) :: v1, v2

    ! ------------------------------------------------------------------
    ! This subroutine swaps the double precision values from variables
    ! v1 and v2.
    ! ------------------------------------------------------------------

    ! LOCAL SCALAR
    real(kind=dp) :: temp

    temp = v1
    v1   = v2
    v2   = temp

  end subroutine swap_dp

  ! --------------------------------------------------------------------

  function norm1(vec)

    ! FUNCTION TYPE
    real(kind=dp) :: norm1

    ! ARRAY ARGUMENTS
    real(kind=dp), dimension(:), intent(in) :: vec

    ! ---------------------------------------------------
    ! This function computes the 1-norm of a vector VEC.
    ! ---------------------------------------------------
    
    if( size(vec) .gt. 0 ) then
       norm1 = sum( abs(vec) )
    else
       norm1 = 0.0d0
    end if

  end function norm1

  ! --------------------------------------------------------------------

  function infnorm_vector(vec)

    ! ARRAY ARGUMENTS
    real(kind=dp), dimension(:), intent(in) :: vec

    ! FUNCTION TYPE
    real(kind=dp) :: infnorm_vector

    ! ----------------------------------------------------------
    ! This function computes the infinity norm of a vector VEC.
    ! ----------------------------------------------------------

    if( size(vec) .gt. 0 ) then
       infnorm_vector = maxval( abs(vec) )
    else
       infnorm_vector = 0.0d0
    end if

  end function infnorm_vector

  ! --------------------------------------------------------------------

  function infnorm_matrix(lines, mnnz, mrow, mval)

    ! SCALAR ARGUMENTS
    integer, intent(in) :: lines, mnnz

    ! ARRAY ARGUMENTS
    integer,       intent(in), dimension(mnnz) :: mrow
    real(kind=dp), intent(in), dimension(mnnz) :: mval

    ! FUNCTION TYPE
    real(kind=dp) :: infnorm_matrix

    ! --------------------------------------------------------
    ! This function computes the infinity norm of a matrix M.
    ! --------------------------------------------------------

    ! LOCAL SCALAR
    integer :: i

    ! LOCAL ARRAY
    real(kind=dp), dimension(lines) :: local

    local = 0.0d0
    do i = 1, mnnz
       local(mrow(i)) = local(mrow(i)) + abs(mval(i))
    end do

    infnorm_matrix = maxval(local)

  end function infnorm_matrix

end module util
