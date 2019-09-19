module user_subs

  use constants, only: dp

  implicit none

  ! Interfaces to user subroutines
  abstract interface
     subroutine evalf(n, x, f, flag)
       import :: dp
       ! SCALAR ARGUMENTS
       integer,       intent(in)  :: n
       integer,       intent(out) :: flag
       real(kind=dp), intent(out) :: f
       ! ARRAY ARGUMENTS
       real(kind=dp), dimension(n), intent(in) :: x
     end subroutine evalf
     
     subroutine evalg(n, x, g, flag)
       import :: dp
       ! SCALAR ARGUMENTS
       integer, intent(in)  :: n
       integer, intent(out) :: flag
       ! ARRAY ARGUMENTS
       real(kind=dp), dimension(n), intent(in)  :: x
       real(kind=dp), dimension(n), intent(out) :: g
     end subroutine evalg
     
     subroutine evalh(n, x, hnnzmax, hrow, hcol, hval, hnnz, flag)
       import :: dp
       ! SCALAR ARGUMENTS
       integer, intent(in)  :: n, hnnzmax
       integer, intent(out) :: flag, hnnz
       ! ARRAY ARGUMENTS
       integer,       dimension(hnnzmax), intent(out) :: hrow, hcol
       real(kind=dp), dimension(n),       intent(in)  :: x
       real(kind=dp), dimension(hnnzmax), intent(out) :: hval
     end subroutine evalh
  end interface

  ! Pointers to user subroutines
  procedure(evalf), pointer :: subf
  procedure(evalg), pointer :: subg
  procedure(evalh), pointer :: subh

end module user_subs
