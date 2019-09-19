program lcminma

  implicit none
  
  ! LOCAL SCALARS
  integer :: n, m, annz, hnnzmax
  integer :: allocstat, flag
  
  ! LOCAL ARRAYS
  integer,      allocatable, dimension(:) :: alin, acol
  real(kind=8), allocatable, dimension(:) :: aval, b, l, lambda, u, x
  
  ! EXTERNAL FUNCIONS
  ! The problem codification must include:
  ! - evalf(n, x, f, flag) : evaluates the objective function 
  !                          at the given point x.
  !
  ! - evalg(n, x, g, flag) : evaluates the gradient of the objective 
  !                          function at the given point x.
  !
  ! - evalh(n, x, hnnzmax, hlin, 
  !         hcol, hval, hnnz, flag) : evaluates the hessian of the
  !                                   objective function at the given
  !                                   point x.
  !
  ! Notice that you can customize the name of these funcions.
  ! They just MUST have this prototype.
  ! Declare the names as external variables and pass as arguments to
  ! lcmin routine
  external :: evalf, evalg, evalh
  
  ! INITIALIZE THE SCALARS
  ! The problem codification must have a routine called INIDIM
  call inidim(n, m, annz, hnnzmax)
  
  allocate(x(n), l(n), u(n), b(m), lambda(m), alin(annz), acol(annz), &
       aval(annz), stat=allocstat)
  if(allocstat /= 0) then
     write (*,*) "ERROR: There was an unexpected problem with memory allocation."
     stop
  end if
  
  ! INITIALIZE PROBLEM DATA
  ! The problem codification must have a routine called INIP
  call inip(n, x, l, u, m, b, lambda, annz, alin, acol, aval, flag)
  
  ! INVOKE THE SOLVER MAIN ROUTINE
  if(flag .eq. 0) then
     call lcmin(n, x, l, u, m, b, lambda, alin, acol, aval, annz, &
          hnnzmax, 1.0d-08, 1.0d-09, evalf, evalg, evalh, flag)
  end if

  deallocate(x, l, u, b, lambda, alin, acol, aval)
    
end program lcminma
