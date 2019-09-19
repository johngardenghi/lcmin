module algsolver

  use constants, only: dp

  implicit none

  save

  ! GLOBAL SCALAR
  integer :: nnz

  ! GLOBAL POINTERS
  integer,       dimension(:), pointer :: irn, jcn
  real(kind=dp), dimension(:), pointer :: a, rhs

contains

  subroutine algsolve(n, x, l, u, m, b, annz, arow, acol, aval, &
       algfeas, innsv, flag)

    use variables, only: epstol, out_file, outunit

    ! SCALAR ARGUMENTS
    character(len=2), intent(in)  :: innsv
    integer,          intent(in)  :: annz, m, n
    integer,          intent(out) :: flag
    real(kind=dp),    intent(in)  :: algfeas

    ! ARRAY ARGUMENTS
    integer,       dimension(annz), intent(in)    :: arow, acol
    real(kind=dp), dimension(annz), intent(in)    :: aval
    real(kind=dp), dimension(m),    intent(in)    :: b
    real(kind=dp), dimension(n),    intent(in)    :: l, u
    real(kind=dp), dimension(n),    intent(inout) :: x

    ! ------------------------------------------------------------
    ! This subroutine solves the feasibility problem
    !
    !                min | Ax - b |
    !                s.a l <= x <= u
    !
    ! in order to find a suitable initial point to LCMIN using
    ! ALGENCAN.
    ! ------------------------------------------------------------

    ! LOCAL SCALARS
    character(len=80) :: specfnm, outputfnm
    logical           :: checkder
    integer           :: allocstat, inform, nvparam
    real(kind=dp)     :: cnorm, efacc, efstain, eoacc, eostain, &
                         f, nlpsupn, snorm

    ! LOCAL ARRAYS
    character(len=80), dimension(4)  :: vparam
    logical,           dimension(11) :: coded(11)
    logical,           dimension(m)  :: equatn, linear
    real(kind=dp),     dimension(m)  :: lambda

    ! Parameters setting
    equatn    = .true.
    linear    = .true.

    checkder  = .false.
    
    efstain   = sqrt( algfeas )
    eostain   = epstol ** 1.5d0
    
    efacc     = sqrt( algfeas )
    eoacc     = sqrt( epstol )
    
    specfnm   = ''
    outputfnm = out_file
    
    nvparam   = 3
    vparam(1) = 'IGNORE-OBJECTIVE-FUNCTION'
    vparam(2) = 'ITERATIONS-OUTPUT-DETAIL 11'
    ! vparam(3) = 'OBJECTIVE-AND-CONSTRAINTS-SCALING-AVOIDED'

    if(innsv .eq. 'tn') then
       vparam(3) = 'TRUNCATED-NEWTON-LINE-SEARCH-INNER-SOLVER'
    elseif(innsv .eq. 'nw') then
       vparam(3) = 'NEWTON-LINE-SEARCH-INNER-SOLVER MA57 MC64'
    else
       nvparam = 2
    end if

    ! Set problem data
    lambda = 0.0d0

    ! Coded subroutines
    coded( 1) = .false. ! evalf
    coded( 2) = .false. ! evalg
    coded( 3) = .false. ! evalh
    coded( 4) = .false. ! evalc
    coded( 5) = .false. ! evaljac
    coded( 6) = .false. ! evalhc
    coded( 7) = .true.  ! evalfc
    coded( 8) = .true.  ! evalgjac
    coded( 9) = .false. ! evalgjacp
    coded(10) = .true.  ! evalhl
    coded(11) = .false. ! evalhlp

    ! Initialize problem data
    call initalg( m, b, annz, arow, acol, aval )

    close( outunit )

    ! Call Algencan solver
    call algencan( algevalf, algevalg, algevalh, algevalc, algevaljac, algevalhc, &
         algevalfc, algevalgjac, algevalgjacp, algevalhl, algevalhlp, annz, 0,    &
         algfeas, epstol, efstain, eostain, efacc, eoacc, outputfnm, specfnm,     &
         nvparam, vparam, n, x, l, u, m, lambda, equatn, linear, coded, checkder, &
         f, cnorm, snorm, nlpsupn, flag )

    open( outunit, file=out_file, status='old', access='append' )

    write (      *, 200)
    write (outunit, 200)

    nullify( a, irn, jcn, rhs )

    ! NONEXECUTABLE STATEMENTS
200 format(/, 1X, "==============================================================================")

  end subroutine algsolve 

  ! ******************************************************************
  ! ******************************************************************
   
  subroutine initalg(m, b, annz, arow, acol, aval)

    ! SCALAR ARGUMENTS
    integer, intent(in) :: annz, m
    
    ! ARRAY ARGUMENTS
    integer,       dimension(annz), intent(in), target :: arow, acol
    real(kind=dp), dimension(annz), intent(in), target :: aval
    real(kind=dp), dimension(m),    intent(in), target :: b

    nnz = annz

    rhs => b
    irn => arow
    jcn => acol
    a   => aval

  end subroutine initalg

  ! ******************************************************************
  ! ******************************************************************
   
  subroutine algevalf(n,x,f,flag)

    ! SCALAR ARGUMENTS
    integer,      intent(in)  :: n
    integer,      intent(out) :: flag
    real(kind=8), intent(out) :: f

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    flag = - 1

  end subroutine algevalf

  ! ******************************************************************
  ! ******************************************************************

  subroutine algevalg(n,x,g,flag)

    ! SCALAR ARGUMENTS
    integer,      intent(in)  :: n
    integer,      intent(out) :: flag

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: g(n)

    flag = - 1

  end subroutine algevalg

  ! ******************************************************************
  ! ******************************************************************

  subroutine algevalh(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)

    ! SCALAR ARGUMENTS
    logical,      intent(out) :: lmem
    integer,      intent(in)  :: lim,n
    integer,      intent(out) :: flag,hnnz

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: hcol(lim),hrow(lim)
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: hval(lim)

    flag = - 1

  end subroutine algevalh

  ! ******************************************************************
  ! ******************************************************************

  subroutine algevalc(n,x,ind,c,flag)

    ! SCALAR ARGUMENTS
    integer,      intent(in)  :: ind,n
    integer,      intent(out) :: flag
    real(kind=8), intent(out) :: c

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: x(n)

    flag = - 1

  end subroutine algevalc

  ! ******************************************************************
  ! ******************************************************************

  subroutine algevaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

    ! SCALAR ARGUMENTS
    logical, intent(out) :: lmem
    integer, intent(in)  :: ind,lim,n
    integer, intent(out) :: flag,jcnnz

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: jcvar(lim)
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: jcval(lim)

    flag = - 1

  end subroutine algevaljac

  ! ******************************************************************
  ! ******************************************************************

  subroutine algevalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

    ! SCALAR ARGUMENTS
    logical, intent(out) :: lmem
    integer, intent(in)  :: ind,lim,n
    integer, intent(out) :: flag,hcnnz

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: hccol(lim),hcrow(lim)
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: hcval(lim)

    flag = - 1

  end subroutine algevalhc

  ! ******************************************************************
  ! ******************************************************************

  subroutine algevalfc(n,x,f,m,c,flag)

    use util, only: smvp

    ! SCALAR ARGUMENTS
    integer,      intent(in)  :: m,n
    integer,      intent(out) :: flag
    real(kind=8), intent(out) :: f

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: c(m)

    integer :: i

    flag = 0

    f = 0.0d0

    c = smvp(m, nnz, irn, jcn, a, x) - rhs

  end subroutine algevalfc

  ! ******************************************************************
  ! ******************************************************************

  subroutine algevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem,flag)

    ! SCALAR ARGUMENTS
    logical,      intent(out) :: lmem
    integer,      intent(in)  :: lim,m,n
    integer,      intent(out) :: flag,jcnnz

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: jcfun(lim),jcvar(lim)
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: g(n),jcval(lim)

    flag = 0
    lmem = .false.

    if(nnz .gt. lim) then
       lmem = .true.
       return
    end if

    g = 0.0d0

    jcfun(1:nnz) = irn
    jcvar(1:nnz) = jcn
    jcval(1:nnz) = a

    jcnnz = nnz

  end subroutine algevalgjac

  ! ******************************************************************
  ! ******************************************************************

  subroutine algevalgjacp(n,x,g,m,p,q,work,gotj,flag)

    ! SCALAR ARGUMENTS
    logical,   intent(inout) :: gotj
    integer,   intent(in)    :: m,n
    integer,   intent(out)   :: flag
    character, intent(in)    :: work

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)    :: x(n)
    real(kind=8), intent(inout) :: p(m),q(n)
    real(kind=8), intent(out)   :: g(n)

    flag = - 1

  end subroutine algevalgjacp

  ! ******************************************************************
  ! ******************************************************************

  subroutine algevalhl(n,x,m,lambda,sf,sc,hlrow,hlcol,hlval,hlnnz,lim,lmem,flag)

    ! SCALAR ARGUMENTS
    logical,      intent(out) :: lmem
    integer,      intent(in)  :: lim,m,n
    integer,      intent(out) :: flag,hlnnz
    real(kind=8), intent(in)  :: sf

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: hlcol(lim),hlrow(lim)
    real(kind=8), intent(in)  :: lambda(m),sc(m),x(n)
    real(kind=8), intent(out) :: hlval(lim)

    flag = 0
    lmem = .false.

    hlnnz = 0

  end subroutine algevalhl

  ! ******************************************************************
  ! ******************************************************************

  subroutine algevalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

    ! SCALAR ARGUMENTS
    logical,      intent(inout) :: goth
    integer,      intent(in)    :: m,n
    integer,      intent(out)   :: flag
    real(kind=8), intent(in)    :: sf

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: lambda(m),p(n),sc(m),x(n)
    real(kind=8), intent(out) :: hp(n)

    flag = - 1

  end subroutine algevalhlp

end module algsolver
