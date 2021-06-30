!     The routine kpsol is the preconditioner solve routine. It must have
!     that specific name be used in order that the c code can find and link
!     to it.  The argument list must also be as illustrated below:

subroutine fkpsol(udata, uscale, fdata, fscale, vv, ftem, ier)
use multicapa
use mkinsol
implicit none

integer neq, i
double precision udata(*), uscale(*), fdata(*), fscale(*)
double precision vv(*), ftem(*)
integer *4 ier ! Kinsol error flag

common /psize/ neq

neq = 3*ntot

do  i = 1, neq
   vv(i) = vv(i) * pp(i)
enddo
ier = 0

return
end

!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!* * *
!c     The routine kpreco is the preconditioner setup routine. It must have
!c     that specific name be used in order that the c code can find and link
!c     to it.  The argument list must also be as illustrated below:

subroutine fkpset(udata, uscale, fdata, fscale,vtemp1,vtemp2, ier)
use multicapa
use mkinsol
implicit none
integer *4 ier ! Kinsol error flag
integer neq, i
double precision udata(*), uscale(*), fdata(*), fscale(*)
double precision vtemp1(*), vtemp2(*)

common /psize/ neq

do i = 1, ntot
   pp(i) = 0.1  / (1.0+exp(-udata(i)))
enddo
do i = ntot+1, ntot*2
   pp(i) = 1.0 !0.1 / (1.0+exp(-udata(i)))
enddo

do i = 2*ntot+1, ntot*3
   pp(i) = 0.1 / (1.0+exp(-udata(i)))
enddo




   ier = 0
return
end

subroutine call_kinsol(x1_old, xg1_old, ier)
use multicapa
implicit none
integer *4 ier ! Kinsol error flag
integer i
real*8 x1(ntot*3), xg1(ntot*3)
real*8 x1_old(ntot*3), xg1_old(ntot*3)
integer*8 iout(15) ! Kinsol additional output information
real*8 rout(2) ! Kinsol additional out information
integer*8 msbpre
real*8 fnormtol, scsteptol
real*8 scale(ntot*3)
real*8 constr(ntot*3)
integer*4  globalstrat, maxl, maxlrst
integer neq ! Kinsol number of equations
integer*4 max_niter
common /psize/ neq ! Kinsol
integer ierr

neq=ntot*3

! INICIA KINSOL

msbpre  = 1 ! maximum number of iterations without prec. setup (?)
fnormtol = 1.0d-8 ! Function-norm stopping tolerance
scsteptol = 1.0d-8 ! Function-norm stopping tolerance

maxl = 5000 ! maximum Krylov subspace dimesion (?!?!?!) ! Esto se usa para el preconditioner
maxlrst = 5 ! maximum number of restarts
max_niter = 2000
globalstrat = 0

call fnvinits(3, neq, ier) ! fnvinits inits NVECTOR module
if (ier .ne. 0) then       ! 3 for Kinsol, neq ecuantion number, ier error flag (0 is OK)
  print*, 'call_kinsol: SUNDIALS_ERROR: FNVINITS returned IER = ', ier
  call MPI_FINALIZE(ierr) ! finaliza MPI
  stop
endif

call fkinmalloc(iout, rout, ier)    ! Allocates memory and output additional information
if (ier .ne. 0) then
   print*, 'call_kinsol: SUNDIALS_ERROR: FKINMALLOC returned IER = ', ier
   call MPI_FINALIZE(ierr) ! finaliza MPI
   stop
 endif

call fkinsetiin('MAX_SETUPS', msbpre, ier)  ! Additional input information
call fkinsetrin('FNORM_TOL', fnormtol, ier)
call fkinsetrin('SSTEP_TOL', scsteptol, ier)
call fkinsetiin('MAX_NITER', max_niter, ier)
call fkinsetiin('MAX_NITERS', max_niter, ier)

do i = 1, ntot  !constraint vector
   constr(i) = 2.0
enddo
do i = ntot+1, 2*ntot  !constraint vector
constr(i) = 0.0        ! cualquier cosa
enddo
do i = 2*ntot+1, 3*ntot  
   constr(i) = 2.0
enddo

call fkinsetvin('CONSTR_VEC', constr, ier) ! constraint vector
! CALL FKINSPTFQMR (MAXL, IER)
call fkinspgmr(maxl, maxlrst, ier) !  Scale Preconditioned GMRES solution of linear system (???)
!call fkinspbcg(maxl, ier) !  Scale Preconditioned BCG

if (ier .ne. 0) then
  print*, 'call_kinsol: SUNDIALS_ERROR: FKINSPGMR returned IER = ', ier
  call fkinfree ! libera memoria
  call MPI_FINALIZE(ierr) ! finaliza MPI
  stop
endif
call fkinspilssetprec(1, ier) ! preconditiones

do i = 1, neq ! scaling vector
  scale(i) = 1.0
enddo

do i = 1, neq ! Initial guess
      x1(i) = x1_old(i)
      xg1(i) = x1(i)  
enddo


call fkinsol(x1, globalstrat, scale, scale, ier)         ! Llama a kinsol

print*, 'IER:', ier

if (ier .lt. 0) then
      print*, 'call_kinsol: SUNDIALS_ERROR: FKINSOL returned IER = ', ier
      print*, 'call_kinsol: Linear Solver returned IER = ', iout(9)
      stop
endif

do i = 1, neq ! output
  x1_old(i) = x1(i)
  xg1_old(i) = x1(i)
enddo

call fkinfree
return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subrutina que llama a kinsol


subroutine call_fkfun(x1_old)
use multicapa
use MPI

integer i

real*8 x1_old(3*ntot)
real*8 x1(3*ntot)
real*8 f(3*ntot)

! MPI

integer tag
parameter(tag = 0)
integer err

x1 = 0.0
do i = 1,3*ntot
  x1(i) = x1_old(i)
enddo

CALL MPI_BCAST(x1, 3*ntot , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)

call fkfun(x1,f, ier) ! todavia no hay solucion => fkfun 
end

