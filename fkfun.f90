subroutine fkfun(x,f,ier2)
use multicapa
use partfunc
use layer
use volume
use bulk
use longs
use MPI
implicit none

integer*4 ier2
real*8 protemp
real*8 x(ntot),f(ntot)
real*8 xh(2*ntot)
real*8 xpot(2*ntot)
real*8 pro(maxcuantas)
integer i,j,k,k1,k2,ii, jj,iz       ! dummy indices
integer err
INTEGER AT
REAL*8 cortepegado
real*8 beta, gammap, gamman
real*16 auxB, auxC
REAL*8 avpolpos(ntot), avpolneg(ntot)
REAL*8 ALGO, ALGO2
integer n
real*8 avpol_tmp(2*ntot), avpol2_tmp(2*ntot)
real*8 avpol_red(ntot)
real*8 avpol_tosend(ntot), avpol2_tosend(ntot)
real*8 nnn

! Jefe
if(rank.eq.0) then ! llama a subordinados y pasa vector x
   flagsolver = 1
   CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
   CALL MPI_BCAST(x, ntot , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)
endif



! common declarations

n = ntot
cortepegado = 0.0
pegado = 2


! avpolall
avpolall = 0
do j = 1, nads ! loop over adsorbed layers
do i = 1, ntot
avpolall(i) = avpolall(i) + avpol(j,i)
!if((iter.eq.1368).and.(i.eq.2))print*, j, avpol(j,i)
end do
end do

! chain parameters

long(1) = long1
long(2) = long2

cuantas(1) = cuantas1
cuantas(2) = cuantas2

do i=1,n                  ! init xh and psi
if(x(i).lt.0.0)x(i)=1d-30
xh(i)=(exp(-x(i)))*(1.0-avpolall(i))             ! solvent density=volume fraction   
xh(i+n)=xsolbulk       ! bulk value for chains protruding to the bulk
enddo

fbound = 0.0

! calculo de xtotal para poor solvent

do i = 1,n
xtotal(i) = 1.0 - xh(i)
xtotal(i+n) = 0.0
enddo

! stoichoimetry

do i = 1, n
avpolpos(i) = 0.0
avpolneg(i) = 0.0
end do

do j = 1, nads ! loop over adsorbed layers
 if (Tcapas(j).eq.1) THEN
  do i = 1, n
   avpolpos(i) = avpolpos(i) + avpol(j, i)
  end do
 else
  do i = 1, n
   avpolneg(i) = avpolneg(i) + avpol(j, i)
  end do
 endif
enddo

! not adsobed

if (Tcapas(nads+1).eq.1) THEN ! adsorbs positive
do i = 1, n
  avpolpos(i) = xtotal(i) - avpolneg(i)
enddo

else  ! adsorbs negative

do i = 1, n
  avpolneg(i) = xtotal(i) - avpolpos(i)
end do
endif

avpolneg(1) = avpolneg(1) + sigma

! maxpol : position of the last layer with complementary polymer

maxpol = 1

if (nads.gt.0) then
 do i = 1,n
  IF(avpol(nads, i).gt.0.0)maxpol=i
 end do
ENDif


!!!
!!! calculation of fbound 
!!!

do i=1,maxpol ! see notes, A = pos = 1, B = neg = 2
auxC = avpolneg(i)/avpolpos(i)
auxB = -1.0 -auxC - 1.0/Kbind0/avpolpos(i)

fbound(1, i) = (-auxB - SQRT(auxB**2 - 4.0*auxC))/2.0
fbound(2, i) = avpolpos(i)*fbound(1, i)/avpolneg(i)
enddo

!if (nads.eq.1) then
!do i = 1, ntot
!print*, i , xh(i), xtotal(i)
!enddo
!stop
!endif

do i=1,n
avpol(nads+1,i)=0.0d0         ! polymer volume fraction
enddo

do i=1,n                 ! polymer to adsorb
avpol2(i)=0.0d0         ! polymer volume fraction
enddo

!! pegado has the last layer with polymer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! LT = type of polymer being adsorbed
! nads = number of layers already adsorbed
!
! avpol => volume fraction of all polimer
! avpol2 => volume fraction for z < pegado



! xpot
AT = Tcapas(nads+1) ! type of layer to add

do i = 1, n
protemp = dlog(xh(i)**(vpol))
if(nads.eq.0)protemp = protemp + (-eps(i))
protemp = protemp-dlog(1.0-fbound(AT, i))

do iz = -Xulimit, Xulimit
if(((iz+i).ge.1).and.(iz+i).le.ntot) then
if(AT.eq.1) then ! pos
protemp=protemp + Xu(1,2,iz)*st/(vpol*vsol)*avpolneg(i+iz)
protemp=protemp + Xu(1,1,iz)*st/(vpol*vsol)*avpolpos(i+iz)
else ! neg
protemp=protemp + Xu(1,1,iz)*st/(vpol*vsol)*avpolneg(i+iz)
protemp=protemp + Xu(1,2,iz)*st/(vpol*vsol)*avpolpos(i+iz)
endif
endif
enddo

xpot(i) = dexp(protemp)
enddo

do i = n+1, 2*n
protemp = dlog(xsolbulk**(vpol))
xpot(i) = dexp(protemp)
enddo

!    probability distribution


q=0.0d0                   ! init q to zero
avpol_tosend = 0.0
avpol_tmp = 0.0
avpol2_tosend = 0.0
avpol2_tmp = 0.0

do ii=1,ntot 
 do i=1,newcuantas(AT)

 pro(i) = expmupol*weight(AT,i)
 nnn = 0.0

 if(ii.le.maxpol) then
    do j=1, maxlayer(AT, i)
     k = j+ii-1
     pro(i)= pro(i) * xpot(k)**in1n(AT, i, j)
     nnn = nnn + in1n(AT,i,j)*fbound(AT,j)
    enddo
 endif

    q=q+pro(i)

    do j=1,maxlayer(AT, i)
     k = j+ii-1
     avpol2_tmp(k)=avpol2_tmp(k)+pro(i)*vpol*in1n(AT,i, j)  ! volume fraction polymer not normed
    enddo

      if(nnn.ge.minn) then
      do j=1,maxlayer(AT, i)
       k = j+ii-1
       avpol_tmp(k)=avpol_tmp(k)+pro(i)*vpol*in1n(AT,i, j) ! only bound polymer!!!
      enddo
    endif
 enddo ! i
enddo   ! ii

avpol_tosend(1:ntot)=avpol_tmp(1:ntot)
avpol2_tosend(1:ntot)=avpol2_tmp(1:ntot)
!------------------ MPI -----------------`-----------------------------
!1. Todos al jefe


avpol_red = 0.0
call MPI_Barrier(MPI_COMM_WORLD, err)

! Jefe
if (rank.eq.0) then
! Junta avpol       
  call MPI_REDUCE(avpol_tosend, avpol_red, ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
  call MPI_REDUCE(avpol2_tosend, avpol2, ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
endif
! Subordinados
if(rank.ne.0) then
! Junta avpol       
  call MPI_REDUCE(avpol_tosend, avpol_red, ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
  call MPI_REDUCE(avpol2_tosend, avpol2, ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
!!!!!!!!!!! IMPORTANTE, LOS SUBORDINADOS TERMINAN ACA... SINO VER !MPI_allreduce!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  goto 3333
endif

avpol(nads+1,:) = avpol_red(:)

!print*, nads
!do i = 1, ntot
!print*, i, avpol(nads+1,i), avpol2(i)
!enddo
!stop

! contruction of f and the volume fractions

do i=1,n
 f(i)=xh(i)-1.0d0
do jj = 1, (nads)
 f(i) = f(i) + avpol(jj, i)
end do
 f(i) = f(i) + avpol2(i)
enddo

iter=iter+1

algo = 0.0
do i = 1, n
 algo = algo + f(i)**2
end do

if(rank.eq.0)PRINT*, iter, algo
norma=algo

3333 continue

return
end
