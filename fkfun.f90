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
real*8 x(2*ntot),f(2*ntot)
real*8 xh(2*ntot)
real*8 xpot(2*ntot,2)
real*8 pro(maxcuantas)
integer i,j,k,k1,k2,ii, jj,iz       ! dummy indices
integer err
INTEGER AT
REAL*8 cortepegado
real*8 beta, gammap, gamman
real*16 auxB, auxC
REAL*8 avpolpos(2*ntot,2), avpolneg(2*ntot,2)
REAL*8 ALGO, ALGO2
integer n
real*8 avpol_tmp(2*ntot,2), avpol2_tmp(2*ntot,2)
real*8 avpol_red(ntot,2)
real*8 avpol_tosend(ntot,2), avpol2_tosend(ntot,2)
real*8 nnn

! Jefe
if(rank.eq.0) then ! llama a subordinados y pasa vector x
   flagsolver = 1
   CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
   CALL MPI_BCAST(x, 2*ntot , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)
endif

! common declarations

n = ntot
cortepegado = 0.0
pegado = 2

! avpolall
avpolall = 0
do j = 1, nads ! loop over adsorbed layers
do i = 1, ntot
avpolall(i,:) = avpolall(i,:) + avpol(j,i,:)
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

xh(i)=(exp(-x(i)))*(1.0-avpolall(i,1)-avpolall(i,2))             ! solvent density=volume fraction   
xh(i+n)=xsolbulk       ! bulk value for chains protruding to the bulk

xtotal(i,1) = (1.0-xh(i))*(exp(-x(i+n))) ! fraction of segments type 1 (sticky!)
xtotal(i,2) = (1.0-xh(i))*(1.0-exp(-x(i+n))) ! fraction of segments type 2 (non-sticky)
enddo

fbound = 0.0

! calculo de xtotal in bulk
AT = Tcapas(nads+1) ! type of layer to add
do i = 1,n
xtotal(i+n,1) = phibulkpol*float(sticky(AT))/float(long(AT))
xtotal(i+n,2) = phibulkpol*(1.0-float(sticky(AT))/float(long(AT)))
enddo

! stoichoimetry
avpolpos = 0.0
avpolneg = 0.0

do j = 1, nads ! loop over adsorbed layers
 if (Tcapas(j).eq.1) THEN
  do i = 1, n
   avpolpos(i,:) = avpolpos(i,:) + avpol(j,i,:)
  end do
 else
  do i = 1, n
   avpolneg(i,:) = avpolneg(i,:) + avpol(j,i,:)
  end do
 endif
enddo

! not adsobed

if (Tcapas(nads+1).eq.1) THEN ! adsorbs positive
do i = 1, 2*n
  avpolpos(i,:) = xtotal(i,:)-avpolneg(i,:)
enddo

else  ! adsorbs negative

do i = 1, 2*n
  avpolneg(i,:) = xtotal(i,:)-avpolpos(i,:)
end do
endif

avpolneg(1,1) = avpolneg(1,1) + sigma

! maxpol : position of the last layer with complementary polymer
! remember avpol(:,:,n) n = 1 for sticky and 2 for non-sticky

maxpol = 1
if (nads.gt.0) then
 do i = 1,n
  IF(avpol(nads,i,1).gt.0.0)maxpol=i
 end do
ENDif


!!!
!!! calculation of fbound 
!!!

do i=1,maxpol ! see notes, A = pos = 1, B = neg = 2
auxC = avpolneg(i,1)/avpolpos(i,1)
auxB = -1.0 -auxC - 1.0/Kbind0/avpolpos(i,1)

fbound(1, i) = (-auxB - SQRT(auxB**2 - 4.0*auxC))/2.0
fbound(2, i) = avpolpos(i,1)*fbound(1, i)/avpolneg(i,1)
enddo

!if (nads.eq.1) then
!do i = 1, ntot
!print*, i , xh(i), xtotal(i)
!enddo
!stop
!endif

avpol(nads+1,:,:)=0.0d0         ! polymer volume fraction
avpol2=0.0d0         ! polymer volume fraction

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

do iz = -Xulimit, Xulimit
if((iz+i).ge.1) then
if(AT.eq.1) then ! pos
protemp=protemp + Xu(1,2,iz)*st/(vpol*vsol)*avpolneg(i+iz,1)
protemp=protemp + Xu(1,2,iz)*st/(vpol*vsol)*avpolneg(i+iz,2)
protemp=protemp + Xu(1,1,iz)*st/(vpol*vsol)*avpolpos(i+iz,1)
protemp=protemp + Xu(1,1,iz)*st/(vpol*vsol)*avpolpos(i+iz,2)
else ! neg
protemp=protemp + Xu(1,1,iz)*st/(vpol*vsol)*avpolneg(i+iz,1)
protemp=protemp + Xu(1,1,iz)*st/(vpol*vsol)*avpolneg(i+iz,2)
protemp=protemp + Xu(1,2,iz)*st/(vpol*vsol)*avpolpos(i+iz,1)
protemp=protemp + Xu(1,2,iz)*st/(vpol*vsol)*avpolpos(i+iz,2)
endif
endif
enddo

xpot(i,2) = dexp(protemp)
protemp = protemp-dlog(1.0-fbound(AT, i))
xpot(i,1) = dexp(protemp)
enddo

do i = n+1, 2*n
xpot(i,:) = xpot(n,:)
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

    do j=1, maxlayer(AT, i)
     k = j+ii-1
     pro(i)= pro(i) * xpot(k,1)**in1n(AT, i, j,1)
     pro(i)= pro(i) * xpot(k,2)**in1n(AT, i, j,2)
     nnn = nnn + in1n(AT,i,j,1)*fbound(AT,k)
    enddo

    q=q+pro(i)

    do j=1,maxlayer(AT, i)
     k = j+ii-1
     avpol2_tmp(k,1)=avpol2_tmp(k,1)+pro(i)*vpol*in1n(AT,i, j,1)  ! volume fraction polymer sticky
     avpol2_tmp(k,2)=avpol2_tmp(k,2)+pro(i)*vpol*in1n(AT,i, j,2)  ! volume fraction polymer non-sticky
    enddo

      if(nnn.ge.minn) then
      do j=1,maxlayer(AT,i)
       k = j+ii-1
       avpol_tmp(k,1)=avpol_tmp(k,1)+pro(i)*vpol*in1n(AT,i,j,1) ! only bound polymer!!!
       avpol_tmp(k,2)=avpol_tmp(k,2)+pro(i)*vpol*in1n(AT,i,j,2) ! only bound polymer!!!
      enddo
    endif
 enddo ! i
enddo   ! ii

avpol_tosend(1:ntot,1)=avpol_tmp(1:ntot,1)
avpol_tosend(1:ntot,2)=avpol_tmp(1:ntot,2)
avpol2_tosend(1:ntot,1)=avpol2_tmp(1:ntot,1)
avpol2_tosend(1:ntot,2)=avpol2_tmp(1:ntot,2)

!------------------ MPI -----------------`-----------------------------
!1. Todos al jefe

avpol_red = 0.0
call MPI_Barrier(MPI_COMM_WORLD, err)

! Jefe
if (rank.eq.0) then
! Junta avpol       
  call MPI_REDUCE(avpol_tosend, avpol_red, 2*ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
  call MPI_REDUCE(avpol2_tosend, avpol2, 2*ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
endif
! Subordinados
if(rank.ne.0) then
! Junta avpol       
  call MPI_REDUCE(avpol_tosend, avpol_red, 2*ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
  call MPI_REDUCE(avpol2_tosend, avpol2, 2*ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
!!!!!!!!!!! IMPORTANTE, LOS SUBORDINADOS TERMINAN ACA... SINO VER !MPI_allreduce!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  goto 3333
endif

avpol(nads+1,:,:) = avpol_red(:,:)


! contruction of f and the volume fractions

! xh = 0.0
do i=1,n

 f(i)=xh(i)-1.0d0
 f(i+n) = 0.0

do jj = 1, (nads)
 f(i) = f(i) + avpol(jj,i,1) + avpol(jj,i,2) 
 f(i+n) = f(i+n) + avpol(jj,i,1)
end do

 f(i) = f(i) + avpol2(i,1)+avpol2(i,2)
 f(i+n) = f(i+n) + avpol2(i,1)
 f(i+n) = f(i+n)/f(i) - xtotal(i,1)/(xtotal(i,1)+xtotal(i,2)) ! fraction of sticky

enddo

! xtotal = avpol2


iter=iter+1

algo = 0.0
do i = 1, 2*n
 algo = algo + f(i)**2
end do

if(rank.eq.0)PRINT*, iter, algo
norma=algo

3333 continue

return
end
