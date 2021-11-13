subroutine fkfun(x,f,ier2)
use multicapa
use partfunc
use layer
use volume
use bulk
use longs
use MPI
use const
implicit none

integer ier2
real*8 protemp
real*8 x(ntot*3),f(ntot*3)
real*8 psi2(0:ntot+1) ! psipsi plus boundaries at z=0 and dimz+1
real*8 xh(ntot)
real*8 xpot(2*ntot)
!real*8 pro(maxcuantas)  !!G
integer i,j,k1,k2,ii, jj,iz       ! dummy indices
integer err
!INTEGER AT
real*8 beta, gammap, gamman
real*8 auxB, auxC
!REAL*8 avpolpos(ntot), avpolneg(ntot)  !G
REAL*8 ALGO, ALGO2
integer n
real*8 avpol_tmp(ntot), avpol2_tmp(ntot)
real*8 avpol_red(ntot)
real*8 avpol_tosend(ntot), avpol2_tosend(ntot),rhopol2_tosend(ntot)
real*8 nnn,rhopol2_tmpvalue
real*8 rhopol2_tmp(ntot) !Nuevo
integer k
real*8 q_tosend(ntot), splp_tosend(ntot)
real*8 q0(ntot), splp0(ntot)
real*8  Check_kbind,check_Kaplus,check_kbmin,check_KANA,check_KBCL
double precision, external :: factorcurv
real*8 alfaA, alfaB, betaA, betaB, gamaA, gamaB

integer lim1,lim2,lim3,lim4 ! limits the regions for fbound calculation for different curvatures


! Jefe
if(rank.eq.0) then ! llama a subordinados y pasa vector x
   flagsolver = 1
   CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
   CALL MPI_BCAST(x, 3*ntot , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)
endif

n = ntot

! common declarations
! avpolall
avpolall = 0
do j = 0, nads ! loop over adsorbed layers
do i = 1, ntot
avpolall(i) = avpolall(i) + avpol(j,i)
end do
end do

! chain parameters

long(1) = long1
long(2) = long2

cuantas(1) = cuantas1
cuantas(2) = cuantas2

do i=1,n                 
xh(i)=x(i)
enddo
fbound = 0.0

do i=1,n
psi(i)=x(i+n)
enddo

do i=1,n
xtotal(i)=x(i+2*n) ! xtotal = volume fraction of adsorbed layer from input
enddo

! Electrostatic potential boundary conditions
psi2(1:n) = psi(1:n)

if (curvature.ge.0) then
  psi2(n+1) = 0.0 ! wall, no charge
  psi2(0) = psi2(1)  ! wall, no charge
else
  psi2(n+1) = psi(n) ! wall, no charge
  psi2(0) = psi2(1)  ! symmetry at channel center dpsi/dr = 0
endif

do i=1,n
   xNaCl(i) = expmuNaCl*(xh(i)**(2.0*vsalt)) ! NaCl
   xpos(i) = expmupos*(xh(i)**vsalt)*exp(-psi2(i)*zpos) ! ion plus volume fraction
   xneg(i) = expmuneg*(xh(i)**vsalt)*exp(-psi2(i)*zneg) ! ion neg volume fraction
   xHplus(i) = expmuHplus*(xh(i))*exp(-psi2(i))         ! H+ volume fraction
   xOHmin(i) = expmuOHmin*(xh(i))*exp(+psi2(i))         ! OH-  volume fraction
enddo

!print*,'xpos,xneg,xhplus,xohmin',xpos(10),xneg(10),xhplus(10),xohmin(10)
!stop

! stoichoimetry

avpolpos = 0.0
avpolneg = 0.0
avpolposcero=0.0
avpolnegcero=0.0
do j = 0, nads ! loop over adsorbed layers

 if (Tcapas(j).eq.1) THEN
  do i = 1, n
   avpolneg(i) = avpolneg(i) + avpol(j, i)
   avpolnegcero(i)= avpolnegcero(i) +avpol(j,i)
  end do
 else
  do i = 1, n
   avpolpos(i) = avpolpos(i) + avpol(j, i)
   avpolposcero(i)=avpolposcero(i) + avpol(j,i)
  end do
 endif
enddo

! not adsobed
if (Tcapas(nads+1).eq.1) THEN ! adsorbs negative
do i = 1, n
  avpolneg(i) = avpolneg(i) + xtotal(i) 
enddo
else  ! adsorbs positive
do i = 1, n
  avpolpos(i) = avpolpos(i) + xtotal(i)
end do
endif

! maxpol : position of the last layer with complementary polymer
maxpol = radio
if (nads.gt.0) then

  if (curvature.lt.0) then
  do i = n,1,-1
   IF(avpol(nads, i).gt.0.0)maxpol=i
  end do

  else

  do i = 1, n
   IF(avpol(nads, i).gt.0.0)maxpol=i
  end do
  endif

endif


!!! calculation of fbound 

! find limits for different curvatures

if (curvature.lt.0) then
  lim1 = maxpol
  lim2 = ntot
  lim3 = 1
  lim4 = maxpol-1
else
  lim1 = radio
  lim2 = maxpol  
  lim3 = maxpol+1
  lim4 = ntot
endif


do i=lim1,lim2

! Escrito de nuevo siguiendo notacion del paper JCP (MT 15-6-2021)
! neg (1) = polimero negativo (A-)
! pos (2) = polimero positivo (B+)

  auxC = avpolpos(i)/avpolneg(i)

  betaA = K0A*xh(i)/xHplus(i)
  betaB = K0B*xh(i)/xOHmin(i)
  alfaA = K0ANa*(xpos(i)/vsalt)/(xh(i)**vsalt)
  alfaB = K0BCl*(xneg(i)/vsalt)/(xh(i)**vsalt)
  gamaA = 1.0+alfaA+1.0/betaA
  gamaB = 1.0+alfaB+1.0/betaB

  auxB = -1.0 -auxC -gamaA*gamaB/Kbind0/(avpolneg(i)/vpol/vsol)
  
  fbound(1, i) = (-auxB - SQRT(auxB**2 - 4.0*auxC))/2.0
  fbound(2, i) = fbound(1,i)/auxC

  fNcharge(1,i) = (1.0 - fbound(1,i))/betaA/gamaA
  fNcharge(2,i) = (1.0 - fbound(2,i))/betaB/gamaB

  fioncharge(1,i) = (1.0 - fbound(1,i))*alfaA/gamaA
  fioncharge(2,i) = (1.0 - fbound(2,i))*alfaB/gamaB



  !Check_Kbind= fbound(1,i)/((1.0-fbound(1,i)-fNcharge(1,i)-fioncharge(1,i))  &
  !  *(1.0-fbound(2,i)-fNcharge(2,i)-fioncharge(2,i))*avpolpos(i)/vpol/vsol) -Kbind0


  !Check_Kbmin=  (xOhmin(i)/xh(i))*(1.0-fbound(2,i)-fNcharge(2,i)-fioncharge(2,i)) & 
  !            /fNcharge(2,i)-K0B

 ! Check_Kaplus= (xHplus(i)/xh(i))*(1.0-fbound(1,i)-fNcharge(1,i)-fioncharge(1,i)) &
 !             /fNcharge(1,i)-K0A

 !! Check_KANa = (xh(i)**vsalt)*fioncharge(1,i)/(1.0-fbound(1,i)-fNcharge(1,i)-fioncharge(1,i)) &
 !             /(xpos(i)/vsalt) - K0ANa

 ! Check_KBCl = (xh(i)**vsalt)*fioncharge(2,i)/(1.0-fbound(2,i)-fNcharge(2,i)-fioncharge(2,i)) &
  !            /(xneg(i)/vsalt) - K0BCl

!! Chequeos
!print*,'checkeo Kbind, KBmin, KAplus, KBCl, KANa',Check_Kbind, Check_Kaplus,Check_Kbmin,Check_KBCl, Check_KANa
!print*,'f asociado', fbound(1, i),fbound(2, i)
!print*,'f nc' ,fNcharge(1,i),fNcharge(2,i)
!print*,'f ioncharge',fioncharge(1,i),fioncharge(2,i)

enddo
!stop
! Solo polimero bulk, sin asociacion
do i=lim3, lim4

  betaA = K0A*xh(i)/xHplus(i)
  betaB = K0B*xh(i)/xOHmin(i)
  alfaA = K0ANa*(xpos(i)/vsalt)/(xh(i)**vsalt)
  alfaB = K0BCl*(xneg(i)/vsalt)/(xh(i)**vsalt)

  fNcharge(1,i) = 1.0/(1.0+betaA+betaA*alfaA)
  fNcharge(2,i) = 1.0/(1.0+betaB+betaB*alfaB)

  fioncharge(1,i) = alfaA*betaA*fNcharge(1,i)
  fioncharge(2,i) = alfaB*betaB*fNcharge(2,i)


!Check_Kbmin=  (xOhmin(i)/xh(i))*(1.0-fbound(2,i)-fNcharge(2,i)-fioncharge(2,i)) & 
!              /fNcharge(2,i)-K0B

!Check_Kaplus= (xHplus(i)/xh(i))*(1.0-fbound(1,i)-fNcharge(1,i)-fioncharge(1,i)) &
!              /fNcharge(1,i)-K0A

!Check_KANa = (xh(i)**vsalt)*fioncharge(1,i)/(1.0-fbound(1,i)-fNcharge(1,i)-fioncharge(1,i)) &
!              /(xpos(i)/vsalt) - K0ANa

!Check_KBCl = (xh(i)**vsalt)*fioncharge(2,i)/(1.0-fbound(2,i)-fNcharge(2,i)-fioncharge(2,i)) &
!              /(xneg(i)/vsalt) - K0BCl



!print*,'checkeo KBmin, KAplus, KBCl, KANa',Check_Kaplus,Check_Kbmin,Check_KBCl, Check_KANa
!print*,'f asociado', fbound(1, i),fbound(2, i)
!print*,'f nc' ,fNcharge(1,i),fNcharge(2,i)
!print*,'f ioncharge',fioncharge(1,i),fioncharge(2,i)
!stop

enddo

avpol(nads+1,:)=0.0d0         ! polymer volume fraction
avpol2=0.0d0         ! polymer volume fraction

!! pegado has the last layer with polymer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LT = type of polymer being adsorbed
! nads = number of layers already adsorbed
! avpol => volume fraction of all polimer
! avpol2 => volume fraction for z < pegado
! xpot

AT = Tcapas(nads+1) ! type of layer to add

do i = 1, ntot

protemp = dlog(xh(i)**(vpol))
protemp = protemp-dlog(1.0-fbound(AT, i)-fNcharge(AT,i)-fioncharge(AT,i))
protemp = protemp-psi2(i)*zpol(AT)


!do iz = -Xulimit, Xulimit
!if((iz+i).ge.1) then
!if(AT.eq.1) then ! pos
!protemp=protemp + Xu(1,2,iz)*st/(vpol*vsol)*avpolneg(i+iz) ! use poor solvent only for curvature = 0
!protemp=protemp + Xu(1,1,iz)*st/(vpol*vsol)*avpolpos(i+iz)
!else ! neg
!protemp=protemp + Xu(1,1,iz)*st/(vpol*vsol)*avpolneg(i+iz)
!protemp=protemp + Xu(1,2,iz)*st/(vpol*vsol)*avpolpos(i+iz)
!endif
!endif
!enddo

xpot(i) = dexp(protemp)


enddo

do i = n+1, 2*n
xpot(i) = xpot(n)
enddo

!    probability distribution
avpol_tosend = 0.0
avpol_tmp = 0.0
avpol2_tosend = 0.0
avpol2_tmp = 0.0
sumprolnpro = 0.0
rhopol2_tosend=0.0
rhopol2_tmp=0.0


splp0 = 0.0
q0 = 0.0

do ii=1,ntot ! position of segment #0 

splp_tosend(ii) = 0.0
q_tosend(ii)=0.0d0                   ! init q to zero (para c/layer)

 do i=1,newcuantas(AT)

 if(weight(AT,i,ii).ne.0.0) then

 pro(i) = 1.0
 nnn = 0.0


! print*, ii, i, minpos(AT,i,ii), maxpos(AT, i, ii)


    do j=minpos(AT,i,ii), maxpos(AT,i,ii) ! posicion dentro del poro

! reflecting boundary condition
     jj = j
     if (jj.gt.n)jj=2*n-jj+1

     k = j-minpos(AT,i,ii)+1
     pro(i)= pro(i) * xpot(jj)**in1n(AT,i,ii,k)

     nnn = nnn + in1n(AT,i,ii,k)*fbound(AT,jj)

    enddo

    q_tosend(ii)=q_tosend(ii)+pro(i) ! single-chain partition function
    splp_tosend(ii) = splp_tosend(ii) + pro(i)*log(pro(i)) 
    rhopol2_tmp(ii)=rhopol2_tmp(ii) + pro(i)*expmupol/vsol*weight(AT,i,ii)
    
   do j=minpos(AT,i,ii), maxpos(AT,i,ii)
     k = j-minpos(AT,i,ii)+1

! reflecting boundary condition
     jj = j
     if (jj.gt.n)jj=2*n-jj+1

     avpol2_tmp(jj)=avpol2_tmp(jj)+pro(i)*expmupol*weight(AT,i,ii)*vpol*in1n(AT,i,ii,k)*factorcurv(ii,jj)  ! volume fraction polymer not normed
    enddo



      if(nnn.ge.minn) then
      do j=minpos(AT,i,ii), maxpos(AT,i,ii)

! reflecting boundary condition
     jj = j
     if (jj.gt.n)jj=2*n-jj+1

       k = j-minpos(AT,i,ii)+1
       avpol_tmp(jj)=avpol_tmp(jj)+pro(i)*expmupol*weight(AT,i,ii)*vpol*in1n(AT,i,ii,k)*factorcurv(ii,jj) ! only bound polymer!!!
      enddo
    endif

 endif ! weight
 
 enddo ! i

enddo   ! ii

avpol_tosend=avpol_tmp
avpol2_tosend=avpol2_tmp

rhopol2_tosend=rhopol2_tmp

!------------------ MPI -----------------`-----------------------------
!1. Todos al jefe


avpol_red = 0.0
call MPI_Barrier(MPI_COMM_WORLD, err)

! Jefe
  call MPI_REDUCE(avpol_tosend, avpol_red, ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
  call MPI_REDUCE(avpol2_tosend, avpol2, ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
  call MPI_REDUCE(rhopol2_tosend, rhopol2, ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
  call MPI_REDUCE(splp_tosend, splp0, ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
  call MPI_REDUCE(q_tosend, q0, ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)



do ii = 1, ntot
 sumprolnpro(ii) = splp0(ii)/q0(ii) - log(q0(ii))
enddo


if(rank.ne.0) then
  goto 3333
endif


avpol(nads+1,:) = avpol_red(:)

! contruction of f and the volume fractions


do i=1,n
 f(i)=xh(i)+ xNaCl(i) + xneg(i) + xpos(i) + xHplus(i) + xOHmin(i) -1.0d0
do jj = 0, nads
 f(i) = f(i) + avpol(jj, i)
end do
 f(i) = f(i) + avpol2(i)

enddo

!print*,'test',xh(10),xneg(10),xpos(10),xhplus(10),xohmin(10)
!stop


do i=1,n
 qtot(i) = (zpos*xpos(i)+zneg*xneg(i))/vsalt + avpolneg(i)*zpol(1)/vpol*(1.0-fbound(1,i)-fNcharge(1,i)-fioncharge(1,i))& !!
+ avpolpos(i)*zpol(2)/vpol*(1.0-fbound(2,i)-fNcharge(2,i)-fioncharge(2,i)) + xHplus(i)-xOHmin(i)                       !
enddo
!do i=1,n
! qtot(i) = (zpos*xpos(i)+zneg*xneg(i))/vsalt + avpolneg(i)*zpol(1)/vpol*(1.0-fbound(1,i)-fNcharge(1,i))& !!
!+ avpolpos(i)*zpol(2)/vpol*(1.0-fbound(2,i)-fNcharge(2,i)) + xHplus(i)-xOHmin(i)                       !
!enddo
!print*,'1est', (zpos*xpos(5)+zneg*xneg(5))/vsalt
!print*,'charg', avpolneg(5)*zpol(1)/vpol*(1.0-fbound(1,5)-fNcharge(1,5)-fioncharge(1,5))
!print*,'qtot',qtot(5),avpolneg(5),fbound(1,5),fncharge(1,5),fioncharge(1,5),avpolpos(5)
!print*,'test',fncharge(2,5),fbound(2,5)
!print*,'qtto',qtot(5)
!stop

!i = 1
!print*, '!', (1.0-fbound(1,1)-fNcharge(1,1)),qtot(1)
!print*, zpos*xpos(i)/vsalt
!print*, zneg*xneg(i)/vsalt 
!print*, avpolneg(i)*zpol(1)/vpol*(1.0-fbound(1,i)-fNcharge(1,i))
!print*, avpolpos(i)*zpol(2)/vpol*(1.0-fbound(2,i)-fNcharge(2,i))
!print*, xHplus(i)
!print*, xOHmin(i)

!i = 1
!print*, 'salt', (zpos*xpos(1)+zneg*xneg(1))/vsalt
!print*, 'polneg', avpolneg(i)*zpol(1)/vpol*(1.0-fbound(1,i)-fNcharge(1,i))
!print*, 'polpos', avpolpos(i)*zpol(2)/vpol*(1.0-fbound(2,i)-fNcharge(2,i))
!print*, 'others', xHplus(i)-xOHmin(i)
!print*,'avpolpos', avpolpos(i)
!print*,'zpol(2)', zpol(2)
!print*,'fractionb', fbound(2,i)
!print*,'fractionN',fNcharge(2,i)
!stop


wperm = 0.114 !water permitivity in units of e^2/kT.nm

do i = 1,n 
   select case (abs(curvature))

    case (0)
     f(i + n)=qtot(i)/vsol &
     +wperm*(psi2(i+1)-2.0*psi2(i)+psi2(i-1))*delta**(-2)

    case(1)
     f(i + n)=qtot(i)/vsol &
     + wperm*(psi2(i+1)-psi2(i))*delta**(-2)/(float(i)-0.5) &
     +wperm*(psi2(i+1)-2.0*psi2(i)+psi2(i-1))*delta**(-2) 

    case(2)
     f(i + n)=qtot(i)/vsol &
     + 2.00*wperm*(psi2(i+1)-psi2(i))*delta**(-2)/(float(i)-0.5) &
     + wperm*(psi2(i+1)-2.0*psi2(i)+psi2(i-1))*delta**(-2) 

    end select

f(n+i)=f(n+i)/(-2.0)

enddo

do i = 1,n
f(i+2*n) = -avpol2(i)+xtotal(i) 
enddo

!print*,'psi2 avpol2', psi2(10),avpol2(10)
!stop
iter=iter+1

algo = 0.0
do i = 1, n*3
 algo = algo + f(i)**2

end do

if(rank.eq.0)PRINT*, iter, algo
norma=algo


3333 continue




ier2 = 0
return
end
