
call initmpi
call read
call allocation
call kais
call solve
end

subroutine solve
use multicapa
use partfunc
use layer
use volume
use bulk
use seed1
use longs
use MPI
use const
implicit none
integer *4 ier ! Kinsol error flag
!real*8 Na               
!parameter (Na=6.02d23)

real*8 alfa, beta

real*8 avpol_red(ntot)

REAL*8 avtotal(ntot)       ! sum over all avpol
!real*8 xsol(ntot)         ! volume fraction solvent

real*8 x1(3*ntot),xg1(3*ntot)   ! density solvent iteration vector
real*8 zc(ntot)           ! z-coordinate layer 

REAL*8 sumrhoz, meanz     ! Espesor medio pesado
!real*8 pro                ! probability distribution function 

real*8 check_KANa,check_KBCl,check_KAplus,check_Kbminbulk

integer n                 ! number of lattice sites
integer itmax             ! maximum number of iteration allowed for 
real*8 xsalt
real*8 fnorm              ! L2 norm of residual vector function fcn
real*8 cHplus, cOHmin
real*8 pOHbulk
real*8 KaA
real*8 KaB

external fcnelect         ! function containing the SCMFT eqs for solver
integer i,j,k,m,ii,flag,c, jj ! dummy indice0s

INTEGER temp
real*8 tempr
real*8 tmp

real*8 min1               ! variable to determine minimal position of chain


integer il,inda,ncha

REAL*8 xfile(3*ntot)                        
real*8 algo, algo2                  

real*8 KaANa!!nw
real*8 KaBCl!!nW

integer*1 in1(maxlong)
real*8 chains(3,maxlong,ncha_max) ! chains(x,i,l)= coordinate x of segement i ,x=2 y=3,z=1
real*8 chainsw(ncha_max), sumweight_tosend
real*8 zp(maxlong)

real*8 sum,sumel          ! auxiliary variable used in free energy computation  
real*8 sumpi,sumrho,sumrhopol, sumrho2, sumrho2mol !suma de la fraccion de polimero

INTEGER LT ! layer type being adsorbed
integer ini,fin,step


! global running files
character*15 meanzfilename
character*15 sigmafilename
character*17 sigmaadfilename

! single layer files
character*18 sysfilename      ! contains value of free energy, input parameter etc
character*26 denssolfilename  ! contains the denisty of the solvent
character*29 denspolfilename(adsmax)
character*28 densendfilename
character*26 densbindfilename(2)
character*26 denschargefilename(2)
CHARACTER*24 totalfilename
CHARACTER*24 qtotfilename
character*27 denspol2filename


character*27 xposfilename
character*27 xnegfilename
character*27 xHplusfilename
character*27 xOHminfilename
character*27 potenfilename
character*27 rho2filename

integer countfile         ! enumerates the outputfiles 
integer countfileuno     ! enumerates the outputfiles para una corrida
integer conf              ! counts number of conformations

integer readsalt          !integer to read salt concentrations


INTEGER cc, ccc

! MPI
integer tag, source
parameter(tag = 0)
integer err
integer ier_tosend
double  precision norma_tosend

! iteration xsolbulk

real*8 xsoliter, xsolerror, xsolnorm

!
seed=435+ 3232*rank               ! seed for random number generator
print*, 'I am', rank, ' and my seed is', seed

if(rank.eq.0)print*, 'Program Multicapa'
if(rank.eq.0)print*, 'GIT Version: ', _VERSION

!     common declarations: used for communciations with other routines
long(1) = long1
long(2) = long2

cuantas(1) = cuantas1/size
cuantas(2) = cuantas2/size

!     initializations of variables 
!pi=dacos(-1.0d0)          ! pi = arccos(-1) 
itmax=200                 ! maximum number of iterations       
n=ntot                    ! size of lattice
conf=0                    ! counter for conformations

vsol=0.030                ! volume solvent molecule in (nm)^3
vpol= ((4.0/3.0)*pi*(0.3)**3)/vsol  ! volume polymer segment in units of vsol r=0.3
!print*,vpol,0.1/vsol

!!!!!GGG!!
zpos = 1.0      ! charge of cation
zneg = -1.0     ! charge of anion
zpol(1) = -1.0      ! charge of polyelectrolyte segment A
zpol(2) = 1.0      ! charge of polyelectrolyte segment B
vsalt=((4.0/3.0)*pi*(0.27)**3)/vsol  ! volume salt in units of vsol 0.27=radius salt  

pKw = 14.0 ! -log10(Kw)
kW=10**(-pKw)
KaA=10**(-pKaA)
KaB=10**(-pKaB)
KaANa=10**(-pKaANa)!!NW
KaBCl=10**(-pKaBCl)!!NW


cHplus = 10**(-pHbulk)    ! concentration H+ in bulk
xHplusbulk = (cHplus*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol
pOHbulk= pKw -pHbulk
cOHmin = 10**(-pOHbulk)   ! concentration OH- in bulk
xOHminbulk = (cOHmin*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol  

xsalt=(csalt*Na/(1.0d24))*(vsalt*vsol)   ! volume fraction salt,csalt in mol/l 

!!GGG!!!

if(rank.eq.0)print*,'nst',nst
do i = 0, adsmax
algo = MOD(i,2)
IF(algo.EQ.0)Tcapas(i) = 2 
IF(algo.NE.0)Tcapas(i) = 1 
end do


! eps

eps(ntot)=eps1
do i=1,ntot-1
eps(i)=0
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CHAIN GENERATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (cadenastype.eq.1) then
   if(rank.eq.1)print*, 'Calling RIS chain generator'
   elseif (cadenastype.eq.2) then
   if(rank.eq.1)print*, 'Calling MK chain generator'
   else
   stop 'Wrong chain generator'
   endif

in1n = 0

do LT = 1,2

   call initcha              ! init matrices for chain generation
   conf=0                    ! counter of number of conformations
   sumweight_tosend = 0.0

   do while (conf.lt.cuantas(LT))

   if (cadenastype.eq.1) then
   call cadenas1(chains,chainsw,ncha,LT)
   elseif (cadenastype.eq.2) then 
   call cadenas_MK(chains,chainsw,ncha,LT)
   endif

   do j=1,ncha

   if(conf.lt.cuantas(LT)) then
   conf=conf+1

   do ii = 1, ntot ! position of first segment
       weight(LT,conf,ii)=chainsw(j)

      minpos(LT,conf,ii) = ntot
      maxpos(LT,conf,ii) = 0 

      in1tmp = 0

      do k=1,long(LT)

      select case (abs(curvature))
      case (2)
        tempr=( (chains(1,k,j)+(float(ii)-0.5)*delta)**2 + chains(2,k,j)**2 +chains(3,k,j)**2 )**(0.5)
        temp=int(tempr/delta)+1  ! put them into the correct layer
      case (1)
        tempr=((chains(1,k,j)+(float(ii)-0.5)*delta)**2+chains(2,k,j)**2)**(0.5)
        temp=int(tempr/delta)+1  ! put them into the correct layer
      case (0) 
        tempr=(chains(1,k,j)+(float(ii)-0.5)*delta)
        temp=int(tempr/delta)+1  ! put them into the correct layer
      endselect

      if(curvature.lt.0) then   ! pore 
        if (temp.le.ntot) then            ! la cadena empieza en el layer 1
            in1tmp(k) = temp
        else
            weight(LT,conf,ii)=0.0 ! out of pore
        endif
      else ! convex
        if (temp.ge.radio) then            ! la cadena empieza en el layer 1
             in1tmp(k) = temp
        else
            weight(LT,conf,ii)=0.0 ! collide with the walls
        endif
      endif

       if(temp.lt.minpos(LT,conf,ii))minpos(LT,conf,ii)=temp
       if(temp.gt.maxpos(LT,conf,ii))maxpos(LT,conf,ii)=temp
 
       enddo ! k


       if(weight(LT,conf,ii).gt.0.0) then
  
       if((maxpos(LT,conf,ii)-minpos(LT,conf,ii)).ge.base) then
       print*,'Rank', rank, 'Increase base'
       call MPI_FINALIZE(ierr) ! finaliza MPI
       stop
       endif

       do k = 1, long(LT)
!       print*, k, in1tmp(k),minpos(LT,conf,ii) 
       temp = in1tmp(k)-minpos(LT,conf,ii)+1 
       in1n(LT,conf,ii,temp) =  in1n(LT,conf,ii,temp) + 1
       enddo

       endif

!   if((conf.eq.3).and.(LT.eq.1).and.(ii.eq.1)) then
!      do jj = 1, ntot
!      print*, 'oo', jj, in1n(LT,conf,ii,jj)
!      enddo
!   endif

   enddo ! ii
   sumweight_tosend = sumweight_tosend +  chainsw(j)
   endif

   enddo ! j
   enddo ! while

   tmp = 0.0
   call MPI_REDUCE(sumweight_tosend, tmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
   CALL MPI_BCAST(tmp, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,err)
   sumweight(LT) = tmp

enddo ! LT

!do LT = 1,2
!do i = 1, cuantas(LT)
!do ii = 1, ntot
!do j = 1, ntot
!print*, LT,i,ii,j,in1n(LT,i,ii,j)
!enddo
!enddo
!enddo
!enddo

!stop

!   do j = 1, ntot
!   print*, j, ii, in1n(1,3,1,j)
!   enddo
!   stop

!     end chains generation 
if(rank.eq.0)print*," chains ready"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     computation starts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     initializations of input depended variables 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


iter=0                    ! iteration counter
newcuantas(1)=cuantas(1)
newcuantas(2)=cuantas(2)

if(rank.eq.0) then
open(unit=533,file='ADS-mol.cm-2.dat')
open(unit=534,file='ADS-cad.nm-2.dat')
open(unit=535,file='meanz.dat')
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  MAIN LOOP OVER LAYERS 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

avpol = 0.0
! surface layer 
if (curvature.lt.0) then ! pore
avpol(0,n) = avpol(0,n) + sigma
else
avpol(0,radio) = avpol(0,radio) + sigma
endif

countfileuno=1           

do nads=preads, adsmax-1 ! number of layers, starts from preads

countfile=1


do cc = 1, nst !loop st
st = sts(cc)

do ccc = 1, nkbind !loop kbind
kbind = kbinds(ccc)

 123 LT = Tcapas(nads+1) ! type of current layer

! inits output files 
write(meanzfilename,'(A6,BZ,I5.5,A4)')'meanz.',countfile,'.dat'
write(sigmafilename,'(A6,BZ,I5.5,A4)')'sigma.',countfile,'.dat'
write(sigmaadfilename,'(A8,BZ,I5.5,A4)')'sigmaad.',countfile,'.dat'

Kbind0 = Kbind/(Na/1.0d24) ! Intrinsic equilibrium constant from uncharged polymers. Definida igual que en el paper JCP

K0A = (KaA*vsol)*(Na/1.0d24)! intrinstic equilibruim constant 
K0B = (Kw/KaB*vsol)*(Na/1.0d24)
K0ANa = (KaANa/vsol)/(Na/1.0d24)! Esta definida igual que en el paper JCP
K0BCl = (KaBCl/vsol)/(Na/1.0d24)! Esta definida igual que en el paper JCP


! Iteration to calculate bulk composition

xsolerror = 1e-6 ! maximum error in xsolbulk
xsoliter = 1.0   ! current iterate
xsolbulk = 100.0


do while (abs(xsoliter-xsolbulk).gt.xsolerror)

xsolbulk = xsoliter

!!!!!!!!! Calculate bulk composition !!!!!!!!!!!!!!!!!!!!!!!!!!

! 1) pH counterions
if(pHbulk.le.7) then  ! pH<= 7
   xposbulk=xsalt/zpos
   xnegbulk= -xsalt/zneg +(xHplusbulk -xOHminbulk) *vsalt ! NaCl+ HCl  
else                  ! pH >7 
   xposbulk=xsalt/zpos +(xOHminbulk -xHplusbulk) *vsalt ! NaCl+ NaOH   
   xnegbulk=-xsalt/zneg
endif

! 2) disociation in bulk  
select case (LT) ! 

case(1) ! polimero negativo

! variables auxiliares
  alfa=xposbulk/vsalt*K0ANa/(xsolbulk**vsalt)
  beta = K0A*xsolbulk/xHplusbulk

! calculo primero fs en bulk con los cationes
  fNchargebulk(1) = 1.0/(1.0+beta*(1.0+alfa))
  fionchargebulk(1)= alfa*beta*fNchargebulk(1)

! ahora, le resto el polimeri a los aniones
  xnegbulk = xnegbulk - (phibulkpol/(vpol*vsol)*(1.0-fNchargebulk(1)-fionchargebulk(1)))*(vsalt*vsol)

! finalmente, calculo el f bulk del polimero positivo con los aniones

  alfa=xnegbulk/vsalt*K0BCl/(xsolbulk**vsalt)
  beta = K0B*xsolbulk/xOHminbulk

  fNchargebulk(2) = 1.0/(1.0+beta*(1.0+alfa))
  fionchargebulk(2)= alfa*beta*fNchargebulk(2)

case(2) ! polimero positivo

! variables auxiliares
  alfa=xnegbulk/vsalt*K0BCl/(xsolbulk**vsalt)
  beta = K0B*xsolbulk/xOHminbulk

! calculo primero fs en bulk con los aniones
  fNchargebulk(2) = 1.0/(1.0+beta*(1.0+alfa))
  fionchargebulk(2)= alfa*beta*fNchargebulk(2)

! ahora, le resto el polimero a los cationes
  xposbulk = xposbulk - (phibulkpol/(vpol*vsol)*(1.0-fNchargebulk(2)-fionchargebulk(2)))*(vsalt*vsol)

! finalmente, calculo el f bulk del polimero negativo con los cationes
  alfa=xposbulk/vsalt*K0ANa/(xsolbulk**vsalt)
  beta = K0A*xsolbulk/xHplusbulk

  fNchargebulk(1) = 1.0/(1.0+beta*(1.0+alfa))
  fionchargebulk(1)= alfa*beta*fNchargebulk(1)

endselect

  xsoliter=1.0 -xHplusbulk -xOHminbulk - xnegbulk -xposbulk -phibulkpol

  if(rank.eq.0)print*,'Error de iteracion', abs(xsoliter-xsolbulk),xsoliter,xsolbulk 

enddo ! iter xsoliter

xsolbulk=xsoliter
! CHEQUEOS

Check_Kbminbulk=xOhminbulk/xsolbulk*(1.0-fNchargebulk(2)-fionchargebulk(2))/fNchargebulk(2)-K0B !!
Check_Kaplus= xHplusbulk/xsolbulk*(1.0-fNchargebulk(1)-fionchargebulk(1))/fNchargebulk(1)-K0A !!

check_KANa=(1.0-fNchargebulk(1)-fionchargebulk(1))*(xposbulk/vsalt)/(fionchargebulk(1)*(xsolbulk**vsalt))-1./K0ANA
check_KBCl=(1.0-fNchargebulk(2)-fionchargebulk(2))*(xnegbulk/vsalt)/(fionchargebulk(2)*(xsolbulk**vsalt))-1./K0BCl

if(rank.eq.0)print*,'Chequeos Equilibrios (Kb,Ka,KbCl,KaNa) ',Check_Kbminbulk,Check_Kaplus,Check_KBCl,Check_KANa

if(rank.eq.0)print*, 'Xsolbulk', xsolbulk
!stop
!!!!!!!! Charge in bulk !!!!!!!!!!!!!!!!!!!

if(rank.eq.0)print*, 'Charge in bulk in q/nm^3'
if(rank.eq.0)print*, '+',xposbulk/(vsol*vsalt)
if(rank.eq.0)print*, '-', xnegbulk/(vsol*vsalt)
if(rank.eq.0)print*, 'H+', xHplusbulk/vsol
if(rank.eq.0)print*, 'OH-', xOHminbulk/vsol
 
if(rank.eq.0) then
  if(LT.eq.1)print*, 'pol-',-phibulkpol/(vpol*vsol)*(1.0-fNchargebulk(1)-fionchargebulk(1))
  if(LT.eq.1)print*, 'f pol-',(1.0-fNchargebulk(1)-fionchargebulk(1))
  if(LT.eq.2)print*, 'pol+',phibulkpol/(vpol*vsol)*(1.0-fNchargebulk(2)-fionchargebulk(2))
  if(LT.eq.2)print*, 'f pol+',(1.0-fNchargebulk(2)-fionchargebulk(2))
  
  print*, 'sum:', xposbulk/(vsol*vsalt)+xHplusbulk/vsol-xnegbulk/(vsol*vsalt)- & 
  xOHminbulk/vsol-phibulkpol/(vpol*vsol)*(1.0-fNchargebulk(1)-fionchargebulk(1))
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
expmupos = xposbulk /xsolbulk**vsalt
expmuneg = xnegbulk /xsolbulk**vsalt
expmuHplus = xHplusbulk /xsolbulk   ! vsol = vHplus 
expmuOHmin = xOHminbulk /xsolbulk   ! vsol = vOHmin 
!!!

expmupol = (phibulkpol*vsol)/(vpol*vsol*long(LT)*xsolbulk**(long(LT)*vpol))        ! exp of bulk value of pol. chem. pot.
expmupol = expmupol/sumweight(LT)
expmupol = expmupol*(1.0-fNchargebulk(LT)-fionchargebulk(LT))**long(LT)


!print*,'expmu',expmupos,expmuneg,expmuhplus,expmuohmin,expmupol
!stop

!expmupol=expmupol/dexp(sumXu11*st/(vpol*vsol)*(phibulkpol)*long(LT))

!!!!
! solver
!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     init guess all 1.0 



if (infile.eq.0) then
do i=1,n
xg1(i)= xsolbulk
x1(i)=xg1(i)
zc(i)= (i-0.5) * delta
enddo

do i=n+1,2*n
xg1(i)=0.0
x1(i)=0.0
enddo

do i=2*n+1,3*n
xg1(i)=phibulkpol
x1(i)=phibulkpol
enddo
endif

!     init guess from files fort.100 (solvent) and fort.200 (potential)                      

if (infile.ge.1) then
do i=1,n
read(100,*)j,xfile(i)   ! solvent
enddo
do i=n+1,2*n
read(200,*)j,xfile(i)   ! potencial
enddo     
do i=2*n+1,3*n
read(300,*)j,xfile(i)   ! polimero
enddo     
endif

if(infile.eq.2) then
do ii=1, preads ! preadsorbed layers
do i = 1, ntot
read(300+ii, *)j, avpol(ii,i)
print*, i, ii, avpol(ii, i)
enddo
enddo
endif

!do i=1,3*n             ! initial gues for x1
!xg1(i)=x1(i)
!enddo

! JEFE
if(rank.eq.0) then ! solo el jefe llama al solver
   iter = 0
   print*, 'solve: Enter solver ', ntot, ' eqs'
   print*, 'KBIND', Kbind
   call call_kinsol(x1, xg1, ier)
   flagsolver = 0
   CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
endif
! Subordinados
if(rank.ne.0) then
  do
     flagsolver = 0
     source = 0
     CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
     if(flagsolver.eq.1) then
        call call_fkfun(x1) ! todavia no hay solucion => fkfun 
     endif ! flagsolver
     if(flagsolver.eq.0) exit ! Detiene el programa para este nodo
   enddo
endif

if (rank.eq.0) then
  avpol_red(:) = avpol(nads+1,:)
  CALL MPI_BCAST(avpol_red, ntot, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,err)
  CALL MPI_BCAST(norma, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,err)
endif
if(rank.ne.0) then
  CALL MPI_BCAST(avpol_red, ntot, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,err)
  CALL MPI_BCAST(norma, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,err)
  avpol(nads+1,:) = avpol_red(:)
endif


do i=1,n
xsol(i)=x1(i)
enddo

infile=-1

!if(norma.gt.error)stop

if(norma.gt.error) then
if(ccc.eq.1) then
Kbind = Kbind/2.0
goto 123
endif
if(rank.eq.0)print*, 'Fail', Kbind,ccc
Kbind=(kbinds(ccc-1)+Kbind)/2.0
if(rank.eq.0)print*, 'Try', Kbind,ccc
!kbinds(ccc-1) = Kbind
goto 123
endif


if(rank.eq.0)print*,"LAYER ", nads+1, " ADSORBED!", st, kbind, norma, error
xg1 = x1 ! work ok, save guess for next iteration

if(kbinds(ccc).ne.Kbind) then
Kbind = kbinds(ccc)
goto 123
endif

infile = 0 ! After the layer is adsorbed, reset initial guess


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determination of adsorbed polymer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! only for last Kbind

if (ccc.eq.nKbind) then

if(rank.eq.0) then
sumrho2 = 0.0
do i=1,n
sumrho2 = sumrho2 + avpol2(i)/vpol
enddo

sumrho2mol = sumrho2/vsol/Na*1.0d21*delta/1.0d7 ! mol.cm-2

sumrho2 = sumrho2/vsol*delta/long(LT)           ! chains by nm2

WRITE(533,*)(nads+1), sumrho2mol  ! mol.cm-2
WRITE(534,*)(nads+1), sumrho2  ! chains by nm2

! weighted z

meanz=0.0
sumrhoz=0.0

do i=1,n
do jj = 1, nads+1
meanz=meanz+(avpol(jj,i))*zc(i)
sumrhoz=sumrhoz+avpol(jj,i)
end do
!meanz=meanz+(avpol2(i))*zc(i) 
!sumrhoz=sumrhoz+avpol2(i)
enddo
meanz=meanz/sumrhoz

endif ! only for last Kbind


write(535,*)(nads+1), meanz
write(sysfilename,'(A7,BZ,I3.3,A1,I3.3,A4)')'system.', countfileuno,'.',countfile,'.dat'
write(denspol2filename,'(A16,BZ,I3.3,A1,I3.3,A4)')'densitypolymerA.',countfileuno,'.',countfile,'.dat'
write(denssolfilename,'(A15,BZ,I3.3,A1,I3.3,A4)')'densitysolvent.', countfileuno,'.',countfile,'.dat'
write(totalfilename,'(A13,BZ,I3.3,A1,I3.3,A4)')'densitytotal.',countfileuno,'.',countfile,'.dat'

write(xposfilename,'(A10,BZ,I3.3,A1,I3.3,A4)')'xpostotal.',countfileuno,'.',countfile,'.dat'
write(xnegfilename,'(A10,BZ,I3.3,A1,I3.3,A4)')'xnegtotal.',countfileuno,'.',countfile,'.dat'
write(xHplusfilename,'(A12,BZ,I3.3,A1,I3.3,A4)')'xHplustotal.',countfileuno,'.',countfile,'.dat'
write(xOHminfilename,'(A12,BZ,I3.3,A1,I3.3,A4)')'xOHmintotal.',countfileuno,'.',countfile,'.dat'
write(potenfilename,'(A6,BZ,I3.3,A1,I3.3,A4)')'poten.',countfileuno,'.',countfile,'.dat'
write(qtotfilename,'(A6,BZ,I3.3,A1,I3.3,A4)')'q_tot.',countfileuno,'.',countfile,'.dat'
write(rho2filename,'(A5,BZ,I3.3,A1,I3.3,A4)')'rho2.',countfileuno,'.',countfile,'.dat'


do jj = 1, nads+1
write(denspolfilename(jj),'(A14,BZ,I2.2,A1,I3.3,A1,I3.3,A4)')'densitypolymer',jj,'.',countfileuno,'.',countfile,'.dat'
end do

do ii = 1, 2
write(densbindfilename(ii),'(A11,BZ,I2.2, A1, I3.3,A1,I3.3,A4)')'densitybind',ii,'.', countfileuno,'.',countfile,'.dat'
end do
do ii = 1, 2
write(denschargefilename(ii),'(A12,BZ,I2.2, A1, I3.3,A1,I3.3,A4)')'densitycharge',ii,'.', countfileuno,'.',countfile,'.dat'
end do


open(unit=310,file=sysfilename)
open(unit=321,file=denspol2filename)
open(unit=323,file=totalfilename)
open(unit=330,file=denssolfilename)

open(unit=331,file=xposfilename)
open(unit=332,file=xnegfilename)
open(unit=333,file=xHplusfilename)
open(unit=334,file=xOHminfilename)
open(unit=335,file=potenfilename)
open(unit=336,file=qtotfilename)
open(unit=337,file=rho2filename)


do ii = 1,2
open(unit=500+ii+5,file=densbindfilename(ii))
open(unit=700+ii+5,file=denschargefilename(ii))
end do

do jj = 1, nads+1
open(unit=600+jj,file=denspolfilename(jj))
end do

do i=1,n
write(321,*)zc(i),avpol2(i)
write(330,*)zc(i),xsol(i)

do ii = 1,2
write(500+ii+5,*)zc(i),fbound(ii, i)
write(700+ii+5,*)zc(i),1.-fNcharge(ii, i)-fbound(ii,i)

end do

do jj = 1, nads+1
write(600+jj,*)zc(i),avpol(jj,i)
avtotal(i) = avpol(jj,i)
end do

WRITE(323,*)zc(i),(avtotal(i))
enddo


do i=1,n
write(331,*)zc(i),xpos(i)
write(332,*)zc(i),xneg(i)
write(333,*)zc(i),xHplus(i)
write(334,*)zc(i),xOHmin(i)
write(335,*)zc(i),psi(i)
write(336,*)zc(i),qtot(i)
write(337,*)zc(i),rhopol2(i)
end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     additional system information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(310,*)'system      = neutral polymer'
write(310,*)'fnorm       = ', norma ! residual size of iteration vector
write(310,*)'error       = ',error
write(310,*)'q           = ',q
write(310,*)'length seg  = ',0.50 ! value see subroutine cadenas
write(310,*)'delta       = ',delta
write(310,*)'vsol        = ',vsol
write(310,*)'vpol        = ',vpol*vsol

write(310,*)'Kbind       = ',Kbind
write(310,*)'Kbind0      = ',Kbind0

write(310,*)'long1       = ',long1
write(310,*)'long2       = ',long2
write(310,*)'radio       = ',radio
write(310,*)'sigma       = ',sigma

write(310,*)'cuantas(1)     = ',cuantas(1)
write(310,*)'cuantas(2)     = ',cuantas(2)

write(310,*)'iterations  = ',iter

WRITE(310,*)'Tipo de capa=', LT

close(310)
close(320)
CLOSE(321)
CLOSE(322)
CLOSE(323)
close(330)
close(350)
close(360)
close(370)
close(380)
close(390)
close(400)
close(410)
close(411)
close(412)
close(420)

do ii=1,2
CLOSE(500+ii)
CLOSE(505+ii)
end do

do jj=1, nads+1
CLOSE(600+ii)
CLOSE(700+ii)
CLOSE(800+ii)
end do

countfile = countfile+1 ! next

endif ! rank
call fe(cc)

END do ! loop de kbind
end do ! loop de st

countfileuno = countfileuno + 1

END do ! loop de numero de paso de adsorcion

close(530)
close(531)
CLOSE(532)
CLOSE(533)
CLOSE(534)
CLOSE(535)
CLOSE(536)
CLOSE(537)

call MPI_FINALIZE(ierr) ! finaliza MPI
stop

end

