!#####################################################################
!
! Este programa calcula los kai para poor-solvent en 1D
!
!#####################################################################

subroutine kais
use multicapa
use MPI
use layer
implicit none

real*8 l ! medio lseg, radio0 del segmento
integer MCsteps ! numero de steps de MC
integer iz
real*8 x,y,z, radio0
integer limit
real*8,allocatable :: matriz11(:)
real*8,allocatable :: matriz12(:)
integer i
real*8 rands,R
real*8 suma11
real*8 suma12
real*8 cutoff
integer jx,jy,jz,dimz,ir
real*8 LJ
real*8, external :: interaction11, interaction12
dimz=1
limit = Xulimit +1
cutoff = float(Xulimit)*delta

ALLOCATE (matriz11(-limit:limit)) ! matriz de kai
ALLOCATE (matriz12(-limit:limit)) ! matriz de kai
if(rank.eq.0)print*,'kais: Kai calculation'
suma11 = 0.0
suma12 = 0.0
do iz = -limit, limit
matriz11(iz) = 0.0
matriz12(iz) = 0.0
enddo

MCsteps = 200
l = lseg 

do jx = 1, MCsteps
do jy = 1, MCsteps
do jz = 1, MCsteps

x = 2.0*cutoff*(((float(jx)-0.5)/float(MCsteps))-0.5) ! number between -cutoff and +cutoff
y = 2.0*cutoff*(((float(jy)-0.5)/float(MCsteps))-0.5)
z = 2.0*cutoff*(((float(jz)-0.5)/float(MCsteps))-0.5)

radio0 = sqrt(x**2 + y**2 + z**2) ! espacio real

   select case (abs(curvature))
   case (0)
   R = abs(x)
   Z = z
   case (1)
   R = sqrt(x**2 + y**2)
   Z = z
   case (2)
   R = sqrt(x**2 + y**2 + z**2)
   if(dimZ.ne.1)stop
   Z = 0.0
   end select

if(radio0.gt.cutoff) cycle ! No esta dentro de la esfera del cut-off   
if(radio0.lt.l) cycle ! esta dentro de la esfera del segmento

 ! celda 
 iR = int(R/delta)+1! 
 iz = int(anint(z/delta))

if(iz.le.ntot)then
 matriz11(iz) = matriz11(iz) + interaction11(radio0)
 matriz12(iz) = matriz12(iz) + interaction12(radio0)
endif
enddo
enddo
enddo

sumXu11 = 0.0
sumXu12 = 0.0
do iz = -Xulimit, Xulimit
 Xu(1,1,iz) = matriz11(iz)/(MCsteps**3)*((2.0*cutoff)**3)
 Xu(2,2,iz) = matriz11(iz)/(MCsteps**3)*((2.0*cutoff)**3)
 Xu(1,2,iz) = matriz12(iz)/(MCsteps**3)*((2.0*cutoff)**3)
 Xu(2,1,iz) = matriz12(iz)/(MCsteps**3)*((2.0*cutoff)**3)
 suma11 = suma11 +  matriz11(iz)/(MCsteps**3)*((2.0*cutoff)**3)
 suma12 = suma12 +  matriz12(iz)/(MCsteps**3)*((2.0*cutoff)**3)
 sumXu11 = sumXu11 + Xu(1,1,iz)
 sumXu12 = sumXu12 + Xu(1,2,iz)
enddo

if(rank.eq.0)print*, 'kais: Sum Xulimit11', suma11
if(rank.eq.0)print*, 'kais: Sum Xulimit12', suma12

suma11 = 0.0
suma12 = 0.0
do iz = -limit, limit
 suma11 = suma11 +  matriz11(iz)/(MCsteps**3)*((2.0*cutoff)**3)
 suma12 = suma12 +  matriz12(iz)/(MCsteps**3)*((2.0*cutoff)**3)
enddo

do iz = -Xulimit, Xulimit
if(rank.eq.0)print*, 'Xu',iz,Xu(1,1,iz),Xu(1,2,iz)
enddo

if(rank.eq.0)print*, 'kais: Total Sum', suma11,suma12

end

double precision function interaction11(d)
use multicapa
implicit none
real*8 d
interaction11 = (lseg/d)**(6.0)
end function

double precision function interaction12(d)
use multicapa
implicit none
real*8 d
interaction12 =  (lseg/d)**(6.0)
end function



