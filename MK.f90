
subroutine cadenas_mk(chains,chainsw,ncha,LT)
! =====================================================================================================
! version 21 may 2011 mk@mat.ethz.ch
! =====================================================================================================

! =====================================================================================================
! cadenas requires ncha = -1 before it is called for the first time! 
! =====================================================================================================

! creates ncha polymers of length N confinded via matriz >= 0
! returns bead positions as cell numbers (chains)
! bond length: lseg, bond diameter: dseg (the latter used in overlap check)
! chains(1..ncha,1..long) = matriz(iR,iZ)

! matriz(iR,iZ) = cell number (outside geometry: < 0)   ! pore.system.h
! indexa(cell number,1) = iR                            ! lookup.h
! indexa(cell number,2) = iZ

use longs
use seed1
use pis
use matrices
use senos
use multicapa
use mpi
use posmk

implicit none

real*8 x(3)
integer ncha
real*8 chains(3,maxlong,ncha_max),chainsw(ncha_max),N
integer ncha_current
common /comncha/ ncha_current
integer N_max
common /compass/ N_max
integer spacer,branch_p
common /com1/ spacer,branch_p
integer LT
real*8 dseg

dseg = lseg

if(LT.eq.1)N=long1
if(LT.eq.2)N=long2


! =====================================================================================================
branch_p      = 5               ! EVENTUALLY SYSTEM AND MODEL SPECIFIC
spacer        = 2               ! BUT q IS AUTOMATICALLY ADJUSTED TO CHOSEN s and p VIA 'TUNING'
! =====================================================================================================

! =====================================================================================================
x(1) = 0.0
x(2) = 0.0
x(3) = 0.0
! =====================================================================================================
! =============================== DO NOT EDIT BELOW THIS LINE =========================================
! =====================================================================================================

ncha_current    = 0
N_max           = N

call walk(1,x,chains,chainsw)
ncha = ncha_current
return
end

! =========================================================================== 

recursive subroutine walk(N,x,chains,chainsw)
use longs
use seed1
use pis
use matrices
use senos
use multicapa
use mpi
use posmk
use layer
implicit none

real*8 x(3),y(3),u(3)
integer N
integer ncha_current
common /comncha/ ncha_current
integer*2 beadno
integer   j,k,ix1,ix2,ix3,iR,iZ,ibranch
real*8    dummy,stretched,power
logical   hit_bead
external  hit_bead
integer N_max
common /compass/ N_max
real*8 chains(3,maxlong,ncha_max), chainsw(ncha_max)
real*8 endtoendtemp(10000)
!        integer*4, allocatable ::  deathcount(:)
integer*4 deathtotal
integer spacer,branch_p
common /com1/ spacer,branch_p
common /crit/ deathtotal
common /endtoend/ endtoendtemp
save power
real*8 tmp, dista
!        save deathcount
real*8, external :: rands
!        if(.not.allocated(deathcount))ALLOCATE(deathcount(maxlong)) 
real*8, external :: interaction11
real*8 dseg
integer k1,k2

dseg = lseg

if (N.eq.1) then
!print*, 3,N,'root',ncha_current
stretched   = (N_max-1)*lseg                           ! EVENTUALLY MODEL SPECIFIC
power = dlog(dble(mcube))/dlog(stretched/dseg)
power = min(1.D0,power)
ncha_current = 0
firstcell    = 0
nextbead     = 0
beadno       = 0
nextbead(N)  = 0
firstcell(0,0,0) = N
!         deathcount   = 0
deathtotal   = 0
! print 1,'long',long
! print 1,'N',N
! print 2,'qprob0',qprob0
! print 1,'spacer',spacer
! print 1,'branch_p',branch_p
! print 1,'nearbonds',nearbonds
! print 1,'ncha_current',ncha_current
do k=1,3; current(1,k) = x(k); enddo
!if (outside(x)) then
!print *,x
!print *,int(sqrt(x(1)**2+x(2)**2)/delta),int(x(3)/delta)
!stop 'anchor is outside!?'
!end if
111      continue
call randomunit(u)
do k=1,3; y(k) = x(k) + lseg*u(k); enddo
!if (outside(y)) goto 111
call walk(2,y,chains,chainsw)          ! calls N=2, next gen is 2 (after some strand ..)
return
end if

if (rands(seed).gt.qprob0) then
! print 3,N,'q prob return',ncha_current
goto 222
return
endif



!if (ncha_current.ge.ncha_max)   return 
!if (outside(x)) then
! print 3,N,'outside return',ncha_current
!goto 222
!return
!endif

dummy = x(1)-current(1,1)
ix1 = nint(dexp(power*dlog(dabs(dummy)/delta))*sign(1.0D0,dummy))
dummy = x(2)-current(1,2)
ix2 = nint(dexp(power*dlog(dabs(dummy)/delta))*sign(1.0D0,dummy))
dummy = x(3)-current(1,3)
ix3 = nint(dexp(power*dlog(dabs(dummy)/delta))*sign(1.0D0,dummy))

!        print*, 'firstcell is', size(firstcell, 1), size(firstcell, 2),
!     & size(firstcell, 3)

if (hit_bead(ix1,ix2,ix3,x,N)) then
!print*, 3,N,'hitbead return',ncha_current
goto 222
end if

beadno                   = firstcell(ix1,ix2,ix3)
nextbead(N)              = beadno
firstcell(ix1,ix2,ix3)   = N

do k=1,3; current(N,k) = x(k); enddo

if (N.eq.N_max) then                   ! complete chain generated successfully

ncha_current = ncha_current + 1


!print*, 3,N,'chain end ---',ncha_current

! determine statistical weight
tmp = 1.0

do k1=1, N
do k2=k1+1,N
dista=(current(k2,1)-current(k1,1))**(2.0)
dista=dista+(current(k2,2)-current(k1,2))**(2.0)
dista=dista+(current(k2,3)-current(k1,3))**(2.0)
dista=sqrt(dista)
tmp=tmp*exp(interaction11(dista))
enddo
enddo

do j=1,N
 chains(1,j,ncha_current)=current(j,1)       ! output 
 chains(2,j,ncha_current)=current(j,2)
 chains(3,j,ncha_current)=current(j,3)
enddo

 chainsw(ncha_current)=tmp

!do j=1,N
!iR = int(sqrt(current(j,1)**2+current(j,2)**2)/delta)+1
!iZ = int(current(j,3)/delta)+1
!if ((iR.gt.dimR).or.(iZ.lt.1).or.(iZ.gt.dimZ)) then
!print*, 'Error in creador'
!print*, 'Increase system size'
!endif
!chains(ncha_current,j) = matriz(iR,iZ)
!enddo

!endtoendtemp(ncha_current) = 
!& ((current(N,1)-current(1,1))**2 +
!& (current(N,2)-current(1,2))**2 +
!& (current(N,3)-current(1,3))**2)**(0.5)

! call checking_actual_config(N)       ! DEACTIVE IN PRODUCTION

if (ncha_current.ge.ncha_max) return

else                                    ! N < N_max

if (mod(N-1,spacer).eq.1) then         ! N=2,2+spacer,...
! print 3,N,'branching',ncha_current
do ibranch=1,branch_p
call randomunit(u)
do k=1,3; y(k) = x(k) + lseg*u(k); enddo
call walk(N+1,y,chains,chainsw)              
if (ncha_current.ge.ncha_max) return
enddo
else
! print 3,N,'linear growth',ncha_current
call randomunit(u)
do k=1,3; y(k) = x(k) + lseg*u(k); enddo
call walk(N+1,y,chains,chainsw)     
if (ncha_current.ge.ncha_max) return
end if

end if

firstcell(ix1,ix2,ix3) = beadno         ! set free
return

1       format("walk ",A30,I10)
2       format("walk ",A30,F10.3)
3       format(I5,1x,A20,I10)

222     continue
!        deathcount(N) = deathcount(N) + 1
deathtotal    = deathtotal + 1

return
end

subroutine randomunit(seg)                      
use seed1
implicit none
real*8, external :: rands
real*8 seg(3),znorm,z(2),mysqrt
znorm=2.D0
do while (znorm.ge.1.D0)
z(1) = 1.D0-2.D0*rands(seed)
z(2) = 1.D0-2.D0*rands(seed)
znorm = z(1)*z(1)+z(2)*z(2)
enddo
mysqrt = dsqrt(1.D0-znorm)
seg(1) = 2.D0*z(1)*mysqrt
seg(2) = 2.D0*z(2)*mysqrt
seg(3) = 1.D0-2.D0*znorm
return
end

function hit_bead(ix1,ix2,ix3,x,N)
use posmk 
use multicapa
implicit none

real*8 tolerance
logical hit_bead
integer ix1,ix2,ix3,k,N,k1,k2,k3
real*8 x(3)
integer*2 beadno
real*8 isd2
real*8 d2, b2

tolerance = 1e-3

d2 = lseg**2 + tolerance
b2 = lseg**2 + tolerance

hit_bead = .false.
do k1=-1,1
do k2=-1,1
do k3=-1,1
!         print*, ix1, ix2, ix3, size(firstcell, 1), 
!     & size(firstcell, 2), size(firstcell, 3)
beadno = firstcell(ix1+k1,ix2+k2,ix3+k3)
if (beadno.ge.N) stop 'BAD'
do while (beadno.gt.0)
isd2 = (current(beadno,1)-x(1))**2+(current(beadno,2)-x(2))**2+(current(beadno,3)-x(3))**2
if (isd2.lt.d2) then                                  ! far away sphere have diameter d
if (abs(beadno-N).le.1) then                         ! treat adjacent as nonoverlapping 23 MAY 2011
hit_bead=.false.
elseif (abs(beadno-N).le.nearbonds) then
if (isd2.lt.b2) then                                ! near spheres have diameter b
hit_bead=.true.;  return
end if
else
hit_bead=.true.;  return
end if
end if
beadno = nextbead(beadno)
enddo
enddo
enddo
enddo
return
end






