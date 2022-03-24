subroutine read
use longs
use multicapa
use bulk
use const
use MPI
use volume
implicit none
integer i

!     reading in of variables from stdin

read(8,*),nada
read(8,*),curvature

read(8,*),nada
read(8,*),radio

read(8,*),nada
read(8,*),adsmax

read(8,*),nada
read(8,*),ntot

read(8,*),nada
read(8,*)cadenastype

read(8,*),nada
read(8,*)cuantas1, cuantas2
maxcuantas=max(cuantas1/size,cuantas2/size)

read(8,*),nada
read(8,*)long1, long2
maxlong = MAX(long1, long2)

read(8,*),nada
read(8,*)minn

read(8,*),nada
read(8,*),npbp

do i =1, npbp
read(8,*),pbps(i)
end do

read(8, *), nada
read(8, *), pKaA    ! pKaA of weak polyacid segments

read(8, *), nada
read(8, *), pKaB    ! pKaB of weak polyacid segments

read(8, *), nada
read(8, *), pKaANa    ! pKaA of weak polyacid segments

read(8, *), nada
read(8, *), pKaBCl   ! pKaB of weak polyacid segments

read(8, *), nada
read(8, *), Ksal

read(8, *), nada
read(8, *), csalt  ! salt concentration in bulk (Molar)

read(8, *), nada
read(8, *), pHbulk ! bulk pH

READ(8,*),nada
read(8,*),nst
do i =1, nst
read(8,*),sts(i)
end do

READ(8,*),nada
READ(8,*), preads

READ(8,*),nada
read(8,*),sigma

READ(8,*),nada
read(8,*),error

READ(8,*),nada
read(8,*),infile

READ(8,*),nada
READ(8,*), kbind

read(8,*)nada
read(8,*)eps1

read(8,*)nada
read(8, *)Xulimit

read(8,*)nada
read(8,*)lseg

read(8,*)nada
read(8,*)AA,BA,CA

read(8,*)nada
read(8,*)vpol0

read(8,*)nada
read(8,*)vsaltpos, vsaltneg



if(curvature.eq.0)radio=1 ! planar
if(curvature.lt.0)ntot=radio ! pore


end
