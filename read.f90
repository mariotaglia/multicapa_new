subroutine read
use longs
use multicapa
use bulk
implicit none
integer i

!     reading in of variables from stdin
read(8,*),nada
read(8,*),adsmax

read(8,*),nada
read(8,*),ntot

read(8,*),nada
read(8,*)cadenastype

read(8,*),nada
read(8,*)cuantas1, cuantas2
maxcuantas = MAX(cuantas1, cuantas2)

read(8,*),nada
read(8,*)long1, long2
maxlong = MAX(long1, long2)

read(8,*),nada
read(8,*)minn

READ(8,*),nada
READ(8,*),phibulkpol

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
READ(8,*), nkbind

do i =1, nkbind
read(8,*),kbinds(i)
end do

read(8,*)nada
read(8,*)eps1

read(8,*)nada
read(8, *)Xulimit

read(8,*)nada
read(8,*)lseg

read(8,*)nada
read(8,*)AA,BA,CA

end
