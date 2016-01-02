subroutine sphere(radius,sph, dc) ! radius: radius of the sph in nm, sphere: vector distribution, must be preallocated
use layer

implicit none
integer dc
real*8 sph(dc)
integer MCsteps ! numero de steps de MC
integer ix, iy , iz
real*8 x,y,z, radio
integer i
real*8 rands
real*8 suma
real*8 radius
integer im
integer iix,iiy,iiz
real*8 vol

print*, 'sphere: init sphere calculation'

suma = 0.0
sph = 0.0

MCsteps = 200

do iix = 1, MCsteps
do iiy = 1, MCsteps
do iiz = 1, MCsteps

x = 2.0*radius*((dfloat(iix-1)/dfloat(MCsteps))-0.5)
y = 2.0*radius*((dfloat(iiy-1)/dfloat(MCsteps))-0.5)
z = 2.0*radius*((dfloat(iiz-1)/dfloat(MCsteps))-0.5)

radio = sqrt(x**2 + y**2 + z**2) ! espacio real

if (radio.gt.radius) cycle ! outside sph

 ! celda 
 iz = int((z+radius)/delta)+1
 sph(iz) = sph(iz) + 1.0
 suma = suma + 1.0

enddo !ix
enddo !iy
enddo !iz
print*, suma, radius, dc, sph
vol = 4.0/3.0*3.14159*(radius**3)
sph = sph/suma*vol ! normalize volume

end



