subroutine allocation
use multicapa
use mkinsol
use posmk
use longs

allocate (Tcapas(adsmax))
allocate (avpol(adsmax,ntot))
allocate (avpol2(ntot))
allocate (avpolall(ntot))
allocate (in1n(2,maxcuantas,ntot,ntot))
!allocate (maxlayer(adsmax,maxcuantas))
allocate (eps(ntot))
allocate (xtotal(ntot))
allocate (fbound(2,ntot))
allocate (pp(ntot))
!allocate (Xu(2,2,-Xulimit:Xulimit))
allocate (weight(2,maxcuantas,ntot))
allocate (current(maxlong,3))
allocate (nextbead(maxlong))
ALLOCATE (firstcell(-mcube:mcube,-mcube:mcube,-mcube:mcube))
end
