subroutine allocation
use multicapa
use mkinsol
use posmk
use longs

allocate (Tcapas(adsmax))
allocate (avpol(adsmax,ntot,2))
allocate (avpol2(ntot,2))
allocate (avpolall(ntot,2))
allocate (in1n(2,maxcuantas,ntot,2))
allocate (maxlayer(adsmax,maxcuantas))
allocate (eps(2*ntot))
allocate (xtotal(2*ntot,2))
allocate (fbound(2,2*ntot))
allocate (pp(2*ntot))
allocate (Xu(2,2,-Xulimit:Xulimit))
allocate (weight(2,maxcuantas))
allocate (current(maxlong,3))
allocate (nextbead(maxlong))
ALLOCATE (firstcell(-mcube:mcube,-mcube:mcube,-mcube:mcube))
end
