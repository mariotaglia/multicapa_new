subroutine allocation
use multicapa
use mkinsol
use posmk
use longs

allocate (Tcapas(adsmax))
allocate (avpol(adsmax,ntot))
allocate (rhopol(adsmax,ntot))
allocate (avpol2(ntot))
allocate (rhopol2(ntot))
allocate (avpolall(ntot))
allocate (in1n(2,maxcuantas,ntot))
allocate (maxlayer(adsmax,maxcuantas))
allocate (eps(2*ntot))
allocate (xtotal(2*ntot))
allocate (fbound(2,2*ntot))
allocate (pp(2*ntot))
allocate (Xu(2,2,-Xulimit:Xulimit))
allocate (weight(2,maxcuantas))
allocate (current(maxlong,3))
allocate (nextbead(maxlong))
ALLOCATE (firstcell(-mcube:mcube,-mcube:mcube,-mcube:mcube))
end
