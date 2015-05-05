subroutine allocation
use multicapa
use mkinsol
allocate (Tcapas(adsmax))
allocate (avpol(adsmax,ntot))
allocate (avpol2(ntot))
allocate (avpolall(ntot))
allocate (in1n(2,maxcuantas,ntot))
allocate (maxlayer(adsmax,maxcuantas))
allocate (eps(2*ntot))
allocate (xtotal(2*ntot))
allocate (fbound(2,2*ntot))
allocate (pp(ntot))
end
