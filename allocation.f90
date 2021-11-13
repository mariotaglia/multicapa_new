subroutine allocation
use multicapa
use mkinsol
use posmk
use longs

allocate (Tcapas(0:adsmax))
allocate (avpol(0:adsmax,ntot))
allocate (avpol2(ntot))
allocate (avpolall(ntot))
allocate (xsol(ntot))  !!!!G
allocate (avpolneg(ntot))!G
allocate (avpolpos(ntot))!G
allocate (rhopol2(ntot))!G
allocate (avpolnegcero(ntot))
allocate (avpolposcero(ntot))
allocate (pro(maxcuantas)) !!!G
allocate (sumprolnpro(ntot))
allocate (in1n(2,maxcuantas,ntot,base))
allocate (in1tmp(maxlong))
allocate (maxpos(2,maxcuantas,2*ntot))
allocate (minpos(2,maxcuantas,2*ntot))
allocate (eps(ntot))
allocate (xtotal(ntot))
allocate (fbound(2,2*ntot))
allocate (fNcharge(2,2*ntot))
allocate (fioncharge(2,2*ntot))
allocate (pp(ntot*3))
allocate (Xu(2,2,-Xulimit:Xulimit))
allocate (weight(2,maxcuantas,ntot))
allocate (current(maxlong,3))
allocate (nextbead(maxlong))
ALLOCATE (firstcell(-mcube:mcube,-mcube:mcube,-mcube:mcube))

allocate(psi(ntot))
allocate(xpos(ntot))
allocate(xNaCl(ntot))
allocate(xneg(ntot))
allocate(xHplus(ntot))
allocate(xOHmin(ntot))
allocate(qtot(ntot))

end
