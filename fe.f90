subroutine fe(cc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! this routine calculates the free energy of the system and chemical potential 
! of chains
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!use results
!use chainsdat
use bulk
use partfunc
use layer
use multicapa
use volume
use longs
use MPI
implicit none

integer cc
integer n
real*8 Free_energy,Free_energy2, Fmupol
real*8 F_mix_s,F_mix_avpolA,F_mix_avpolb,F_conf,F_EQ,Fpro
integer i,ii,j,jj
real*8 sumas,sumrho,sumrhopol,sumpi,sum
real*8 Fact_rhobulk
integer iz,ir,dimz
double precision, external :: jacobian
integer sumnewcuantas
integer err

call MPI_REDUCE(newcuantas(AT), sumnewcuantas, 1, MPI_INTEGER, MPI_SUM,0, MPI_COMM_WORLD, err)

dimz=1

F_mix_s= 0.0
n=ntot

!!!!!!!!! F mix s !!!!!!!!!!!!!

Fact_rhobulk=phibulkpol/(float(long(AT))*vpol*vsol)
do iR =1, n
 F_mix_s=F_mix_s+xsol(iR)*(log(xsol(iR))-1.0)*jacobian(iR)
 F_mix_s=F_mix_s-xsolbulk*(log(xsolbulk)-1.0)*jacobian(ir)
enddo

F_mix_s=F_mix_s *delta /vsol 
if(rank.eq.0)print*,'fmixs',F_mix_s

Free_energy= Free_energy +F_mix_s

!!!!!!! F mix A !!!!!!!!!!!!!!!

F_mix_avpolA=0.0
do ir=1,n
if (rhopol2(ir).ne.0.0)then
 F_mix_avpolA=F_mix_avpolA+(rhopol2(ir))*(log(rhopol2(ir)*vsol)-1.0) *jacobian(ir)
endif
 F_mix_avpolA=F_mix_avpolA-(Fact_rhobulk)*(log(Fact_rhobulk*vsol)-1.0) *jacobian(ir)
enddo
F_mix_avpolA=F_mix_avpolA *delta

Free_energy= Free_energy +F_mix_avpolA

if(rank.eq.0)print*,'fmixA',F_mix_avpolA

!!!!!  F mupol !!!!!!!!!!!!!!!!!!

Fmupol=0.0
do ir=1,n
 Fmupol=Fmupol-rhopol2(ir)*log(expmupol)*jacobian(ir)
 Fmupol=Fmupol+Fact_rhobulk*log(expmupol)*jacobian(ir)
enddo
Fmupol=Fmupol*delta !-delta

Free_energy= Free_energy +Fmupol

if(rank.eq.0)print*,'fmupol',Fmupol

!!!!!! Fconf !!!!!!!!!!!!!!!!!!!!!

F_Conf=0
!!!!!!
do iR=1,n
if (rhopol2(ir).ne.0.0)then
 F_Conf=F_Conf+rhopol2(iR)*sumprolnpro(iR)*jacobian(iR)
endif
 F_Conf=F_Conf-Fact_rhobulk*log(1.0/float(sumnewcuantas)) *jacobian(iR)
enddo
!Linea  82 es porque  tenemos en el bulk pol?
F_conf=F_conf*delta  
if(rank.eq.0)print*,'F_conf',F_conf
!!!!!!

Free_energy=Free_energy+F_conf

F_EQ=0.
!!!!   avpolpos==<na> y avpolneg == <nb> 
do iR=1, n

F_EQ=F_EQ+avpolpos(iR)*fbound(1,iR)*(-log(Kbind0)) *jacobian(iR)
!F_EQ=F_EQ+avpolpos(i)*fbound(1,i)*(-log( fbound(2,i)/(1.0-fbound(1,i))/(1.0-fbound(2,i))/(avpolpos(i)/vsol/vpol) )) !chequear signo
F_EQ=F_EQ+avpolpos(iR)*(1.0-fbound(1,iR))*(log(1.0-fbound(1,iR))) *jacobian(iR)

if (fbound(1,iR)>0.0)then
  F_eq=F_eq+avpolpos(iR)*(fbound(1,iR))*(log(fbound(1,iR))) *jacobian(iR)
  F_eq=F_eq-avpolpos(iR)*fbound(1,iR)*(log(avpolpos(iR)/vpol/vsol*fbound(1,iR))-1.0)*jacobian(iR) !vab =1 <= no, ojo, aca va vpol
endif
enddo

if(rank.eq.0)print*,'Feq',F_eq*delta/(vpol*vsol) 

do iR = 1,n
if (fbound(2,iR)>0.)then
 F_eq=F_eq + avpolneg(iR)*fbound(2,iR)*log(fbound(2,iR)) *jacobian(iR)
 F_eq=F_eq + avpolneg(iR)*(1.-fbound(2,iR))*log(1.-fbound(2,iR)) *jacobian(iR) !<= ojo, habia error de parentesis
endif
enddo

if(rank.eq.0)print*,'Feq',F_eq*delta/(vpol*vsol) 


F_EQ=F_EQ*delta/(vpol*vsol) 


Free_energy=Free_energy+F_EQ

Free_energy2=0.0
sumpi=0.0
sumrho=0.0
sumrhopol=0.0
sumas=0.0
do iR=1, n
! if (xsol(i)>0.)then
 sumpi= sumpi + log(xsol(iR))*(1.0-avpolnegcero(iR)-avpolposcero(iR)) *jacobian(iR)
 sumpi= sumpi - log(xsolbulk)*jacobian(iR)
 sumrho= sumrho+(-xsol(iR))*jacobian(iR)
 sumrho= sumrho-(-xsolbulk)*jacobian(iR)

 sumrhopol=sumrhopol-(rhopol2(iR)) *jacobian(iR)
 sumrhopol=sumrhopol-(-Fact_rhobulk)*jacobian(iR)
if (AT.eq.1) then
 sumas=sumas+avpolpos(iR)*fbound(1,iR)*jacobian(iR)
 sumas=sumas+avpolneg(iR)*log(1.-fbound(2,iR))*jacobian(iR)
 sumas=sumas+avpolposcero(iR)*log(1.-fbound(1,iR))*jacobian(iR)
else
 sumas=sumas+avpolneg(iR)*fbound(2,iR)*jacobian(iR)
 sumas=sumas+avpolpos(iR)*log(1.-fbound(1,iR))*jacobian(iR)
 sumas=sumas+avpolnegcero(iR)*log(1.-fbound(2,iR))*jacobian(iR)
endif
enddo

sumpi=sumpi*delta/vsol!
sumrho=sumrho*delta/vsol!
sumrhopol=sumrhopol*delta!
sumas=sumas*delta/(vpol*vsol)!



if(rank.eq.0)print*,'sumpi',sumpi
if(rank.eq.0)print*,'sumrho',sumrho
if(rank.eq.0)print*,'sumrhopol',sumrhopol
if(rank.eq.0)print*,'sumas',sumas

sum=sumpi+sumrho+sumrhopol+sumas

Free_energy2=sum

if(rank.eq.0)print*, 'FREE_energy :', Free_energy, Free_energy2

if ((abs(Free_energy-Free_energy2)>1.).and.(rank.eq.0))then
stop
endif
!fesum=fesum+Free_energy
!fesum2=fesum2+Free_energy2
!print*, 'fesum',  fesum,  fesum2
return


end subroutine


