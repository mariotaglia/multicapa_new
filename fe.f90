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
use const
implicit none

integer cc
integer n
real*8 Free_energy,Free_energy2, Fmupol
real*8 F_mix_s,F_mix_avpolA,F_mix_avpolb,F_conf,F_EQ,Fpro,F_mix_pos,F_mix_neg
real*8 F_mix_OHmin, F_mix_Hplus,F_electro
integer i,ii,j,jj
real*8 sumas,sumrho,sumrhopol,sumpi,sum,sumel
real*8 Fact_rhobulk
integer iz,ir,dimz
double precision, external :: jacobian
integer sumnewcuantas
integer err
integer maxR, minR

if(curvature.ge.0) then
 maxR = ntot
 minR = radio
else 
 maxR = ntot
 minR = 1
endif

Free_energy=0.0
Free_energy2=0.0

if(rank.eq.0) then
print*, 'Starting free energy calculation...'
! open files
open(unit=301, file='F_tot.dat')
!open(unit=302, file='F_mixs.dat')
!open(unit=303, file='F_mixpos.dat')
!open(unit=304, file='F_mixneg.dat')
!open(unit=305, file='F_mixH.dat')
!open(unit=306, file='F_mixOH.dat')
!open(unit=307, file='F_conf.dat')
!open(unit=308, file='F_eq.dat')

!do is=1,Npoorsv
!do js=1,Npoorsv
!write(F_vdWfilename(is,js),'(A6,BZ,I3.3,A1,I3.3,A4)')'F_vdW.',is,'.',js,'.dat'
!open(unit=10000*is+js, file=F_vdWfilename(is,js) )
!enddo
!enddo

open(unit=312, file='F_tot2.dat')
!open(unit=313, file='F_mixp.dat')
!open(unit=314, file='F_Uchain.dat')
endif




call MPI_REDUCE(newcuantas(AT), sumnewcuantas, 1, MPI_INTEGER, MPI_SUM,0, MPI_COMM_WORLD, err)

dimz=1

F_mix_s= 0.0
n=ntot

!!!!!!!!! F mix s !!!!!!!!!!!!!

Fact_rhobulk=phibulkpol/(float(long(AT))*vpol*vsol)
do iR =minR, maxR
 F_mix_s=F_mix_s+xsol(iR)*(log(xsol(iR))-1.0)*jacobian(iR)
 F_mix_s=F_mix_s-xsolbulk*(log(xsolbulk)-1.0)*jacobian(ir)
enddo

F_mix_s=F_mix_s *delta /vsol 
if(rank.eq.0)print*,'fmixs',F_mix_s

Free_energy= Free_energy +F_mix_s



! 2. cations entropy

F_Mix_pos = 0.0

do iR = minR, maxR
  F_Mix_pos = F_Mix_pos + xpos(ir)*(log(xpos(ir)/vsalt)-1.0 - log(expmupos) + dlog(vsalt))*jacobian(ir)
  F_Mix_pos = F_Mix_pos - xposbulk*(log(xposbulk/vsalt)-1.0 - log(expmupos) + dlog(vsalt))*jacobian(ir)
enddo

F_Mix_pos = F_Mix_pos * delta/vsol/vsalt
Free_Energy = Free_Energy + F_Mix_pos


!
! 3. anions entropy

F_Mix_neg = 0.0

do iR = minR, maxR
  F_Mix_neg = F_Mix_neg + xneg(ir)*(log(xneg(ir)/vsalt)-1.0 - log(expmuneg) + dlog(vsalt))*jacobian(ir)
  F_Mix_neg = F_Mix_neg - xnegbulk*(log(xnegbulk/vsalt)-1.0 - log(expmuneg) + dlog(vsalt))*jacobian(ir)
enddo

F_Mix_neg = F_Mix_neg * delta/vsol/vsalt
Free_Energy = Free_Energy + F_Mix_neg
!


! 4. proton entropy

F_Mix_Hplus = 0.0

do iR = minR,maxR
  F_Mix_Hplus = F_Mix_Hplus + xHplus(ir)*(log(xHplus(ir))-1.0 - log(expmuHplus))*jacobian(ir)
  F_Mix_Hplus = F_Mix_Hplus - xHplusbulk*(log(xHplusbulk)-1.0 - log(expmuHplus))*jacobian(ir)
enddo

F_Mix_Hplus = F_Mix_Hplus * delta/vsol
Free_Energy = Free_Energy + F_Mix_Hplus
!


! 5. hydroxyl ions entropy

F_Mix_OHmin = 0.0

do iR = minR, maxR
  F_Mix_OHmin = F_Mix_OHmin + xOHmin(ir)*(log(xOHmin(ir))-1.0 - log(expmuOHmin))*jacobian(ir)
  F_Mix_OHmin = F_Mix_OHmin - xOHminbulk*(log(xOHminbulk)-1.0 - log(expmuOHmin))*jacobian(ir)
enddo

F_Mix_OHmin = F_Mix_OHmin * delta/vsol
Free_Energy = Free_Energy + F_Mix_OHmin


!!!!!!! F mix A !!!!!!!!!!!!!!! polymer in solution

F_mix_avpolA=0.0
do iR=minR,maxR
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
do iR=minR,maxR
 Fmupol=Fmupol-rhopol2(ir)*log(expmupol)*jacobian(ir)
 Fmupol=Fmupol+Fact_rhobulk*log(expmupol)*jacobian(ir)
enddo
Fmupol=Fmupol*delta !-delta

Free_energy= Free_energy +Fmupol

if(rank.eq.0)print*,'fmupol',Fmupol

!!!!!! Fconf !!!!!!!!!!!!!!!!!!!!!

F_Conf=0
!!!!!!
do iR=minR,maxR
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


!!!!!!!!! F eq !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


F_EQ=0.
!!!!   avpolpos==<na> y avpolneg == <nb> 

do iR=minR, maxR ! Negativo

F_EQ=F_EQ+avpolneg(iR)*fbound(1,iR)*(-log(Kbind0)) *jacobian(iR)
F_EQ=F_EQ+avpolneg(iR)*(1.0-fbound(1,iR)-fNcharge(1,ir))*(log(1.0-fbound(1,iR)-fNcharge(1,ir))) *jacobian(iR)
F_Eq=F_eq+avpolneg(iR)*fNcharge(1,iR)*log(fNcharge(1,iR))*jacobian(iR)
F_Eq=F_eq+avpolneg(iR)*fNcharge(1,iR)*log(K0A/xsolbulk)*jacobian(iR)
F_Eq=F_eq-avpolneg(iR)*fNcharge(1,iR)*(log(expmuHplus))*jacobian(iR)

if (fbound(1,iR).gt.0.0)then
  F_eq=F_eq+avpolpos(iR)*(fbound(1,iR))*(log(fbound(1,iR))) *jacobian(iR)
  F_eq=F_eq-avpolpos(iR)*fbound(1,iR)*(log(avpolpos(iR)/vpol/vsol*fbound(1,iR))-1.0)*jacobian(iR) !vab =1 <= no, ojo, aca va vpol
endif
enddo

if(rank.eq.0)print*,'Feq',F_eq*delta/(vpol*vsol) 

do iR = minR,maxR ! Positivo

if (fbound(2,iR).gt.0.0) then
 F_eq=F_eq + avpolpos(iR)*fbound(2,iR)*log(fbound(2,iR)) *jacobian(iR)
endif
 F_eq=F_eq + avpolpos(iR)*(1.-fbound(2,iR)-fNcharge(2,ir))*log(1.-fbound(2,iR)-fNcharge(2,ir)) *jacobian(iR) !<= ojo, habia error de parentesis
 F_eq=F_Eq + avpolpos(iR)*(fncharge(2,iR)*log(fncharge(2,iR)))*jacobian(iR)
 F_eq=F_eq + avpolpos(iR)*fNcharge(2,iR)*log(k0B/xsolbulk)*jacobian(iR)
 F_Eq=F_eq - avpolpos(iR)*fNcharge(2,iR)*log(expmuOHmin)*jacobian(iR)
enddo

if(rank.eq.0)print*,'Feq',F_eq*delta/(vpol*vsol) 

!!!! Substract bulk

if (AT.eq.1) then 
  do iR=minR, maxR ! Negativo
     F_EQ=F_EQ-phibulkpol*(1.0-fNchargebulk(AT))*(log(1.0-fNchargebulk(AT))) *jacobian(iR)
     F_Eq=F_eq-phibulkpol*fNchargebulk(AT)*log(fNchargebulk(AT))*jacobian(iR)
     F_Eq=F_eq-phibulkpol*fNchargebulk(AT)*log(K0A/xsolbulk)*jacobian(iR)
     F_Eq=F_eq+phibulkpol*fNchargebulk(AT)*(log(expmuHplus))*jacobian(iR)
 enddo

else if (AT.eq.2) then

   do iR = minR,maxR ! Positivo
     F_EQ=F_EQ-phibulkpol*(1.0-fNchargebulk(AT))*(log(1.0-fNchargebulk(AT))) *jacobian(iR)
     F_Eq=F_eq-phibulkpol*fNchargebulk(AT)*log(fNchargebulk(AT))*jacobian(iR)
     F_eq=F_eq-phibulkpol*fNchargebulk(AT)*log(k0B/xsolbulk)*jacobian(iR)
     F_Eq=F_eq+phibulkpol*fNchargebulk(AT)*log(expmuOHmin)*jacobian(iR)
  enddo
endif

if(rank.eq.0)print*,'Feq',F_eq*delta/(vpol*vsol) 

F_EQ=F_EQ*delta/(vpol*vsol) 
Free_energy=Free_energy+F_EQ

! 9. Electrostatic !

F_electro = 0.0

do iR  = minR, maxR
  F_electro = F_electro + delta*psi(iR)*qtot(iR)/2.0/vsol *jacobian(iR)
enddo

Free_Energy = Free_Energy + F_electro
if(rank.eq.0)print*,'Felec',F_electro



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Free energy at extrema 

Free_energy2=0.0
sumpi=0.0
sumrho=0.0
sumrhopol=0.0
sumas=0.0
do iR=minR, maxR
! if (xsol(i)>0.)then
 sumpi= sumpi + log(xsol(iR))*(1.0-avpolnegcero(iR)-avpolposcero(iR)) *jacobian(iR)
 sumpi= sumpi - log(xsolbulk)*jacobian(iR)

 sumrho= sumrho+(-xsol(iR) -xHplus(ir) -xOHmin(ir)-(xpos(ir)+xneg(ir))/vsalt)*jacobian(iR) !!
 sumrho= sumrho-(-xsolbulk-xHplusbulk -xOHminbulk-(xposbulk+xnegbulk)/vsalt)*jacobian(iR)!!

 sumel = sumel - qtot(ir)*psi(ir)/2.0 * jacobian(iR)
 sumel = sumel + avpolnegcero(iR)*zpol(1)/vpol*psi(iR)*jacobian(iR)    ! electrostatic part free en
 sumel = sumel + avpolposcero(iR)*zpol(2)/vpol*psi(iR)*jacobian(iR)    ! electrostatic part free energy  

! sumel = sumel + (avpolnegcero(iR)*zpol(1)+avpolposcero(iR)*zpol(2))/vpol*psi(iR)/2.0*jacobian(iR)    ! electrostatic part free energy  


 sumrhopol=sumrhopol-(rhopol2(iR)) *jacobian(iR)
 sumrhopol=sumrhopol-(-Fact_rhobulk)*jacobian(iR)


if (AT.eq.1) then
 sumas=sumas+avpolneg(iR)*fbound(1,iR)*jacobian(iR)

 sumas=sumas+avpolnegcero(iR)*log(1.-fbound(1,iR)-fncharge(1,iR))*jacobian(iR)

! sumas=sumas+avpolnegcero(iR)*fNcharge(1,iR)*log(fNcharge(1,iR))*jacobian(iR)
! sumas=sumas+avpolnegcero(iR)*fNcharge(1,iR)*log(K0A/xsolbulk)*jacobian(iR)
! sumas=sumas-avpolnegcero(iR)*fNcharge(1,iR)*(log(expmuHplus))*jacobian(iR)

 sumas=sumas+avpolposcero(iR)*log(1.-fbound(2,iR)-fncharge(2,iR))*jacobian(iR)

! sumas=sumas+avpolposcero(iR)*(fncharge(2,iR)*log(fncharge(2,iR)))*jacobian(iR)
! sumas=sumas+avpolposcero(iR)*fNcharge(2,iR)*log(k0B/xsolbulk)*jacobian(iR)
! sumas=sumas-avpolposcero(iR)*fNcharge(2,iR)*log(expmuOHmin)*jacobian(iR)

else

 sumas=sumas+avpolpos(iR)*fbound(2,iR)*jacobian(iR)

 sumas=sumas+avpolposcero(iR)*log(1.-fbound(2,iR)-fncharge(2,iR))*jacobian(iR)

! sumas=sumas+avpolposcero(iR)*(fncharge(2,iR)*log(fncharge(2,iR)))*jacobian(iR)
! sumas=sumas+avpolposcero(iR)*fNcharge(2,iR)*log(k0B/xsolbulk)*jacobian(iR)
! sumas=sumas-avpolposcero(iR)*fNcharge(2,iR)*log(expmuOHmin)*jacobian(iR)

 sumas=sumas+avpolnegcero(iR)*log(1.-fbound(1,iR)-fNcharge(1,iR))*jacobian(iR)

! sumas=sumas+avpolnegcero(iR)*fNcharge(1,iR)*log(fNcharge(1,iR))*jacobian(iR)
! sumas=sumas+avpolnegcero(iR)*fNcharge(1,iR)*log(K0A/xsolbulk)*jacobian(iR)
! sumas=sumas-avpolnegcero(iR)*fNcharge(1,iR)*(log(expmuHplus))*jacobian(iR)

endif
enddo

sumpi=sumpi*delta/vsol!
sumrho=sumrho*delta/vsol!
sumrhopol=sumrhopol*delta!
sumel = sumel*delta/vsol  
sumas=sumas*delta/(vpol*vsol)!



if(rank.eq.0)print*,'sumel',sumel
if(rank.eq.0)print*,'sumpi',sumpi
if(rank.eq.0)print*,'sumrho',sumrho
if(rank.eq.0)print*,'sumrhopol',sumrhopol
if(rank.eq.0)print*,'sumas',sumas

sum=sumpi+sumrho+sumrhopol+sumel+sumas

Free_energy2=sum

if(rank.eq.0)print*, 'FREE_energy :', Free_energy, Free_energy2

!if ((abs(Free_energy-Free_energy2)>1.).and.(rank.eq.0))then
!stop
!endif


if(rank.eq.0) then
  write(301,*) nads+1,cc,  Free_energy
!  write(302,*)  F_Mix_s
!  write(303,*)  F_Mix_avpolA
!  write(304,*)  F_mupol
 ! write(305,*) npol, F_Mix_Hplus/npol
!  write(306,*) npol, F_Mix_OHmin/npol
 ! write(307,*) npol, F_Conf/npol
 ! write(308,*) npol, F_Eq/npol
! write(313,*)counter, counter2, F_Eq_P                              
 ! do is=1,Npoorsv
 !   do js=1,Npoorsv
 !      write(10000*is+js,*) npol, F_vdW(is,js)/npol
 !   enddo
 ! enddo
  write(312,*)cc,  Free_energy2
 ! write(313,*) npol, F_Mix_p/npol
 ! write(314,*) npol, F_Uchain/npol
endif

return



end subroutine


