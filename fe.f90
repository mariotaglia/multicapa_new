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
implicit none

integer cc
integer n
real*8 Free_energy,Free_energy2
real*8 F_mix_s,F_mix_avpolA,F_mix_avpolb,F_conf,F_EQ,Fpro
integer i,ii,j,jj
real*8 sumas,sumrho,sumrhopol,sumpi,sum
Free_energy = 0.0

F_mix_s= 0.0
n=ntot


do i =1,n
!if (xsol(i)>0.0)then
 F_mix_s=F_mix_s+xsol(i)*(log(xsol(i))-1.0)  
 F_mix_s=F_mix_s-xsolbulk*(log(xsolbulk)-1.0)
!endif
enddo

F_mix_s=F_mix_s *delta /vsol
print*,'fmixs',F_mix_s

Free_energy= Free_energy +F_mix_s


F_mix_avpolA=0.0
!!Consulta   avpol2(i) == rho_a*vp ? Lo estoy tomando asi por ahora
do i=1,n
!if (avpol2(i)>0.) then
 F_mix_avpolA=F_mix_avpolA+(avpol2(i))*(log((avpol2(i))/vpol)-1.0 ) 
 F_mix_avpolA=F_mix_avpolA-(phibulkpol)*(log((phibulkpol)/vpol)-1.0)
!endif
enddo


F_mix_avpolA=F_mix_avpolA *delta/(vsol*vpol) 

Free_energy= Free_energy +F_mix_avpolA

print*,'fmixA',F_mix_avpolA

F_Conf=0
!!!!!!
do i=1,n
do j=1,newcuantas(AT)
!if (avpol2(i)>0.)then
F_Conf=F_Conf+(avpol2(i))*(pro(j)/q)*log(pro(j)/q) 
!endif

enddo
enddo

F_conf=F_conf*delta/(vpol*vsol)  !! TEstear si es vpol o vpol *M
print*,'F_conf',F_conf
!!!!!!

Free_energy=Free_energy+F_conf

F_EQ=0.
!!!!   avpolpos==<na> y avpolneg == <nb> 
do i=1, n
F_EQ=F_EQ+avpolpos(i)*fbound(1,i)*(-log(Kbind0)) !chequear signo
F_EQ=F_EQ+avpolpos(i)*(1.-fbound(1,i))*(log(1-fbound(1,i)))
if (fbound(1,i)>0.)then
F_EQ=F_EQ+avpolpos(i)*(fbound(1,i))*(log(fbound(1,i)))
if (avpolpos(i)>0)then
F_eq=F_eq-avpolpos(i)*fbound(1,i)*(log(avpolpos(i)*fbound(1,i))-1.0 ) !vab =1
endif
endif
enddo

print*,'Feq',F_eq


do i = 1,n
if (fbound(2,i)>0.)then
 F_eq=F_eq + avpolneg(i)*(fbound(2,i)*log(fbound(2,i)))	 !! 
 F_eq=F_eq + avpolneg(i)*(1.-(fbound(2,i))*log(1.-fbound(2,i))) !
endif
enddo


F_EQ=F_EQ*delta/(vpol*vsol)

print*,'Feq',F_eq

Free_energy=Free_energy+F_EQ

Free_energy2=0.0
sumpi=0.0
sumrho=0.0
sumrhopol=0.0
sumas=0.0
do i=1, n
! if (xsol(i)>0.)then
 sumpi= sumpi + log(xsol(i)) 
 sumpi= sumpi - log(xsolbulk)
 sumpi= sumpi*(1.0-avpolneg(i)*vpol*vsol-avpolposcero(i)*vpol*vsol)! 
! sumpi= sumpi*(1.0-avpolneg(i)-avpolposcero(i))! 


 sumrho= sumrho+(-xsol(i)) 
 sumrho= sumrho-(-xsolbulk)
!endif
! if(avpol2(i)>0.)then
 sumrhopol=sumrhopol-(avpol2(i))
 sumrhopol=sumrhopol-(-phibulkpol)
! endif

 sumas=sumas+avpolpos(i)*fbound(1,i)
 sumas=sumas+avpolneg(i)*log(1.-fbound(2,i))
 sumas=sumas+avpolposcero(i)*log(1.-fbound(1,i))

enddo

sumpi=sumpi*delta/vsol
sumrho=sumrho*delta/vsol
sumrhopol=sumrhopol*delta/(vpol*vsol)
sumas=sumas*delta/(vpol*vsol)



print*,'sumpi',sumpi
print*,'sumrho',sumrho
print*,'sumrhopol',sumrhopol
print*,'sumas',sumas

sum=sumpi+sumrho+sumrhopol+sumas



Free_energy2=sum

print*, 'FREE_energy :', Free_energy, Free_energy2
!fesum=fesum+Free_energy
!fesum2=fesum2+Free_energy2
!print*, 'fesum',  fesum,  fesum2
return

end subroutine


