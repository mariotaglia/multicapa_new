
module mkinsol
double precision, allocatable :: pp(:)
endmodule

module multicapa

integer, parameter :: base = 25

integer cadenastype
integer, parameter :: ncha_max = 700
real*8, allocatable :: weight(:,:,:)
real*8 sumweight(2)
real*8 lseg
real*8  vsalt
integer kaitype
real*8, allocatable :: Xu(:,:,:)
real*8 sumXu11, sumXu12
real*8 AA, BA, CA
integer Xulimit
real*8 minn
REAL*8 sts(100), kbinds(100)
INTEGER nst, nkbind
real*8 error              ! error imposed accuaracy
real*8 infile             ! inputfile control variable for reading input files  value 0,1
CHARACTER nada
INTEGER preads ! number of pre-adsorbed layers
real*8  fesum
real*8 fesum2
real*8 norma
INTEGER adsmax
integer AT !G
integer cuantas(2)          ! number of polymer configuration or  bound sequences
INTEGER nads ! layers already assembled
INTEGER,allocatable :: Tcapas(:) ! layers already assembled
real*8 splp
real*8, allocatable :: sumprolnpro(:)

integer ntot ! lattice sites
real*8, allocatable :: avpol(:,:) ! volume fraction polymers already adsorbed
real*8, allocatable :: avpolneg(:) !G 
real*8, allocatable :: avpolpos(:)  !G
real*8, allocatable :: rhopol2(:)
real*8, allocatable :: avpolnegcero(:)
real*8, allocatable :: avpolposcero(:)
real*8, allocatable :: avpol2(:) ! volume fraction polymer in solution
real*8, allocatable :: avpolall(:)
real*8, allocatable :: xsol(:) !!!G
real*8, allocatable :: pro(:)  !!!G
real*16, allocatable :: psi(:)
real*16, allocatable :: xpos(:)
real*16, allocatable :: xneg(:)
real*16, allocatable :: xHplus(:)
real*16, allocatable :: xOHmin(:)
real*8, allocatable :: qtot(:)

INTEGER cuantas1, cuantas2
INTEGER maxcuantas

integer curvature, radio

integer*1, allocatable :: in1n(:,:,:,:)
integer*4, allocatable :: in1tmp(:)
integer, allocatable ::  maxpos(:,:,:)
integer, allocatable ::  minpos(:,:,:)

real*8 sigma

real*8, allocatable :: eps(:)
real*8 eps1
integer newcuantas(2)       ! number of polymer configuration of allowed conformation

integer iter              ! counts number of iterations

REAL*8, allocatable ::  xtotal(:)
real*8 st
integer pegado

REAL*8 Kbind, Kbind0
REAL*8,allocatable :: fbound(:,:)
real*8, allocatable ::  fNcharge(:,:)
integer maxpol

endmodule

module partfunc
real*8 q
endmodule

module layer
real*8, parameter :: delta = 0.5
!real*8, parameter :: deltaR= 0.5!!TEstear
!real*8, parameter :: deltaZ= 1.0!!TEstear
endmodule

module volume
real*8 vpol, vsol
endmodule


module bulk
REAL*8 expmupol
real*8 xsolbulk, phibulkpol           ! volume fraction of solvent in bulk
real*16 expmupos, expmuneg, expmuHplus, expmuOHmin  ! exp(-beta*mu)*(bulk volume fraction), where mu is the chemical potential
real*8 csalt
real*16 pHbulk
real*16 xHplusbulk, xOHminbulk ! bulk volume fraction of H+ and OH-
real*16 xposbulk, xnegbulk  
endmodule

module seed1
integer seed              ! seed for random number generator
endmodule


module longs
integer long(2)            ! length of polymer

INTEGER long1,long2,maxlong
endmodule

module const
real*8, parameter :: Na = 6.02d23 ! Avogadro's number
real*8, parameter :: lb = 0.714   ! bjerrum lenght in water in nm
real*8 zpos, zneg, zpol(2) 
real*8 constq
real*8 wperm
real*8, parameter :: pi=dacos(-1.0d0)          ! pi = arccos(-1) 
real*16 K0A, K0B ,K0ANa,K0BCl, K0Eo!K0
real*16 pKaA
real*16 pkaB
real*16 Kw
real*16 pKw

endmodule

module matrices
real*8 tt(3,3),tp(3,3),tm(3,3)
endmodule

module senos
real*8 sitheta,cotheta,siphip,cophip
endmodule

module MPI
include 'mpif.h' ! librerias MPI
integer rank, size, ierr
integer flagsolver
endmodule

module posmk
real*8, allocatable :: current(:,:)
integer*2, allocatable :: nextbead(:)
integer*2, allocatable :: firstcell(:,:,:)
integer, parameter :: mcube = 100
integer, parameter :: calq = 0
real*8, parameter :: qprob0 = 0.6933
integer, parameter :: nearbonds = 5
endmodule posmk

