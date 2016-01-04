module colloids
integer nc1, nc2, nc(2), maxnc
real*8 rc1, rc2, rc(2)
real*8 vc(2)
integer dc(2), maxdc
real*8, allocatable :: sph(:,:)
real*8, allocatable :: sphs(:,:)

endmodule colloids

module mkinsol
double precision, allocatable :: pp(:)
endmodule

module multicapa
integer cadenastype
integer, parameter :: ncha_max = 700
real*8, allocatable :: weight(:,:)
real*8 sumweight(2)
real*8 lseg
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

real*8 norma
INTEGER adsmax
integer cuantas(2)          ! number of polymer configuration or  bound sequences
INTEGER nads ! layers already assembled
INTEGER,allocatable :: Tcapas(:) ! layers already assembled

integer ntot ! lattice sites
real*8, allocatable :: avpol(:,:) ! volume fraction polymers already adsorbed
real*8, allocatable :: rhopol(:,:) ! volume fraction polymers already adsorbed
real*8, allocatable :: avpol2(:) ! volume fraction polymer in solution
real*8, allocatable :: rhopol2(:) ! volume fraction polymer in solution
real*8, allocatable :: avpolall(:)

INTEGER cuantas1, cuantas2
INTEGER maxcuantas

integer*4, allocatable :: in1n(:,:,:)
integer, allocatable ::  maxlayer(:,:)

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
integer maxpol

endmodule

module partfunc
real*8 q
endmodule

module layer
real*8, parameter :: delta = 0.5
endmodule

module volume
real*8 vpol, vsol
endmodule


module bulk
REAL*8 expmupol
real*8 xsolbulk, phibulkpol           ! volume fraction of solvent in bulk
endmodule

module seed1
integer seed              ! seed for random number generator
endmodule


module longs
integer long(2)            ! length of polymer

INTEGER long1,long2,maxlong
endmodule

module pis
real*8 pi
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

