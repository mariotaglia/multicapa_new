salt="0.01" # 1 0.1#
npKEo="5"
phb="5"
dimz="120"
infile="0"
delta="0.5"
pkaA="5"
pkaB="9"
long="40"
curvature="0"
radius="6"
adsmax="16"
bulkpol="1.d-3"
sigma="0.1"
pkion="-0.2"
nsal="4"
lseg="0.35"

for j in $long ; do
echo $j

mkdir testcsal_longjul3_$j
cd testcsal_longjul3_$j

for ph in $phb ; do
echo $ph

mkdir ph_$ph
cd ph_$ph

for sal in $salt ; do
echo $sal

mkdir csalt_$sal
cd csalt_$sal

for pkE in $npKEo ; do
echo $pkE

mkdir pkEo_$pkE
cd pkEo_$pkE

for pAion in $pkion ; do
echo $pAion

mkdir pkaANa_$pAion
cd pkaANa_$pAion

for pA in $pkaA ; do
echo $pA

mkdir pkaA_$pA
cd pkaA_$pA

for pB in $pkaB ; do
echo $pB

mkdir pkaB_$pB
cd pkaB_$pB

for phibulkpol in $bulkpol ; do
echo $phibulkpol

mkdir phipol_$phibulkpol
cd phipol_$phibulkpol


for sig in $sigma ; do
echo $sig

mkdir sigma_$sig
cd sigma_$sig

for z in $dimz ; do
echo $z
 
mkdir dimz_$z
cd dimz_$z
 
for l in $lseg ; do
echo $l

mkdir lseg_$l
cd lseg_$l

echo "
# curvature #
$curvature
#radius
$radius
#adsmax
$adsmax
# dimz #
$z
# chaintype
2
# cuantas1 cuantas 2 #
10000 10000
#long1 long2#
$j $j
#minn corte#
1.0
#phibulkpol#
$phibulkpol
#pKaA
$pA
#pKaB
$pB
#pKaNA
$pAion
#pKaBCl
$pAion
#nsal
$nsal #csal 
1. 
0.5 #0.01 #0.05  0.01 0.001 0.05 $sal0.1 
0.1
$sal
#pH
$ph
#scale interations intermolecule#
1
0.0
#preads#
0
#sigma#
$sig
#error#
1.0d-6
#infile#
$infile
#Kbind# 5 0.001  1. 10. 100.0
125.   #1000.0 1000.0
#eps#
0.0
#Xulimit#
2
#lseg#
$l
#epsilonLJ YukawaA YukawaL#
0.0 0.0 1.0" > fort.8
~/bin/multicapa_ioncsalt

cd ..
done

cd ..
done

cd ..
done

cd ..
done

cd ..
done

cd ..
done

cd ..
done

cd ..
done

cd ..
done

cd ..
done

cd ..
done
