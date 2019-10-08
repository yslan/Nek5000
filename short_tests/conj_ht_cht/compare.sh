
str=$1 # user profide tag
usr=conj_ht
rea=conj_ht
np=2

logT="logfile_"$str"_IFCHAR_T"
logF="logfile_"$str"_IFCHAR_F"

vtkT="vtk_"$logT
vtkF="vtk_"$logF

# Compile
if [ -f "./makenek" ]; then
  ./makenek $usr | tee tttt
  stest=`grep "Compilation successful!" tttt`
else
  make -j  | tee tttt
  stest=`grep "Compilation successful!" tttt`
fi
if [ -z "$stest" ]; then
  echo "Compile fail"
  exit 1
fi
rm tttt

# neklmpi
if [ ! -f "./neklmpi" ]; then
  echo "copy neklmpi from previus nek_v19-rc1/bin/neklmpi !"
  exit 1
fi

# backup .par
cp conj_ht.par conj_ht_bak.par

# IFCHAR=T
echo "Run 1/2: IFCHAR=T"
cp conj_ht_IFCHAR_T.par conj_ht.par
./neklmpi $rea $np
cp logfile $logT
  mkdir -p $vtkT
  visnek $rea
  mv *0.f0* $vtkT
  mv $rea.nek5000 $vtkT

# IFCHAR=T
echo "Run 2/2: IFCHAR=F"
cp conj_ht_IFCHAR_F.par conj_ht.par
./neklmpi $rea $np
cp logfile $logF
  mkdir -p "vtk_"$logF
  visnek $rea
  mv *0.f0* $vtkF
  mv $rea.nek5000 $vtkF

# recover .par
cp conj_ht_bak.par conj_ht.par

# show comparison

echo $logT
tac $logT|grep -m10 Step |tac
tac $logT|grep -m10 tmax$ |tac

echo $logF
tac $logF|grep -m10 Step |tac
tac $logF|grep -m10 tmax$ |tac

