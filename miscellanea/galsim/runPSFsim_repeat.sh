txt='.txt'
fits='.fits'
NX=100
NY=100
NCPU=32
filename='rand_dlst_dls_depth_m19to27_8000by8000pix'
ext='_MoffatPSF_set'
randum='/home/inchani/galsim/randnum.txt'

try=1
while [ $try -lt 21 ];do
num=$(sed -n $try"p" $randum)
echo Try Num $try
echo random number : $num 
cp makePSFs0.yaml makePSFs.yaml
sed -i -e "20s/NX/$NX/g" makePSFs.yaml    
sed -i -e "21s/NY/$NY/g" makePSFs.yaml    
sed -i -e "27s/REPLACEIT/{ type : Catalog , col : 9 }/g" makePSFs.yaml    
sed -i -e "29s/NCPU/$NCPU/g" makePSFs.yaml    
sed -i -e "34s/RADOMNUMBER/$num/g" makePSFs.yaml    
sed -i -e "40s/FILENAME/$filename$txt/g" makePSFs.yaml    
sed -i -e "48s/FILENAME/$filename$ext$try$fits/g" makePSFs.yaml    
galsim makePSFs.yaml    
try=$(($try+1))
done
