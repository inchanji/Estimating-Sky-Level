txt='.txt'
fits='.fits'
NCPU=32
filename='rand_dist_dls_depth_m27to29_2000by2000pix'
ext='_set'
randum='/home/inchani/galsim/randnum.txt'  # make your own a list of random numbers
try=1
while [ $try -lt 21 ];do
num=$(sed -n $try"p" $randum)
echo Try Num $try
echo random number : $num 
cp makegalaxies_in_one0Bndry.yaml makegalaxies_in_one.yaml 
sed -i -e "31s/RADOMNUMBER/$num/g" makegalaxies_in_one.yaml
sed -i -e "32s/NCPU/$NCPU/g" makegalaxies_in_one.yaml 
sed -i -e "42s/FILENAME/$filename$txt/g" makegalaxies_in_one.yaml 
sed -i -e "50s/FILENAME/$filename$ext$try$fits/g" makegalaxies_in_one.yaml 
galsim makegalaxies_in_one.yaml 
try=$(($try+1))
done
