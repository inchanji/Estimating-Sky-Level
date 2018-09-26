#!/bin/sh
imgname='SET6_variousPSF_DLSlikeNew_sky3240_stack'
_conv='.conv'
_fits='.fits'

galsimdir='/home/inchani/galsim/sex/'                    
imgdir='/home/inchani/galsim/output_yaml/'              

#THRESH='0.6000'
#THRESH='0.6000'
THRESH='0.5000'
FILTER='gauss_3.0_5x5'
#MINAREA=' 22'
MINAREA=' 6 '
PIXSCALE='0.257'
STARFWHM='0.9'

galsimfile0='myconfig0.sex'
galsimfile='myconfig.sex'
cp $galsimfile0 $galsimfile

OUTPUT=$imgname'_'$THRESH'_'$FILTER$_ext'.sexcat'
sed -i -e "7s/REPLACEIT/$OUTPUT/g" $galsimdir$galsimfile
sed -i -e "18s/REPLACEIT/$MINAREA/g" $galsimdir$galsimfile
sed -i -e "21s/REPLACEIT/$THRESH/g" $galsimdir$galsimfile
sed -i -e "22s/REPLACEIT/$THRESH/g" $galsimdir$galsimfile
sed -i -e "25s/REPLACEIT/$FILTER$_conv/g" $galsimdir$galsimfile
sed -i -e "67s/REPLACEIT/$PIXSCALE/g" $galsimdir$galsimfile
sed -i -e "71s/REPLACEIT/$STARFWHM/g" $galsimdir$galsimfile
sex $imgdir$imgname$_fits -c $galsimfile
mv -f  $OUTPUT $imgdir
mv -f rms.fits $imgdir$imgname"_rms.fits"
mv -f back.fits $imgdir$imgname"_back.fits"
mv -f seg.fits $imgdir$imgname"_seg.fits"
