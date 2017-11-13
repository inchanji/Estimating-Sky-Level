# Subroutine for SkyEstimation 

You can find all data discussed here in example.tar.gz file

To run our example in your computer, you need to specify the locations of:

  1. directory (All.Path)
  2. input mag file (All.inputmagsfile)
  3. mag measurements by SExtracotr (All.sexdatafile)
  4. image file (All.imgfilename)
  5. segmentation map (All.segfilename)
  
in Global 'All' variable.

To estimate sky for the (detected) objects in the entire field, 
you specify the locations of image file (All.imgfilename) and segmentation map(All.segfilename)


When incorporating our subroutines to yours, the only thing you do is to import all subroutines and global variables, then
run estimate_sky(xc,yc,Re) where (xc,yc) is the center of position and Re is the half-light-radius estimated by SExtractor.
Sometines you need to redefine the maximum position of an image (All.MaxPos) or zero-point mag (All.m0).
