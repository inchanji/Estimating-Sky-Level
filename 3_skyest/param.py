import numpy as np
from astropy.io import fits

dirpath 		= "../data/"

imgfilename		= 'exImg.fits'
#segfilename		= "exImg_seg.fits"
newsegfilename  = "exImg_seg2.fits"
sexcatfilename  = "exImg.sexcat"

NpixMin 		= 4000  # min number of pixels of image cutout for sky estimation
LpixMin 		= 80 	# min width/height of image cutout for sky estimation 



nsexcat_buf 	= 28; iXc = 1; iYc = 2; iRe = 34;
ReFac = 15; npix_iter = 5;

path_img  		= dirpath + imgfilename 
#path_seg  		= dirpath + segfilename 
path_newseg 	= dirpath + newsegfilename
path_sexcat 	= dirpath + sexcatfilename

hdu 		= fits.open(path_img)
img 		= hdu[0].data
Nrow, Ncol  = np.shape(img)
#hdu 		= fits.open(path_seg)
#seg 		= hdu[0].data
hdu 		= fits.open(path_newseg)
newseg 		= hdu[0].data

sexcatdata  = np.array([line.rstrip('\n') for line in open(path_sexcat)])[nsexcat_buf:]

