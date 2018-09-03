import numpy as np
from astropy.io import fits

dirpath 		= "../data/"

#imgfilename		= 'exImg.fits'
segfilename		= "exImg_seg.fits"
newsegfilename  = "exImg_seg2.fits"
sexcatfilename  = "exImg.sexcat"

m0 				= 32.73 # zero point for DLS
mu_lim 			= 0.8 # for ~30 mag /arcsec^2 

nsexcat_buf 	= 28

#path_img  		= dirpath + imgfilename 
path_seg  		= dirpath + segfilename 
path_newseg 	= dirpath + newsegfilename
path_sexcat 	= dirpath + sexcatfilename



#hdu = fits.open(path_img)
#img = hdu[0].data
hdu = fits.open(path_seg)
seg = hdu[0].data

i_seg_max 	= np.max(seg) 	# Max sexcat ID of indentified object
#img_ravel 	= img.ravel()
#seg_ravel 	= seg.ravel()
Nrow, Ncol  = seg.shape

sexcatdata  = np.array([line.rstrip('\n') for line in open(path_sexcat)])[nsexcat_buf:]

