import numpy as np
from astropy.io import fits

dirpath 		= "../data/" 
imgfilename		= "exImg.fits"
segfilename		= "exImg_seg.fits"
rmsfilename  	= "exImg_rms.fits"
sexcatfilename  = "exImg.sexcat"

path_img  		= dirpath + imgfilename 
path_seg  		= dirpath + segfilename 
path_rms 		= dirpath + rmsfilename
path_sexcat 	= dirpath + sexcatfilename


	




class sExcatCol:
	nsexcat_buf = 28; 
	ID = 0; X_IMAGE = 1; Y_IMAGE = 2; XMIN_IMAGE = 3; YMIN_IMAGE = 4; XMAX_IMAGE = 5; YMAX_IMAGE = 6;
	MAG_AUTO = 9; BACKGROUND = 30; FLUX_RADIUS = 34; FLUX_AUTO = 7; ELLIPTICITY = 41; THETA_IMAGE = 40;
	CXX_IMAGE = 35; CYY_IMAGE = 36; CXY_IMAGE = 37; A_IMAGE = 38; 
	#X_galfit = 1; Y_galfit = 2; MAG_galfit = 3; Re_galfit = 4; n_galfit = 6; PA = 9;

class galfit:
	galfitProfType 	= "sersic"
	plateSacle 		= 0.257 	# arcsec / pixel
	initsersicn 	= 2.5 		# init geuss for sersic index
	facImgCutout 	= 1.0		# free param for size of image cutout 
	useGalfitSky	= False

	m0				= 25.35		# zero mag for Galfit run (see userguide)
	

	# CCD parameters
	GAIN 			= 3.00 
	EXPTIME			= 900.  	# in [sec]
	NCOMBINE        = 20 	


class image:
	m0				= 32.73		# zero mag of CCD
	mu_lim 			= 0.8 		# for ~30 mag /arcsec^2 

	hdu 			= fits.open(path_img)
	img 			= hdu[0].data.copy()
	Nrow, Ncol  	= np.shape(img)
	hdu 			= fits.open(path_seg)
	seg 			= hdu[0].data.copy()
	hdu 			= fits.open(path_rms)
	rms 			= hdu[0].data.copy()

class data:
	sexcatdata  	= np.array([line.rstrip('\n') for line in open(path_sexcat)])[sExcatCol.nsexcat_buf:]



