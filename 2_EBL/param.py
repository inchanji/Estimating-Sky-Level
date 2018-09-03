import numpy as np
from astropy.io import fits

dirpath 		= "../data/entirefield/" 

imgfilename		= "entirefield.fits"
segfilename		= "entirefield_seg.fits"
newsegfilename  = "entirefield_seg2.fits"
sexcatfilename  = "entirefield.sexcat"

m0 				= 	32.73  # zero point for DLS
mbrightest 		= 	19 	   # the brightest magnitude of number count [number / mag]
mcut			= 	25.5   # bright end of number count for EBL calculation
mfaintest		= 	29.    # faint end of number count for EBL calculation
mslope 			= 	0.4	   # slope of power-law number count [number / mag]
dmag			= 	0.002  # flux resolution ( ~ 500 ADU) at bright end for EBL flux integration
normconst 		=   6.4347 # normalization constant for galaxy number count 


nsexcat_buf 	= 28; iFlux = 7; iMag = 9;
path_img  		= dirpath + imgfilename 
path_seg  		= dirpath + segfilename 
path_newseg 	= dirpath + newsegfilename
path_sexcat 	= dirpath + sexcatfilename

hdu 		= fits.open(path_img)
img 		= hdu[0].data.ravel()
Ntotpix		= len(img)
hdu 		= fits.open(path_seg)
seg 		= hdu[0].data.ravel()
hdu 		= fits.open(path_newseg)
newseg 		= hdu[0].data.ravel()

sexcatdata  = np.array([line.rstrip('\n') for line in open(path_sexcat)])[nsexcat_buf:]

SExfluxlist = []; SExmaglist = []
for i in range(len(sexcatdata)):
	SExfluxlist.append(float(sexcatdata[i].split()[iFlux]))
	SExmaglist.append(float(sexcatdata[i].split()[iMag]))

SExfluxlist = np.array(SExfluxlist); 
SExmaglist = np.array(SExmaglist)