import numpy as np 
import param as pm
from scipy import special
from astropy.io import fits
import scipy.ndimage as ndimage
import os

def concentration(n):
	return 1.9992 * n - 0.3271

def rlim(n, Ie, mu_lim):
	bn = concentration(n)
	return (1. - np.log(mu_lim/Ie)/bn)**n

def sersic(x, p):
	bn = concentration(p[2])
	return p[0] * np.exp(bn) * np.exp( -bn * (x/p[1])**(1./p[2]) )

def Ieff(m, n, hlr) :
	L = 10.**((pm.image.m0-m)/2.5)
	bn = 1.9992*n-0.3271
	Ie = L * bn**(2.*n) / (2.*np.pi * hlr**2. * n * np.exp(bn) * special.gamma(2*n))
	return Ie	


def growMaskMore(MaskMap, fac = 2, val = 1):
	Ny, Nx = MaskMap.shape
	x = np.linspace(-1.*fac,1.*fac,fac*2+1)
	X,Y = np.meshgrid(x, x)
	Z   = X**2. + Y**2.
	Z[Z <= 1.*fac**2.] = np.nan; 
	MaskMap0 = np.zeros((Ny, Nx)) * 1.
	MaskMap0[:] = MaskMap[:]
	ind = np.where(MaskMap[:] == val)   

	MaskMap0[ind] = np.nan
	if np.sum(ind) > 0:
		i_row = ind[0]; i_col = ind[1]
		N = len(i_col)   
		for i in range(N):
			x = i_col[i]; y = i_row[i]
			x0 = x-fac; y0 = y-fac; x1 = x+fac; y1 = y+fac
			x00 = y00 = x11 = y11 = 0
			if x0 < 0: 
				x0 = 0; x00 = fac - x; 
			if y0 < 0: 
				y0 = 0; y00 = fac - y
			if y1 > Ny-1:
				y11 = y1 - Ny +1; y1 = Ny-1; 
			if x1 > Nx-1:
				x11 = x1 - Nx +1; x1 = Nx-1; 
			MaskMap0[y0:y1+1,x0:x1+1] += Z[y00:2*fac+1-y11,x00:2*fac+1-x11]
	MaskMap[np.isnan(MaskMap0)] = val
	del MaskMap0
	return MaskMap

def getSizeImgCutout(data, minsize = 5):

	Re 			= float(data[pm.sExcatCol.FLUX_RADIUS])
	A_IMAGE 	= float(data[pm.sExcatCol.A_IMAGE])
	Mag 		= float(data[pm.sExcatCol.MAG_AUTO])
	axisratio 	= 1. - float(data[pm.sExcatCol.ELLIPTICITY])
	PA 			= float(data[pm.sExcatCol.THETA_IMAGE])
	PA 			*= np.pi / 180.
	Ie  		= Ieff(Mag, 4, Re) # ADU / pixel^2
	R30 		= np.min([rlim(4,Ie,pm.image.mu_lim), 20.]) # in the unit of Re

	if np.isnan(R30):  R30 = 5. 
	if Re * R30 < 5: R30 = 5. / Re

	RlimImgCutout  	= R30 * pm.galfit.facImgCutout	

	height 	= np.max([RlimImgCutout * Re * (np.abs(np.sin(PA)) + axisratio * np.abs(np.cos(PA)) ), minsize])	
	width  	= np.max([RlimImgCutout * Re * (np.abs(np.cos(PA)) + axisratio * np.abs(np.sin(PA)) ), minsize])
	return width, height


def getImgRange(data, width, height):
	xc = float(data[pm.sExcatCol.X_IMAGE])
	yc = float(data[pm.sExcatCol.Y_IMAGE])
	
	col0 = int(xc - width-1);  col1 = int(xc + width-1);
	row0 = int(yc - height-1); row1 = int(yc + height-1);	
	if col0 < 0: col0 = 0
	if row0 < 0: row0 = 0
	if col1 > pm.image.Ncol-1: col1 = pm.image.Ncol-1
	if row1 > pm.image.Nrow-1: row1 = pm.image.Nrow-1

	cen_col = xc - col0
	cen_row = yc - row0

	return [cen_row, cen_col], [row0,row1,col0,col1]


def make_objimg(data, imgrange, sky, useGalfitSky):
	FileOut = pm.dirpath + "galfit/objimg/ObjImg_" + data[pm.sExcatCol.ID] + ".fits"

	if useGalfitSky: 
		img 	= pm.image.img[imgrange[0]:imgrange[1]+1, imgrange[2]:imgrange[3]+1].copy()
	else:  
		img 	= pm.image.img[imgrange[0]:imgrange[1]+1, imgrange[2]:imgrange[3]+1].copy() - sky

	hduout 	= fits.PrimaryHDU(img)

	prihdr = hduout.header

	prihdr.set("CREATED", "username")
	prihdr.set("GAIN", pm.galfit.GAIN, "[e-/ADU]")
	prihdr.set("EXPTIME", pm.galfit.EXPTIME, "[s] exposure time (of exposure 1)")
	prihdr.set("NCOMBINE", pm.galfit.NCOMBINE, "number of exposures")

	hduout.writeto(FileOut, overwrite = "true")
	
	return FileOut	



def make_mask_img(ID, imgrange, useGalfitSky):
	FileOut  	= pm.dirpath + "galfit/segmentation/MaskImg_" + str(ID) + ".fits"
	Img      	= pm.image.img[imgrange[0]:imgrange[1]+1, imgrange[2]:imgrange[3]+1].copy()
	MaskMap0 	= pm.image.seg[imgrange[0]:imgrange[1]+1, imgrange[2]:imgrange[3]+1].copy()
	MaskMap 	= MaskMap0.copy()
	SegImgRavel = MaskMap.ravel()

	IDsSeg = []
	for i in range(1, np.max(SegImgRavel)+1):
		if np.sum(i == SegImgRavel) > 0: IDsSeg.append(i)
	IDsSeg = np.array(IDsSeg)


	if not(useGalfitSky):
		for i in range(len(IDsSeg)):
			data 	= pm.data.sexcatdata[IDsSeg[i]-1].split()
			x   	= int(float(data[pm.sExcatCol.X_IMAGE])); 
			y  		= int(float(data[pm.sExcatCol.Y_IMAGE])); 
			cxx 	= float(data[pm.sExcatCol.CXX_IMAGE]);
			cyy 	= float(data[pm.sExcatCol.CYY_IMAGE]);
			cxy 	= float(data[pm.sExcatCol.CXY_IMAGE]);
			A   	= float(data[pm.sExcatCol.A_IMAGE]);
			hlr 	= float(data[pm.sExcatCol.FLUX_RADIUS]);
			Mag 	= float(data[pm.sExcatCol.MAG_AUTO]);

			Ie  	= Ieff(Mag, 4, hlr) 		# ADU / pixel^2
			fac 	= np.min([rlim(4, Ie, pm.image.mu_lim), 20.]) 

			if np.isnan(fac):  fac = 5. 
			if hlr * fac < 5: fac = 5. / hlr

			lim = hlr * fac / A	

			y1d 		= np.arange(imgrange[1]-imgrange[0]+1) * 1. + imgrange[0] - y
			x1d 		= np.arange(imgrange[3]-imgrange[2]+1) * 1. + imgrange[2] - x

			x2d1, y2d1 	= np.meshgrid(x1d,y1d) 

			if not(IDsSeg[i] == ID):
				#MaskMap 		  	= growMaskMore(MaskMap, 2, IDsSeg[i])
				ind_inc			  	= cxx * (x2d1)**2. + cyy * (y2d1)**2. + cxy * (x2d1) * (y2d1) <  lim**2. 
				MaskMap[ind_inc] 	= IDsSeg[i]			

	ind 		= (MaskMap == 0) | (MaskMap == ID)
	indn 		= np.logical_not(ind)

	MaskMap[ind] = 0; MaskMap[indn]= 1
	hduout 	 = fits.PrimaryHDU(MaskMap)
	hduout.writeto(FileOut, overwrite="true")

	return FileOut, MaskMap0


def make_sigma_img(MaskMap, ID, imgrange, useGalfitSky):
	FileOut  	= pm.dirpath + "galfit/sigma/SigmaImg_" + str(ID) + ".fits"

	Img      	= pm.image.img[imgrange[0]:imgrange[1]+1, imgrange[2]:imgrange[3]+1].copy()
	RmsMap  	= pm.image.rms[imgrange[0]:imgrange[1]+1, imgrange[2]:imgrange[3]+1].copy()
	
	#rms1 		= np.median(RmsMap[MaskMap == ID])
	#rms2 		= np.median(RmsMap)
	#RmsMap[:] 	= rms2

	Sigma = np.sqrt(Img / pm.galfit.GAIN +  RmsMap **2.) 

	#################################################################
	if not(useGalfitSky):
		Sigma 		= ndimage.gaussian_filter(Sigma, sigma=0.5, order=0)
		ind 		= np.logical_not((MaskMap == ID) | (MaskMap == 0))
		#ind 		= np.logical_not(MaskMap == 0)
		Sigma[ind] 	= 1e10
	#################################################################

	hduout = fits.PrimaryHDU(Sigma)
	hduout.writeto(FileOut, overwrite="true")

	del hduout, Img, RmsMap, Sigma
	return FileOut

def writeFeedmeFile(FileNames, tgtpos, imgrange, data, useGalfitSky):
	f = open(FileNames[0],'w')

	f.write("A)\t"+ FileNames[1] + "\t# Input data image (FITS file)\n")
	f.write("B)\t"+ FileNames[2] + "\t# Output data image block\n")
	f.write("C)\t"+ FileNames[3] + "\t# Sigma image name (made from data if blank or 'none') \n")
	#f.write("C)\t none\t# Sigma image name (made from data if blank or "none") \n")
	f.write("D)\t"+ FileNames[4] + "\t# Input PSF image and (optional) diffusion kernel\n")
	f.write("E)\t1\t# PSF fine sampling factor relative to data \n")
	f.write("F)\t"+ FileNames[5] + "\t# Bad pixel mask (FITS image or ASCII coord list)\n")
	f.write("G)\t"+ FileNames[6] +"\t# File with parameter constraints (ASCII file) \n")

	f.write("H)\t"+ str(1)+"\t"+str(imgrange[3]-imgrange[2]+1)+"\t"+str(1)+"\t"+str(imgrange[1]-imgrange[0]+1)+"\t# Image region to fit (xmin xmax ymin ymax)\n")
	f.write("I)\t"+ str(int((imgrange[3]-imgrange[2]+1)*4.))+"\t"+str(int((imgrange[1]-imgrange[0]+1)*4.))+"\t# Size of the convolution box (x y)\n")

	f.write("J)\t"+ str(pm.galfit.m0)+ "\t# Magnitude photometric zeropoint \n")  # If I adjust it to 2E+04... it will make Mags comparable to SExtractor
	f.write("K)\t"+ str(pm.galfit.plateSacle) + "\t" + str(pm.galfit.plateSacle) + "\t# Plate scale (dx dy)    [arcsec per pixel]\n")
	f.write("O)\tregular\t# Display type (regular, curses, both)\n")
	f.write("P)\t0\t# Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n")
	f.write("\n")

	f.write(" 0) "+ str(pm.galfit.galfitProfType) 				+ " \t# Object type\n")
	f.write(" 1) "+ str(tgtpos[1])+" "+str(tgtpos[0])			+ " 1 1 # position x, y [pixel]\n")
	f.write(" 3) "+ data[pm.sExcatCol.MAG_AUTO]					+"  1       # total magnitude\n")
	f.write(" 4) "+ data[pm.sExcatCol.FLUX_RADIUS]				+"  1       # R_e [Pixels]\n")
	f.write(" 5) "+ str(pm.galfit.initsersicn)      			+"  1       # Sersic exponent (deVauc=4, expdisk=1)\n")
	f.write(" 9) "+ str(1-float(data[pm.sExcatCol.ELLIPTICITY]))+"  1       # axis ratio (b/a)\n")
	f.write("10) "+ data[pm.sExcatCol.THETA_IMAGE]	   			+"  1       # position angle (PA)  [Degrees: Up=0, Left=90]\n")
	f.write(" Z) "+ "0" 	    								+"          #  Skip this model in output image?  (yes=1, no=0)\n")

	
	f.write("\n")

	if useGalfitSky:
		f.write(" 0) sky             		# Object type\n")
		################################################################################################
		f.write(" 1) "+data[pm.sExcatCol.BACKGROUND]   +"      1 #  SExtractor sky at center of fitting region [ADUs]\n") ##
		###############################################################################################s#
		f.write(" 2) 0.0000      1          #  dsky/dx (sky gradient in x)\n")
		f.write(" 3) 0.0000      1          #  dsky/dy (sky gradient in y)\n")
		f.write(" Z) 0                      #  output option (0 = resid., 1 = Don\"t subtract)\n")

	f.close()


def writeConstraint(FileOut): 
	f = open(FileOut,'w')

	f.write('1 \t re \t -2 2\n') 		# Re (or HLR): init val +/- 2 mag
	f.write('1 \t mag \t  -2 2\n') 		# mag
	f.write('1 \t x \t -1 1\n') 		# x pos
	f.write('1 \t y \t -1 1\n') 		# y pos		
	f.write('1 \t n \t 0.2 to 8\n') 	# sersic n : n = 0.2 ~ 8
	f.write('1 \t 9 \t  0.1 to 1\n') 	# axis ratio

	f.close()


def sbpFitGal(IDgal, sky, useGalfitSky):
	OutputFileName  		= pm.dirpath + "galfit/GalfitResult_ObjID_" + str(IDgal) + ".fits"
	ConstraintsFileName 	= pm.dirpath + "galfit/constraint/galfit_ID_" + str(IDgal) + ".constraint"
	FeedmeFileName			= pm.dirpath + "galfit/feedme/galfit_ID_"+str(IDgal)+".feedme"	
	PSFFileName 			= pm.dirpath + "galfit/PSF/PSFid_"+str(IDgal)+".fits"	

	SExCatData 				= pm.data.sexcatdata[IDgal-1].split()

	wImgCutout, hImgCutout 	= getSizeImgCutout(SExCatData)
	TgtPos, ImgRange		= getImgRange(SExCatData, wImgCutout, hImgCutout)
	print(TgtPos, ImgRange)
	ImgFileName				= make_objimg(SExCatData, ImgRange, sky, useGalfitSky)
	MaskFileName, MaskMap	= make_mask_img(IDgal, ImgRange, useGalfitSky)
	SigmaFileName			= make_sigma_img(MaskMap, IDgal, ImgRange, useGalfitSky)
	
	FileNames = [FeedmeFileName, ImgFileName, OutputFileName, SigmaFileName, PSFFileName, MaskFileName, ConstraintsFileName]
	writeFeedmeFile(FileNames,TgtPos, ImgRange, SExCatData, useGalfitSky)
	writeConstraint(ConstraintsFileName)

	command = ("galfit " + FeedmeFileName)
	os.system(command)




	return 




