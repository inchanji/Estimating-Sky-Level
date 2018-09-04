import sys
import numpy as np
import param as pm
import function as fn





def main(argv):
	if len(argv) == 0: 
		iGal  			= 1 		# default sexcat ID. the true seric index of the target gal. is ~ 2.4.
		skyval 			= 0.8541	# estimated sky background for ID = 1
		#skyval 		= 1.727129  # Sextractor

		pm.galfit.useGalfitSky = False 	# True: use Galfit's sky, False: use estimated sky value
	else: 
		iGal   			= int(argv[0])
		skyval 			= float(argv[1])
		if int(argv[2]) == 1: pm.galfit.useGalfitSky = True
		else: pm.galfit.useGalfitSky = False

	print("SBP fitting for SExtractor ID = {}".format(iGal))
	
	fn.sbpFitGal(iGal, skyval, pm.galfit.useGalfitSky)

	print("Fitting done.")

if __name__ == '__main__':
	main(sys.argv[1:])
