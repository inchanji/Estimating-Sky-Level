import sys
import numpy as np
import param as pm
import function as fn





def main(argv):
	if len(argv) == 0: iGal = 1 # default sexcat ID
	else: iGal = int(argv[0])

	print("Estimating sky background SExtractor ID = {}".format(iGal))
	
	skyval, skynoise, Npix, meanImg, medianImg = fn.estimate_sky(iGal)

	print("Estimated Sky Background: {} [ADU]".format(skyval))
	print("Nosie: {} [ADU]".format(skynoise))
	print("Number of pixels used: {} [pix]".format(Npix))
	print("Mean pixel value of outside masks: {} [ADU]".format(meanImg))
	print("Median pixel value of outside masks: {} [ADU]".format(medianImg))
	print("Estimation done.")

if __name__ == '__main__':
	main(sys.argv[1:])
