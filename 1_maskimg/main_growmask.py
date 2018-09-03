import sys
import numpy as np
import param as pm
import function as fn





def main(argv):
	print("Start masking, i.e., spatial filtering.")

	fn.growMask()

	print("Masking complete. \nNew segmentation map is saved at:{}".format(pm.path_newseg))

if __name__ == '__main__':
	main(sys.argv[1:])
