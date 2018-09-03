import sys
import numpy as np
import param as pm
import function as fn





def main(argv):
	print("Start estimating extragalactic background light (EBL).")

	fn.calEBL()

	print("Estimation done.")

if __name__ == '__main__':
	main(sys.argv[1:])
