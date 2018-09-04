# Pipeline for Galaxy Surface Brightness Profile Fitting (Ji et al. 2018).

These subroutines are used in Ji et al. (http://adsabs.harvard.edu/abs/2018PASP..130h4504J) for sky background estimation and galaxy surface brightness profile (SBP) fit. The original codes used in the paper were written in Python 2, but the codes uploaded here are modified to be operational in Python 3.

The steps to fit the SBP of galaxy are as follows:
1. Spatial masking (1_maskimg)
2. Extragalactic Background Light (EBL) estimation (2_EBL): This task should be performed using a large field of CCD image to avoid sample vaiance. Because the size of the entire image (i.e., "random distribution simulation" in the paper) is so huge that it is impossible to upload to this repository. 
