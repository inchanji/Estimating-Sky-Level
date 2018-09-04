# Pipeline for Galaxy Surface Brightness Profile Fitting.

These subroutines are used in Ji et al. (2018; see http://adsabs.harvard.edu/abs/2018PASP..130h4504J) for sky background estimation and galaxy surface brightness profile (SBP) fit. The original codes used in the paper were written in Python 2, but the codes uploaded here are modified to be operational in Python 3.

The steps to fit the SBP of galaxy are as follows:
1. Spatial filtering (1_maskimg)
2. Extragalactic Background Light (EBL) estimation (2_EBL): this task should be performed using a large field of CCD image to avoid sample vaiance. Because the size of the entire image (i.e., "random distribution simulation" in the paper) is so huge that it is impossible to upload to this repository. 
3. Statistical filtering (3_skyest): combined with spatial filter, this step yields local sky background (i.e., sky background without EBL correction).
4. SBP fit (4_sbpFit): This subroutine runs GALFIT (ver. 3) to estimate the best-fit parameters of Sersic profile. Parameters (MAG_AUTO, X_IMAGE, Y_IMAGE, A_IMAGE, etc.) and images (image cutout and rms image) described in "Feedme" file for GALFIT run are from SExtractor's measurement. 

For more details, please see Ji et al. (2018).
