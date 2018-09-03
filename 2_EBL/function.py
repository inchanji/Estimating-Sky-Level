import numpy as np
import param as pm


# PDF of galaxy number count
def pdf(x, l, mfaint):
	return 1./l * np.log10(x) + mfaint

# assumed galaxy number count 
def galNum(R):
	return 10.**(pm.mslope*R - pm.normconst)

# Caulculate the total number of galaxies using power law
def NtotGal(dmag = 0.2): 
	Ntot = 0; Num = []; mbright = pm.mbrightest - dmag / 2
	Niter = int((pm.mfaintest - pm.mbrightest) / dmag)
	for i in range(Niter):
		R = mbright + dmag * (i+1) # Rmag: 19 - 29
		Num.append(galNum(R))
		Ntot += Num[i]
	Ntot = int(Ntot)    
	print("The total number of galaxies assuming a power-law function : ", Ntot)
	return Ntot

# Estimate distribution of magnitudes
def distMag(Ntot):
	Rmags = []; n = 0;
	while n < Ntot:
		mag = pdf(np.random.uniform(), pm.mslope, pm.mfaintest)
		if (mag > pm.mbrightest) & (mag < pm.mfaintest):
			n += 1
			Rmags.append(mag)
	return np.array(Rmags)


def calEBL():
	Ntot 	= NtotGal()
	Rmags   = distMag(Ntot)

	F_EBL = 0.; m = pm.mcut; N_EBL = 0;

	while m < pm.mfaintest:
		ind = (pm.SExmaglist < m+pm.dmag) & (pm.SExmaglist > m)
		N_SEx = np.sum(ind)
		ind = np.where((Rmags < m+pm.dmag) & (Rmags > m))[0]

		mags_bin = Rmags[ind]
		np.random.shuffle(mags_bin)
		mags_bin2 = np.array(mags_bin)[N_SEx:]

		EBL_bin = np.sum( 10.**( (pm.m0 - mags_bin2) / 2.5 )) 

		N_EBL += len(mags_bin2)
		F_EBL += EBL_bin
		m 	  += pm.dmag
	

	Npix_outside_mask0 =  np.sum(pm.seg==0)
	Npix_outside_mask  =  np.sum(pm.newseg==0)
	data = pm.img[(pm.newseg == 0)].copy()

	print('The total flux of the EBL in the image is :', F_EBL) 
	print('The total number of galaxies undetected : ',  N_EBL)

	print ('Estimated EBL [ADU / pix] is :', F_EBL / Npix_outside_mask0)
	print ('The number of pixels outside grwon masks is :', Npix_outside_mask)
	print ('Mean flux per pixel outside grown masks : ', np.mean(data) )
	print ('Median flux per pixel outside grown masks : ', np.median(data) )
	
	

