import numpy as np
import param as pm
#from mpfit import mpfit
from scipy import special
from astropy.io import fits

def concentration(n):
	return 1.9992 * n - 0.3271

def rlim(n, Ie, mu_lim):
	bn = concentration(n)
	return (1. - np.log(mu_lim/Ie)/bn)**n


def sersic(x,p):
	bn = concentration(p[2])
	return p[0] * np.exp(bn) * np.exp( -bn * (x/p[1])**(1./p[2]) )

'''
def costfn(p, fjac = None, x = None, y = None, yerr = None):
	status  = 0        
	model   = sersic(x,p)
	deviation = (y - model) / yerr
	return [ status, deviation ] 
'''

def Ieff(m,n,hlr) :
	L = 10.**((pm.m0-m)/2.5)
	bn = 1.9992*n-0.3271
	Ie = L * bn**(2.*n) / (2.*np.pi * hlr**2. * n * np.exp(bn) * special.gamma(2*n))
	return Ie
'''
def sersicfit( x = [], y = [], yerr = [], param = [], Quiet = 1):
	nparam = len(param)
	fa = {'x':x, 'y':y ,'yerr':yerr }
	par_info = [{'value': 0., 'fixed': 0, 'limited': [0, 0], 'lsmits' : [0., 0.], 'tied' : ''} \
		for i in range(nparam)]

	for i in range(nparam): 
		if param[i] < 0.: param[i] = np.abs(param[i])
		par_info[i]['value']      =  param[i]
		par_info[i]['limited']    =  [1, 1]

	par_info[0]['limits'] = [ param[0] * 1e-10     , param[0]  * 1e2    ] #  amplitude
	par_info[1]['limits'] = [ 0.05 , 3. ]  
	par_info[1]['fixed'] = 1
	par_info[2]['limits'] = [ 0.5      , 5.    ]  # standard devidation

	m     =     mpfit( costfn, parinfo = par_info, functkw = fa, quiet = Quiet )
	return m.params, m.fnorm, m.perror
'''


def growMask():
	index 	= np.arange(pm.i_seg_max) + 1; 
	seg 	= pm.seg.copy()

	for i in range(len(pm.sexcatdata)):
		ID =  index[i]
		if ID != int(pm.sexcatdata[ID-1].split()[0]):
			print("No matching... exit")

		x   	= int(float(pm.sexcatdata[ID-1].split()[1])); 
		y  		= int(float(pm.sexcatdata[ID-1].split()[2]));
		hlr 	= float(pm.sexcatdata[ID-1].split()[34]);
		cxx 	= float(pm.sexcatdata[ID-1].split()[35]); 	
		cyy 	= float(pm.sexcatdata[ID-1].split()[36]);
		cxy 	= float(pm.sexcatdata[ID-1].split()[37]);
		A  		= float(pm.sexcatdata[ID-1].split()[38]);
		KronR 	= float(pm.sexcatdata[ID-1].split()[43]);
		mag 	= float(pm.sexcatdata[ID-1].split()[9])	;


		Ie  = Ieff(mag, 4, hlr)
		fac = np.min([rlim(4,Ie,pm.mu_lim), 20]) # in units of hlr
		#fac = np.max([3.,fac])

		if np.isnan(fac): fac = 5. 

		cut = hlr * fac * 1.1
		lim = hlr * fac / A

		x0 = np.max([int(x - cut),0]);  x1 = np.min([int(x + cut + 1),pm.Ncol])
		y0 = np.max([int(y - cut),0]);  y1 = np.min([int(y + cut + 1),pm.Nrow])

		
		seg1 = seg[y0:y1,x0:x1].copy()

		y1d = np.arange(y1-y0) * 1. + y0  - y
		x1d = np.arange(x1-x0) * 1. + x0  - x

		x2d1, y2d1 = np.meshgrid(x1d,y1d)  
		
		ind1 = cxx * (x2d1)**2. + cyy * (y2d1)**2. + cxy * (x2d1) * (y2d1) <  lim**2.
		seg1[ind1] = ID
		
		seg[y0:y1,x0:x1] = seg1
		
	
	hdu = fits.PrimaryHDU(seg)
	hdu.writeto(pm.path_newseg, overwrite=True)



