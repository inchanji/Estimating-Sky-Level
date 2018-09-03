import numpy as np 
import param as pm
from mpfit import mpfit

def startingParams(igal):
	i = igal - 1
	xc = int(round(float(pm.sexcatdata[i].split()[pm.iXc])))
	yc = int(round(float(pm.sexcatdata[i].split()[pm.iYc])))
	Re = float(pm.sexcatdata[i].split()[pm.iRe])

	return xc, yc, Re


def imgseg_cutout(xc, yc, Re):
	dx = dy = int(Re * pm.ReFac) 
	iterate = True 
	while iterate:
		x0 = int(xc - dx); x1 = int(xc + dx); y0 = int(yc - dy); y1 = int(yc + dy);
		if x0 < 0: x0 = 0; 
		if x1 > pm.Ncol: x1 = pm.Ncol;
		if y0 < 0: y0 = 0; 
		if y1 > pm.Nrow: y1 = pm.Now;
		img = pm.img[y0:y1+1, x0:x1+1].copy()
		seg = pm.newseg[y0:y1+1, x0:x1+1].copy()

		ind_not	 = np.logical_not(seg == 0) 
		seg[ind_not] = -1

		if (np.sum(seg != -1) > pm.NpixMin) and ((x1 - x0 + 1  > pm.LpixMin) or  (y1 - y0 + 1  > pm.LpixMin)): 
			iterate = False
		else: 
			dx += pm.npix_iter; dy += pm.npix_iter;

	return img, seg

def skyInitGuess(data, fac = 0.2):
	# Initial params are estimated using SExtractors' method
	sky = np.median(data); sig = np.std(data); N0  = len(data);
	ind = np.where(np.abs(data-sky) < 3. * sig)[0]
	data = data[ind].copy()
	N   = len(data); sig_ori = np.std(data)
	N0  = N + 1 
	Crowd = False
	while N != N0:
		N0  = N; sky = np.median(data); sig = np.std(data);
		ind = np.where(np.abs(data-sky) < 3. * sig)[0]
		data = data[ind].copy()
		sig = np.std(data); N = len(ind);
	if np.abs(sig-sig_ori) / sig_ori > fac: Crowd = True        
	if Crowd: return 2.5 * np.median(data) - 1.5 * np.mean(data), sig, True
	else:     return np.mean(data), sig, False	

def gauss(x, p):
	return p[0] * np.exp( -(x-p[1])**2. / (2.*p[2]**2.) )

def costfn(p, fjac = None, x = None, y = None, yerr = None):
	status = 0
	model  = gauss(x,p)
	return [ status, (y - model) / yerr  ] 	

def gaussFit( x = [], y = [], param = [], Quiet = 1):
	nparam  = len(param)
	fa = {'x':x, 'y':y, 'yerr':np.sqrt(y)}
	par_info = [{'value': 0., 'fixed': 0, 'limited': [0, 0], 'limits' : [0., 0.], 'tied' : ''} for i in range(nparam)]

	for i in range(nparam): 
		par_info[i]['value']      =  param[i]
		par_info[i]['limited']    =  [1, 1]

	par_info[0]['limits'] = [ param[0] * 1e-10      , param[0] / 0.5      ]  # amplitude
	par_info[1]['limits'] = [ param[1] - param[2] , param[1] + param[2] ]  # mean +/- sigma
	if param[2] > 200.: par_info[2]['value'] = 100.
	par_info[2]['limits'] = [ 1.      , 200.    ]  # standard devidation

	m     =     mpfit( costfn, parinfo = par_info, functkw = fa, quiet = Quiet )
	return m.params, m.fnorm, m.perror	


def estimate_sky(igal):
	xc, yc, Re  = startingParams(igal)
	img, seg 	= imgseg_cutout(xc, yc, Re)
	#print(np.shape(img))

	img[seg == -1] = np.nan
	data = img[np.isfinite(img)].ravel() 

	skySEx, skysigSEx, crowd =  skyInitGuess(data)

	print("Initial Guess: ",skySEx, skysigSEx, crowd)

	minval = np.max([skySEx - 5. * skysigSEx, -150]);   
	binsize =  3.5 * np.std(data) *  len(data) **(-0.33333333) * 1. # Bin size by Scott rule

	minval = round(minval*10.) / 10.; binsize = round(binsize*10.) / 10.; 
	maxval = skySEx + 5.* skysigSEx    
	maxval = minval + round((maxval-minval)/binsize) * binsize


	bins, edges = np.histogram(data, int((maxval - minval)/binsize), range = [minval,maxval])
	left, right = edges[:-1],edges[1:]
	x_hist 		= np.array([left,right],dtype = 'float').T.flatten(); 
	y_hist 		= np.array([bins,bins], dtype = 'float').T.flatten();    
	x = x_hist[::2] + 0.5 * (x_hist[1]-x_hist[0]); 
	y = y_hist[::2];  

	ind = np.where(y > 0)[0]
	x   = x[ind]; y = y[ind]

	param, xsqr, error =  gaussFit( x = x, y = y, param = [np.max(y), skySEx, skysigSEx] )

	return param[1], param[2], len(data) , np.mean(data), np.median(data)



