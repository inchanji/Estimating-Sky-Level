import numpy as np
from astropy.io import fits
from mpfit import mpfit
import warnings
import time

from multiprocessing import Process
from multiprocessing.queues import Queue



class All:
	Path = '/Users/inchanji/Research/galsim/practice/ctr_set/github/'

	imgfilename   = 'SET6_variousPSF_DLSlike_sky3240_stack.fits'
	segfilename   = 'SET6_variousPSF_DLSlike_sky3240_stack_seg.fits'    
	catfilename   = 'SET6.cat'
	skyresultfile = 'SET6_SkyResult.txt'

	Ncpu 		  = 2

	MaxPos = 7999
	NCOMBINE = 20
	m0       = 32.73
	Modellist = ['mygalR24_avg.fits', 'mygalR25_avg.fits', \
				'mygalR26_avg.fits','mygalR27_avg.fits',\
				'mygalR28_avg.fits','mygalR29_avg.fits']
	Rmaglist = [ 26.1, 26.3, 26.5,  26.7,  26.9,  27.1,  27.3,  27.5,  \
				27.7,  27.9,  28.1,  28.3,  28.5,   28.7, 28.9]
	Ngallist = [ 1865, 4878, 9218, 14440, 19863, 25296, 30381, 36731, \
				44176, 53079, 64187, 76422, 92257, 111001, 133692]    

	HistFac  = (80./8000.)**2.
	Rmaglist = np.array(Rmaglist); Ngallist = np.array(Ngallist,dtype='float')
	Ngallist *=  HistFac    

	HistUGalsPix0 = []; 
	for i in range(len(Rmaglist)):
		Rmag = int(round(Rmaglist[i]))
		filename = 'mygalR'+str(Rmag)+'_avg.fits'
		hdu = fits.open(Path+filename)
		prof0 = hdu[0].data.ravel() 

		Ftot = 10.**((m0-Rmaglist[i])/2.5)
		prof0 = prof0[np.argsort(prof0)[::-1]]
		prof0 *= Ftot            
		prof0 = prof0[prof0 > 0.81] # 30 mag /arcsec^2
		for j in range(int(Ngallist[i])): 
			HistUGalsPix0 += list(prof0)
	HistUGalsPix0 = np.array(HistUGalsPix0)  


class Map:
	print 'Reading fits images ',
	hdu  = fits.open(All.Path + All.imgfilename)
	print '.',
	Img  = hdu[0].data.copy() # image
	print '.',
	hdu = fits.open(All.Path + All.segfilename)
	print '.',
	Seg = hdu[0].data.copy() # seg map
	print '.',
	del hdu
	print 'done'


def GalNum(R):
	return 10.**(0.4*R - 6.4347)

def return_mask_fac(totalmag):
	if totalmag < 20.:  fac = 15
	elif totalmag < 21: fac = 10
	elif totalmag < 22: fac = 8
	elif totalmag < 23: fac = 6
	else:               fac = 4
	return fac


def grow_mask(MaskMap, fac = 2, val = 1):
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

def clipping3sig(data):
	data = data.ravel()
	mu  = np.mean(data)
	sig = np.std(data)
	ind0 = data > mu + 3. * sig
	while np.sum(ind0):
		ind  = np.where( data < mu + 3. * sig )[0]
		data = data[ind].copy()
		mu   = np.mean(data)
		sig  = np.std(data)
		ind0 = data > mu + 3. * sig
	return mu, sig

def est_sky_SEx(data, fac = 0.2):
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

def EBL(x, X, Y, binsize, sky, noise):
	y = np.zeros(len(x))
	for i in range(len(X)):
		sig = np.sqrt(noise**2. + 1. * X[i]**2. / All.NCOMBINE)    
		amp = Y[i] * binsize / np.sqrt(2. * np.pi) / sig 
		if i == 0:
			y = amp * np.exp(-0.5 * (x - sky - X[i])**2./sig/sig)
		else:
			y += amp * np.exp(-0.5 * (x - sky - X[i])**2./sig/sig)    
	return y

def gauss(x, p):
	return p[0] * np.exp( -(x-p[1])**2. / (2.*p[2]**2.) )

def myskyfit(p, fjac = None, x = None, y = None, x2 = None, y2 = None, binsize= None):
	status  = 0        
	model   = gauss(x,p[:3]) + EBL(x, x2,y2,binsize,p[1], p[2]) * p[3]

	deviation = (y - model) / np.sqrt(y) 
	ind = y > 0.
	return [ status, deviation[ind] ] 

def skyfit( x = [], y = [],  x2 = [], y2 = [],param = [], binsize = 3.,Quiet = 1):
	nparam  = len(param)
	fa = {'x':x, 'y':y ,'x2':x2, 'y2':y2, 'binsize': binsize }
	par_info = [{'value': 0., 'fixed': 0, 'limited': [0, 0], 'limits' : [0., 0.], 'tied' : ''} \
			for i in range(nparam)]

	for i in range(nparam): 
		par_info[i]['value']      =  param[i]
		par_info[i]['limited']    =  [1, 1]

	par_info[0]['limits'] = [ param[0] * 1e-10      , param[0]  * 1.5    ]  # amplitude
	par_info[1]['limits'] = [ param[1] - param[2] , param[1] + param[2] ]  # mean +/- sigma
	if param[2] > 200.: par_info[2]['value'] = 100.
	par_info[2]['limits'] = [ 1.      , 200.    ]  # standard devidation
	par_info[3]['limits'] = [ param[3]*0.9999, param[3]*1.00001]

	m     =     mpfit( myskyfit, parinfo = par_info, functkw = fa, quiet = Quiet )
	return m.params, m.fnorm, m.perror


def estimate_sky(xc,yc,Re):
	dx = dy = int(Re*20.) 
	iterate = True 

	while iterate:
		x0 = xc - dx; x1 = xc + dx; y0 = yc - dy; y1 = yc + dy;
		if x0 < 0: x0 = 0 ; 
		if x1 > All.MaxPos: x1 = All.MaxPos;
		if y0 < 0: y0 = 0 ; 
		if y1 > All.MaxPos: y1 = All.MaxPos;
		img = Map.Img[y0:y1+1, x0:x1+1].copy()
		seg = Map.Seg[y0:y1+1, x0:x1+1].copy()

		Nmax = np.max(seg.ravel())
		for i in range(1, Nmax+1):
			if np.sum(seg.ravel() == i) == 0: continue
			flux = np.sum(img.ravel()[seg.ravel() == i])
			mag = All.m0 - 2.5*np.log10(flux)
			fac = return_mask_fac(mag) 
			seg = grow_mask(seg, fac, i)

		ind_not	 = np.logical_not(seg == 0) 
		seg[ind_not] = -1
		seg = grow_mask(seg, fac, -1)
		if (np.sum(seg != -1) > 4000) and ((x1 - x0 + 1  > 79) \
		    or  (y1 - y0 + 1  > 79)): iterate = False
		else: 
			dx += 5; dy += 5;

	img[seg == -1] = np.nan
	data = img[np.isfinite(img)].copy().ravel()
	mu_SEx, sig_SEx, Crowd =  est_sky_SEx(data.copy())
	minval = np.max([mu_SEx - 3. * sig_SEx, -150]);   

	data = img[np.isfinite(img)].copy().ravel()
	mu_SEx, sig_SEx, Crowd =  est_sky_SEx(data.copy())
	minval = np.max([mu_SEx - 3. * sig_SEx, -150]);   
	binsize =  3.5 * np.std(data) *  len(data) **(-0.33333333) * 1. # Bin size by Scott rule
	maxval = mu_SEx + 3.* sig_SEx    
	maxval = minval + round((maxval-minval)/binsize) * binsize

	bins, edges = np.histogram(data, int((maxval - minval)/binsize), \
	                           range = [minval,maxval])
	left, right = edges[:-1],edges[1:]
	x_hist = np.array([left,right],dtype = 'float').T.flatten(); 
	y_hist = np.array([bins,bins], dtype = 'float').T.flatten();    
	x = x_hist[::2] + 0.5 * (x_hist[1]-x_hist[0]); y = y_hist[::2];   


	bins, edges = np.histogram(All.HistUGalsPix0, int((maxval - minval)/binsize), \
	                           range = [minval,maxval])
	left, right = edges[:-1],edges[1:]
	x_hist = np.array([left,right],dtype = 'float').T.flatten(); 
	y_hist = np.array([bins,bins], dtype = 'float').T.flatten();    
	x2 = x_hist[::2] + 0.5 * (x_hist[1]-x_hist[0]); y2 = y_hist[::2];   
	amp_ratio = 1. * np.sum(np.isfinite(data)) / (80.)**2. 

	param, xsqr, error =  skyfit( x = x, y = y, x2 = x2, y2 = y2, \
	                    param = [np.max(y), mu_SEx, sig_SEx,  amp_ratio] )    
	return param[1],param[2]


def doFindSky(alldata, iproc, Nproc, index = [], queue = 0):
	i_st = 0
	i_end = len(index) - 1
	data    = alldata[index].copy()
	output = []
	print '> CPU ID = ', iproc , ', Index = [ ', i_st,' ', i_end,' ]'		
	dim 	= i_end - i_st + 1	
	for i in range(dim):
		
		ID = data[i].split()[0]
		xc = float(data[i].split()[1])
		yc = float(data[i].split()[2])
		Re = float(data[i].split()[5])
		sky, noise = estimate_sky(xc,yc,Re)
		
		output.append(ID+'\t'+str(sky) + '\t' + str(noise))
		print output[i]
	return queue.put(np.array(output))


def main():
	alldata   = np.array([line.rstrip('\n') for line in open(All.Path+All.catfilename)])[1:]
	Ndata  = len(alldata)

	if All.Ncpu > Ndata:
		print 'The number of cores > the total number of objects... '
		All.Ncpu = Ndata	

	indarr = [] ; index = np.arange(Ndata); 
	#np.random.shuffle(index)
	N = int( round(1. * Ndata / All.Ncpu));

	for i in range(All.Ncpu-1):
		ind 	= index[:N].copy()
		ind 	= ind[np.argsort(ind)].copy()
		indarr.append(ind)
		index 	= index[N:].copy()
	indarr.append(index[np.argsort(index)])	


	t_0 = time.time()
	queues = [Queue() for i in range(All.Ncpu)]
	args = [(alldata, i, All.Ncpu, indarr[i], queues[i]) for i in range(All.Ncpu)]
	jobs = [Process(target = doFindSky, args=(a)) for a in args]

	for j in jobs: j.start()

	data = queues[0].get()
	for q in queues[1:]: 
		data = np.concatenate(( data, q.get() ))

	for j in jobs: j.join()
	t_1 = time.time()

	print 'fitting is done.'
	print 'time :', t_1 - t_0, ' seconds, ', (t_1 - t_0) / 60. , ' minutes. ', All.Ncpu, ' CPUs'
	print '      ', (t_1 - t_0) * All.Ncpu / Ndata, 'seconds / obj'

	f = open(All.Path + All.skyresultfile,'w')
	for i in range(len(data)):
		f.write(data[i]+'\n')
	f.close()	


if __name__ == '__main__':
	main()
