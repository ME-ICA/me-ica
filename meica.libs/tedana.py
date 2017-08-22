#!/usr/bin/env python
__version__="v3.2 beta1"
welcome_block="""
# Multi-Echo ICA, Version %s
#
# Kundu, P., Brenowitz, N.D., Voon, V., Worbe, Y., Vertes, P.E., Inati, S.J., Saad, Z.S., 
# Bandettini, P.A. & Bullmore, E.T. Integrated strategy for improving functional 
# connectivity mapping using multiecho fMRI. PNAS (2013).
#
# Kundu, P., Inati, S.J., Evans, J.W., Luh, W.M. & Bandettini, P.A. Differentiating 
#   BOLD and non-BOLD signals in fMRI time series using multi-echo EPI. NeuroImage (2011).
# http://dx.doi.org/10.1016/j.neuroimage.2011.12.028
#
# PROCEDURE 2 : Computes ME-PCA and ME-ICA
# -Computes T2* map 
# -Computes PCA of concatenated ME data, then computes TE-dependence of PCs
# -Computes ICA of TE-dependence PCs
# -Identifies TE-dependent ICs, outputs high-\kappa (BOLD) component
#    and denoised time series
# -or- Computes TE-dependence of each component of a general linear model
#    specified by input (includes MELODIC FastICA mixing matrix)
""" % (__version__)

import os
import sys
from optparse import OptionParser
import numpy as np
import nibabel as nib
from sys import stdout,argv
import scipy.stats as stats
import time
import datetime
if __name__=='__main__':
	selfuncfile='%s/select_model_fft20e.py' % os.path.dirname(argv[0])
	execfile(selfuncfile)
import cPickle as pickle
import gzip
import pdb

if 'DEBUG' in argv: 
	import ipdb
	debug_mode = True
else: debug_mode = False

F_MAX=500
Z_MAX = 8

def _interpolate(a, b, fractionepsmap):
    """Returns the point at the given fraction between a and b, where
    'fraction' must be between 0 and 1.
    """
    return a + (b - a)*fraction;

def scoreatpercentile(a, per, limit=(), interpolation_method='lower'):
    """
    This function is grabbed from scipy

    """
    values = np.sort(a, axis=0)
    if limit:
        values = values[(limit[0] <= values) & (values <= limit[1])]

    idx = per /100. * (values.shape[0] - 1)
    if (idx % 1 == 0):
        score = values[idx]
    else:
        if interpolation_method == 'fraction':
            score = _interpolate(values[int(idx)], values[int(idx) + 1],
                                 idx % 1)
        elif interpolation_method == 'lower':
            score = values[int(np.floor(idx))]
        elif interpolation_method == 'higher':
            score = values[int(np.ceil(idx))]
        else:
            raise ValueError("interpolation_method can only be 'fraction', " \
                             "'lower' or 'higher'")
    return score

def MAD(data, axis=None):
    return np.median(np.abs(data - np.median(data, axis)), axis)

def dice(A,B):
	denom = np.array(A!=0,dtype=np.int).sum(0)+(np.array(B!=0,dtype=np.int).sum(0))
	if denom!=0:
		AB_un = andb([A!=0,B!=0])==2
		numer = np.array(AB_un,dtype=np.int).sum(0)
		return 2.*numer/denom
	else:
		return 0.

def spatclust(data,mask,csize,thr,header,aff,infile=None,dindex=0,tindex=0):
	if infile==None:
		data = data.copy()
		data[data<thr] = 0
		niwrite(unmask(data,mask),aff,'__clin.nii.gz',header)
		infile='__clin.nii.gz'
	addopts=""
	if data!=None and len(np.squeeze(data).shape)>1 and dindex+tindex==0: addopts="-doall"
	else: addopts="-1dindex %s -1tindex %s" % (str(dindex),str(tindex))
	os.system('3dmerge -overwrite %s -dxyz=1  -1clust 1 %i -1thresh %.02f -prefix __clout.nii.gz %s' % (addopts,int(csize),float(thr),infile))
	clustered = fmask(nib.load('__clout.nii.gz').get_data(),mask)!=0
	return clustered

def rankvec(vals):
	asort = np.argsort(vals)
	ranks = np.zeros(vals.shape[0])
	ranks[asort]=np.arange(vals.shape[0])+1
	return ranks

def niwrite(data,affine, name , header=None):
	stdout.write(" + Writing file: %s ...." % name) 
	
	thishead = header
	if thishead == None:
		thishead = head.copy()
		thishead.set_data_shape(list(data.shape))

	outni = nib.Nifti1Image(data,affine,header=thishead)
	outni.to_filename(name)
	print 'done.'

def cat2echos(data,Ne):
	"""
	cat2echos(data,Ne)

	Input:
	data shape is (nx,ny,Ne*nz,nt)
	"""
	nx,ny = data.shape[0:2]
	nz = data.shape[2]/Ne
	if len(data.shape) >3:
		nt = data.shape[3]
	else:
		nt = 1
	return np.reshape(data,(nx,ny,nz,Ne,nt),order='F')

def uncat2echos(data,Ne):
	"""
	uncat2echos(data,Ne)

	Input:
	data shape is (nx,ny,Ne,nz,nt)
	"""
    	nx,ny = data.shape[0:2]
	nz = data.shape[2]*Ne
	if len(data.shape) >4:
		nt = data.shape[4]
	else:
		nt = 1
	return np.reshape(data,(nx,ny,nz,nt),order='F')

def makemask(cdat):

	nx,ny,nz,Ne,nt = cdat.shape

	mask = np.ones((nx,ny,nz),dtype=np.bool)

	for i in range(Ne):
		tmpmask = (cdat[:,:,:,i,:] != 0).prod(axis=-1,dtype=np.bool)
		mask = mask & tmpmask

	return mask

def makeadmask(cdat,min=True,getsum=False):

	nx,ny,nz,Ne,nt = cdat.shape

	mask = np.ones((nx,ny,nz),dtype=np.bool)

	if min:
		mask = cdat[:,:,:,:,:].prod(axis=-1).prod(-1)!=0
		return mask
	else:
		#Make a map of longest echo that a voxel can be sampled with,
		#with minimum value of map as X value of voxel that has median
		#value in the 1st echo. N.b. larger factor leads to bias to lower TEs
		emeans = cdat[:,:,:,:,:].mean(-1)
		medv = emeans[:,:,:,0] == stats.scoreatpercentile(emeans[:,:,:,0][emeans[:,:,:,0]!=0],33,interpolation_method='higher')
		lthrs = np.squeeze(np.array([ emeans[:,:,:,ee][medv]/3 for ee in range(Ne) ]))
		if len(lthrs.shape)==1: lthrs = np.atleast_2d(lthrs).T
		lthrs = lthrs[:,lthrs.sum(0).argmax()]
		mthr = np.ones([nx,ny,nz,ne])
		for ee in range(Ne): mthr[:,:,:,ee]*=lthrs[ee]
		mthr = np.abs(emeans[:,:,:,:])>mthr
		masksum = np.array(mthr,dtype=np.int).sum(-1)
		mask = masksum!=0
		if getsum: return mask,masksum
		else: return mask

def fmask(data,mask):
	"""
	fmask(data,mask)

	Input:
	data shape is (nx,ny,nz,...)
	mask shape is (nx,ny,nz)

	Output:
	out shape is (Nm,...)
	"""

	s = data.shape
	sm = mask.shape

	N = s[0]*s[1]*s[2]
	news = []
	news.append(N)

	if len(s) >3:
		news.extend(s[3:])

	tmp1 = np.reshape(data,news)
	fdata = tmp1.compress((mask > 0 ).ravel(),axis=0)

	return fdata.squeeze()

def unmask (data,mask):
	"""
	unmask (data,mask)

	Input:

	data has shape (Nm,nt)
	mask has shape (nx,ny,nz)

	"""
	M = (mask != 0).ravel()
	Nm = M.sum()

	nx,ny,nz = mask.shape

	if len(data.shape) > 1:
		nt = data.shape[1]
	else:
		nt = 1

	out = np.zeros((nx*ny*nz,nt),dtype=data.dtype)
	out[M,:] = np.reshape(data,(Nm,nt))

	return np.squeeze(np.reshape(out,(nx,ny,nz,nt)))

def t2smap(catd,mask,tes):
	"""
	t2smap(catd,mask,tes)

	Input:

	catd  has shape (nx,ny,nz,Ne,nt)
	mask  has shape (nx,ny,nz)
	tes   is a 1d numpy array
	"""
	nx,ny,nz,Ne,nt = catd.shape
	N = nx*ny*nz

	echodata = fmask(catd,mask)
	Nm = echodata.shape[0]

	#Do Log Linear fit
	B = np.reshape(np.abs(echodata[:,:ne])+1, (Nm,(ne)*nt)).transpose()
	B = np.log(B)
	x = np.array([np.ones(Ne),-tes])
	X = np.tile(x,(1,nt))
	X = np.sort(X)[:,::-1].transpose()

	beta,res,rank,sing = np.linalg.lstsq(X,B)
	t2s = 1/beta[1,:].transpose()
	s0  = np.exp(beta[0,:]).transpose()

	#Goodness of fit
	alpha = (np.abs(B)**2).sum(axis=0)
	t2s_fit = blah = (alpha - res)/(2*res)
	
	out = np.squeeze(unmask(t2s,mask)),np.squeeze(unmask(s0,mask)),unmask(t2s_fit,mask)

	return out

def t2sadmap(catd,mask,tes):
	"""
	t2smap(catd,mask,tes)

	Input:

	catd  has shape (nx,ny,nz,Ne,nt)
	mask  has shape (nx,ny,nz)
	tes   is a 1d numpy array
	"""
	nx,ny,nz,Ne,nt = catd.shape
	N = nx*ny*nz

	echodata = fmask(catd,mask)
	Nm = echodata.shape[0]

	t2ss = np.zeros([nx,ny,nz,Ne-1])
	s0vs = np.zeros([nx,ny,nz,Ne-1])

	for ne in range(1,Ne+1):

		#Do Log Linear fit
		B = np.reshape(np.abs(echodata[:,:ne])+1, (Nm,(ne)*nt)).transpose()
		B = np.log(B)
		x = np.array([np.ones(ne),-tes[:ne] ])
		X = np.tile(x,(1,nt))
		X = np.sort(X)[:,::-1].transpose()

		beta,res,rank,sing = np.linalg.lstsq(X,B)
		t2s = 1/beta[1,:].transpose()
		s0  = np.exp(beta[0,:]).transpose()

		t2s[np.isinf(t2s)] = 500.
		s0[np.isnan(s0)] = 0.

		t2ss[:,:,:,ne-2] = np.squeeze(unmask(t2s,mask))
		s0vs[:,:,:,ne-2] = np.squeeze(unmask(s0,mask))

	#Limited T2* and S0 maps
	fl = np.zeros([nx,ny,nz,len(tes)-2+1])
	for ne in range(Ne-1):
		fl_ = np.squeeze(fl[:,:,:,ne])
		fl_[masksum==ne+2] = True
		fl[:,:,:,ne] = fl_
	fl = np.array(fl,dtype=bool)
	t2sa = np.squeeze(unmask(t2ss[fl],masksum>1))
	s0va = np.squeeze(unmask(s0vs[fl],masksum>1))

	#Full T2* maps with S0 estimation errors
	t2saf = t2sa.copy()
	s0vaf = s0va.copy()
	t2saf[masksum==1] = t2ss[masksum==1,0]
	s0vaf[masksum==1] = s0vs[masksum==1,0]

	return t2sa,s0va,t2ss,s0vs,t2saf,s0vaf

def get_coeffs(data,mask,X,add_const=False):
	"""
	get_coeffs(data,X)

	Input:

	data has shape (nx,ny,nz,nt)
	mask has shape (nx,ny,nz)
	X    has shape (nt,nc)

	Output:

	out  has shape (nx,ny,nz,nc)
	""" 
	mdata = fmask(data,mask).transpose()
        
        X=np.atleast_2d(X)
        if X.shape[0]==1: X=X.T
        Xones = np.atleast_2d(np.ones(np.min(mdata.shape))).T
        if add_const: X = np.hstack([X,Xones])

	tmpbetas = np.linalg.lstsq(X,mdata)[0].transpose()
        if add_const: tmpbetas = tmpbetas[:,:-1]
	out = unmask(tmpbetas,mask)

	return out


def andb(arrs):
	result = np.zeros(arrs[0].shape)
	for aa in arrs: result+=np.array(aa,dtype=np.int)
	return result


def optcom(data,t2s,tes,mask,useG=True):
	"""
	out = optcom(data,t2s)


	Input:

	data.shape = (nx,ny,nz,Ne,Nt)
	t2s.shape  = (nx,ny,nz)
	tes.shape  = (Ne,)

	Output:

	out.shape = (nx,ny,nz,Nt)
	"""
	nx,ny,nz,Ne,Nt = data.shape 

	if useG:
		fdat = fmask(data,mask)
		ft2s = fmask(t2sG,mask)

	else:
		fdat = fmask(data,mask)
		ft2s = fmask(t2s,mask)
	
	tes = tes[np.newaxis,:]
	ft2s = ft2s[:,np.newaxis]
	
	if options.combmode == 'ste':
		alpha = fdat.mean(-1)*tes
	else: 
		alpha = tes * np.exp(-tes /ft2s)
	
	alpha = np.tile(alpha[:,:,np.newaxis],(1,1,Nt))

	fout  = np.average(fdat,axis = 1,weights=alpha)
	out = unmask(fout,mask)
	print 'Out shape is ', out.shape
	return out

def getelbow(ks,val=False):
	#Elbow using linear projection method - moderate
	ks = np.sort(ks)[::-1]
	nc = ks.shape[0]
	coords = np.array([np.arange(nc),ks])
	p  = coords - np.tile(np.reshape(coords[:,0],(2,1)),(1,nc))
	b  = p[:,-1] 
	b_hat = np.reshape(b/np.sqrt((b**2).sum()),(2,1))
	proj_p_b = p - np.dot(b_hat.T,p)*np.tile(b_hat,(1,nc))
	d = np.sqrt((proj_p_b**2).sum(axis=0))
	k_min_ind = d.argmax()
	k_min  = ks[k_min_ind]
	if val: return ks[k_min_ind]
	else: return k_min_ind

def getelbow2(ks,val=False):
	#Elbow using mean/variance method - conservative
	ks = np.sort(ks)[::-1]
	nk = len(ks)
	ds = np.array([  (ks[nk-5-ii-1]>ks[nk-5-ii:nk].mean()+2*ks[nk-5-ii:nk].std()) for ii in range(nk-5) ][::-1],dtype=np.int)
	dsum = []
	c_ = 0
	for d_ in ds: 
		c_=(c_+d_)*d_
		dsum.append(c_)
	e2=np.argmax(np.array(dsum))
	elind = np.max([getelbow(ks),e2])
	if val: return ks[elind]
	else: return elind

def getelbow3(ks,val=False):
	#Elbow using curvature - aggressive
	ks = np.sort(ks)[::-1]
	dKdt = ks[:-1]-ks[1:]
	dKdt2 = dKdt[:-1]-dKdt[1:]		
	curv = np.abs((dKdt2/(1+dKdt[:-1]**2.)**(3./2.)))
	curv[np.isnan(curv)]=-1*10**6
	maxcurv = np.argmax(curv)+2
	if val: return(ks[maxcurv])
	else:return maxcurv

def getfbounds(ne):
	F05s=[None,None,18.5,10.1,7.7,6.6,6.0,5.6,5.3,5.1,5.0]
	F025s=[None,None,38.5,17.4,12.2,10,8.8,8.1,7.6,7.2,6.9]
	F01s=[None,None,98.5,34.1,21.2,16.2,13.8,12.2,11.3,10.7,10.]
	return F05s[ne-1],F025s[ne-1],F01s[ne-1]

def eimask(dd,ees=None):
	if ees==None: ees=range(dd.shape[1])
	imask = np.zeros([dd.shape[0],len(ees)])
	for ee in ees:
		print ee
		lthr = 0.001*scoreatpercentile(dd[:,ee,:].flatten(),98)
		hthr = 5*scoreatpercentile(dd[:,ee,:].flatten(),98)
		print lthr,hthr
		imask[dd[:,ee,:].mean(1) > lthr,ee]=1
		imask[dd[:,ee,:].mean(1) > hthr,ee]=0
	return imask

def dwtmat(mmix):
    #ipdb.set_trace()
    print "++Wavelet transforming data"
    llt = len(np.hstack(pywt.dwt(mmix[0],'db2')))
    mmix_wt = np.zeros([mmix.shape[0],llt])
    for ii in range(mmix_wt.shape[0]):
        wtx = pywt.dwt(mmix[ii],'db2')
        #print len(wtx[0]),len(wtx[1])
        cAlen = len(wtx[0])
        mmix_wt[ii] = np.hstack(wtx)
    return mmix_wt,cAlen

def idwtmat(mmix_wt,cAl):
    print "++Inverse wavelet transforming"
    lt = len(pywt.idwt(mmix_wt[0,:cAl],mmix_wt[0,cAl:],'db2',correct_size=True))
    mmix_iwt = np.zeros([mmix_wt.shape[0],lt])
    for ii in range(mmix_iwt.shape[0]):
        mmix_iwt[ii] = pywt.idwt(mmix_wt[ii,:cAl],mmix_wt[ii,cAl:],'db2',correct_size=True)
    return mmix_iwt

def tedpca(ste=0,mlepca=True): 
	nx,ny,nz,ne,nt = catd.shape
	ste = np.array([int(ee) for ee in str(ste).split(',')])
	cAl = None
	if len(ste) == 1 and ste[0]==-1:
		print "-Computing PCA of optimally combined multi-echo data"
		OCmask = makemask(OCcatd[:,:,:,np.newaxis,:])
		d = fmask(OCcatd,OCmask)
		eim = eimask(d[:,np.newaxis,:])
		eim = eim[:,0]==1
		d = d[eim,:]
		#ipdb.set_trace()
	elif len(ste) == 1 and ste[0]==0:
		print "-Computing PCA of spatially concatenated multi-echo data"
		ste = np.arange(ne)
		d = np.float64(fmask(catd,mask))
		eim = eimask(d)==1
		d = d[eim]
	else:
		print "-Computing PCA of TE #%s" % ','.join([str(ee) for ee in ste])
		d = np.float64(np.concatenate([fmask(catd[:,:,:,ee,:],mask)[:,np.newaxis,:] for ee in ste-1],axis=1))
		eim = eimask(d)==1
		eim = np.squeeze(eim)
		d = np.squeeze(d[eim])

	dz = ((d.T-d.T.mean(0))/d.T.std(0)).T #Variance normalize timeseries
	dz = (dz-dz.mean())/dz.std() #Variance normalize everything

	#if wvpca: 
	#	import ipdb
	#	ipdb.set_trace()
	#	print "++Transforming time series from time domain to wavelet domain."
	#	dz,cAl = dwtmat(dz)

	pcastate_fn = 'pcastate.pklgz'

	if not os.path.exists(pcastate_fn):
		##Do PC dimension selection
		#Get eigenvalue cutoff

		if mlepca:
			from sklearn.decomposition import PCA
			ppca = PCA(n_components='mle')
			ppca.fit(dz)
			v = ppca.components_
			s = ppca.explained_variance_
			u = np.dot(np.dot(dz,v.T),np.diag(1./s))
		else:
			u,s,v = np.linalg.svd(dz,full_matrices=0)
		
		sp = s/s.sum()
		eigelb = sp[getelbow(sp)]

		spdif = np.abs(sp[1:]-sp[:-1])
		spdifh = spdif[spdif.shape[0]/2:]
		spdmin = spdif.min()
		spdthr = np.mean([spdifh.max(),spdmin])
		spmin = sp[(spdif.shape[0]/2)+(np.arange(spdifh.shape[0])[spdifh>=spdthr][0])+1]
		spcum = []
		spcumv = 0
		for sss in sp:
			spcumv+=sss
			spcum.append(spcumv)
		spcum = np.array(spcum)
			
		#Compute K and Rho for PCA comps
		
		#ipdb.set_trace()

		eimum = np.atleast_2d(eim)
		eimum = np.transpose(eimum,np.argsort(np.atleast_2d(eim).shape)[::-1])
		eimum = np.array(np.squeeze(unmask(eimum.prod(1),mask)),dtype=np.bool)
		vTmix = v.T
		vTmixN =((vTmix.T-vTmix.T.mean(0))/vTmix.T.std(0)).T
		#ctb,KRd,betasv,v_T = fitmodels2(catd,v.T,eimum,t2s,tes,mmixN=vTmixN)
		none,ctb,betasv,v_T = fitmodels_direct(catd,v.T,eimum,t2s,tes,mmixN=vTmixN,full_sel=False)
		ctb = ctb[ctb[:,0].argsort(),:]
		ctb = np.vstack([ctb.T[0:3],sp]).T

		#Save state
		print "Saving PCA"
		pcastate = {'u':u,'s':s,'v':v,'ctb':ctb,'eigelb':eigelb,'spmin':spmin,'spcum':spcum}
		try:		
			pcastate_f = gzip.open(pcastate_fn,'wb')
			pickle.dump(pcastate,pcastate_f)
			pcastate_f.close()
		except:
			print "Could not save PCA solution!"

	else:
		print "Loading PCA"
		pcastate_f = gzip.open(pcastate_fn,'rb')
		pcastate = pickle.load(pcastate_f)
		for key,val in pcastate.items(): exec(key + '=val')

	np.savetxt('comp_table_pca.txt',ctb[ctb[:,1].argsort(),:][::-1])
	np.savetxt('mepca_mix.1D',v[ctb[:,1].argsort()[::-1],:].T)
	
	kappas = ctb[ctb[:,1].argsort(),1]
	rhos = ctb[ctb[:,2].argsort(),2]
	fmin,fmid,fmax = getfbounds(ne)
	kappa_thr = np.average(sorted([fmin,getelbow(kappas,True)/2,fmid]),weights=[kdaw,1,1])
	rho_thr = np.average(sorted([fmin,getelbow2(rhos,True)/2,fmid]),weights=[rdaw,1,1])
	if int(kdaw)==-1:
		kappas_lim = kappas[andb([kappas<fmid,kappas>fmin])==2]
		#kappas_lim = kappas[andb([kappas<kappas[getelbow(kappas)],kappas>fmin])==2]
		kappa_thr = kappas_lim[getelbow(kappas_lim)]
		rhos_lim = rhos[andb([rhos<fmid,rhos>fmin])==2]
		rho_thr = rhos_lim[getelbow(rhos_lim)]
		options.stabilize=True
	if int(kdaw)!=-1 and int(rdaw)==-1:
		rhos_lim = rhos[andb([rhos<fmid,rhos>fmin])==2]
		rho_thr = rhos_lim[getelbow(rhos_lim)]
	if options.stabilize:
		pcscore = (np.array(ctb[:,1]>kappa_thr,dtype=np.int)+np.array(ctb[:,2]>rho_thr,dtype=np.int)+np.array(ctb[:,3]>eigelb,dtype=np.int))*np.array(ctb[:,3]>spmin,dtype=np.int)*np.array(spcum<0.95,dtype=np.int)*np.array(ctb[:,2]>fmin,dtype=np.int)*np.array(ctb[:,1]>fmin,dtype=np.int)*np.array(ctb[:,1]!=F_MAX,dtype=np.int)*np.array(ctb[:,2]!=F_MAX,dtype=np.int) 
	else:
		pcscore = (np.array(ctb[:,1]>kappa_thr,dtype=np.int)+np.array(ctb[:,2]>rho_thr,dtype=np.int)+np.array(ctb[:,3]>eigelb,dtype=np.int))*np.array(ctb[:,3]>spmin,dtype=np.int)*np.array(ctb[:,1]!=F_MAX,dtype=np.int)*np.array(ctb[:,2]!=F_MAX,dtype=np.int)
	pcsel = pcscore > 0 
	pcrej = np.array(pcscore==0,dtype=np.int)*np.array(ctb[:,3]>spmin,dtype=np.int) > 0

	"""
kdaw=15
rdaw=5
kappa_thr = np.average(sorted([fmin,getelbow(kappas,True)/2,fmid]),weights=[kdaw,1,1])
rho_thr = np.average(sorted([fmin,getelbow2(rhos,True)/2,fmid]),weights=[rdaw,1,1])
print kappa_thr, rho_thr
pcscore = (np.array(ctb[:,1]>kappa_thr,dtype=np.int)+np.array(ctb[:,2]>rho_thr,dtype=np.int)+np.array(ctb[:,3]>eigelb,dtype=np.int))*np.array(ctb[:,3]>spmin,dtype=np.int)*np.array(ctb[:,1]!=F_MAX,dtype=np.int)*np.array(ctb[:,2]!=F_MAX,dtype=np.int)
np.array(pcscore>0).sum()
	"""

	if 'DEBUG' in sys.argv:
		import ipdb
		ipdb.set_trace()

	dd = u.dot(np.diag(s*np.array(pcsel,dtype=np.int))).dot(v)
	
	#if wvpca: 
	#	print "++Transforming PCA solution from wavelet domain to time domain"
	#	dd = idwtmat(dd,cAl)

	nc = s[pcsel].shape[0]
	print pcsel
	print "--Selected %i components. Minimum Kappa=%0.2f Rho=%0.2f" % (nc,kappa_thr,rho_thr)

	dd = ((dd.T-dd.T.mean(0))/dd.T.std(0)).T #Variance normalize timeseries
	dd = (dd-dd.mean())/dd.std() #Variance normalize everything

	return nc,dd

def tedica(dd,cost):
	"""
	Input is dimensionally reduced spatially concatenated multi-echo time series dataset from tedpca()
	Output is comptable, mmix, smaps from ICA, and betas from fitting catd to mmix
	"""
	#Do ICA
	climit = float("%s" % options.conv)
	#icanode = mdp.nodes.FastICANode(white_comp=nc, white_parm={'svd':True},approach='symm', g=cost, fine_g=options.finalcost, limit=climit, verbose=True)
	icanode = mdp.nodes.FastICANode(white_comp=nc,approach='symm', g=cost, fine_g=options.finalcost, primary_limit=climit*100, limit=climit, verbose=True)
	icanode.train(dd)
	smaps = icanode.execute(dd)
	mmix = icanode.get_recmatrix().T
	mmix = (mmix-mmix.mean(0))/mmix.std(0)
	return mmix

def write_split_ts(data,comptable,mmix,suffix=''):
	mdata = fmask(data,mask)
	betas = fmask(get_coeffs(unmask((mdata.T-mdata.T.mean(0)).T,mask),mask,mmix),mask)
	dmdata = mdata.T-mdata.T.mean(0)
	varexpl = (1-((dmdata.T-betas.dot(mmix.T))**2.).sum()/(dmdata**2.).sum())*100
	print 'Variance explained: ', varexpl , '%'
	midkts = betas[:,midk].dot(mmix.T[midk,:])
	lowkts = betas[:,rej].dot(mmix.T[rej,:])
	if len(acc)!=0:
		niwrite(unmask(betas[:,acc].dot(mmix.T[acc,:]),mask),aff,'_'.join(['hik_ts',suffix])+'.nii')
	if len(midk)!=0: 
		niwrite(unmask(midkts,mask),aff,'_'.join(['midk_ts',suffix])+'.nii')
	if len(rej)!=0: 
		niwrite(unmask(lowkts,mask),aff,'_'.join(['lowk_ts',suffix])+'.nii')
	niwrite(unmask(fmask(data,mask)-lowkts-midkts,mask),aff,'_'.join(['dn_ts',suffix])+'.nii')
	return varexpl

def split_ts(data,comptable,mmix):
	cbetas = get_coeffs(data-data.mean(-1)[:,:,:,np.newaxis],mask,mmix)
	betas = fmask(cbetas,mask)
	if len(acc)!=0:
		hikts=unmask(betas[:,acc].dot(mmix.T[acc,:]),mask)
	else:
		hikts = None
	return hikts,data-hikts

def writefeats(cbetas,comptable,mmix,suffix=''):
	#Write signal changes (dS)
	niwrite(cbetas[:,:,:,:],aff,'_'.join(['betas',suffix])+'.nii')
	niwrite(cbetas[:,:,:,acc],aff,'_'.join(['betas_hik',suffix])+'.nii')
	#Compute features (dS/S)
	if options.e2d==None: e2d=np.floor(ne/2)+1
	edm = fmask(catd[:,:,:,e2d-1,:],mask)
	edms = edm/edm.std(-1)[:,np.newaxis]
	edms[edm<1]=0
	hik,noise = split_ts(unmask(edms,mask),comptable,mmix)
	noise = noise-noise.mean(-1)[:,:,:,np.newaxis]

	zfac = 1./(mmix.shape[0]-len(acc)-1)*(noise**2).sum(-1) #noise scaling
	niwrite(zfac,aff,'zfac.nii')
	
	cbetam = fmask(cbetas[:,:,:,acc],mask)
	cbetam = (cbetam-cbetam.mean(0))/cbetam.std(0)
	cbetam = cbetam/fmask(zfac,mask)[:,np.newaxis]
	cbetam[edm.mean(-1)<1,:] = 0
	
	niwrite(unmask(cbetam,mask),aff,'_'.join(['feats',suffix])+'.nii')
	
def computefeats2(data,mmix,mask,normalize=True):
	#Write feature versions of components 
	data = data[mask]
	data_vn = (data-data.mean(axis=-1)[:,np.newaxis])/data.std(axis=-1)[:,np.newaxis]
	data_R = get_coeffs(unmask(data_vn,mask),mask,mmix)[mask]
	data_R[data_R<-.999] = -0.999
	data_R[data_R>.999] = .999
	data_Z = np.arctanh(data_R)
	if len(data_Z.shape)==1: data_Z = np.atleast_2d(data_Z).T
	if normalize:
		#data_Z2 = ((data_Z.T-data_Z.mean(0)[:,np.newaxis])/data_Z.std(0)[:,np.newaxis]).T
		data_Z = (((data_Z.T-data_Z.mean(0)[:,np.newaxis])/data_Z.std(0)[:,np.newaxis])  + (data_Z.mean(0)/data_Z.std(0))[:,np.newaxis]).T
	return data_Z

def writefeats2(data,mmix,mask,suffix=''):
	#Write feature versions of components 
	feats = computefeats2(data,mmix,mask)
	niwrite(unmask(feats,mask),aff,'_'.join(['feats',suffix])+'.nii')

def writect(comptable,ctname='',varexpl='-1',classarr=[]):
	global acc,rej,midk,empty
	if len(classarr)!=0:
		acc,rej,midk,empty = classarr
	nc = comptable.shape[0]
	ts = time.time()
	st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
	sortab = comptable[comptable[:,1].argsort()[::-1],:]
	if ctname=='': ctname = 'comp_table.txt'
	open('accepted.txt','w').write(','.join([str(int(cc)) for cc in acc]))
	open('rejected.txt','w').write(','.join([str(int(cc)) for cc in rej]))
	open('midk_rejected.txt','w').write(','.join([str(int(cc)) for cc in midk]))
	with open(ctname,'w') as f:
		f.write("#\n#ME-ICA Component statistics table for: %s \n#Run on %s \n#\n" % (os.path.abspath(os.path.curdir),st) )
		f.write("#Dataset variance explained by ICA (VEx): %.02f \n" %  ( varexpl ) )
		f.write("#Total components generated by decomposition (TCo): %i \n" %  ( nc ) )
		f.write("#No. accepted BOLD-like components, i.e. effective degrees of freedom for correlation (lower bound; DFe): %i\n" %  ( len(acc) ) )
		f.write("#Total number of rejected components (RJn): %i\n" %  (len(midk)+len(rej)) )
		f.write("#Nominal degress of freedom in denoised time series (..._medn.nii.gz; DFn): %i \n" %  (nt-len(midk)-len(rej)) )
		f.write("#ACC %s \t#Accepted BOLD-like components\n" % ','.join([str(int(cc)) for cc in acc]) )
		f.write("#REJ %s \t#Rejected non-BOLD components\n" % ','.join([str(int(cc)) for cc in rej]) )
		f.write("#MID %s \t#Rejected R2*-weighted artifacts\n" % ','.join([str(int(cc)) for cc in midk]) )
		f.write("#IGN %s \t#Ignored components (kept in denoised time series)\n" % ','.join([str(int(cc)) for cc in empty]) )
		f.write("#VEx	TCo	DFe	RJn	DFn	\n")
		f.write("##%.02f	%i	%i	%i	%i \n" % (varexpl,nc,len(acc),len(midk)+len(rej),nt-len(midk)-len(rej)))
		f.write("#	comp	Kappa	Rho	%%Var	%%Var(norm)	\n")
		for i in range(nc):
			f.write('%d\t%f\t%f\t%.2f\t%.2f\n'%(sortab[i,0],sortab[i,1],sortab[i,2],sortab[i,3],sortab[i,4]))

def gscontrol_raw(dtrank=4):
	"""
	This function uses the spatial global signal estimation approach to modify catd (global variable) to 
	removal global signal out of individual echo time series datasets. The spatial global signal is estimated
	from the optimally combined data after detrending with a Legendre polynomial basis of order=0 and degree=dtrank.
	"""

	print "++ Applying amplitude-based T1 equilibration correction"
    
	#Legendre polynomial basis for denoising
	from scipy.special import lpmv
	Lmix = np.array([lpmv(0,vv,np.linspace(-1,1,OCcatd.shape[-1])) for vv in range(dtrank)]).T

	#Compute mean, std, mask local to this function - inefficient, but makes this function a bit more modular
	Gmu = OCcatd.mean(-1)
	Gstd = OCcatd.std(-1)
	Gmask = Gmu!=0

	#Find spatial global signal
	dat = OCcatd[Gmask] - Gmu[Gmask][:,np.newaxis]
	sol = np.linalg.lstsq(Lmix,dat.T) #Legendre basis for detrending
	detr = dat - np.dot(sol[0].T,Lmix.T)[0]
	sphis = (detr).min(1)
	sphis -= sphis.mean()
	niwrite(unmask(sphis,Gmask),aff,'T1gs.nii',head)
    
	#Find time course of the spatial global signal, make basis with the Legendre basis
	glsig = np.linalg.lstsq(np.atleast_2d(sphis).T,dat)[0]
	glsig = (glsig-glsig.mean() ) / glsig.std()
	np.savetxt('glsig.1D',glsig)
	glbase = np.hstack([Lmix,glsig.T])

	#Project global signal out of optimally combined data
	sol = np.linalg.lstsq(np.atleast_2d(glbase),dat.T)
	tsoc_nogs = dat - np.dot(np.atleast_2d(sol[0][dtrank]).T,np.atleast_2d(glbase.T[dtrank])) + Gmu[Gmask][:,np.newaxis]

	global OCcatd, sphis_raw
	sphis_raw = sphis
	niwrite(OCcatd,aff,'tsoc_orig.nii',head)
	OCcatd = unmask(tsoc_nogs,Gmask)
	niwrite(OCcatd,aff,'tsoc_nogs.nii',head)
	
	#Project glbase out of each echo
	for ii in range(Ne):
		dat = catd[:,:,:,ii,:][Gmask]
		sol = np.linalg.lstsq(np.atleast_2d(glbase),dat.T)
		e_nogs = dat - np.dot(np.atleast_2d(sol[0][dtrank]).T,np.atleast_2d(glbase.T[dtrank]))
		catd[:,:,:,ii,:] = unmask(e_nogs,Gmask)
		if 'DEBUG' in sys.argv:
			ipdb.set_trace()

def gscontrol_mmix():

	Gmu = OCcatd.mean(-1)
	Gstd = OCcatd.std(-1)
	Gmask = Gmu!=0

	"""
	Compute temporal regression
	"""
	dat = (OCcatd[Gmask] - Gmu[Gmask][:,np.newaxis]) / Gstd[mask][:,np.newaxis]
	solG = np.linalg.lstsq(mmix,dat.T)
	resid = dat - np.dot(solG[0].T,mmix.T)

	"""
	Build BOLD time series without amplitudes, and save T1-like effect
	"""
	bold_ts = np.dot(solG[0].T[:,acc],mmix[:,acc].T)
	sphis = bold_ts.min(-1)
	sphis -= sphis.mean()
	print sphis.shape
	niwrite(unmask(sphis,mask),aff,'sphis_hik.nii',header=head)

	"""
	Find the global signal based on the T1-like effect
	"""
	sol = np.linalg.lstsq(np.atleast_2d(sphis).T,dat)
	glsig = sol[0]

	"""
	T1 correct time series by regression
	"""
	bold_noT1gs = bold_ts - np.dot(np.linalg.lstsq(glsig.T,bold_ts.T)[0].T,glsig)
	niwrite(unmask(bold_noT1gs*Gstd[mask][:,np.newaxis],mask),aff,'hik_ts_OC_T1c.nii',header=head)

	"""
	Make medn version of T1 corrected time series
	"""
	niwrite(Gmu[:,:,:,np.newaxis]+unmask((bold_noT1gs+resid)*Gstd[mask][:,np.newaxis] ,mask),aff,'dn_ts_OC_T1c.nii',header=head)

	"""
	Orthogonalize mixing matrix w.r.t. T1-GS
	"""
	mmixnogs = mmix.T - np.dot(np.linalg.lstsq(glsig.T,mmix)[0].T,glsig)
	mmixnogs_mu = mmixnogs.mean(-1)
	mmixnogs_std = mmixnogs.std(-1)
	mmixnogs_norm = (mmixnogs-mmixnogs_mu[:,np.newaxis])/mmixnogs_std[:,np.newaxis]
	mmixnogs_norm = np.vstack([np.atleast_2d(np.ones(max(glsig.shape))),glsig,mmixnogs_norm])
	
	"""
	Write T1-GS corrected components and mixing matrix
	"""
	sol = np.linalg.lstsq(mmixnogs_norm.T,dat.T)
	niwrite(unmask(sol[0].T[:,2:],mask),aff,'betas_hik_OC_T1c.nii',header=head)
	np.savetxt('meica_mix_T1c.1D',mmixnogs)


def dwtmat(mmix):
	if 'DEBUG' in sys.argv:
		import ipdb
		ipdb.set_trace()
	print "++Wavelet transforming data"
	llt = len(np.hstack(pywt.dwt(mmix[0],'db2')))
	mmix_wt = np.zeros([mmix.shape[0],llt])
	for ii in range(mmix_wt.shape[0]):
		wtx = pywt.dwt(mmix[ii],'db2')
		#print len(wtx[0]),len(wtx[1])
		cAlen = len(wtx[0])
		mmix_wt[ii] = np.hstack(wtx)
	return mmix_wt,cAlen

def idwtmat(mmix_wt,cAl):
	print "++Inverse wavelet transforming"
	lt = len(pywt.idwt(mmix_wt[0,:cAl],mmix_wt[0,cAl:],'db2',correct_size=True))
	mmix_iwt = np.zeros([mmix_wt.shape[0],lt])
	for ii in range(mmix_iwt.shape[0]):
		mmix_iwt[ii] = pywt.idwt(mmix_wt[ii,:cAl],mmix_wt[ii,cAl:],'db2',correct_size=True)
	return mmix_iwt

def dwtcatd():
	if 'DEBUG' in sys.argv:
		import ipdb
		ipdb.set_trace()
	stackmask = tile(mask,[1,1,ne])
	maskin = uncat2echos(catd,ne)[stackmask]
	trans_mask_catd,newlen = dwtmat(uncat2echos(catd,ne)[stackmask])
	newcatd = unmask(trans_mask_catd,stackmask)

def writeresults():
	print "++ Writing optimally combined time series"
	ts = OCcatd
	niwrite(ts,aff,'ts_OC.nii')
	print "++ Writing Kappa-filtered optimally combined timeseries"
	varexpl = write_split_ts(ts,comptable,mmix,'OC')
	print "++ Writing signal versions of components"
	ts_B = get_coeffs(ts,mask,mmix)
	niwrite(ts_B[:,:,:,:],aff,'_'.join(['betas','OC'])+'.nii')
	if len(acc)!=0:
		niwrite(ts_B[:,:,:,acc],aff,'_'.join(['betas_hik','OC'])+'.nii')
		print "++ Writing optimally combined high-Kappa features"
		writefeats2(split_ts(ts,comptable,mmix)[0],mmix[:,acc],mask,'OC2')
	print "++ Writing component table"
	writect(comptable,'comp_table.txt',varexpl)

def writeresults_echoes():
	for ii in range(ne):
		print "++ Writing Kappa-filtered TE#%i timeseries" % (ii+1)
		write_split_ts(catd[:,:,:,ii,:],comptable,mmix,'e%i' % (ii+1))

def ctabsel(ctabfile):
	ctlines = open(ctabfile).readlines()
	class_tags = ['#ACC','#REJ','#MID','#IGN']
	class_dict = {}
	for ii,ll in enumerate(ctlines):
		for kk in class_tags:
			if ll[:4]==kk and ll[4:].strip() != '':
				class_dict[kk] = ll[4:].split('#')[0].split(',')
	return tuple([np.array(class_dict[kk],dtype=int) for kk in class_tags])


###################################################################################################
# 						Begin Main
###################################################################################################

if __name__=='__main__':

	parser=OptionParser()
	parser.add_option('-d',"--orig_data",dest='data',help="Spatially Concatenated Multi-Echo Dataset",default=None)
	parser.add_option('-e',"--TEs",dest='tes',help="Echo times (in ms) ex: 15,39,63",default=None)
	parser.add_option('',"--mix",dest='mixm',help="Mixing matrix. If not provided, ME-PCA & ME-ICA (MDP) is done.",default=None)
	parser.add_option('',"--ctab",dest='ctab',help="Component table extract pre-computed classifications from.",default=None)
	parser.add_option('',"--manacc",dest='manacc',help="Comma separated list of manually accepted components",default=None)
	parser.add_option('',"--strict",dest='strict',action='store_true',help="Ignore low-variance ambiguous components",default=False)
	#parser.add_option('',"--wav",dest='wav',help="Perform wavelet PCA, default False",default=False)
	parser.add_option('',"--no_gscontrol",dest='no_gscontrol',action='store_true',help="Control global signal using spatial approach",default=False)
	parser.add_option('',"--kdaw",dest='kdaw',help="Dimensionality augmentation weight (Kappa). Default 10. -1 for low-dimensional ICA",default=10.)
	parser.add_option('',"--rdaw",dest='rdaw',help="Dimensionality augmentation weight (Rho). Default 1. -1 for low-dimensional ICA",default=1.)
	parser.add_option('',"--conv",dest='conv',help="Convergence limit. Default 2.5e-5",default='2.5e-5')
	parser.add_option('',"--sourceTEs",dest='ste',help="Source TEs for models. ex: -ste 2,3 ; -ste 0 for all, -1 for opt. com. Default -1.",default=0-1)	
	parser.add_option('',"--combmode",dest='combmode',help="Combination scheme for TEs: t2s (Posse 1999, default),ste(Poser)",default='t2s')	
	parser.add_option('',"--denoiseTEs",dest='dne',action='store_true',help="Denoise each TE dataset separately",default=False)	
	parser.add_option('',"--initcost",dest='initcost',help="Initial cost func. for ICA: pow3,tanh(default),gaus,skew",default='tanh')
	parser.add_option('',"--finalcost",dest='finalcost',help="Final cost func, same opts. as initial",default='tanh')	
	parser.add_option('',"--stabilize",dest='stabilize',action='store_true',help="Stabilize convergence by reducing dimensionality, for low quality data",default=False)
	parser.add_option('',"--fout",dest='fout',help="Output TE-dependence Kappa/Rho SPMs",action="store_true",default=False)
	parser.add_option('',"--filecsdata",dest='filecsdata',help="Save component selection data",action="store_true",default=False)
	parser.add_option('',"--label",dest='label',help="Label for output directory.",default=None)

	(options,args) = parser.parse_args()
	args = list(set(args) - set(['DEBUG']))

	print "-- ME-PCA/ME-ICA Component for ME-ICA %s--" % __version__

	if options.tes==None or options.data==None: 
		print "*+ Need at least data and TEs, use -h for help."		
		sys.exit()

	print "++ Loading Data"
	tes = np.fromstring(options.tes,sep=',',dtype=np.float32)
	ne = tes.shape[0]
	catim  = nib.load(options.data)	
	head   = catim.get_header()
	head.extensions = []
	head.set_sform(head.get_sform(),code=1)
	aff = catim.get_affine()
	catd = cat2echos(catim.get_data(),ne)
	nx,ny,nz,Ne,nt = catd.shape
	mu  = catd.mean(axis=-1)
	sig  = catd.std(axis=-1)

	"""Parse options, prepare output directory"""
	if options.fout: options.fout = aff
	else: options.fout=None
	kdaw = float(options.kdaw)
	rdaw = float(options.rdaw)
	if options.label!=None: dirname='%s' % '.'.join(['TED',options.label])
	else: dirname='TED'
	os.system('mkdir %s' % dirname)
	if options.mixm!=None: 
		try: os.system('cp %s %s/meica_mix.1D; cp %s %s/%s' % (options.mixm,dirname,options.mixm,dirname,os.path.basename(options.mixm)))
		except: pass
	if options.ctab!=None: 
		try: os.system('cp %s %s/comp_table.txt; cp %s %s/%s' % (options.mixm,dirname,options.mixm,dirname,os.path.basename(options.mixm)))
		except: pass

	os.chdir(dirname)
	
	print "++ Computing Mask"
	mask,masksum = makeadmask(catd,min=False,getsum=True)

	print "++ Computing T2* map"
	t2s,s0,t2ss,s0s,t2sG,s0G = t2sadmap(catd,mask,tes)

	#Condition values
	cap_t2s = scoreatpercentile(t2s.flatten(),99.5)
	t2s[t2s>cap_t2s*10]=cap_t2s
	niwrite(s0,aff,'s0v.nii')
	niwrite(t2s,aff,'t2sv.nii')
	niwrite(t2ss,aff,'t2ss.nii')
	niwrite(s0s,aff,'s0vs.nii')
	niwrite(s0G,aff,'s0vG.nii')
	niwrite(t2sG,aff,'t2svG.nii')

	#Optimally combine data
	OCcatd = optcom(catd,t2s,tes,mask)

	if not options.no_gscontrol: 
		gscontrol_raw()

	if options.mixm == None:
		print "++ Doing ME-PCA and ME-ICA with MDP"
		import mdp
		if 'DEBUG' in sys.argv:
			import ipdb
			ipdb.set_trace()
		
		nc,dd = tedpca(options.ste)
		mmix_orig = tedica(dd,cost=options.initcost)
		np.savetxt('__meica_mix.1D',mmix_orig)
		seldict,comptable,betas,mmix = fitmodels_direct(catd,mmix_orig,mask,t2s,tes,options.fout,reindex=True)
		np.savetxt('meica_mix.1D',mmix)
		if 'GROUP0' in sys.argv:
			group0_flag = True
		else: group0_flag = False
		acc,rej,midk,empty = selcomps(seldict,knobargs=args,debug=debug_mode,group0_only=group0_flag,strict_mode=options.strict)
		del dd
	else:
		mmix_orig = np.loadtxt('meica_mix.1D')
		eim = eimask(np.float64(fmask(catd,mask)))==1
		eimum = np.array(np.squeeze(unmask(np.array(eim,dtype=np.int).prod(1),mask)),dtype=np.bool)
		seldict,comptable,betas,mmix = fitmodels_direct(catd,mmix_orig,mask,t2s,tes,options.fout)
		if options.ctab == None: acc,rej,midk,empty = selcomps(seldict,knobargs=args,debug=debug_mode,strict_mode=options.strict)
		else: acc,rej,midk,empty = ctabsel(ctabfile)

	if len(acc)==0:
		print "\n** WARNING! No BOLD components detected!!! Please check data and results!\n"
    
	writeresults()
	gscontrol_mmix()
	if options.dne: writeresults_echoes()    
