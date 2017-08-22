#!/usr/bin/env python
__version__="v3.1 beta1"
welcome_block="""
# Multi-Echo ICA, Version %s
# See http://dx.doi.org/10.1016/j.neuroimage.2011.12.028
# Kundu, P., Inati, S.J., Evans, J.W., Luh, W.M. & Bandettini, P.A. Differentiating 
#	BOLD and non-BOLD signals in fMRI time series using multi-echo EPI. NeuroImage (2011).
#
# Kundu, P., Inati, S.J., Evans, J.W., Luh, W.M. & Bandettini, P.A. Differentiating 
#   BOLD and non-BOLD signals in fMRI time series using multi-echo EPI. NeuroImage (2011).
# http://dx.doi.org/10.1016/j.neuroimage.2011.12.028
#
# t2smap.py version %s 	(c) 2014 Prantik Kundu, Noah Brenowitz, Souheil Inati
#
#Computes T2* map
""" % (__version__,__version__)

import os
from optparse import OptionParser
import numpy as np
import nibabel as nib
from sys import stdout,argv

import scipy.stats as stats

if 'DEBUG' in argv: 
	import ipdb
	debug_mode = True
else: debug_mode = False

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
            score = values[np.floor(idx)]
        elif interpolation_method == 'higher':
            score = values[np.ceil(idx)]
        else:
            raise ValueError("interpolation_method can only be 'fraction', " \
                             "'lower' or 'higher'")
    return score

def niwrite(data,affine, name , header=None):
	data[np.isnan(data)]=0
	stdout.write(" + Writing file: %s ...." % name) 
	
	thishead = header
	if thishead == None:
		thishead = head.copy()
		thishead.set_data_shape(list(data.shape))

	outni = nib.Nifti1Image(data,affine,header=thishead)
	outni.set_data_dtype('float64')
	outni.to_filename(name)
	

	print 'done.'

	return outni

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

def makemask(cdat,min=True,getsum=False):

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

	return np.reshape(out,(nx,ny,nz,nt))

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
	B = np.reshape(np.abs(echodata)+1, (Nm,Ne*nt)).transpose()
	B = np.log(B)
	x = np.array([np.ones(Ne),-tes])
	X = np.tile(x,(1,nt))
	X = np.sort(X)[:,::-1].transpose()

	beta,res,rank,sing = np.linalg.lstsq(X,B)
	t2s = 1/beta[1,:].transpose()
	s0  = np.exp(beta[0,:]).transpose()

	out = unmask(t2s,mask),unmask(s0,mask)
	out[0][np.isnan(out[0])]=0.

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

	for ne in range(2,Ne+1):
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

	fl = np.zeros([nx,ny,nz,len(tes)-2+1])
	for ne in range(Ne-1):
		fl_ = np.squeeze(fl[:,:,:,ne])
		fl_[masksum==ne+2] = True
		fl[:,:,:,ne] = fl_
	fl = np.array(fl,dtype=bool)

	if debug_mode : ipdb.set_trace()

	t2sa = unmask(t2ss[fl],masksum>1)
	s0va = unmask(s0vs[fl],masksum>1)
	
	return t2sa,s0va,t2ss,s0vs

def optcom(data,t2s,tes,mask):
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

###################################################################################################
# 						Begin Main
###################################################################################################

if __name__=='__main__':

	parser=OptionParser()
	parser.add_option('-d',"--orig_data",dest='data',help="Spatially Concatenated Multi-Echo Dataset",default=None)
	parser.add_option('-c',"--combmode",dest='combmode',help="Combination scheme for TEs: t2s (Posse 1999),ste(Poser,2006 default)",default='ste')	
	parser.add_option('-l',"--label",dest='label',help="Optional label to tag output files with",default=None)
	parser.add_option('-e',"--TEs",dest='tes',help="Echo times (in ms) ex: 15,39,63",default=None)

	(options,args) = parser.parse_args()

	print "-- T2* Map Component for ME-ICA %s --" % (__version__)

	if options.tes==None or options.data==None: 
		print "*+ Need at least data and TEs, use -h for help."		
		sys.exit()

	if options.label!=None:
		suf='_%s' % str(options.label)
	else:
		suf=''

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
	
	print "++ Computing Mask"
	mask,masksum  = makemask(catd,min=False,getsum=True)

	print "++ Computing Adaptive T2* map"
	t2s,s0,t2ss,s0vs   = t2sadmap(catd,mask,tes)
	niwrite(masksum,aff,'masksum%s.nii' % suf )
	niwrite(t2ss,aff,'t2ss%s.nii' % suf )
	niwrite(s0vs,aff,'s0vs%s.nii' % suf )

	print "++ Computing optimal combination"
	tsoc = np.array(optcom(catd,t2s,tes,mask),dtype=float)

	#Clean up numerical errors
	t2sm = t2s.copy()
	tsoc[np.isnan(tsoc)]=0
	s0[np.isnan(s0)]=0
	s0[s0<0]=0
	t2s[np.isnan(t2s)]=0
	t2s[t2s<0]=0
	t2sm[np.isnan(t2sm)]=0
	t2sm[t2sm<0]=0

	niwrite(tsoc,aff,'ocv%s.nii' % suf)
	niwrite(s0,aff,'s0v%s.nii' % suf)
	niwrite(t2s,aff,'t2sv%s.nii' % suf )
	niwrite(t2sm,aff,'t2svm%s.nii' % suf )
	
	

	


	




