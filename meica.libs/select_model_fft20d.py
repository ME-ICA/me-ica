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
# PROCEDURE 2a: Model fitting and component selection routines
"""

import numpy as np
import scipy.stats as stats
import scipy.signal as SS
from numpy import random
from sklearn import svm
import scipy.optimize

def fitmodels_direct(catd,mmix,mask,t2s,tes,fout=None,reindex=False,mmixN=None,full_sel=True,debugout=False):
	"""
   	Usage:
   	
   	fitmodels_direct(fout)
	
   	Input:
   	fout is flag for output of per-component TE-dependence maps
   	t2s is a (nx,ny,nz) ndarray
   	tes is a 1d array
   	"""

   	#Compute opt. com. raw data
	tsoc = np.array(optcom(catd,t2s,tes,mask),dtype=float)[mask]
	tsoc_mean = tsoc.mean(axis=-1)
	tsoc_dm = tsoc-tsoc_mean[:,np.newaxis]
	
	#Compute un-normalized weight dataset (features)
	if mmixN == None: mmixN=mmix
	#WTS = computefeats2(unmask(unmask(tsoc,mask)[t2s!=0],t2s!=0),mmixN,t2s!=0,normalize=False)
	WTS = computefeats2(unmask(tsoc,mask),mmixN,mask,normalize=False)

	#Compute PSC dataset - shouldn't have to refit data
	tsoc_B = get_coeffs(unmask(tsoc_dm,mask),mask,mmix)[mask]
	tsoc_Babs = np.abs(tsoc_B)
	PSC = tsoc_B/tsoc.mean(axis=-1)[:,np.newaxis]*100

	#Compute skews to determine signs based on unnormalized weights, correct mmix & WTS signs based on spatial distribution tails
	from scipy.stats import skew
	signs = skew(WTS,axis=0)
	signs /= np.abs(signs)
	mmix = mmix.copy()
	mmix*=signs
	WTS*=signs
	PSC*=signs
	totvar = (tsoc_B**2).sum()
	totvar_norm = (WTS**2).sum()

	#Compute Betas and means over TEs for TE-dependence analysis
	Ne = tes.shape[0]
	betas = cat2echos(get_coeffs(uncat2echos(catd,Ne),np.tile(mask,(1,1,Ne)),mmix),Ne)
	nx,ny,nz,Ne,nc = betas.shape
	Nm = mask.sum()
	NmD = (t2s!=0).sum()
	mu = catd.mean(axis=-1)
	tes = np.reshape(tes,(Ne,1))
	fmin,fmid,fmax = getfbounds(ne)

	#Mask arrays
	mumask   = fmask(mu,t2s!=0)
	#t2smask  = fmask(t2s,mask)
	t2smask  = fmask(t2s,t2s!=0)
	betamask = fmask(betas,t2s!=0)

	if debugout: fout=aff

	#Setup Xmats
	#Model 1
	X1 = mumask.transpose()
	
	#Model 2
	X2 = np.tile(tes,(1,NmD))*mumask.transpose()/t2smask.transpose()
	
	#Tables for component selection
	Kappas = np.zeros([nc])
	Rhos = np.zeros([nc])
	varex = np.zeros([nc])
	varex_norm = np.zeros([nc])
	Z_maps = np.zeros([Nm,nc])
	F_R2_maps = np.zeros([NmD,nc])
	F_S0_maps = np.zeros([NmD,nc])
	Z_clmaps = np.zeros([Nm,nc])
	F_R2_clmaps = np.zeros([NmD,nc])
	F_S0_clmaps = np.zeros([NmD,nc])
	Br_clmaps_R2 = np.zeros([Nm,nc])
	Br_clmaps_S0 = np.zeros([Nm,nc])

	for i in range(nc):

		#size of B is (nc, nx*ny*nz)
		B = np.atleast_3d(betamask)[:,:,i].transpose()
		alpha = (np.abs(B)**2).sum(axis=0)
		varex[i] = (tsoc_B[:,i]**2).sum()/totvar*100.
		varex_norm[i] = (unmask(WTS,mask)[t2s!=0][:,i]**2).sum()/totvar_norm*100.

		#S0 Model
		coeffs_S0 = (B*X1).sum(axis=0)/(X1**2).sum(axis=0)
		SSE_S0 = (B - X1*np.tile(coeffs_S0,(Ne,1)))**2
		SSE_S0 = SSE_S0.sum(axis=0)
		F_S0 = (alpha - SSE_S0)*2/(SSE_S0)
		F_S0_maps[:,i] = F_S0

		#R2 Model
		coeffs_R2 = (B*X2).sum(axis=0)/(X2**2).sum(axis=0)
		SSE_R2 = (B - X2*np.tile(coeffs_R2,(Ne,1)))**2
		SSE_R2 = SSE_R2.sum(axis=0)
		F_R2 = (alpha - SSE_R2)*2/(SSE_R2)
		F_R2_maps[:,i] = F_R2

		#Compute weights as Z-values
		wtsZ=(WTS[:,i]-WTS[:,i].mean())/WTS[:,i].std()
		wtsZ[np.abs(wtsZ)>Z_MAX]=(Z_MAX*(np.abs(wtsZ)/wtsZ))[np.abs(wtsZ)>Z_MAX]
		Z_maps[:,i] = wtsZ 

		#Compute Kappa and Rho
		F_S0[F_S0>F_MAX] = F_MAX
		F_R2[F_R2>F_MAX] = F_MAX
		Kappas[i] = np.average(F_R2,weights=np.abs(np.squeeze(unmask(wtsZ,mask)[t2s!=0]**2.)))
		Rhos[i] = np.average(F_S0,weights=np.abs(np.squeeze(unmask(wtsZ,mask)[t2s!=0]**2.)))

	#Tabulate component values
	comptab_pre = np.vstack([np.arange(nc),Kappas,Rhos,varex,varex_norm]).T
	if reindex:
		#Re-index all components in Kappa order
		comptab = comptab_pre[comptab_pre[:,1].argsort()[::-1],:]
		Kappas = comptab[:,1]; Rhos = comptab[:,2]; varex = comptab[:,3]; varex_norm = comptab[:,4]
		nnc = np.array(comptab[:,0],dtype=np.int)
		mmix_new = mmix[:,nnc]
		F_S0_maps = F_S0_maps[:,nnc]; F_R2_maps = F_R2_maps[:,nnc]; Z_maps = Z_maps[:,nnc]
		WTS = WTS[:,nnc]; PSC=PSC[:,nnc]; tsoc_B=tsoc_B[:,nnc]; tsoc_Babs=tsoc_Babs[:,nnc]
		comptab[:,0] = np.arange(comptab.shape[0])
	else:
		comptab = comptab_pre
		mmix_new = mmix

	#Full selection including clustering criteria
	seldict=None
	if full_sel: 
		for i in range(nc):	

			#Save out files
			out = np.zeros((nx,ny,nz,4))
			if fout!=None: 
				ccname = "cc%.3d.nii" % i
			else: ccname = ".cc_temp.nii.gz"

			out[:,:,:,0] = np.squeeze(unmask(PSC[:,i],mask))
			out[:,:,:,1] = np.squeeze(unmask(F_R2_maps[:,i],t2s!=0))
			out[:,:,:,2] = np.squeeze(unmask(F_S0_maps[:,i],t2s!=0))
			out[:,:,:,3] = np.squeeze(unmask(Z_maps[:,i],mask))

			#import ipdb; ipdb.set_trace()

			niwrite(out,fout,ccname)
			os.system('3drefit -sublabel 0 PSC -sublabel 1 F_R2  -sublabel 2 F_SO -sublabel 3 Z_sn %s 2> /dev/null > /dev/null'%ccname)

			csize = np.max([int(Nm*0.0005)+5,20])
			#csize = 10

			#Do simple clustering on F
			os.system("3dcalc -overwrite -a %s[1..2] -expr 'a*step(a-%i)' -prefix .fcl_in.nii.gz -overwrite" % (ccname,fmin))
			os.system('3dmerge -overwrite -dxyz=1 -1clust 1 %i -doall -prefix .fcl_out.nii.gz .fcl_in.nii.gz' % (csize))
			sel = fmask(nib.load('.fcl_out.nii.gz').get_data(),t2s!=0)!=0
			sel = np.array(sel,dtype=np.int)
			F_R2_clmaps[:,i] = sel[:,0]
			F_S0_clmaps[:,i] = sel[:,1]

			#Do simple clustering on Z at p<0.05
			sel = spatclust(None,mask,csize,1.95,head,aff,infile=ccname,dindex=3,tindex=3)
			Z_clmaps[:,i] = sel

			#Do simple clustering on ranked signal-change map
			countsigFR2 = F_R2_clmaps[:,i].sum()
			countsigFS0 = F_S0_clmaps[:,i].sum()
			Br_clmaps_R2[:,i] = spatclust(rankvec(tsoc_Babs[:,i]),mask,csize,max(tsoc_Babs.shape)-countsigFR2,head,aff)
			Br_clmaps_S0[:,i] = spatclust(rankvec(tsoc_Babs[:,i]),mask,csize,max(tsoc_Babs.shape)-countsigFS0,head,aff)

		seldict = {}
		selvars = ['Kappas','Rhos','WTS','varex','Z_maps','F_R2_maps','F_S0_maps',\
			'Z_clmaps','F_R2_clmaps','F_S0_clmaps','tsoc_B','Br_clmaps_R2','Br_clmaps_S0','PSC']
		for vv in selvars:
			seldict[vv] = eval(vv)
		
		if debugout or ('DEBUGOUT' in args):
			#Package for debug
			import cPickle as cP
			import zlib
			try: os.system('mkdir compsel.debug')
			except: pass
			selvars = ['Kappas','Rhos','WTS','varex','Z_maps','Z_clmaps','F_R2_clmaps','F_S0_clmaps','Br_clmaps_R2','Br_clmaps_S0','PSC']
			for vv in selvars:
				with open('compsel.debug/%s.pkl.gz' % vv,'wb') as ofh:
					print "Writing debug output: compsel.debug/%s.pkl.gz" % vv
					ofh.write(zlib.compress(cP.dumps(eval(vv))))
					ofh.close()

	return seldict,comptab,betas,mmix_new

def do_svm(train_set,train_labs,test_set,svmtype=0):
	if svmtype==2: probability=True
	else: probability = False
	clf = svm.SVC(kernel='linear',probability=probability)
	if svmtype==1: clf = svm.LinearSVC(loss='squared_hinge',penalty='l1',dual=False)
	clf.fit(train_set,train_labs)
	return clf.predict(test_set),clf	

def fft_variance(fproj_arr,fproj_arr_val,A,B):
	fproj_sel_T = stats.ttest_ind(fproj_arr[:,A].T,fproj_arr[:,B].T)
	fproj_sel_A = (andb([fproj_sel_T[0]>0,fproj_sel_T[1]<0.05])==2).reshape(mask.shape[0:2])
	fproj_sel_B = (andb([fproj_sel_T[0]<0,fproj_sel_T[1]<0.05])==2).reshape(mask.shape[0:2])
	return fproj_arr_val[fproj_sel_A.flatten()].sum(0),fproj_arr_val[fproj_sel_B.flatten()].sum(0)

def gaussian(height, center_x, center_y, width_x, width_y):
	"""Returns a gaussian function with the given parameters"""
	width_x = float(width_x)
	width_y = float(width_y)
	return lambda x,y: height*np.exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
	"""Returns (height, x, y, width_x, width_y)
	the gaussian parameters of a 2D distribution by calculating its
	moments """
	total = data.sum()
	X, Y = np.indices(data.shape)
	x = (X*data).sum()/total
	y = (Y*data).sum()/total
	col = data[:, int(y)]
	width_x = np.sqrt(abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
	row = data[int(x), :]
	width_y = np.sqrt(abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
	height = data.max()
	return height, x, y, width_x, width_y

def fitgaussian(data):
	"""Returns (height, x, y, width_x, width_y)
	the gaussian parameters of a 2D distribution found by a fit"""
	params = moments(data)
	errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) - data)
	p, success = scipy.optimize.leastsq(errorfunction, params)
	return p

def selcomps(seldict,debug=False,olevel=2,oversion=99,knobargs='',filecsdata=False,savecsdiag=True,group0_only=False,strict_mode=False):

	selmodelversion='fft20c.051517'

	#import ipdb
	import numpy.fft as fft
	from sklearn import svm
	from sklearn.cluster import DBSCAN

	try:
		if options.filecsdata: filecsdata=True
	except: 
		pass
	
	if filecsdata: 
		import cPickle as pickle
		import bz2
		if seldict!=None:
			print "Saving component selection data"
			csstate_f = bz2.BZ2File('compseldata.pklbz','wb')
			pickle.dump(seldict,csstate_f)
			csstate_f.close()
		else:
			try:
				csstate_f = bz2.BZ2File('compseldata.pklbz','rb')
				seldict = pickle.load(csstate_f)
				csstate_f.close()
			except:
				print "No component data found!"
				return None

	#Dump dictionary into variable names
	for key in seldict.keys(): exec("%s=seldict['%s']" % (key,key))

	#List of components
	midk = []
	ign = []
	nc = np.arange(len(Kappas))
	ncl = np.arange(len(Kappas))

	#If user has specified components to accept manually
	try:
		if options.manacc:
			acc = sorted([int(vv) for vv in options.manacc.split(',')])
			midk = []
			rej = sorted(np.setdiff1d(ncl,acc))
			return acc,rej,midk,[] #Add string for ign
	except: 
		pass

	"""
	Set knobs
	"""
	if knobargs!='': 
		for knobarg in ''.join(knobargs).split(','): exec(knobarg)

	"""
	Do some tallies for no. of significant voxels
	"""
	countsigZ = Z_clmaps.sum(0)
	countsigFS0 = F_S0_clmaps.sum(0)
	countsigFR2 = F_R2_clmaps.sum(0)
	countnoise = np.zeros(len(nc))

	"""
	Make table of dice values
	"""
	dice_table = np.zeros([nc.shape[0],2])
	csize = np.max([int(mask.sum()*0.0005)+5,20])
	for ii in ncl:
		dice_FR2 = dice(unmask(Br_clmaps_R2[:,ii],mask)[t2s!=0],F_R2_clmaps[:,ii])
		dice_FS0 = dice(unmask(Br_clmaps_S0[:,ii],mask)[t2s!=0],F_S0_clmaps[:,ii])
		dice_table[ii,:] = [dice_FR2,dice_FS0] #step 3a here and above
	dice_table[np.isnan(dice_table)]=0

	if debug:
		import pdb
		pdb.set_trace()
		#import IPython
		#from IPython.core.debugger import Tracer; Tracer()()

	"""
	Make table of noise gain
	"""
	tt_table = np.zeros([len(nc),4])
	counts_FR2_Z = np.zeros([len(nc),2])
	for ii in nc:
		comp_noise_sel = andb([np.abs(Z_maps[:,ii])>1.95,Z_clmaps[:,ii]==0])==2
		countnoise[ii] = np.array(comp_noise_sel,dtype=np.int).sum()
		noise_FR2_Z = np.log10(np.unique(F_R2_maps[unmask(comp_noise_sel,mask)[t2s!=0],ii]))
		signal_FR2_Z  = np.log10(np.unique(F_R2_maps[unmask(Z_clmaps[:,ii],mask)[t2s!=0]==1,ii]))
		counts_FR2_Z[ii,:] = [len(signal_FR2_Z),len(noise_FR2_Z)]
		try: 
			ttest = stats.ttest_ind(signal_FR2_Z,noise_FR2_Z,equal_var=True)
			mwu = stats.norm.ppf(stats.mannwhitneyu(signal_FR2_Z,noise_FR2_Z)[1])
			tt_table[ii,0] = np.abs(mwu)*ttest[0]/np.abs(ttest[0])
			tt_table[ii,1] = ttest[1]
		except: pass
	tt_table[np.isnan(tt_table)]=0
	#import pdb; pdb.set_trace()
	tt_table[np.isinf(tt_table[:,0]),0]=np.percentile(tt_table[~np.isinf(tt_table[:,0]),0],98) 

	#Time series derivative kurtosis
	mmix_dt = (mmix[:-1]-mmix[1:])
	mmix_dt2 = (mmix_dt[:-1]-mmix_dt[1:])
	mmix_kurt = stats.kurtosis(mmix_dt)

	#Polynomial detrend of mmix
	p0base = np.array([(np.arange(mmix.shape[0])-np.mean(np.arange(mmix.shape[0])))/np.std(np.arange(mmix.shape[0])),np.ones(mmix.shape[0])])
	mmixp0 = mmix-np.dot(np.linalg.lstsq(mmix,p0base.T)[0],p0base).T
	mmixp0_dt = (mmixp0[:-1]-mmixp0[1:])
	mmixp0_dt2 = (mmixp0_dt[:-1]-mmixp0_dt[1:])
	mmixp0_kurt = stats.kurtosis(mmixp0_dt)

	if debug:
		import ipdb
		ipdb.set_trace()

	"""
	Step 1: Reject anything that's obviously an artifact
	a. Estimate a null variance
	"""
	rej = ncl[andb([Rhos>Kappas,countsigFS0>countsigFR2])>0]
	ncl = np.setdiff1d(ncl,rej)

	if debug:
		import ipdb
		ipdb.set_trace()

	"""
	Step 2: Compute 3-D spatial FFT of Beta maps to detect high-spatial frequency artifacts
	"""
	fproj_arr = np.zeros([np.prod(mask.shape[0:2]),len(nc)])
	fproj_arr_val = np.zeros([np.prod(mask.shape[0:2]),len(nc)])
	spr = []
	fdist = []
	for ii in nc: 
		fproj = np.fft.fftshift(np.abs(np.fft.rfftn(unmask(seldict['PSC'],mask)[:,:,:,ii])))
		fproj_z = fproj.max(2)
		fproj[fproj==fproj.max()] = 0
		fproj_arr[:,ii] = rankvec(fproj_z.flatten())
		fproj_arr_val[:,ii] = fproj_z.flatten()
		spr.append(np.array(fproj_z>fproj_z.max()/4,dtype=np.int).sum())
		fprojr = np.array([fproj,fproj[:,:,::-1]]).max(0)
		fdist.append(np.max([ fitgaussian(fproj.max(jj))[3:].max() for jj in range(len(fprojr.shape)) ]))
	fdist = np.array(fdist)
	spr = np.array(spr)

	if debug:
		import ipdb
		ipdb.set_trace()

	"""
	Step 3: Create feature space of component properties
	"""	
	fdist_pre = fdist.copy()
	fdist_pre[fdist>np.median(fdist)*3] = np.median(fdist)*3
	fdist_z = (fdist_pre - np.median(fdist_pre) ) / fdist_pre.std()
	spz = (spr-spr.mean())/spr.std()
	Tz = (tt_table[:,0]-tt_table[:,0].mean())/tt_table[:,0].std()
	varex_ = np.log(varex)
	Vz = (varex_-varex_.mean())/varex_.std()
	Kz = (Kappas-Kappas.mean())/Kappas.std()
	Rz = (Rhos-Rhos.mean())/Rhos.std()
	Ktz = np.log(Kappas)/2
	Ktz = (Ktz-Ktz.mean())/Ktz.std()
	Rtz = np.log(Rhos)/2
	Rtz = (Rtz-Rtz.mean())/Rtz.std()
	KRr = stats.zscore(np.log(Kappas)/np.log(Rhos))
	cnz = (countnoise-countnoise.mean())/countnoise.std()
	Dz = stats.zscore(np.arctanh(dice_table[:,0]+0.001))
	fz = np.array([Tz,Vz,Ktz,KRr,cnz,Rz,mmix_kurt,fdist_z])

	if debug:
		import ipdb
		ipdb.set_trace()

	"""
	Step 3: Make initial guess of where BOLD components are and use DBSCAN to exclude noise components and find a sample set of 'good' components
	"""
	#epsmap is [index,level of overlap with dicemask,number of high Rho components]
	F05,F025,F01 = getfbounds(ne)
	epsmap = []
	Rhos_sorted = np.array(sorted(Rhos))[::-1]
	#Make an initial guess as to number of good components based on consensus of control points across Rhos and Kappas
	KRcutguesses = [getelbow(Rhos),getelbow2(Rhos),getelbow3(Rhos),getelbow(Kappas),getelbow2(Kappas),getelbow3(Kappas)]
	Kelbowval = np.median([getelbow(Kappas,True),getelbow2(Kappas,True),getelbow3(Kappas,True)]+list(getfbounds(ne)))
	Khighelbowval = stats.scoreatpercentile([getelbow(Kappas,True),getelbow2(Kappas,True),getelbow3(Kappas,True)]+list(getfbounds(ne)),75)
	KRcut = np.median(KRcutguesses)
	#only use exclusive when inclusive is extremely inclusive - double KRcut
	if getelbow2(Kappas) > KRcut*2 and getelbow(Kappas,True)<F01: Kcut = getelbow(Kappas,True)
	else: Kcut = getelbow2(Kappas,True)
	#only use inclusive when exclusive is extremely exclusive - half KRcut (remember for Rho inclusive is higher, so want both Kappa and Rho to defaut to lower) 
	if getelbow2(Rhos) > KRcut*2 : Rcut = getelbow(Rhos,True) #consider something like min([getelbow(Rhos,True),sorted(Rhos)[::-1][KRguess] ])
	else: Rcut = getelbow2(Rhos,True)
	if Rcut > Kcut: Kcut = Rcut #Rcut should never be higher than Kcut
	KRelbow = andb([Kappas>Kcut,Rhos<Rcut ] )
	#Make guess of Kundu et al 2011 plus remove high frequencies, generally high variance, and high variance given low Kappa 
	tt_lim = scoreatpercentile(tt_table[tt_table[:,0]>0,0],75)/3
	KRguess = np.setdiff1d(np.setdiff1d(nc[KRelbow==2],rej),np.union1d(nc[tt_table[:,0]<tt_lim],np.union1d(np.union1d(nc[spz>1],nc[Vz>2]),nc[andb([varex>0.5*sorted(varex)[::-1][int(KRcut)],Kappas<2*Kcut])==2])))
	guessmask = np.zeros(len(nc))
	guessmask[KRguess] = 1

	#Throw lower-risk bad components out
	rejB = ncl[andb([tt_table[ncl,0]<0,varex[ncl]>np.median(varex),ncl > KRcut])==3]
	rej = np.union1d(rej,rejB)
	ncl = np.setdiff1d(ncl,rej)

	if debug:
		import ipdb
		ipdb.set_trace()

	for ii in range(20000):
		db = DBSCAN(eps=.005+ii*.005, min_samples=3).fit(fz.T)
		if db.labels_.max() > 1 and db.labels_.max() < len(nc)/6 and np.intersect1d(rej,nc[db.labels_==0]).shape[0]==0 and np.array(db.labels_==-1,dtype=int).sum()/float(len(nc))<.5:
			epsmap.append([ii, dice(guessmask,db.labels_==0),np.intersect1d(nc[db.labels_==0],nc[Rhos>getelbow(Rhos_sorted,True)]).shape[0]   ])
			if debug: print "found solution", ii, db.labels_
		db = None

	if debug:
		import pdb
		pdb.set_trace()

	epsmap = np.array(epsmap)
	group0 = []
	dbscanfailed=False
	if len(epsmap)!=0 :
		#Select index that maximizes Dice with guessmask but first minimizes number of higher Rho components
		ii = epsmap[np.argmax(epsmap[epsmap[:,2]==np.min(epsmap[:,2]),1],0),0]
		print 'Component selection tuning: ' , epsmap[:,1].max()
		db = DBSCAN(eps=.005+ii*.005, min_samples=3).fit(fz.T)
		ncl = nc[db.labels_==0]
		ncl = np.setdiff1d(ncl,rej)
		ncl = np.setdiff1d(ncl,ncl[ncl>len(nc)-len(rej)])
		group0 = ncl.copy()
		group_n1 = nc[db.labels_==-1]
		to_clf = np.setdiff1d(nc,np.union1d(ncl,rej))
	if len(group0)==0 or len(group0) < len(KRguess)*.5:
		dbscanfailed=True
		print "DBSCAN bassed guess failed. Using elbow guess method."
		ncl = np.setdiff1d(np.setdiff1d(nc[KRelbow==2],rej),np.union1d(nc[tt_table[:,0]<tt_lim],np.union1d(np.union1d(nc[spz>1],nc[Vz>2]),nc[andb([varex>0.5*sorted(varex)[::-1][int(KRcut)],Kappas<2*Kcut])==2])))
		group0 = ncl.copy()
		group_n1 = []
		to_clf = np.setdiff1d(nc,np.union1d(group0,rej))
	if len(group0)<2 or (len(group0)<4 and float(len(rej))/len(group0)>3):
		print "WARNING: Extremely limited reliable BOLD signal space. Not filtering further into midk etc."
		midkfailed = True
		min_acc = np.array([])
		if len(group0)!=0:		
			toacc_hi = np.setdiff1d(nc [andb([ fdist <= np.max(fdist[group0]), Rhos<F025, Vz>-2 ])==3  ],np.union1d(group0,rej)) #For extremes, building in a 20% tolerance
			min_acc = np.union1d(group0,toacc_hi)
			to_clf = np.setdiff1d(nc , np.union1d(min_acc,rej) )		
		diagstepkeys=['rej','KRcut','Kcut','Rcut','dbscanfailed','midkfailed','KRguess','group0','min_acc','toacc_hi']
		diagstepout=[]
		for ddk in diagstepkeys: diagstepout.append("%s: %s" %  (ddk,eval('str(%s)' % ddk) ) )
		with open('csstepdata.txt','w') as ofh:
			ofh.write('\n'.join(diagstepout))
		ofh.close()
		return list(sorted(min_acc)),list(sorted(rej)),[],list(sorted(to_clf))

	if group0_only: return list(sorted(group0)),list(sorted(rej)),[],list(sorted(to_clf))

	if debug:
		import ipdb
		ipdb.set_trace()

	#Find additional components to reject based on Dice - doing this here since Dice is a little unstable, need to reference group0
	rej_supp = []
	dice_rej = False
	if not dbscanfailed and len(rej)+len(group0)<0.75*len(nc):
		dice_rej = True
		rej_supp = np.setdiff1d(np.setdiff1d(np.union1d(rej,nc[dice_table[nc,0]<=dice_table[nc,1]]  ),group0),group_n1)
		rej = np.union1d(rej,rej_supp)

	if debug:
		import ipdb
		ipdb.set_trace()


	#Temporal features
	mmix_kurt_z = (mmix_kurt-mmix_kurt[group0].mean())/mmix_kurt[group0].std()
	mmixp0_kurt_z = (mmixp0_kurt-mmixp0_kurt[group0].mean())/mmixp0_kurt[group0].std()
	mmix_kurt_z_max  = np.max([mmix_kurt_z,mmixp0_kurt_z],0)
	
	"""
	Step 2: Classifiy midk and ignore using separte SVMs for difference variance regimes
	#To render hyperplane:
	min_x = np.min(spz2);max_x=np.max(spz2)
   	# plotting separating hyperplane
    	ww = clf_.coef_[0]
    	aa = -ww[0] / ww[1]
    	xx = np.linspace(min_x - 2, max_x + 2)  # make sure the line is long enough
    	yy = aa * xx - (clf_.intercept_[0]) / ww[1]
    	plt.plot(xx, yy, '-')	
	"""
	if debug:
		import pdb
		pdb.set_trace()

	toacc_hi = np.setdiff1d(nc [andb([ fdist <= np.max(fdist[group0]), Rhos<F025, Vz>-2 ])==3  ],np.union1d(group0,rej))  #Tried getting rid of accepting based on SVM altogether, now using only rejecting
	toacc_lo = np.intersect1d(to_clf,nc[andb([spz<1,Rz<0,mmix_kurt_z_max<5,Dz>-1,Tz>-1,Vz<0,Kappas>=F025,fdist<3*np.percentile(fdist[group0],98)])==8])		
	midk_clf,clf_ = do_svm(fproj_arr_val[:,np.union1d(group0,rej)].T,[0]*len(group0) + [1]*len(rej),fproj_arr_val[:,to_clf].T,svmtype=2)
	midk = np.setdiff1d(to_clf[andb([midk_clf==1,varex[to_clf]>np.median(varex[group0])  ])==2],np.union1d(toacc_hi,toacc_lo))
	if len(np.intersect1d(to_clf[andb([midk_clf==1,Vz[to_clf]>0 ])==2],toacc_hi))==0:
		svm_acc_fail = True
		toacc_hi = np.union1d(toacc_hi,to_clf[midk_clf==0])  #only use SVM to augment toacc_hi only if toacc_hi isn't already conflicting with SVM choice
	else: svm_acc_fail = False
	
	"""
	Step 3: Compute variance associated with low T2* areas (e.g. draining veins and low T2* areas) 
	#To write out veinmask
	veinout = np.zeros(t2s.shape)
	veinout[t2s!=0] = veinmaskf
	niwrite(veinout,aff,'veinmaskf.nii',header=head)
	veinBout = unmask(veinmaskB,mask)
	niwrite(veinBout,aff,'veins50.nii',header=head)
	"""
	tsoc_B_Zcl = np.zeros(tsoc_B.shape)
	tsoc_B_Zcl[Z_clmaps!=0] = np.abs(tsoc_B)[Z_clmaps!=0]
 	sig_B = [ stats.scoreatpercentile(tsoc_B_Zcl[tsoc_B_Zcl[:,ii]!=0,ii],25) if len(tsoc_B_Zcl[tsoc_B_Zcl[:,ii]!=0,ii]) != 0  else 0 for ii in nc   ] 
 	sig_B = np.abs(tsoc_B)>np.tile(sig_B,[tsoc_B.shape[0],1])
	veinmask = andb([t2s<scoreatpercentile(t2s[t2s!=0],15),t2s!=0])==2
	veinmaskf = veinmask[t2s!=0]
	veinR = np.array(sig_B[veinmaskf].sum(0),dtype=float)/sig_B[~veinmaskf].sum(0)
	veinR[np.isnan(veinR)] = 0
	veinc = np.union1d(rej,midk)
	rej_veinRZ = ((veinR-veinR[veinc].mean())/veinR[veinc].std())[veinc]
	rej_veinRZ[rej_veinRZ<0] = 0
	rej_veinRZ[ countsigFR2[veinc] > np.array(veinmaskf,dtype=int).sum()] =0 
	t2s_lim = [stats.scoreatpercentile(t2s[t2s!=0],50),stats.scoreatpercentile(t2s[t2s!=0],80)/2]
	phys_var_zs = []
	for t2sl_i in range(len(t2s_lim)):
		t2sl = t2s_lim[t2sl_i]
		veinW = sig_B[:,veinc]*np.tile(rej_veinRZ,[sig_B.shape[0],1])
		veincand = fmask(unmask(andb([s0[t2s!=0]<np.median(s0[t2s!=0]),t2s[t2s!=0]<t2sl])>=1,t2s!=0),mask)
		veinW[~veincand]=0
		invein = veinW.sum(1)[fmask(unmask(veinmaskf,t2s!=0)*unmask(veinW.sum(1)>1,mask),mask)]
		minW = 10*(np.log10(invein).mean())-1*10**(np.log10(invein).std())
		veinmaskB = veinW.sum(1)>minW
		tsoc_Bp = tsoc_B.copy()
		tsoc_Bp[tsoc_Bp<0]=0
		sig_Bp = sig_B*tsoc_Bp>0
		vvex = np.array([(tsoc_Bp[veinmaskB,ii]**2.).sum()/(tsoc_Bp[:,ii]**2.).sum() for ii in nc])
		group0_res = np.intersect1d(KRguess,group0)
		phys_var_zs.append( (vvex-vvex[group0_res].mean())/vvex[group0_res].std() )
		veinBout = unmask(veinmaskB,mask)
		niwrite(veinBout,aff,'veins_l%i.nii' % t2sl_i,header=head)
	#Mask to sample veins
	phys_var_z = np.array(phys_var_zs).max(0)
	Vz2 = (varex_ - varex_[group0].mean())/varex_[group0].std()

	"""
	Step 4: Learn joint TE-dependence spatial and temporal models to move remaining artifacts to ignore class
	"""
	if debug:
		import ipdb
		ipdb.set_trace()

	to_ign = []
	
	minK_ign = np.max([F05,getelbow2(Kappas,True)])
	newcest = len(group0)+len(toacc_hi[ Kappas[toacc_hi]>minK_ign ])
	phys_art = np.setdiff1d(nc[andb([phys_var_z>3.5,Kappas<minK_ign])==2],group0)
	phys_art = np.union1d(np.setdiff1d(nc[andb([phys_var_z>2,rankvec(phys_var_z)-rankvec(Kappas)>newcest/2,Vz2>-1])==3],group0),phys_art)
	#Want to replace field_art with an acf/SVM based approach instead of a kurtosis/filter one
	field_art = np.setdiff1d(nc[andb([mmix_kurt_z_max>5,Kappas<minK_ign])==2],group0)
	field_art = np.union1d(np.setdiff1d(nc[andb([mmix_kurt_z_max>2,rankvec(mmix_kurt_z_max)-rankvec(Kappas)>newcest/2,Vz2>1,Kappas<F01])==4],group0),field_art)
	misc_art = np.setdiff1d(nc[andb([(rankvec(Vz)-rankvec(Ktz))>newcest/2,Kappas<Khighelbowval])==2],group0)
	ign_cand = np.unique(list(field_art)+list(phys_art)+list(misc_art))
	g0_red = np.setdiff1d(group0,ign_cand)
	midkrej = np.union1d(midk,rej)
	to_ign = np.setdiff1d(list(ign_cand),midkrej)
	toacc = np.union1d(toacc_hi,toacc_lo)
	ncl = np.setdiff1d(np.union1d(ncl,toacc),np.union1d(to_ign,midkrej))
	ign = np.setdiff1d(nc,list(ncl)+list(midk)+list(rej))
	orphan = np.setdiff1d(nc,list(ncl)+list(to_ign)+list(midk)+list(rej))

	#Last ditch effort to save some transient components
	if not strict_mode:
		Vz3 = (varex_ - varex_[ncl].mean())/varex_[ncl].std()
		ncl = np.union1d(ncl,np.intersect1d(orphan,nc[andb([Kappas>F05,Rhos<F025,Kappas>Rhos,Vz3<=-1,Vz3>-3,mmix_kurt_z_max<2.5])==6]))
		ign = np.setdiff1d(nc,list(ncl)+list(midk)+list(rej))
		orphan = np.setdiff1d(nc,list(ncl)+list(to_ign)+list(midk)+list(rej))

	if debug:
		import pdb
		pdb.set_trace()

	if savecsdiag:
		diagstepkeys=['selmodelversion','rej','KRcut','Kcut','Rcut','dbscanfailed','KRguess','group0','dice_rej','rej_supp','to_clf','midk', 'svm_acc_fail', 'toacc_hi','toacc_lo','field_art','phys_art','misc_art','ncl','ign']
		diagstepout=[]
		for ddk in diagstepkeys: diagstepout.append("%s: %s" %  (ddk,eval('str(%s)' % ddk) ) )
		with open('csstepdata.txt','w') as ofh:
			ofh.write('\n'.join(diagstepout))
		allfz = np.array([Tz,Vz,Ktz,KRr,cnz,Rz,mmix_kurt,fdist_z])
		np.savetxt('csdata.txt',allfz)

	return list(sorted(ncl)),list(sorted(rej)),list(sorted(midk)),list(sorted(ign))
