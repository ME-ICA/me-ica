import os, re, sys
import commands
from gooey import Gooey, GooeyParser
from os import system,getcwd,mkdir,chdir,popen
from string import rstrip,split

startdir=os.path.abspath(os.path.curdir)
meicadir=os.path.dirname(os.path.abspath(os.path.expanduser(sys.argv[0])))


#Filename parser for NIFTI and AFNI files
def dsprefix(idn):
	def prefix(datasetname):
		return split(datasetname,'+')[0]
	if len(split(idn,'.'))!=0:
		if split(idn,'.')[-1]=='HEAD' or split(idn,'.')[-1]=='BRIK' or split(idn,'.')[-2:]==['BRIK','gz']:
			return prefix(idn)
		elif split(idn,'.')[-1]=='nii' and not split(idn,'.')[-1]=='nii.gz':
			return '.'.join(split(idn,'.')[:-1])
		elif split(idn,'.')[-2:]==['nii','gz']:
			return '.'.join(split(idn,'.')[:-2])
		else:
			return prefix(idn)
	else:
		return prefix(idn)

def dssuffix(idna):
	suffix = idna.split(dsprefix(idna))[-1]
	#print suffix
	spl_suffix=suffix.split('.')
	#print spl_suffix
	if len(spl_suffix[0])!=0 and spl_suffix[0][0] == '+': return spl_suffix[0]
	else: return suffix

def version_checker(cur_ver,ref_ver):
	"""
	Checks version in major/minor format of a current version versus a reference version.
	Supports 2 or more level versioning, only 2 levels checked)

	Input:
	cur_ver: float or string of version to check
	ref_ver: float or string of reference version

	Returns:
	bool for pass or fail
	"""
	cur_ver = str(cur_ver)
	ref_ver = str(ref_ver)
	cur_Varr = [int(vvv) for vvv in cur_ver.split('.')[0:2]]
	cur_V = cur_Varr[0]*10**2 + cur_Varr[1]%100
	ref_Varr = [int(vvv) for vvv in ref_ver.split('.')[0:2]]
	ref_V = ref_Varr[0]*10**2 + ref_Varr[1]%100
	if cur_V > ref_V: return True
	return False

#Run dependency check
def dep_check():
	print '++ Checking system for dependencies...'
	fails=0
	numpy_installed = 0
	scipy_installed = 0
	sklearn_installed = 0
	python_version_ok = 0
	global grayweight_ok
	grayweight_ok = 0
	print " + Python version: %s" % ('.'.join([str(v) for v in sys.version_info[0:3]]))
	if sys.version_info < (2,6) or sys.version_info > (3,0):
		print "*+ Python 2.x is <2.6, please upgrade to Python 2.x >= 2.6 & Numpy >= 1.5.x."
		fails+=1
	else:
		python_version_ok = 1
	try:
		import numpy
		numpy_installed = 1
	except:
		print "*+ Can't import Numpy! Please check Numpy installation for this Python environment."
		fails+=1
	try:
		import scipy
		scipy_installed = 1
	except:
		print "*+ Can't import Scipy! Please check Scipy installation for this Python environment."
		fails+=1
	try:
		import sklearn
		sklearn_installed = 1
	except:
		print "*+ Can't import scikit-learn! Please install scikit-learn version >0.15.0 for this Python environment."
		fails+=1

	if numpy_installed:
		print " + Numpy version: %s" % (numpy.__version__)
		if version_checker(numpy.__version__,1.5)==False:
			fails+=1
			print "*+ Numpy version is too old! Please upgrade to Numpy >=1.5.x!"
		import numpy.__config__ as nc
		if nc.blas_opt_info == {}:
			fails+=1
			print "*+ Numpy is not linked to BLAS! Please check Numpy installation."
	if scipy_installed:
		print " + Scipy version: %s" % (scipy.__version__)
		if version_checker(scipy.__version__,0.11)==False:
			fails+=1
			print "*+ Scipy version is too old! Please upgrade to Scipy >=0.11.x!"
	if sklearn_installed:
		print " + scikit-learn version: %s" % (sklearn.__version__)
		if version_checker(sklearn.__version__,0.15)==False:
			fails+=1
			print "*+ scikit-learn version is too old! Please upgrade to Scipy >=0.15.x!"
	afnicheck = commands.getstatusoutput("3dinfo")
	afnisegcheck = commands.getstatusoutput("3dSeg -help")
	if afnicheck[0]!=0:
		print "*+ Can't run AFNI binaries. Make sure AFNI is on the path!"
		fails+=1
	elif not afnicheck[1].__contains__('Alternate Alternative Usage') and len(afnisegcheck[1]) < 1000:
		print "*+ This seems like an old version of AFNI. Please upgrade to latest version of AFNI."
		fails+=1
	if afnisegcheck[0]==0 and len(afnisegcheck[1]) >= 1000:
		print " + Using AFNI 3dSeg for gray matter weighted anatomical-functional coregistration"
		grayweight_ok = 1
	if grayweight_ok==0:
		print "*+ WARNING: AFNI 3dSeg not available for gray matter weighted coregistration. See README."
	if fails==0:
		print " + Dependencies OK."
	else:
		print "*+ ABORTING. Please see error messages."
		return False
	return True

def getdsname(e_ii,prefixonly=False):
	if shorthand_dsin: dsname = '%s%s%s%s' % (prefix,datasets[e_ii],trailing,isf)
	else: dsname =  datasets_in[e_ii]
	if prefixonly: return dsprefix(dsname)
	else: return dsname

running = True

@Gooey(advanced=True,required_cols=1, optional_cols=1, program_name="Multi-Echo ICA v3.2", default_size=(800, 650),monospace_display=True)
def main():
	#Configure options and help dialog
	parser=GooeyParser(description='Enter dataset configuration for basic ME-ICA processing.')
	#opssubparsers = parser.add_subparsers()
	#mainopts = opssubparsers.add_parser("Main Processing Options")
	#ex: mainopts.add_argument('--echoes',"--echoes",dest='TEs',help="ex: -e 14.5,38.5,62.5  Echo times (in ms)",default='')
	#parser.add_argument('--subid',"--subid",dest='Subject ID',help="Subject ID",default='S1')
	parser.add_argument('--outdir',"--outdir",dest='Results Directory',help="Output directory",default=startdir,widget="DirChooser")
	parser.add_argument('--echoes',"--echoes",dest='TEs',help="Comma separated TEs in units of ms. Ex: 14.5,38.5,62.5",default='')
	parser.add_argument('--data',"--data",dest='ME-EPI Datasets',help="Choose multiple: e1.nii.gz:e2.nii.gz:e3.nii.gz",nargs='+',widget="MultiFileChooser",default='')
	parser.add_argument('--anat',"--anat",dest='Anatomical',help="Example:   anat.nii.gz",widget="FileChooser",default='')
	parser.add_argument('--basetime',"--basetime",dest='Equilibration Time',help="Example (s for sec, v for vol):   10s OR 10v",default='10s')
	parser.add_argument('--MNI',"--MNI",dest='MNI Normalization',choices=['None', 'Affine', 'Nonlinear (QWarp)'],help="Warp to MNI space using high-resolution template",default='None')
	parser.add_argument('-f',"--fres",dest='EPI Upsampling',help="In mm. Leave blank for none.  Example:   2.5mm", default='2.5mm')
	parser.add_argument('-S',"--no_skullstrip",action="store_true",dest='Disable Automatic Anatomical Skullstrip',help="Anatomical is already skull-stripped and intensity-normalized and skull-stripped",default=False)
	parser.add_argument('-R',"--report",action="store_true",dest='Generate ME-ICA Report',help="Generate HTML report of ME-ICA components and denoising",default=False)
	parser.add_argument('--native',"--native",action="store_true",dest='Native Space Outputs',help="Output result datasets in native space in addition to standard space",default=False)
	args = parser.parse_args()
	return args

def sysrun(thecmd):
	print 'Executing command: ', thecmd
	os.system(thecmd)


if __name__ == '__main__':

	args = main()
	print "CHECKING DEPENDENCIES..."
	ready = dep_check()
	if not ready:  
		print "Dependencies are not OK. Try downloading AFNI and setting on the path, and installing Anaconda Python distribution."
	print

	theTEs = []
	theEPIs = []
	theAnat = None
	outdir = vars(args)['Results Directory']

	#make links to desired start directory, and create links to source files in desired start directory

	#Parse TE string
	if ',' not in args.TEs: 
		print 'ERROR: Please check the syntax of the comma-separated list of TEs'
		print
	args_TEs = args.TEs
	theTEs = args_TEs.replace(' ','').replace('\t','').split(',')
	#print theTEs

	#Parse file string
	theEPIs_in = sorted(vars(args)['ME-EPI Datasets'])
	theEPIs = sorted(list(set([os.path.basename(epi) for epi in theEPIs_in])))
	#print 'the datasets are:', theEPIs_in
	if len(theEPIs) < 3:  
		print 'ERROR: Please select a ME-EPI dataset with >=3 image time series with unique filenames using the file chooser'
		print

	if len(theTEs) != len(theEPIs) and len(theEPIs)<3:
		#raise proper errors here
		print 'ERROR: The number of TEs and the number of ME-EPI datasets do not match!'
		print

	if not os.path.exists(outdir): os.mkdir(outdir)
	print "PREPARING INPUT DATASETS..."
	
	if theAnat!=None:  
		sysrun('ln -s %s %s' % (args.Anatomical,outdir ))
		theAnat = os.path.basename(args.Anatomical)
	for epi in theEPIs_in:
		sysrun('ln -s %s %s' % (epi,outdir ))
	print

	btarg = ''
	mniarg = ''
	fresarg = ''
	nossarg = ''
	natarg = ''
	anatarg = ''
	if vars(args)['Equilibration Time']!='': btarg = '-b %s' % vars(args)['Equilibration Time']
	if theAnat != None:
		if vars(args)['MNI Normalization']=='Affine':   mniarg = '--MNI'
		if vars(args)['MNI Normalization']=='Nonlinear (QWarp)':   mniarg = '--MNI --qwarp'
		if vars(args)['Disable Automatic Anatomical Skullstrip']:  nossarg = '--no_skullstrip'
		if vars(args)['Native Space Outputs']:  natarg = '--native'
		anatarg = '-a %s' % theAnat
	if vars(args)['EPI Upsampling']!='': fresarg = '--fres %s' % vars(args)['EPI Upsampling']

	print "RUNNING ME-ICA:"
	meica_cmd = "cd %s; %s %s/meica.py -e %s -d \"%s\" %s %s %s %s %s %s --OVERWRITE "   %  (outdir, sys.executable,meicadir,','.join(theTEs),','.join(theEPIs),anatarg,btarg,mniarg,fresarg,nossarg,natarg)
	sysrun(meica_cmd)
	print


