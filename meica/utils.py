import numpy as np
import nibabel as nib
from scipy.stats import scoreatpercentile


def niwrite(data, affine, name, head, header=None):
    data[np.isnan(data)]=0
    if header == None:
        this_header = head.copy()
        this_header.set_data_shape(list(data.shape))
    else:
        this_header = header

    outni = nib.Nifti1Image(data,affine,header=this_header)
    outni.set_data_dtype('float64')
    outni.to_filename(name)


def cat2echos(data,Ne):
    """
    cat2echos(data,Ne)

    Input:
    data shape is (nx,ny,Ne*nz,nt)
    """
    nx,ny = data.shape[0:2]
    nz = data.shape[2]//Ne
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
        medv = emeans[:,:,:,0] == scoreatpercentile(emeans[:,:,:,0][emeans[:,:,:,0]!=0],33,interpolation_method='higher')
        lthrs = np.squeeze(np.array([ emeans[:,:,:,ee][medv]/3 for ee in range(Ne) ]))

        if len(lthrs.shape)==1:
            lthrs = np.atleast_2d(lthrs).T
        lthrs = lthrs[:,lthrs.sum(0).argmax()]

        mthr = np.ones([nx,ny,nz,Ne])
        for ee in range(Ne): mthr[:,:,:,ee]*=lthrs[ee]
        mthr = np.abs(emeans[:,:,:,:])>mthr
        masksum = np.array(mthr,dtype=np.int).sum(-1)
        mask = masksum!=0
        if getsum:
            return mask,masksum
        else:
            return mask


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
