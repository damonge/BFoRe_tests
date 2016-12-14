import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import sys as sys

fname_read=sys.argv[1]

def read_debug(fname) :
    f=open(fname,"rb")
    fsize,ipix,nside,nside_spec,n_sub,n_nu,n_pol,n_comp,n_spec_vary,n_samples,nb_samples=np.fromfile(f,dtype=np.int32,count=11)
    if fsize==4 :
        print "is float"
        dt=np.float32
    else :
        print "is double"
        dt=np.float64
    n_side_sub=nside/nside_spec
    if n_side_sub**2!=n_sub :
        print "shit"
        exit(1)
    amin2perpix=4*np.pi*(180*60/np.pi)**2/hp.nside2npix(nside)

    chain=(np.fromfile(f,dtype=dt,count=n_spec_vary*(n_samples+nb_samples))).reshape(n_samples+nb_samples,n_spec_vary)
    
    indices_mean=np.fromfile(f,dtype=dt,count=n_spec_vary)
    indices_covar=(np.fromfile(f,dtype=dt,count=n_spec_vary*n_spec_vary)).reshape(n_spec_vary,n_spec_vary)
    amp_mean=(np.fromfile(f,dtype=dt,count=n_comp*n_pol*n_sub)).reshape(n_side_sub,n_side_sub,n_pol,n_comp)
    amp_covar=(np.fromfile(f,dtype=dt,count=n_comp*n_comp*n_pol*n_sub)).reshape(n_side_sub,n_side_sub,n_pol,n_comp,n_comp)
    input_data=(np.fromfile(f,dtype=dt,count=n_nu*n_pol*n_sub)).reshape(n_side_sub,n_side_sub,n_pol,n_nu)
    input_noise=(np.fromfile(f,dtype=dt,count=n_nu*n_pol*n_sub)).reshape(n_side_sub,n_side_sub,n_pol,n_nu)
    
    for i in np.arange(n_spec_vary) :
        plt.hist(chain[nb_samples:,i],bins=100); plt.show()
        plt.plot(chain[nb_samples:,i]); plt.xlabel('#sample',fontsize=16); plt.ylabel('$\\beta_s$',fontsize=16); plt.show()
        for j in np.arange(n_spec_vary-i-1)+i+1 :
            plt.plot(chain[:nb_samples,i],chain[:nb_samples,j],'r-',markersize=0.1);
            plt.plot(chain[nb_samples:,i],chain[nb_samples:,j],'b.',markersize=1); plt.show()

    f.close()

    print "Mean and rms:", indices_mean, np.sqrt(np.diag(indices_covar))
    print "Covar:"
    print indices_covar
  
    print "Mean amplitudes"
    for icomp in np.arange(n_comp) :
        for ipol in np.arange(n_pol) :
            plt.title("%d "%icomp+"%d"%ipol)
            plt.imshow(amp_mean[:,:,ipol,icomp],interpolation='none',origin='lower')
            plt.show()

    print "RMS amplitudes"
    for icomp in np.arange(n_comp) :
        for ipol in np.arange(n_pol) :
            plt.title("%d "%icomp+"%d"%ipol)
            plt.imshow(amp_covar[:,:,ipol,icomp,icomp],interpolation='none',origin='lower')
            print np.sqrt(np.mean(amp_covar[:,:,ipol,icomp,icomp])*amin2perpix)
            plt.show()

read_debug(fname_read)
