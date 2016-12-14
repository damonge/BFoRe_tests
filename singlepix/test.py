import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import scipy.linalg as lng
import common as com

theta_patch=90.
phi_patch=50.
plot_stuff=True
n_samples=10000
n_samples_keep=n_samples*3/4
n_check=1000

nside_spec=16
nside=256
fname_nulist="characteristics_full.txt"
include_polarization=True
beta_s_free=True
beta_d_free=True
temp_d_free=False
independent_polarization=False
include_synchrotron=True
include_dust=True
beta_d_0=1.536735
temp_d_0=20.92867
beta_s_0=-1.
sigma_beta_s=0.01
sigma_beta_d=0.01
sigma_temp_d=0.1

par=com.ParamClean(nside_spec,nside,fname_nulist,
                   include_polarization,independent_polarization,
                   include_synchrotron,include_dust,
                   beta_s_free,beta_d_free,temp_d_free,
                   beta_s_0,beta_d_0,temp_d_0,
                   sigma_beta_s,sigma_beta_d,sigma_temp_d)
print par.n_dof

print "Will vary %d non-linear parameters"%par.n_spec_vary, par.x_spec_0

map_cmb_true=np.zeros([par.n_pix,par.n_pol])
for ipol in np.arange(par.n_pol) :
    map_cmb_true[:,ipol]=hp.read_map("../../r0p00s4321/cmb_r0p00_ns256s4321.fits",field=ipol)

maps_obs=np.zeros([par.n_pix,par.n_pol,par.n_nu])
for inu in np.arange(par.n_nu) :
    for ipol in np.arange(par.n_pol) :
        maps_obs[:,ipol,inu]=hp.read_map("../../r0p00s4321/obs_r0p00_nu%03d.fits"%(inu+1),field=ipol)

amin2_per_pix=4*np.pi*(180*60/np.pi)**2/par.n_pix
sigma2_per_pix=par.noise_list**2/amin2_per_pix
maps_s2_noise=np.ones([par.n_pix,par.n_pol,par.n_nu])*sigma2_per_pix[np.newaxis,np.newaxis,:]
if par.include_polarization :
    maps_s2_noise[:,1:,:]*=2
maps_noise_weights=1./maps_s2_noise

npix_spec=hp.nside2npix(nside_spec)
ipix0=hp.ring2nest(par.nside_spec,hp.ang2pix(par.nside_spec,theta_patch*np.pi/180,phi_patch*np.pi/180))
print "pixel", ipix0
map_mean=np.zeros([par.n_pix,par.n_pol,par.n_comp])
map_sigma=np.zeros([par.n_pix,par.n_pol,par.n_comp])
map_xspec_mean=np.zeros([npix_spec,par.n_spec_vary])
map_xspec_sigma=np.zeros([npix_spec,par.n_spec_vary])

#for ipix in [ipix0] :
for ipix in np.arange(npix_spec) :
    if ipix!=ipix0 :
        continue
    ipix_list=hp.nest2ring(par.nside,ipix*par.n_sub+np.arange(par.n_sub))

    if plot_stuff :
        map_show=np.zeros(par.n_pix); map_show[ipix_list]=1.0; hp.mollview(map_show); plt.show()
        
    patch_cmb_true=map_cmb_true[ipix_list,:]
    patch_obs=maps_obs[ipix_list,:,:]
    patch_noise_weights=maps_noise_weights[ipix_list,:,:]

    if plot_stuff :
        for ipol in np.arange(par.n_pol) :
            plt.title("True CMB %d"%ipol)
            plt.imshow(np.reshape(patch_cmb_true[:,ipol],(par.n_side_sub,par.n_side_sub)),interpolation='none',origin='lower');
            plt.show()

    if plot_stuff==2 :
        for inu in np.arange(par.n_nu) :
            for ipol in np.arange(par.n_pol) :
                plt.title("%d "%inu+"%d"%ipol)
                plt.imshow(np.reshape(patchs_obs[:,ipol,inu],(par.n_side_sub,par.n_side_sub)),interpolation='none',origin='lower');
                plt.show()

    t_mean=np.zeros([par.n_sub,par.n_pol,par.n_comp])
    t_sigma=np.zeros([par.n_sub,par.n_pol,par.n_comp])
    x_spec_samples=np.zeros([n_samples,par.n_param_max])
    x_spec_old=par.x_spec_0.copy()
    x_spec_mean=np.zeros_like(par.x_spec_0)
    x_spec_sigma=np.zeros_like(par.x_spec_0)
    n_accepted=0.
    sigma_vary=par.sigma_spec.copy()
    rate_target=0.2
    for i in np.arange(n_samples) :
        if i%n_check==0 :
            print "Step # %d :"%i
            toprint=True
        else :
            toprint=False
    #Gibbs sampling
        t_new=com.sample_amplitudes(patch_obs,patch_noise_weights,x_spec_old,par) #A_{n+1}(b_{n})

        accept,x_spec_new=com.sample_indices(patch_obs,patch_noise_weights,t_new,x_spec_old,sigma_vary,
                                             toprint,par) #b_{n+1}(A_{n+1})

        n_accepted+=accept
        x_spec_old[:]=x_spec_new[:]

        if i%n_check==0 :
            rate=float(n_accepted)/n_check
            print "   Acceptance rate: %lE"%rate
            print "   Current step size: ", sigma_vary
            n_accepted=0.
#            if i<(n_samples-n_samples_keep) : #Still burning
#                if rate>1.0E-5 :
#                    sigma_vary*=rate/rate_target
#                else :
#                    sigma_vary*=0.5

        x_spec_samples[i,:]=x_spec_old
        if i>=(n_samples-n_samples_keep) :
            x_spec_mean+=x_spec_old
            x_spec_sigma+=x_spec_old**2
            t_mean+=t_new
            t_sigma+=t_new**2

    t_mean/=n_samples_keep
    x_spec_mean/=n_samples_keep
    t_sigma=np.sqrt(t_sigma/n_samples_keep-t_mean**2)
    x_spec_sigma=np.sqrt(x_spec_sigma/n_samples_keep-x_spec_mean**2)
    print x_spec_mean, x_spec_sigma
    map_mean[ipix_list,:,:]=t_mean
    map_sigma[ipix_list,:,:]=t_sigma
    map_xspec_mean[ipix,:]=x_spec_mean[:par.n_spec_vary]
    map_xspec_sigma[ipix,:]=x_spec_sigma[:par.n_spec_vary]
            
    np.savetxt("xspec.txt"%ipix,x_spec_samples)

    print com.compute_chi2(patch_obs,patch_noise_weights,t_mean,x_spec_mean,par)/par.n_dof, par.n_dof
    print np.std(patch_cmb_true-t_mean[:,:,par.index_cmb],axis=0)*np.sqrt(amin2_per_pix)
    print np.sqrt(np.mean(t_sigma[:,:,par.index_cmb]**2,axis=0))*np.sqrt(amin2_per_pix)
