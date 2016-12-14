import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import scipy.linalg as lng

class ParamClean:
    nside=256
    nside_spec=16
    fname_nulist="characteristics.txt"
    nu_list=None
    noise_list=None
    include_polarization=False
    beta_s_free=False
    beta_d_free=False
    temp_d_free=False
    independent_polarization=False
    include_synchrotron=False
    include_dust=False
    beta_s_0=-1.
    beta_d_0=1.536735
    temp_d_0=20.92867
    sigma_beta_s=0.01
    sigma_beta_d=0.01
    sigma_temp_d=0.3
    n_pol=1
    n_comp=1
    index_synchrotron=None
    index_dust=None
    n_param_max=None
    n_spec_vary=None
    index_beta_s_t=None
    index_beta_s_p=None
    index_beta_d_t=None
    index_beta_d_p=None
    index_temp_d_t=None
    index_temp_d_p=None
    n_dof=None
    n_side_sub=16
    n_sub=256
    n_pix=786432
    n_nu=None
    x_spec_0=None
    sigma_spec=None

    def __init__(self,nside_spec,nside,fname_nulist,
                 include_polarization,independent_polarization,
                 include_synchrotron,include_dust,
                 beta_s_free,beta_d_free,temp_d_free,
                 beta_s_0,beta_d_0,temp_d_0,
                 sigma_beta_s,sigma_beta_d,sigma_temp_d) :
        self.nside_spec=nside_spec
        self.nside=nside
        self.fname_nulist=fname_nulist
        self.nu_list,self.noise_list=np.loadtxt(fname_nulist,unpack=True)[[1,3]]
        self.include_polarization=include_polarization
        self.beta_s_free=beta_s_free
        self.beta_d_free=beta_d_free
        self.temp_d_free=temp_d_free
        self.independent_polarization=independent_polarization
        self.include_synchrotron=include_synchrotron
        self.include_dust=include_dust
        self.beta_s_0=beta_s_0
        self.beta_d_0=beta_d_0
        self.temp_d_0=temp_d_0
        self.sigma_beta_s=sigma_beta_s
        self.sigma_beta_d=sigma_beta_d
        self.sigma_temp_d=sigma_temp_d

        if self.include_polarization :
            self.n_pol=3
        else :
            self.n_pol=1

        self.n_comp=0
        self.index_cmb=self.n_comp
        self.n_comp+=1
        if self.include_synchrotron :
            self.index_synchrotron=self.n_comp
            self.n_comp+=1
        if self.include_dust :
            self.index_dust=self.n_comp
            self.n_comp+=1

        self.n_param_max=3
        if self.independent_polarization :
            self.n_param_max*=2
            
        self.n_spec_vary=0
        i_spec_novary=self.n_param_max
        if self.include_synchrotron and self.beta_s_free :
            self.index_beta_s_t=self.n_spec_vary; self.n_spec_vary+=1
        else :
            i_spec_novary-=1; self.index_beta_s_t=i_spec_novary
        if self.include_dust and self.beta_d_free :
            self.index_beta_d_t=self.n_spec_vary; self.n_spec_vary+=1
        else :
            i_spec_novary-=1; self.index_beta_d_t=i_spec_novary
        if self.include_dust and self.temp_d_free :
            self.index_temp_d_t=self.n_spec_vary; self.n_spec_vary+=1
        else :
            i_spec_novary-=1; self.index_temp_d_t=i_spec_novary
        if independent_polarization :
            if self.include_synchrotron and self.beta_s_free :
                self.index_beta_s_p=self.n_spec_vary; self.n_spec_vary+=1
            else :
                i_spec_novary-=1; self.index_beta_s_p=i_spec_novary
            if self.include_dust and self.beta_d_free :
                self.index_beta_d_p=self.n_spec_vary; self.n_spec_vary+=1
            else :
                i_spec_novary-=1; self.index_beta_d_p=i_spec_novary
            if self.include_dust and self.temp_d_free :
                self.index_temp_d_p=self.n_spec_vary; self.n_spec_vary+=1
            else :
                i_spec_novary-=1; self.index_temp_d_p=i_spec_novary
        else :
            self.index_beta_s_p=self.index_beta_s_t
            self.index_beta_d_p=self.index_beta_d_t
            self.index_temp_d_p=self.index_temp_d_t

        self.n_side_sub=(self.nside/self.nside_spec)
        self.n_sub=self.n_side_sub**2
        self.n_pix=hp.nside2npix(self.nside)
        self.n_nu=len(self.nu_list)
        self.n_dof=self.n_sub*self.n_pol*(self.n_nu-self.n_comp)-self.n_spec_vary

        self.x_spec_0=np.zeros(self.n_param_max)
        self.sigma_spec=np.zeros(self.n_param_max)
        self.x_spec_0[self.index_beta_s_t]=self.beta_s_0
        self.x_spec_0[self.index_beta_s_p]=self.beta_s_0
        self.x_spec_0[self.index_beta_d_t]=self.beta_d_0
        self.x_spec_0[self.index_beta_d_p]=self.beta_d_0
        self.x_spec_0[self.index_temp_d_t]=self.temp_d_0
        self.x_spec_0[self.index_temp_d_p]=self.temp_d_0
        self.sigma_spec[self.index_beta_s_t]=self.sigma_beta_s
        self.sigma_spec[self.index_beta_s_p]=self.sigma_beta_s
        self.sigma_spec[self.index_beta_d_t]=self.sigma_beta_d
        self.sigma_spec[self.index_beta_d_p]=self.sigma_beta_d
        self.sigma_spec[self.index_temp_d_t]=self.sigma_temp_d
        self.sigma_spec[self.index_temp_d_p]=self.sigma_temp_d

def freq_evolve(spec_type,nu_0,beta,temp,nu) :
    if spec_type=="BB" : #CMB
        x=0.017611907*nu
        ex=np.exp(x)
        return ex*(x/(ex-1))**2
    elif spec_type=="PL" : #Synch
        return (nu/nu_0)**(beta-2.)
    elif spec_type=="mBB" : #Dust
        x_to=0.0479924466*nu/temp
        x_from=0.0479924466*nu_0/temp
        return (nu/nu_0)**(1+beta)*(np.exp(x_from)-1)/(np.exp(x_to)-1)

def get_evolution_matrices(x_spec,par) :
    f_mat=np.zeros([par.n_pol,par.n_comp,par.n_nu])
    #CMB
    f_mat[:,par.index_cmb,:]=freq_evolve("BB",None,None,None,par.nu_list)[np.newaxis,:]
    #Synchrotron
    if par.include_synchrotron :
        f_mat[0,par.index_synchrotron,:]=freq_evolve("PL",23.,x_spec[par.index_beta_s_t],None,par.nu_list)
        if par.include_polarization :
            f_mat[1:,par.index_synchrotron,:]=freq_evolve("PL",23.,x_spec[par.index_beta_s_p],None,par.nu_list)[np.newaxis,:]
    #Dust
    if par.include_dust :
        f_mat[0 ,par.index_dust,:]=freq_evolve("mBB",353.,x_spec[par.index_beta_d_t],x_spec[par.index_temp_d_t],par.nu_list)
        if par.include_polarization :
            f_mat[1:,par.index_dust,:]=freq_evolve("mBB",353.,x_spec[par.index_beta_d_p],
                                                   x_spec[par.index_temp_d_p],par.nu_list)[np.newaxis,:]

    return f_mat

def compute_chi2(map_data,map_w2,map_amp,x_spec,par) :
    f_mat=get_evolution_matrices(x_spec,par) #[npol,ncomp,nnu]
                  #[npol,ncomp,nnu] [npix,npol,ncomp] -> [npix,npol,nnu]
    res=map_data-np.sum(f_mat[np.newaxis,:,:,:]*map_amp[:,:,:,np.newaxis],axis=2)
    
    return np.sum(res**2*map_w2)

def analyze_linear_chi2(map_data,map_w2,x_spec,par) :
    f_mat=get_evolution_matrices(x_spec,par)
    #                           [npol,ncomp,ncomp,nnu]                                     [npix,npol,nnu]   = [npix,npol,ncomp,ncomp]
    c_inv=np.sum((f_mat[:,np.newaxis,:,:]*f_mat[:,:,np.newaxis,:])[np.newaxis,:,:,:,:]*map_w2[:,:,np.newaxis,np.newaxis,:],axis=4)
    #          [npix,npol,nnu]               [npol,ncomp,nnu]
    vec_0=np.sum((map_data*map_w2)[:,:,np.newaxis,:]*f_mat[np.newaxis,:,:,:],axis=3) #[npix,npol,ncomp]
    t_mean=np.linalg.solve(c_inv,vec_0)

    return c_inv,t_mean

def sample_amplitudes(map_data,map_w2,x_spec,par) :
    c_inv,t_mean=analyze_linear_chi2(map_data,map_w2,x_spec,par) #Mean and inverse covariance of linear system
    t_extra=np.linalg.solve(np.transpose(np.linalg.cholesky(c_inv),[0,1,3,2]),np.random.randn(par.n_sub,par.n_pol,par.n_comp))
    #    u=np.random.randn(n_sub,n_pol,n_comp) #[npix,npol,ncomp]
    #    t_extra=np.sum(np.linalg.cholesky(np.linalg.inv(c_inv))*u[:,:,np.newaxis,:],axis=3)

    return t_mean+t_extra

def sample_indices(map_data,map_w2,map_amp,x_spec_old,sigma_vary,do_print,par) :
    x_spec_new=x_spec_old.copy()
    x_spec_new[:par.n_spec_vary]+=sigma_vary[:par.n_spec_vary]*np.random.randn(par.n_spec_vary)
    chi2_old=compute_chi2(map_data,map_w2,map_amp,x_spec_old,par)
    chi2_new=compute_chi2(map_data,map_w2,map_amp,x_spec_new,par)
    r=np.exp(-0.5*(chi2_new-chi2_old))
#    if do_print :
#    print x_spec_old[:par.n_spec_vary], x_spec_new[:par.n_spec_vary], chi2_old, chi2_new, r
    replaced=1.
    if r<1 :
        if np.random.rand(1)>r :
            x_spec_new[:]=x_spec_old[:]
            replaced=0.

    return replaced,x_spec_new
