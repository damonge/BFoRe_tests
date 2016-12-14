#include "common.h"

static void decouple_covariance(ParamBFoRe *par,gsl_matrix *cov)
{
  if(par->flag_include_polarization && par->flag_independent_polarization) {
    if(par->flag_beta_s_free) {
      gsl_matrix_set(cov,par->index_beta_s_t,par->index_beta_s_p,0); //bst-bsp
      gsl_matrix_set(cov,par->index_beta_s_p,par->index_beta_s_t,0); //bsp-bst
      if(par->flag_beta_d_free) {
	gsl_matrix_set(cov,par->index_beta_s_t,par->index_beta_d_p,0); //bst-bdp
	gsl_matrix_set(cov,par->index_beta_d_p,par->index_beta_s_t,0); //bdp-bst
	gsl_matrix_set(cov,par->index_beta_s_p,par->index_beta_d_t,0); //bsp-bdt
	gsl_matrix_set(cov,par->index_beta_d_t,par->index_beta_s_p,0); //bdt-bsp
      }
      if(par->flag_temp_d_free) {
	gsl_matrix_set(cov,par->index_beta_s_t,par->index_temp_d_p,0); //bst,tdp
	gsl_matrix_set(cov,par->index_temp_d_p,par->index_beta_s_t,0); //tdp,bst
	gsl_matrix_set(cov,par->index_beta_s_p,par->index_temp_d_t,0); //bsp,tdt
	gsl_matrix_set(cov,par->index_temp_d_t,par->index_beta_s_p,0); //tdt,bsp
      }
    }
    if(par->flag_beta_d_free) {
      gsl_matrix_set(cov,par->index_beta_d_t,par->index_beta_d_p,0); //bdt-bdp
      gsl_matrix_set(cov,par->index_beta_d_p,par->index_beta_d_t,0); //bdp-bdt
      if(par->flag_temp_d_free) {
	gsl_matrix_set(cov,par->index_beta_d_t,par->index_temp_d_p,0); //bdt-tdp
	gsl_matrix_set(cov,par->index_temp_d_p,par->index_beta_d_t,0); //tdp,bdt
	gsl_matrix_set(cov,par->index_beta_d_p,par->index_temp_d_t,0); //bdp,tdt
	gsl_matrix_set(cov,par->index_temp_d_t,par->index_beta_d_p,0); //tdt,bdp
      }
    }
    if(par->flag_temp_d_free) {
      gsl_matrix_set(cov,par->index_temp_d_t,par->index_temp_d_p,0); //tdt-tdp
      gsl_matrix_set(cov,par->index_temp_d_p,par->index_temp_d_t,0); //tdp-tdt
    }
  }
}

static void init_priors(ParamBFoRe *par,PixelState *pst,int ipix_big)
{
  pst->prior_mean[par->index_beta_s_t]=par->map_prior_centres[par->index_beta_s_t+par->n_param_max*ipix_big];
  if(par->flag_beta_s_free) {
    pst->prior_isigma[par->index_beta_s_t]=1./par->map_prior_widths[par->index_beta_s_t+
								    par->n_param_max*ipix_big];
  }    
  pst->prior_mean[par->index_beta_d_t]=par->map_prior_centres[par->index_beta_d_t+par->n_param_max*ipix_big];
  if(par->flag_beta_d_free) {
    pst->prior_isigma[par->index_beta_d_t]=1./par->map_prior_widths[par->index_beta_d_t+
								    par->n_param_max*ipix_big];
  }    
  pst->prior_mean[par->index_temp_d_t]=par->map_prior_centres[par->index_temp_d_t+par->n_param_max*ipix_big];
  if(par->flag_temp_d_free) {
    pst->prior_isigma[par->index_temp_d_t]=1./par->map_prior_widths[par->index_temp_d_t+
								    par->n_param_max*ipix_big];
  }    
  if(par->flag_include_polarization && par->flag_independent_polarization) {
    pst->prior_mean[par->index_beta_s_p]=par->map_prior_centres[par->index_beta_s_p+par->n_param_max*ipix_big];
    if(par->flag_beta_s_free) {
      pst->prior_isigma[par->index_beta_s_p]=1./par->map_prior_widths[par->index_beta_s_p+
								      par->n_param_max*ipix_big];
    }    
    pst->prior_mean[par->index_beta_d_p]=par->map_prior_centres[par->index_beta_d_p+par->n_param_max*ipix_big];
    if(par->flag_beta_d_free) {
      pst->prior_isigma[par->index_beta_d_p]=1./par->map_prior_widths[par->index_beta_d_p+
								      par->n_param_max*ipix_big];
    }    
    pst->prior_mean[par->index_temp_d_p]=par->map_prior_centres[par->index_temp_d_p+par->n_param_max*ipix_big];
    if(par->flag_temp_d_free) {
      pst->prior_isigma[par->index_temp_d_p]=1./par->map_prior_widths[par->index_temp_d_p+
								      par->n_param_max*ipix_big];
    }
  }
}

PixelState *pixel_state_new(ParamBFoRe *par)
{
  int ip;
  PixelState *pst=my_malloc(sizeof(PixelState));
  pst->f_matrix=my_malloc(par->n_pol*par->n_nu*par->n_comp*sizeof(flouble));
  pst->cov_inv=my_malloc(par->n_sub*par->n_pol*sizeof(gsl_matrix *));
  pst->vec_mean=my_malloc(par->n_sub*par->n_pol*sizeof(gsl_vector *));
  pst->prior_mean=my_malloc(par->n_param_max*sizeof(flouble));
  pst->prior_isigma=my_malloc(par->n_param_max*sizeof(flouble));
  for(ip=0;ip<par->n_sub*par->n_pol;ip++) {
    pst->cov_inv[ip]=gsl_matrix_alloc(par->n_comp,par->n_comp);
    pst->vec_mean[ip]=gsl_vector_alloc(par->n_comp);
  }
  pst->rand_spec=my_malloc(par->n_spec_vary*sizeof(flouble));
  pst->vaux=gsl_vector_alloc(par->n_comp);

  return pst;
}

void pixel_state_free(PixelState *pst,ParamBFoRe *par)
{
  int ip;
  free(pst->f_matrix);
  for(ip=0;ip<par->n_sub*par->n_pol;ip++) {
    gsl_matrix_free(pst->cov_inv[ip]);
    gsl_vector_free(pst->vec_mean[ip]);
  }
  free(pst->cov_inv);
  free(pst->vec_mean);
  free(pst->rand_spec);
  free(pst->prior_mean);
  free(pst->prior_isigma);
  gsl_vector_free(pst->vaux);
  free(pst);
}

static flouble freq_evolve(int spec_type,double nu_0,double beta,double temp,double nu)
{
  flouble x_to,x_from,ex;
  switch(spec_type)  {
  case 0 :
    x_to=0.017611907*nu;
    ex=exp(x_to);
    x_to=x_to/(ex-1);
    return ex*x_to*x_to;
    break;
  case 1 :
    return pow(nu/nu_0,beta-2.);
    break;
  case 2 :
    x_to=0.0479924466*nu/temp; //DAM: possible optimization, use 1/T instead of T
    x_from=0.0479924466*nu_0/temp;
    return pow(nu/nu_0,beta+1.)*(exp(x_from)-1)/(exp(x_to)-1);
  }
  return -1;
}

static void compute_f_matrix(ParamBFoRe *par,flouble *x_spec,flouble *f_matrix)
{
  int inu;
  for(inu=0;inu<par->n_nu;inu++) {
    int ipol;
    flouble nu=par->freqs[inu];
    if(par->flag_include_cmb) {
      for(ipol=0;ipol<par->n_pol;ipol++) {
	f_matrix[par->index_cmb+par->n_comp*(inu+ipol*par->n_nu)]=
	  freq_evolve(0,-1,-1,-1,nu);
      }
    }
    if(par->flag_include_synchrotron) {
      f_matrix[par->index_synchrotron+par->n_comp*(inu+0*par->n_nu)]=
	freq_evolve(1,par->nu0_s,x_spec[par->index_beta_s_t],-1,nu);
      for(ipol=1;ipol<par->n_pol;ipol++) {
	f_matrix[par->index_synchrotron+par->n_comp*(inu+ipol*par->n_nu)]=
	  freq_evolve(1,par->nu0_s,x_spec[par->index_beta_s_p],-1,nu);
      }
    }
    if(par->flag_include_dust) {
      f_matrix[par->index_dust+par->n_comp*(inu+0*par->n_nu)]=
	freq_evolve(2,par->nu0_d,x_spec[par->index_beta_d_t],
		    x_spec[par->index_temp_d_t],nu);
      for(ipol=1;ipol<par->n_pol;ipol++) {
	f_matrix[par->index_dust+par->n_comp*(inu+ipol*par->n_nu)]=
	  freq_evolve(2,par->nu0_d,x_spec[par->index_beta_d_p],
		      x_spec[par->index_temp_d_p],nu);
      }
    }
  }
}

static flouble chi2_prior_correctvolume(ParamBFoRe *par,PixelState *pst,
					flouble *noise_w,flouble *x_spec)
{
  int ipix;
  double chi2=0;

  compute_f_matrix(par,x_spec,pst->f_matrix);

  for(ipix=0;ipix<par->n_sub;ipix++) {
    int ipol;
    for(ipol=0;ipol<par->n_pol;ipol++) {
      int inu,ic1;
      int index_pix=ipol+par->n_pol*ipix;
      gsl_matrix *mat_here=pst->cov_inv[index_pix];
      gsl_matrix_set_zero(mat_here);
      for(inu=0;inu<par->n_nu;inu++) {
	flouble invsigma2=noise_w[inu+par->n_nu*index_pix];
	for(ic1=0;ic1<par->n_comp;ic1++) {
	  int ic2;
	  flouble fm1=pst->f_matrix[ic1+par->n_comp*(inu+par->n_nu*ipol)];
	  for(ic2=0;ic2<=ic1;ic2++) {
	    flouble fm2=pst->f_matrix[ic2+par->n_comp*(inu+par->n_nu*ipol)];
	    gsl_matrix_set(mat_here,ic1,ic2,gsl_matrix_get(mat_here,ic1,ic2)+fm1*fm2*invsigma2);
	  }
	}
      }
      
      gsl_linalg_cholesky_decomp(mat_here);
      for(ic1=0;ic1<par->n_comp;ic1++)
	chi2+=log(gsl_matrix_get(mat_here,ic1,ic1));
    }
  }

  return -2*chi2;
}

static flouble chi2_prior(ParamBFoRe *par,PixelState *pst,flouble *x_spec)
{
  int ipar;
  flouble x,chi2=0;

  for(ipar=0;ipar<par->n_spec_vary;ipar++) {
    if(pst->prior_isigma[ipar]>0) {
      x=(x_spec[ipar]-pst->prior_mean[ipar])*pst->prior_isigma[ipar];
      chi2+=x*x;
    }
  }

  return chi2;
}

static double compute_chi2(ParamBFoRe *par,flouble *data,flouble *noise_w,
			   flouble *amps,flouble *x_spec,PixelState *pst)
{
  int ipix;
  double chi2=0;

  compute_f_matrix(par,x_spec,pst->f_matrix);

  for(ipix=0;ipix<par->n_sub;ipix++) {
    int ipol;
    for(ipol=0;ipol<par->n_pol;ipol++) {
      int inu;
      int index_pix=ipol+par->n_pol*ipix;
      for(inu=0;inu<par->n_nu;inu++) {
	int icomp;
	double res=(double)(data[inu+par->n_nu*index_pix]);
	for(icomp=0;icomp<par->n_comp;icomp++) {
	  res-=pst->f_matrix[icomp+par->n_comp*(inu+par->n_nu*ipol)]*
	    amps[icomp+par->n_comp*index_pix];
	}
	chi2+=res*res*noise_w[inu+par->n_nu*index_pix];
      }
    }
  }

  return chi2;
}

static void compute_marginalized_chi2(ParamBFoRe *par,flouble *data,flouble *noise_w,
				      flouble *x_spec,PixelState *pst,flouble subtract)
{
  int ipix;
  flouble offset=subtract/(par->n_sub*par->n_pol*par->n_comp);

  pst->chi2=chi2_prior(par,pst,x_spec);
  compute_f_matrix(par,x_spec,pst->f_matrix);

  for(ipix=0;ipix<par->n_sub;ipix++) {
    int ipol;
    for(ipol=0;ipol<par->n_pol;ipol++) {
      int ic1,inu;
      int index_pix=ipol+par->n_pol*ipix;
      gsl_matrix *mat_here=pst->cov_inv[index_pix];
      gsl_vector *vec_here=pst->vec_mean[index_pix];
      gsl_matrix_set_zero(mat_here);
      gsl_vector_set_zero(pst->vaux);
      for(inu=0;inu<par->n_nu;inu++) {
	int index_f=par->n_comp*(inu+par->n_nu*ipol);
	flouble invsigma2=noise_w[inu+par->n_nu*index_pix];
	flouble data_here=data[inu+par->n_nu*index_pix];
	for(ic1=0;ic1<par->n_comp;ic1++) {
	  int ic2;
	  flouble fm1=pst->f_matrix[index_f+ic1];
	  gsl_vector_set(pst->vaux,ic1,gsl_vector_get(pst->vaux,ic1)+fm1*data_here*invsigma2);
	  for(ic2=0;ic2<=ic1;ic2++) {
	    flouble fm2=pst->f_matrix[index_f+ic2];
	    gsl_matrix_set(mat_here,ic1,ic2,gsl_matrix_get(mat_here,ic1,ic2)+fm1*fm2*invsigma2);
	  }
	}
      }
      gsl_linalg_cholesky_decomp(mat_here);
      if(par->flag_include_volume_prior==0) {
	for(ic1=0;ic1<par->n_comp;ic1++)
	  pst->chi2+=2*log(gsl_matrix_get(mat_here,ic1,ic1)); //0.5*log(det((F^T N^-1 F)^-1))
      }
      gsl_linalg_cholesky_invert(mat_here);
      gsl_blas_dgemv(CblasNoTrans,1.,mat_here,pst->vaux,0,vec_here); //v_mean = (F^T N^-1 F)^-1 F^T N^-1 d
      for(ic1=0;ic1<par->n_comp;ic1++)
	pst->chi2-=gsl_vector_get(pst->vaux,ic1)*gsl_vector_get(vec_here,ic1)+offset; //chi2= (F^T N^-1 d)^T (F^T N^-1 F)^-1 (F^T N^-1 d)
    }
  }
}

static void analyze_linear_chi2(ParamBFoRe *par,flouble *data,flouble *noise_w,
				flouble *x_spec,PixelState *pst)
{
  int ipix;

  compute_f_matrix(par,x_spec,pst->f_matrix);

  for(ipix=0;ipix<par->n_sub;ipix++) {
    int ipol;
    for(ipol=0;ipol<par->n_pol;ipol++) {
      int ic1;
      int index_pix=ipol+par->n_pol*ipix;
      for(ic1=0;ic1<par->n_comp;ic1++) {
	int ic2,inu;
	flouble vec=0;
	for(inu=0;inu<par->n_nu;inu++) {
	  vec+=
	    pst->f_matrix[ic1+par->n_comp*(inu+par->n_nu*ipol)]*
	    data[inu+par->n_nu*index_pix]*noise_w[inu+par->n_nu*index_pix]; //v = F^T N^-1 d
	}
	gsl_vector_set(pst->vec_mean[index_pix],ic1,vec);
	for(ic2=0;ic2<par->n_comp;ic2++) {
	  flouble icov=0;
	  for(inu=0;inu<par->n_nu;inu++) {
	    icov+=
	      pst->f_matrix[ic1+par->n_comp*(inu+par->n_nu*ipol)]*
	      pst->f_matrix[ic2+par->n_comp*(inu+par->n_nu*ipol)]*
	      noise_w[inu+par->n_nu*index_pix]; //C^-1 = F^T N^-1 F
	  }
	  gsl_matrix_set(pst->cov_inv[index_pix],ic1,ic2,icov);
	}
      }
      gsl_linalg_cholesky_decomp(pst->cov_inv[index_pix]);
      gsl_linalg_cholesky_svx(pst->cov_inv[index_pix],pst->vec_mean[index_pix]); //v = (F^T N^-1 F)^-1 F^T N^-1 d
    }
  }
}

//static void solve_lower_triangular(int n,gsl_matrix *mat,flouble *v)
//{
//  int i;
//  for(i=0;i<n;i++) {
//    int j;
//    flouble res=v[i];
//    for(j=0;j<i;j++)
//      res-=v[j]*gsl_matrix_get(mat,i,j);
//    v[i]=res/gsl_matrix_get(mat,i,i);
//  }
//}

static void solve_upper_triangular(int n,gsl_matrix *mat,flouble *v)
{
  int i;
  for(i=n-1;i>=0;i--) {
    int j;
    flouble res=v[i];
    for(j=n-1;j>i;j--)
      res-=gsl_matrix_get(mat,i,j)*v[j];
    v[i]=res/gsl_matrix_get(mat,i,i);
  }
}

static void draw_amplitudes(ParamBFoRe *par,Rng *rng,flouble *data,flouble *noise_w,
			    flouble *x_spec,PixelState *pst,flouble *amps)
{
  int ipix;

  analyze_linear_chi2(par,data,noise_w,x_spec,pst);
  
  for(ipix=0;ipix<par->n_sub;ipix++) {
    int ipol;
    for(ipol=0;ipol<par->n_pol;ipol++) {
      int icomp;
      int index_pix=ipol+par->n_pol*ipix;
      flouble *amp_here=&(amps[par->n_comp*index_pix]);
      for(icomp=0;icomp<par->n_comp;icomp++)
	amp_here[icomp]=rand_gauss(rng);
      solve_upper_triangular(par->n_comp,pst->cov_inv[index_pix],amp_here); //t_rand = L^T^-1 u
      for(icomp=0;icomp<par->n_comp;icomp++)
	amp_here[icomp]+=gsl_vector_get(pst->vec_mean[index_pix],icomp); //t = t_mean + t_rand
    }
  }
}

static int draw_spectral_indices(ParamBFoRe *par,Rng *rng,flouble *data,flouble *noise_w,
				 flouble *amps,flouble *x_spec_old,PixelState *pst,
				 gsl_matrix *mat_step,flouble *x_spec_new)
{ //DAM: possible optimization: demote mat_step to flouble *
  int ipar;
  double chi2_old,chi2_new,ratio;
  memcpy(x_spec_new,x_spec_old,par->n_param_max*sizeof(flouble));
  for(ipar=0;ipar<par->n_spec_vary;ipar++)
    pst->rand_spec[ipar]=rand_gauss(rng);
  for(ipar=0;ipar<par->n_spec_vary;ipar++) {
    int ipar2;
    for(ipar2=0;ipar2<=ipar;ipar2++)
      x_spec_new[ipar]+=gsl_matrix_get(mat_step,ipar,ipar2)*pst->rand_spec[ipar2];
  }

  chi2_new=compute_chi2(par,data,noise_w,amps,x_spec_new,pst)+chi2_prior(par,pst,x_spec_new);
  chi2_old=compute_chi2(par,data,noise_w,amps,x_spec_old,pst)+chi2_prior(par,pst,x_spec_old);
  if(par->flag_include_volume_prior) {
    chi2_new+=chi2_prior_correctvolume(par,pst,noise_w,x_spec_new);
    chi2_old+=chi2_prior_correctvolume(par,pst,noise_w,x_spec_old);
  }

  ratio=exp(-0.5*(chi2_new-chi2_old));

  if(ratio<1) {
    if(rand_real01(rng)>ratio)
      return 0;
  }

  return 1;
}

//Sets x_spec to prior mean
//Sets step matrix to diagonal with m[i,i]= step_i*factor
static void restart_mcmc(ParamBFoRe *par,PixelState *pst,flouble *x_spec,
			 gsl_matrix *mat_step,flouble factor)
{
  x_spec[par->index_beta_s_t]=pst->prior_mean[par->index_beta_s_t];
  x_spec[par->index_beta_d_t]=pst->prior_mean[par->index_beta_d_t];
  x_spec[par->index_temp_d_t]=pst->prior_mean[par->index_temp_d_t];
  if(par->n_spec_vary>0) {
    gsl_matrix_set_zero(mat_step);
    if(par->flag_beta_s_free)
      gsl_matrix_set(mat_step,par->index_beta_s_t,par->index_beta_s_t,par->beta_s_step*factor);
    if(par->flag_beta_d_free)
      gsl_matrix_set(mat_step,par->index_beta_d_t,par->index_beta_d_t,par->beta_d_step*factor);
    if(par->flag_temp_d_free)
      gsl_matrix_set(mat_step,par->index_temp_d_t,par->index_temp_d_t,par->temp_d_step*factor);
  }
  if(par->flag_include_polarization && par->flag_independent_polarization) {
    x_spec[par->index_beta_s_p]=pst->prior_mean[par->index_beta_s_p];
    x_spec[par->index_beta_d_p]=pst->prior_mean[par->index_beta_d_p];
    x_spec[par->index_temp_d_p]=pst->prior_mean[par->index_temp_d_p];
    if(par->n_spec_vary>0) {
      if(par->flag_beta_s_free)
	gsl_matrix_set(mat_step,par->index_beta_s_p,par->index_beta_s_p,par->beta_s_step*factor);
      if(par->flag_beta_d_free)
	gsl_matrix_set(mat_step,par->index_beta_d_p,par->index_beta_d_p,par->beta_d_step*factor);
      if(par->flag_temp_d_free)
	gsl_matrix_set(mat_step,par->index_temp_d_p,par->index_temp_d_p,par->temp_d_step*factor);
    }
  }
}

static void compute_covariance_sampling(ParamBFoRe *par,Rng *rng,flouble *data,flouble *noise_w,
					flouble *x_spec,PixelState *pst_new,PixelState *pst_old,
					gsl_matrix *cov_out,int do_print)
{
  int ic1,ic2,i_sample;
  flouble stepping_factor=1.;
  int nsamples=par->n_spec_vary*par->n_spec_vary*50;
  flouble ratio_accepted=0;
  flouble *x_spec_old=my_malloc(par->n_param_max*sizeof(flouble));
  flouble *x_spec_new=my_malloc(par->n_param_max*sizeof(flouble));
  flouble *mean=my_calloc(par->n_spec_vary,sizeof(flouble));
  gsl_matrix *mat_step=gsl_matrix_alloc(par->n_spec_vary,par->n_spec_vary);
  flouble chi0=pst_old->chi2;

  gsl_matrix_set_zero(cov_out);
  restart_mcmc(par,pst_old,x_spec_old,mat_step,stepping_factor);
  memcpy(x_spec_old,x_spec,par->n_param_max*sizeof(flouble));

  for(i_sample=0;i_sample<nsamples;i_sample++) {
    int accepted=draw_spectral_indices_marginal(par,rng,data,noise_w,x_spec_old,pst_old,mat_step,x_spec_new,pst_new);
    if(accepted) {
      PixelState *tmp=pst_old;
      memcpy(x_spec_old,x_spec_new,par->n_param_max*sizeof(flouble));
      pst_old=pst_new;
      pst_new=tmp;
    }

    for(ic1=0;ic1<par->n_spec_vary;ic1++) {
      mean[ic1]+=x_spec_old[ic1];
      for(ic2=0;ic2<par->n_spec_vary;ic2++) {
	gsl_matrix_set(cov_out,ic1,ic2,gsl_matrix_get(cov_out,ic1,ic2)+
		       x_spec_old[ic1]*x_spec_old[ic2]);
      }
    }

    ratio_accepted+=accepted;
    if(i_sample%500==499) {
      if(ratio_accepted<50) {
	stepping_factor*=0.5;
	dbg_printf(do_print,"Not enough samples, restarting with smaller step size %lf\n",stepping_factor);
	restart_mcmc(par,pst_old,x_spec_old,mat_step,stepping_factor);
	memcpy(x_spec_old,x_spec,par->n_param_max*sizeof(flouble));
	compute_marginalized_chi2(par,data,noise_w,x_spec_old,pst_old,0);
	i_sample=-1;
	gsl_matrix_set_zero(cov_out);
	memset(mean,0,par->n_spec_vary*sizeof(flouble));
      }
      else {
	ratio_accepted/=500.;
	dbg_printf(do_print,"%d Acceptance ratio %.2lf\n",i_sample,ratio_accepted);
	ratio_accepted=0;
      }
    }
  }

  dbg_printf(do_print,"Mean : ");
  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
    mean[ic1]/=nsamples;
    dbg_printf(do_print,"%lf ",mean[ic1]);
  }
  dbg_printf(do_print,"\n");

  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
    for(ic2=0;ic2<par->n_spec_vary;ic2++) {
      flouble cov=gsl_matrix_get(cov_out,ic1,ic2)/nsamples-mean[ic1]*mean[ic2];
      gsl_matrix_set(cov_out,ic1,ic2,cov);
    }
  }
  decouple_covariance(par,cov_out);

  pst_old->chi2=chi0;

  free(mean);
  free(x_spec_old);
  free(x_spec_new);
  gsl_matrix_free(mat_step);
}

static void compute_covariance_numerical(ParamBFoRe *par,flouble *data,flouble *noise_w,
					 flouble *x_spec,PixelState *pst,gsl_matrix *cov_out)
{
  int ic1;
  flouble chi00,chipp,chimm,chipm,chimp,chi_off;
  flouble *diag_save=my_malloc(par->n_spec_vary*sizeof(flouble));
  flouble *stepsize=my_malloc(par->n_spec_vary*sizeof(flouble));
  flouble *x_dum=my_malloc(par->n_param_max*sizeof(flouble));

  for(ic1=0;ic1<par->n_param_max;ic1++)
    printf("%lE ",x_spec[ic1]);
  printf("\n");

  memcpy(x_dum,x_spec,par->n_param_max*sizeof(flouble));
  compute_marginalized_chi2(par,data,noise_w,x_dum,pst,0); chi00=pst->chi2;
  chi_off=chi00; chi00-=chi_off;

  //Compute diagonal first
  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
    flouble h,cov_elem;
    int converged=0;
    char fn[256];

    memcpy(x_dum,x_spec,par->n_param_max*sizeof(flouble));
    h=2.0/pst->prior_isigma[ic1];

    x_dum[ic1]=x_spec[ic1]+h*0.5;
    compute_marginalized_chi2(par,data,noise_w,x_dum,pst,chi_off); chipp=pst->chi2;
    x_dum[ic1]=x_spec[ic1]-h*0.5;
    compute_marginalized_chi2(par,data,noise_w,x_dum,pst,chi_off); chimm=pst->chi2;
    cov_elem=2*(chipp-2*chi00+chimm)/(h*0.5*h*0.5);
    int iter=0;
    while((!converged) && (iter<100)) {
      flouble cov_elem_new;
      h*=0.5;
      iter++;

      memcpy(x_dum,x_spec,par->n_param_max*sizeof(flouble));
      x_dum[ic1]=x_spec[ic1]+h;
      compute_marginalized_chi2(par,data,noise_w,x_dum,pst,chi_off); chipp=pst->chi2;
      x_dum[ic1]=x_spec[ic1]-h;
      compute_marginalized_chi2(par,data,noise_w,x_dum,pst,chi_off); chimm=pst->chi2;
      cov_elem_new=(chipp-2*chi00+chimm)/(h*h);

      if((cov_elem_new>0) && (fabs(cov_elem_new/cov_elem-1)<1E-2))
	converged=1;
      cov_elem=cov_elem_new;
    }

    gsl_matrix_set(cov_out,ic1,ic1,cov_elem);
    stepsize[ic1]=h;
  }

  //Compute off diagonal
  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
    int ic2;
    for(ic2=ic1+1;ic2<par->n_spec_vary;ic2++) {
      flouble cov_elem;

      memcpy(x_dum,x_spec,par->n_param_max*sizeof(flouble));

      x_dum[ic1]=x_spec[ic1]+stepsize[ic1];
      x_dum[ic2]=x_spec[ic2]+stepsize[ic2];
      compute_marginalized_chi2(par,data,noise_w,x_dum,pst,chi_off); chipp=pst->chi2;

      x_dum[ic1]=x_spec[ic1]+stepsize[ic1];
      x_dum[ic2]=x_spec[ic2]-stepsize[ic2];
      compute_marginalized_chi2(par,data,noise_w,x_dum,pst,chi_off); chipm=pst->chi2;

      x_dum[ic1]=x_spec[ic1]-stepsize[ic1];
      x_dum[ic2]=x_spec[ic2]+stepsize[ic2];
      compute_marginalized_chi2(par,data,noise_w,x_dum,pst,chi_off); chimp=pst->chi2;

      x_dum[ic1]=x_spec[ic1]-stepsize[ic1];
      x_dum[ic2]=x_spec[ic2]-stepsize[ic2];
      compute_marginalized_chi2(par,data,noise_w,x_dum,pst,chi_off); chimm=pst->chi2;

      cov_elem=(chipp-chipm-chimp+chimm)/(4*stepsize[ic1]*stepsize[ic2]);
      gsl_matrix_set(cov_out,ic1,ic2,cov_elem);
      gsl_matrix_set(cov_out,ic2,ic1,cov_elem);
    }
  }

  //Divide by two
  decouple_covariance(par,cov_out);
  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
    int ic2;
    for(ic2=0;ic2<par->n_spec_vary;ic2++) {
      gsl_matrix_set(cov_out,ic1,ic2,gsl_matrix_get(cov_out,ic1,ic2)*0.5);
    }
    diag_save[ic1]=gsl_matrix_get(cov_out,ic1,ic1);
    if(diag_save[ic1]<=0)
      report_error(1,"Bad covariance\n");
  }

  //Invert
  int err=gsl_linalg_cholesky_decomp(cov_out);
  if(err) {
    gsl_matrix_set_zero(cov_out);
    for(ic1=0;ic1<par->n_spec_vary;ic1++)
      gsl_matrix_set(cov_out,ic1,ic1,1./diag_save[ic1]);
  }
  else
    gsl_linalg_cholesky_invert(cov_out);
  free(diag_save);

  pst->chi2=chi_off;
  free(x_dum); free(stepsize);
}

static int draw_spectral_indices_marginal(ParamBFoRe *par,Rng *rng,flouble *data,flouble *noise_w,
					  flouble *x_spec_old,PixelState *pst_old,gsl_matrix *mat_step,
					  flouble *x_spec_new,PixelState *pst_new)
{ //DAM: possible optimization: demote mat_step to flouble *
  int ipar;
  double ratio;
  memcpy(x_spec_new,x_spec_old,par->n_param_max*sizeof(flouble));
  for(ipar=0;ipar<par->n_spec_vary;ipar++)
    pst_new->rand_spec[ipar]=rand_gauss(rng);
  for(ipar=0;ipar<par->n_spec_vary;ipar++) {
    int ipar2;
    for(ipar2=0;ipar2<=ipar;ipar2++)
      x_spec_new[ipar]+=gsl_matrix_get(mat_step,ipar,ipar2)*pst_new->rand_spec[ipar2];
  }
  compute_marginalized_chi2(par,data,noise_w,x_spec_new,pst_new,0);

  ratio=exp(-0.5*(pst_new->chi2-pst_old->chi2));

  if(ratio<1) {
    if(rand_real01(rng)>ratio)
      return 0;
  }

  return 1;
}

static void close_region_file(FILE *fo,int ipix)
{
  my_fwrite(&ipix,sizeof(ipix),1,fo);
  fclose(fo);
}

static FILE *open_region_file(ParamBFoRe *par,int ipix)
{
  FILE *fo;
  int size_float=sizeof(flouble);
  char fname[256];
  sprintf(fname,"%s_pix%d.dat",par->output_prefix,ipix);
  fo=my_fopen(fname,"wb");

  my_fwrite(&ipix,sizeof(ipix),1,fo);
  my_fwrite(&NodeThis,sizeof(NodeThis),1,fo);
  my_fwrite(&size_float,sizeof(size_float),1,fo);
  my_fwrite(&(par->nside),sizeof(par->nside),1,fo);
  my_fwrite(&(par->nside_spec),sizeof(par->nside_spec),1,fo);
  my_fwrite(&(par->n_sub),sizeof(par->n_sub),1,fo);
  my_fwrite(&(par->n_nu),sizeof(par->n_nu),1,fo);
  my_fwrite(&(par->n_pol),sizeof(par->n_pol),1,fo);
  my_fwrite(&(par->n_comp),sizeof(par->n_comp),1,fo);
  my_fwrite(&(par->n_spec_vary),sizeof(par->n_spec_vary),1,fo);
  my_fwrite(&(par->n_samples),sizeof(par->n_samples),1,fo);
  my_fwrite(&(par->n_samples_burn),sizeof(par->n_samples_burn),1,fo);
  my_fwrite(&(par->n_output_rate),sizeof(par->n_output_rate),1,fo);

  return fo;
}

typedef struct {
  ParamBFoRe *par;
  PixelState *pst;
  flouble *data;
  flouble *noise_w;
  flouble *x_spec;
  flouble chi0;
  int ipar;
  int ider;
} ParamChi2;

static flouble chi2_marg_func(flouble *x,void *pars)
{
  ParamChi2 *p=(ParamChi2 *)pars;
  memcpy(p->x_spec,x,p->par->n_spec_vary*sizeof(flouble));
  compute_marginalized_chi2(p->par,p->data,p->noise_w,p->x_spec,p->pst,0);
  return p->pst->chi2;
}

static flouble chi2_marg_func_1p(flouble x,void *pars)
{
  flouble chi2_out;
  ParamChi2 *p=(ParamChi2 *)pars;
  flouble x_save=p->x_spec[p->ipar];
  flouble chi2_save=p->pst->chi2;

  p->x_spec[p->ipar]=x;
  compute_marginalized_chi2(p->par,p->data,p->noise_w,p->x_spec,p->pst,p->chi0);

  p->x_spec[p->ipar]=x_save;
  chi2_out=p->pst->chi2;
  p->pst->chi2=chi2_save;
  return chi2_out;
}

static flouble dchi2_marg_func_1p(flouble x,void *pars)
{
  ParamChi2 *p=(ParamChi2 *)pars;
  ParamChi2 p2;
  flouble x_save=p->x_spec[p->ipar];
  flouble chi2_save=p->pst->chi2;
  flouble der_out,err_out;

  p->x_spec[p->ipar]=x;
  
  p2.par=p->par;
  p2.pst=p->pst;
  p2.data=p->data;
  p2.noise_w=p->noise_w;
  p2.x_spec=p->x_spec;
  p2.chi0=p->chi0;
  p2.ipar=p->ider;

  gsl_function F;
  F.function=&chi2_marg_func_1p;
  F.params=&p2;

  gsl_deriv_central(&F,p->x_spec[p->ider],1./p->pst->prior_isigma[p->ider],&der_out,&err_out);

  p->x_spec[p->ipar]=x_save;
  p->pst->chi2=chi2_save;
  return der_out;
}

static void get_ml_marginal(ParamBFoRe *par,flouble *data,flouble *noise_w,
			    flouble *x_spec,PixelState *pst)
{
  int ii;
  ParamChi2 pc2;
  PowellParams *par_pow;

  pc2.par=par;
  pc2.data=data;
  pc2.noise_w=noise_w;
  pc2.pst=pst;
  pc2.x_spec=x_spec;

  par_pow=powell_params_new(par->n_spec_vary,x_spec,&chi2_marg_func,&pc2,100,1E-7);
  powell(par_pow);

  for(ii=0;ii<par->n_spec_vary;ii++)
    x_spec[ii]=par_pow->p[ii];
  free_powell_params(par_pow);

  compute_marginalized_chi2(par,data,noise_w,x_spec,pst,0);
}

static void compute_covariance_numerical_b(ParamBFoRe *par,flouble *data,flouble *noise_w,
					   flouble *x_spec,PixelState *pst,gsl_matrix *cov_out)
{
  int ic1;
  flouble chi00;
  ParamChi2 pc2;
  flouble *diag_save=my_malloc(par->n_spec_vary*sizeof(flouble));
  compute_marginalized_chi2(par,data,noise_w,x_spec,pst,0); chi00=pst->chi2;
  pc2.par=par;
  pc2.data=data;
  pc2.noise_w=noise_w;
  pc2.pst=pst;
  pc2.x_spec=x_spec;
  pc2.chi0=chi00;

  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
    int ic2;
    for(ic2=ic1;ic2<par->n_spec_vary;ic2++) {
      flouble result,error;
      pc2.ipar=ic2;
      pc2.ider=ic1;

      gsl_function F;
      F.function=&dchi2_marg_func_1p;
      F.params=&pc2;
      gsl_deriv_central(&F,x_spec[ic2],1./pst->prior_isigma[ic2],&result,&error);
      gsl_matrix_set(cov_out,ic1,ic2,result);
      if(ic1!=ic2)
	gsl_matrix_set(cov_out,ic2,ic1,result);
    }
  }

  decouple_covariance(par,cov_out);
  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
    int ic2;
    for(ic2=0;ic2<par->n_spec_vary;ic2++) {
      gsl_matrix_set(cov_out,ic1,ic2,gsl_matrix_get(cov_out,ic1,ic2)*0.5);
    }
    diag_save[ic1]=gsl_matrix_get(cov_out,ic1,ic1);
    if(diag_save[ic1]<=0)
      report_error(1,"Bad covariance\n");
  }

  //Invert
  int err=gsl_linalg_cholesky_decomp(cov_out);
  if(err) {
    gsl_matrix_set_zero(cov_out);
    for(ic1=0;ic1<par->n_spec_vary;ic1++)
      gsl_matrix_set(cov_out,ic1,ic1,1./diag_save[ic1]);
  }
  else
    gsl_linalg_cholesky_invert(cov_out);
  free(diag_save);
}


void clean_pixel(ParamBFoRe *par,Rng *rng,PixelState *pst,int ipix_big)
{
  int i_sample,ic1,ic2,ipix,n_updated,err,accepted;
  FILE *fo;
  flouble ratio_accepted=0;
  flouble stepping_factor=1.;
  int ip_spc=par->n_spec_vary*ipix_big;
  int id_cell=ipix_big*par->n_sub*par->n_pol;
  flouble *data=&(par->maps_data[id_cell*par->n_nu]);
  flouble *noise_w=&(par->maps_noise_weight[id_cell*par->n_nu]);
  flouble *amps_mean=&(par->map_components_mean[id_cell*par->n_comp]);
  flouble *amps_covar=&(par->map_components_covar[id_cell*par->n_comp*par->n_comp]);
  flouble *amps_dum=my_malloc(par->n_sub*par->n_pol*par->n_comp*sizeof(flouble));
  flouble *x_spec_old=my_malloc(par->n_param_max*sizeof(flouble));
  flouble *x_spec_ml=my_malloc(par->n_param_max*sizeof(flouble));
  flouble *x_spec_new=my_malloc(par->n_param_max*sizeof(flouble));
  gsl_matrix *mat_step=gsl_matrix_alloc(par->n_spec_vary,par->n_spec_vary);
  gsl_matrix *cov_save=gsl_matrix_alloc(par->n_spec_vary,par->n_spec_vary);
  gsl_matrix *cov_spec=gsl_matrix_alloc(par->n_spec_vary,par->n_spec_vary);
  flouble *mean_spec=my_calloc(par->n_spec_vary,sizeof(flouble));
  flouble factor_rescale=2.4/sqrt((double)(par->n_spec_vary));
  int do_print=(ipix_big==par->dbg_ipix);

  if(par->flag_write_samples)
    fo=open_region_file(par,ipix_big);
  else
    fo=NULL;

  init_priors(par,pst,ipix_big);
  memset(amps_mean,0,par->n_sub*par->n_pol*par->n_comp*sizeof(flouble));
  memset(amps_covar,0,par->n_sub*par->n_pol*par->n_comp*par->n_comp*sizeof(flouble));
  restart_mcmc(par,pst,x_spec_old,mat_step,stepping_factor);
  get_ml_marginal(par,data,noise_w,x_spec_old,pst);
  memcpy(x_spec_ml,x_spec_old,par->n_param_max*sizeof(flouble));

  if(par->n_spec_vary>0) {
    dbg_printf(do_print,"Burning\n");
    ratio_accepted=0;
    gsl_matrix_set_zero(cov_spec);
    gsl_matrix_set_zero(cov_save);
    n_updated=0;
    for(i_sample=0;i_sample<par->n_samples_burn;i_sample++) {
      if(i_sample%par->n_update_covar==0)
	draw_amplitudes(par,rng,data,noise_w,x_spec_old,pst,amps_dum); //A_{n+1}(b_n)
      accepted=draw_spectral_indices(par,rng,data,noise_w,amps_dum,
				     x_spec_old,pst,mat_step,x_spec_new); //b_{n+1}(A_{n+1})
      if(accepted)
	memcpy(x_spec_old,x_spec_new,par->n_param_max*sizeof(flouble));
#ifdef _DEBUG
      if(ipix_big==par->dbg_ipix)
	memcpy(&(par->dbg_extra[i_sample*par->n_spec_vary]),x_spec_old,par->n_spec_vary*sizeof(flouble));
#endif //_DEBUG
      for(ic1=0;ic1<par->n_spec_vary;ic1++) {
	mean_spec[ic1]+=x_spec_old[ic1];
	for(ic2=0;ic2<par->n_spec_vary;ic2++) {
	  flouble cov_plus=x_spec_old[ic1]*x_spec_old[ic2];
	  gsl_matrix_set(cov_spec,ic1,ic2,gsl_matrix_get(cov_spec,ic1,ic2)+cov_plus);
	}
      }
      
      ratio_accepted+=accepted;
      if(i_sample%par->n_update_covar==par->n_update_covar-1) { //Update covariance
	if(ratio_accepted<=5) { //Check if too few samples were accepted
	  dbg_printf(do_print,"No samples were accepted! Restarting with smaller step size\n");
	  stepping_factor*=0.5;
	  restart_mcmc(par,pst,x_spec_old,mat_step,stepping_factor);
	  memcpy(x_spec_old,x_spec_ml,par->n_param_max*sizeof(flouble));
	  gsl_matrix_set_zero(cov_save);
	  n_updated=0;
	  i_sample=-1;
	}
	else {
	  //Compute mean in this batch
	  ratio_accepted/=par->n_update_covar;
	  dbg_printf(do_print,"%d Acceptance ratio %.2lf\n",i_sample,ratio_accepted);
	  
	  dbg_printf(do_print,"Current mean :");
	  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
	    mean_spec[ic1]/=par->n_update_covar;
	    dbg_printf(do_print," %lf",mean_spec[ic1]);
	  }
	  dbg_printf(do_print,"\n");
	  //Compute covariance in this batch
	  dbg_printf(do_print,"Current covariance :\n");
	  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
	    dbg_printf(do_print,"   |");
	    for(ic2=0;ic2<par->n_spec_vary;ic2++) {
	      flouble cov=gsl_matrix_get(cov_spec,ic1,ic2)/par->n_update_covar-
		mean_spec[ic1]*mean_spec[ic2];
	      gsl_matrix_set(cov_spec,ic1,ic2,cov);
	      dbg_printf(do_print," %lE",gsl_matrix_get(cov_spec,ic1,ic2));
	    }
	    dbg_printf(do_print,"|\n");
	  }
	  //Add to mean covariance
	  if(i_sample>par->n_samples_burn/2) {
	    n_updated++;
	    for(ic1=0;ic1<par->n_spec_vary;ic1++) {
	      for(ic2=0;ic2<par->n_spec_vary;ic2++)
		gsl_matrix_set(cov_save,ic1,ic2,gsl_matrix_get(cov_save,ic1,ic2)+
			       gsl_matrix_get(cov_spec,ic1,ic2));
	    }
	  }
	  
	  //Save diagonal
	  for(ic1=0;ic1<par->n_spec_vary;ic1++)
	    mean_spec[ic1]=sqrt(gsl_matrix_get(cov_spec,ic1,ic1));
	  //Cholesky decomposition of covariance
	  err=gsl_linalg_cholesky_decomp(cov_spec);
	  gsl_matrix_set_zero(mat_step);
	  if(err==GSL_EDOM) { //If covariance is not positive definite just save standard deviations
	    dbg_printf(do_print,"Covariance is not positive definite\n");
	    for(ic1=0;ic1<par->n_spec_vary;ic1++)
	      gsl_matrix_set(mat_step,ic1,ic1,factor_rescale*mean_spec[ic1]);
	  }
	  else { //Store cholesky decomposition in stepping function
	    for(ic1=0;ic1<par->n_spec_vary;ic1++) {
	      for(ic2=0;ic2<=ic1;ic2++)
		gsl_matrix_set(mat_step,ic1,ic2,factor_rescale*gsl_matrix_get(cov_spec,ic1,ic2));
	    }
	  }
	}
	
	gsl_matrix_set_zero(cov_spec);
	for(ic1=0;ic1<par->n_spec_vary;ic1++)
	  mean_spec[ic1]=0;
	ratio_accepted=0;
      }
    }
    dbg_printf(do_print,"Final covariance :\n");
    for(ic1=0;ic1<par->n_spec_vary;ic1++) { //Compute covariance
      dbg_printf(do_print,"   |");
      for(ic2=0;ic2<par->n_spec_vary;ic2++) {
	flouble cov=gsl_matrix_get(cov_save,ic1,ic2)/n_updated;
	gsl_matrix_set(cov_save,ic1,ic2,cov);
	dbg_printf(do_print," %lE",cov);
      }
      dbg_printf(do_print,"|\n");
      mean_spec[ic1]=sqrt(gsl_matrix_get(cov_save,ic1,ic1)); //Save diagonal
    }
    err=gsl_linalg_cholesky_decomp(cov_save); //Cholesky decomposition of covariance
    gsl_matrix_set_zero(mat_step);
    if(err==GSL_EDOM) { //If covariance is not positive definite just save standard deviations
      dbg_printf(do_print,"Final covariance is not positive definite!!!!\n");
      for(ic1=0;ic1<par->n_spec_vary;ic1++)
	gsl_matrix_set(mat_step,ic1,ic1,factor_rescale*mean_spec[ic1]);
    }
    else { //Store cholesky decomposition in stepping function
      for(ic1=0;ic1<par->n_spec_vary;ic1++) {
	for(ic2=0;ic2<=ic1;ic2++)
	  gsl_matrix_set(mat_step,ic1,ic2,factor_rescale*gsl_matrix_get(cov_save,ic1,ic2));
      }
    }
  }

  dbg_printf(do_print,"Starting actual sampling\n");
#ifdef _DEBUG
  if(ipix_big==par->dbg_ipix)
    ratio_accepted=0;
#endif //_DEBUG
  for(i_sample=0;i_sample<par->n_samples;i_sample++) {
    if(i_sample%par->n_spec_resample==0)
      draw_amplitudes(par,rng,data,noise_w,x_spec_old,pst,amps_dum); //A_{n+1}(b_n)
    if(par->n_spec_vary>0) {
      accepted=draw_spectral_indices(par,rng,data,noise_w,amps_dum,
				     x_spec_old,pst,mat_step,x_spec_new); //b_{n+1}(A_{n+1})
      if(accepted)
	memcpy(x_spec_old,x_spec_new,par->n_param_max*sizeof(flouble));
    }
    if((par->flag_write_samples) && (i_sample%par->n_output_rate==0)) {
      int i_sample_write=i_sample/par->n_output_rate;
      my_fwrite(&i_sample_write,sizeof(int),1,fo);
      my_fwrite(x_spec_old,sizeof(flouble),par->n_spec_vary,fo);
      my_fwrite(amps_dum,sizeof(flouble),par->n_sub*par->n_pol*par->n_comp,fo);
      my_fwrite(&i_sample_write,sizeof(int),1,fo);
    }

#ifdef _DEBUG
    if(par->n_spec_vary>0) {
      if(ipix_big==par->dbg_ipix) {
	memcpy(&(par->dbg_extra[(i_sample+par->n_samples_burn)*par->n_spec_vary]),
	       x_spec_old,par->n_spec_vary*sizeof(flouble));
	ratio_accepted+=accepted;
	if(i_sample && (i_sample%par->n_update_covar==par->n_update_covar-1)) {
	  ratio_accepted/=par->n_update_covar;
	  dbg_printf(do_print,"%d Acceptance ratio %.2lf\n",i_sample,ratio_accepted);
	  ratio_accepted=0;
	}
      }
    }
#endif //_DEBUG

    //Add amplitudes
    for(ipix=0;ipix<par->n_sub*par->n_pol;ipix++) {
      int index_here=ipix*par->n_comp;
      for(ic1=0;ic1<par->n_comp;ic1++) {
	amps_mean[ic1+index_here]+=amps_dum[ic1+index_here];
	for(ic2=0;ic2<par->n_comp;ic2++)
	  amps_covar[ic2+par->n_comp*(ic1+index_here)]+=amps_dum[ic1+index_here]*amps_dum[ic2+index_here];
      }
    }
    //Add indices
    for(ic1=0;ic1<par->n_spec_vary;ic1++) {
      par->map_indices_mean[ic1+ip_spc]+=x_spec_old[ic1];
      for(ic2=0;ic2<par->n_spec_vary;ic2++)
	par->map_indices_covar[ic1+par->n_spec_vary*(ic2+ip_spc)]+=x_spec_old[ic1]*x_spec_old[ic2];
    }
  }

  //Compute averages
  for(ipix=0;ipix<par->n_sub*par->n_pol;ipix++) {
    int index_here=ipix*par->n_comp;
    for(ic1=0;ic1<par->n_comp;ic1++)
      amps_mean[ic1+index_here]/=par->n_samples;
  }
  for(ipix=0;ipix<par->n_sub*par->n_pol;ipix++) {
    int index_here=ipix*par->n_comp;
    for(ic1=0;ic1<par->n_comp;ic1++) {
      for(ic2=0;ic2<par->n_comp;ic2++) {
	amps_covar[ic2+par->n_comp*(ic1+index_here)]/=par->n_samples;
	amps_covar[ic2+par->n_comp*(ic1+index_here)]-=amps_mean[ic1+index_here]*amps_mean[ic2+index_here];
      }
    }
  }
  for(ic1=0;ic1<par->n_spec_vary;ic1++)
    par->map_indices_mean[ic1+ip_spc]/=par->n_samples;
  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
    for(ic2=0;ic2<par->n_spec_vary;ic2++) {
      par->map_indices_covar[ic1+par->n_spec_vary*(ic2+ip_spc)]/=par->n_samples;
      par->map_indices_covar[ic1+par->n_spec_vary*(ic2+ip_spc)]-=
	par->map_indices_mean[ic1+ip_spc]*par->map_indices_mean[ic2+ip_spc];
    }
  }

  if(par->flag_write_samples)
    close_region_file(fo,ipix_big);
  free(amps_dum);
  free(x_spec_old);
  free(x_spec_new);
  free(mean_spec);
  gsl_matrix_free(mat_step);
  gsl_matrix_free(cov_save);
  gsl_matrix_free(cov_spec);
}

#define COV_TYPE_NUMER_GSL 0
#define COV_TYPE_NUMER_BFR 1
#define COV_TYPE_SAMPLING 2

static void compute_covariance_wrap(ParamBFoRe *par,Rng *rng,flouble *data,flouble *noise_w,
				    flouble *x_spec,PixelState *pst_new,PixelState *pst_old,
				    gsl_matrix *cov_out,int do_print,int cov_type)
{
  if(cov_type==COV_TYPE_NUMER_GSL)
    compute_covariance_numerical_b(par,data,noise_w,x_spec,pst_old,cov_out);
  else if(cov_type==COV_TYPE_NUMER_BFR)
    compute_covariance_numerical(par,data,noise_w,x_spec,pst_old,cov_out);
  else if(cov_type==COV_TYPE_SAMPLING)
    compute_covariance_sampling(par,rng,data,noise_w,x_spec,pst_new,pst_old,cov_out,do_print);
  else 
    compute_covariance_numerical_b(par,data,noise_w,x_spec,pst_old,cov_out);
}

void clean_pixel_from_marginal(ParamBFoRe *par,Rng *rng,PixelState *pst_old,
			       PixelState *pst_new,int ipix_big)
{
  int i_sample,ic1,ic2,ipix,n_updated,err,accepted;
  flouble ratio_accepted;
  int cov_type=COV_TYPE_NUMER_BFR;
  int ip_spc=par->n_spec_vary*ipix_big;
  int id_cell=ipix_big*par->n_sub*par->n_pol;
  flouble *data=&(par->maps_data[id_cell*par->n_nu]);
  flouble *noise_w=&(par->maps_noise_weight[id_cell*par->n_nu]);
  flouble *amps_mean=&(par->map_components_mean[id_cell*par->n_comp]);
  flouble *amps_covar=&(par->map_components_covar[id_cell*par->n_comp*par->n_comp]);
  flouble *x_spec_old=my_malloc(par->n_param_max*sizeof(flouble));
  flouble *x_spec_mean=my_calloc(par->n_param_max,sizeof(flouble));
  flouble *x_spec_new=my_malloc(par->n_param_max*sizeof(flouble));
  gsl_matrix *mat_step=gsl_matrix_alloc(par->n_spec_vary,par->n_spec_vary);
  gsl_matrix *cov_save=gsl_matrix_alloc(par->n_spec_vary,par->n_spec_vary);
  gsl_matrix *cov_spec=gsl_matrix_alloc(par->n_spec_vary,par->n_spec_vary);
  flouble *mean_spec=my_calloc(par->n_spec_vary,sizeof(flouble));
  flouble factor_rescale=2.4/sqrt((double)(par->n_spec_vary));
  int do_print=(ipix_big==par->dbg_ipix);

  init_priors(par,pst_old,ipix_big);
  init_priors(par,pst_new,ipix_big);
  memset(amps_mean,0,par->n_sub*par->n_pol*par->n_comp*sizeof(flouble));
  memset(amps_covar,0,par->n_sub*par->n_pol*par->n_comp*par->n_comp*sizeof(flouble));

  //Compute ML point
  restart_mcmc(par,pst_old,x_spec_old,mat_step,1.);
  get_ml_marginal(par,data,noise_w,x_spec_old,pst_old);
  memcpy(x_spec_mean,x_spec_old,par->n_param_max*sizeof(flouble));

  dbg_printf(do_print,"Mean: ");
  for(ic1=0;ic1<par->n_spec_vary;ic1++)
    dbg_printf(do_print,"%lf ",x_spec_old[ic1]);
  dbg_printf(do_print,"\n");

  //Compute Covariance around ML
  dbg_printf(do_print,"Computing numerical covariance\n");
  compute_covariance_wrap(par,rng,data,noise_w,x_spec_old,pst_new,pst_old,cov_spec,do_print,cov_type);

  dbg_printf(do_print,"Numerical covariance :\n");
  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
    dbg_printf(do_print,"   |");
    for(ic2=0;ic2<par->n_spec_vary;ic2++)
      dbg_printf(do_print," %lE",gsl_matrix_get(cov_spec,ic1,ic2));
    dbg_printf(do_print,"|\n");
  }

  err=0;
  //Save diagonal
  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
    flouble cov=gsl_matrix_get(cov_spec,ic1,ic1);
    if(cov<=0)
      err=GSL_EDOM;
    mean_spec[ic1]=sqrt(cov);
  }
  if(err==GSL_EDOM) {
    dbg_printf(do_print,"Something is wrong iwth the covariance matrix\n");
    report_error(1,"Exiting\n");
  }

  //Cholesky decomposition of covariance
  err=gsl_linalg_cholesky_decomp(cov_spec);
  gsl_matrix_set_zero(mat_step);
  if(err==GSL_EDOM) { //If covariance is not positive definite just save standard deviations
    dbg_printf(do_print,"Covariance is not positive definite\n");
    for(ic1=0;ic1<par->n_spec_vary;ic1++)
      gsl_matrix_set(mat_step,ic1,ic1,factor_rescale*mean_spec[ic1]);
  }
  else { //Store cholesky decomposition in stepping function
    for(ic1=0;ic1<par->n_spec_vary;ic1++) {
      for(ic2=0;ic2<=ic1;ic2++)
	gsl_matrix_set(mat_step,ic1,ic2,factor_rescale*gsl_matrix_get(cov_spec,ic1,ic2));
    }
  }

  if(par->n_spec_vary>0) {
    dbg_printf(do_print,"Burning\n");

    ratio_accepted=0;
    for(ic1=0;ic1<par->n_spec_vary;ic1++)
      x_spec_mean[ic1]=0;

    for(i_sample=0;i_sample<par->n_samples_burn;i_sample++) {
      accepted=draw_spectral_indices_marginal(par,rng,data,noise_w,x_spec_old,pst_old,
					      mat_step,x_spec_new,pst_new); //b_{n+1}(A_{n+1})
      if(accepted) {
	PixelState *tmp=pst_old;
	memcpy(x_spec_old,x_spec_new,par->n_param_max*sizeof(flouble));
	pst_old=pst_new;
	pst_new=tmp;
      }

      for(ic1=0;ic1<par->n_spec_vary;ic1++)
	x_spec_mean[ic1]+=x_spec_old[ic1];

#ifdef _DEBUG
      if(ipix_big==par->dbg_ipix)
	memcpy(&(par->dbg_extra[i_sample*par->n_spec_vary]),x_spec_old,par->n_spec_vary*sizeof(flouble));
#endif //_DEBUG

      ratio_accepted+=accepted;
      if(i_sample%par->n_update_covar==par->n_update_covar-1) { //Update covariance
	//Compute mean in this batch
	ratio_accepted/=par->n_update_covar;
	dbg_printf(do_print,"%d Acceptance ratio %.2lf\n",i_sample,ratio_accepted);
	ratio_accepted=0;
      }
    }

    for(ic1=0;ic1<par->n_spec_vary;ic1++)
      x_spec_mean[ic1]/=par->n_samples_burn;
  }

  dbg_printf(do_print,"Current: ");
  for(ic1=0;ic1<par->n_spec_vary;ic1++)
    dbg_printf(do_print,"%lf ",x_spec_old[ic1]);
  dbg_printf(do_print,"\n");

  dbg_printf(do_print,"Mean: ");
  for(ic1=0;ic1<par->n_spec_vary;ic1++)
    dbg_printf(do_print,"%lf ",x_spec_mean[ic1]);
  dbg_printf(do_print,"\n");

  //Compute Covariance around ML
  dbg_printf(do_print,"Computing numerical covariance\n");
  compute_covariance_wrap(par,rng,data,noise_w,x_spec_old,pst_new,pst_old,cov_spec,do_print,cov_type);

  dbg_printf(do_print,"Numerical covariance :\n");
  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
    dbg_printf(do_print,"   |");
    for(ic2=0;ic2<par->n_spec_vary;ic2++)
      dbg_printf(do_print," %lE",gsl_matrix_get(cov_spec,ic1,ic2));
    dbg_printf(do_print,"|\n");
  }

  err=0;
  //Save diagonal
  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
    flouble cov=gsl_matrix_get(cov_spec,ic1,ic1);
    if(cov<=0)
      err=GSL_EDOM;
    mean_spec[ic1]=sqrt(cov);
  }
  if(err==GSL_EDOM) {
    dbg_printf(do_print,"Something is wrong iwth the covariance matrix\n");
    report_error(1,"Exiting\n");
  }

  //Cholesky decomposition of covariance
  err=gsl_linalg_cholesky_decomp(cov_spec);
  gsl_matrix_set_zero(mat_step);
  if(err==GSL_EDOM) { //If covariance is not positive definite just save standard deviations
    dbg_printf(do_print,"Covariance is not positive definite\n");
    for(ic1=0;ic1<par->n_spec_vary;ic1++)
      gsl_matrix_set(mat_step,ic1,ic1,factor_rescale*mean_spec[ic1]);
  }
  else { //Store cholesky decomposition in stepping function
    for(ic1=0;ic1<par->n_spec_vary;ic1++) {
      for(ic2=0;ic2<=ic1;ic2++)
	gsl_matrix_set(mat_step,ic1,ic2,factor_rescale*gsl_matrix_get(cov_spec,ic1,ic2));
    }
  }

  dbg_printf(do_print,"Starting actual sampling\n");
#ifdef _DEBUG
  if(ipix_big==par->dbg_ipix)
    ratio_accepted=0;
#endif //_DEBUG
  int n_samples_real=par->n_samples;
  if(par->n_spec_vary==0)
    n_samples_real=1;
  for(i_sample=0;i_sample<n_samples_real;i_sample++) {
    if(par->n_spec_vary>0) {
      accepted=draw_spectral_indices_marginal(par,rng,data,noise_w,x_spec_old,pst_old,
					      mat_step,x_spec_new,pst_new); //b_{n+1}(A_{n+1})
      if(accepted) {
	PixelState *tmp=pst_old;
	memcpy(x_spec_old,x_spec_new,par->n_param_max*sizeof(flouble));
	pst_old=pst_new;
	pst_new=tmp;
      }
    }

#ifdef _DEBUG
    if(par->n_spec_vary>0) {
      if(ipix_big==par->dbg_ipix) {
	memcpy(&(par->dbg_extra[(i_sample+par->n_samples_burn)*par->n_spec_vary]),
	       x_spec_old,par->n_spec_vary*sizeof(flouble));
	ratio_accepted+=accepted;
	if(i_sample && (i_sample%par->n_update_covar==par->n_update_covar-1)) {
	  ratio_accepted/=par->n_update_covar;
	  dbg_printf(do_print,"%d Acceptance ratio %.2lf\n",i_sample,ratio_accepted);
	  ratio_accepted=0;
	}
      }
    }
#endif //_DEBUG

    //Add amplitudes
    for(ipix=0;ipix<par->n_sub*par->n_pol;ipix++) {
      int index_here=ipix*par->n_comp;
      for(ic1=0;ic1<par->n_comp;ic1++) {
	amps_mean[ic1+index_here]+=gsl_vector_get(pst_old->vec_mean[ipix],ic1);
	for(ic2=0;ic2<par->n_comp;ic2++)
	  amps_covar[ic2+par->n_comp*(ic1+index_here)]+=gsl_matrix_get(pst_old->cov_inv[ipix],ic1,ic2);
      }
    }
    //Add indices
    for(ic1=0;ic1<par->n_spec_vary;ic1++) {
      par->map_indices_mean[ic1+ip_spc]+=x_spec_old[ic1];
      for(ic2=0;ic2<par->n_spec_vary;ic2++)
	par->map_indices_covar[ic1+par->n_spec_vary*(ic2+ip_spc)]+=x_spec_old[ic1]*x_spec_old[ic2];
    }
  }

  //Compute averages
  for(ipix=0;ipix<par->n_sub*par->n_pol;ipix++) {
    int index_here=ipix*par->n_comp;
    for(ic1=0;ic1<par->n_comp;ic1++) {
      amps_mean[ic1+index_here]/=n_samples_real;
      for(ic2=0;ic2<par->n_comp;ic2++)
	amps_covar[ic2+par->n_comp*(ic1+index_here)]/=n_samples_real;
    }
  }
  for(ic1=0;ic1<par->n_spec_vary;ic1++) 
    par->map_indices_mean[ic1+ip_spc]/=n_samples_real;

  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
    for(ic2=0;ic2<par->n_spec_vary;ic2++) {
      par->map_indices_covar[ic1+par->n_spec_vary*(ic2+ip_spc)]/=n_samples_real;
      par->map_indices_covar[ic1+par->n_spec_vary*(ic2+ip_spc)]-=
	par->map_indices_mean[ic1+ip_spc]*par->map_indices_mean[ic2+ip_spc];
    }
  }

  dbg_printf(do_print,"Mean: ");
  for(ic1=0;ic1<par->n_spec_vary;ic1++)
    dbg_printf(do_print,"%lf ",par->map_indices_mean[ic1+ip_spc]);
  dbg_printf(do_print,"\n");

  dbg_printf(do_print,"Output covariance :\n");
  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
    dbg_printf(do_print,"   |");
    for(ic2=0;ic2<par->n_spec_vary;ic2++)
      dbg_printf(do_print," %lE",par->map_indices_covar[ic1+par->n_spec_vary*(ic2+ip_spc)]);
    dbg_printf(do_print,"|\n");
  }

  free(x_spec_old);
  free(x_spec_new);
  free(x_spec_mean);
  free(mean_spec);
  gsl_matrix_free(mat_step);
  gsl_matrix_free(cov_save);
  gsl_matrix_free(cov_spec);
}

/*
void clean_pixel_from_marginal(ParamBFoRe *par,Rng *rng,PixelState *pst_old,
			       PixelState *pst_new,int ipix_big)
{
  int i_sample,ic1,ic2,ipix,n_updated,err,accepted;
  flouble ratio_accepted;
  flouble stepping_factor=1.;
  int ip_spc=par->n_spec_vary*ipix_big;
  int id_cell=ipix_big*par->n_sub*par->n_pol;
  flouble *data=&(par->maps_data[id_cell*par->n_nu]);
  flouble *noise_w=&(par->maps_noise_weight[id_cell*par->n_nu]);
  flouble *amps_mean=&(par->map_components_mean[id_cell*par->n_comp]);
  flouble *amps_covar=&(par->map_components_covar[id_cell*par->n_comp*par->n_comp]);
  flouble *x_spec_old=my_malloc(par->n_param_max*sizeof(flouble));
  flouble *x_spec_ml=my_malloc(par->n_param_max*sizeof(flouble));
  flouble *x_spec_new=my_malloc(par->n_param_max*sizeof(flouble));
  gsl_matrix *mat_step=gsl_matrix_alloc(par->n_spec_vary,par->n_spec_vary);
  gsl_matrix *cov_save=gsl_matrix_alloc(par->n_spec_vary,par->n_spec_vary);
  gsl_matrix *cov_spec=gsl_matrix_alloc(par->n_spec_vary,par->n_spec_vary);
  gsl_matrix *cov_numerical=gsl_matrix_alloc(par->n_spec_vary,par->n_spec_vary);
  gsl_matrix *cov_numerical_b=gsl_matrix_alloc(par->n_spec_vary,par->n_spec_vary);
  gsl_matrix *cov_sampling=gsl_matrix_alloc(par->n_spec_vary,par->n_spec_vary);
  flouble *mean_spec=my_calloc(par->n_spec_vary,sizeof(flouble));
  flouble factor_rescale=2.4/sqrt((double)(par->n_spec_vary));
  int do_print=(ipix_big==par->dbg_ipix);

  init_priors(par,pst_old,ipix_big);
  init_priors(par,pst_new,ipix_big);
  memset(amps_mean,0,par->n_sub*par->n_pol*par->n_comp*sizeof(flouble));
  memset(amps_covar,0,par->n_sub*par->n_pol*par->n_comp*par->n_comp*sizeof(flouble));
  restart_mcmc(par,pst_old,x_spec_old,mat_step,stepping_factor);
  get_ml_marginal(par,data,noise_w,x_spec_old,pst_old);
  memcpy(x_spec_ml,x_spec_old,par->n_param_max*sizeof(flouble));

  for(ic1=0;ic1<par->n_param_max;ic1++)
    printf("%lE ",x_spec_old[ic1]);
  printf("\n");

  dbg_printf(do_print,"Computing numerical covariance\n");
  compute_covariance_numerical(par,data,noise_w,x_spec_old,pst_old,cov_numerical);
  dbg_printf(do_print,"Numerical covariance :\n");
  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
    dbg_printf(do_print,"   |");
    for(ic2=0;ic2<par->n_spec_vary;ic2++)
      dbg_printf(do_print," %lE",gsl_matrix_get(cov_numerical,ic1,ic2));
    dbg_printf(do_print,"|\n");
  }
  compute_covariance_numerical_b(par,data,noise_w,x_spec_old,pst_old,cov_numerical_b);
  dbg_printf(do_print,"Numerical covariance :\n");
  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
    dbg_printf(do_print,"   |");
    for(ic2=0;ic2<par->n_spec_vary;ic2++)
      dbg_printf(do_print," %lE",gsl_matrix_get(cov_numerical_b,ic1,ic2));
    dbg_printf(do_print,"|\n");
  }
  compute_covariance_sampling(par,rng,data,noise_w,x_spec_ml,pst_new,pst_old,cov_sampling,do_print);
  dbg_printf(do_print,"Sampling covariance :\n");
  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
    dbg_printf(do_print,"   |");
    for(ic2=0;ic2<par->n_spec_vary;ic2++)
      dbg_printf(do_print," %lE",gsl_matrix_get(cov_sampling,ic1,ic2));
    dbg_printf(do_print,"|\n");
  }
  exit(1);

  if(par->n_spec_vary>0) {
    dbg_printf(do_print,"Burning\n");
    n_updated=0;
    ratio_accepted=0;
    gsl_matrix_set_zero(cov_spec);
    gsl_matrix_set_zero(cov_save);
    for(i_sample=0;i_sample<par->n_samples_burn;i_sample++) {
      accepted=draw_spectral_indices_marginal(par,rng,data,noise_w,x_spec_old,pst_old,
					      mat_step,x_spec_new,pst_new); //b_{n+1}(A_{n+1})
      if(accepted) {
	PixelState *tmp=pst_old;
	memcpy(x_spec_old,x_spec_new,par->n_param_max*sizeof(flouble));
	pst_old=pst_new;
	pst_new=tmp;
      }
#ifdef _DEBUG
      if(ipix_big==par->dbg_ipix)
	memcpy(&(par->dbg_extra[i_sample*par->n_spec_vary]),x_spec_old,par->n_spec_vary*sizeof(flouble));
#endif //_DEBUG
      for(ic1=0;ic1<par->n_spec_vary;ic1++) {
	mean_spec[ic1]+=x_spec_old[ic1];
	for(ic2=0;ic2<par->n_spec_vary;ic2++) {
	  gsl_matrix_set(cov_spec,ic1,ic2,gsl_matrix_get(cov_spec,ic1,ic2)+
			 x_spec_old[ic1]*x_spec_old[ic2]);
	}
      }

      ratio_accepted+=accepted;
      if(i_sample%par->n_update_covar==par->n_update_covar-1) { //Update covariance
	if(ratio_accepted<=5) { //Check if too few samples were accepted
	  dbg_printf(do_print,"No samples were accepted! Restarting with smaller step size\n");
	  stepping_factor*=0.5;
	  restart_mcmc(par,pst_old,x_spec_old,mat_step,stepping_factor);
	  memcpy(x_spec_old,x_spec_ml,par->n_param_max*sizeof(flouble));
	  compute_marginalized_chi2(par,data,noise_w,x_spec_old,pst_old,0);
	  gsl_matrix_set_zero(cov_save);
	  n_updated=0;
	  i_sample=-1;
	}
	else {
	  //Compute mean in this batch
	  ratio_accepted/=par->n_update_covar;
	  dbg_printf(do_print,"%d Acceptance ratio %.2lf\n",i_sample,ratio_accepted);

	  dbg_printf(do_print,"Current mean :");
	  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
	    mean_spec[ic1]/=par->n_update_covar;
	    dbg_printf(do_print," %lf",mean_spec[ic1]);
	  }
	  dbg_printf(do_print,"\n");
	  //Compute covariance in this batch
	  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
	    for(ic2=0;ic2<par->n_spec_vary;ic2++) {
	      flouble cov=gsl_matrix_get(cov_spec,ic1,ic2)/par->n_update_covar-
		mean_spec[ic1]*mean_spec[ic2];
	      gsl_matrix_set(cov_spec,ic1,ic2,cov);
	    }
	  }
	  decouple_covariance(par,cov_spec);
	  dbg_printf(do_print,"Current covariance :\n");
	  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
	    dbg_printf(do_print,"   |");
	    for(ic2=0;ic2<par->n_spec_vary;ic2++)
	      dbg_printf(do_print," %lE",gsl_matrix_get(cov_spec,ic1,ic2));
	    dbg_printf(do_print,"|\n");
	  }
	  //Add to mean covariance
	  if(i_sample>par->n_samples_burn/2) {
	    n_updated++;
	    for(ic1=0;ic1<par->n_spec_vary;ic1++) {
	      for(ic2=0;ic2<par->n_spec_vary;ic2++)
		gsl_matrix_set(cov_save,ic1,ic2,gsl_matrix_get(cov_save,ic1,ic2)+
			       gsl_matrix_get(cov_spec,ic1,ic2));
	    }
	  }

	  //Save diagonal
	  for(ic1=0;ic1<par->n_spec_vary;ic1++)
	    mean_spec[ic1]=sqrt(gsl_matrix_get(cov_spec,ic1,ic1));
	  //Cholesky decomposition of covariance
	  err=gsl_linalg_cholesky_decomp(cov_spec);
	  gsl_matrix_set_zero(mat_step);
	  if(err==GSL_EDOM) { //If covariance is not positive definite just save standard deviations
	    dbg_printf(do_print,"Covariance is not positive definite\n");
	    for(ic1=0;ic1<par->n_spec_vary;ic1++)
	      gsl_matrix_set(mat_step,ic1,ic1,factor_rescale*mean_spec[ic1]);
	  }
	  else { //Store cholesky decomposition in stepping function
	    for(ic1=0;ic1<par->n_spec_vary;ic1++) {
	      for(ic2=0;ic2<=ic1;ic2++)
		gsl_matrix_set(mat_step,ic1,ic2,factor_rescale*gsl_matrix_get(cov_spec,ic1,ic2));
	    }
	  }
	}

	gsl_matrix_set_zero(cov_spec);
	for(ic1=0;ic1<par->n_spec_vary;ic1++)
	  mean_spec[ic1]=0;
	ratio_accepted=0;
      }
    }
    dbg_printf(do_print,"Final covariance :\n");
    for(ic1=0;ic1<par->n_spec_vary;ic1++) { //Compute covariance
      dbg_printf(do_print,"   |");
      for(ic2=0;ic2<par->n_spec_vary;ic2++) {
	flouble cov=gsl_matrix_get(cov_save,ic1,ic2)/n_updated;
	gsl_matrix_set(cov_save,ic1,ic2,cov);
	dbg_printf(do_print," %lE",cov);
      }
      dbg_printf(do_print,"|\n");
      mean_spec[ic1]=sqrt(gsl_matrix_get(cov_save,ic1,ic1)); //Save diagonal
    }
    err=gsl_linalg_cholesky_decomp(cov_save); //Cholesky decomposition of covariance
    gsl_matrix_set_zero(mat_step);
    if(err==GSL_EDOM) { //If covariance is not positive definite just save standard deviations
      dbg_printf(do_print,"Final covariance is not positive definite!!!!\n");
      for(ic1=0;ic1<par->n_spec_vary;ic1++)
	gsl_matrix_set(mat_step,ic1,ic1,factor_rescale*mean_spec[ic1]);
    }
    else { //Store cholesky decomposition in stepping function
      for(ic1=0;ic1<par->n_spec_vary;ic1++) {
	for(ic2=0;ic2<=ic1;ic2++)
	  gsl_matrix_set(mat_step,ic1,ic2,factor_rescale*gsl_matrix_get(cov_save,ic1,ic2));
      }
    }
  }

  dbg_printf(do_print,"Starting actual sampling\n");
#ifdef _DEBUG
  if(ipix_big==par->dbg_ipix)
    ratio_accepted=0;
#endif //_DEBUG
  int n_samples_real=par->n_samples;
  if(par->n_spec_vary==0)
    n_samples_real=1;
  for(i_sample=0;i_sample<n_samples_real;i_sample++) {
    if(par->n_spec_vary>0) {
      accepted=draw_spectral_indices_marginal(par,rng,data,noise_w,x_spec_old,pst_old,
					      mat_step,x_spec_new,pst_new); //b_{n+1}(A_{n+1})
      if(accepted) {
	PixelState *tmp=pst_old;
	memcpy(x_spec_old,x_spec_new,par->n_param_max*sizeof(flouble));
	pst_old=pst_new;
	pst_new=tmp;
      }
    }

#ifdef _DEBUG
    if(par->n_spec_vary>0) {
      if(ipix_big==par->dbg_ipix) {
	memcpy(&(par->dbg_extra[(i_sample+par->n_samples_burn)*par->n_spec_vary]),
	       x_spec_old,par->n_spec_vary*sizeof(flouble));
	ratio_accepted+=accepted;
	if(i_sample && (i_sample%par->n_update_covar==par->n_update_covar-1)) {
	  ratio_accepted/=par->n_update_covar;
	  dbg_printf(do_print,"%d Acceptance ratio %.2lf\n",i_sample,ratio_accepted);
	  ratio_accepted=0;
	}
      }
    }
#endif //_DEBUG

    //Add amplitudes
    for(ipix=0;ipix<par->n_sub*par->n_pol;ipix++) {
      int index_here=ipix*par->n_comp;
      for(ic1=0;ic1<par->n_comp;ic1++) {
	amps_mean[ic1+index_here]+=gsl_vector_get(pst_old->vec_mean[ipix],ic1);
	for(ic2=0;ic2<par->n_comp;ic2++)
	  amps_covar[ic2+par->n_comp*(ic1+index_here)]+=gsl_matrix_get(pst_old->cov_inv[ipix],ic1,ic2);
      }
    }
    //Add indices
    for(ic1=0;ic1<par->n_spec_vary;ic1++) {
      par->map_indices_mean[ic1+ip_spc]+=x_spec_old[ic1];
      for(ic2=0;ic2<par->n_spec_vary;ic2++)
	par->map_indices_covar[ic1+par->n_spec_vary*(ic2+ip_spc)]+=x_spec_old[ic1]*x_spec_old[ic2];
    }
  }

  //Compute averages
  for(ipix=0;ipix<par->n_sub*par->n_pol;ipix++) {
    int index_here=ipix*par->n_comp;
    for(ic1=0;ic1<par->n_comp;ic1++) {
      amps_mean[ic1+index_here]/=n_samples_real;
      for(ic2=0;ic2<par->n_comp;ic2++)
	amps_covar[ic2+par->n_comp*(ic1+index_here)]/=n_samples_real;
    }
  }
  for(ic1=0;ic1<par->n_spec_vary;ic1++)
    par->map_indices_mean[ic1+ip_spc]/=n_samples_real;

  for(ic1=0;ic1<par->n_spec_vary;ic1++) {
    for(ic2=0;ic2<par->n_spec_vary;ic2++) {
      par->map_indices_covar[ic1+par->n_spec_vary*(ic2+ip_spc)]/=n_samples_real;
      par->map_indices_covar[ic1+par->n_spec_vary*(ic2+ip_spc)]-=
	par->map_indices_mean[ic1+ip_spc]*par->map_indices_mean[ic2+ip_spc];
    }
  }

  free(x_spec_old);
  free(x_spec_new);
  free(x_spec_ml);
  free(mean_spec);
  gsl_matrix_free(mat_step);
  gsl_matrix_free(cov_save);
  gsl_matrix_free(cov_spec);
}
*/

  /*
#ifdef _DEBUG_SINGLEPIX
#define NGRID_BETA_S 256
#define NGRID_BETA_D 128
#define BETA_S_MIN -1.5
#define BETA_S_MAX -0.5
#define BETA_D_MIN 1.5
#define BETA_D_MAX 1.75
  if(ipix_big==par->dbg_ipix) {
    int isam,ip_s,ip_d;
    char fffname[256];
    FILE *fff;
    double chiret[3];
    flouble *x_spec_here=my_malloc(par->n_param_max*sizeof(flouble));

    printf("Computing brute-force marginalized likelihood\n");
    sprintf(fffname,"test_chi2_pix%d.txt",ipix_big);
    fff=my_fopen(fffname,"w");
    memcpy(x_spec_here,x_spec_old,par->n_param_max*sizeof(flouble));
    for(ip_s=0;ip_s<NGRID_BETA_S;ip_s++) {
      x_spec_here[par->index_beta_s_t]=BETA_S_MIN+
	(BETA_S_MAX-BETA_S_MIN)*(ip_s+0.5)/NGRID_BETA_S;
      for(ip_d=0;ip_d<NGRID_BETA_D;ip_d++) {
	x_spec_here[par->index_beta_d_t]=BETA_D_MIN+
	  (BETA_D_MAX-BETA_D_MIN)*(ip_d+0.5)/NGRID_BETA_D;
	compute_marginalized_lpdf(par,data,noise_w,x_spec_here,pst);
	fprintf(fff,"%lE \n",pst->chi2);
      }
    }
    fclose(fff);
    free(x_spec_here);
  }
#endif //_DEBUG_SINGLEPIX
  */
