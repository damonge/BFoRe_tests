#ifndef _COMMON_BFORE
#define _COMMON_BFORE

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#ifdef _WITH_OMP
#include <omp.h>
#endif //_WITH_OMP
#ifdef _WITH_MPI
#include <mpi.h>
#endif //_WITH_MPI
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multimin.h>

#ifdef _LONGIDS
typedef long lint;
#else //_LONGIDS
typedef int lint;
#endif //_LONGIDS

#ifdef _SPREC
typedef float flouble;
typedef float complex fcomplex;
#ifdef _WITH_MPI
#define FLOUBLE_MPI MPI_FLOAT
#endif //_WITH_MPI
#else //_SPREC
typedef double flouble;
typedef double complex fcomplex;
#ifdef _WITH_MPI
#define FLOUBLE_MPI MPI_DOUBLE
#endif //_WITH_MPI
#endif //_SPREC

extern int NNodes;
extern int NodeThis;
extern int IThread0;

typedef struct {
  int nside; //Map nside
  int nside_spec; //Nside for spectral indices
  int n_side_sub; //Number of sub_pixels per side in each spectral index pixel
  int n_sub; //Number of sub_pixels in each spectral index pixel
  int n_pix; //Number of pixels per map
  int n_pix_spec; //Number of spectral index pixels
  int n_pix_spec_unmasked; //Number of unmasked spectral index pixels

  char input_data_prefix[256]; //Prefix pointing to the data maps
  flouble *maps_data; //Data frequency maps

  char input_noise_prefix[256]; //Prefix pointing to the noixe variance maps
  flouble *maps_noise_weight; //Noise weight maps (basically 1/noise_variance_per_pixel)

  char input_beta_s_t_prior[256]; //Path to prior map on beta_s (temperature)
  char input_beta_s_p_prior[256]; //Path to prior map on beta_s (polarization)
  char input_beta_d_t_prior[256]; //Path to prior map on beta_d (temperature)
  char input_beta_d_p_prior[256]; //Path to prior map on beta_d (polarization)
  char input_temp_d_t_prior[256]; //Path to prior map on temp_d (temperature)
  char input_temp_d_p_prior[256]; //Path to prior map on temp_d (polarization)
  flouble *map_prior_centres; //Maps of the spectral index prior mean
  flouble *map_prior_widths; //Maps of the spectral index prior width

  char input_mask_fname[256]; //Mask filename
  int *ipix_unmasked; //Indices of unmasked pixels

  char output_prefix[256]; //Output prefix
  int flag_write_samples; //Do we want to output samples?
  flouble *map_components_mean; //Mean of amplitudes
  flouble *map_components_covar; //Covariance of amplitudes
  flouble *map_indices_mean; //Mean of spectral indices
  flouble *map_indices_covar; //Covariance of spectral indices
  flouble *map_chi2; //Chi^2 map

  char fname_nulist[256]; //File containing frequencies
  int n_nu; //Number of frequencies
  flouble *freqs; //Frequencies

  int flag_include_polarization; //Do we have polarization?
  int n_pol; //Number of polarization channels (1-> T, 3-> T,Q,U)

  int flag_include_cmb; //Include CMB in sky model?
  int flag_include_synchrotron; //Incude synchrotron in sky model?
  int flag_include_dust; //Include dust in sky model?
  int flag_include_volume_prior; //Use volume (Jeffeys) prior?
  int flag_use_marginal; //Sample spectral indices from marginal distribution?
  int n_comp; //Number of components (up to 3)
  int index_cmb; //Index for CMB component
  int index_synchrotron; //Index for synchrotron component
  int index_dust; //Index for dust component

  int flag_independent_polarization; //Assume independent spectral indices in polarization?
  int flag_beta_s_free; //Is beta_s free?
  int flag_beta_d_free; //Is beta_d free?
  int flag_temp_d_free; //Is temp_d free?
  int n_param_max; //Maximum number of parameters to sample
  int n_spec_vary; //Number of free spectral indices
  int n_dof_pix; //Number of degrees of freedom per spectral index pixel
  int index_beta_s_t; //Index for beta_s (temperature)
  int index_beta_s_p; //Index for beta_s (polarization)
  int index_beta_d_t; //Index for beta_d (temperature)
  int index_beta_d_p; //Index for beta_d (polarization)
  int index_temp_d_t; //Index for temp_d (temperature)
  int index_temp_d_p; //Index for temp_d (polarization)
  flouble beta_s_step; //Initial step size in beta_s
  flouble beta_d_step; //Initial step size in beta_d
  flouble temp_d_step; //Initial step size in temp_d
  flouble nu0_s; //Reference frequency for synchrotron
  flouble nu0_d; //Reference frequency for dust

  unsigned long seed; //Seed
  int n_samples; //Total number of samples per spectral index pixel
  int n_output_rate; //Output sample every so many samples
  flouble frac_samples_burn; //Fraction of samples used for burning
  int n_samples_burn; //Number of samples used for burning
  int n_update_covar; //Number of samples used to compute the initial covariance matrix
  int n_spec_resample; //Number of spectral index samples takeng for each amplitudes sample

  int ipix_0; //First pixel corresponding to this node
  int ipix_f; //Last pixel for this node (this one actually corresponds to the next node)

  int dbg_ipix; //Pixel index used for debugging
  flouble *dbg_extra; //Debugging data
} ParamBFoRe;

//Defined in common.c
int my_linecount(FILE *f);
void report_error(int level,char *fmt,...);
void *my_malloc(size_t size);
void *my_calloc(size_t nmemb,size_t size);
FILE *my_fopen(const char *path,const char *mode);
size_t my_fwrite(const void *ptr, size_t size, size_t nmemb,FILE *stream);
size_t my_fread(void *ptr,size_t size,size_t count,FILE *stream);
void param_bfore_free(ParamBFoRe *par);
ParamBFoRe *read_params(char *fname);
void write_output(ParamBFoRe *par);
void dbg_printf(int do_print,char *fmt,...);

//Defined in healpix_extra.c
#ifdef _WITH_SHT
long he_nalms(int lmax);
long he_indexlm(int l,int m,int lmax);
void he_alm2map(int nside,int lmax,int ntrans,flouble **maps,fcomplex **alms);
void he_map2alm(int nside,int lmax,int ntrans,flouble **maps,fcomplex **alms);
void he_anafast(flouble **maps_1,flouble **maps_2,
		int nmaps_1,int nmaps_2,
		int pol_1,int pol_2,
		flouble **cls,int nside,int lmax);
double *he_generate_beam_window(int lmax,double fwhm_amin);
void he_alter_alm(int lmax,double fwhm_amin,fcomplex *alms,double *window);
#endif //_WITH_SHT
void he_write_healpix_map(flouble **tmap,int nfields,long nside,char *fname);
flouble *he_read_healpix_map(char *fname,long *nside,int nfield);
int he_ring_num(long nside,double z);
void he_ring2nest_inplace(flouble *map_in,long nside);
void he_nest2ring_inplace(flouble *map_in,long nside);
void he_udgrade(flouble *map_in,long nside_in,
		flouble *map_out,long nside_out,
		int nest);

//Defined in rng.c
#define RNG_NRAN 624
#define RNG_MRAN 397
#define RNG_MATRIX_A 0x9908b0df
#define RNG_UPPER_MASK 0x80000000UL
#define RNG_LOWER_MASK 0x7fffffffUL
typedef struct {
  unsigned long mt[RNG_NRAN];
  int mti;
  int calc_gauss;
  double u;
  double phi;
} Rng;
Rng *init_rng(unsigned long seed);
void end_rng(Rng *rng);
unsigned long rand_ulong(Rng *rng);
double rand_real01(Rng *rng);
double rand_gauss(Rng *rng);

//Defined in powell.c
#define TINY 1.0E-25
#define ZEPS 1.0E-10
#define GLIMIT 100.0
#define GOLD 1.618034
#define CGOLD 0.3819660
#define LIN_TOL 2.0E-4
#define LIN_ITMAX 100
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);  

typedef struct {
  int n;
  double *p;
  double **xi;
  int iter;
  int max_iter;
  double fret;
  double ftol;
  double *xdum;
  double *xdir;
  double (*fun)(double *,void *);
  void *params;
} PowellParams;
void free_powell_params(PowellParams *par);
PowellParams *powell_params_new(int n,flouble *p,flouble (*fun)(flouble *,void *),
				void *params,int max_iter,flouble ftol);
void powell(PowellParams *par);

//Defined in bfore.c
typedef struct {
  flouble *f_matrix; //Frequency evolution matrix
  gsl_matrix **cov_inv; //Set of covariance matrices for component amplitudes (one per pixel)
  gsl_vector **vec_mean; //Set of vectors of mean component amplitudes (one per pixel)
  flouble *prior_mean; //Prior mean (one per spectral param)
  flouble *prior_isigma; //1/sigma of prior (one per spectral param)
  flouble *rand_spec; //Dummy array to fill with random numbers (one element per spectral param)
  gsl_vector *vaux; //Dummy vector (one element per component)
  double chi2; //Pixel chi2
} PixelState;
PixelState *pixel_state_new(ParamBFoRe *par);
void pixel_state_free(PixelState *pst,ParamBFoRe *par);
void clean_pixel(ParamBFoRe *par,Rng *rng,PixelState *pst,int ipix_big);
void clean_pixel_from_marginal(ParamBFoRe *par,Rng *rng,PixelState *pst_old,
			       PixelState *pst_new,int ipix_big);

#endif //_COMMON_BFORE
