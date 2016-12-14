#include "common.h"

void mpi_init(int *p_argc,char ***p_argv)
{
#ifdef _WITH_MPI
  int ii,nthreads_this;
  int *nthreads_all;
  MPI_Init(p_argc,p_argv);

  MPI_Comm_size(MPI_COMM_WORLD,&NNodes);
  MPI_Comm_rank(MPI_COMM_WORLD,&NodeThis);

  nthreads_all=my_malloc(NNodes*sizeof(int));
#ifdef _WITH_OMP
  nthreads_this=omp_get_max_threads();
#else //_WITH_OMP
  nthreads_this=1;
#endif //_WITH_OMP
  MPI_Allgather(&nthreads_this,1,MPI_INT,nthreads_all,1,MPI_INT,MPI_COMM_WORLD);
  if(NodeThis==0) {
    for(ii=0;ii<NNodes;ii++)
      printf("Node %d has %d threads\n",ii,nthreads_all[ii]);
  }
  IThread0=0;
  for(ii=0;ii<NodeThis;ii++)
    IThread0+=nthreads_all[ii];
#ifdef _DEBUG
  printf("Node %d, thread count starts at %d\n",NodeThis,IThread0);
#endif //_DEBUG
  free(nthreads_all);

#else //_WITH_MPI

  NNodes=1;
  NodeThis=0;
  IThread0=0;
#endif //_WITH_MPI
}

int main(int argc,char **argv)
{ 
  char fname_init[256];
  if(argc!=2) {
    printf("Usage: fg_rm.x param_file\n");
    exit(0);
  }
  sprintf(fname_init,"%s",argv[1]);

  mpi_init(&argc,&argv);
  gsl_set_error_handler_off();
  ParamBFoRe *par=read_params(fname_init);

  int ii,n_threads;
#ifdef _WITH_OMP
  n_threads=omp_get_max_threads();
#else //_WITH_OMP
  n_threads=1;
#endif //_WITH_OMP
  PixelState **pst_old,**pst_new=NULL;

  pst_old=my_malloc(n_threads*sizeof(PixelState *));
  for(ii=0;ii<n_threads;ii++)
    pst_old[ii]=pixel_state_new(par);
  if(par->flag_use_marginal) {
    pst_new=my_malloc(n_threads*sizeof(PixelState *));
    for(ii=0;ii<n_threads;ii++)
      pst_new[ii]=pixel_state_new(par);
  }

#ifdef _WITH_OMP
#ifndef _DEBUG_SINGLEPIX
#pragma omp parallel default(none) shared(par,IThread0,NodeThis,pst_old,pst_new)
#endif //_DEBUG_SINGLEPIX
#endif //_WITH_OMP
  {
    int ipix_big,ithr;
    unsigned long seed_thr;
    Rng *rng;

    ithr=0;
#ifdef _WITH_OMP
#ifndef _DEBUG_SINGLEPIX
    ithr=omp_get_thread_num();
#endif //_DEBUG_SINGLEPIX
#endif //_WITH_OMP
    seed_thr=par->seed+IThread0+ithr;
    rng=init_rng(seed_thr);

#ifdef _DEBUG_SINGLEPIX
    ipix_big=par->dbg_ipix;
#else //_DEBUG_SINGLEPIX
#ifdef _WITH_OMP
#pragma omp for
#endif //_WITH_OMP
    for(ipix_big=par->ipix_0;ipix_big<par->ipix_f;ipix_big++)
#endif //_DEBUG_SINGLEPIX
      {
	int ip=par->ipix_unmasked[ipix_big];
	printf("Node %d, thread %d, pixel %d\n",NodeThis,ithr,ip);
	if(par->flag_use_marginal)
	  clean_pixel_from_marginal(par,rng,pst_old[ithr],pst_new[ithr],ip); //TODO: This needs checking
	else
	  clean_pixel(par,rng,pst_old[ithr],ip); //TODO: This needs checking
      }//end omp for
    end_rng(rng);
  }//end omp parallel

  for(ii=0;ii<n_threads;ii++)
    pixel_state_free(pst_old[ii],par);
  free(pst_old);
  if(par->flag_use_marginal) {
    for(ii=0;ii<n_threads;ii++)
      pixel_state_free(pst_new[ii],par);
    free(pst_new);
  }

  if(NodeThis==0)
    printf("Writing output\n");
  write_output(par);

  param_bfore_free(par);

#ifdef _WITH_MPI
  MPI_Finalize();
#endif //_WITH_MPI

  return 0;
}
