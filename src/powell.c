#include "common.h"

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

static flouble get_powell_1d(PowellParams *par,flouble x)
{
  int j;
  for(j=0;j<par->n;j++)
    par->xdum[j]=par->p[j]+x*par->xdir[j];
  
  return (*(par->fun))(par->xdum,par->params);
}

void free_powell_params(PowellParams *par)
{
  int i;
  for(i=0;i<par->n;i++) 
    free(par->xi[i]);
  free(par->p);
  free(par->xi);
  free(par->xdum);
  free(par->xdir);
  free(par);
}

PowellParams *powell_params_new(int n,flouble *p,flouble (*fun)(flouble *,void *),
				void *params,int max_iter,flouble ftol)
{
  int ii;
  PowellParams *par=(PowellParams *)my_malloc(sizeof(PowellParams));
  par->n=n;
  par->p=(flouble *)my_malloc(n*sizeof(flouble));
  for(ii=0;ii<n;ii++)
    par->p[ii]=p[ii];
  par->xi=(flouble **)my_malloc(n*sizeof(flouble *));
  for(ii=0;ii<n;ii++) {
    par->xi[ii]=(flouble *)my_calloc(n,sizeof(flouble));
    par->xi[ii][ii]=1;
  }
  par->xdum=(flouble *)my_malloc(n*sizeof(flouble));
  par->xdir=(flouble *)my_malloc(n*sizeof(flouble));
  par->iter=0;
  par->max_iter=max_iter;
  par->fun=fun;
  par->fret=(*fun)(p,params);
  par->params=params;
  par->ftol=ftol;

  return par;
}

static void mnbrak(PowellParams *par,flouble *ax,flouble *bx,flouble *cx,
		   flouble *fa,flouble *fb,flouble *fc)
{
  flouble ulim,u,r,q,fu,dum;
  *fa=get_powell_1d(par,*ax);
  *fb=get_powell_1d(par,*bx);
  if(*fb > *fa) {
    SHFT(dum,*ax,*bx,dum);
    SHFT(dum,*fb,*fa,dum);
  }
  *cx=(*bx)+GOLD*(*bx-*ax);
  *fc=get_powell_1d(par,*cx);
  while(*fb>*fc) {
    r=(*bx-*ax)*(*fb-*fc);
    q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
      (2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
    ulim=(*bx)+GLIMIT*(*cx-*bx);

    if((*bx-u)*(u-*cx) > 0) {
      fu=get_powell_1d(par,u);
      if(fu<*fc) {
	*ax=(*bx);
	*bx=u;
	*fa=(*fb);
	*fa=fu;
	return;
      }
      else if(fu>*fb) {
	*cx=u;
	*fc=fu;
	return;
      }
      u=(*cx)+GOLD*(*cx-*bx);
      fu=get_powell_1d(par,u);
    }
    else if((*cx-u)*(u-ulim) > 0.0) {
      fu=get_powell_1d(par,u);
      if(fu<*fc) {
	SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx));
	SHFT(*fb,*fc,fu,get_powell_1d(par,u));
      }
    }
    else if((u-ulim)*(ulim-*cx) > 0.0) {
      u=ulim;
      fu=get_powell_1d(par,u);
    }
    else {
      u=(*cx)+GOLD*(*cx-*bx);
      fu=get_powell_1d(par,u);
    }
    SHFT(*ax,*bx,*cx,u);
    SHFT(*fa,*fb,*fc,fu);
  }
}

static flouble brent(PowellParams *par,flouble ax,flouble bx,
		    flouble cx,flouble tol,flouble *xmin)
{
  int iter;
  flouble a,b,d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  flouble e=0.0;

  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=get_powell_1d(par,x);
  for(iter=0;iter<LIN_ITMAX;iter++) {
    xm=0.5*(a+b);
    tol2=2.*(tol1=tol*fabs(x)+ZEPS);
    if(fabs(x-xm)<=tol2-0.5*(b-a)) {
      *xmin=x;
      return fx;
    }
    if(fabs(e)>tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.*(q-r);
      if(q>0.) p=-p;
      q=fabs(q);
      etemp=e;
      e=d;
      if(fabs(p)>fabs(0.5*q*etemp) || p<=q*(a-x) || p>=q*(b-x))
	d=CGOLD*(e=(x>=xm ? a-x : b-x));
      else {
	d=p/q;
	u=x+d;
	if(u-a<tol2 || b-u<tol2)
	  d=SIGN(tol1,xm-x);
      }
    }
    else {
      d=CGOLD*(e=(x>=xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=get_powell_1d(par,u);
    
    if(fu<=fx) {
      if(u>=x) a=x; else b=x;
      SHFT(v,w,x,u);
      SHFT(fv,fw,fx,fu);
    }
    else {
      if(u<x) a=u; else b=u;
      if(fu<=fw || w==x) {
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      }
      else if(fu<=fv || v==x || v==w) {
	v=u;
	fv=fu;
      }
    }
  }
  report_error(0,"Brent didn't converge after %d iterations\n",LIN_ITMAX);
  *xmin=x;
  return fx;
}

static void linmin(PowellParams *par,flouble *xdir)
{
  int j;
  flouble xx,xmin,fx,fb,fa,bx,ax;

  for(j=0;j<par->n;j++) par->xdir[j]=xdir[j];

  ax=0.;
  xx=1.;
  mnbrak(par,&ax,&xx,&bx,&fa,&fx,&fb);
  par->fret=brent(par,ax,xx,bx,LIN_TOL,&xmin);

  for(j=0;j<par->n;j++) {
    xdir[j]*=xmin;
    par->p[j]+=xdir[j];
  }
}

void powell(PowellParams *par)
{
  int i,ibig,j;
  flouble del,fp,fptt,t,*pt,*ptt,*xit;

  pt=(flouble *)my_malloc(par->n*sizeof(flouble));
  ptt=(flouble *)my_malloc(par->n*sizeof(flouble));
  xit=(flouble *)my_malloc(par->n*sizeof(flouble));

  for(j=0;j<par->n;j++) pt[j]=par->p[j]; //Save initial point;

  for(par->iter=1;par->iter<par->max_iter;par->iter++) {
    fp=par->fret;
    ibig=0;
    del=0; //Will become the biggest function decrease
    for(i=0;i<par->n;i++) { 
      for(j=0;j<par->n;j++) xit[j]=par->xi[j][i];
      fptt=par->fret;
      linmin(par,xit);
      if(fptt-par->fret > del) { //Record largest decrease
	del=fptt-par->fret;
	ibig=i;
      }
    }
    if(2.*(fp-par->fret) < par->ftol*(fabs(fp)+fabs(par->fret))+TINY) {
      free(pt);
      free(ptt);
      free(xit);

      return;
    }
    if(par->iter==par->max_iter) {
      report_error(0,"Powell didn't converge after %d iterations\n",
		   par->max_iter);
      return;
    }
    for(j=0;j<par->n;j++) {
      ptt[j]=2.*par->p[j]-pt[j];
      xit[j]=par->p[j]-pt[j];
      pt[j]=par->p[j];
    }
    fptt=(*(par->fun))(ptt,par->params);
    if(fptt<fp) {
      flouble x=fp-par->fret-del;
      t=2.*(fp-2*par->fret+fptt)*x*x-del*(fp-fptt)*(fp-fptt);
      if(t<0.0) {
	linmin(par,xit);
	for(j=0;j<par->n;j++) {
	  par->xi[j][ibig]=par->xi[j][par->n-1];
	  par->xi[j][par->n-1]=xit[j];
	}
      }
    }
  }
}
