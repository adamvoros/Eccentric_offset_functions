/*****************************************************************************/
/* optsn.c								     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <malloc.h>

#include "downhill.h"
#include "lmfit.h"

/*****************************************************************************/

double get_gaussian(double mean,double stddev)
{
 double u,v,r;
 u=sqrt(-2.0*log(drand48()));
 v=2*M_PI*drand48();
 r=u*cos(v);
 return(r*stddev+mean);
}

/*****************************************************************************/

double solve_kepler_equ(double m,double e)
{
 double	s,d,d0;
 int	n;

 s=sin(m);

 if ( s==0.0 )
	return(m);

 else if ( e>=0.8 && s*s<0.1 && cos(m)>0.0 )
  {	if ( s>0 )	d=+pow(+6*s,1.0/3.0);
	else		d=-pow(-6*s,1.0/3.0);
  }
 else
  {	d=e*s/(1-sin(m+e)+s);
	if ( d<-e )	d=-e;
	else if ( d>e )	d=+e;
  }

 for ( n=15,d0=0.0 ; d != d0 && n>0 ; n-- )
  {	d0=d;
	d=d-(d-e*sin(m+d))/(1-e*cos(m+d));
  };

 return(m+d);
}

double eoq(double lambda,double k,double h)
{
 double	e,E;
 e=sqrt(h*h+k*k);
 if ( e<=0.0 )
	return(0.0); 
 else
  {	E=solve_kepler_equ(lambda-atan2(h,k),e);
	return(e*cos(E));
  }
}

double eop(double lambda,double k,double h)
{
 double	e,E;
 e=sqrt(h*h+k*k);
 if ( e<=0.0 )
	return(0.0); 
 else
  {	E=solve_kepler_equ(lambda-atan2(h,k),e);
	return(e*sin(E));
  }
 
}

int eoqp(double lambda,double k,double h,double *rq,double *rp)
{
 double	e,E;
 e=sqrt(h*h+k*k);
 if ( e<=0.0 )
  {	*rq=0.0;
	*rp=0.0;
	return(0);
  }
 else
  {	E=solve_kepler_equ(lambda-atan2(h,k),e);
	*rq=e*cos(E);
	*rp=e*sin(E);
	return(0);
  }
}

double arg(double x,double y)
{
 return(atan2(y,x));
}

double lamtran(double l0,double k,double h)
{
 double l,c,s,j;

 c=cos(l0);
 s=sin(l0);
 j=sqrt(1-k*k-h*h);

 l=arg(k+c+h*(k*s-h*c)/(1+j),h+s-k*(k*s-h*c)/(1+j))-(k*s-h*c)*j/(1+k*c+h*s);

 return(l);
}
double dlamtrandk(double l0,double k,double h)
{
 double c,s,j,dldk,ep;

 c=cos(l0);
 s=sin(l0);
 j=sqrt(1-k*k-h*h);
 ep=k*c+h*s;
 
 dldk=-h/(1+j)-j*(h+(2+ep)*s)/((1+ep)*(1+ep));

 return(dldk);
}
double dlamtrandh(double l0,double k,double h)
{
 double c,s,j,dldh,ep;

 c=cos(l0);
 s=sin(l0);
 j=sqrt(1-k*k-h*h);
 ep=k*c+h*s;
 
 dldh=+k/(1+j)+j*(k+(2+ep)*c)/((1+ep)*(1+ep));

 return(dldh);
}

double vy0(double l,double k,double h)
{
 double	p,q,c,vy,j;

 eoqp(l,k,h,&q,&p);

 c=cos(l+p);
 j=sqrt(1-k*k-h*h);
 
 vy=(c-k*q/(1+j))/(1-q);

 return(vy);
}

double dvy0dk(double l,double k,double h)
{
 double	p,q,c,s,dv,j;

 eoqp(l,k,h,&q,&p);

 c=cos(l+p);
 s=sin(l+p);
 j=sqrt(1-k*k-h*h);

 dv=(c*c-s*s+q*s*s-c*k-k*(c-k)/(1+j))/((1-q)*(1-q)*(1-q))-q*(1+k*k/(j*(1+j)))/((1-q)*(1+j));

 return(dv);
}

double dvy0dh(double l,double k,double h)
{
 double	p,q,c,s,dv,j;

 eoqp(l,k,h,&q,&p);

 c=cos(l+p);
 s=sin(l+p);
 j=sqrt(1-k*k-h*h);

 dv=(2*s*c-q*s*c-c*h-k*(s-h)/(1+j))/((1-q)*(1-q)*(1-q))-q*k*h/((1-q)*(1+j)*(1+j)*j);

 return(dv);
}
double dvy0dl(double l,double k,double h)
{
 double	p,q,c,s,dv,j;

 eoqp(l,k,h,&q,&p);

 c=cos(l+p);
 s=sin(l+p);
 j=sqrt(1-k*k-h*h);
 
 dv=(-s+h+p/(1+j)*k)/((1-q)*(1-q)*(1-q));

 return(dv);
}

double get_sn_ratio(double a,double b,double k,double h,double *ls,int nl)
{
 int	n,nvar;
 double	**amatrix,*bvector,*fvars,w,sn,l,y,l0;
 int	i,j;
 
 /* vy = b + a*vy0(l,k,h) */
 /* a, b, k, h */
 nvar=4;
 amatrix=matrix_alloc(nvar);
 bvector=vector_alloc(nvar);
 fvars  =vector_alloc(nvar);

 for ( i=0 ; i<nvar ; i++ )
  {	for ( j=0 ; j<nvar ; j++ )
	 {	amatrix[i][j]=0.0;		}
	bvector[i]=0.0;
  }

 for ( n=0 ; n<nl ; n++ )
  {	l0=lamtran(M_PI/2.0,k,h);
	l=ls[n]+l0;

	fvars[0]=vy0(l,k,h);
	fvars[1]=1;
	/*
	fvars[2]=a*(dvy0dk(l,k,h));
	fvars[3]=a*(dvy0dh(l,k,h));
	*/
	fvars[2]=a*(dvy0dk(l,k,h)+dvy0dl(l,k,h)*dlamtrandk(M_PI/2.0,k,h));
	fvars[3]=a*(dvy0dh(l,k,h)+dvy0dl(l,k,h)*dlamtrandh(M_PI/2.0,k,h));
	

	y=b+a*fvars[0];

	w=1.0;
	for ( i=0 ; i<nvar ; i++ )
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	amatrix[i][j] += w*fvars[i]*fvars[j];	}
		bvector[i] += w*y*fvars[i];
	 }
  }

 if ( invert_gauss(amatrix,nvar) )
	sn=0.0;
 else
  {	sn=amatrix[2][2]*amatrix[3][3]-amatrix[2][3]*amatrix[3][2]; 
	/* sn=amatrix[0][0]; */
	if ( sn <= 0.0 )
		sn=0.0;
	else
		sn=1.0/sn;
  }

 vector_free(fvars);
 vector_free(bvector);
 matrix_free(amatrix);

 return(sn);
}

int compare_double(const void *v1,const void *v2)
{
 if ( *(double *)v1 < *(double *)v2 )
	return(-1);
 else if ( *(double *)v1 > *(double *)v2 )
	return(+1);
 else
	return(0);
}

typedef struct
 {	double	a,b,k,h;
	int	nlong;
 } khsnparam;

double dh_callback(void *param,double *arr)
{
 double		sn;
 khsnparam	*p=(khsnparam *)param;

 sn=-get_sn_ratio(p->a,p->b,p->k,p->h,arr,p->nlong);

 return(sn);
}

int main(int argc,char *argv[])
{
 double	a,b,k,h;
 double	*phase,*longs,*ophss,dp;
 double	*fixphases;
 int	i,nfixphase,nlong,seed,niter,no_sort;
 int	method;

 k=0.0,h=0.0;

 seed=0;
 nlong=4;
 niter=10000;
 no_sort=0;

 method=2;
 fixphases=NULL;
 nfixphase=0;

 for ( i=1 ; i<argc ; i++ )
  {	if ( strcmp(argv[i],"-h")==0 || strcmp(argv[i],"--help")==0 )
	 {	fprintf(stderr,"Usage:\toptsn [-h|--help]\n");
		fprintf(stderr,"\t[-1|-2] [-e|--eccentricity <k>,<h>] [-s|--seed <seed>]\n");
		fprintf(stderr,"\t[-p|--phases <list-of-fixed-phases>]\n");
		fprintf(stderr,"\t[-n|--nlong <nlong>] [-i|--iterations <iterations>] [-q|--no-sort]\n");
		return(0);
	 }
	else if ( (strcmp(argv[i],"-e")==0 || strcmp(argv[i],"--eccentricity")==0) && i<argc-1 )
	 {	i++;
		sscanf(argv[i],"%lg,%lg",&k,&h);
	 }
	else if ( (strcmp(argv[i],"-s")==0 || strcmp(argv[i],"--seed")==0) && i<argc-1 )
	 {	i++;
		sscanf(argv[i],"%d",&seed);
	 }
	else if ( (strcmp(argv[i],"-n")==0 || strcmp(argv[i],"--nlong")==0) && i<argc-1 )
	 {	i++;
		sscanf(argv[i],"%d",&nlong);
	 }
	else if ( (strcmp(argv[i],"-i")==0 || strcmp(argv[i],"--iterations")==0) && i<argc-1 )
	 {	i++;
		sscanf(argv[i],"%d",&niter);
	 }
	else if ( (strcmp(argv[i],"-p")==0 || strcmp(argv[i],"--phase")==0 || strcmp(argv[i],"--phases")==0) && i<argc-1 )
	 {	char	*p;
		i++;
		p=argv[i];
		nfixphase=0;
		fixphases=NULL;
		while ( *p )
		 {	double	d;
			int	k;
			if ( sscanf(p,"%lg",&d)<1 )
				break;
			fixphases=(double *)realloc(fixphases,sizeof(double)*(nfixphase+1));
			fixphases[nfixphase]=d-floor(d);
			nfixphase++;
			while ( *p )
			 {	k=0;
				while ( *p==',' || *p==' ' )
				 {	p++;k++;			}
				if ( k )	
					break;
				p++;
			 }
		 }
	 }
	else if ( (strcmp(argv[i],"-q")==0 || strcmp(argv[i],"--no-sort")==0) )
		no_sort=1;
	else if ( strcmp(argv[i],"--sort")==0 )
		no_sort=0;
 	/* else if ( strcmp(argv[i],"-0")==0 )	method=0; */
	else if ( strcmp(argv[i],"-1")==0 )	method=1;
	else if ( strcmp(argv[i],"-2")==0 )	method=2;
	else
	 {	fprintf(stderr,"Error: invalid command line argument near '%s'.\n",argv[i]);
		return(1);
	 }
  }

 a=1.0;
 b=0.0;

 phase=(double *)malloc(sizeof(double)*nlong);
 longs=(double *)malloc(sizeof(double)*nlong);
 ophss=(double *)malloc(sizeof(double)*nlong);


/*
 if ( method == 0 )
  {	dp=0.02;
	for ( i=0 ; i<nlong ; i++ )
	 {	phase[i]=dp*(double)i;		}

	while ( 1 )
	 {	for ( i=0 ; i<nlong ; i++ )
		 {	longs[i]=phase[i]*2*M_PI;		}
		sn=get_sn_ratio(a,b,k,h,longs,nlong);

		fprintf(stdout,"%4.2f %4.2f %4.2f %4.2f %12g\n",phase[0],phase[1],phase[2],phase[3],sn);

		for ( i=nlong-1,dm=1-0.5*dp ; i>=0 ; i--,dm-=dp )
		 {	phase[i]+=dp;
			if ( phase[i]<dm )
				break;
		 }
		if ( i<0 )
			break;
		dm=phase[i]+dp;
		for ( i++ ; i<nlong ; i++ )
		 {	phase[i]=dm;
			dm+=dp;
		 }
	 }
  }
*/

 if ( method==1 )
  {	double		**sphases,sn0;
	int		i,j;
	FILE		*fw;
	khsnparam	p;

	fw=stdout;

	sphases=(double **)matrix_alloc_gen(nlong+1,nlong+1);

	srand48(seed);

	/*
	for ( i=0 ; i<nlong ; i++ )
	 {	phase[i]=2.0*M_PI*((double)i+0.5)/(double)nlong;		}
	*/

	for ( i=0 ; i<nlong ; i++ )
	 {	phase[i]=2.0*M_PI*drand48();				}
	for ( i=0 ; i<nfixphase ; i++ )
	 {	phase[i]=2.0*M_PI*fixphases[i];				}
			
	for ( i=0 ; i<nlong ; i++ )
	 {	sphases[0][i]=phase[i];		}
	for ( i=0 ; i<nlong-nfixphase ; i++ )
	 {	for ( j=0 ; j<nlong ; j++ )
		 {	sphases[i+1][j]=phase[j]+((i+nfixphase)==j?0.01:0.00);	}
	 }

	p.a=a,p.b=b;
	p.k=k,p.h=h;
	p.nlong=nlong;
	downhill_simplex(sphases,nlong,nlong-nfixphase,dh_callback,&p,1e-10,2000);

	sn0=dh_callback(&p,sphases[0]);

	for ( i=0 ; i<nlong ; i++ )
	 {	phase[i]=sphases[0][i]/(2*M_PI);
		phase[i]-=floor(phase[i]);
		ophss[i]=phase[i];
	 }

	qsort(ophss,nlong,sizeof(double),compare_double);
	if ( no_sort )
	 {	for ( i=0 ; i<nlong ; i++ )
	 	 {	fprintf(fw,"%7.5f ",phase[i]);	}
	 }
	else
	 {	for ( i=0 ; i<nlong ; i++ )
	 	 {	fprintf(fw,"%7.5f ",ophss[i]);	}
	 }

	fprintf(fw,"# %12g\n",-sn0);

  }
 else if ( method==2 )
  {	int	acnt,tcnt;
	double	sn,sn0,*prtph;
	double	vx,alpha,u;
	FILE	*fw;

	srand48(seed);

	dp=0.01;
	/*
	for ( i=0 ; i<nlong ; i++ )
	 {	phase[i]=((double)i+0.5)/(double)nlong;		}
	*/
	for ( i=0 ; i<nlong ; i++ )
	 {	phase[i]=drand48();				}
	/*
	phase[0]=0.20, 
	phase[1]=0.42, 
	phase[2]=0.58, 
	phase[3]=0.80;
	*/

	acnt=tcnt=0;

	fw=stdout;

	for ( i=0 ; i<nlong ; i++ )
	 {	longs[i]=phase[i]*2*M_PI;		}
	sn0=get_sn_ratio(a,b,k,h,longs,nlong);

	prtph=(double *)malloc(sizeof(double)*nlong);


	while ( acnt<niter )
	 {	for ( i=0 ; i<nlong ; i++ )
		 {	prtph[i]=get_gaussian(phase[i],dp);
			prtph[i]-=floor(prtph[i]);
			longs[i]=prtph[i]*2*M_PI;
		 }
		sn=get_sn_ratio(a,b,k,h,longs,nlong);
		vx=exp(0.5*(sn-sn0)*10.0);
		alpha=(vx<1.0?vx:1.0);
		if ( alpha>0.0 && (u=drand48()) <= alpha )
		 {	for ( i=0 ; i<nlong ; i++ )
			 {	ophss[i]=phase[i]=prtph[i];		}
			qsort(ophss,nlong,sizeof(double),compare_double);
			if ( no_sort )
			 {	for ( i=0 ; i<nlong ; i++ )
			 	 {	fprintf(fw,"%6.4f ",phase[i]);	}
			 }
			else
			 {	for ( i=0 ; i<nlong ; i++ )
			 	 {	fprintf(fw,"%6.4f ",ophss[i]);	}
			 }
			fprintf(fw,"\n");
			acnt++;
			sn0=sn;
		 }
		tcnt++;
	 }
	fprintf(stderr,"# %d/%d (%.6f)\n",acnt,tcnt,(double)acnt/(double)tcnt);
  }
 else
  {	fprintf(stderr,"%s: error: invalid method code: %d.\n",argv[0],method);
	return(1);
  }
	
 return(0);
}
