#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "complex.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


void fft(complex float* x, complex float* X, int N){
  complex float W;
  complex float** A1;
  int j1,j0,k1,k0,r1,r2,i,m;

  complex float** cp1;

  r1=(int) sqrt(N);

  while(N%r1!=0)
    r1--;

  /*fprintf(stderr,"%d %d %d\n",N,N%r1,r1);*/
  
  r2=N/r1;

  A1=(complex float**) malloc(r1*sizeof(complex float*));

  for(i=0; i<N; ++i){
    X[i]=0;
  }

  for(i=0; i<r1; ++i){
    A1[i]=(complex float*) malloc(r2*sizeof(complex float));
    for(m=0; m<r2; ++m){
      A1[i][m]=0;
    }
  }

  W=cexp(-2*M_PI*I/N);

  cp1=(complex float**) malloc(r1*sizeof(complex float*));
  for(j0=0; j0<r1; ++j0){
    cp1[j0]=(complex float*) malloc(r1*sizeof(complex float));
    for(k1=0; k1<r1; ++k1){
      cp1[j0][k1]=cpow(W,j0*k1*r2);
    }
  }
      
  /*fprintf(stderr,"%d %d\n",r1,r2);*/

  for(j0=0; j0<r1; ++j0){
    for(k0=0; k0<r2; ++k0){
      for(k1=0; k1<r1; ++k1){
	/*A1[j0][k0]+=x[k1*r2+k0]*cpow(W,j0*k1*r2);*/
	A1[j0][k0]+=x[k1*r2+k0]*cp1[j0][k1];
      }
    }  
  }
  for(j1=0; j1<r2; ++j1){
    for(j0=0; j0<r1; ++j0){
      for(k0=0; k0<r2; ++k0){
	X[j1*r1+j0]+=A1[j0][k0]*cpow(W,(j1*r1+j0)*k0);
      }     
    }
  }
  
  
  return;
}
  

void ifft(complex float* X, complex float* x, int N){
  complex float* xr;
  int i;
  xr=(complex float*) malloc(sizeof(complex float)*N);
  fft(X, xr, N);
  x[0]=xr[0];
  for(i=1; i<N; ++i){
    x[i]=xr[N-i];
  }
  for(i=0; i<N; ++i){
    x[i]=x[i]/N;
  }
  free(xr);
}

