#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "complex.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void fft(complex float* x, complex float* X, int N){
  complex float* xp;
  complex float* Xp;
  int depth;
  int i;
  int Np;
  
  complex float* coeffs;
  complex float s;
  complex float s1;
  complex float z;
  int j,k;
  int p2d,N_over_p2d;
  
  for(i=0; pow(2,i)<N; ++i);
  Np=pow(2,i);
  xp=(complex float*) malloc(Np*sizeof(complex float));
  Xp=(complex float*) malloc(Np*sizeof(complex float));
  for(i=0; i<N; ++i)
    xp[i]=x[i];
  for(i=N; i<Np; ++i)
    xp[i]=0;
  depth=(int) (log(Np)/log(2))/2;

  p2d=(int) pow(2,depth);
  N_over_p2d=Np/p2d;
  coeffs=(complex float*) malloc(sizeof(complex float)*Np/pow(2,depth));
  for(k=0; k<N; ++k){
    z=-2*M_PI*I*k/N;
    for(i=0; i<Np/pow(2,depth); ++i){
      coeffs[i]=cexp(z*i*p2d);
    }
    
    s=0;
    for(i=0; i<p2d; ++i){
      s1=0;
      for(j=0; j<N_over_p2d; ++j){
	s1+=coeffs[j]*xp[ (int) (j*p2d+i)];
      }
      s+=s1*cexp(z*i);
    }
    Xp[k]=s;
  }
  
  for(i=0; i<N; ++i)
    X[i]=Xp[i];
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
}

