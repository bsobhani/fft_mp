#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "complex.h"

/*
void dft(complex float* x, complex float* X, int N){
  int n,k;
  complex float s,z;
  z=0;
  for(k=0; k<N; ++k){
    //fprintf(stderr,"%d %d\n",k,N);
    s=0;
    z=-2*M_PI*I*k/N;
    for(n=0; n<N; ++n){
      s+=x[n]*cexp(z*n);
    }
    X[k]=s;
  }
  return;
}
*/

void fft_depth(complex float* x, complex float* X, int N, int depth){
  complex float* coeffs;
  complex float s;
  complex float s1;
  complex float z;
  int i,j,k;
  int p2d,N_over_p2d;
  p2d=(int) pow(2,depth);
  N_over_p2d=N/p2d;
  coeffs=(complex float*) malloc(sizeof(complex float)*N/pow(2,depth));
  for(k=0; k<N; ++k){
    //fprintf(stderr,"%d %d\n",k,N);
    z=-2*M_PI*I*k/N;
    //z=-2*M_PI*I*k*p2d/N;
    for(i=0; i<N/pow(2,depth); ++i){
      //coeffs[i]=cexp(z*i*pow(2,depth));
      coeffs[i]=cexp(z*i*p2d);
    }
    
    s=0;
    for(i=0; i<p2d; ++i){
      s1=0;
      for(j=0; j<N_over_p2d; ++j){
	//s1+=coeffs[j]*x[ (int) (j*pow(2,depth)+i)];
	//s1+=cexp(z*j*p2d)*x[ (int) (j*p2d+i)];
	s1+=coeffs[j]*x[ (int) (j*p2d+i)];
      }
      s+=s1*cexp(z*i);
    }
    X[k]=s;
  }
  
  return;
}

void fft(complex float* x, complex float* X, int N){
  int depth=(int) (log(N)/log(2))/2;
  fft_depth(x, X, N, depth);
}

void ifft(complex float* x, complex float* X, int N){
  complex float* xr;
  int depth=(int) (log(N)/log(2))/2;
  int i;
  xr=(complex float*) malloc(sizeof(complex float)*N);
  fft_depth(X, xr, N, depth);
  for(i=0; i<N; ++i){
    x[i]=xr[N-i-1];
  }
  for(i=0; i<N; ++i){
    x[i]=x[i]/N;
  }
}



int test(){
  complex float buffer[8];
  complex float out[8];
  int i;
  buffer[0]=1; buffer[1]=2; buffer[2]=3; buffer[3]=4; buffer[4]=5; buffer[5]=6; buffer[6]=7; buffer[7]=8;

  //dft(buffer,out,8);
  
  fft_depth(buffer,out,8,2);

  
  for(i=0; i<8; ++i){
    printf("%f+j%f\n",creal(out[i]),cimag(out[i]));
  }
  return 0;
}
