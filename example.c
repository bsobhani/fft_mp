#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "complex.h"

void fft(complex float* x, complex float* X, int N);
void ifft(complex float* X, complex float* x, int N);


int main(void){
  complex float x[5];
  complex float X[5];
  int k;
  x[0]=1; x[1]=2; x[2]=3; x[3]=4; x[4]=5;

  printf("Original signal:\n");

  for(k=0; k<5; ++k)
    printf("%f + %fj\n",creal(x[k]),cimag(x[k]));
  
  fft(x,X,5); /* FFT */

  printf("\nFFT of signal:\n");

  for(k=0; k<5; ++k)
    printf("%f + %fj\n",creal(X[k]),cimag(X[k]));
  
  ifft(X,x,5); /* Inverse FFT */
  
  printf("\nInverse FFT of FFT of signal:\n");
  for(k=0; k<5; ++k)
    printf("%f + %fj\n",creal(x[k]),cimag(x[k]));
  return 0;
}
