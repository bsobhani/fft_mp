#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "complex.h"

void fft(complex float* x, complex float* X, int N);
void ifft(complex float* X, complex float* x, int N);

#define NSIG 15
#define OUTPUTS_ON 1

int main(void){
  complex float x[160000];
  complex float X[160000];
  int k;
  /*x[0]=1; x[1]=2; x[2]=3; x[3]=4; x[4]=5; x[5]=6; x[6]=7; x[7]=8;
  x[8]=1; x[9]=2; x[10]=3; x[11]=4; x[12]=5; x[13]=6; x[14]=7; x[15]=8;
  */

  for(k=0; k<NSIG; ++k){
    x[k]=k;
  }
  
  printf("Original signal:\n");

  if(OUTPUTS_ON)
  for(k=0; k<NSIG; ++k)
    printf("%f + %fj\n",creal(x[k]),cimag(x[k]));
  
  fft(x,X,NSIG); /* FFT */

  printf("\nFFT of signal:\n");

  if(OUTPUTS_ON)
  for(k=0; k<NSIG; ++k)
     printf("%f + %fj\n",creal(X[k]),cimag(X[k]));
  
  ifft(X,x,NSIG); /* Inverse FFT */
  
  printf("\nInverse FFT of FFT of signal:\n");
  if(OUTPUTS_ON)
  for(k=0; k<NSIG; ++k)
     printf("%f + %fj\n",creal(x[k]),cimag(x[k]));
  return 0;
}
