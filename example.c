#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "complex.h"

void fft(complex float* x, complex float* X, int N);

int main(){
  complex float x[8];
  complex float X[8];
  int k;
  x[0]=1; x[1]=2; x[2]=3; x[3]=4; x[4]=5; x[5]=6; x[6]=7; x[7]=8;
  fft(x,X,8);
  for(k=0; k<8; ++k)
    printf("%f j%f\n",creal(X[k]),cimag(X[k]));
  return 0;
}
