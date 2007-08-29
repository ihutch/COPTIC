#include <stdlib.h>
#include <stdio.h>
/*    int rand(void); 
    void srand(unsigned int seed);
      Here is a code fragment to convert the
      random ints to random floats in the range 0. to 1.:
*/

float rand_()
{
  long i;
  float x;
  static double xfac;
  static long il=0;
  if(!il){
    il=(long)RAND_MAX;
    xfac=1./( ((double) il) + (((double) il)/1000000.));
  }
  i = rand(); 
  x = ((double) i)*xfac ;
  if(x >= 1. || x < 0.) {
    printf("RAND Error: x=%f, i=%d, il=%d\n",x,i,il);
  }
  return x;
}

void srand_(long *iseed)
{
  srand( (unsigned int) *iseed); 
}

