#include <stdlib.h>
#include <stdio.h>

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
    printf("RAND Error: x=%f, i=%ld, il=%ld\n",x,i,il);
  }
  return x;
}

void srand_(long *iseed)
{
  srand( (unsigned int) *iseed); 
}

/*******************************************************************/
float random1_()
{
  long i;
  float x;
  static double xfac;
  static long il=0;
  if(!il){
    il=(long)RAND_MAX;
    xfac=1./( ((double) il) + (((double) il)/1000000.));
  }
  i = random(); 
  x = ((double) i)*xfac ;
  if(x >= 1. || x < 0.) {
    printf("RAND Error: x=%f, i=%ld, il=%ld\n",x,i,il);
  }
  return x;
}

/* Fortran call to this should pass state as integer*1 and statelen as
   its length (or integer*2 and half its length)*/
int initrand_(int *seed, char *state, int *statelen)
{
  if(initstate((unsigned int) *seed,state,*statelen)==NULL){
    printf("Initstate Error: seed=%d\n",*seed);
  }
}

int setrand_(char *state){
  setstate(state);
}
