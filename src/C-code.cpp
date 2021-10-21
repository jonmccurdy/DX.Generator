#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <R_ext/Random.h>
using namespace Rcpp;

#define K_MAX	3500
#define PP 2147483647  //2^31-1
#define IPP 4.6566129e-10
static long long B_X1 = 536869888; // will change based on K, S
static long long B_X2 = 65011712;
static long long B_X3 = 67633152;
static long long B_X4 = 67108736;

static int I_X;           /* running index */

unsigned long MODP(unsigned long z) {return ((((z)&PP)+((z)>>31)) &PP);}

static int K_X = 47;
static int S_X = 1;
static int nseed = K_X+2;
static Int32 seed[K_MAX];
static double res;
int set=0;

int i, K12, K13, K23;

// [[Rcpp::export]]
void changeK(int K, int S){
  K_X = K;
  S_X = S;
  // B_X1 change
  nseed = K_X+2;
  set=1;  /* determines if changeK has been used for printing message purpose*/
}

void dx_1(){
  int II0 = I_X;
  if(++I_X >= K_X+2)  I_X = 2;     /*wrap around running index */
  seed[I_X] = MODP(B_X1 * seed[I_X] + seed[II0]);
  res = (double) seed[I_X] * IPP;
}

void dx_2(){
  int II0 = I_X;
  if(++I_X >= K_X+2)  I_X = 2;     /*wrap around running index */
  seed[I_X] = MODP(B_X2 * (seed[I_X] + seed[II0]));
  res = (double) seed[I_X] * IPP;
}

void dx_3(){
  int II0 = I_X;
  if(++I_X >= K_X+2)  I_X = 2;     /*wrap around running index */
  if(++K12 >= K_X+2) K12 = 2;    /*wrap around K12*/
  seed[I_X] = MODP(B_X3 * (seed[I_X] + seed[K12] + seed[II0]));
  res= (double) seed[I_X] * IPP;
}

void dx_4(){
  int II0 = I_X;
  if(++I_X >= K_X+2)   I_X = 2;    /*wrap around running index */
  if(++K13 >= K_X+2) K13 = 2;    /*wrap around K13*/
  if(++K23 >= K_X+2) K23 = 2;    /*wrap around K23*/
  seed[I_X] = MODP(B_X4 * (seed[I_X]+seed[K13]+seed[K23]+seed[II0]));
  res= (double) seed[I_X] * IPP;
}

void (*dx_gen)();

void generator_type(){
  switch(S_X){
  case 1: dx_gen=&dx_1;
    break;
  case 2: dx_gen=&dx_2;
    break;
  case 3: dx_gen=&dx_3;
    break;
  case 4: dx_gen=&dx_4;
    break;
  default: dx_gen=&dx_1;
  Rcerr << "The value of 'S' that was chosen is not compatible with this package. \nBy default, a DX-" << K_X<< "-1 generator was chosen.\n";
  S_X=1;
  break;
  }
}

double * user_unif_rand(){
  dx_gen();
  return &res;
}

void user_unif_init(Int32 seed_in) { 
  generator_type();
  seed[0] = K_X;
  seed[1] = S_X;
  srand(seed_in);
  for (i=2; i< nseed; i++) seed[i] = rand() & PP;
  I_X = K_X+1;
  K12 = K_X/2+1;
  K13 = 2*K_X/3+1;   
  K23 = K_X/3+1;
  if(set==1){
    Rcout << "You are currently using a DX-" << K_X << "-" << S_X << " Generator. \nTo alter it use the set.dx() function.\n";
  }
}

int * user_unif_nseed() {return &nseed; }
int * user_unif_seedloc() { return (int *) &seed; }