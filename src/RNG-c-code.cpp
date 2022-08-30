#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <R_ext/Random.h>
using namespace Rcpp;

#define K_MAX	50000
#define PP 2147483647  //2^31-1
#define IPP 4.6566129e-10
#define I32 2.328306437e-10;
#define ROTR(x,n) ((((unsigned long) x)>>(n))|(((unsigned long) x)<<(32-(n))))

#define TSIZE 256 		// actual table size used
#define maxT 1024 		//maximum table size
#define maxB 15
#define B_X 32896  		
#define B_Y 32776
#define T_MOD TSIZE-1
#define MOD2_32 4294967295  

static long long B_X1 = 536869888; // will change based on K, S- can be done more efficiently
static long long B_X2 = 65011712;
static long long B_X3 = 67633152;
static long long B_X4 = 67108736;

static int I_X;           /* running index */

unsigned long MODP(unsigned long z) {return ((((z)&PP)+((z)>>31)) &PP);}

static int K_X = 1597;
static int S_X = 1;
static int nseed = 3;
static int seed[K_MAX]={K_X,S_X};
static int XX[K_MAX];
static double res;
int set=0;
int number;

int i, K12, K13, K23;


int get_random(){
  number = (2147483629*number + 2147483587) & PP;
  return number;
}

void dx_1(){
  K_X = seed[0];
  S_X = seed[1];
  int II0 = I_X;
  if(++I_X >= K_X)  I_X = 0;     /*wrap around running index */
  XX[I_X] = MODP(B_X1 * XX[I_X] + XX[II0]);
  res = (double) XX[I_X] * IPP;
}

void dx_2(){
  K_X = seed[0];
  S_X = seed[1];
  int II0 = I_X;
  if(++I_X >= K_X)  I_X = 0;     /*wrap around running index */
  XX[I_X] = MODP(B_X2 * (XX[I_X] + XX[II0]));
  res = (double) XX[I_X] * IPP;
}

void dx_3(){
  K_X = seed[0];
  S_X = seed[1];
  int II0 = I_X;
  if(++I_X >= K_X)  I_X = 0;     /*wrap around running index */
  if(++K12 >= K_X) K12 = 0;    /*wrap around K12*/
  XX[I_X] = MODP(B_X3 * (XX[I_X] + XX[K12] + XX[II0]));
  res= (double) XX[I_X] * IPP;
}

void dx_4(){
  K_X = seed[0];
  S_X = seed[1];
  int II0 = I_X;
  if(++I_X >= K_X)   I_X = 0;    /*wrap around running index */
  if(++K13 >= K_X) K13 = 0;    /*wrap around K13*/
  if(++K23 >= K_X) K23 = 0;    /*wrap around K23*/
  XX[I_X] = MODP(B_X4 * (XX[I_X]+XX[K13]+XX[K23]+XX[II0]));
  res= (double) XX[I_X] * IPP;
}

void dx_init(){
  for (i=0; i< K_MAX; i++){
    XX[i] = get_random();
  }
}

void (*dx_gen)();

void generator_type(){
  switch(S_X){
  case 1: dx_gen=&dx_1;
    dx_init();
    break;
  case 2: dx_gen=&dx_2;
    dx_init();
    break;
  case 3: dx_gen=&dx_3;
    dx_init();
    break;
  case 4: dx_gen=&dx_4;
    dx_init();
    break;
  default: dx_gen=&dx_1;
  dx_init();
  S_X=1;
  seed[1]=S_X;
  Rcerr << "The value of 'S' that was chosen is not compatible with this package. \nBy default, a DX-" << K_X<< "-1 generator was chosen.\n";
  break;
  }
}

double * user_unif_rand(){
  dx_gen();
  return &res;
}

void user_unif_init(Int32 seed_in) {
  K_X = seed[0];
  S_X = seed[1];

  number=seed_in;
  seed[2] =  get_random();
  
  if(K_X<5){
    K_X=1597;
    seed[0]=K_X;
    Rcerr << "The value of 'K' that was chosen is not compatible with this package. \nBy default, a DX-1597-" << S_X<< " generator was chosen.\n";
  }
  
  if(K_X>50000){
    K_X=1597;
    seed[0]=K_X;
    Rcerr << "The value of 'K' that was chosen is not compatible with this package. \nBy default, a DX-1597-" << S_X<< " generator was chosen.\n";
  }
  generator_type();
  
  I_X = seed[0]-1;
  K12 = seed[0]/2-1;
  K13 = 2*seed[0]/3-1;   
  K23 = seed[0]/3-1;
  
  if(set!=0){
    Rcout << "# You are currently using a DX-" << K_X << "-" << S_X << " Generator.\n# To alter it use the dx_init() function.\n";
  }
  set=1;
}

int * user_unif_nseed() {return &nseed; }
int * user_unif_seedloc() { return (int *) &seed; }
