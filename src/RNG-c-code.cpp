#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <R_ext/Random.h>
using namespace Rcpp;

#define K_MAX	35000
#define PP 2147483647  //2^31-1
#define IPP 4.6566129e-10
#define I32 2.328306437e-10;
#define ROTR(x,n) ((((unsigned long) x)>>(n))|(((unsigned long) x)<<(32-(n))))

#define P_mod 2147483647 //2^31 - 82845

static long long B_X1 = 2097101; //536869888; // will change based on K, S- can be done more efficiently
static long long B_X2 = 65011712;
static long long B_X3 = 67633152;
static long long B_X4 = 67108736;

static int I_X;           /* running index */

unsigned long MODP(unsigned long z) {return ((((z)&PP)+((z)>>31)) &PP);}

#define MUL20(x) ( ((x)>>11) + ((x)<<20)&PP ) 
#define MUL9(x) ( ((x)>>22) + ((x)<<9) &PP )

static int K_X = 9;
static int S_X = 1;
static int nseed = 3;
static int seed[K_MAX]={K_X,S_X};
static int XX[K_MAX];
static double res;
int set=0;
long long number;

int i, K12, K13, K23, K_X_old, S_X_old;
unsigned long S;

int get_random(){
  number = (1664525*number + 1013904223) % PP;
  return number;
}

void dx_1(){
  K_X = seed[0];
  S_X = seed[1];                  /*makes sure we stay with this S-value*/
  int II0 = I_X;
  if(++I_X >= K_X)  I_X = 0;     /*wrap around running index */
  XX[I_X] = MODP(B_X1 * XX[I_X] + XX[II0]);
  res = (double) XX[I_X] * IPP;
}

void dx_2(){
  K_X = seed[0];
  S_X = seed[1];                 /*makes sure we stay with this S-value*/
  int II0 = I_X;
  if(++I_X >= K_X)  I_X = 0;     /*wrap around running index */
  XX[I_X] = MODP(B_X2 * (XX[I_X] + XX[II0]));
  res = (double) XX[I_X] * IPP;
}

void dx_3(){
  K_X = seed[0];
  S_X = seed[1];               /*makes sure we stay with this S-value*/
  int II0 = I_X;
  if(++I_X >= K_X)  I_X = 0;     /*wrap around running index */
  if(++K12 >= K_X) K12 = 0;    /*wrap around K12*/
  XX[I_X] = MODP(B_X3 * (XX[I_X] + XX[K12] + XX[II0]));
  res= (double) XX[I_X] * IPP;
}

void dx_4(){
  K_X = seed[0];
  S_X = seed[1];                /*makes sure we stay with this S-value*/
  int II0 = I_X;
  if(++I_X >= K_X)   I_X = 0;    /*wrap around running index */
  if(++K13 >= K_X) K13 = 0;    /*wrap around K13*/
  if(++K23 >= K_X) K23 = 0;    /*wrap around K23*/
  XX[I_X] = MODP(B_X4 * (XX[I_X]+XX[K13]+XX[K23]+XX[II0]));
  res= (double) XX[I_X] * IPP;
}
/*
void dx_B2_perc(){
  K_X = seed[0];
  S_X = seed[1];                      //makes sure we stay with this S-value
  int II0 = I_X;
  if(++I_X >= K_X)  I_X = 0;     //wrap around running index
  S = MODP(XX[I_X] + XX[II0]); 
  XX[I_X] = MODP( MUL20(S) + MUL9(S));
  res = (double) XX[I_X] * IPP;
}

void dx_B2_MODP(){
  K_X = seed[0];
  S_X = seed[1];                      //makes sure we stay with this S-value
  int II0 = I_X;
  if(++I_X >= K_X)  I_X = 0;     //wrap around running index
  S = MODP(XX[I_X] + XX[II0]) ; 
  XX[I_X] = MODP(MUL20(S) + MUL9(S));
  res= (double) XX[I_X] * IPP;
}

void dx_B3_perc(){
  K_X = seed[0];
  S_X = seed[1]; //makes sure we stay with this S-value
  int II0 = I_X;
  if(++I_X >= K_X)  I_X = 0;     //wrap around running index 
  if(++K12 >= K_X) K12 = 0;    //wrap around K12
  S = (XX[I_X] + XX[K12] + XX[II0]) % PP ; 
  XX[I_X] = (MUL20(S) + MUL9(S)) % PP;
  res = (double) XX[I_X] * IPP;
}

void dx_B3_MODP(){
  K_X = seed[0];
  S_X = seed[1];                      //makes sure we stay with this S-value
  int II0 = I_X;
  if(++I_X >= K_X)  I_X = 0;     //wrap around running index 
  if(++K12 >= K_X) K12 = 0;    //wrap around K12  
  S = MODP(XX[I_X] + XX[K12] + XX[II0]) ; 
  XX[I_X] = MODP(MUL20(S) + MUL9(S));
  res= (double) XX[I_X] * IPP;
}


void dx_B4_perc(){
  K_X = seed[0];
  S_X = seed[1]; //makes sure we stay with this S-value
  int II0 = I_X;
  if(++I_X >= K_X)   I_X = 0;    //wrap around running index 
  if(++K13 >= K_X) K13 = 0;    //wrap around K13
  if(++K23 >= K_X) K23 = 0;    //wrap around K23
  S = (XX[I_X]+XX[K13]+XX[K23]+XX[II0]) % PP ; 
  XX[I_X] = (MUL20(S) + MUL9(S)) % PP;
  res = (double) XX[I_X] * IPP;
}

void dx_B4_MODP(){
  K_X = seed[0];
  S_X = seed[1];                      //makes sure we stay with this S-value
  int II0 = I_X;
  if(++I_X >= K_X)   I_X = 0;    //wrap around running index 
  if(++K13 >= K_X) K13 = 0;    //wrap around K13
  if(++K23 >= K_X) K23 = 0;    //wrap around K23
  S = MODP(XX[I_X]+XX[K13]+XX[K23]+XX[II0]) ; 
  XX[I_X] = MODP(MUL20(S) + MUL9(S));
  res= (double) XX[I_X] * IPP;
}

void dx_P1(){
  K_X = seed[0];
  S_X = seed[1];                  //makes sure we stay with this S-value
  int II0 = I_X;
  if(++I_X >= K_X)  I_X = 0;     //wrap around running index
  XX[I_X] = (B_X1 * XX[I_X] + XX[II0]) % P_mod;
  res = (double) XX[I_X] * IPP;
}

void dx_P2(){
  K_X = seed[0];
  S_X = seed[1];                 //makes sure we stay with this S-value
  int II0 = I_X;
  if(++I_X >= K_X)  I_X = 0;     //wrap around running index
  XX[I_X] = (B_X2 * (XX[I_X] + XX[II0])) % P_mod;
  res = (double) XX[I_X] * IPP;
}

void dx_P3(){
  K_X = seed[0];
  S_X = seed[1];               //makes sure we stay with this S-value
  int II0 = I_X;
  if(++I_X >= K_X)  I_X = 0;     //wrap around running index 
  if(++K12 >= K_X) K12 = 0;    //wrap around K12
  XX[I_X] = (B_X3 * (XX[I_X] + XX[K12] + XX[II0])) % P_mod;
  res= (double) XX[I_X] * IPP;
}

void dx_P4(){
  K_X = seed[0];
  S_X = seed[1];                //makes sure we stay with this S-value
  int II0 = I_X;
  if(++I_X >= K_X)   I_X = 0;    //wrap around running index 
  if(++K13 >= K_X) K13 = 0;    //wrap around K13
  if(++K23 >= K_X) K23 = 0;    //wrap around K23
  XX[I_X] = (B_X4 * (XX[I_X]+XX[K13]+XX[K23]+XX[II0])) % P_mod;
  res= (double) XX[I_X] * IPP;
}
*/

void dx_init(){
  for (i=0; i< K_X; i++){
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
  /*
  case 5: dx_gen=&dx_P1;
    dx_init();
    break;
  case 6: dx_gen=&dx_P2;
    dx_init();
    break;
  case 7: dx_gen=&dx_P3;
    dx_init();
    break;
  case 8: dx_gen=&dx_P4;
    dx_init();
    break;
  case 9: dx_gen=&dx_B2_perc;
    dx_init();
    break;
  case 10: dx_gen=&dx_B2_MODP;
    dx_init();
    break;
  case 11: dx_gen=&dx_B3_perc;
    dx_init();
    break;
  case 12: dx_gen=&dx_B3_MODP;
    dx_init();
    break;
  case 13: dx_gen=&dx_B4_perc;
    dx_init();
    break;
  case 14: dx_gen=&dx_B4_MODP;
    dx_init();
    break;
   */
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
  seed[2] = get_random();
  
  
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
  
  if(set==0||K_X_old != K_X || S_X_old != S_X){
     Rcout << "# You are currently using a DX-" << K_X << "-" << S_X << " Generator.\n# To alter it use the dx_init() function.\n";
   }
  K_X_old = K_X;
  S_X_old = S_X;
  set=1;
}

int * user_unif_nseed() {return &nseed; }
int * user_unif_seedloc() { return (int *) &seed; }
