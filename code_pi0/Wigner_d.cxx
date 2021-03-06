#include "TROOT.h"
#include <cmath>
#include <iostream>
#include <complex>
#include <string>
#include <cctype>
#include <vector>
#include <map>

 using namespace std; 

#define MAX(x,y) (x>y ? x : y)
#define MIN(x,y) (x<y ? x : y)

inline Int_t factorial(Int_t __i) {
  Int_t f = 1;
  if((__i == 0)||(__i == 1)) f = 1;  
  else{
    while(__i > 0){
      f = f*__i;
      __i--;
    }
  }
  return f;
} 

Float_t Wigner_d(const Float_t __j,const Float_t __m,const Float_t __n,
		Float_t __beta){

  Int_t J = (Int_t)(2.*__j);
  Int_t M = (Int_t)(2.*__m);
  Int_t N = (Int_t)(2.*__n);
  Int_t temp_M, k, k_low, k_hi;
  Float_t const_term = 0.0, sum_term = 0.0, d = 1.0;
  Int_t m_p_n, j_p_m, j_p_n, j_m_m, j_m_n;
  Int_t kmn1, kmn2, jmnk, jmk, jnk;
  Float_t kk;

  if (J < 0 || abs (M) > J || abs (N) > J) {
    cerr << endl;
    cerr << "d: you have entered an illegal number for J, M, N." << endl;
    cerr << "Must follow these rules: J >= 0, abs(M) <= J, and abs(N) <= J." 
	 << endl;
    cerr << "J = " << J <<  " M = " << M <<  " N = " << N << endl;
    return 0.;
  }
  
  if (__beta < 0) {
    __beta = fabs (__beta);
    temp_M = M;
    M = N;
    N = temp_M;
  }

  m_p_n = (M + N) / 2;
  j_p_m = (J + M) / 2;
  j_m_m = (J - M) / 2;
  j_p_n = (J + N) / 2;
  j_m_n = (J - N) / 2;
  
  kk = (Float_t)factorial(j_p_m)*(Float_t)factorial(j_m_m)
    *(Float_t)factorial(j_p_n) * (Float_t)factorial(j_m_n) ;
  const_term = pow((-1.0),(j_p_m)) * sqrt(kk);	
  
  k_low = MAX(0, m_p_n);
  k_hi = MIN(j_p_m, j_p_n);

  for (k = k_low; k <= k_hi; k++) {
    
    kmn1 = 2 * k - (M + N) / 2;
    jmnk = J + (M + N) / 2 - 2 * k;
    jmk = (J + M) / 2 - k;
    jnk = (J + N) / 2 - k;
    kmn2 = k - (M + N) / 2;
	
    sum_term += pow ((-1.0), (k)) *
      ((pow (cos (__beta / 2.0), kmn1)) * (pow (sin (__beta / 2.0), jmnk))) /
      (factorial (k) * factorial (jmk) * factorial (jnk) * factorial (kmn2));
  }

  d = const_term * sum_term;
  return d;
}
