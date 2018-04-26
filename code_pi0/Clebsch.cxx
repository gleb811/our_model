#include "TROOT.h"
#include <cmath>
#include <iostream>
#include <complex>
#include <string>
#include <cctype>
#include <vector>
#include <map>

/// Utility function used by Clebsch 
inline Double_t dfact(const Double_t &__x){
  if((__x < 0.00001) && (__x >= 0.0)) return 1.;
  if(__x < 0) return 0.;
  return __x*dfact(__x - 1.);
}

void Clebsch(Float_t  __j1,Float_t  __j2,Float_t  __J,Float_t  __m1,Float_t  __m2,Float_t  __M, Double_t &A){

/*8888
float  __j1 = 3.;
float  __m1 = 0.;
float  __j2 = 1.5;
float  __m2 = 0.5;
float  __J = 1.5;
float  __M = 0.5;
*/
//cout << __j1 << " " << __j2 << " " << __J << " " << __m1 << " " << __m2 << " " << __M << endl;



  // convert to pure integers (each 2*spin)
  Int_t j1 = (Int_t)(2.*__j1);
  Int_t m1 = (Int_t)(2.*__m1);
  Int_t j2 = (Int_t)(2.*__j2);
  Int_t m2 = (Int_t)(2.*__m2);
  Int_t J = (Int_t)(2.*__J);
  Int_t M = (Int_t)(2.*__M);

  Double_t n0,n1,n2,n3,n4,n5,d0,d1,d2,d3,d4,exp;
  Int_t nu = 0;
  
  Double_t sum = 0;
  while(((d3=(j1-j2-M)/2+nu) < 0)||((n2=(j1-m1)/2+nu) < 0 )) { nu++;}
  while (((d1=(J-j1+j2)/2-nu) >= 0) && ((d2=(J+M)/2-nu) >= 0) 
	 &&((n1=(j2+J+m1)/2-nu) >= 0 )){
    d3=((j1-j2-M)/2+nu);
    n2=((j1-m1)/2+nu);
    d0=dfact((Double_t) nu);
    exp=nu+(j2+m2)/2;
    n0 = (Double_t) pow(-1.,exp);
    sum += ((n0*dfact(n1)*dfact(n2))/(d0*dfact(d1)*dfact(d2)*dfact(d3)));
    nu++;
  }

//  if (sum == 0) return 0;

  n0 = J+1;
  n1 = dfact((Double_t) (J+j1-j2)/2);
  n2 = dfact((Double_t) (J-j1+j2)/2);
  n3 = dfact((Double_t) (j1+j2-J)/2);
  n4 = dfact((Double_t) (J+M)/2);
  n5 = dfact((J-M)/2);
  
  d0 = dfact((Double_t) (j1+j2+J)/2+1);
  d1 = dfact((Double_t) (j1-m1)/2);
  d2 = dfact((Double_t) (j1+m1)/2);
  d3 = dfact((Double_t) (j2-m2)/2);
  d4 = dfact((Double_t) (j2+m2)/2);
  
  A = ((Double_t) (n0*n1*n2*n3*n4*n5))/((Double_t) (d0*d1*d2*d3*d4));
  A = sqrt(A)*sum;  

  if((__m1 + __m2) != __M) A=0.;
  if(abs(__m1) > __j1) A=0;
  if(abs(__m2) > __j2) A=.0;
 
	
  	
}


