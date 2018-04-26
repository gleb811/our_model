#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"
#include "res_param.h"
#include "Math/SpecFunc.h"
#include "TF1.h"

using namespace std;











Float_t DSIMPS(Float_t *Y, Float_t *X, Int_t Npoints){
Float_t Integral, h;
Integral = 0.;
for(Int_t i=2;i<Npoints;i += 2){

h = X[i] - X[i-2];
Integral += (Y[i-2] + 4.*Y[i-1] + Y[i])*h; 
};

Integral /= 6.0;

return Integral;
};



Bool_t rho_prot_orb_mom(Float_t res_j, Float_t res_parity,Int_t &max_val) {

/*
Int_t max_val_tmp = Int_t(res_j + 1. + 0.5);
if ((max_val_tmp%2 == 0) && (res_parity < 0.)) max_val = Int_t(max_val_tmp);
if ((max_val_tmp%2 == 1) && (res_parity < 0.)) max_val = Int_t(max_val_tmp-1);
if ((max_val_tmp%2 == 0) && (res_parity > 0.)) max_val = Int_t(max_val_tmp-1);
if ((max_val_tmp%2 == 1) && (res_parity > 0.)) max_val = Int_t(max_val_tmp);

*/

// An alternative method to calculate max_val

max_val = Int_t(res_j + 1. + 0.5) - Int_t(res_j + res_parity/2.)%2;

if (res_parity < 0.) return false;
if (res_parity > 0.) return true;
};



Float_t avrg_rho_mom (Int_t res_num) {

Float_t A,B;
Float_t P_DISTR[27],SRHO_run[27];
Float_t ERHO_run,MRHO_run,PRHO_run;


MP = 0.938;
MRHO = 0.769;
MPIP = 0.138;


A = (2.*MPIP)*(2.*MPIP);
B= (MRHO+4.*WIDTH_RHO)*(MRHO+4.*WIDTH_RHO);

for (Int_t i=0; i<=26; i++) {
SRHO_run[i] = A + (B-A)*i/26.;
ERHO_run = (res_mass[res_num]*res_mass[res_num]+SRHO_run[i]-MP*MP)/(2.*res_mass[res_num]);
MRHO_run = sqrt(SRHO_run[i]);

if (ERHO_run > MRHO_run) {
PRHO_run = sqrt(ERHO_run*ERHO_run-SRHO_run[i]);
P_DISTR[i] = PRHO_run/M_PI*WIDTH_RHO*MRHO/((SRHO_run[i]-MRHO*MRHO)*(SRHO_run[i]-MRHO*MRHO)+MRHO*MRHO*WIDTH_RHO*WIDTH_RHO);
} else {
P_DISTR[i]=0.;
};


}; 

return DSIMPS(P_DISTR,SRHO_run,27);

};



Float_t avrg_rho_mom_no_res (Float_t W) {

Float_t A,B;
Float_t P_DISTR[27],SRHO_run[27];
Float_t ERHO_run,MRHO_run,PRHO_run;


MP = 0.938;
MRHO = 0.769;
MPIP = 0.138;


A = (2.*MPIP)*(2.*MPIP);
B= (MRHO+4.*WIDTH_RHO)*(MRHO+4.*WIDTH_RHO);

for (Int_t i=0; i<=26; i++) {
SRHO_run[i] = A + (B-A)*i/26.;
ERHO_run = (W*W+SRHO_run[i]-MP*MP)/(2.*W);
MRHO_run = sqrt(SRHO_run[i]);

if (ERHO_run > MRHO_run) {
PRHO_run = sqrt(ERHO_run*ERHO_run-SRHO_run[i]);
P_DISTR[i] = PRHO_run/M_PI*WIDTH_RHO*MRHO/((SRHO_run[i]-MRHO*MRHO)*(SRHO_run[i]-MRHO*MRHO)+MRHO*MRHO*WIDTH_RHO*WIDTH_RHO);
} else {
P_DISTR[i]=0.;
};


}; 

return DSIMPS(P_DISTR,SRHO_run,27);

};





void corr_centrifug_bar_rho(Float_t W, TLorentzVector P4_PIP, TLorentzVector P4_PIM) {

for (Int_t res_num=0;res_num<Number_of_res;res_num++) {
for (Int_t i=0;i<3;i++) {
for (Int_t j=0;j<2;j++) {
factor_centrifugal[res_num][i][j] = 0.;
};
};
};

//cout << "avrg_rho_mom = " << avrg_rho_mom(1);

TF1 *Bessel = new TF1("bessel0","ROOT::Math::cyl_bessel_j([0],x)",0.,10.);
TF1 *Neumann = new TF1("neuman0","ROOT::Math::cyl_neumann([0],x)",0.,10.);

Int_t max;

for (Int_t res_num=0;res_num<Number_of_res;res_num++) {

for (Int_t RHO_NU_MOM = rho_prot_orb_mom(res_j[res_num],res_parity[res_num],max); RHO_NU_MOM <= max;  RHO_NU_MOM+=2)  {

TLorentzVector P4_RHO;
Float_t R_INTERACTION = 10.102;
Float_t Bess_rho_mom_res,Bess_rho_mom;
Float_t Neum_rho_mom_res,Neum_rho_mom;

Float_t factor;
Float_t tmp;

P4_RHO = P4_PIP + P4_PIM;
// qc = P4_RHO.Vect().Mag()
//cout << P4_RHO.Vect().Mag() << endl;

Float_t E_RHO_res,P_RHO_mag_res;

E_RHO_res = (res_mass[res_num]*res_mass[res_num] + P4_RHO.Mag2() - MP*MP)/2./res_mass[res_num];

if ((E_RHO_res*E_RHO_res - P4_RHO.Mag2()) > 0.) {
P_RHO_mag_res = sqrt(E_RHO_res*E_RHO_res - P4_RHO.Mag2());

// qr = P_RHO_mag_res
//cout << E_RHO_res*E_RHO_res - P4_RHO.Mag2() << endl;

//TF1 *Bessel = new TF1("bessel0","ROOT::Math::cyl_bessel_j([0],x)",0.,10.);
Bessel->SetParameters(RHO_NU_MOM+0.5,0);


//TF1 *Neumann = new TF1("neuman0","ROOT::Math::cyl_neumann([0],x)",0.,10.);

Neumann->SetParameters(RHO_NU_MOM+0.5,0);


Bess_rho_mom = (Bessel->Eval((P4_RHO.Vect().Mag())*R_INTERACTION));
Neum_rho_mom = (Neumann->Eval((P4_RHO.Vect().Mag())*R_INTERACTION));
Bess_rho_mom_res = (Bessel->Eval(P_RHO_mag_res*R_INTERACTION));
Neum_rho_mom_res = (Neumann->Eval(P_RHO_mag_res*R_INTERACTION));

factor = (Bess_rho_mom_res*Bess_rho_mom_res + Neum_rho_mom_res*Neum_rho_mom_res);

factor = factor/(Bess_rho_mom*Bess_rho_mom + Neum_rho_mom*Neum_rho_mom);

factor = res_mass[res_num]*factor/W;

if (abs(factor)>1.5)factor = factor/abs(factor)*1.5; // Taken from bw_correction1_rhon.f

//if (factor > 10.) cout << "factor = " << factor << "P4_RHOs = " << P_RHO_mag_res*R_INTERACTION << endl;
//factor = 1.;

} else {

factor = 0.;

};
//cout << res_name[res_num] << "  " << RHO_NU_MOM << "  " << (Int_t(RHO_NU_MOM-res_parity[res_num]))/2 << "  " << factor << endl;

factor = 1.;

factor_centrifugal[res_num][(RHO_NU_MOM-(max%2))/2][0] = factor;
factor_centrifugal[res_num][(RHO_NU_MOM-(max%2))/2][1] = factor;

/*
BR_RHO_LS[res_num][(RHO_NU_MOM-(max%2))/2][0] = BR_RHO_LS[res_num][(RHO_NU_MOM-(max%2))/2][0]*factor;
BR_RHO_LS[res_num][(RHO_NU_MOM-(max%2))/2][1] = BR_RHO_LS[res_num][(RHO_NU_MOM-(max%2))/2][1]*factor;
*/


//cout << res_name[res_num] << "  " << RHO_NU_MOM << "  " << (Int_t(RHO_NU_MOM-res_parity[res_num]))/2 << "  " << BR_RHO_LS[res_num][(Int_t(RHO_NU_MOM-res_parity[res_num]))/2][0] << "  " << BR_RHO_LS[res_num][(Int_t(RHO_NU_MOM-res_parity[res_num]))/2][1] << endl;

}; //end loop over RHO_NU_MOM
}; // end loop over resonances



};





