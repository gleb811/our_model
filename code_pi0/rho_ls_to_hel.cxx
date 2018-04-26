#include "TROOT.h"
#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <iostream>
#include "global.h"
#include "res_param.h"
#include "corr_centrifug_bar_rho.h"
#include "Clebsch.h"
#include "Wigner_d.h"
#include "TComplex.h"
#include "TParticle.h"

using namespace std;

void rho_ls_to_hel(TComplex *amp_dec_to_rho_p,TLorentzVector P4_proton,Float_t M12){

Int_t max;

TComplex *amp_dec_to_rho_p_pointer;

Float_t S[2]={0.5,1.5};
Float_t HRO[3]={-1.,0.,1.};
Float_t HNU[2]={0.5,-0.5};
Float_t HRE[4]={-1.5,-0.5,0.5,1.5};

Float_t ERHO_res;


//for (Int_t i_hre = 0; i_hre < 4; i_hre++) {
//HRE[i_hre] = i_hre - (max_proj_j-1.)/2.;

//cout << HRE[i_hre] << endl;
//};

Float_t PHI_proton = P4_proton.Phi();
if (PHI_proton<0.)PHI_proton = 2.*M_PI+PHI_proton;

//cout << "THETA PROTON = " << P4_proton.Theta() << " PHI PROTON = " << PHI_proton << endl;

Double_t A,B;

for (Int_t res_num=0;res_num<Number_of_res;res_num++) {
for (Int_t i_s = 0; i_s < 2; i_s++) {
for (Int_t i_hro = 0; i_hro < 3; i_hro++) {
for (Int_t i_hnu = 0; i_hnu < 2; i_hnu++) {
rho_width_hel_sqrt[res_num][i_hro][i_hnu] = 0.;
};
};
};
};



for (Int_t res_num=0;res_num<Number_of_res;res_num++) {

for (Int_t RHO_NU_MOM = rho_prot_orb_mom(res_j[res_num],res_parity[res_num],max); RHO_NU_MOM <= max;  RHO_NU_MOM+=2)  {

for (Int_t i_s = 0; i_s < 2; i_s++) {
//for (Int_t i_s = 0; i_s <= Int_t(abs(max/2.-RHO_NU_MOM/2.-1.)); i_s++) {
for (Int_t i_hro = 0; i_hro < 3; i_hro++) {
for (Int_t i_hnu = 0; i_hnu < 2; i_hnu++) {

Clebsch(RHO_NU_MOM,S[i_s],res_j[res_num],0.,HNU[i_hnu]-HRO[i_hro],HNU[i_hnu]-HRO[i_hro],A);
Clebsch(0.5,1.,S[i_s],HNU[i_hnu],-HRO[i_hro],HNU[i_hnu]-HRO[i_hro],B);


rho_width_hel_sqrt[res_num][i_hro][i_hnu] = rho_width_hel_sqrt[res_num][i_hro][i_hnu] + A*B*sqrt((2.*RHO_NU_MOM+1.)/(2.*res_j[res_num]+1.))*sqrt(res_width[res_num]*BR_RHO_LS[res_num][(RHO_NU_MOM-(max%2))/2][i_s]);
//if(res_num == 5)cout << RHO_NU_MOM << "  " << res_j[res_num] << "  " << A*B << "  " << sqrt(res_width[res_num]*BR_RHO_LS[res_num][(RHO_NU_MOM-(max%2))/2][i_s]) << "  " << rho_width_hel_sqrt[res_num][i_hro][i_hnu] << endl;
};// end of loop over i_hnu (helicity of the nucleon)
}; // end of loop over i_hro (helicity of the rho)
}; // end of loop over s (the total spin of the rho and nucleon)

}; //end loop over RHO_NU_MOM
}; // end loop over resonances



for (Int_t res_num = 0; res_num < Number_of_res; res_num++) {

// Energy of the rho meson at the resonant point

ERHO_res = (res_mass[res_num]*res_mass[res_num]+M12*M12-MN*MN)/(2.*res_mass[res_num]);

for (Int_t i_hnu = 0; i_hnu < 2; i_hnu++) {
for (Int_t i_hro = 0; i_hro < 3; i_hro++) {

//amp_dec_to_rho_p_pointer = amp_dec_to_rho_p + res_num*3*2 + i_hnu*3 + i_hro;

for (Int_t i_hre = 0; i_hre < 4; i_hre++) {


amp_dec_to_rho_p_pointer = amp_dec_to_rho_p + res_num*3*2*4+ i_hnu*3*4 + i_hro*4 + i_hre;

/*
if ((res_num == 11)&&(i_hnu == 0)&&(i_hro == 1)) {
amp_dec_to_rho_p_pointer[0] = TComplex(1.,2.); 

} else {

amp_dec_to_rho_p_pointer[0] = TComplex(3.,4.); 

};
*/
// The condition (ERHO_res*ERHO_res-M12*M12)>0. checks if the production of the rho meson via the decay of given resonance is possible

if ((abs(HNU[i_hnu]-HRO[i_hro]) <= res_j[res_num])&&(abs(HRE[i_hre]) <= res_j[res_num])&&((ERHO_res*ERHO_res-M12*M12)>0.)) {

//if (res_num == 1) {
//cout << "res_num_j = " << res_j[res_num] << " HNU = " << HNU[i_hnu] << " hro = " << HRO[i_hro] << endl;
//};

//if (i_hre < Int_t(2.*res_j[res_num]+1.)) {
amp_dec_to_rho_p_pointer[0] = TComplex(rho_width_hel_sqrt[res_num][i_hro][i_hnu],0.);
amp_dec_to_rho_p_pointer[0] = amp_dec_to_rho_p_pointer[0]*Double_t(res_mass[res_num]*sqrt(8.*M_PI*(2.*res_j[res_num]+1.)));
amp_dec_to_rho_p_pointer[0] = amp_dec_to_rho_p_pointer[0]*Double_t(sqrt(1./avrg_rho_mom(res_num)));
amp_dec_to_rho_p_pointer[0] = amp_dec_to_rho_p_pointer[0]*Double_t(Wigner_d(res_j[res_num],HRE[i_hre],HNU[i_hnu]-HRO[i_hro],P4_proton.Theta()));


amp_dec_to_rho_p_pointer[0] = amp_dec_to_rho_p_pointer[0]*TComplex::Exp(TComplex(0.,1.)*Double_t(HRE[i_hre]*PHI_proton));


//if (res_num == 2) { cout << "res_j[res_num] = " << res_j[res_num] << " HRE = " << HRE[i_hre]<< "  " <<  Double_t(Wigner_d(res_j[res_num],HRE[i_hre],HNU[i_hnu]-HRO[i_hro],P4_proton.Theta())) << "   " <<  TComplex::Exp(TComplex(0.,1.)*Double_t(HRE[i_hre+Int_t((max_proj_j -(2.*res_j[res_num]+1.))/2.)]*PHI_proton)) << " res_num = " << res_num   << endl; };
//cout << P4_proton.Theta() << "   " << P4_proton.Phi() << endl;
} else {
amp_dec_to_rho_p_pointer[0] = TComplex(0.,0.);
};


};
//cout << "   " << endl;
//cout <<  i_hro << "  " << i_hnu << "  " << rho_width_hel_sqrt[res_num][i_hro][i_hnu] << endl;

};
};
};


//cout << Wigner_d(1.5,-1.5,-0.5,0.63596630096435547) << endl;
};
