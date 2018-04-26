#include "TROOT.h"
#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <iostream>
#include "global.h"
#include "res_param.h"
#include "TComplex.h"
#include "TParticle.h"
#include "corr_centrifug_bar_rho.h"

void amp_elmag(TComplex *amp_gam_nu_to_res,Float_t W,Float_t Q2, TLorentzVector P4_PIP, TLorentzVector P4_PIM){

TComplex *amp_gam_nu_to_res_pointer;
Float_t HGAM[3]={-1.,0.,1.};
Float_t HNUI[2]={0.5,-0.5};

Float_t E_gamma_res_cms,q_gamma_res_cms;
Float_t E_gamma_cms,q_gamma_cms;
Float_t symfact;


E_gamma_cms = (W*W-Q2-MP*MP)/2./W;

q_gamma_cms = sqrt(E_gamma_cms*E_gamma_cms+Q2);

//cout << "aaaa " << (P4_PIP + P4_PIM).Vect().Mag() << "   " << avrg_rho_mom(1) << endl;

for (Int_t res_num = 0; res_num < Number_of_res; res_num++)  {

E_gamma_res_cms = (res_mass[res_num]*res_mass[res_num]-Q2-MP*MP)/2./res_mass[res_num];

q_gamma_res_cms = sqrt(E_gamma_res_cms*E_gamma_res_cms+Q2);

symfact = -1.*res_parity[res_num]*pow(-1.,res_j[res_num]-1.5);

for (Int_t i_hnui = 0; i_hnui < 2; i_hnui++) {
for (Int_t i_hgam = 0; i_hgam < 3; i_hgam++) {


//if(abs(HNUI[i_hnui]-HGAM[i_hgam]) <= res_j[res_num]) { 

amp_gam_nu_to_res_pointer = amp_gam_nu_to_res +res_num*3*2+ i_hnui*3 + i_hgam;


if ((i_hnui == 1)&&(i_hgam == 2)) {
amp_gam_nu_to_res_pointer[0] = TComplex(Double_t(A32[res_num]));
};
if ((i_hnui == 0)&&(i_hgam == 2)) {
amp_gam_nu_to_res_pointer[0] = TComplex(Double_t(A12[res_num]));
};
if ((i_hnui == 1)&&(i_hgam == 1)) {
amp_gam_nu_to_res_pointer[0] = TComplex(Double_t(S12[res_num]*sqrt(2.)));
};

if ((i_hnui == 0)&&(i_hgam == 0)) {
amp_gam_nu_to_res_pointer[0] = TComplex(Double_t(symfact*A32[res_num]));
};
if ((i_hnui == 1)&&(i_hgam == 0)) {
amp_gam_nu_to_res_pointer[0] = TComplex(Double_t(symfact*A12[res_num]));
};
if ((i_hnui == 0)&&(i_hgam == 1)) {
amp_gam_nu_to_res_pointer[0] = TComplex(Double_t(symfact*S12[res_num]*sqrt(2.)));
};


amp_gam_nu_to_res_pointer[0] =amp_gam_nu_to_res_pointer[0]* Double_t(W/res_mass[res_num]*sqrt((8.*MP*res_mass[res_num]*q_gamma_res_cms)/(4.*M_PI*alph_const)));





amp_gam_nu_to_res_pointer[0] = amp_gam_nu_to_res_pointer[0]*Double_t(sqrt(q_gamma_res_cms/q_gamma_cms)/1000.);

amp_gam_nu_to_res_pointer[0] = amp_gam_nu_to_res_pointer[0]*Double_t(sqrt(avrg_rho_mom(res_num)/((P4_PIP + P4_PIM).Vect().Mag())));

//cout << " QQQ1 = " << ((P4_PIP + P4_PIM).Vect().Mag()) << endl;
//} else {
//amp_gam_nu_to_res_pointer[0] = TComplex(0.,0.);
//};




};
};
};


};
