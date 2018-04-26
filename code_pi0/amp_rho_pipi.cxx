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

void amp_rho_pipi(TComplex amp[], TLorentzVector P4_PIM, TLorentzVector P4_PIP) {

Float_t theta_rho = (P4_PIM+P4_PIP).Theta();
Float_t phi_rho = (P4_PIM+P4_PIP).Phi();

TComplex epsilon[3][4];
Double_t E_rho = (P4_PIM+P4_PIP)[3];
Double_t M_rho_run = (P4_PIM+P4_PIP).Mag();
Double_t px_rho = (P4_PIM+P4_PIP)[0];
Double_t py_rho = (P4_PIM+P4_PIP)[1];
Double_t pz_rho = (P4_PIM+P4_PIP)[2];

Double_t G = 5.724*1.056; // Three-level coupling constant
Double_t FF;

/*
TComplex epsilon_tmp[3][3];
TComplex epsilon_tmp2[3][4];

Double_t x = 1./(E_rho+M_rho_run);
TComplex epm;

epsilon_tmp[0][0] = TComplex(cos(theta_rho)*cos(phi_rho)/sqrt(2.),sin(phi_rho)/sqrt(2.));

epsilon_tmp[0][1] = TComplex(cos(theta_rho)*sin(phi_rho)/sqrt(2.),-1.*cos(phi_rho)/sqrt(2.));

epsilon_tmp[0][2] = TComplex(-1.*sin(theta_rho)/sqrt(2.),0.);

epsilon_tmp[1][0] = TComplex(sin(theta_rho)*cos(phi_rho),0.);

epsilon_tmp[1][1] = TComplex(sin(theta_rho)*sin(phi_rho),0.);

epsilon_tmp[1][2] = TComplex(cos(theta_rho),0.);

epsilon_tmp[2][0] = TComplex(-1.*cos(theta_rho)*cos(phi_rho)/sqrt(2.),sin(phi_rho)/sqrt(2.));

epsilon_tmp[2][1] = TComplex(-1.*cos(theta_rho)*sin(phi_rho)/sqrt(2.),-1.*cos(phi_rho)/sqrt(2.));

epsilon_tmp[2][2] = TComplex(sin(theta_rho)/sqrt(2.),0.);


for (Int_t i=0; i<3; i++){

epm = (epsilon_tmp[i][0]*px_rho +
       epsilon_tmp[i][1]*py_rho +
       epsilon_tmp[i][2]*pz_rho)/M_rho_run; 

epsilon_tmp2[i][0] = epsilon_tmp[i][0] + TComplex(px_rho*x,0.)*epm;
epsilon_tmp2[i][1] = epsilon_tmp[i][1] + TComplex(py_rho*x,0.)*epm;
epsilon_tmp2[i][2] = epsilon_tmp[i][2] + TComplex(pz_rho*x,0.)*epm;
epsilon_tmp2[i][3] = TComplex(epm,0.);
};
*/

//-------------------------------------------------------------------------

TLorentzVector eps_1,eps_2,eps_0;

/* 
eps_1 and eps_2 are the polarization vectors needed to construct eps_r and eps_l as
eps_r = -1.*sqrt(1./2.)*(eps_1 + i*eps_2) for \lambda_{rho} = 1
eps_l = 1.*sqrt(1./2.)*(eps_1 - i*eps_2) for \lambda_{rho} = -1

eps_0 is the polarization vector that corresponds \lambda_{rho} = 0
They are in the system with the z-axis along rho-momentum. 
*/

eps_1.SetXYZT(1.,0.,0.,0.);
eps_2.SetXYZT(0.,1.,0.,0.);
eps_0.SetXYZT(0.,0.,E_rho/M_rho_run,sqrt(px_rho*px_rho+py_rho*py_rho+pz_rho*pz_rho)/M_rho_run);

/*
Rotation of eps_1, eps_2, and eps_0 to the system with the z-axis along
the virtual photon.
*/

eps_1.RotateY(theta_rho);
eps_1.RotateZ(phi_rho);

eps_2.RotateY(theta_rho);
eps_2.RotateZ(phi_rho);

eps_0.RotateY(theta_rho);
eps_0.RotateZ(phi_rho);

/*
Construction of the polarization matrix from eps_l, eps_0, eps_2
*/

for (Int_t i=0; i<4; i++){
epsilon[0][i] = TComplex(eps_1[i],-1.*eps_2[i])*Double_t(1./sqrt(2.)); 
epsilon[1][i] = TComplex(eps_0[i],0.); 
epsilon[2][i] = TComplex(eps_1[i],eps_2[i])*Double_t(-1./sqrt(2.)); 
};




Double_t metrg[4] = {-1.,-1.,-1.,1.} ;



for (Int_t i_hro=0; i_hro<3; i_hro++){
amp[i_hro] = TComplex(0.,0.);
for (Int_t j=0; j<4; j++){
amp[i_hro] = amp[i_hro] + epsilon[i_hro][j]*Double_t(P4_PIP[j]-P4_PIM[j])*metrg[j];
};
//cout << "M = " << amp[i] << endl;
};


//Rho decay strong form fector
//      R.S.Longarce J.Dolbeau Nucl.Phys.B122 (1977) 493



Double_t p_pi_rho_rest_run;
Double_t p_pi_rho_rest;

p_pi_rho_rest_run = sqrt(M_rho_run*M_rho_run/4.-(MPIP+MPI0)*(MPIP+MPI0)/4.);
p_pi_rho_rest = sqrt(MRHO*MRHO/4.-(MPIP+MPI0)*(MPIP+MPI0)/4.);


Double_t cut_rho = 0.4;

Double_t FFC;
FFC = sqrt(p_pi_rho_rest_run*p_pi_rho_rest_run/(cut_rho*cut_rho+p_pi_rho_rest_run*p_pi_rho_rest_run))*p_pi_rho_rest/p_pi_rho_rest_run;

Double_t FFR;
FFR = sqrt(p_pi_rho_rest*p_pi_rho_rest/(cut_rho*cut_rho+p_pi_rho_rest*p_pi_rho_rest));

FF = FFC/FFR;

for (Int_t i_hro=0; i_hro<3; i_hro++){


amp[i_hro] = amp[i_hro]*G*FF;

};



};
