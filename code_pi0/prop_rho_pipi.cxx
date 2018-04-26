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
#include "amp_rho_pipi.h"
#include "corr_centrifug_bar_rho.h"

void prop_rho_pipi(TComplex &prop, TLorentzVector P4_PIM, TLorentzVector P4_PIP) {
TLorentzVector P4_PIM_rho_rest,P4_PIP_rho_rest;
TComplex amp_dec_rho_pipi_rho_rest[3];
//Double_t FF;
Float_t p_rho = (P4_PIM+P4_PIP).Vect().Mag();
Float_t e_rho = (P4_PIM+P4_PIP)[3];
Float_t M_rho_run = (P4_PIM+P4_PIP).Mag();
Float_t p_pi = sqrt((M_rho_run/2.)*(M_rho_run/2.)-(MPIP+MPI0)*(MPIP+MPI0)/4.);

//cout <<"p_rho = "<< p_rho <<" e_rho = "<<e_rho<<"p_pi = "<<p_pi<<endl;

P4_PIP_rho_rest.SetXYZT(0.,0.,p_pi,M_rho_run/2.);
P4_PIM_rho_rest.SetXYZT(0.,0.,-1.*p_pi,M_rho_run/2.);

amp_rho_pipi(amp_dec_rho_pipi_rho_rest,P4_PIM_rho_rest,P4_PIP_rho_rest);

//cout <<" rrrr " << amp_dec_rho_pipi_rho_rest[1] << endl;

Double_t csum = amp_dec_rho_pipi_rho_rest[1]*(amp_dec_rho_pipi_rho_rest[1].Conjugate(TComplex(amp_dec_rho_pipi_rho_rest[1])));

Double_t Gamma_rho_dec_rest;

Gamma_rho_dec_rest = p_pi/(24.*MRHO*MRHO*M_PI)*csum;

Float_t beta = p_rho/e_rho;
Float_t gamma_lorentz = 1./(sqrt(1.-beta*beta));

//Float_t G = 5.724*1.056;
Float_t Gamma_rho_dec_cms;

//Gamma_rho_dec_cms = Gamma_rho_dec_rest*G*G*FF*FF/gamma_lorentz;
Gamma_rho_dec_cms = Gamma_rho_dec_rest/gamma_lorentz;

prop = TComplex(1,0);
prop = prop/TComplex(M_rho_run*M_rho_run-MRHO*MRHO,MRHO*Gamma_rho_dec_cms);
//cout <<  Gamma_rho_dec_rest*G*G << "  " << M_rho_run <<" iiii\n";
//cout << "eee " << prop << endl;

};
