#include "TROOT.h"
#include "THnSparse.h"
#include "TFile.h"
#include "TH1.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include "global.h"
#include "res_param.h"
#include "corr_centrifug_bar_rho.h"
#include "TParticle.h"
#include "anti_rot.h"
#include "rho_ls_to_hel.h"
#include "TComplex.h"
#include "Clebsch.h"
#include "amp_elmag.h"
#include "amp_rho_pipi.h"
#include "prop_rho_pipi.h"

 using namespace std; 
 
 
 TComplex BWP(Float_t W, Int_t res_num, Float_t width_rho) {

TComplex BWP;


BWP = TComplex(1.,0)/(TComplex(Double_t(res_mass[res_num]*res_mass[res_num]-W*W),-1.*Double_t(res_mass[res_num]*width_rho)));

//cout << "width_rho = " << width_rho << " BWP = " << BWP << " W = " << W << endl;

//BWP = TComplex(res_mass[res_num]*res_mass[res_num]-W*W,-1.*res_mass[res_num]*res_width[res_num]);
//BWP = 1./BWP;

return BWP;

};

 
Float_t G_BYCKLING(Float_t x, Float_t y, Float_t z, Float_t u, Float_t v, Float_t w) {
    return x*x*y+x*y*y+z*z*u+z*u*u+v*v*w+v*w*w+x*z*w+x*u*v+y*z*v+y*u*w-x*y*(z+u+v+w)-z*u*(x+y+v+w)-v*w*(x+y+z+u);
    }; 
 
  
  
int main(int argc, char** argv) {

global();
res_param();


THnSparseD *topology_1;
TH1D *h_alpha;
TH1D *h_theta;
TH1D *h_dcostheta;
TH1D *h_spippim;

Int_t bins[4];
Double_t xmin[4];
Double_t xmax[4];

Float_t sigma_int;

Float_t sigma_tot;
Float_t sigma_t,sigma_l;
Float_t W = 1.4375;
Float_t Q2 = 0.65;
//Float_t phi_el = 0.; 
Float_t E_beam = 2.445;
Float_t M12;
Float_t M23;
Float_t theta_hadr;
Float_t alpha_hadr;
Float_t phi_hadr  = 0.;
Float_t m1 = MPIM;
Float_t m2 = MPIP;
Float_t m3 = MP;
TLorentzVector P4_1,P4_2,P4_3;
TLorentzVector P4_RHO;

Float_t nu;
Float_t theta_e_prime;
Float_t eps_t,eps_l;

Float_t S[2]={0.5,1.5};
Float_t HRO[3]={-1.,0.,1.};
Float_t HNUF[2]={0.5,-0.5};
Float_t HNUI[2]={0.5,-0.5};
Float_t HGAM[3]={-1.,0.,1};
Float_t HRE[4]={-1.5,-0.5,0.5,1.5};

Float_t E_RHO_res,P_RHO_mag_res;

TComplex amp_dec_to_rho_p[Number_of_res][2][3][4];
TComplex amp_gam_nu_to_res[Number_of_res][2][3];
TComplex amp_gam_nu_res_rho_p[Number_of_res][2][3][2][3];
TComplex amp_dec_rho_pipi[3];
TComplex amp_n_rho_tot_res[2][3][2];
TComplex tmp;

Double_t Isospin_Clebsch[Number_of_res];

for (Int_t res_num = 0; res_num < Number_of_res; res_num++) {
//(I_proton,I_rho,I_res,Iz_proton,Iz_rho,Iz_res=Iz_proton+Iz_rho,Isospin_Clebsch)

Clebsch(0.5,1.,res_i[res_num],0.5,0.,0.5,Isospin_Clebsch[res_num]);
//Clebsch(0.5,1.,res_i[res_num],-0.5,1.,0.5,Isospin_Clebsch[res_num]); for gamma p -> n pi+ pi0 channel

};

Float_t tot_width_rho[Number_of_res];

for (Int_t res_num = 0; res_num < Number_of_res; res_num++) {
tot_width_rho[res_num] = 0.;
for (Int_t i_l_rho_nu = 0; i_l_rho_nu < 3; i_l_rho_nu++) {
for (Int_t i_s_rho_nu = 0; i_s_rho_nu < 2; i_s_rho_nu++) {
tot_width_rho[res_num] =  tot_width_rho[res_num] + BR_RHO_LS[res_num][i_l_rho_nu][i_s_rho_nu]*res_width[res_num];
};
};
};

Int_t Ns12=12;
Float_t s12min = (m1+m2)*(m1+m2);


//h_spippim = new TH1D("h_spippim","h_spippim",Ns12,(s12min - s12step/2.),(s12max + s12step/2.));

Int_t Ns23=12;
Float_t s23min = (m2+m3)*(m2+m3);


Int_t Ntheta=6;
Float_t thetamin = 0.01;
Float_t thetamax = M_PI-0.01;
Float_t thetastep=(thetamax-thetamin)/(Ntheta-1);

h_dcostheta = new TH1D("h_dcostheta","h_dcostheta",Ntheta,(thetamin - thetastep/2.)*180./M_PI,(thetamax + thetastep/2.)*180./M_PI);

h_theta = new TH1D("h_theta","h_theta",Ntheta,(thetamin - thetastep/2.)*180./M_PI,(thetamax + thetastep/2.)*180./M_PI);

for (Int_t itheta=1; itheta<=Ntheta; itheta++) {
//h_dcostheta->Fill((thetamin+thetastep*(itheta-1))*180./M_PI,cos(thetamin+thetastep*(itheta-1.5))-cos(thetamin+thetastep*(itheta-0.5)));
h_dcostheta->Fill((thetamin+thetastep*(itheta-1))*180./M_PI,sin(thetamin+thetastep*(itheta-1))*thetastep);
};

Int_t Nphi=6;
Float_t phimin = 0. + 0.01;
Float_t phimax = 2.*M_PI - 0.01;
Float_t phistep=(phimax-phimin)/(Nphi-1);
Float_t alpha[Nphi];
for (Int_t iphi=1; iphi<=Nphi; iphi++) {alpha[iphi-1] = 0.;};

bins[0] = Ns12;
bins[1] = Ns23;
bins[2] = Ntheta;
bins[3] = Nphi;


xmin[2] = (thetamin - thetastep/2.)*180./M_PI;
xmin[3] = (phimin - phistep/2.)*180./M_PI;


xmax[2] = (thetamax + thetastep/2.)*180./M_PI;
xmax[3] = (phimax + phistep/2.)*180./M_PI;

topology_1 = new THnSparseD("topology_1","topology_1",4,bins,xmin,xmax);


// Start loop over W

ostringstream out_file_name;

for (Int_t iw=0; iw<=16; iw++) {
W = 1.4125 + 0.025*iw;

Float_t s12max = (W-m3)*(W-m3);
Float_t s12step=(s12max-s12min)/(Ns12-1);

Float_t s23max = (W-m1)*(W-m1);
Float_t s23step=(s23max-s23min)/(Ns23-1);

xmin[0] = s12min - s12step/2.;
xmin[1] = s23min - s23step/2.;

xmax[0] = s12max + s12step/2.;
xmax[1] = s23max + s23step/2.;

out_file_name.str("");
out_file_name << "sigma_diff_" << W*10000 << ".dat";
//cout << out_file_name.str() << endl;

std::ofstream ofs (out_file_name.str().c_str(), std::ofstream::out);

sigma_int = 0.;

for (Int_t is23=1; is23<=Ns23; is23++) {
for (Int_t is12=1; is12<=Ns12; is12++) {
for (Int_t itheta=1; itheta<=Ntheta; itheta++) {
for (Int_t iphi=1; iphi<=Nphi; iphi++) {
//for (Int_t is23=9; is23<=9; is23++) {
//for (Int_t is12=10; is12<=10; is12++) {
//for (Int_t itheta=6; itheta<=6; itheta++) {
//for (Int_t iphi=2; iphi<=2; iphi++) {


M12 = sqrt(s12min+s12step*(is12-1));
M23 = sqrt(s23min+s23step*(is23-1));
theta_hadr = thetamin+thetastep*(itheta-1);
alpha_hadr = phimin+phistep*(iphi-1);



Double_t Var_1[4];
Var_1[0] = M12*M12;
Var_1[1] = M23*M23;
Var_1[2] = theta_hadr*180./M_PI;
Var_1[3] = alpha_hadr*180./M_PI;

//cout << Var_1[0] << "  " << Var_1[1] << "  " << Var_1[2]*M_PI/180. << "  " << Var_1[3]*M_PI/180.<< "  " << W << endl;

//theta_hadr = 9.9999997764825821E-003;
//alpha_hadr = 9.9999997764825821E-003;

//cout << "Bycling " << G_BYCKLING(M12*M12, M23*M23, W*W, m2*m2, m1*m1, m3*m3) << endl;

if (G_BYCKLING(M12*M12, M23*M23, W*W, m2*m2, m1*m1, m3*m3) < 0.) {

anti_rot(W,Q2,M12,M23,theta_hadr,alpha_hadr,phi_hadr,0.13956789672374725,0.13956789672374725,0.93827229738235474,P4_1,P4_2,P4_3);



amp_elmag((TComplex *) amp_gam_nu_to_res,W,Q2,P4_1,P4_2);
corr_centrifug_bar_rho(W,P4_1,P4_2);


/*
for (Int_t i_hnui = 0; i_hnui < 2; i_hnui++) {
for (Int_t i_hgam = 0; i_hgam < 3; i_hgam++) {
cout << " Hnui = " <<  HNUI[i_hnui] << " Hgam = " << HGAM[i_hgam] << amp_gam_nu_to_res[4][i_hnui][i_hgam] << endl;
};
};
*/

rho_ls_to_hel((TComplex *) amp_dec_to_rho_p,P4_3,M12);

for (Int_t res_num = 0; res_num < Number_of_res; res_num++) {
//for (Int_t res_num = 0; res_num <= 12; res_num++) {




for (Int_t i_hnuf = 0; i_hnuf < 2; i_hnuf++) {
for (Int_t i_hro = 0; i_hro < 3; i_hro++) {
for (Int_t i_hre = 0; i_hre < 4; i_hre++) {
amp_dec_to_rho_p[res_num][i_hnuf][i_hro][i_hre] = amp_dec_to_rho_p[res_num][i_hnuf][i_hro][i_hre]*Isospin_Clebsch[res_num];


if (res_num == 5) {

//cout <<" Hnuf = " <<  HNUF[i_hnuf] << " Hro = " << HRO[i_hro] << " Hre = " << HRE[i_hre] << amp_dec_to_rho_p[res_num][i_hnuf][i_hro][i_hre]/Isospin_Clebsch[res_num] << endl;
//cout << "BWP = " << BWP(W,res_num,tot_width_rho[res_num])  << " W = " << W << endl;
};

};
};
};


for (Int_t i_hnui = 0; i_hnui < 2; i_hnui++) {
for (Int_t i_hgam = 0; i_hgam < 3; i_hgam++) {
for (Int_t i_hnuf = 0; i_hnuf < 2; i_hnuf++) {
for (Int_t i_hro = 0; i_hro < 3; i_hro++) {

Int_t i_hre = 1.5 + HGAM[i_hgam] - HNUI[i_hnui];

amp_gam_nu_res_rho_p[res_num][i_hnui][i_hgam][i_hnuf][i_hro] = amp_dec_to_rho_p[res_num][i_hnuf][i_hro][i_hre]*amp_gam_nu_to_res[res_num][i_hnui][i_hgam]*BWP(W,res_num,res_width[res_num]);


};   // hro
};  // hnuf
};  //hgam
}; // hnui



};// num_res





amp_rho_pipi(amp_dec_rho_pipi,P4_1,P4_2);

TComplex prop;
prop_rho_pipi(prop,P4_1,P4_2);



for (Int_t i_hgam = 0; i_hgam < 3; i_hgam++) {
for (Int_t i_hnui = 0; i_hnui < 2; i_hnui++) {
for (Int_t i_hnuf = 0; i_hnuf < 2; i_hnuf++) {
amp_n_rho_tot_res[i_hnuf][i_hgam][i_hnui] = TComplex(0.,0.);
//tmp = TComplex(0.,0.);

for (Int_t i_hro = 0; i_hro < 3; i_hro++) {


//tmp = tmp + amp_dec_rho_pipi[i_hro]*prop*amp_gam_nu_res_rho_p[res_num][i_hnui][i_hgam][i_hnuf][i_hro];

for (Int_t res_num = 0; res_num < Number_of_res; res_num++) {
amp_n_rho_tot_res[i_hnuf][i_hgam][i_hnui] = amp_n_rho_tot_res[i_hnuf][i_hgam][i_hnui] + amp_dec_rho_pipi[i_hro]*prop*amp_gam_nu_res_rho_p[res_num][i_hnui][i_hgam][i_hnuf][i_hro];
}; // res_num



}; // i_hro

//amp_n_rho_tot_res[i_hnuf][i_hgam][i_hnui] = amp_n_rho_tot_res[i_hnuf][i_hgam][i_hnui] + tmp;

}; // i_hnuf
}; // i_hnui
}; // i_hgam






//cout << " prop = " << prop << endl;

TComplex Mx = TComplex(0.,0.);
TComplex My = TComplex(0.,0.);
TComplex Mz = TComplex(0.,0.);

Float_t Mxx =0.;
Float_t Myy =0.;
Float_t Mzz =0.;

for (Int_t i_hnui = 0; i_hnui < 2; i_hnui++) {
for (Int_t i_hnuf = 0; i_hnuf < 2; i_hnuf++) {

Mx = Double_t(-1./sqrt(2.))*(amp_n_rho_tot_res[i_hnuf][2][i_hnui]-amp_n_rho_tot_res[i_hnuf][0][i_hnui]);

My = Double_t(1./sqrt(2.))*TComplex(0.,1.)*(amp_n_rho_tot_res[i_hnuf][2][i_hnui]+amp_n_rho_tot_res[i_hnuf][0][i_hnui]);

Mz = amp_n_rho_tot_res[i_hnuf][1][i_hnui];

//cout << "Mx = " << Mx << " My = " << My << " Mz = " << Mz <<endl;

Mxx = Mxx + 0.5*(Mx*Mx.Conjugate(Mx)).Re();
Myy = Myy + 0.5*(My*My.Conjugate(My)).Re();
Mzz = Mzz + 0.5*(Mz*Mz.Conjugate(Mz)).Re();

}; // i_hnuf
}; // i_hnui

sigma_t = (Mxx+Myy);
sigma_l = 2.*Mzz;

//cout << "Mxx = " << Mxx << " Myy = " << Myy << " Mzz = " << Mzz <<endl;

// The photon energy in the laboratory system at Q^2 = 0 (i.e. the equivalent energy of a real photon).
Float_t K=(W*W-MP*MP)/(2.*MP); 
// 4.*K*MP is the invariant flux of photons that is equivalent to the invariant fluz of real photons 
Float_t F = 1./(4.*K*MP); 
// Phase space kinematic factor that corresponds to three-body final state
Float_t PHF = 1./(32.*W*W*pow((2.*M_PI),5.)); 

Float_t factor = F*alph_const*4.*M_PI*0.5*PHF*389.379; // the factor 389.379 is (hc)^2 in units GeV^2*mub

sigma_t = sigma_t*factor;
sigma_l = sigma_l*factor;

//cout << "S12 = " << Var_1[0] << " S23 = " << Var_1[1] << " theta = " << Var_1[2]*M_PI/180. << " alpha = " << Var_1[3]*M_PI/180.<< "\n";
//cout << "sigma_t = " << sigma_t << " sigma_l = " << sigma_l <<endl;

nu = (W*W+Q2-MP*MP)/(2.*MP); 
theta_e_prime = acos(1.-Q2/E_beam/(E_beam-nu)/2.);

eps_t = 1./(1.+2.*(1.+nu*nu/Q2)*tan(theta_e_prime/2.)*tan(theta_e_prime/2.));
eps_l = Q2/nu/nu*eps_t;

sigma_tot = sigma_t + eps_l*sigma_l;

 } else { // Bycling function >0
 
sigma_tot = 0.; 
sigma_t = 0.; 
sigma_l = 0.; 
 
 }; // end if Bycling function < 0
 
//cout << "sigma_tot = " << sigma_tot << " S12 = " << M12*M12 << " S23 = " << M23*M23 << " theta = " << theta_hadr*180./M_PI << " alpha = " << alpha_hadr*180./M_PI <<endl; 





// Output 5diff cr.sec. to the text file


ofs << Var_1[0] << "\n";
ofs << Var_1[1] << "\n";
ofs << Var_1[2]*M_PI/180. << "\n";
ofs << Var_1[3]*M_PI/180. << "\n";
ofs << sigma_t << "\n";
ofs << sigma_l << "\n";
ofs << 0. << "\n";
ofs << 0. << "\n";
ofs << 0. << "\n";
ofs << 0. << "\n";
ofs << eps_l << "\n";
ofs << W << "\n";



// Fill out 5D histogram

topology_1->Fill(Var_1,sigma_tot*sin(theta_hadr)*thetastep); 
alpha[iphi-1] = alpha[iphi-1] + sigma_tot*sin(theta_hadr)*thetastep*s12step*s23step*2.*M_PI;

if ((is12 == 1)||(is12 == Ns12))sigma_tot = sigma_tot/2.;
if ((is23 == 1)||(is23 == Ns23))sigma_tot = sigma_tot/2.;
if ((itheta == 1)||(itheta == Ntheta))sigma_tot = sigma_tot/2.;
if ((iphi == 1)||(iphi == Nphi))sigma_tot = sigma_tot/2.;
sigma_int = sigma_int +sigma_tot*sin(theta_hadr)*thetastep*s12step*s23step*2.*M_PI*phistep;

}; // iphi
}; // itheta
}; // is23
}; // is12



h_alpha=topology_1->Projection(3,"");
h_alpha->Scale(2.*M_PI);
h_alpha->Scale(s12step);
h_alpha->Scale(s23step);
h_alpha->SetName("h_alpha");
h_alpha->SetTitle("h_alpha");

h_theta=topology_1->Projection(2,"");
h_theta->Scale(2.*M_PI);
h_theta->Scale(s12step);
h_theta->Scale(s12step);
h_theta->Scale(phistep);
h_theta->Divide(h_dcostheta);
h_theta->SetName("h_theta");
h_theta->SetTitle("h_theta");

/*
h_spippim=topology_1->Projection(0,"");
h_spippim->Scale(2.*M_PI);
h_spippim->Scale((s23max-s23min));
h_spippim->SetName("h_spippim");
h_spippim->SetTitle("h_spippim");
*/

TFile *outFile = new TFile("output.root","recreate");

outFile->cd();
topology_1->Write();
h_alpha->Write();
h_theta->Write();
//h_spippim->Write();
h_dcostheta->Write();
outFile->Close();
//for (Int_t iphi=1; iphi<=Nphi; iphi++) {cout << alpha[iphi-1] << endl;};

cout << "W = " << W << " sigma_int = " << sigma_int << endl;


ofs.close();

}; // iw
};

