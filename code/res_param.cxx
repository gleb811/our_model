#include <stdio.h>
#include <stdlib.h>  
#include "TROOT.h"
#include <RQ_OBJECT.h>
#include <iostream>
#include <string> 
#include <dlfcn.h>
#include <sstream>
#include <fstream>
#include <algorithm>



using namespace std;

const Int_t Number_of_res_const = 13;

Int_t Number_of_res = Number_of_res_const;

Int_t max_proj_j;


Float_t res_mass[Number_of_res_const],res_j[Number_of_res_const],res_i[Number_of_res_const],res_l[Number_of_res_const],res_parity[Number_of_res_const],res_width[Number_of_res_const];
Float_t A12[Number_of_res_const],A32[Number_of_res_const],S12[Number_of_res_const],BR_RHO_LS[Number_of_res_const][3][2];

Float_t factor_centrifugal[Number_of_res_const][3][2];

/*
BR_RHO_LS[Number_of_res_const][3][2] is the array of branching fractions for decay to the rho proton channel for different resonances.

The index [3] in the second brackets corresponds to the rho orbital momentum (L). For different resonances L can be 0,2... (for that case index 0 corresponds to L=0 (S), 1 to L=2 (D) ...) or 1,3... (for that case index 0 corresponds to L=1 (P), 1 to L=3 (F) ...).

The index [2] in the third brackets corresponds to the total spin of the rho and final state nucleon (index 0 corresponds to the total spin 1/2, while index 1 to the total spin 3/2).
*/

string res_name[Number_of_res_const];
Float_t rho_width_hel_sqrt[Number_of_res_const][3][2];

int res_param(){


res_name[0] = "P11_1440"; // N(1440)1/2+
res_mass[0] = 1.44;
res_j[0] = 0.5;
res_i[0] = 0.5;
res_l[0] = 1.0;
res_parity[0] = +1;
res_width[0] = 0.350;
A12[0]=21.8;
A32[0]=0.;
S12[0]=0.;
BR_RHO_LS[0][0][0] = 0.;
BR_RHO_LS[0][1][0] = 0.;
BR_RHO_LS[0][2][0] = 0.;
BR_RHO_LS[0][0][1] = 0.;
BR_RHO_LS[0][1][1] = 0.;
BR_RHO_LS[0][2][1] = 0.;


res_name[1] = "D13_1520"; // N(1520)3/2-
res_mass[1] = 1.520;
res_j[1] = 1.5;
res_i[1] = 0.5;
res_l[1] = 2.0;
res_parity[1] = -1;
res_width[1] = 0.120;
A12[1]=-62.2;
A32[1]=48.1;
S12[1]=-33.1*1.375;
BR_RHO_LS[1][0][0] = 0.;
BR_RHO_LS[1][1][0] = 0.;
BR_RHO_LS[1][2][0] = 0.;
BR_RHO_LS[1][0][1] = 0.21; // (rho_{3}N)_{S}
BR_RHO_LS[1][1][1] = 0.;
BR_RHO_LS[1][2][1] = 0.;


res_name[2] = "S11_1535"; // N(1535)1/2-
res_mass[2] = 1.535;
res_j[2] = 0.5;
res_i[2] = 0.5;
res_l[2] = 0.0;
res_parity[2] = -1;
res_width[2] = 0.150;
A12[2]=91.;
A32[2]=0.;
S12[2]=-13.6*1.406;
BR_RHO_LS[2][0][0] = 0.02; // (rho_{1}N)_{S}
BR_RHO_LS[2][1][0] = 0.;
BR_RHO_LS[2][2][0] = 0.;
BR_RHO_LS[2][0][1] = 0.;
BR_RHO_LS[2][1][1] = 0.01; // (rho_{3}N)_{D}
BR_RHO_LS[2][2][1] = 0.;



res_name[3] = "S11_1650"; // N(1650)1/2-
res_mass[3] = 1.650;
res_j[3] = 0.5;
res_i[3] = 0.5;
res_l[3] = 0.0;
res_parity[3] = -1;
res_width[3] = 0.150;
A12[3]=65.;
A32[3]=0.;
S12[3]=-12.8*1.65;
BR_RHO_LS[3][0][0] = 0.0; // (rho_{1}N)_{S}
BR_RHO_LS[3][1][0] = 0.;
BR_RHO_LS[3][2][0] = 0.;
BR_RHO_LS[3][0][1] = 0.;
BR_RHO_LS[3][1][1] = 0.03; // (rho_{3}N)_{D}
BR_RHO_LS[3][2][1] = 0.;


res_name[4] = "D15_1675"; // N(1675)5/2-
res_mass[4] = 1.675;
res_j[4] = 2.5;
res_i[4] = 0.5;
res_l[4] = 2.0;
res_parity[4] = -1;
res_width[4] = 0.150;
A12[4]=-0.143;
A32[4]=3.588;
S12[4]=0.;
BR_RHO_LS[4][0][0] = 0.;
BR_RHO_LS[4][1][0] = 0.004; // (rho_{1}N)_{D}
BR_RHO_LS[4][2][0] = 0.;
BR_RHO_LS[4][0][1] = 0.;
BR_RHO_LS[4][1][1] = 0.002; // (rho_{3}N)_{D}
BR_RHO_LS[4][2][1] = 0.;



res_name[5] = "F15_1680"; // N(1680)5/2+
res_mass[5] = 1.680;
res_j[5] = 2.5;
res_i[5] = 0.5;
res_l[5] = 3.0;
res_parity[5] = +1;
res_width[5] = 0.130;
A12[5]=-27.;
A32[5]=55.;
S12[5]=-17.*1.72;
BR_RHO_LS[5][0][0] = 0.;
BR_RHO_LS[5][1][0] = 0.; 
BR_RHO_LS[5][2][0] = 0.;
BR_RHO_LS[5][0][1] = 0.05; // (rho_{3}N)_{P}
BR_RHO_LS[5][1][1] = 0.02; // (rho_{3}N)_{F}
BR_RHO_LS[5][2][1] = 0.;



res_name[6] = "D13_1700"; // N(1700)3/2-
res_mass[6] = 1.700;
res_j[6] = 1.5;
res_i[6] = 0.5;
res_l[6] = 2.0;
res_parity[6] = -1;
res_width[6] = 0.100;
A12[6]=-2.;
A32[6]=2.;
S12[6]=-3.*1.76;
BR_RHO_LS[6][0][0] = 0.;
BR_RHO_LS[6][1][0] = 0.; 
BR_RHO_LS[6][2][0] = 0.;
BR_RHO_LS[6][0][1] = 0.13; // (rho_{3}N)_{S}
BR_RHO_LS[6][1][1] = 0.; 
BR_RHO_LS[6][2][1] = 0.;



res_name[7] = "P13_1720"; // N(1720)3/2+
res_mass[7] = 1.720;
res_j[7] = 1.5;
res_i[7] = 0.5;
res_l[7] = 1.0;
res_parity[7] = +1;
res_width[7] = 0.150;
A12[7]=43.;
A32[7]=-39.;
S12[7]=-37.*1.80;
BR_RHO_LS[7][0][0] = 0.87; // (rho_{1}N)_{P}
BR_RHO_LS[7][1][0] = 0.; 
BR_RHO_LS[7][2][0] = 0.;
BR_RHO_LS[7][0][1] = 0.;
BR_RHO_LS[7][1][1] = 0.; 
BR_RHO_LS[7][2][1] = 0.;



res_name[8] = "S31_1620"; // Delta(1620)1/2-
res_mass[8] = 1.620;
res_j[8] = 0.5;
res_i[8] = 1.5;
res_l[8] = 0.0;
res_parity[8] = -1;
res_width[8] = 0.150;
A12[8]=9.;
A32[8]=0.;
S12[8]=-44.*1.58;
BR_RHO_LS[8][0][0] = 0.25; // (rho_{1}N)_{S}
BR_RHO_LS[8][1][0] = 0.; 
BR_RHO_LS[8][2][0] = 0.;
BR_RHO_LS[8][0][1] = 0.;
BR_RHO_LS[8][1][1] = 0.04; // (rho_{3}N)_{D} 
BR_RHO_LS[8][2][1] = 0.;


res_name[9] = "D33_1700"; // Delta(1700)3/2-
res_mass[9] = 1.700;
res_j[9] = 1.5;
res_i[9] = 1.5;
res_l[9] = 2.0;
res_parity[9] = -1;
res_width[9] = 0.300;
A12[9]=35.;
A32[9]=38.;
S12[9]=17.*1.76;
BR_RHO_LS[9][0][0] = 0.;
BR_RHO_LS[9][1][0] = 0.; 
BR_RHO_LS[9][2][0] = 0.;
BR_RHO_LS[9][0][1] = 0.08;// (rho_{3}N)_{S} 
BR_RHO_LS[9][1][1] = 0.; 
BR_RHO_LS[9][2][1] = 0.;



res_name[10] = "F35_1905"; // Delta(1905)5/2+
res_mass[10] = 1.905;
res_j[10] = 2.5;
res_i[10] = 1.5;
res_l[10] = 3.0;
res_parity[10] = +1;
res_width[10] = 0.350;
A12[10]=-25.859;
A32[10]=-53.879;
S12[10]=0;
BR_RHO_LS[10][0][0] = 0.;
BR_RHO_LS[10][1][0] = 0.; 
BR_RHO_LS[10][2][0] = 0.;
BR_RHO_LS[10][0][1] = 0.86;// (rho_{3}N)_{P} 
BR_RHO_LS[10][1][1] = 0.; 
BR_RHO_LS[10][2][1] = 0.;


res_name[11] = "F37_1950"; // Delta(1950)7/2+
res_mass[11] = 1.950;
res_j[11] = 3.5;
res_i[11] = 1.5;
res_l[11] = 3.0;
res_parity[11] = +1;
res_width[11] = 0.300;
A12[11]=-23.830;
A32[11]=-68.789;
S12[11]=0;
BR_RHO_LS[11][0][0] = 0.;
BR_RHO_LS[11][1][0] = 0.; 
BR_RHO_LS[11][2][0] = 0.;
BR_RHO_LS[11][0][1] = 0.;
BR_RHO_LS[11][1][1] = 0.43;// (rho_{3}N)_{F}  
BR_RHO_LS[11][2][1] = 0.;



res_name[12] = "P33_1600"; // Delta(1600)3/2+
res_mass[12] = 1.600;
res_j[12] = 1.5;
res_i[12] = 1.5;
res_l[12] = 1.0;
res_parity[12] = +1;
res_width[12] = 0.350;
A12[12]=0.;
A32[12]=0.;
S12[12]=0.;
BR_RHO_LS[12][0][0] = 0.;
BR_RHO_LS[12][1][0] = 0.; 
BR_RHO_LS[12][2][0] = 0.;
BR_RHO_LS[12][0][1] = 0.;
BR_RHO_LS[12][1][1] = 0.;
BR_RHO_LS[12][2][1] = 0.;

max_proj_j = Int_t(2. * *std::max_element(res_j,res_j+Number_of_res) + 1.);

};
