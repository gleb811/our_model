void mok_sigma_t_l_etc() {

ostringstream qqq;

for (Int_t ii = 0; ii < 17; ii++) {
qqq.str("");
qqq << "sigma_diff_" << 14125 + ii*250 << ".dat";
//cout << qqq.str() << endl;
Double_t S12_val, S23_val, TH_val, ALP_val;
Double_t s12_min, s12_max,s23_min, s23_max, th_min, th_max, alp_min, alp_max;
Double_t ds12_tmp,ds23_tmp, dalpha_tmp;
Double_t dalpha, ds12, ds23;
Double_t th_l, th_r, dtheta, dtheta_tmp;

Double_t mpr = 0.938272;
Double_t mpi = 0.13957;

//Double_t true_min_s12, true_max_s12, true_min_s23, true_max_s23;


string dummy,xsect;
Double_t cross_sect, cross_sect_t, cross_sect_l,cross_sect_c2f,cross_sect_s2f,cross_sect_cf,cross_sect_sf, epsl,W,Xsect_int, sigma_t_int, sigma_l_int, sigma_c2f_int, sigma_s2f_int, sigma_cf_int, sigma_sf_int;

Xsect_int = 0.;
sigma_t_int = 0.;
sigma_l_int  = 0.;
sigma_c2f_int = 0.;
sigma_s2f_int = 0.;
sigma_cf_int = 0.;
sigma_sf_int = 0.;





//Read from file minimum and maximum values of all variables
ifstream input(qqq.str().c_str());
if(input.is_open()){

for (Int_t is23 = 1; is23 <=12; is23++) {
for (Int_t is12 = 1; is12 <=12; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {
getline(input,dummy);
S12_val = atof(dummy.c_str());
getline(input,dummy);
S23_val = atof(dummy.c_str());
getline(input,dummy);
TH_val = atof(dummy.c_str());
getline(input,dummy);
ALP_val = atof(dummy.c_str());

getline(input,xsect);
cross_sect_t = atof(xsect.c_str());
getline(input,xsect);
cross_sect_l = atof(xsect.c_str());

getline(input,xsect);
cross_sect_c2f = atof(xsect.c_str());
getline(input,xsect);
cross_sect_s2f = atof(xsect.c_str());
getline(input,xsect);
cross_sect_cf = atof(xsect.c_str());
getline(input,xsect);
cross_sect_sf = atof(xsect.c_str());


getline(input,dummy);
epsl = atof(dummy.c_str());
getline(input,dummy);
W = atof(dummy.c_str());



if ((is23==1)&&(itheta==1)&&(ialpha==1)&&(is12==1)) s12_min = S12_val;
if ((is23==1)&&(itheta==1)&&(ialpha==1)&&(is12==12)) s12_max = S12_val;

if ((is12==1)&&(itheta==1)&&(ialpha==1)&&(is23==1)) s23_min = S23_val;
if ((is12==1)&&(itheta==1)&&(ialpha==1)&&(is23==12)) s23_max = S23_val;

if ((is12==1)&&(itheta==1)&&(ialpha==1)&&(is23==1)) th_min = TH_val;
if ((is12==1)&&(itheta==6)&&(ialpha==1)&&(is23==1)) th_max = TH_val;

if ((is12==1)&&(itheta==1)&&(ialpha==1)&&(is23==1)) alp_min = ALP_val;
if ((is12==1)&&(itheta==1)&&(ialpha==6)&&(is23==1)) alp_max = ALP_val;

};
};
};
};
};
input.close();


//Determine the with of the bin over all variables 
ds12 = (s12_max - s12_min)/11;
ds23 = (s23_max - s23_min)/11;
dalpha = (alp_max - alp_min)/5;
dtheta = (th_max - th_min)/5;

ds12_tmp = ds12; 
ds23_tmp = ds23; 
dalpha_tmp = dalpha; 
dtheta_tmp = dtheta;
//cout << "\n";
//cout << "s12_min = "<< s12_min<<", s12_max = "<< s12_max <<", ds12 = " << ds12<< "\n";
//cout << "s23_min = "<< s23_min<<", s23_max = "<< s23_max <<", ds23 = " << ds23<< "\n";
//cout << "th_min = "<< th_min<<", th_max = "<< th_max <<", dtheta = " << dtheta<< "\n";
//cout << "alp_min = "<< alp_min<<", alp_max = "<< alp_max <<", dalpha = " << dalpha<< "\n";

//cout << "W = "<< W << " GeV \n";


//true_min_s12 = (mpi+mpi)*(mpi+mpi);
//true_max_s12 = (W - mpr)*(W - mpr);

//true_min_s23 = (mpi+mpr)*(mpi+mpr);
//true_max_s23 = (W - mpi)*(W - mpi);


//cout << true_min_s12-s12_min  << " WWWWWWW " << true_max_s12-s12_max << "\n";

ifstream inp1(qqq.str().c_str());
if(inp1.is_open()){


for (Int_t is23 = 1; is23 <=12; is23++) {
for (Int_t is12 = 1; is12 <=12; is12++) {

for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {
//I am doing this to force this variables renew each time loops run, beceuse they are determied outside the loop and sometimes change inside the loop
ds12=ds12_tmp;
ds23=ds23_tmp;
dalpha = dalpha_tmp;
dtheta = dtheta_tmp;

getline(inp1,dummy);
getline(inp1,dummy);
getline(inp1,dummy);
getline(inp1,dummy);

getline(inp1,xsect);
cross_sect_t = atof(xsect.c_str());
getline(inp1,xsect);
cross_sect_l = atof(xsect.c_str());

getline(inp1,xsect);
cross_sect_c2f = atof(xsect.c_str());
getline(inp1,xsect);
cross_sect_s2f = atof(xsect.c_str());
getline(inp1,xsect);
cross_sect_cf = atof(xsect.c_str());
getline(inp1,xsect);
cross_sect_sf = atof(xsect.c_str());

getline(inp1,dummy);
epsl = atof(dummy.c_str());
getline(inp1,dummy);
W = atof(dummy.c_str());
//cout << W << "\n";getline(inp1,xsect);


//if the point is the first or the last then the width of the bin is smaller 
if ((is12==1)||(is12==12)) ds12=ds12_tmp/2;
//if ((is12==16)) cout << ds12 << "\n";

if ((is23==1)||(is23==12)) ds23=ds23_tmp/2;
//if ((is23==14)) cout << ds23 << "\n";

if ((ialpha==1)||(ialpha==6)) dalpha=dalpha_tmp/2+0.01;
//if ((ialpha==2)) cout << dalpha << "\n";


//left and right edge for theta if the point is not first or last
th_l = th_min + dtheta/2.+(itheta - 2.)*dtheta;
th_r = th_min + dtheta/2.+(itheta -1.)*dtheta;

//left and right edge for theta if the point is first or last
if (itheta==1) th_l = th_min-0.01;
if (itheta==1) th_r = th_min+dtheta/2.;

if (itheta==6) th_l = th_max - dtheta/2.;
if (itheta==6) th_r = th_max+0.01;




//if ((is12==2)&&(itheta==2)&&(ialpha==2)) cout << is23<< "    "<< "ds23 = " <<ds23 << "\n";

//if ((is23==2)&&(itheta==2)&&(ialpha==2)) cout << is12<< "    "<< "ds12 = " <<ds12 << "\n";

//if ((is23==2)&&(is12==2)&&(itheta==2)) cout << ialpha<< "    "<< "dalpha = " <<dalpha << "\n";

//if ((is12==2)&&(is23==2)&&(ialpha==2)) cout << itheta<< "    "<<"th_l = "<< th_l <<", th_r = "<<th_r<<"  " << "   "<< ", th_r - th_l = "<<th_r-th_l<<  "\n";



cross_sect = cross_sect_t+epsl*cross_sect_l; 
//+ epsl*cross_sect_l;
//cross_sect = epsl*cross_sect_l;
Xsect_int = Xsect_int + cross_sect*(cos(th_l)-cos(th_r))*ds12*ds23*dalpha; 

sigma_t_int = sigma_t_int + cross_sect_t*(cos(th_l)-cos(th_r))*ds12*ds23*dalpha; 
sigma_l_int = sigma_l_int +  cross_sect_l*(cos(th_l)-cos(th_r))*ds12*ds23*dalpha; 
sigma_c2f_int = sigma_c2f_int + cross_sect_c2f*(cos(th_l)-cos(th_r))*ds12*ds23*dalpha; 
sigma_s2f_int = sigma_s2f_int + cross_sect_s2f*(cos(th_l)-cos(th_r))*ds12*ds23*dalpha; 
sigma_cf_int  = sigma_cf_int + cross_sect_cf*(cos(th_l)-cos(th_r))*ds12*ds23*dalpha; 
sigma_sf_int = sigma_sf_int +  cross_sect_sf*(cos(th_l)-cos(th_r))*ds12*ds23*dalpha; 
//cout<<cross_sect_sf <<"\n";

};
};
};
};
};



Xsect_int = Xsect_int*(2*3.1415);

sigma_t_int = sigma_t_int*(2*3.1415);
sigma_l_int = sigma_l_int*(2*3.1415);
sigma_c2f_int = sigma_c2f_int*(2*3.1415);
sigma_s2f_int = sigma_s2f_int*(2*3.1415);
sigma_cf_int  = sigma_cf_int*(2*3.1415); 
sigma_sf_int = sigma_sf_int*(2*3.1415);

//std::ofstream ofs ("qwqwqw.txt", std::ofstream::out);

//ofs << W << "   " << Xsect_int << "\n";


//cout << W <<"  " << Xsect_int << "\n";

if (ii < 16 ) {
cout << " W = " << W << "  " << Xsect_int <<  endl;
//cout <<  Xsect_int<< ", ";
} else {
cout << " W = " << W << "  " << Xsect_int << endl;
//cout <<  Xsect_int<< endl;
};

//cout << sigma_t_int+ epsl*sigma_l_int << "\n";

inp1.close();

};


};
