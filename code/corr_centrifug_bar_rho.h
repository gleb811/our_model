#include <TLorentzVector.h>


Float_t DSIMPS(Float_t *Y, Float_t *X, Int_t Npoints);

Bool_t rho_prot_orb_mom(Float_t res_j, Float_t res_parity,Int_t &max_val);

Float_t avrg_rho_mom (Int_t res_num);

Float_t avrg_rho_mom_no_res (Float_t W);

void corr_centrifug_bar_rho(Float_t W, TLorentzVector P4_PIP, TLorentzVector P4_PIM);



  

