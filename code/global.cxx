#include <stdio.h>
#include <stdlib.h>  
#include "TROOT.h"
#include <RQ_OBJECT.h>
#include <iostream>
#include "TDatabasePDG.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include <string> 
#include <dlfcn.h>
#include <sstream>
#include <fstream>
#include "TTree.h"



using namespace std;

Float_t MP,MRHO,MPIP,MPIM,Me;

Float_t WIDTH_RHO;

Float_t alph_const; 

int global(){
//Starting valuesconst char *HOME_ROOT;

const char *HOME_ROOT;

    
//HOME_ROOT = getenv("ROOT");
system("root_home=`root-config --etcdir`");
HOME_ROOT = getenv("root_home");
ostringstream ROOT_DIR;
ROOT_DIR << HOME_ROOT << "/pdg_table.txt";

TDatabasePDG *pdg = new TDatabasePDG();
pdg->ReadPDGTable(ROOT_DIR.str().c_str());
TParticlePDG *part1 = new TParticlePDG();

part1 = pdg->GetParticle("proton");
MP= part1->Mass();  

part1 = pdg->GetParticle("rho0");
MRHO = part1->Mass();
WIDTH_RHO = part1->Width();


part1 = pdg->GetParticle("pi+");
MPIP= part1->Mass(); 
 

part1 = pdg->GetParticle("pi-");
MPIM= part1->Mass();

alph_const = 1./137.035; 

part1 = pdg->GetParticle("e-");
Me= part1->Mass(); 
};
