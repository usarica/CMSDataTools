#ifndef EXT_CCONST_H
#define EXT_CCONST_H

#include <cmath>
#include "TMath.h"

using namespace std;


float getDVBF2jetsConstant(float ZZMass){
  float par[9]={
    1.876,
    -55.488,
    403.32,
    0.3906,
    80.8,
    27.7,
    -0.06,
    54.97,
    309.96
  };
  float kappa =
    pow(1.-atan((ZZMass-par[1])/par[2])*2./TMath::Pi(), par[0])
    + par[3]*exp(-pow((ZZMass-par[4])/par[5], 2))
    + par[6]*exp(-pow((ZZMass-par[7])/par[8], 2));
  float constant = kappa/(1.-kappa);
  return constant;
}

float getDVBF1jetConstant(float ZZMass){
  float par[8]={
    0.395,
    -0.07,
    85.,
    30.,
    -0.691,
    -5659.47,
    5734.37,
    0.75
  };
  float kappa =
    par[0]
    + par[1]*exp(-pow((ZZMass-par[2])/par[3], 2))
    + par[4]*pow(log((ZZMass-par[5])/par[6]), par[7])*(ZZMass>=(par[5]+par[6]));
  float constant = kappa/(1.-kappa);
  return constant;
}

// [0]+[1]*exp(-pow((x-[2])/[3], 2))+[4]*exp(-pow((x-[5])/[6], 2))+[7]*(exp(-pow((x-[8])/[9], 2))*(x<[8])+exp(-pow((x-[8])/[10], 2))*(x>=[8]))+[11]*exp(-pow((x-[12])/[13], 2))
float getDbkgkinConstant(int ZZflav, float ZZMass){ // ZZflav==id1*id2*id3*id4
  float par[14]={
    0.775,
    -0.565,
    70.,
    5.90,
    -0.235,
    130.1,
    13.25,
    -0.33,
    191.04,
    16.05,
    187.47,
    -0.21,
    1700.,
    400.
  };
  if (abs(ZZflav)==121*121 || abs(ZZflav)==121*242 || abs(ZZflav)==242*242){ par[11]=-0.3; par[13]=500.; } // 4e
  else if (abs(ZZflav)==169*121 || abs(ZZflav)==169*242) par[11]=-0.2; // 2e2mu
  else if (abs(ZZflav)==169*169) par[11]=-0.16; // 4mu
  float kappa =
    par[0]
    +par[1]*exp(-pow(((ZZMass-par[2])/par[3]), 2))
    +par[4]*exp(-pow(((ZZMass-par[5])/par[6]), 2))
    +par[7]*(
    exp(-pow(((ZZMass-par[8])/par[9]), 2))*(ZZMass<par[8])
    + exp(-pow(((ZZMass-par[8])/par[10]), 2))*(ZZMass>=par[8])
    );
  if (abs(ZZflav)==121*121 || abs(ZZflav)==121*242 || abs(ZZflav)==242*242) kappa += par[11]*exp(-pow(((ZZMass<par[12] ? (ZZMass-par[12]) : 0.)/par[13]), 2));
  else kappa += par[11]*exp(-pow((ZZMass-par[12]) /par[13], 2));

  float constant = kappa/(1.-kappa);
  return constant;
}

float getDbkgConstant(int ZZflav, float ZZMass){
  float cbkgkin = getDbkgkinConstant(ZZflav, ZZMass);
  if (abs(ZZflav)==121*121 || abs(ZZflav)==121*242 || abs(ZZflav)==242*242) return cbkgkin*35.6; // 4e
  else if (abs(ZZflav)==169*169) return cbkgkin*22.8; // 4mu
  else if (abs(ZZflav)==121*169 || abs(ZZflav)==242*169) return cbkgkin*41.8; // 2e2mu
  else return 1.;
}

#endif

