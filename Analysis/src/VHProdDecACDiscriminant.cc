#include "VHProdDecACDiscriminant.h"
#include "HelperFunctions.h"
#include <cassert>


using namespace std;
using namespace HelperFunctions;


VHProdDecACDiscriminant::VHProdDecACDiscriminant(
  const TString cfilename, const TString splinename,
  const TString gfilename, const TString gsplinename,
  const float gscale_
) : Discriminant(cfilename, splinename, gfilename, gsplinename, gscale_){}

void VHProdDecACDiscriminant::eval(const std::vector<float>& vars, const float& valReco){
  enum{
    iZHSM=0,
    iWHSM,
    iDecSM,
    iZHBSM,
    iWHBSM,
    iDecBSM,
    iZHConst,
    iWHConst,
    nvarsreq
  };
  enum{
    iGVH, // There ought to be a dedicated g for WH+ZH
    iGDec,
    nGs
  };
  this->resetVal();
  if (checkNanInf(vars) && checkNonNegative(vars) && vars.size()==(unsigned int) nvarsreq && theG.size()==(unsigned int)nGs){
    float pZHSM = vars[iZHSM]/vars[iZHConst];
    float pWHSM = vars[iWHSM]/vars[iWHConst];
    float pZHBSM = vars[iZHBSM]/vars[iZHConst];
    float pWHBSM = vars[iWHBSM]/vars[iWHConst];
    float pVHdecSM = (pZHSM + pWHSM)*vars[iDecSM];
    float pVHdecBSM = (pZHBSM + pWHBSM)*vars[iDecBSM];

    float gCommon = pow((theG.at(0).second->Eval(valReco))*(theG.at(1).second->Eval(valReco))*gscale, 2);

    val = pVHdecSM / (pVHdecSM + gCommon*pVHdecBSM);
  }
}
