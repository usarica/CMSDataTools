#include "VHProdIntTrigPhaseACDiscriminant.h"
#include "HelperFunctions.h"
#include <cassert>


using namespace std;
using namespace HelperFunctions;


VHProdIntTrigPhaseACDiscriminant::VHProdIntTrigPhaseACDiscriminant() : Discriminant("", "", "", "", 1){}

void VHProdIntTrigPhaseACDiscriminant::eval(const std::vector<float>& vars, const float& valReco){
  enum{
    iZHSM=0,
    iZHBSM,
    iZHBSMSMInt, // Should already be subtracted
    iZHConst,
    iWHSM,
    iWHBSM,
    iWHBSMSMInt, // Should already be subtracted
    iWHConst,
    nvarsreq
  };
  this->resetVal();
  if (checkNanInf(vars) && checkNonNegative(vars) && vars.size()==(unsigned int) nvarsreq){
    float pZHSM = vars[iZHSM]/vars[iZHConst];
    float pWHSM = vars[iWHSM]/vars[iWHConst];
    float pZHBSM = vars[iZHBSM]/vars[iZHConst];
    float pWHBSM = vars[iWHBSM]/vars[iWHConst];
    float pZHBSMSMInt = vars[iZHBSMSMInt]/vars[iZHConst];
    float pWHBSMSMInt = vars[iWHBSMSMInt]/vars[iWHConst];
    float pVHSM = (pZHSM + pWHSM);
    float pVHBSM = (pZHBSM + pWHBSM);
    float pVHBSMSMInt = (pZHBSMSMInt + pWHBSMSMInt);

    val = pVHBSMSMInt / sqrt(pVHSM*pVHBSM);
  }
}
