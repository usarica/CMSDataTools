#include "VHProdIntTrigPhaseACDiscriminant.h"
#include "HelperFunctions.h"
#include <cassert>


using namespace std;
using namespace HelperFunctions;


VHProdIntTrigPhaseACDiscriminant::VHProdIntTrigPhaseACDiscriminant() : Discriminant("", "", "", "", 1){}

void VHProdIntTrigPhaseACDiscriminant::eval(const std::vector<float>& vars, const float& /*valReco*/){
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
  if (
    checkNanInf(vars)
    &&
    checkNonNegative(vars, -1, (int) iZHBSMSMInt) && checkNonNegative(vars, (int) iZHConst, (int) iWHBSMSMInt) && checkNonNegative(vars, (int) iWHConst, -1)
    &&
    vars.size()==(unsigned int) nvarsreq
    ){
    val=0;
    float pZHSM = vars[iZHSM]/vars[iZHConst];
    float pZHBSM = vars[iZHBSM]/vars[iZHConst];
    float pZHBSMSMInt = vars[iZHBSMSMInt]/vars[iZHConst];
    float pWHSM = vars[iWHSM]/vars[iWHConst];
    float pWHBSM = vars[iWHBSM]/vars[iWHConst];
    float pWHBSMSMInt = vars[iWHBSMSMInt]/vars[iWHConst];

    unsigned int count=0;
    if (pZHSM!=0. && pZHBSM!=0.){ val += pZHBSMSMInt / (2.*sqrt(pZHSM*pZHBSM)); count++; }
    if (pWHSM!=0. && pWHBSM!=0.){ val += pWHBSMSMInt / (2.*sqrt(pWHSM*pWHBSM)); count++; }
    if (count>0) val /= static_cast<float>(count);
  }
}
