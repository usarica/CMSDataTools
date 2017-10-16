#include "SimpleInterferenceDiscriminant.h"
#include "HelperFunctions.h"
#include <cassert>


using namespace std;
using namespace HelperFunctions;


SimpleInterferenceDiscriminant::SimpleInterferenceDiscriminant(const TString cfilename, const TString splinename) : Discriminant(cfilename, splinename){}

void SimpleInterferenceDiscriminant::eval(const std::vector<float>& vars, const float& valReco){
  const unsigned int nvarsreq=3;
  const unsigned int nvarsreq_withC=5;
  bool isWithC=(vars.size()==nvarsreq_withC);
  assert((checkNonNegative(vars) && isWithC) || vars.size()==nvarsreq);
  if (!checkNanInf(vars)) val = -999;
  else{
    float constant = getCval(valReco);
    if (!isWithC) val = vars[2]*sqrt(constant)/(vars[0]+constant*vars[1]);
    else val = (vars[2]*(1./vars[3]+1./vars[4])-vars[0]/vars[3]-vars[1]/vars[4])*sqrt(vars[3]*vars[4]*constant)/(vars[0]+constant*vars[1]);
  }
}
