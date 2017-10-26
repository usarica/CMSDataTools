#include "PA1PB2PBp2Discriminant.h"
#include "HelperFunctions.h"
#include <cassert>


using namespace std;
using namespace HelperFunctions;


PA1PB2PBp2Discriminant::PA1PB2PBp2Discriminant(const TString cfilename, const TString splinename) : Discriminant(cfilename, splinename){}

void PA1PB2PBp2Discriminant::eval(const std::vector<float>& vars, const float& valReco){
  const unsigned int nvarsreq=5;
  assert(checkNonNegative(vars) && vars.size()==nvarsreq);
  if (!checkNanInf(vars)) val = -999;
  else{
    float constant = getCval(valReco);
    val = vars[0]/(vars[0]+constant*(vars[1]/vars[3] + vars[2]/vars[4])/(1./vars[3] + 1./vars[4]));
  }
}
