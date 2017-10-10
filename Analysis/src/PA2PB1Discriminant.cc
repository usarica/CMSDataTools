#include "PA2PB1Discriminant.h"
#include "HelperFunctions.h"
#include <cassert>


using namespace std;
using namespace HelperFunctions;


PA2PB1Discriminant::PA2PB1Discriminant(const TString cfilename, const TString splinename) : Discriminant(cfilename, splinename){}

void PA2PB1Discriminant::eval(const std::vector<float>& vars, const float& valReco){
  const unsigned int nvarsreq=3;
  assert(checkNonZero(vars) && vars.size()==nvarsreq);
  if (!checkNanInf(vars)) val = -999;
  else{
    float constant = getCval(valReco);
    val = vars[0]*vars[1]/(vars[0]*vars[1]+constant*vars[2]);
  }
}
