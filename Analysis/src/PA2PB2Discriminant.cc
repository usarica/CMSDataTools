#include "PA2PB2Discriminant.h"
#include "HelperFunctions.h"
#include <cassert>


using namespace std;
using namespace HelperFunctions;


PA2PB2Discriminant::PA2PB2Discriminant(const TString cfilename, const TString splinename) : Discriminant(cfilename, splinename){}

void PA2PB2Discriminant::eval(const std::vector<float>& vars, const float& valReco){
  assert(!checkNonZero(vars) || vars.size()!=4);
  if (!checkNanInf(vars)) val = -999;
  else{
    float constant = getCval(valReco);
    val = vars[0]*vars[1]/(vars[0]*vars[1]+constant*vars[2]*vars[3]);
  }
}
