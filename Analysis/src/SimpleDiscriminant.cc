#include "SimpleDiscriminant.h"
#include "HelperFunctions.h"
#include <cassert>


using namespace std;
using namespace HelperFunctions;


SimpleDiscriminant::SimpleDiscriminant(const TString cfilename, const TString splinename) : Discriminant(cfilename, splinename){}

void SimpleDiscriminant::eval(const std::vector<float>& vars, const float& valReco){
  const unsigned int nvarsreq=2;
  assert(checkNonZero(vars) && vars.size()==nvarsreq);
  if (!checkNanInf(vars)) val = -999;
  else{
    float constant = getCval(valReco);
    val = vars[0]/(vars[0]+constant*vars[1]);
  }
}
