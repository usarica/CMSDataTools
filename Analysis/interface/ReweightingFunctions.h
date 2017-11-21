#ifndef REWEIGHTINGFUNCTIONS_H
#define REWEIGHTINGFUNCTIONS_H

#include "CJLSTTree.h"
#include <vector>
#include "TString.h"
#include "HelperFunctions.h"


namespace ReweightingFunctions{
  std::vector<float*> getWeightRefs(CJLSTTree* tree, const std::vector<TString>& strWeights);
  float* getWeightRef(CJLSTTree* tree, const TString& strWeight);

  float getSimpleWeight(CJLSTTree* tree, const std::vector<float*>& vals);

}


#endif
