#ifndef REWEIGHTINGFUNCTIONS_H
#define REWEIGHTINGFUNCTIONS_H

#include "CJLSTTree.h"
#include <vector>
#include "TString.h"
#include "HelperFunctions.h"


namespace ReweightingFunctions{
  std::vector<float> getWeightValues(CJLSTTree* tree, const std::vector<TString>& strWeights);

  float getSimpleWeight(CJLSTTree* tree, const std::vector<TString>& strWeights);

}


#endif
