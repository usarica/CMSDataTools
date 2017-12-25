#ifndef REWEIGHTINGFUNCTIONS_H
#define REWEIGHTINGFUNCTIONS_H

#include "CJLSTTree.h"
#include <vector>
#include "TString.h"
#include "HelperFunctions.h"


namespace ReweightingFunctions{
  std::vector<float*> getWeightRefs(CJLSTTree* tree, const std::vector<TString>& strWeights);
  float* getWeightRef(CJLSTTree* tree, const TString& strWeight);

  typedef float(*ReweightingFunction_t)(CJLSTTree*, const std::vector<float*>&);
  float getSimpleWeight(CJLSTTree* tree, const std::vector<float*>& vals);
  float getA1PlusB1Weight(CJLSTTree* tree, const std::vector<float*>& vals);
  float getA1MinusB1Weight(CJLSTTree* tree, const std::vector<float*>& vals);
  float getOnePlusB1OverA1Weight(CJLSTTree* tree, const std::vector<float*>& vals);
  float getOneMinusB1OverA1Weight(CJLSTTree* tree, const std::vector<float*>& vals);
  float getA1OverB1Weight(CJLSTTree* tree, const std::vector<float*>& vals);

}


#endif
