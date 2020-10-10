#ifndef REWEIGHTINGFUNCTIONS_H
#define REWEIGHTINGFUNCTIONS_H

#include <vector>
#include "BaseTree.h"
#include "HelperFunctions.h"
#include "TVar.hh"


namespace ReweightingFunctions{
  std::vector<float*> getWeightRefs(BaseTree* tree, std::vector<TString> const& strWeights);
  float* getWeightRef(BaseTree* tree, TString const& strWeight);

  typedef float(*ReweightingFunction_t)(BaseTree*, std::vector<float*> const&);
  float getSimpleWeight(BaseTree* tree, std::vector<float*> const& vals); // Product of weights
  float getA1PlusB1Weight(BaseTree* tree, std::vector<float*> const& vals); // wgt = wA + wB
  float getA1MinusB1Weight(BaseTree* tree, std::vector<float*> const& vals); // wgt = wA - wB
  float getA1MinusB1MinusC1Weight(BaseTree* tree, std::vector<float*> const& vals); // wgt = wA - wB - wC
  float getOnePlusB1OverA1Weight(BaseTree* tree, std::vector<float*> const& vals); // wgt = 1 + wB/wA
  float getOneMinusB1OverA1Weight(BaseTree* tree, std::vector<float*> const& vals); // wgt = 1 - wB/wA
  float getA1OverB1Weight(BaseTree* tree, std::vector<float*> const& vals); // wgt = wA/wB

  float getAbsWeightThresholdByNeff(BaseTree* tree, std::vector<float*> const& vals, ReweightingFunction_t rule, double thr_Neff, TVar::VerbosityLevel verbosity=TVar::ERROR);

}


#endif
