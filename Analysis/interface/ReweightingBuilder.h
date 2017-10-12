#ifndef REWEIGHTINGBUILDER_H
#define REWEIGHTINGBUILDER_H

#include "CJLSTTree.h"


class ReweightingBuilder{
protected:
  float(*rule)(CJLSTTree*, const std::vector<TString>&);
  std::vector<TString> strWeights;

  std::unordered_map<CJLSTTree*, std::vector<float>> weightThresholds;

public:
  ReweightingBuilder(TString inStrWeight, float(*infcn)(CJLSTTree*, const std::vector<TString>&));
  ReweightingBuilder(std::vector<TString> inStrWeights, float(*infcn)(CJLSTTree*, const std::vector<TString>&));

  virtual float eval(CJLSTTree* theTree) const;
  void setupWeightVariables(CJLSTTree* theTree, std::vector<std::pair<float, float>>& boundaries, TString strOrderingVal, float fractionRequirement=0.999);
  std::vector<float> getWeightThresholds(CJLSTTree* theTree);

};


#endif
