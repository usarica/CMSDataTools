#ifndef REWEIGHTINGBUILDER_H
#define REWEIGHTINGBUILDER_H

#include "CJLSTTree.h"


class ReweightingBuilder{
public:
  enum BWScheme{
    kDoNotConsider,
    kDivideSampleBW,
    kMultiplySampleBW
  };

protected:
  const BWScheme theScheme;
  std::vector<TString> strWeights;

  std::unordered_map<CJLSTTree*, std::vector<float>> weightThresholds;

public:
  ReweightingBuilder(TString inStrWeight, ReweightingBuilder::BWScheme inscheme = ReweightingBuilder::kDoNotConsider);
  ReweightingBuilder(std::vector<TString> inStrWeights, ReweightingBuilder::BWScheme inscheme = ReweightingBuilder::kDoNotConsider);

  virtual float eval(CJLSTTree* theTree) const;
  std::vector<float> determineWeightThresholds(CJLSTTree* theTree, std::vector<std::pair<float, float>>& boundaries, float fractionRequirement=0.999, TString orderingVal="GenHMass");
  std::vector<float> getWeightThresholds(CJLSTTree* theTree);

};


#endif
