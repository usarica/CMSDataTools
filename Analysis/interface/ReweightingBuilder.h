#ifndef REWEIGHTINGBUILDER_H
#define REWEIGHTINGBUILDER_H

#include "CJLSTTree.h"
#include "ExtendedBinning.h"


class ReweightingBuilder{
protected:
  float(*rule)(CJLSTTree*, const std::vector<TString>&);
  std::vector<TString> strWeights;

  std::unordered_map<CJLSTTree*, ExtendedBinning> weightBinning;
  std::unordered_map<CJLSTTree*, std::vector<float>> weightThresholds;
  std::unordered_map<CJLSTTree*, std::vector<float>> sumPostThrWeights;
  std::unordered_map<CJLSTTree*, std::vector<unsigned int>> sumEvents;
  std::unordered_map<CJLSTTree*, std::vector<unsigned int>> sumNonZeroWgtEvents;

public:
  ReweightingBuilder(TString inStrWeight, float(*infcn)(CJLSTTree*, const std::vector<TString>&));
  ReweightingBuilder(std::vector<TString> inStrWeights, float(*infcn)(CJLSTTree*, const std::vector<TString>&));

  virtual float eval(CJLSTTree* theTree) const;
  float getPostThresholdWeight(CJLSTTree* theTree) const;
  void setupWeightVariables(CJLSTTree* theTree, const ExtendedBinning& binning, float fractionRequirement=0.999);
  std::vector<float> getWeightThresholds(CJLSTTree* theTree);

};


#endif
