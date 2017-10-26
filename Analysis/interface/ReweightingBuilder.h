#ifndef REWEIGHTINGBUILDER_H
#define REWEIGHTINGBUILDER_H

#include "CJLSTTree.h"
#include "ExtendedBinning.h"


class ReweightingBuilder{
protected:
  float(*rule)(CJLSTTree*, const std::vector<TString>&);
  std::vector<TString> strWeights;

  ExtendedBinning weightBinning;
  std::unordered_map<CJLSTTree*, std::vector<float>> weightThresholds;
  std::unordered_map<CJLSTTree*, std::vector<float>> sumPostThrWeights;
  std::unordered_map<CJLSTTree*, std::vector<unsigned int>> sumEvents;
  std::unordered_map<CJLSTTree*, std::vector<unsigned int>> sumNonZeroWgtEvents;

public:
  ReweightingBuilder(TString inStrWeight, float(*infcn)(CJLSTTree*, const std::vector<TString>&));
  ReweightingBuilder(std::vector<TString> inStrWeights, float(*infcn)(CJLSTTree*, const std::vector<TString>&));

  virtual float eval(CJLSTTree* theTree) const;

  std::vector<float> getWeightThresholds(CJLSTTree* theTree) const;
  float getPostThresholdWeight(CJLSTTree* theTree) const;
  float getPostThresholdSumWeights(CJLSTTree* theTree) const;
  unsigned int getSumEvents(CJLSTTree* theTree) const;
  unsigned int getSumNonZeroWgtEvents(CJLSTTree* theTree) const;
  int findBin(CJLSTTree* theTree) const;

  void setWeightBinning(const ExtendedBinning& binning);
  void setupWeightVariables(CJLSTTree* theTree, float fractionRequirement=0.999);

};


#endif
