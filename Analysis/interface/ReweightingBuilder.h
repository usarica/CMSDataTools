#ifndef REWEIGHTINGBUILDER_H
#define REWEIGHTINGBUILDER_H

#include "CJLSTTree.h"
#include "ExtendedBinning.h"


class ReweightingBuilder{
protected:
  bool divideByNSample;
  float(*rule)(CJLSTTree*, const std::vector<float*>&);
  std::vector<TString> strWeights;

  ExtendedBinning weightBinning;
  std::unordered_map<CJLSTTree*, float*> binningVarRefs;
  std::unordered_map<CJLSTTree*, std::vector<float*>> componentRefs;
  std::unordered_map<CJLSTTree*, std::vector<float>> weightThresholds;
  std::unordered_map<CJLSTTree*, std::vector<float>> sumPostThrWeights;
  std::unordered_map<CJLSTTree*, std::vector<unsigned int>> sumEvents;
  std::unordered_map<CJLSTTree*, std::vector<unsigned int>> sumNonZeroWgtEvents;

public:
  ReweightingBuilder(TString inStrWeight, float(*infcn)(CJLSTTree*, const std::vector<float*>&));
  ReweightingBuilder(std::vector<TString> inStrWeights, float(*infcn)(CJLSTTree*, const std::vector<float*>&));

  std::vector<TString> const& getWeightVariables() const;

  virtual float eval(CJLSTTree* theTree) const;

  std::vector<float> getWeightThresholds(CJLSTTree* theTree) const;
  float getPostThresholdWeight(CJLSTTree* theTree) const;
  float getSumPostThresholdWeights(CJLSTTree* theTree) const;
  unsigned int getSumEvents(CJLSTTree* theTree) const;
  unsigned int getSumNonZeroWgtEvents(CJLSTTree* theTree) const;
  int findBin(CJLSTTree* theTree) const;

  void setDivideByNSample(const bool flag);
  void setWeightBinning(const ExtendedBinning& binning);
  void setupWeightVariables(CJLSTTree* theTree, float fractionRequirement=0.999, unsigned int minimumNevents=0);
  std::vector<CJLSTTree*> getRegisteredTrees() const;

  float getSumAllPostThresholdWeights(CJLSTTree* theTree) const; // Tree is passed here to find the bin
  float getSumAllPostThresholdWeights(int bin) const;

  unsigned int getSumAllEvents(CJLSTTree* theTree) const; // Tree is passed here to find the bin
  unsigned int getSumAllEvents(int bin) const;

  unsigned int getSumAllNonZeroWgtEvents(CJLSTTree* theTree) const; // Tree is passed here to find the bin
  unsigned int getSumAllNonZeroWgtEvents(int bin) const;

  float getNormComponent(CJLSTTree* theTree) const; // Tree is passed here to find the bin
  float getNormComponent(int bin) const;
  float getNorm() const;
};


#endif
