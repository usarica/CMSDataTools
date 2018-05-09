#ifndef REWEIGHTINGBUILDER_H
#define REWEIGHTINGBUILDER_H

#include "CJLSTTree.h"
#include "ExtendedBinning.h"
#include "ReweightingFunctions.h"


class ReweightingBuilder{
public:
  constexpr static bool useNeffInNormComponent=true;

protected:
  bool allowNegativeWeights;
  bool divideByNSample;
  ReweightingFunctions::ReweightingFunction_t rule;
  std::vector<TString> strWeights;

  ExtendedBinning weightBinning;
  std::unordered_map<CJLSTTree*, float*> binningVarRefs;
  std::unordered_map<CJLSTTree*, std::vector<float*>> componentRefs;
  std::unordered_map<CJLSTTree*, std::vector<float>> weightThresholds;
  std::unordered_map<CJLSTTree*, std::vector<float>> sumPostThrWeights;
  std::unordered_map<CJLSTTree*, std::vector<float>> sumPostThrSqWeights;
  std::unordered_map<CJLSTTree*, std::vector<unsigned int>> sumEvents;
  std::unordered_map<CJLSTTree*, std::vector<unsigned int>> sumNonZeroWgtEvents;

public:
  ReweightingBuilder(TString inStrWeight, ReweightingFunctions::ReweightingFunction_t infcn);
  ReweightingBuilder(std::vector<TString> inStrWeights, ReweightingFunctions::ReweightingFunction_t infcn);

  std::vector<TString> const& getWeightVariables() const;

  virtual float eval(CJLSTTree* theTree) const;

  std::vector<float> getWeightThresholds(CJLSTTree* theTree) const;
  float getPostThresholdWeight(CJLSTTree* theTree) const;
  float getPostThresholdSqWeight(CJLSTTree* theTree) const;
  float getPostThresholdSqWeightInv(CJLSTTree* theTree) const;
  float getSumPostThresholdWeights(CJLSTTree* theTree) const;
  float getSumPostThresholdSqWeights(CJLSTTree* theTree) const;
  float getSumPostThresholdSqWeightInvs(CJLSTTree* theTree) const;
  float getSumPostThresholdNeffs(CJLSTTree* theTree) const;
  unsigned int getSumEvents(CJLSTTree* theTree) const;
  unsigned int getSumNonZeroWgtEvents(CJLSTTree* theTree) const;
  int findBin(CJLSTTree* theTree) const;

  void rejectNegativeWeights(const bool flag);
  void setDivideByNSample(const bool flag);
  void setWeightBinning(const ExtendedBinning& binning);
  void setupWeightVariables(CJLSTTree* theTree, float fractionRequirement=0.999, unsigned int minimumNevents=0);
  std::vector<CJLSTTree*> getRegisteredTrees() const;

  float getSumAllPostThresholdWeights(CJLSTTree* theTree) const; // Tree is passed here to find the bin
  float getSumAllPostThresholdWeights(int bin) const;

  float getSumAllPostThresholdSqWeights(CJLSTTree* theTree) const; // Tree is passed here to find the bin
  float getSumAllPostThresholdSqWeights(int bin) const;

  float getSumAllPostThresholdSqWeightInvs(CJLSTTree* theTree) const; // Tree is passed here to find the bin
  float getSumAllPostThresholdSqWeightInvs(int bin) const;

  float getSumAllPostThresholdNeffs(CJLSTTree* theTree) const; // Tree is passed here to find the bin
  float getSumAllPostThresholdNeffs(int bin) const;

  unsigned int getSumAllEvents(CJLSTTree* theTree) const; // Tree is passed here to find the bin
  unsigned int getSumAllEvents(int bin) const;

  unsigned int getSumAllNonZeroWgtEvents(CJLSTTree* theTree) const; // Tree is passed here to find the bin
  unsigned int getSumAllNonZeroWgtEvents(int bin) const;

  float getNormComponent(CJLSTTree* theTree) const; // Tree is passed here to find the bin
  float getNormComponent(int bin) const;
  float getNorm() const;

  static ExtendedBinning getTrueMassBinning(std::vector<CJLSTTree*> const& trees, float forceMass=-1);
};


#endif
