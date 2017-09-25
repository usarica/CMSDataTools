#ifndef CJLSTSET_H
#define CJLSTSET_H

#include <iostream>
#include <vector>
#include "CJLSTTree.h"


class CJLSTSet{
protected:
  std::vector<CJLSTTree*> treeList;
  std::unordered_map<CJLSTTree*, float> permanentWeights;

public:
  CJLSTSet(const TString& strname);
  CJLSTSet(const std::vector<TString>& strlist);
  ~CJLSTSet();

  bool addCJLSTTree(const TString& strname);

  CJLSTTree* getCJLSTTree(TString sampleid);

  // Call this to book xsec, does not change per event
  void bookXS();
  // Call this to book "reco-dependent" event weights (x gen weights), basically product of weights that can change per event
  void bookOverallEventWgt();
  // Call this to store the immutable tree weights (those that depend on xsec and NGenTrue)
  void getPermanentWeights(bool useXS, bool useNormPerMass, bool useNgen, bool renormalizeWeights);

  // Get the overall "reco" weight (no xsec)
  float getOverallEventWgt(CJLSTTree* sample);
  float getOverallEventWgt(TString sampleid);

  // Get the permanent, immutable weights of the different trees
  float getPermanentWeight(CJLSTTree* sample);
  float getPermanentWeight(TString sampleid);

};


#endif
