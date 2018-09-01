#ifndef CJLSTSET_H
#define CJLSTSET_H

#include <iostream>
#include <vector>
#include "CJLSTTree.h"


class CJLSTSet{
public:
  enum NormScheme{
    NormScheme_None,
    NormScheme_OneOverNgen,
    NormScheme_OneOverNgen_RenormBySumOneOverNgen,
    NormScheme_OneOverNgen_RelRenormToSumNgen,
    NormScheme_NgenOverNgenWPU,
    NormScheme_XsecOnly,
    NormScheme_XsecOverNgen,
    NormScheme_XsecOverNgen_RenormBySumXsecOverNgen,
    NormScheme_XsecOverNgen_RelRenormToSumNgen
  };

protected:
  std::vector<CJLSTTree*> treeList;
  std::unordered_map<CJLSTTree*, float> permanentWeights;

public:
  CJLSTSet(const TString& strname, const TString treename=TREE_NAME, const TString failedtreename=TREE_FAILED_NAME, const TString countersname=COUNTERS_NAME);
  CJLSTSet(const std::vector<TString>& strlist, const TString treename=TREE_NAME, const TString failedtreename=TREE_FAILED_NAME, const TString countersname=COUNTERS_NAME);
  ~CJLSTSet();

  bool addCJLSTTree(const TString& strname, const TString treename=TREE_NAME, const TString failedtreename=TREE_FAILED_NAME, const TString countersname=COUNTERS_NAME);
  bool addCJLSTTreeList(const std::vector<TString>& strlist, const TString treename=TREE_NAME, const TString failedtreename=TREE_FAILED_NAME, const TString countersname=COUNTERS_NAME);
  bool dissociateCJLSTTree(CJLSTTree*& tree);
  bool associateCJLSTTree(CJLSTTree*& tree);

  CJLSTTree* getCJLSTTree(TString sampleid) const;
  const std::vector<CJLSTTree*>& getCJLSTTreeList() const;
  std::vector<CJLSTTree*>& getCJLSTTreeList();

  // Call this to book xsec, does not change per event
  void bookXS();
  // Call this to book "reco-dependent" event weights (x gen weights), basically product of weights that can change per event
  void bookOverallEventWgt();
  // Call this to store the immutable tree weights (those that depend on xsec and NGenTrue)
  void setPermanentWeights(const CJLSTSet::NormScheme scheme, const bool useNormPerMass, const bool useNgenWPU=false);

  // Get the overall "reco" weight (no xsec)
  float getOverallEventWgt(CJLSTTree* const sample) const;
  float getOverallEventWgt(TString sampleid) const;

  // Get the permanent, immutable weights of the different trees
  float getPermanentWeight(CJLSTTree* const sample) const;
  float getPermanentWeight(TString sampleid) const;

  CJLSTTree* getSelectedEvent(const int evid);
  CJLSTTree* getFailedEvent(const int evid);

};


#endif
