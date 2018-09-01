#ifndef CJLSTTREE_H
#define CJLSTTREE_H

#include "Samples.h"
#include "BaseTree.h"


// Forward declarations
class CJLSTSet;


class CJLSTTree : public BaseTree{
protected:
  CJLSTSet* associatedSet;

public:
  const TString sampleIdentifier;
  float MHVal;

  static TString constructCJLSTSamplePath(TString strsample);

  CJLSTTree(TString strsample, const TString treename=TREE_NAME, const TString failedtreename=TREE_FAILED_NAME, const TString countersname=COUNTERS_NAME);
  ~CJLSTTree(){}

  unsigned int getNGenNoPU();
  float getNGenWithPU();

  void setAssociatedSet(CJLSTSet* inSet){ associatedSet = inSet; }
  CJLSTSet* getAssociatedSet(){ return associatedSet; }
  CJLSTSet const* getAssociatedSet() const{ return associatedSet; }

};


#endif
