#ifndef CJLSTTREE_H
#define CJLSTTREE_H

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

  CJLSTTree(TString strsample);
  ~CJLSTTree(){}

  unsigned int getNGenNoPU();
  float getNGenWithPU();

  void setAssociatedSet(CJLSTSet* inSet){ associatedSet = inSet; }
  CJLSTSet* getAssociatedSet(){ return associatedSet; }
  CJLSTSet const* getAssociatedSet() const{ return associatedSet; }

};


#endif
