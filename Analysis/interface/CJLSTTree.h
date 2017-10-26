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
  float GHVal;

  static TString constructCJLSTSamplePath(TString strsample);

  CJLSTTree(TString strsample);
  ~CJLSTTree(){}

  float getNGenWithPU();
  float getTrueBW(float const* overrideMH=nullptr);

  void setAssociatedSet(CJLSTSet* inSet);
  CJLSTSet* getAssociatedSet(){ return associatedSet; }
  CJLSTSet const* getAssociatedSet() const{ return associatedSet; }

};


#endif
