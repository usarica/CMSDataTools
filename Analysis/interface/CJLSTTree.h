#ifndef CJLSTTREE_H
#define CJLSTTREE_H

#include "BaseTree.h"


class CJLSTTree : public BaseTree{
public:
  const TString sampleIdentifier;
  float MHVal;
  float GHVal;

  static TString constructCJLSTSamplePath(TString strsample);

  CJLSTTree(TString strsample);
  ~CJLSTTree(){}

  float getNGenWithPU();

  float getTrueBW(float const* overrideMH=nullptr);

};


#endif
