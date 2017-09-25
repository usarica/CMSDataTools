#ifndef CJLSTTREE_H
#define CJLSTTREE_H

#include "BaseTree.h"


class CJLSTTree : public BaseTree{
public:
  const TString sampleIdentifier;

  static TString constructCJLSTSamplePath(TString strsample);

  CJLSTTree(TString strsample);
  ~CJLSTTree(){}

  float getNGenWithPU();

};


#endif
