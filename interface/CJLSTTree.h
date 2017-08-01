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


CJLSTTree::CJLSTTree(TString strsample) :
BaseTree(CJLSTTree::constructCJLSTSamplePath(strsample), TREE_NAME, TREE_FAILED_NAME, COUNTERS_NAME),
sampleIdentifier(strsample)
{}


TString CJLSTTree::constructCJLSTSamplePath(TString strsample){
  TString res = CJLSTsamplesdir + "/" + strsample + "/ZZ4lAnalysis.root";
  return res;
}

float CJLSTTree::getNGenWithPU(){ return (hCounters ? hCounters->GetBinContent(40): 0.); }

#endif
