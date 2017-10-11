#include "CJLSTTree.h"


CJLSTTree::CJLSTTree(TString strsample) :
BaseTree(CJLSTTree::constructCJLSTSamplePath(strsample), TREE_NAME, TREE_FAILED_NAME, COUNTERS_NAME),
sampleIdentifier(strsample),
MHVal(-1),
GHVal(0)
{
  if (valid){
    SampleHelpers::makeGlobalMELA(theSqrts);
    MHVal = SampleHelpers::findPoleMass(sampleIdentifier);
    GHVal = SampleHelpers::GlobalMELA->getHiggsWidthAtPoleMass(MHVal);
  }
}


TString CJLSTTree::constructCJLSTSamplePath(TString strsample){
  TString res = CJLSTsamplesdir + "/" + strsample + "/ZZ4lAnalysis.root";
  return res;
}

float CJLSTTree::getNGenWithPU(){ return (hCounters ? hCounters->GetBinContent(40): 0.); }

