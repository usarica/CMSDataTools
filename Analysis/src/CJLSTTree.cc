#include "CJLSTTree.h"


CJLSTTree::CJLSTTree(TString strsample) :
  BaseTree(CJLSTTree::constructCJLSTSamplePath(strsample), TREE_NAME, TREE_FAILED_NAME, COUNTERS_NAME),
  associatedSet(nullptr),
  sampleIdentifier(strsample),
  MHVal(-1)
{
  if (valid) MHVal = SampleHelpers::findPoleMass(sampleIdentifier);
}

TString CJLSTTree::constructCJLSTSamplePath(TString strsample){
  TString res = CJLSTsamplesdir + "/" + strsample + "/ZZ4lAnalysis.root";
  return res;
}

unsigned int CJLSTTree::getNGenNoPU(){ return (hCounters ? hCounters->GetBinContent(1) : 0.); }
float CJLSTTree::getNGenWithPU(){ return (hCounters ? hCounters->GetBinContent(40): 0.); }
