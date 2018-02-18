#include "CJLSTTree.h"


CJLSTTree::CJLSTTree(TString strsample, const TString treename, const TString failedtreename, const TString countersname) :
  BaseTree(CJLSTTree::constructCJLSTSamplePath(strsample), treename, failedtreename, countersname),
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
