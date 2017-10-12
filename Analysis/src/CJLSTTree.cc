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
float CJLSTTree::getTrueBW(){
  if (MHVal<0. || GHVal==0.) return 1.;
  else{
    if (valfloats.find("GenHMass")==valfloats.end()){ // Don't spend time trying to find the branch
      float valdef=-1;
      this->bookBranch("GenHMass", valdef);
      this->refreshCurrentEvent();
    }
    float& GenHMass = valfloats["GenHMass"]->first;
    if (GenHMass>=0.){
      float invProp = pow(pow(GenHMass, 2)-pow(MHVal, 2), 2) + pow(MHVal*GHVal, 2);
      return 1./invProp;
    }
    else return 1.;
  }
}

