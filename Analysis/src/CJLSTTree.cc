#include "CJLSTTree.h"


CJLSTTree::CJLSTTree(TString strsample) :
  BaseTree(CJLSTTree::constructCJLSTSamplePath(strsample), TREE_NAME, TREE_FAILED_NAME, COUNTERS_NAME),
  associatedSet(nullptr),
  sampleIdentifier(strsample),
  MHVal(-1),
  GHVal(0)
{
  if (valid){
    //SampleHelpers::makeGlobalMELA(theSqrts);
    MHVal = SampleHelpers::findPoleMass(sampleIdentifier);
    //GHVal = SampleHelpers::GlobalMELA->getHiggsWidthAtPoleMass(MHVal);
  }
}

TString CJLSTTree::constructCJLSTSamplePath(TString strsample){
  TString res = CJLSTsamplesdir + "/" + strsample + "/ZZ4lAnalysis.root";
  return res;
}

unsigned int CJLSTTree::getNGenNoPU(){ return (hCounters ? hCounters->GetBinContent(1) : 0.); }
float CJLSTTree::getNGenWithPU(){ return (hCounters ? hCounters->GetBinContent(40): 0.); }

float CJLSTTree::getTrueBW(float const* overrideMH){
  if ((!overrideMH && (MHVal<0. || GHVal==0.)) || (overrideMH && (*overrideMH)<0.)) return 1.;
  else{
    if (valfloats.find("GenHMass")==valfloats.end()){ // Don't spend time trying to find the branch
      float valdef=-1;
      this->bookBranch("GenHMass", valdef);
      this->refreshCurrentEvent();
    }
    float& GenHMass = valfloats["GenHMass"]->first;
    if (GenHMass>=0.){
      float invProp;
      if (overrideMH) invProp = pow(pow(GenHMass, 2)-pow(*overrideMH, 2), 2) + pow((*overrideMH)*SampleHelpers::GlobalMELA->getHiggsWidthAtPoleMass(*overrideMH), 2);
      else invProp = pow(pow(GenHMass, 2)-pow(MHVal, 2), 2) + pow(MHVal*GHVal, 2);
      return 1./invProp;
    }
    else return 1.;
  }
}

void CJLSTTree::setAssociatedSet(CJLSTSet* inSet){ associatedSet = inSet; }

