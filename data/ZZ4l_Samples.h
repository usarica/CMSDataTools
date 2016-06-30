#include <iostream>
#include <string>
#include "TMath.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TCanvas.h"

const float PI_VAL = TMath::Pi();

enum sample{
  kSig=0,
  kBkg=1,
  kBSI=2,
  kNumSamples=3
};
enum channel{
  k4mu,
  k4e,
  k2e2mu,
  nChannels
};
const TString sampleName[kNumSamples] ={
	"Sig",
	"Bkg",
	"BSI"
};
const TString sample_suffix[kNumSamples]={
	"Sig",
	"Bkg",
	"BSI"
};
const TString sample_suffix_MCFM[kNumSamples]={
  "ggZZ",
  "ggZZ",
  "ggZZ"
};
const TString sample_suffix_Phantom[kNumSamples]={
  "",
  "",
  ""
};
const TString user_folder[nChannels]={
  "4mu",
  "4e",
  "2e2mu"
};

const float luminosity[1]={ 10 };
const double MCFM_widthrescale = 0.9771908764; // Signal g1=sqrt(MCFM_widthrescale)

//Directories
const TString user_gg2VV_location="/work-zfs/lhc/heshy/CJLSTtrees/160121/";
const TString user_dir="/work-zfs/lhc/usarica/CMS-related/Analysis/MassWidth_2016/CMSSW_8_0_12/src/Analysis/ICHEP2016/";
const TString user_TemplateswithTrees_dir="/work-zfs/lhc/usarica/CMS-related/Analysis/MassWidth_2016/CMSSW_8_0_12/src/Analysis/ICHEP2016/Templates/28062016/";

const TString user_treefile = "ZZ4lAnalysis.root";
const TString user_treename = "ZZTree/candTree";
const TString user_countersname = "ZZTree/Counters";

