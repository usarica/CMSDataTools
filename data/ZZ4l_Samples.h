#ifndef ZZ4L_SAMPLES_H
#define ZZ4L_SAMPLES_H


#include <string>
#include "TMath.h"
#include "TString.h"

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
enum CJLST_FinalState{
  MMMM=0,
  EEEE=1,
  EEMM=2,
  ZZ=3,
  LLTT=4, // tau final states (L=e,mu,tau); obsolete.
  llTT=27,// tau final states (L=e,mu)
  TTTT=28,// 4tau
  ZLL=18, // Generic label for a CR (used in job configuration)
  ZZOnShell = 26,
  CRMMMMss=5, // AA CRs (ss loose letptons; os ones are for x-check)
  CRMMMMos=6,
  CREEEEss=7,
  CREEEEos=8,
  CREEMMss=9,
  CREEMMos=10,
  CRMMEEss=11,
  CRMMEEos=12,
  CRZLLss=21,  //merged AA ss CRs
  CRZLLos_2P2F=22,  //CR for Z+2l opposite sign, 2P2F
  CRZLLos_3P1F=23,  //CR for Z+2l opposite sign, 3P1F
  CRZLLos_2P2F_ZZOnShell=24,  //CR for Z+2l opposite sign, 2P2F, 2 on-shell Zs
  CRZLLos_3P1F_ZZOnShell=25,  //CR for Z+2l opposite sign, 3P1F, 2 on-shell Zs
  ZL=13, // Fake rate CR (Z+loose lepton)
  CRZLLHiSIP=14,   // Old inverted SIP CRs
  CRZLLHiSIPMM=15,
  CRZLLHiSIPKin=16,
  CRZLL=17,   // Old CR for Z2 with no SIP
  CRZ2mLL=19, // A CR: Z->mumu + l+l- (with at least 1F)
  CRZ2eLL=20, // A CR: Z->ee   + l+l- (with at least 1F)
  NONE = 99,
  BUGGY=666
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
const TString sample_prefix_MCFM[kNumSamples]={
  "ggTo",
  "ggTo",
  "ggTo"
};
const TString sample_suffix_MCFM[kNumSamples]={
  "_Contin_MCFM701",
  "_Contin_MCFM701",
  "_Contin_MCFM701"
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
//const TString user_gg2VV_location="/scratch0/hep/usarical/CJLST/LHC_13TeV/4l/160714/";
//const TString user_gg2VV_location="/scratch0/hep/usarical/CJLST/LHC_13TeV/4l/160720/";
//const TString user_dir="/scratch0/hep/usarical/CJLST/Analysis/ICHEP2016_mainstream/";
//const TString user_TemplateswithTrees_dir="/scratch0/hep/usarical/CJLST/Analysis/ICHEP2016_mainstream/Templates/";
const TString user_gg2VV_location="/work-zfs/lhc/usarica/hep/CJLST/LHC_13TeV/4l/160720/";
const TString user_dir="/work-zfs/lhc/usarica/CMS-related/Analysis/MassWidth_2016/ICHEP2016_mainstream/";
const TString user_TemplateswithTrees_dir="/work-zfs/lhc/usarica/CMS-related/Analysis/MassWidth_2016/ICHEP2016_mainstream/Templates/";

const TString user_treefile = "ZZ4lAnalysis.root";

const TString user_treename = "ZZTree/candTree";
const TString user_countersname = "ZZTree/Counters";

const TString user_CRtreename = "CRZLLTree/candTree";
const TString user_CRcountersname = "CRZLLTree/Counters";

const TString user_RSE_treename = "ZZTreelooseEle/candTree";

const TString user_RSE_CRtreename = "CRZLLTreelooseEle/candTree";

const TString user_TLE_treename = "ZZTreetle/candTree";

const TString user_TLE_CRtreename = "CRZLLTreetle/candTree";


#endif

