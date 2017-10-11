#ifndef ZZ4L_SAMPLES_H
#define ZZ4L_SAMPLES_H


#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include "TMath.h"
#include "TString.h"

const float PI_VAL = TMath::Pi();

enum sample{
  kSig=0,
  kBkg=1,
  kBSI=2,
  kNumSamples=3
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

enum Channel{
  k4mu,
  k4e,
  k2e2mu,
  nChannels
};
TString channame(const int ichan){
  if (ichan==(const int)Channel::k4mu) return TString("4mu");
  else if (ichan==(const int)Channel::k4e) return TString("4e");
  else if (ichan==(const int)Channel::k2e2mu) return TString("2e2mu");
  else return TString("4l");
}

enum SampleType{ // Unless specified otherwise, the samples are from POWHEG
  gg_Sig,
  gg_Sig_PythiaScaleDn,
  gg_Sig_PythiaScaleUp,
  gg_Sig_PythiaTuneDn,
  gg_Sig_PythiaTuneUp,

  gg_Sig_MiNLO,
  gg_Sig_PythiaScaleDn_MiNLO,
  gg_Sig_PythiaScaleUp_MiNLO,
  gg_Sig_PythiaTuneDn_MiNLO,
  gg_Sig_PythiaTuneUp_MiNLO,

  VBF_Sig,
  VBF_Sig_PythiaScaleDn,
  VBF_Sig_PythiaScaleUp,
  VBF_Sig_PythiaTuneDn,
  VBF_Sig_PythiaTuneUp,

  WH_Sig,
  WH_Sig_PythiaScaleDn,
  WH_Sig_PythiaScaleUp,
  WH_Sig_PythiaTuneDn,
  WH_Sig_PythiaTuneUp,

  ZH_Sig,
  ZH_Sig_PythiaScaleDn,
  ZH_Sig_PythiaScaleUp,
  ZH_Sig_PythiaTuneDn,
  ZH_Sig_PythiaTuneUp,

  tt_Sig,
  tt_Sig_PythiaScaleDn,
  tt_Sig_PythiaScaleUp,
  tt_Sig_PythiaTuneDn,
  tt_Sig_PythiaTuneUp,

  gg_Bkg_MCFM,

  VBF_Sig_Phantom,
  VBF_Bkg_Phantom,

  qq_Bkg,

  nSampleTypes
};
TString nameSampleType(int is){
  if (is==(int)SampleType::gg_Sig) return TString("gg_Sig");
  else if (is==(int)SampleType::gg_Sig_PythiaScaleDn) return TString("gg_Sig_PythiaScaleDn");
  else if (is==(int)SampleType::gg_Sig_PythiaScaleUp) return TString("gg_Sig_PythiaScaleUp");
  else if (is==(int)SampleType::gg_Sig_PythiaTuneDn) return TString("gg_Sig_PythiaTuneDn");
  else if (is==(int)SampleType::gg_Sig_PythiaTuneUp) return TString("gg_Sig_PythiaTuneUp");

  else if (is==(int)SampleType::gg_Sig_MiNLO) return TString("gg_Sig_MiNLO");
  else if (is==(int)SampleType::gg_Sig_PythiaScaleDn_MiNLO) return TString("gg_Sig_PythiaScaleDn_MiNLO");
  else if (is==(int)SampleType::gg_Sig_PythiaScaleUp_MiNLO) return TString("gg_Sig_PythiaScaleUp_MiNLO");
  else if (is==(int)SampleType::gg_Sig_PythiaTuneDn_MiNLO) return TString("gg_Sig_PythiaTuneDn_MiNLO");
  else if (is==(int)SampleType::gg_Sig_PythiaTuneUp_MiNLO) return TString("gg_Sig_PythiaTuneUp_MiNLO");

  else if (is==(int)SampleType::VBF_Sig) return TString("VBF_Sig");
  else if (is==(int)SampleType::VBF_Sig_PythiaScaleDn) return TString("VBF_Sig_PythiaScaleDn");
  else if (is==(int)SampleType::VBF_Sig_PythiaScaleUp) return TString("VBF_Sig_PythiaScaleUp");
  else if (is==(int)SampleType::VBF_Sig_PythiaTuneDn) return TString("VBF_Sig_PythiaTuneDn");
  else if (is==(int)SampleType::VBF_Sig_PythiaTuneUp) return TString("VBF_Sig_PythiaTuneUp");

  else if (is==(int)SampleType::WH_Sig) return TString("WH_Sig");
  else if (is==(int)SampleType::WH_Sig_PythiaScaleDn) return TString("WH_Sig_PythiaScaleDn");
  else if (is==(int)SampleType::WH_Sig_PythiaScaleUp) return TString("WH_Sig_PythiaScaleUp");
  else if (is==(int)SampleType::WH_Sig_PythiaTuneDn) return TString("WH_Sig_PythiaTuneDn");
  else if (is==(int)SampleType::WH_Sig_PythiaTuneUp) return TString("WH_Sig_PythiaTuneUp");

  else if (is==(int)SampleType::ZH_Sig) return TString("ZH_Sig");
  else if (is==(int)SampleType::ZH_Sig_PythiaScaleDn) return TString("ZH_Sig_PythiaScaleDn");
  else if (is==(int)SampleType::ZH_Sig_PythiaScaleUp) return TString("ZH_Sig_PythiaScaleUp");
  else if (is==(int)SampleType::ZH_Sig_PythiaTuneDn) return TString("ZH_Sig_PythiaTuneDn");
  else if (is==(int)SampleType::ZH_Sig_PythiaTuneUp) return TString("ZH_Sig_PythiaTuneUp");

  else if (is==(int)SampleType::tt_Sig) return TString("tt_Sig");
  else if (is==(int)SampleType::tt_Sig_PythiaScaleDn) return TString("tt_Sig_PythiaScaleDn");
  else if (is==(int)SampleType::tt_Sig_PythiaScaleUp) return TString("tt_Sig_PythiaScaleUp");
  else if (is==(int)SampleType::tt_Sig_PythiaTuneDn) return TString("tt_Sig_PythiaTuneDn");
  else if (is==(int)SampleType::tt_Sig_PythiaTuneUp) return TString("tt_Sig_PythiaTuneUp");

  else if (is==(int)SampleType::gg_Bkg_MCFM) return TString("gg_Bkg_MCFM");

  else if (is==(int)SampleType::VBF_Sig_Phantom) return TString("VBF_Sig_Phantom");
  else if (is==(int)SampleType::VBF_Bkg_Phantom) return TString("VBF_Bkg_Phantom");

  else if (is==(int)SampleType::qq_Bkg) return TString("qq_Bkg");

  else return TString("nSampleTypes");
}
std::vector<std::pair<const TString, const int>> sampleTypeList;
void constructSampleTypeList(){
  sampleTypeList.clear();

  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_gg_Sig_POWHEG.txt", gg_Sig));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_gg_Sig_PythiaScaleDn_POWHEG.txt", gg_Sig_PythiaScaleDn));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_gg_Sig_PythiaScaleUp_POWHEG.txt", gg_Sig_PythiaScaleUp));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_gg_Sig_PythiaTuneDn_POWHEG.txt", gg_Sig_PythiaTuneDn));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_gg_Sig_PythiaTuneUp_POWHEG.txt", gg_Sig_PythiaTuneUp));

  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_gg_Sig_MiNLO_POWHEG.txt", gg_Sig_MiNLO));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_gg_Sig_PythiaScaleDn_MiNLO_POWHEG.txt", gg_Sig_PythiaScaleDn_MiNLO));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_gg_Sig_PythiaScaleUp_MiNLO_POWHEG.txt", gg_Sig_PythiaScaleUp_MiNLO));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_gg_Sig_PythiaTuneDn_MiNLO_POWHEG.txt", gg_Sig_PythiaTuneDn_MiNLO));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_gg_Sig_PythiaTuneUp_MiNLO_POWHEG.txt", gg_Sig_PythiaTuneUp_MiNLO));

  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_VBF_Sig_POWHEG.txt", VBF_Sig));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_VBF_Sig_PythiaScaleDn_POWHEG.txt", VBF_Sig_PythiaScaleDn));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_VBF_Sig_PythiaScaleUp_POWHEG.txt", VBF_Sig_PythiaScaleUp));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_VBF_Sig_PythiaTuneDn_POWHEG.txt", VBF_Sig_PythiaTuneDn));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_VBF_Sig_PythiaTuneUp_POWHEG.txt", VBF_Sig_PythiaTuneUp));

  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_WH_Sig_POWHEG.txt", WH_Sig));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_WH_Sig_PythiaScaleDn_POWHEG.txt", WH_Sig_PythiaScaleDn));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_WH_Sig_PythiaScaleUp_POWHEG.txt", WH_Sig_PythiaScaleUp));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_WH_Sig_PythiaTuneDn_POWHEG.txt", WH_Sig_PythiaTuneDn));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_WH_Sig_PythiaTuneUp_POWHEG.txt", WH_Sig_PythiaTuneUp));

  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_ZH_Sig_POWHEG.txt", ZH_Sig));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_ZH_Sig_PythiaScaleDn_POWHEG.txt", ZH_Sig_PythiaScaleDn));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_ZH_Sig_PythiaScaleUp_POWHEG.txt", ZH_Sig_PythiaScaleUp));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_ZH_Sig_PythiaTuneDn_POWHEG.txt", ZH_Sig_PythiaTuneDn));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_ZH_Sig_PythiaTuneUp_POWHEG.txt", ZH_Sig_PythiaTuneUp));

  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_tt_Sig_POWHEG.txt", tt_Sig));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_tt_Sig_PythiaScaleDn_POWHEG.txt", tt_Sig_PythiaScaleDn));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_tt_Sig_PythiaScaleUp_POWHEG.txt", tt_Sig_PythiaScaleUp));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_tt_Sig_PythiaTuneDn_POWHEG.txt", tt_Sig_PythiaTuneDn));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_tt_Sig_PythiaTuneUp_POWHEG.txt", tt_Sig_PythiaTuneUp));

  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_gg_Bkg_MCFM.txt", gg_Bkg_MCFM));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_VBF_Sig_Phantom.txt", VBF_Sig_Phantom));
  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_VBF_Bkg_Phantom.txt", VBF_Bkg_Phantom));

  sampleTypeList.push_back(std::pair<const TString, const int>("Samples_qq_Bkg_POWHEG.txt", qq_Bkg));
}

const float luminosity[1]={ 1 };
const double MCFM_widthrescale = 0.9771908764; // Signal g1=sqrt(MCFM_widthrescale)

//Directories
const TString package_dir="/work-zfs/lhc/usarica/CMS-related/Analysis/MassWidth_2016/CMSSW_8_0_12/src/Analysis/HiggsMassWidth/";
const TString user_input_dir="/work-zfs/lhc/usarica/CMS-related/CJLSTProduction/170222/";
const TString user_output_dir="/work-zfs/lhc/usarica/CMS-related/Analysis/MassWidth_2016/Moriond2017_mainstream/";
const TString user_TemplateswithTrees_dir="/work-zfs/lhc/usarica/CMS-related/Analysis/MassWidth_2016/Moriond2017_mainstream/Templates/";

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

