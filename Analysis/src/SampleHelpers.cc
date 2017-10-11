#include "SampleHelpers.h"


using namespace std;


namespace SampleHelpers{
  shared_ptr<Mela> GlobalMELA;
}

void SampleHelpers::makeGlobalMELA(int CoM){ if (!GlobalMELA) GlobalMELA.reset(new Mela(CoM, 125, TVar::ERROR)); }

float SampleHelpers::findPoleMass(const TString samplename){
  float mass = -1;
  if (samplename=="") return mass;
  std::string strtmp = samplename.Data();
  std::size_t extpos = strtmp.find(".root");
  if (extpos!=std::string::npos) strtmp.erase(extpos, 5);
  std::vector<std::string> strsplit;
  HelperFunctions::splitOptionRecursive(strtmp, strsplit, 'H');
  if (strsplit.size()>1){
    std::string strmass = strsplit.at(1);
    strsplit.clear();
    HelperFunctions::splitOptionRecursive(strmass, strsplit, '_');
    strmass = strsplit.at(0);
    mass = std::stod(strmass);
  }
  return mass;
}
TTree* SampleHelpers::findTree(std::vector<TTree*> treeList, int evid){
  int ev = evid;
  for (unsigned int t=0; t<treeList.size(); t++){
    TTree* tree = treeList.at(t);
    int nevts = tree->GetEntries();
    if (ev<nevts) return tree;
    else ev -= nevts;
    if (ev<0) cerr << "findTree::ERROR: Could not find the event " << evid << endl;
  }
  return 0;
}
void SampleHelpers::getEntry(std::vector<TTree*> treeList, int evid){
  int ev = evid;
  for (unsigned int t=0; t<treeList.size(); t++){
    TTree* tree = treeList.at(t);
    int nevts = tree->GetEntries();
    if (ev<nevts){ tree->GetEntry(ev); break; }
    else ev -= nevts;
    if (ev<0) cerr << "getEntry::ERROR: Could not find the event " << evid << endl;
  }
}
float SampleHelpers::getEntry(std::vector<std::pair<TTree*, TH1F*>> treeList, int evid){
  int ev = evid;
  for (unsigned int t=0; t<treeList.size(); t++){
    TTree* tree = treeList.at(t).first;
    int nevts = tree->GetEntries();
    if (ev<nevts){ tree->GetEntry(ev); return float(1./treeList.at(t).second->GetBinContent(40)); }
    else ev -= nevts;
    if (ev<0) cerr << "getEntry::ERROR: Could not find the event " << evid << endl;
  }
  return 0;
}

TString SampleHelpers::getChannelName(const SampleHelpers::Channel chan){
  switch (chan){
  case k4mu:
    return "4mu";
  case k4e:
    return "4e";
  case k2e2mu:
    return "2e2mu";
  default:
    return "";
  }
}

std::vector<TString> SampleHelpers::constructSamplesList(TString strsample, float sqrts){
  std::vector<TString> samples;
  if (strsample=="VBF_Sig_POWHEG"){
    if (sqrts==13){
      samples.push_back("VBFH115");
      samples.push_back("VBFH120");
      samples.push_back("VBFH124");
      samples.push_back("VBFH125");
      samples.push_back("VBFH126");
      samples.push_back("VBFH130");
      samples.push_back("VBFH135");
      samples.push_back("VBFH140");
      samples.push_back("VBFH150");
      samples.push_back("VBFH155");
      samples.push_back("VBFH160");
      samples.push_back("VBFH165");
      samples.push_back("VBFH170");
      samples.push_back("VBFH175");
      samples.push_back("VBFH180");
      samples.push_back("VBFH190");
      samples.push_back("VBFH210");
      samples.push_back("VBFH230");
      samples.push_back("VBFH250");
      samples.push_back("VBFH270");
      samples.push_back("VBFH300");
      samples.push_back("VBFH350");
      samples.push_back("VBFH450");
      samples.push_back("VBFH500");
      samples.push_back("VBFH550");
      samples.push_back("VBFH600");
      samples.push_back("VBFH700");
      samples.push_back("VBFH750");
      samples.push_back("VBFH800");
      samples.push_back("VBFH900");
      samples.push_back("VBFH1000");
      samples.push_back("VBFH2000");
      samples.push_back("VBFH2500");
      samples.push_back("VBFH3000");
    }
  }
  else if (strsample=="WH_Sig_POWHEG"){
    if (sqrts==13){
      samples.push_back("WminusH115"); samples.push_back("WplusH115");
      samples.push_back("WminusH120"); samples.push_back("WplusH120");
      samples.push_back("WminusH124"); samples.push_back("WplusH124");
      samples.push_back("WminusH125"); samples.push_back("WplusH125");
      samples.push_back("WminusH126"); samples.push_back("WplusH126");
      samples.push_back("WminusH130"); samples.push_back("WplusH130");
      samples.push_back("WminusH135"); samples.push_back("WplusH135");
      samples.push_back("WminusH140"); samples.push_back("WplusH140");
      samples.push_back("WminusH145"); samples.push_back("WplusH145");
      samples.push_back("WminusH150"); samples.push_back("WplusH150");
      samples.push_back("WminusH155"); samples.push_back("WplusH155");
      samples.push_back("WminusH160"); samples.push_back("WplusH160");
      samples.push_back("WminusH165"); samples.push_back("WplusH165");
      samples.push_back("WminusH170"); samples.push_back("WplusH170");
      samples.push_back("WminusH175"); samples.push_back("WplusH175");
      samples.push_back("WminusH180"); samples.push_back("WplusH180");
      samples.push_back("WminusH190"); samples.push_back("WplusH190");
      samples.push_back("WminusH200"); samples.push_back("WplusH200");
      samples.push_back("WminusH210"); samples.push_back("WplusH210");
      samples.push_back("WminusH230"); samples.push_back("WplusH230");
    }
  }
  else if (strsample=="ZH_Sig_POWHEG"){
    if (sqrts==13){
      samples.push_back("ZH115");
      samples.push_back("ZH120");
      samples.push_back("ZH124");
      samples.push_back("ZH125");
      samples.push_back("ZH126");
      samples.push_back("ZH130");
      samples.push_back("ZH135");
      samples.push_back("ZH140");
      samples.push_back("ZH145");
      samples.push_back("ZH150");
      samples.push_back("ZH155");
      samples.push_back("ZH160");
      samples.push_back("ZH165");
      samples.push_back("ZH170");
      samples.push_back("ZH175");
      samples.push_back("ZH180");
      samples.push_back("ZH190");
      samples.push_back("ZH200");
      samples.push_back("ZH210");
      samples.push_back("ZH230");
    }
  }
  else if (strsample=="gg_Sig_POWHEG"){
    if (sqrts==13){
      samples.push_back("ggH115");
      samples.push_back("ggH120");
      samples.push_back("ggH124");
      samples.push_back("ggH125");
      samples.push_back("ggH126");
      samples.push_back("ggH130");
      samples.push_back("ggH135");
      samples.push_back("ggH145");
      samples.push_back("ggH150");
      samples.push_back("ggH155");
      samples.push_back("ggH160");
      samples.push_back("ggH165");
      samples.push_back("ggH170");
      samples.push_back("ggH175");
      samples.push_back("ggH180");
      samples.push_back("ggH190");
      samples.push_back("ggH200");
      samples.push_back("ggH210");
      samples.push_back("ggH230");
      samples.push_back("ggH250");
      samples.push_back("ggH270");
      samples.push_back("ggH300");
      samples.push_back("ggH350");
      samples.push_back("ggH400");
      samples.push_back("ggH450");
      samples.push_back("ggH500");
      samples.push_back("ggH550");
      samples.push_back("ggH600");
      samples.push_back("ggH700");
      samples.push_back("ggH750");
      samples.push_back("ggH800");
      samples.push_back("ggH900");
      samples.push_back("ggH1000");
      samples.push_back("ggH1500");
      samples.push_back("ggH2000");
      samples.push_back("ggH2500");
      samples.push_back("ggH3000");
    }
  }
  else if (strsample=="gg_Sig_SM_MCFM"){
    if (sqrts==13){
      samples.push_back("ggTo4mu_0PMH125_MCFM701");
      samples.push_back("ggTo4e_0PMH125_MCFM701");
      samples.push_back("ggTo2e2mu_0PMH125_MCFM701");
      samples.push_back("ggTo2e2tau_0PMH125_MCFM701");
      samples.push_back("ggTo2mu2tau_0PMH125_MCFM701");
      samples.push_back("ggTo4tau_0PMH125_MCFM701");
    }
  }
  else if (strsample=="gg_Bkg_MCFM"){
    if (sqrts==13){
      samples.push_back("ggTo4mu_Contin_MCFM701");
      samples.push_back("ggTo4e_Contin_MCFM701");
      samples.push_back("ggTo2e2mu_Contin_MCFM701");
      samples.push_back("ggTo2e2tau_Contin_MCFM701");
      samples.push_back("ggTo2mu2tau_Contin_MCFM701");
      samples.push_back("ggTo4tau_Contin_MCFM701");
    }
  }
  else if (strsample=="VV_Sig_Phantom"){
    if (sqrts==13){
      samples.push_back("VBFTo2e2muJJ_0PMH125_phantom128");
      samples.push_back("VBFTo4muJJ_0PMH125_phantom128");
      samples.push_back("VBFTo4eJJ_0PMH125_phantom128");
    }
  }
  else if (strsample=="VV_Bkg_Phantom"){
    if (sqrts==13){
      samples.push_back("VBFTo2e2muJJ_Contin_phantom128");
      samples.push_back("VBFTo4muJJ_Contin_phantom128");
      samples.push_back("VBFTo4eJJ_Contin_phantom128");
    }
  }
  else if (strsample=="qq_Bkg"){
    if (sqrts==13){
      samples.push_back("ZZTo4l");
      samples.push_back("ZZTo4l_ext");
    }
  }
  else if (strsample=="qq_Bkg_Combined"){ // ZZTo4l + ZZTo4l_ext (use this one)
    if (sqrts==13){
      samples.push_back("ZZTo4lCombined");
    }
  }
  return samples;
}
void SampleHelpers::getSamplesList(float sqrts, std::vector<TString> s, std::vector<TString>& vs){
  for (auto& ss : s){
    vector<TString> dumappend = constructSamplesList(ss, sqrts);
    HelperFunctions::appendVector<TString>(vs, dumappend);
  }
}
void SampleHelpers::getSamplePairs(float sqrts, std::vector<TString> s1, std::vector<TString> s2, std::vector<TString>& vs1, std::vector<TString>& vs2){
  getSamplesList(sqrts, s1, vs1);
  getSamplesList(sqrts, s2, vs2);
}


