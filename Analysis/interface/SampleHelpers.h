#ifndef SAMPLEHELPERS_H
#define SAMPLEHELPERS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <memory>
#include <utility>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <unordered_map>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TList.h"
#include "HelperFunctions.h"


namespace SampleHelpers{
  enum Channel{
    k4mu,
    k4e,
    k2e2mu,
    k4l,
    k2l2l,
    NChannels
  };

  extern shared_ptr<Mela> GlobalMELA;

  void makeGlobalMELA(int CoM);

  float findPoleMass(const TString samplename);
  TTree* findTree(std::vector<TTree*> const& treeList, int evid);
  void getEntry(std::vector<TTree*>& treeList, int evid);
  float getEntry(std::vector<std::pair<TTree*, TH1F*>>& treeList, int evid);

  TString getChannelName(const SampleHelpers::Channel chan);
  SampleHelpers::Channel getChannelFromName(const TString channame);
  bool testChannel(SampleHelpers::Channel const& targetChannel, short const& Z1Flav, short const& Z2Flav);

  std::vector<TString> constructSamplesList(TString strsample, float sqrts);
  void getSamplesList(float sqrts, std::vector<TString> const& s, std::vector<TString>& vs);
  void getSamplePairs(float sqrts, std::vector<TString> const& s1, std::vector<TString> const& s2, std::vector<TString>& vs1, std::vector<TString>& vs2);

  template <typename T> std::vector<std::pair<T*, T*>> getZXFR_SS();

  template<typename T> void bookBranch(TTree* tree, TString strname, T* var);
  template<typename T> void putBranch(TTree* tree, TString strname, T& var);
}

// Pairs are (e-EB, e-EE), (mu-EB, mu-EE)
template <typename T> std::vector<std::pair<T*, T*>> SampleHelpers::getZXFR_SS(){
  TString hname[2][2]={
    { "FR_SS_electron_EB", "FR_SS_electron_EE" },
    { "FR_SS_muon_EB", "FR_SS_muon_EE" }
  };
  TFile* finput;
  finput = TFile::Open("../data/FakeRate_SS_Moriond368.root", "read");
  std::vector<std::pair<T*, T*>> result;
  for (int f=0; f<2; f++){
    if (!(finput!=0 && finput->IsOpen())){ std::cerr << "getZXFR_SS: File is not open!" << std::endl; return result; }
    else std::cout << "getZXFR_SS: File opened" << std::endl;
    T* htmp[2];
    gROOT->cd();
    for (unsigned int t=0; t<2; t++) htmp[t] = (T*) finput->Get(hname[f][t]);
    result.push_back(std::pair<T*, T*>((T*) htmp[0]->Clone(Form("%s_copy", hname[f][0].Data())), (T*) htmp[1]->Clone(Form("%s_copy", hname[f][1].Data()))));
  }
  if (finput!=0 && finput->IsOpen()) finput->Close();
  return result;
}


template<typename T> void SampleHelpers::bookBranch(TTree* tree, TString strname, T* var){
  if (tree!=nullptr){
    bool found=false;
    const TList* blist = (const TList*)tree->GetListOfBranches();
    for (int ib=0; ib<blist->GetSize(); ib++){
      TString bname = blist->At(ib)->GetName();
      if (strname==bname){ found=true; break; }
    }
    if (found){
      tree->SetBranchStatus(strname, 1);
      tree->SetBranchAddress(strname, var);
    }
  }
}
template<typename T> void SampleHelpers::putBranch(TTree* tree, TString strname, T& var){
  if (tree!=nullptr) tree->Branch(strname, &var);
}

template std::vector<std::pair<TH1F*, TH1F*>> SampleHelpers::getZXFR_SS<TH1F>();
template std::vector<std::pair<TGraphErrors*, TGraphErrors*>> SampleHelpers::getZXFR_SS<TGraphErrors>();


template void SampleHelpers::bookBranch<short>(TTree* tree, TString strname, short* var);
template void SampleHelpers::bookBranch<unsigned int>(TTree* tree, TString strname, unsigned int* var);
template void SampleHelpers::bookBranch<int>(TTree* tree, TString strname, int* var);
template void SampleHelpers::bookBranch<long>(TTree* tree, TString strname, long* var);
template void SampleHelpers::bookBranch<float>(TTree* tree, TString strname, float* var);
template void SampleHelpers::bookBranch<double>(TTree* tree, TString strname, double* var);
template void SampleHelpers::bookBranch<std::vector<unsigned int>*>(TTree* tree, TString strname, std::vector<unsigned int>** var);
template void SampleHelpers::bookBranch<std::vector<int>*>(TTree* tree, TString strname, std::vector<int>** var);
template void SampleHelpers::bookBranch<std::vector<long>*>(TTree* tree, TString strname, std::vector<long>** var);
template void SampleHelpers::bookBranch<std::vector<float>*>(TTree* tree, TString strname, std::vector<float>** var);
template void SampleHelpers::bookBranch<std::vector<double>*>(TTree* tree, TString strname, std::vector<double>** var);

template void SampleHelpers::putBranch<short>(TTree* tree, TString strname, short& var);
template void SampleHelpers::putBranch<unsigned int>(TTree* tree, TString strname, unsigned int& var);
template void SampleHelpers::putBranch<int>(TTree* tree, TString strname, int& var);
template void SampleHelpers::putBranch<long>(TTree* tree, TString strname, long& var);
template void SampleHelpers::putBranch<float>(TTree* tree, TString strname, float& var);
template void SampleHelpers::putBranch<double>(TTree* tree, TString strname, double& var);
template void SampleHelpers::putBranch<std::vector<short>>(TTree* tree, TString strname, std::vector<short>& var);
template void SampleHelpers::putBranch<std::vector<unsigned int>>(TTree* tree, TString strname, std::vector<unsigned int>& var);
template void SampleHelpers::putBranch<std::vector<int>>(TTree* tree, TString strname, std::vector<int>& var);
template void SampleHelpers::putBranch<std::vector<long>>(TTree* tree, TString strname, std::vector<long>& var);
template void SampleHelpers::putBranch<std::vector<float>>(TTree* tree, TString strname, std::vector<float>& var);
template void SampleHelpers::putBranch<std::vector<double>>(TTree* tree, TString strname, std::vector<double>& var);


#endif
