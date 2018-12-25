#ifndef SAMPLEHELPERSCORE_H
#define SAMPLEHELPERSCORE_H

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
#include "MELAStreamHelpers.hh"
#include "CMSLorentzVector.h"


namespace SampleHelpers{
  extern std::shared_ptr<Mela> GlobalMELA;

  void makeGlobalMELA(int CoM, TVar::VerbosityLevel verbosity=TVar::ERROR);

  std::vector<TString> lsdir(TString indir);

  float findPoleMass(const TString samplename);
  TTree* findTree(std::vector<TTree*> const& treeList, int evid);
  bool branchExists(TTree* tree, TString strname);
  void getEntry(std::vector<TTree*>& treeList, int evid);
  float getEntry(std::vector<std::pair<TTree*, TH1F*>>& treeList, int evid);

  template<typename T> void bookBranch(TTree* tree, TString strname, T* var);
  template<typename T> void putBranch(TTree* tree, TString strname, T& var);
}

template<typename T> void SampleHelpers::bookBranch(TTree* tree, TString strname, T* var){
  if (tree){
    TString strbname = strname;
    // First search in aliases and replace strbname
    const TList* aliasList = (const TList*) tree->GetListOfAliases();
    for (int ia=0; ia<aliasList->GetSize(); ia++){
      TObject* alias = aliasList->At(ia);
      if (!alias) continue;
      TString aname = alias->GetName();
      if (aname==strbname){
        strbname = alias->GetTitle();
        MELAStreamHelpers::MELAout << "SampleHelpers::bookBranch: Using branch name " << strbname << " instead of the alias " << strname << "." << std::endl;
        break;
      }
    }
    if (SampleHelpers::branchExists(tree, strbname)){
      tree->SetBranchStatus(strbname, 1);
      tree->SetBranchAddress(strbname, var);
    }
    else MELAStreamHelpers::MELAout << "SampleHelpers::bookBranch: Branch " << strbname << " does not exist in tree " << tree->GetName() << "!" << std::endl;
  }
}
template<typename T> void SampleHelpers::putBranch(TTree* tree, TString strname, T& var){
  if (tree){
    // Do not check for branch alias
    if (!SampleHelpers::branchExists(tree, strname)) tree->Branch(strname, &var);
    else{
      tree->SetBranchStatus(strname, 1);
      tree->SetBranchAddress(strname, &var);
    }
  }
}

template void SampleHelpers::bookBranch<bool>(TTree* tree, TString strname, bool* var);
template void SampleHelpers::bookBranch<short>(TTree* tree, TString strname, short* var);
template void SampleHelpers::bookBranch<unsigned int>(TTree* tree, TString strname, unsigned int* var);
template void SampleHelpers::bookBranch<int>(TTree* tree, TString strname, int* var);
template void SampleHelpers::bookBranch<unsigned long>(TTree* tree, TString strname, unsigned long* var);
template void SampleHelpers::bookBranch<long>(TTree* tree, TString strname, long* var);
template void SampleHelpers::bookBranch<long long>(TTree* tree, TString strname, long long* var);
template void SampleHelpers::bookBranch<float>(TTree* tree, TString strname, float* var);
template void SampleHelpers::bookBranch<double>(TTree* tree, TString strname, double* var);
template void SampleHelpers::bookBranch<CMSLorentzVector>(TTree* tree, TString strname, CMSLorentzVector* var);
template void SampleHelpers::bookBranch<std::vector<bool>*>(TTree* tree, TString strname, std::vector<bool>** var);
template void SampleHelpers::bookBranch<std::vector<short>*>(TTree* tree, TString strname, std::vector<short>** var);
template void SampleHelpers::bookBranch<std::vector<unsigned int>*>(TTree* tree, TString strname, std::vector<unsigned int>** var);
template void SampleHelpers::bookBranch<std::vector<int>*>(TTree* tree, TString strname, std::vector<int>** var);
template void SampleHelpers::bookBranch<std::vector<unsigned long>*>(TTree* tree, TString strname, std::vector<unsigned long>** var);
template void SampleHelpers::bookBranch<std::vector<long>*>(TTree* tree, TString strname, std::vector<long>** var);
template void SampleHelpers::bookBranch<std::vector<long long>*>(TTree* tree, TString strname, std::vector<long long>** var);
template void SampleHelpers::bookBranch<std::vector<float>*>(TTree* tree, TString strname, std::vector<float>** var);
template void SampleHelpers::bookBranch<std::vector<double>*>(TTree* tree, TString strname, std::vector<double>** var);
template void SampleHelpers::bookBranch<std::vector<CMSLorentzVector>*>(TTree* tree, TString strname, std::vector<CMSLorentzVector>** var);

template void SampleHelpers::putBranch<bool>(TTree* tree, TString strname, bool& var);
template void SampleHelpers::putBranch<short>(TTree* tree, TString strname, short& var);
template void SampleHelpers::putBranch<unsigned int>(TTree* tree, TString strname, unsigned int& var);
template void SampleHelpers::putBranch<int>(TTree* tree, TString strname, int& var);
template void SampleHelpers::putBranch<unsigned long>(TTree* tree, TString strname, unsigned long& var);
template void SampleHelpers::putBranch<long>(TTree* tree, TString strname, long& var);
template void SampleHelpers::putBranch<long long>(TTree* tree, TString strname, long long& var);
template void SampleHelpers::putBranch<float>(TTree* tree, TString strname, float& var);
template void SampleHelpers::putBranch<double>(TTree* tree, TString strname, double& var);
template void SampleHelpers::putBranch<CMSLorentzVector>(TTree* tree, TString strname, CMSLorentzVector& var);
template void SampleHelpers::putBranch<std::vector<bool>>(TTree* tree, TString strname, std::vector<bool>& var);
template void SampleHelpers::putBranch<std::vector<short>>(TTree* tree, TString strname, std::vector<short>& var);
template void SampleHelpers::putBranch<std::vector<unsigned int>>(TTree* tree, TString strname, std::vector<unsigned int>& var);
template void SampleHelpers::putBranch<std::vector<int>>(TTree* tree, TString strname, std::vector<int>& var);
template void SampleHelpers::putBranch<std::vector<unsigned long>>(TTree* tree, TString strname, std::vector<unsigned long>& var);
template void SampleHelpers::putBranch<std::vector<long>>(TTree* tree, TString strname, std::vector<long>& var);
template void SampleHelpers::putBranch<std::vector<long long>>(TTree* tree, TString strname, std::vector<long long>& var);
template void SampleHelpers::putBranch<std::vector<float>>(TTree* tree, TString strname, std::vector<float>& var);
template void SampleHelpers::putBranch<std::vector<double>>(TTree* tree, TString strname, std::vector<double>& var);
template void SampleHelpers::putBranch<std::vector<CMSLorentzVector>>(TTree* tree, TString strname, std::vector<CMSLorentzVector>& var);


#endif
