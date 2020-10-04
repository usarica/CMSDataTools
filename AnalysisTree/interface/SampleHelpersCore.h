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
#include <csignal>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TList.h"
#include "HostHelpersCore.h"
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"
#include "CMSEDMWrapper.h"
#include "CMSLorentzVector.h"


namespace SampleHelpers{
  extern std::shared_ptr<Mela> GlobalMELA;
  extern TDirectory* const rootTDirectory;
  extern volatile sig_atomic_t doSignalInterrupt;

  void makeGlobalMELA(int CoM, TVar::VerbosityLevel verbosity=TVar::ERROR);

  void setSignalInterrupt(int snum);

  std::vector<TString> lsdir(TString indir, HostHelpers::Hosts const* target_host=nullptr);

  float findPoleMass(const TString samplename);
  TTree* findTree(std::vector<TTree*> const& treeList, int evid);
  bool branchExists(TTree* tree, TString strname);
  bool aliasExists(TTree* tree, TString strname);
  void getEntry(std::vector<TTree*>& treeList, int evid);
  float getEntry(std::vector<std::pair<TTree*, TH1F*>>& treeList, int evid);

  template<typename T> void bookBranch(TTree* tree, TString strname, T* var);
  template<typename T> void bookEDMBranch(TTree* tree, TString strname, T* var);
  template<typename T> void putBranch(TTree* tree, TString strname, T& var);
}

template<typename T> void SampleHelpers::bookBranch(TTree* tree, TString strname, T* var){
  if (tree){
    TString strbname = strname;
    // First search in aliases and replace strbname
    const TList* aliasList = (const TList*) tree->GetListOfAliases();
    if (aliasList){
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
    }
    if (SampleHelpers::branchExists(tree, strbname)){
      TString bnamegen="";
      if (strbname.Contains(".")){
        std::vector<TString> tmplist;
        HelperFunctions::splitOptionRecursive(strbname, tmplist, '.');
        if (!tmplist.empty()) bnamegen = tmplist.front() + "*";
      }
      if (bnamegen!="") tree->SetBranchStatus(bnamegen, 1);
      else tree->SetBranchStatus(strbname, 1);
      tree->SetBranchAddress(strbname, var);
      MELAStreamHelpers::MELAout << "SampleHelpers::bookBranch: Booked " << strbname << " to address " << var << "." << std::endl;
      if (!tree->GetBranchStatus(strbname)) MELAStreamHelpers::MELAerr << "SampleHelpers::bookBranch: Booked branch " << strbname << " has status 0!" << std::endl;
    }
    else MELAStreamHelpers::MELAout << "SampleHelpers::bookBranch: Branch " << strbname << " does not exist in tree " << tree->GetName() << "!" << std::endl;
  }
}
template<typename T> void SampleHelpers::bookEDMBranch(TTree* tree, TString strname, T* var){
  if (tree){
    TString strbname = strname;
    // First search in aliases and replace strbname
    const TList* aliasList = (const TList*) tree->GetListOfAliases();
    if (aliasList){
      for (int ia=0; ia<aliasList->GetSize(); ia++){
        TObject* alias = aliasList->At(ia);
        if (!alias) continue;
        TString aname = alias->GetName();
        if (aname==strbname){
          strbname = alias->GetTitle();
          break;
        }
      }
    }
    // Ensure that the branch name is "edm_product_name."
    if (strbname.Contains(".")){
      std::vector<TString> tmplist;
      HelperFunctions::splitOptionRecursive(strbname, tmplist, '.');
      if (!tmplist.empty()) strbname = tmplist.front() + ".";
    }
    SampleHelpers::bookBranch(tree, strbname, var);
  }
}
template<typename T> void SampleHelpers::putBranch(TTree* tree, TString strname, T& var){
  if (tree){
    // Do not check for branch alias
    if (!SampleHelpers::branchExists(tree, strname)){
      if (!tree->Branch(strname, &var)) MELAStreamHelpers::MELAerr << "SampleHelpers::putBranch: Did not succeed in creating a new branch for " << strname << " in tree " << tree->GetName() << std::endl;
    }
  }
}

template void SampleHelpers::bookBranch<bool>(TTree* tree, TString strname, bool* var);
template void SampleHelpers::bookBranch<short>(TTree* tree, TString strname, short* var);
template void SampleHelpers::bookBranch<unsigned int>(TTree* tree, TString strname, unsigned int* var);
template void SampleHelpers::bookBranch<int>(TTree* tree, TString strname, int* var);
template void SampleHelpers::bookBranch<unsigned long>(TTree* tree, TString strname, unsigned long* var);
template void SampleHelpers::bookBranch<long>(TTree* tree, TString strname, long* var);
template void SampleHelpers::bookBranch<unsigned long long>(TTree* tree, TString strname, unsigned long long* var);
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
template void SampleHelpers::bookBranch<std::vector<unsigned long long>*>(TTree* tree, TString strname, std::vector<unsigned long long>** var);
template void SampleHelpers::bookBranch<std::vector<long long>*>(TTree* tree, TString strname, std::vector<long long>** var);
template void SampleHelpers::bookBranch<std::vector<float>*>(TTree* tree, TString strname, std::vector<float>** var);
template void SampleHelpers::bookBranch<std::vector<double>*>(TTree* tree, TString strname, std::vector<double>** var);
template void SampleHelpers::bookBranch<std::vector<CMSLorentzVector>*>(TTree* tree, TString strname, std::vector<CMSLorentzVector>** var);
// CMSEDMWrapper equivalents of fundamental types
template void SampleHelpers::bookEDMBranch<CMSEDMWrapper<bool>*>(TTree* tree, TString strname, CMSEDMWrapper<bool>** var);
template void SampleHelpers::bookEDMBranch<CMSEDMWrapper<short>*>(TTree* tree, TString strname, CMSEDMWrapper<short>** var);
template void SampleHelpers::bookEDMBranch<CMSEDMWrapper<unsigned int>*>(TTree* tree, TString strname, CMSEDMWrapper<unsigned int>** var);
template void SampleHelpers::bookEDMBranch<CMSEDMWrapper<int>*>(TTree* tree, TString strname, CMSEDMWrapper<int>** var);
template void SampleHelpers::bookEDMBranch<CMSEDMWrapper<unsigned long>*>(TTree* tree, TString strname, CMSEDMWrapper<unsigned long>** var);
template void SampleHelpers::bookEDMBranch<CMSEDMWrapper<long>*>(TTree* tree, TString strname, CMSEDMWrapper<long>** var);
template void SampleHelpers::bookEDMBranch<CMSEDMWrapper<unsigned long long>*>(TTree* tree, TString strname, CMSEDMWrapper<unsigned long long>** var);
template void SampleHelpers::bookEDMBranch<CMSEDMWrapper<long long>*>(TTree* tree, TString strname, CMSEDMWrapper<long long>** var);
template void SampleHelpers::bookEDMBranch<CMSEDMWrapper<float>*>(TTree* tree, TString strname, CMSEDMWrapper<float>** var);
template void SampleHelpers::bookEDMBranch<CMSEDMWrapper<double>*>(TTree* tree, TString strname, CMSEDMWrapper<double>** var);
template void SampleHelpers::bookEDMBranch<CMSEDMWrapper<CMSLorentzVector>*>(TTree* tree, TString strname, CMSEDMWrapper<CMSLorentzVector>** var);
template void SampleHelpers::bookEDMBranch<CMSEDMWrapper<std::vector<bool>>*>(TTree* tree, TString strname, CMSEDMWrapper<std::vector<bool>>** var);
template void SampleHelpers::bookEDMBranch<CMSEDMWrapper<std::vector<short>>*>(TTree* tree, TString strname, CMSEDMWrapper<std::vector<short>>** var);
template void SampleHelpers::bookEDMBranch<CMSEDMWrapper<std::vector<unsigned int>>*>(TTree* tree, TString strname, CMSEDMWrapper<std::vector<unsigned int>>** var);
template void SampleHelpers::bookEDMBranch<CMSEDMWrapper<std::vector<int>>*>(TTree* tree, TString strname, CMSEDMWrapper<std::vector<int>>** var);
template void SampleHelpers::bookEDMBranch<CMSEDMWrapper<std::vector<unsigned long>>*>(TTree* tree, TString strname, CMSEDMWrapper<std::vector<unsigned long>>** var);
template void SampleHelpers::bookEDMBranch<CMSEDMWrapper<std::vector<long>>*>(TTree* tree, TString strname, CMSEDMWrapper<std::vector<long>>** var);
template void SampleHelpers::bookEDMBranch<CMSEDMWrapper<std::vector<unsigned long long>>*>(TTree* tree, TString strname, CMSEDMWrapper<std::vector<unsigned long long>>** var);
template void SampleHelpers::bookEDMBranch<CMSEDMWrapper<std::vector<long long>>*>(TTree* tree, TString strname, CMSEDMWrapper<std::vector<long long>>** var);
template void SampleHelpers::bookEDMBranch<CMSEDMWrapper<std::vector<float>>*>(TTree* tree, TString strname, CMSEDMWrapper<std::vector<float>>** var);
template void SampleHelpers::bookEDMBranch<CMSEDMWrapper<std::vector<double>>*>(TTree* tree, TString strname, CMSEDMWrapper<std::vector<double>>** var);
template void SampleHelpers::bookEDMBranch<CMSEDMWrapper<std::vector<CMSLorentzVector>>*>(TTree* tree, TString strname, CMSEDMWrapper<std::vector<CMSLorentzVector>>** var);

// This framework does NOT 'put' CMSEDMWrappers: It is not a practical analysis tool!
template void SampleHelpers::putBranch<bool>(TTree* tree, TString strname, bool& var);
template void SampleHelpers::putBranch<short>(TTree* tree, TString strname, short& var);
template void SampleHelpers::putBranch<unsigned int>(TTree* tree, TString strname, unsigned int& var);
template void SampleHelpers::putBranch<int>(TTree* tree, TString strname, int& var);
template void SampleHelpers::putBranch<unsigned long>(TTree* tree, TString strname, unsigned long& var);
template void SampleHelpers::putBranch<long>(TTree* tree, TString strname, long& var);
template void SampleHelpers::putBranch<unsigned long long>(TTree* tree, TString strname, unsigned long long& var);
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
template void SampleHelpers::putBranch<std::vector<unsigned long long>>(TTree* tree, TString strname, std::vector<unsigned long long>& var);
template void SampleHelpers::putBranch<std::vector<long long>>(TTree* tree, TString strname, std::vector<long long>& var);
template void SampleHelpers::putBranch<std::vector<float>>(TTree* tree, TString strname, std::vector<float>& var);
template void SampleHelpers::putBranch<std::vector<double>>(TTree* tree, TString strname, std::vector<double>& var);
template void SampleHelpers::putBranch<std::vector<CMSLorentzVector>>(TTree* tree, TString strname, std::vector<CMSLorentzVector>& var);


#endif
