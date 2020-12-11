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
      //MELAStreamHelpers::MELAout << "SampleHelpers::putBranch: Branching " << strname << std::endl;
      if (!tree->Branch(strname, &var)) MELAStreamHelpers::MELAerr << "SampleHelpers::putBranch: Did not succeed in creating a new branch for " << strname << " in tree " << tree->GetName() << std::endl;
    }
  }
}


#endif
