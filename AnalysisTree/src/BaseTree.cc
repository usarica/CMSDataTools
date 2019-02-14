#include "TSystem.h"
#include "TDirectory.h"
#include "BaseTree.h"
#include "BaseTree.hpp"
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;


BaseTree::BaseTree() :
  sampleIdentifier(""),
  finput(nullptr),
  tree(nullptr),
  failedtree(nullptr),
  hCounters(nullptr),
  valid(false),
  receiver(true),
  currentEvent(-1),
  currentTree(nullptr)
{}
BaseTree::BaseTree(const TString cinput, const TString treename, const TString failedtreename, const TString countersname) :
  sampleIdentifier(""), // Sample identifier is supposed to be overwritten by the daughter class
  finput(nullptr),
  tree(nullptr),
  failedtree(nullptr),
  hCounters(nullptr),
  valid(false),
  receiver(true),
  currentEvent(-1),
  currentTree(nullptr)
{
  TDirectory* curdir = gDirectory; // Save current directory to return back to it later
  if (!gSystem->AccessPathName(cinput)){
    finput = TFile::Open(cinput, "read");
    if (finput){
      if (finput->IsOpen() && !finput->IsZombie()){
        finput->cd();
        tree = (TTree*) finput->Get(treename);
        valid = (tree!=nullptr);
        if (!valid){ finput->Close(); finput=nullptr; }
        else{
          if (failedtreename!=""){
            failedtree = (TTree*) finput->Get(failedtreename);
            if (!failedtree) cout << "BaseTree::BaseTree(" << cinput << ") does not contain " << failedtreename << endl;
          }
          if (countersname!="") hCounters = (TH1F*) finput->Get(countersname);
        }
      }
      else if (finput->IsOpen()){ finput->Close(); finput=nullptr; }
    }
  }
  curdir->cd(); // Return back to the directory before opening the input file
}
BaseTree::BaseTree(const TString treename) :
  sampleIdentifier(""),
  finput(nullptr),
  tree(new TTree(treename, "")),
  failedtree(nullptr),
  hCounters(nullptr),
  valid(true),
  receiver(false),
  currentEvent(-1),
  currentTree(nullptr)
{}

BaseTree::~BaseTree(){
  HelperFunctions::cleanUnorderedMap(valTBitss);
  HelperFunctions::cleanUnorderedMap(valbools);
  HelperFunctions::cleanUnorderedMap(valshorts);
  HelperFunctions::cleanUnorderedMap(valuints);
  HelperFunctions::cleanUnorderedMap(valints);
  HelperFunctions::cleanUnorderedMap(valulongs);
  HelperFunctions::cleanUnorderedMap(vallongs);
  HelperFunctions::cleanUnorderedMap(valulonglongs);
  HelperFunctions::cleanUnorderedMap(vallonglongs);
  HelperFunctions::cleanUnorderedMap(valfloats);
  HelperFunctions::cleanUnorderedMap(valdoubles);
  HelperFunctions::cleanUnorderedMap(valstrings);
  HelperFunctions::cleanUnorderedMap(valTStrings);
  HelperFunctions::cleanUnorderedMap(valCMSLorentzVectors);
  if (!receiver){
    HelperFunctions::cleanUnorderedMap(valVbools);
    HelperFunctions::cleanUnorderedMap(valVshorts);
    HelperFunctions::cleanUnorderedMap(valVuints);
    HelperFunctions::cleanUnorderedMap(valVints);
    HelperFunctions::cleanUnorderedMap(valVulongs);
    HelperFunctions::cleanUnorderedMap(valVlongs);
    HelperFunctions::cleanUnorderedMap(valVulonglongs);
    HelperFunctions::cleanUnorderedMap(valVlonglongs);
    HelperFunctions::cleanUnorderedMap(valVfloats);
    HelperFunctions::cleanUnorderedMap(valVdoubles);
    HelperFunctions::cleanUnorderedMap(valVstrings);
    HelperFunctions::cleanUnorderedMap(valVTStrings);
    HelperFunctions::cleanUnorderedMap(valVCMSLorentzVectors);

    HelperFunctions::cleanUnorderedMap(valVVbools);
    HelperFunctions::cleanUnorderedMap(valVVshorts);
    HelperFunctions::cleanUnorderedMap(valVVuints);
    HelperFunctions::cleanUnorderedMap(valVVints);
    HelperFunctions::cleanUnorderedMap(valVVulongs);
    HelperFunctions::cleanUnorderedMap(valVVlongs);
    HelperFunctions::cleanUnorderedMap(valVVulonglongs);
    HelperFunctions::cleanUnorderedMap(valVVlonglongs);
    HelperFunctions::cleanUnorderedMap(valVVfloats);
    HelperFunctions::cleanUnorderedMap(valVVdoubles);

    delete hCounters;
    delete failedtree;
    delete tree;
  }
}

BaseTree::BranchType BaseTree::searchBranchType(TString branchname) const{
  if (valTBitss.find(branchname)!=valTBitss.cend()) return BranchType_TBits_t;
  else if (valbools.find(branchname)!=valbools.cend()) return BranchType_bool_t;
  else if (valshorts.find(branchname)!=valshorts.cend()) return BranchType_short_t;
  else if (valuints.find(branchname)!=valuints.cend()) return BranchType_uint_t;
  else if (valints.find(branchname)!=valints.cend()) return BranchType_int_t;
  else if (valulongs.find(branchname)!=valulongs.cend()) return BranchType_ulong_t;
  else if (vallongs.find(branchname)!=vallongs.cend()) return BranchType_long_t;
  else if (valulonglongs.find(branchname)!=valulonglongs.cend()) return BranchType_ulonglong_t;
  else if (vallonglongs.find(branchname)!=vallonglongs.cend()) return BranchType_longlong_t;
  else if (valfloats.find(branchname)!=valfloats.cend()) return BranchType_float_t;
  else if (valdoubles.find(branchname)!=valdoubles.cend()) return BranchType_double_t;
  else if (valstrings.find(branchname)!=valstrings.cend()) return BranchType_string_t;
  else if (valTStrings.find(branchname)!=valTStrings.cend()) return BranchType_TString_t;
  else if (valCMSLorentzVectors.find(branchname)!=valCMSLorentzVectors.cend()) return BranchType_CMSLorentzVector_t;

  else if (valVbools.find(branchname)!=valVbools.cend()) return BranchType_vbool_t;
  else if (valVshorts.find(branchname)!=valVshorts.cend()) return BranchType_vshort_t;
  else if (valVuints.find(branchname)!=valVuints.cend()) return BranchType_vuint_t;
  else if (valVints.find(branchname)!=valVints.cend()) return BranchType_vint_t;
  else if (valVulongs.find(branchname)!=valVulongs.cend()) return BranchType_vulong_t;
  else if (valVlongs.find(branchname)!=valVlongs.cend()) return BranchType_vlong_t;
  else if (valVulonglongs.find(branchname)!=valVulonglongs.cend()) return BranchType_vulonglong_t;
  else if (valVlonglongs.find(branchname)!=valVlonglongs.cend()) return BranchType_vlonglong_t;
  else if (valVfloats.find(branchname)!=valVfloats.cend()) return BranchType_vfloat_t;
  else if (valVdoubles.find(branchname)!=valVdoubles.cend()) return BranchType_vdouble_t;
  else if (valVstrings.find(branchname)!=valVstrings.cend()) return BranchType_vstring_t;
  else if (valVTStrings.find(branchname)!=valVTStrings.cend()) return BranchType_vTString_t;
  else if (valVCMSLorentzVectors.find(branchname)!=valVCMSLorentzVectors.cend()) return BranchType_vCMSLorentzVector_t;

  else if (valVVbools.find(branchname)!=valVVbools.cend()) return BranchType_vvbool_t;
  else if (valVVshorts.find(branchname)!=valVVshorts.cend()) return BranchType_vvshort_t;
  else if (valVVuints.find(branchname)!=valVVuints.cend()) return BranchType_vvuint_t;
  else if (valVVints.find(branchname)!=valVVints.cend()) return BranchType_vvint_t;
  else if (valVVulongs.find(branchname)!=valVVulongs.cend()) return BranchType_vvulong_t;
  else if (valVVlongs.find(branchname)!=valVVlongs.cend()) return BranchType_vvlong_t;
  else if (valVVulonglongs.find(branchname)!=valVVulonglongs.cend()) return BranchType_vvulonglong_t;
  else if (valVVlonglongs.find(branchname)!=valVVlonglongs.cend()) return BranchType_vvlonglong_t;
  else if (valVVfloats.find(branchname)!=valVVfloats.cend()) return BranchType_vvfloat_t;
  else if (valVVdoubles.find(branchname)!=valVVdoubles.cend()) return BranchType_vvdouble_t;

  else return BranchType_unknown_t;
}

TFile* BaseTree::getInputFile(){ return finput; }
TTree* BaseTree::getSelectedTree(){ return tree; }
TTree* BaseTree::getFailedTree(){ return failedtree; }
TFile const* BaseTree::getInputFile() const{ return finput; }
TTree const* BaseTree::getSelectedTree() const{ return tree; }
TTree const* BaseTree::getFailedTree() const{ return failedtree; }

bool BaseTree::getSelectedEvent(int ev){
  this->resetBranches();
  bool result=false;
  if (tree && ev<tree->GetEntries()) result = (tree->GetEntry(ev)>0);
  if (result){
    currentEvent = ev;
    currentTree = tree;
  }
  return result;
}
bool BaseTree::getFailedEvent(int ev){
  this->resetBranches();
  bool result=false;
  if (failedtree && ev<failedtree->GetEntries()) result = (failedtree->GetEntry(ev)>0);
  if (result){
    currentEvent = ev;
    currentTree = failedtree;
  }
  return result;
}
bool BaseTree::getEvent(int ev){
  if (ev<this->getSelectedNEvents()) return this->getSelectedEvent(ev);
  else return this->getFailedEvent(ev-this->getSelectedNEvents());
}
void BaseTree::refreshCurrentEvent(){
  TTree* tmpTree = currentTree;
  int tmpEv = currentEvent;
  this->resetBranches();
  if (tmpTree){
    tmpTree->GetEntry(tmpEv);
    currentEvent = tmpEv;
    currentTree = tmpTree;
  }
}

int BaseTree::getSelectedNEvents() const{ return (tree ? tree->GetEntries() : 0); }
int BaseTree::getFailedNEvents() const{ return (failedtree ? failedtree->GetEntries() : 0); }
int BaseTree::getNEvents() const{ return (this->getSelectedNEvents() + this->getFailedNEvents()); }

bool BaseTree::isValid() const{ return valid; }
bool BaseTree::branchExists(TString branchname, BranchType* type){
  BranchType theType = searchBranchType(branchname);
  if (type) *type=theType;
  return (theType!=BranchType_unknown_t);
}

void BaseTree::print() const{
  for (auto const& it:valTBitss){ if (it.second){ MELAout << "\t- " << it.first << " value: "; for (unsigned int ibit=0; ibit<it.second->first.GetNbits(); ibit++){ MELAout << it.second->first.TestBitNumber(ibit); } MELAout << " (address: " << &(it.second->first) << ")" << endl; } }
  for (auto const& it:valbools){ if (it.second){ MELAout << "\t- " << it.first << " value: " << it.second->first << " (address: " << &(it.second->first) << ")" << endl; } }
  for (auto const& it:valshorts){ if (it.second){ MELAout << "\t- " << it.first << " value: " << it.second->first << " (address: " << &(it.second->first) << ")" << endl; } }
  for (auto const& it:valuints){ if (it.second){ MELAout << "\t- " << it.first << " value: " << it.second->first << " (address: " << &(it.second->first) << ")" << endl; } }
  for (auto const& it:valints){ if (it.second){ MELAout << "\t- " << it.first << " value: " << it.second->first << " (address: " << &(it.second->first) << ")" << endl; } }
  for (auto const& it:valulongs){ if (it.second){ MELAout << "\t- " << it.first << " value: " << it.second->first << " (address: " << &(it.second->first) << ")" << endl; } }
  for (auto const& it:vallongs){ if (it.second){ MELAout << "\t- " << it.first << " value: " << it.second->first << " (address: " << &(it.second->first) << ")" << endl; } }
  for (auto const& it:valulonglongs){ if (it.second){ MELAout << "\t- " << it.first << " value: " << it.second->first << " (address: " << &(it.second->first) << ")" << endl; } }
  for (auto const& it:vallonglongs){ if (it.second){ MELAout << "\t- " << it.first << " value: " << it.second->first << " (address: " << &(it.second->first) << ")" << endl; } }
  for (auto const& it:valfloats){ if (it.second){ MELAout << "\t- " << it.first << " value: " << it.second->first << " (address: " << &(it.second->first) << ")" << endl; } }
  for (auto const& it:valdoubles){ if (it.second){ MELAout << "\t- " << it.first << " value: " << it.second->first << " (address: " << &(it.second->first) << ")" << endl; } }
  for (auto const& it:valstrings){ if (it.second){ MELAout << "\t- " << it.first << " value: " << it.second->first << " (address: " << &(it.second->first) << ")" << endl; } }
  for (auto const& it:valTStrings){ if (it.second){ MELAout << "\t- " << it.first << " value: " << it.second->first << " (address: " << &(it.second->first) << ")" << endl; } }
  for (auto const& it:valCMSLorentzVectors){ if (it.second){ MELAout << "\t- " << it.first << " value: " << it.second->first << " (address: " << &(it.second->first) << ")" << endl; } }

  for (auto const& it:valVbools){ MELAout << "\t- " << it.first << " value: "; if (it.second) MELAout << *(it.second); else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
  for (auto const& it:valVshorts){ MELAout << "\t- " << it.first << " value: "; if (it.second) MELAout << *(it.second); else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
  for (auto const& it:valVuints){ MELAout << "\t- " << it.first << " value: "; if (it.second) MELAout << *(it.second); else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
  for (auto const& it:valVints){ MELAout << "\t- " << it.first << " value: "; if (it.second) MELAout << *(it.second); else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
  for (auto const& it:valVulongs){ MELAout << "\t- " << it.first << " value: "; if (it.second) MELAout << *(it.second); else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
  for (auto const& it:valVlongs){ MELAout << "\t- " << it.first << " value: "; if (it.second) MELAout << *(it.second); else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
  for (auto const& it:valVulonglongs){ MELAout << "\t- " << it.first << " value: "; if (it.second) MELAout << *(it.second); else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
  for (auto const& it:valVlonglongs){ MELAout << "\t- " << it.first << " value: "; if (it.second) MELAout << *(it.second); else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
  for (auto const& it:valVfloats){ MELAout << "\t- " << it.first << " value: "; if (it.second) MELAout << *(it.second); else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
  for (auto const& it:valVdoubles){ MELAout << "\t- " << it.first << " value: "; if (it.second) MELAout << *(it.second); else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
  for (auto const& it:valVstrings){ MELAout << "\t- " << it.first << " value: "; if (it.second) MELAout << *(it.second); else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
  for (auto const& it:valVTStrings){ MELAout << "\t- " << it.first << " value: "; if (it.second) MELAout << *(it.second); else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
  for (auto const& it:valVCMSLorentzVectors){ MELAout << "\t- " << it.first << " value: "; if (it.second) MELAout << *(it.second); else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }

  for (auto const& it:valVVbools){ MELAout << "\t- " << it.first << " value: "; if (it.second){ for (auto const& v:*(it.second)){ MELAout << "{ " << v << " }"; } } else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
  for (auto const& it:valVVshorts){ MELAout << "\t- " << it.first << " value: "; if (it.second){ for (auto const& v:*(it.second)){ MELAout << "{ " << v << " }"; } } else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
  for (auto const& it:valVVuints){ MELAout << "\t- " << it.first << " value: "; if (it.second){ for (auto const& v:*(it.second)){ MELAout << "{ " << v << " }"; } } else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
  for (auto const& it:valVVints){ MELAout << "\t- " << it.first << " value: "; if (it.second){ for (auto const& v:*(it.second)){ MELAout << "{ " << v << " }"; } } else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
  for (auto const& it:valVVulongs){ MELAout << "\t- " << it.first << " value: "; if (it.second){ for (auto const& v:*(it.second)){ MELAout << "{ " << v << " }"; } } else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
  for (auto const& it:valVVlongs){ MELAout << "\t- " << it.first << " value: "; if (it.second){ for (auto const& v:*(it.second)){ MELAout << "{ " << v << " }"; } } else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
  for (auto const& it:valVVulonglongs){ MELAout << "\t- " << it.first << " value: "; if (it.second){ for (auto const& v:*(it.second)){ MELAout << "{ " << v << " }"; } } else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
  for (auto const& it:valVVlonglongs){ MELAout << "\t- " << it.first << " value: "; if (it.second){ for (auto const& v:*(it.second)){ MELAout << "{ " << v << " }"; } } else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
  for (auto const& it:valVVfloats){ MELAout << "\t- " << it.first << " value: "; if (it.second){ for (auto const& v:*(it.second)){ MELAout << "{ " << v << " }"; } } else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
  for (auto const& it:valVVdoubles){ MELAout << "\t- " << it.first << " value: "; if (it.second){ for (auto const& v:*(it.second)){ MELAout << "{ " << v << " }"; } } else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
}

void BaseTree::resetBranches(){
  currentEvent = -1;
  currentTree = nullptr;

  this->resetBranch<BaseTree::BranchType_TBits_t>();
  this->resetBranch<BaseTree::BranchType_bool_t>();
  this->resetBranch<BaseTree::BranchType_short_t>();
  this->resetBranch<BaseTree::BranchType_uint_t>();
  this->resetBranch<BaseTree::BranchType_int_t>();
  this->resetBranch<BaseTree::BranchType_ulong_t>();
  this->resetBranch<BaseTree::BranchType_long_t>();
  this->resetBranch<BaseTree::BranchType_ulonglong_t>();
  this->resetBranch<BaseTree::BranchType_longlong_t>();
  this->resetBranch<BaseTree::BranchType_float_t>();
  this->resetBranch<BaseTree::BranchType_double_t>();
  this->resetBranch<BaseTree::BranchType_string_t>();
  this->resetBranch<BaseTree::BranchType_TString_t>();
  this->resetBranch<BaseTree::BranchType_CMSLorentzVector_t>();
  if (!receiver){
    this->resetBranch<BaseTree::BranchType_vbool_t>();
    this->resetBranch<BaseTree::BranchType_vshort_t>();
    this->resetBranch<BaseTree::BranchType_vuint_t>();
    this->resetBranch<BaseTree::BranchType_vint_t>();
    this->resetBranch<BaseTree::BranchType_vulong_t>();
    this->resetBranch<BaseTree::BranchType_vlong_t>();
    this->resetBranch<BaseTree::BranchType_vulonglong_t>();
    this->resetBranch<BaseTree::BranchType_vlonglong_t>();
    this->resetBranch<BaseTree::BranchType_vfloat_t>();
    this->resetBranch<BaseTree::BranchType_vdouble_t>();
    this->resetBranch<BaseTree::BranchType_vstring_t>();
    this->resetBranch<BaseTree::BranchType_vTString_t>();
    this->resetBranch<BaseTree::BranchType_vCMSLorentzVector_t>();

    this->resetBranch<BaseTree::BranchType_vvbool_t>();
    this->resetBranch<BaseTree::BranchType_vvshort_t>();
    this->resetBranch<BaseTree::BranchType_vvuint_t>();
    this->resetBranch<BaseTree::BranchType_vvint_t>();
    this->resetBranch<BaseTree::BranchType_vvulong_t>();
    this->resetBranch<BaseTree::BranchType_vvlong_t>();
    this->resetBranch<BaseTree::BranchType_vvulonglong_t>();
    this->resetBranch<BaseTree::BranchType_vvlonglong_t>();
    this->resetBranch<BaseTree::BranchType_vvfloat_t>();
    this->resetBranch<BaseTree::BranchType_vvdouble_t>();
  }
}

void BaseTree::getValidBranchNamesWithoutAlias(TTree* t, std::vector<TString>& res) const{
  if (!t) return;

  const TList* alist = (const TList*) t->GetListOfAliases();
  const TList* llist = (const TList*) t->GetListOfLeaves();
  const TList* blist = (const TList*) t->GetListOfBranches();
  // First check all aliases and record the proper names
  for (int ib=0; ib<alist->GetSize(); ib++){
    TString bname = alist->At(ib)->GetName();
    TString bnameproper = t->GetAlias(bname);
    TString bnamegen="";
    if (bnameproper.Contains(".")){
      std::vector<TString> tmplist;
      splitOptionRecursive(bnameproper, tmplist, '.');
      if (!tmplist.empty()) bnamegen = tmplist.front() + "*";
    }
    if (searchBranchType(bname)!=BranchType_unknown_t){
      if (bnamegen!="" && !checkListVariable(res, bnamegen)) res.push_back(bnamegen);
      else if (!checkListVariable(res, bnameproper)) res.push_back(bnameproper);
    }
  }
  // Then check all leaves
  for (int ib=0; ib<llist->GetSize(); ib++){
    TString bname = llist->At(ib)->GetName();
    TString bnamegen="";
    if (bname.Contains(".")){
      std::vector<TString> tmplist;
      splitOptionRecursive(bname, tmplist, '.');
      if (!tmplist.empty()) bnamegen = tmplist.front() + "*";
    }
    if (searchBranchType(bname)!=BranchType_unknown_t){
      if (bnamegen!="" && !checkListVariable(res, bnamegen)) res.push_back(bnamegen);
      else if (!checkListVariable(res, bname)) res.push_back(bname);
    }
  }
  // Then check all branches
  for (int ib=0; ib<blist->GetSize(); ib++){
    TString bname = blist->At(ib)->GetName();
    TString bnamegen="";
    if (bname.Contains(".")){
      std::vector<TString> tmplist;
      splitOptionRecursive(bname, tmplist, '.');
      if (!tmplist.empty()) bnamegen = tmplist.front() + "*";
    }
    if (searchBranchType(bname)!=BranchType_unknown_t){
      if (bnamegen!="" && !checkListVariable(res, bnamegen)) res.push_back(bnamegen);
      else if (!checkListVariable(res, bname)) res.push_back(bname);
    }
  }
}
void BaseTree::silenceUnused(){
  constexpr unsigned int ntrees = 2;
  TTree* trees[ntrees]={
    tree,
    failedtree
  };
  for (unsigned int it=0; it<ntrees; it++){
    if (!trees[it]) continue;
    trees[it]->SetBranchStatus("*", 0);

    std::vector<TString> currentBranchList;
    this->getValidBranchNamesWithoutAlias(trees[it], currentBranchList);

    for (TString const& bname:currentBranchList){
      trees[it]->SetBranchStatus(bname, 1);
      //cout << "Unmuting branch " << bname << endl;
    }
  }
}
void BaseTree::releaseBranch(TString branchname){
  const BranchType btype = searchBranchType(branchname);
  if (btype!=BranchType_unknown_t){
    constexpr unsigned int ntrees = 2;
    TTree* trees[ntrees]={
      tree,
      failedtree
    };
    for (unsigned int it=0; it<ntrees; it++){
      if (!trees[it]) continue;
      trees[it]->ResetBranchAddress(trees[it]->GetBranch(branchname));
      trees[it]->SetBranchStatus(branchname, 0);
    }
  }

  switch (btype){
  case BranchType_TBits_t:
    this->removeBranch<BranchType_TBits_t>(branchname);
    break;
  case BranchType_bool_t:
    this->removeBranch<BranchType_bool_t>(branchname);
    break;
  case BranchType_short_t:
    this->removeBranch<BranchType_short_t>(branchname);
    break;
  case BranchType_uint_t:
    this->removeBranch<BranchType_uint_t>(branchname);
    break;
  case BranchType_int_t:
    this->removeBranch<BranchType_int_t>(branchname);
    break;
  case BranchType_ulong_t:
    this->removeBranch<BranchType_ulong_t>(branchname);
    break;
  case BranchType_long_t:
    this->removeBranch<BranchType_long_t>(branchname);
    break;
  case BranchType_ulonglong_t:
    this->removeBranch<BranchType_ulonglong_t>(branchname);
    break;
  case BranchType_longlong_t:
    this->removeBranch<BranchType_longlong_t>(branchname);
    break;
  case BranchType_float_t:
    this->removeBranch<BranchType_float_t>(branchname);
    break;
  case BranchType_double_t:
    this->removeBranch<BranchType_double_t>(branchname);
    break;
  case BranchType_string_t:
    this->removeBranch<BranchType_string_t>(branchname);
    break;
  case BranchType_TString_t:
    this->removeBranch<BranchType_TString_t>(branchname);
    break;
  case BranchType_CMSLorentzVector_t:
    this->removeBranch<BranchType_CMSLorentzVector_t>(branchname);
    break;

  case BranchType_vbool_t:
    this->removeBranch<BranchType_vbool_t>(branchname);
    break;
  case BranchType_vshort_t:
    this->removeBranch<BranchType_vshort_t>(branchname);
    break;
  case BranchType_vuint_t:
    this->removeBranch<BranchType_vuint_t>(branchname);
    break;
  case BranchType_vint_t:
    this->removeBranch<BranchType_vint_t>(branchname);
    break;
  case BranchType_vulong_t:
    this->removeBranch<BranchType_vulong_t>(branchname);
    break;
  case BranchType_vlong_t:
    this->removeBranch<BranchType_vlong_t>(branchname);
    break;
  case BranchType_vulonglong_t:
    this->removeBranch<BranchType_vulonglong_t>(branchname);
    break;
  case BranchType_vlonglong_t:
    this->removeBranch<BranchType_vlonglong_t>(branchname);
    break;
  case BranchType_vfloat_t:
    this->removeBranch<BranchType_vfloat_t>(branchname);
    break;
  case BranchType_vdouble_t:
    this->removeBranch<BranchType_vdouble_t>(branchname);
    break;
  case BranchType_vstring_t:
    this->removeBranch<BranchType_vstring_t>(branchname);
    break;
  case BranchType_vTString_t:
    this->removeBranch<BranchType_vTString_t>(branchname);
    break;
  case BranchType_vCMSLorentzVector_t:
    this->removeBranch<BranchType_vCMSLorentzVector_t>(branchname);
    break;

  case BranchType_vvbool_t:
    this->removeBranch<BranchType_vvbool_t>(branchname);
    break;
  case BranchType_vvshort_t:
    this->removeBranch<BranchType_vvshort_t>(branchname);
    break;
  case BranchType_vvuint_t:
    this->removeBranch<BranchType_vvuint_t>(branchname);
    break;
  case BranchType_vvint_t:
    this->removeBranch<BranchType_vvint_t>(branchname);
    break;
  case BranchType_vvulong_t:
    this->removeBranch<BranchType_vvulong_t>(branchname);
    break;
  case BranchType_vvlong_t:
    this->removeBranch<BranchType_vvlong_t>(branchname);
    break;
  case BranchType_vvulonglong_t:
    this->removeBranch<BranchType_vvulonglong_t>(branchname);
    break;
  case BranchType_vvlonglong_t:
    this->removeBranch<BranchType_vvlonglong_t>(branchname);
    break;
  case BranchType_vvfloat_t:
    this->removeBranch<BranchType_vvfloat_t>(branchname);
    break;
  case BranchType_vvdouble_t:
    this->removeBranch<BranchType_vvdouble_t>(branchname);
    break;

  default:
    break;
  }
}
void BaseTree::setAutoSave(Long64_t fsave){
  if (receiver || !tree) return;
  tree->SetAutoSave(fsave);
}

bool BaseTree::isValidEvent() const{ return BaseTree::isValid(); } // To be overloaded in the daughter tree

void BaseTree::fill(){
  if (receiver || !tree) return;
  tree->Fill();
}
void BaseTree::writeToFile(TFile* file){
  if (receiver || !tree || !file || !(file->IsOpen() && !file->IsZombie())) return;
  file->WriteTObject(tree);
}

void BaseTree::writeSimpleEntries(std::vector<SimpleEntry>::iterator const& vecBegin, std::vector<SimpleEntry>::iterator const& vecEnd, BaseTree* const& tree){
  if (!tree) return;
  for (std::vector<SimpleEntry>::iterator it=vecBegin; it!=vecEnd; it++){
    SimpleEntry& entry = *it;
    if (it==vecBegin){
      for (auto itb=entry.namedbools.begin(); itb!=entry.namedbools.end(); itb++) tree->putBranch(itb->first, itb->second);
      for (auto itb=entry.namedshorts.begin(); itb!=entry.namedshorts.end(); itb++) tree->putBranch(itb->first, itb->second);
      for (auto itb=entry.nameduints.begin(); itb!=entry.nameduints.end(); itb++) tree->putBranch(itb->first, itb->second);
      for (auto itb=entry.namedints.begin(); itb!=entry.namedints.end(); itb++) tree->putBranch(itb->first, itb->second);
      for (auto itb=entry.namedulongs.begin(); itb!=entry.namedulongs.end(); itb++) tree->putBranch(itb->first, itb->second);
      for (auto itb=entry.namedlongs.begin(); itb!=entry.namedlongs.end(); itb++) tree->putBranch(itb->first, itb->second);
      for (auto itb=entry.namedulonglongs.begin(); itb!=entry.namedulonglongs.end(); itb++) tree->putBranch(itb->first, itb->second);
      for (auto itb=entry.namedlonglongs.begin(); itb!=entry.namedlonglongs.end(); itb++) tree->putBranch(itb->first, itb->second);
      for (auto itb=entry.namedfloats.begin(); itb!=entry.namedfloats.end(); itb++) tree->putBranch(itb->first, itb->second);
      for (auto itb=entry.nameddoubles.begin(); itb!=entry.nameddoubles.end(); itb++) tree->putBranch(itb->first, itb->second);
      for (auto itb=entry.namedCMSLorentzVectors.begin(); itb!=entry.namedCMSLorentzVectors.end(); itb++) tree->putBranch(itb->first, itb->second);
      for (auto itb=entry.namedVbools.begin(); itb!=entry.namedVbools.end(); itb++) tree->putBranch(itb->first, &(itb->second));
      for (auto itb=entry.namedVshorts.begin(); itb!=entry.namedVshorts.end(); itb++) tree->putBranch(itb->first, &(itb->second));
      for (auto itb=entry.namedVuints.begin(); itb!=entry.namedVuints.end(); itb++) tree->putBranch(itb->first, &(itb->second));
      for (auto itb=entry.namedVints.begin(); itb!=entry.namedVints.end(); itb++) tree->putBranch(itb->first, &(itb->second));
      for (auto itb=entry.namedVulongs.begin(); itb!=entry.namedVulongs.end(); itb++) tree->putBranch(itb->first, &(itb->second));
      for (auto itb=entry.namedVlongs.begin(); itb!=entry.namedVlongs.end(); itb++) tree->putBranch(itb->first, &(itb->second));
      for (auto itb=entry.namedVulonglongs.begin(); itb!=entry.namedVulonglongs.end(); itb++) tree->putBranch(itb->first, &(itb->second));
      for (auto itb=entry.namedVlonglongs.begin(); itb!=entry.namedVlonglongs.end(); itb++) tree->putBranch(itb->first, &(itb->second));
      for (auto itb=entry.namedVfloats.begin(); itb!=entry.namedVfloats.end(); itb++) tree->putBranch(itb->first, &(itb->second));
      for (auto itb=entry.namedVdoubles.begin(); itb!=entry.namedVdoubles.end(); itb++) tree->putBranch(itb->first, &(itb->second));
      for (auto itb=entry.namedVCMSLorentzVectors.begin(); itb!=entry.namedVCMSLorentzVectors.end(); itb++) tree->putBranch(itb->first, &(itb->second));
    }
    for (auto itb=entry.namedbools.begin(); itb!=entry.namedbools.end(); itb++) tree->setVal(itb->first, itb->second);
    for (auto itb=entry.namedshorts.begin(); itb!=entry.namedshorts.end(); itb++) tree->setVal(itb->first, itb->second);
    for (auto itb=entry.nameduints.begin(); itb!=entry.nameduints.end(); itb++) tree->setVal(itb->first, itb->second);
    for (auto itb=entry.namedints.begin(); itb!=entry.namedints.end(); itb++) tree->setVal(itb->first, itb->second);
    for (auto itb=entry.namedulongs.begin(); itb!=entry.namedulongs.end(); itb++) tree->setVal(itb->first, itb->second);
    for (auto itb=entry.namedlongs.begin(); itb!=entry.namedlongs.end(); itb++) tree->setVal(itb->first, itb->second);
    for (auto itb=entry.namedulonglongs.begin(); itb!=entry.namedulonglongs.end(); itb++) tree->setVal(itb->first, itb->second);
    for (auto itb=entry.namedlonglongs.begin(); itb!=entry.namedlonglongs.end(); itb++) tree->setVal(itb->first, itb->second);
    for (auto itb=entry.namedfloats.begin(); itb!=entry.namedfloats.end(); itb++) tree->setVal(itb->first, itb->second);
    for (auto itb=entry.nameddoubles.begin(); itb!=entry.nameddoubles.end(); itb++) tree->setVal(itb->first, itb->second);
    for (auto itb=entry.namedCMSLorentzVectors.begin(); itb!=entry.namedCMSLorentzVectors.end(); itb++) tree->setVal(itb->first, itb->second);
    for (auto itb=entry.namedVbools.begin(); itb!=entry.namedVbools.end(); itb++) tree->setVal(itb->first, &(itb->second));
    for (auto itb=entry.namedVshorts.begin(); itb!=entry.namedVshorts.end(); itb++) tree->setVal(itb->first, &(itb->second));
    for (auto itb=entry.namedVuints.begin(); itb!=entry.namedVuints.end(); itb++) tree->setVal(itb->first, &(itb->second));
    for (auto itb=entry.namedVints.begin(); itb!=entry.namedVints.end(); itb++) tree->setVal(itb->first, &(itb->second));
    for (auto itb=entry.namedVulongs.begin(); itb!=entry.namedVulongs.end(); itb++) tree->setVal(itb->first, &(itb->second));
    for (auto itb=entry.namedVlongs.begin(); itb!=entry.namedVlongs.end(); itb++) tree->setVal(itb->first, &(itb->second));
    for (auto itb=entry.namedVulonglongs.begin(); itb!=entry.namedVulonglongs.end(); itb++) tree->setVal(itb->first, &(itb->second));
    for (auto itb=entry.namedVlonglongs.begin(); itb!=entry.namedVlonglongs.end(); itb++) tree->setVal(itb->first, &(itb->second));
    for (auto itb=entry.namedVfloats.begin(); itb!=entry.namedVfloats.end(); itb++) tree->setVal(itb->first, &(itb->second));
    for (auto itb=entry.namedVdoubles.begin(); itb!=entry.namedVdoubles.end(); itb++) tree->setVal(itb->first, &(itb->second));
    for (auto itb=entry.namedVCMSLorentzVectors.begin(); itb!=entry.namedVCMSLorentzVectors.end(); itb++) tree->setVal(itb->first, &(itb->second));
    tree->fill();
  }
}
