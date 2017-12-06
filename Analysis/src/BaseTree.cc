#include "BaseTree.h"
#include "BaseTree.hpp"


using namespace std;


BaseTree::BaseTree() :
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
          failedtree = (TTree*) finput->Get(failedtreename);
          if (!failedtree) cout << "BaseTree::BaseTree(" << cinput << ") does not contain " << failedtreename << endl;
          hCounters = (TH1F*) finput->Get(countersname);
        }
      }
      else if (finput->IsOpen()){ finput->Close(); finput=nullptr; }
    }
  }
  curdir->cd(); // Return back to the directory before opening the input file
}
BaseTree::BaseTree(const TString treename) :
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
  HelperFunctions::cleanUnorderedMap(valshorts);
  HelperFunctions::cleanUnorderedMap(valuints);
  HelperFunctions::cleanUnorderedMap(valints);
  HelperFunctions::cleanUnorderedMap(valfloats);
  HelperFunctions::cleanUnorderedMap(valdoubles);
  if (!receiver){
    HelperFunctions::cleanUnorderedMap(valVshorts);
    HelperFunctions::cleanUnorderedMap(valVuints);
    HelperFunctions::cleanUnorderedMap(valVints);
    HelperFunctions::cleanUnorderedMap(valVfloats);
    HelperFunctions::cleanUnorderedMap(valVdoubles);

    delete hCounters;
    delete failedtree;
    delete tree;
  }
}

BaseTree::BranchType BaseTree::searchBranchType(TString branchname) const{
  if (valshorts.find(branchname)!=valshorts.cend()) return BranchType_short_t;
  else if (valuints.find(branchname)!=valuints.cend()) return BranchType_uint_t;
  else if (valints.find(branchname)!=valints.cend()) return BranchType_int_t;
  else if (valfloats.find(branchname)!=valfloats.cend()) return BranchType_float_t;
  else if (valdoubles.find(branchname)!=valdoubles.cend()) return BranchType_double_t;

  else if (valVshorts.find(branchname)!=valVshorts.cend()) return BranchType_vshort_t;
  else if (valVuints.find(branchname)!=valVuints.cend()) return BranchType_vuint_t;
  else if (valVints.find(branchname)!=valVints.cend()) return BranchType_vint_t;
  else if (valVfloats.find(branchname)!=valVfloats.cend()) return BranchType_vfloat_t;
  else if (valVdoubles.find(branchname)!=valVdoubles.cend()) return BranchType_vdouble_t;

  else return BranchType_unknown_t;
}

bool BaseTree::getSelectedEvent(int ev){
  resetBranches();
  bool result=false;
  if (tree && ev<tree->GetEntries()) result = (tree->GetEntry(ev)>0);
  if (result){
    currentEvent = ev;
    currentTree = tree;
  }
  return result;
}
bool BaseTree::getFailedEvent(int ev){
  resetBranches();
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
  resetBranches();
  if (tmpTree){
    tmpTree->GetEntry(tmpEv);
    currentEvent = tmpEv;
    currentTree = tmpTree;
  }
}

int BaseTree::getSelectedNEvents(){ return (tree ? tree->GetEntries() : 0); }
int BaseTree::getFailedNEvents(){ return (failedtree ? failedtree->GetEntries() : 0); }
int BaseTree::getNEvents(){ return (this->getSelectedNEvents() + this->getFailedNEvents()); }

bool BaseTree::isValid() const{ return valid; }
bool BaseTree::branchExists(TString branchname, BranchType* type){
  BranchType theType = searchBranchType(branchname);
  if (type) *type=theType;
  return (theType!=BranchType_unknown_t);
}

void BaseTree::resetBranches(){
  currentEvent = -1;
  currentTree = nullptr;

  this->resetBranch<BaseTree::BranchType_short_t>();
  this->resetBranch<BaseTree::BranchType_uint_t>();
  this->resetBranch<BaseTree::BranchType_int_t>();
  this->resetBranch<BaseTree::BranchType_float_t>();
  this->resetBranch<BaseTree::BranchType_double_t>();
  if (!receiver){
    this->resetBranch<BaseTree::BranchType_vshort_t>();
    this->resetBranch<BaseTree::BranchType_vuint_t>();
    this->resetBranch<BaseTree::BranchType_vint_t>();
    this->resetBranch<BaseTree::BranchType_vfloat_t>();
    this->resetBranch<BaseTree::BranchType_vdouble_t>();
  }
}

void BaseTree::silenceUnused(){
  const unsigned int ntrees = 2;
  TTree* trees[ntrees]={
    tree,
    failedtree
  };
  for (unsigned int it=0; it<ntrees; it++){
    if (!trees[it]) continue;
    const TList* blist = (const TList*)trees[it]->GetListOfBranches();
    for (int ib=0; ib<blist->GetSize(); ib++){
      TString bname = blist->At(ib)->GetName();
      if (searchBranchType(bname)==BranchType_unknown_t) trees[it]->SetBranchStatus(bname, 0);
    }
  }
}
void BaseTree::releaseBranch(TString branchname){
  const BranchType btype = searchBranchType(branchname);
  if (btype!=BranchType_unknown_t){
    const unsigned int ntrees = 2;
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
  case BranchType_short_t:
    this->removeBranch<BranchType_short_t>(branchname);
    break;
  case BranchType_uint_t:
    this->removeBranch<BranchType_uint_t>(branchname);
    break;
  case BranchType_int_t:
    this->removeBranch<BranchType_int_t>(branchname);
    break;
  case BranchType_float_t:
    this->removeBranch<BranchType_float_t>(branchname);
    break;
  case BranchType_double_t:
    this->removeBranch<BranchType_double_t>(branchname);
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
  case BranchType_vfloat_t:
    this->removeBranch<BranchType_vfloat_t>(branchname);
    break;
  case BranchType_vdouble_t:
    this->removeBranch<BranchType_vdouble_t>(branchname);
    break;
  default:
    break;
  }
}

bool BaseTree::isValidEvent() const{ return true; } // To be overloaded in CJLSTTree to check for POWHEG mH<300 GeV, ZZMass < mH

