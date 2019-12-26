#include "TDirectory.h"
#include "HostHelpersCore.h"
#include "BaseTree.h"
#include "BaseTree.hpp"
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;


bool BaseTree::robustSaveWrite = false;
BaseTree::BaseTree() :
  sampleIdentifier(""),
  finput(nullptr),
  tree(nullptr),
  failedtree(nullptr),
  hCounters(nullptr),
  valid(false),
  receiver(true),
  acquireTreePossession(!receiver),
  isTChain(false),
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
  acquireTreePossession(!receiver),
  isTChain(false),
  currentEvent(-1),
  currentTree(nullptr)
{
  TDirectory* curdir = gDirectory; // Save current directory to return back to it later
  if (cinput.Contains("*")){ // Use a TChain
    isTChain=true;
    if (treename!=""){
      tree = new TChain(treename);
      valid = (tree->Add(cinput) != 0);
      if (!valid){
        delete tree;
        tree=nullptr;
      }
    }
    if (failedtreename!=""){
      failedtree = new TChain(failedtreename);
      valid = (failedtree->Add(cinput) != 0);
      if (!valid){
        delete failedtree;
        failedtree=nullptr;
      }
    }
    if (countersname!=""){
      MELAerr << "BaseTree::BaseTree: Cannot add histograms in chain mode." << endl;
      assert(0);
    }
  }
  else if (HostHelpers::FileReadable(cinput.Data())){
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
  acquireTreePossession(!receiver),
  isTChain(false),
  currentEvent(-1),
  currentTree(nullptr)
{}
BaseTree::BaseTree(TFile* finput_, TTree* tree_, TTree* failedtree_, TH1F* hCounters_, bool receiver_override) :
  sampleIdentifier(""),
  finput(finput_),
  tree(tree_),
  failedtree(failedtree_),
  hCounters(hCounters_),
  valid(false),
  receiver(receiver_override || finput!=nullptr),
  acquireTreePossession(!receiver),
  isTChain(false),
  currentEvent(-1),
  currentTree(nullptr)
{
  if (finput){
    if (finput->IsOpen() && !finput->IsZombie()){
      if (tree || failedtree) valid = true;
    }
    else if (finput->IsOpen()){ finput->Close(); finput=nullptr; tree=nullptr; failedtree=nullptr; }
  }
  else{
    if (tree || failedtree) valid = true;
  }
}

BaseTree::~BaseTree(){
#define SIMPLE_DATA_INPUT_DIRECTIVE(name, type, default_value) HelperFunctions::cleanUnorderedMap(val##name##s);
#define VECTOR_DATA_INPUT_DIRECTIVE(name, type) HelperFunctions::cleanUnorderedMap(valV##name##s);
#define DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(name, type) HelperFunctions::cleanUnorderedMap(valVV##name##s);
  SIMPLE_DATA_INPUT_DIRECTIVES
  if (!receiver){
    VECTOR_DATA_INPUT_DIRECTIVES
    DOUBLEVECTOR_DATA_INPUT_DIRECTIVES
    if (acquireTreePossession){
      delete hCounters;
      delete failedtree;
      delete tree;
    }
  }
  else if (isTChain){
    delete hCounters;
    delete failedtree;
    delete tree;
  }
#undef SIMPLE_DATA_INPUT_DIRECTIVE
#undef VECTOR_DATA_INPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_INPUT_DIRECTIVE

  if (finput && finput->IsOpen()) finput->Close();
}

BaseTree::BranchType BaseTree::searchBranchType(TString branchname) const{
#define SIMPLE_DATA_INPUT_DIRECTIVE(name, type, default_value) else if (val##name##s.find(branchname)!=val##name##s.cend()) return BranchType_##name##_t;
#define VECTOR_DATA_INPUT_DIRECTIVE(name, type) else if (valV##name##s.find(branchname)!=valV##name##s.cend()) return BranchType_v##name##_t;
#define DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(name, type) else if (valVV##name##s.find(branchname)!=valVV##name##s.cend()) return BranchType_vv##name##_t;
  if (false) return BranchType_unknown_t;
  SIMPLE_DATA_INPUT_DIRECTIVES
  VECTOR_DATA_INPUT_DIRECTIVES
  DOUBLEVECTOR_DATA_INPUT_DIRECTIVES
  else return BranchType_unknown_t;
#undef SIMPLE_DATA_INPUT_DIRECTIVE
#undef VECTOR_DATA_INPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_INPUT_DIRECTIVE
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

// Overloads of getNEvents
unsigned int BaseTree::getNGenNoPU(){ int res = this->getNEvents(); return std::max(res, 0); }
float BaseTree::getNGenWithPU(){ return this->getNEvents(); }

bool BaseTree::isValid() const{ return valid; }
bool BaseTree::branchExists(TString branchname, BranchType* type){
  BranchType theType = searchBranchType(branchname);
  if (type) *type=theType;
  return (theType!=BranchType_unknown_t);
}

void BaseTree::print() const{
#define SIMPLE_DATA_INPUT_DIRECTIVE(name, type, default_value) for (auto const& it:val##name##s){ if (it.second){ MELAout << "\t- " << it.first << " value: " << it.second->first << " (address: " << &(it.second->first) << ")" << endl; } }
#define VECTOR_DATA_INPUT_DIRECTIVE(name, type) for (auto const& it:valV##name##s){ MELAout << "\t- " << it.first << " value: "; if (it.second) MELAout << *(it.second); else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
#define DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(name, type) for (auto const& it:valVV##name##s){ MELAout << "\t- " << it.first << " value: "; if (it.second){ for (auto const& v:*(it.second)){ MELAout << "{ " << v << " }"; } } else MELAout << "null"; MELAout << " (address: " << it.second << ")" << endl; }
  SIMPLE_DATA_INPUT_DIRECTIVES
  VECTOR_DATA_INPUT_DIRECTIVES
  DOUBLEVECTOR_DATA_INPUT_DIRECTIVES
#undef SIMPLE_DATA_INPUT_DIRECTIVE
#undef VECTOR_DATA_INPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_INPUT_DIRECTIVE
}

void BaseTree::resetBranches(){
  currentEvent = -1;
  currentTree = nullptr;

#define SIMPLE_DATA_INPUT_DIRECTIVE(name, type, default_value) this->resetBranch<BaseTree::BranchType_##name##_t>();
#define VECTOR_DATA_INPUT_DIRECTIVE(name, type) this->resetBranch<BaseTree::BranchType_v##name##_t>();
#define DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(name, type) this->resetBranch<BaseTree::BranchType_vv##name##_t>();
  SIMPLE_DATA_INPUT_DIRECTIVES
  if (!receiver){
    VECTOR_DATA_INPUT_DIRECTIVES
    DOUBLEVECTOR_DATA_INPUT_DIRECTIVES
  }
#undef SIMPLE_DATA_INPUT_DIRECTIVE
#undef VECTOR_DATA_INPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_INPUT_DIRECTIVE
}

void BaseTree::getValidBranchNamesWithoutAlias(TTree* t, std::vector<TString>& res, bool check_linked) const{
  if (!t) return;

  const TList* alist = (const TList*) t->GetListOfAliases();
  const TList* llist = (const TList*) t->GetListOfLeaves();
  const TList* blist = (const TList*) t->GetListOfBranches();
  // First check all aliases and record the proper names
  if (alist){
    for (int ib=0; ib<alist->GetSize(); ib++){
      TString bname = alist->At(ib)->GetName();
      TString bnameproper = t->GetAlias(bname);
      TString bnamegen="";
      if (bnameproper.Contains(".")){
        std::vector<TString> tmplist;
        splitOptionRecursive(bnameproper, tmplist, '.');
        if (!tmplist.empty()) bnamegen = tmplist.front() + "*";
      }
      if (searchBranchType(bname)!=BranchType_unknown_t || !check_linked){
        if (bnamegen!="" && !checkListVariable(res, bnamegen)) res.push_back(bnamegen);
        else if (!checkListVariable(res, bnameproper)) res.push_back(bnameproper);
      }
    }
  }
  // Then check all leaves
  if (llist){
    for (int ib=0; ib<llist->GetSize(); ib++){
      TString bname = llist->At(ib)->GetName();
      TString bnamegen="";
      if (bname.Contains(".")){
        std::vector<TString> tmplist;
        splitOptionRecursive(bname, tmplist, '.');
        if (!tmplist.empty()) bnamegen = tmplist.front() + "*";
      }
      if (searchBranchType(bname)!=BranchType_unknown_t || !check_linked){
        if (bnamegen!="" && !checkListVariable(res, bnamegen)) res.push_back(bnamegen);
        else if (!checkListVariable(res, bname)) res.push_back(bname);
      }
    }
  }
  if (blist){
    // Then check all branches
    for (int ib=0; ib<blist->GetSize(); ib++){
      TString bname = blist->At(ib)->GetName();
      TString bnamegen="";
      if (bname.Contains(".")){
        std::vector<TString> tmplist;
        splitOptionRecursive(bname, tmplist, '.');
        if (!tmplist.empty()) bnamegen = tmplist.front() + "*";
      }
      if (searchBranchType(bname)!=BranchType_unknown_t || !check_linked){
        if (bnamegen!="" && !checkListVariable(res, bnamegen)) res.push_back(bnamegen);
        else if (!checkListVariable(res, bname)) res.push_back(bname);
      }
    }
  }
}

void BaseTree::getValidBranchNamesWithoutAlias(std::vector<TString>& res, bool check_linked) const{
  constexpr unsigned int ntrees = 2;
  TTree* trees[ntrees]={
    tree,
    failedtree
  };
  for (unsigned int it=0; it<ntrees; it++) this->getValidBranchNamesWithoutAlias(trees[it], res, check_linked);
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
    this->getValidBranchNamesWithoutAlias(trees[it], currentBranchList, true);

    for (TString const& bname:currentBranchList){
      trees[it]->SetBranchStatus(bname, 1);
      //cout << "Unmuting branch " << bname << endl;
    }
  }
}
void BaseTree::unmuteAllBranches(){
  constexpr unsigned int ntrees = 2;
  TTree* trees[ntrees]={
    tree,
    failedtree
  };
  for (unsigned int it=0; it<ntrees; it++){
    if (!trees[it]) continue;
    trees[it]->SetBranchStatus("*", 1);
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

#define SIMPLE_DATA_INPUT_DIRECTIVE(name, type, default_value) \
  case BranchType_##name##_t: \
  this->removeBranch<BranchType_##name##_t>(branchname); \
  break;
#define VECTOR_DATA_INPUT_DIRECTIVE(name, type) \
  case BranchType_v##name##_t: \
  this->removeBranch<BranchType_v##name##_t>(branchname); \
  break;
#define DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(name, type) \
  case BranchType_vv##name##_t: \
  this->removeBranch<BranchType_vv##name##_t>(branchname); \
  break;

  switch (btype){
  SIMPLE_DATA_INPUT_DIRECTIVES
  VECTOR_DATA_INPUT_DIRECTIVES
  DOUBLEVECTOR_DATA_INPUT_DIRECTIVES

  default:
    break;
  }

#undef SIMPLE_DATA_INPUT_DIRECTIVE
#undef VECTOR_DATA_INPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_INPUT_DIRECTIVE
}

void BaseTree::setAcquireTreePossession(bool flag){ acquireTreePossession = flag; }

void BaseTree::setAutoSave(Long64_t fsave){
  if (receiver || !tree) return;
  tree->SetAutoSave(fsave);
}
long long BaseTree::doAutoSave(const char* opts){
  if (receiver || !tree) return 0;
  return tree->AutoSave(opts);
}

void BaseTree::setAutoFlush(Long64_t fflush){
  if (receiver || !tree) return;
  tree->SetAutoFlush(fflush);
}
int BaseTree::doFlushBaskets(){
  if (receiver || !tree) return 0;
  return tree->FlushBaskets();
}

bool BaseTree::isValidEvent() const{ return BaseTree::isValid(); } // To be overloaded in the daughter tree

void BaseTree::fill(){
  if (receiver || !tree) return;
  tree->Fill();
}
void BaseTree::writeToFile(TFile* file){
  if (receiver || !tree || !file || !file->IsOpen() || file->IsZombie()) return;
  if (robustSaveWrite) file->WriteTObject(tree, nullptr, "WriteDelete");
  else file->WriteTObject(tree);
}

void BaseTree::setRobustSaveWrite(bool flag){ BaseTree::robustSaveWrite = flag; }
void BaseTree::writeSimpleEntries(std::vector<SimpleEntry>::iterator const& vecBegin, std::vector<SimpleEntry>::iterator const& vecEnd, BaseTree* const& tree){
  if (!tree) return;
  for (std::vector<SimpleEntry>::iterator it=vecBegin; it!=vecEnd; it++){
    SimpleEntry& entry = *it;
    if (it==vecBegin){
#define SIMPLE_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=entry.named##name_t##s.begin(); itb!=entry.named##name_t##s.end(); itb++) tree->putBranch(itb->first, itb->second);
#define VECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=entry.namedV##name_t##s.begin(); itb!=entry.namedV##name_t##s.end(); itb++) tree->putBranch(itb->first, &(itb->second));
#define DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=entry.namedVV##name_t##s.begin(); itb!=entry.namedVV##name_t##s.end(); itb++) tree->putBranch(itb->first, &(itb->second));
      SIMPLE_DATA_OUTPUT_DIRECTIVES
      VECTOR_DATA_OUTPUT_DIRECTIVES
      DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVES
#undef SIMPLE_DATA_OUTPUT_DIRECTIVE
#undef VECTOR_DATA_OUTPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE
    }
#define SIMPLE_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=entry.named##name_t##s.begin(); itb!=entry.named##name_t##s.end(); itb++) tree->setVal(itb->first, itb->second);
#define VECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=entry.namedV##name_t##s.begin(); itb!=entry.namedV##name_t##s.end(); itb++) tree->setVal(itb->first, &(itb->second));
#define DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=entry.namedVV##name_t##s.begin(); itb!=entry.namedVV##name_t##s.end(); itb++) tree->setVal(itb->first, &(itb->second));
    SIMPLE_DATA_OUTPUT_DIRECTIVES
    VECTOR_DATA_OUTPUT_DIRECTIVES
    DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVES
#undef SIMPLE_DATA_OUTPUT_DIRECTIVE
#undef VECTOR_DATA_OUTPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE
    tree->fill();
  }

  // Save if this flag is specified
  if (robustSaveWrite){
    tree->doAutoSave("FlushBaskets");
  }
}
