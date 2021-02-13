#include "TDirectory.h"
#include "TRegexp.h"
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
    isTChain = true;
    valid = true;

    std::vector<TString> treenamelist; treenamelist.reserve(2);
    if (treename!="") treenamelist.push_back(treename);
    if (failedtreename!="") treenamelist.push_back(failedtreename);
    std::vector< std::vector<TString> > valid_files;
    valid = getValidFilesForTreeList(cinput, treenamelist, valid_files);

    if (valid){
      if (treename!=""){
        TChain* tc = new TChain(treename);
        for (auto const& fname:valid_files.front()) tc->Add(fname);
        tree = tc;
      }
      if (failedtreename!=""){
        TChain* tc = new TChain(failedtreename);
        auto const& vfiles = (treename!="" ? valid_files.back() : valid_files.front());
        for (auto const& fname:vfiles) tc->Add(fname);
        failedtree = tc;
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
        valid = true;
        finput->cd();
        if (treename!=""){
          tree = (TTree*) finput->Get(treename);
          if (!tree) cout << "BaseTree::BaseTree(" << cinput << ") does not contain " << treename << endl;
          valid &= (tree!=nullptr);
        }
        if (failedtreename!=""){
          failedtree = (TTree*) finput->Get(failedtreename);
          if (!failedtree) cout << "BaseTree::BaseTree(" << cinput << ") does not contain " << failedtreename << endl;
          valid &= (failedtree!=nullptr);
        }
        if (countersname!=""){
          hCounters = (TH1F*) finput->Get(countersname);
          if (!hCounters) cout << "BaseTree::BaseTree(" << cinput << ") does not contain " << countersname << endl;
          valid &= (hCounters!=nullptr);
        }
        if (!valid){
          tree = nullptr;
          failedtree = nullptr;
          hCounters = nullptr;
          finput->Close();
          finput = nullptr;
        }
      }
      else{
        if (finput->IsOpen()) finput->Close();
        else delete finput;
        finput=nullptr;
      }
    }
  }

  treelist.reserve(2);
  if (tree) treelist.push_back(tree);
  if (failedtree) treelist.push_back(failedtree);

  curdir->cd(); // Return back to the directory before opening the input file
}
BaseTree::BaseTree(const TString cinput, std::vector<TString> const& treenames, const TString countersname) :
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
    isTChain = true;
    valid = true;
    treelist.reserve(treenames.size());
    std::vector< std::vector<TString> > valid_files;
    valid = getValidFilesForTreeList(cinput, treenames, valid_files);
    if (valid){
      unsigned int itree=0;
      for (auto const& treename:treenames){
        if (treename!=""){
          TChain* tc = new TChain(treename);
          for (auto const& fname:valid_files.at(itree)) tc->Add(fname);
          treelist.push_back(tc);
        }
        itree++;
      }
    }
    if (!valid){
      for (auto& tt:treelist) delete tt;
      treelist.clear();
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
        valid = true;
        finput->cd();
        treelist.reserve(treenames.size());
        for (auto const& treename:treenames){
          if (treename!=""){
            TTree* tt = (TTree*) finput->Get(treename);
            if (!tt) cout << "BaseTree::BaseTree(" << cinput << ") does not contain " << treename << endl;
            else treelist.push_back(tt);
            valid &= (tt!=nullptr);
          }
        }
        if (countersname!=""){
          hCounters = (TH1F*) finput->Get(countersname);
          if (!hCounters) cout << "BaseTree::BaseTree(" << cinput << ") does not contain " << countersname << endl;
          valid &= (hCounters!=nullptr);
        }
        if (!valid){
          for (auto& tt:treelist) delete tt;
          treelist.clear();
          hCounters = nullptr;
          finput->Close();
          finput = nullptr;
        }
      }
      else{
        if (finput->IsOpen()) finput->Close();
        else delete finput;
        finput=nullptr;
      }
    }
  }

  if (!treelist.empty()) tree = treelist.front();
  if (treelist.size()>1) failedtree = treelist.back();

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
{
  treelist.reserve(1);
  treelist.push_back(tree);
}
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
    else{
      if (finput->IsOpen()) finput->Close();
      else delete finput;
      finput=nullptr;
      tree=nullptr;
      failedtree=nullptr;
      hCounters=nullptr;
    }
  }
  else{
    if (tree || failedtree) valid = true;
  }

  treelist.reserve(2);
  if (tree) treelist.push_back(tree);
  if (failedtree) treelist.push_back(failedtree);
}
BaseTree::BaseTree(TFile* finput_, std::vector<TTree*> const& treelist_, TH1F* hCounters_, bool receiver_override) :
  sampleIdentifier(""),
  finput(finput_),
  tree(nullptr),
  failedtree(nullptr),
  hCounters(hCounters_),
  valid(false),
  receiver(receiver_override || finput!=nullptr),
  acquireTreePossession(!receiver),
  isTChain(false),
  currentEvent(-1),
  currentTree(nullptr)
{
  treelist.reserve(treelist_.size());
  if (finput){
    if (finput->IsOpen() && !finput->IsZombie()){
      for (auto const& tt:treelist_){
        if (tt) treelist.push_back(tt);
      }
    }
  }
  else{
    for (auto const& tt:treelist_){
      if (tt) treelist.push_back(tt);
    }
  }

  valid = (!treelist.empty() && treelist.size()==treelist_.size());
  if (!valid){
    tree=nullptr;
    failedtree=nullptr;
    hCounters=nullptr;
    if (finput){
      if (finput->IsOpen()) finput->Close();
      else delete finput;
      finput=nullptr;
    }
  }

  if (!treelist.empty()) tree = treelist.front();
  if (treelist.size()>1) failedtree = treelist.back();
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
      for (auto& tt:treelist) delete tt;
    }
  }
  else if (isTChain){
    delete hCounters;
    for (auto& tt:treelist) delete tt;
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
std::vector<TTree*>& BaseTree::getValidTrees(){ return treelist; }
TFile const* BaseTree::getInputFile() const{ return finput; }
TTree const* BaseTree::getSelectedTree() const{ return tree; }
TTree const* BaseTree::getFailedTree() const{ return failedtree; }
std::vector<TTree*> const& BaseTree::getValidTrees() const{ return treelist; }

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
  int n_acc=0;
  for (auto& tt:treelist){
    int nEntries = tt->GetEntries();
    if (ev<n_acc+nEntries){
      this->resetBranches();
      bool result = (tt->GetEntry(ev-n_acc)>0);
      if (result){
        currentEvent = ev-n_acc;
        currentTree = tt;
      }
      return result;
    }
    n_acc += nEntries;
  }
  return false;
}
bool BaseTree::updateBranch(int ev, TString const& bname, bool check_linked){
  if (check_linked){
    if (searchBranchType(bname)==BranchType_unknown_t) MELAerr << "BaseTree::updateBranch: Branch " << bname << " is not linked." << endl;
    return false;
  }
  int n_acc=0;
  for (auto& tt:treelist){
    int nEntries = tt->GetEntries();
    if (ev<n_acc+nEntries){
      int ev_tree = ev-n_acc;
      TBranch* tbr = nullptr;
      if (isTChain){
        // If the tree is a TChain, the entry index and the total number of entries
        // need to correspond to the loaded tree of the TChain.
        TChain* tc = dynamic_cast<TChain*>(tt);
        ev_tree = tc->LoadTree(ev_tree);
        TTree* tc_tree = tc->GetTree();
        nEntries = (tc_tree ? tc_tree->GetEntries() : -1);
        // It is important to get the branch AFTER loading the tree
        // because one might get the branch of another tree accidentally.
        tbr = (tc_tree ? tc_tree->GetBranch(bname) : nullptr);
      }
      else tbr = tt->GetBranch(bname);
      if (tbr && tbr->GetEntries()!=nEntries){
        MELAerr << "BaseTree::updateBranch: Branch " << bname << " has " << tbr->GetEntries() << " entries != " << nEntries << "." << endl;
        assert(0);
      }
      return (tbr && tbr->GetEntries()==nEntries && ev_tree>=0 && tbr->GetEntry(ev_tree)>0);
    }
    n_acc += nEntries;
  }
  return false;
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
bool BaseTree::getCurrentEventInfo(TTree*& currentTree_, int& currentEvent_) const{
  currentTree_ = this->currentTree;
  currentEvent_ = this->currentEvent;
  return (currentTree_!=nullptr && currentEvent_>=0);
}
bool BaseTree::isSameEvent(TTree* const& currentTree_, int const& currentEvent_) const{
  return (
    currentTree_ && currentEvent_>=0
    &&
    this->currentTree == currentTree_
    &&
    this->currentEvent == currentEvent_
    );
}

int BaseTree::getSelectedNEvents() const{ return (tree ? tree->GetEntries() : 0); }
int BaseTree::getFailedNEvents() const{ return (failedtree ? failedtree->GetEntries() : 0); }
int BaseTree::getNEvents() const{
  int n_acc=0;
  for (auto const& tt:treelist) n_acc += tt->GetEntries();
  return n_acc;
}

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
      auto const& bmem = alist->At(ib);
      if (!bmem) continue;
      TString bname = bmem->GetName();
      TString bnameproper = t->GetAlias(bname);
      TString bnamegen="";
      if (bnameproper.Contains(".")){
        std::vector<TString> tmplist;
        splitOptionRecursive(bnameproper, tmplist, '.');
        if (!tmplist.empty()) bnamegen = tmplist.front() + "*";
      }
      if (!check_linked || searchBranchType(bname)!=BranchType_unknown_t){
        if (bnamegen!="" && !checkListVariable(res, bnamegen)) res.push_back(bnamegen);
        else if (!checkListVariable(res, bnameproper)) res.push_back(bnameproper);
      }
    }
  }
  // Then check all leaves
  if (llist){
    for (int ib=0; ib<llist->GetSize(); ib++){
      auto const& bmem = llist->At(ib);
      if (!bmem) continue;
      TString bname = bmem->GetName();
      TString bnamegen="";
      if (bname.Contains(".")){
        std::vector<TString> tmplist;
        splitOptionRecursive(bname, tmplist, '.');
        if (!tmplist.empty()) bnamegen = tmplist.front() + "*";
      }
      if (!check_linked || searchBranchType(bname)!=BranchType_unknown_t){
        if (bnamegen!="" && !checkListVariable(res, bnamegen)) res.push_back(bnamegen);
        else if (!checkListVariable(res, bname)) res.push_back(bname);
      }
    }
  }
  // Then check all branches
  if (blist){
    for (int ib=0; ib<blist->GetSize(); ib++){
      auto const& bmem = blist->At(ib);
      if (!bmem) continue;
      TString bname = bmem->GetName();
      TString bnamegen="";
      if (bname.Contains(".")){
        std::vector<TString> tmplist;
        splitOptionRecursive(bname, tmplist, '.');
        if (!tmplist.empty()) bnamegen = tmplist.front() + "*";
      }
      if (!check_linked || searchBranchType(bname)!=BranchType_unknown_t){
        if (bnamegen!="" && !checkListVariable(res, bnamegen)) res.push_back(bnamegen);
        else if (!checkListVariable(res, bname)) res.push_back(bname);
      }
    }
  }
}

void BaseTree::getValidBranchNamesWithoutAlias(std::vector<TString>& res, bool check_linked) const{
  for (auto const& tt:treelist) this->getValidBranchNamesWithoutAlias(tt, res, check_linked);
}

void BaseTree::silenceUnused(){
  for (auto const& tt:treelist){
    tt->SetBranchStatus("*", 0);

    std::vector<TString> currentBranchList;
    this->getValidBranchNamesWithoutAlias(tt, currentBranchList, true);

    for (TString const& bname:currentBranchList){
      tt->SetBranchStatus(bname, 1);
      //cout << "Unmuting branch " << bname << endl;
    }
  }
}
void BaseTree::unmuteAllBranches(){
  for (auto const& tt:treelist) tt->SetBranchStatus("*", 1);
}
void BaseTree::releaseBranch(TString branchname){
  const BranchType btype = searchBranchType(branchname);
  if (btype!=BranchType_unknown_t){
    for (auto& tt:treelist){
      tt->ResetBranchAddress(tt->GetBranch(branchname));
      tt->SetBranchStatus(branchname, 0);
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
void BaseTree::muteAllBranchesExcept(std::vector<TString> const& bnames_excepted){
  for (auto& tt:treelist){
    tt->SetBranchStatus("*", 0);

    std::vector<TString> currentBranchList;
    this->getValidBranchNamesWithoutAlias(tt, currentBranchList, true);

    for (TString const& bname:bnames_excepted){
      if (HelperFunctions::checkListVariable(currentBranchList, bname)) tt->SetBranchStatus(bname, 1);
      else MELAerr << "BaseTree::muteAllBranchesExcept: Branch " << bname << " was not booked. Please book this branch and re-run this function." << endl;
    }
  }
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
void BaseTree::writeToDirectory(TDirectory* dir){
  if (receiver || !tree || !dir) return;
  dir->WriteTObject(tree);
}

bool BaseTree::getValidFilesForTreeList(TString cinput, std::vector<TString> const& treenames, std::vector< std::vector<TString> >& res) const{
  TDirectory* curdir = gDirectory; // Save current directory to return back to it later

  HostHelpers::ExpandEnvironmentVariables(cinput);

  std::vector<TString> fnames;
  if (cinput.Contains("*")){
    TString cinputcore, fpattern;
    size_t ipos = cinput.Last('/');
    if (ipos>0){
      cinputcore = cinput(0, ipos+1);
      fpattern = cinput(ipos+1, cinput.Length());
    }
    else if (ipos==0){
      MELAerr << "BaseTree::getValidFilesForTreeList: Invalid pattern " << cinput << endl;
      return false;
    }
    else{
      fpattern = cinput;
      cinputcore = "./";
    }
    std::vector<TString> fnames_all = SampleHelpers::lsdir(cinputcore);
    for (auto const& fname:fnames_all){
      if (fname.Index(TRegexp(fpattern.Data(), true))>=0) fnames.push_back(cinputcore + fname);
    }
  }
  else fnames.push_back(cinput);

  res.assign(treenames.size(), std::vector<TString>()); for (auto& vv:res) vv.reserve(fnames.size());
  for (auto& fname:fnames){
    TFile* ftmp = TFile::Open(fname, "read");
    if (ftmp){
      if (ftmp->IsOpen() && !ftmp->IsZombie()){
        ftmp->cd();
        unsigned int it=0;
        for (auto const& treename:treenames){
          if (treename!=""){
            TTree* tt = (TTree*) ftmp->Get(treename);
            if (!tt){
              MELAerr << "BaseTree::getValidFilesForTreeList: " << fname << " does not contain " << treename << endl;
              continue;
            }
            else if (tt->GetEntries()==0){
              std::vector<TString> bnames;
              BaseTree::getValidBranchNamesWithoutAlias(tt, bnames, false);
              if (bnames.empty()){
                MELAerr << "BaseTree::getValidFilesForTreeList: " << fname << " contains " << treename << ", but the tree has no branches." << endl;
                continue;
              }
            }
            res.at(it).push_back(fname);
          }
          it++;
        }
        ftmp->Close();
      }
      else{
        if (ftmp->IsOpen()) ftmp->Close();
        else delete ftmp;

        MELAerr << "BaseTree::getValidFilesForTreeList: File " << fname << " is not readable and was set to a zombie state! Aborting operation..." << endl;
        curdir->cd();
        return false;
      }
    }
    else{
      MELAerr << "BaseTree::getValidFilesForTreeList: File " << fname << " is not readable and could not be opened! Aborting operation..." << endl;
      curdir->cd();
      return false;
    }

    curdir->cd();
  }

  {
    unsigned int it=0;
    for (auto& vv:res){
      if (treenames.at(it)!=""){
        if (vv.empty()) return false;
        else if (vv.size() == fnames.size()){
          vv.clear();
          vv.push_back(cinput);
        }
      }
      it++;
    }
  }

  return true;
}

void BaseTree::setRobustSaveWrite(bool flag){ BaseTree::robustSaveWrite = flag; }
void BaseTree::writeSimpleEntries(std::vector<SimpleEntry>::const_iterator const& vecBegin, std::vector<SimpleEntry>::const_iterator const& vecEnd, BaseTree* const& tree_, bool createBranches){
  if (!tree_) return;
  for (std::vector<SimpleEntry>::const_iterator it=vecBegin; it!=vecEnd; it++){
    SimpleEntry const& entry = *it;
    if (createBranches && it==vecBegin){
#define SIMPLE_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=entry.named##name_t##s.cbegin(); itb!=entry.named##name_t##s.cend(); itb++) tree_->putBranch(itb->first, itb->second);
#define VECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=entry.namedV##name_t##s.cbegin(); itb!=entry.namedV##name_t##s.cend(); itb++) tree_->putBranch(itb->first, &(itb->second));
#define DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=entry.namedVV##name_t##s.cbegin(); itb!=entry.namedVV##name_t##s.cend(); itb++) tree_->putBranch(itb->first, &(itb->second));
      SIMPLE_DATA_OUTPUT_DIRECTIVES
      VECTOR_DATA_OUTPUT_DIRECTIVES
      DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVES
#undef SIMPLE_DATA_OUTPUT_DIRECTIVE
#undef VECTOR_DATA_OUTPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE
    }
#define SIMPLE_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=entry.named##name_t##s.cbegin(); itb!=entry.named##name_t##s.cend(); itb++) tree_->setVal(itb->first, itb->second);
#define VECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=entry.namedV##name_t##s.cbegin(); itb!=entry.namedV##name_t##s.cend(); itb++) tree_->setVal(itb->first, &(itb->second));
#define DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=entry.namedVV##name_t##s.cbegin(); itb!=entry.namedVV##name_t##s.cend(); itb++) tree_->setVal(itb->first, &(itb->second));
    SIMPLE_DATA_OUTPUT_DIRECTIVES
    VECTOR_DATA_OUTPUT_DIRECTIVES
    DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVES
#undef SIMPLE_DATA_OUTPUT_DIRECTIVE
#undef VECTOR_DATA_OUTPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE
    tree_->fill();
  }

  // Save if this flag is specified
  if (robustSaveWrite){
    tree_->doAutoSave("FlushBaskets");
  }
}
