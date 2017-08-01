#ifndef BASETREE_H
#define BASETREE_H

#include "SampleHelpers.h"
#include "Samples.h"
#include "TSystem.h"
#include "TDirectory.h"


class BaseTree{
public:
  enum BranchType{
    BranchType_short_t,
    BranchType_uint_t,
    BranchType_int_t,
    BranchType_float_t,
    BranchType_double_t,

    BranchType_vshort_t,
    BranchType_vuint_t,
    BranchType_vint_t,
    BranchType_vfloat_t,
    BranchType_vdouble_t,

    BranchType_unknown_t
  };

protected:
  TFile* finput;
  TTree* tree;
  TTree* failedtree;
  TH1F* hCounters;
  bool valid;
  const bool receiver;

  unordered_map<TString, pair<short, short>*> valshorts;
  unordered_map<TString, pair<unsigned int, unsigned int>*> valuints;
  unordered_map<TString, pair<int, int>*> valints;
  unordered_map<TString, pair<float, float>*> valfloats;
  unordered_map<TString, pair<double, double>*> valdoubles;

  unordered_map<TString, vector<short>*> valVshorts;
  unordered_map<TString, vector<unsigned int>*> valVuints;
  unordered_map<TString, vector<int>*> valVints;
  unordered_map<TString, vector<float>*> valVfloats;
  unordered_map<TString, vector<double>*> valVdoubles;

  BranchType searchBranchType(TString branchname);

  template<BranchType T> void resetBranch();
  void resetBranches();

public:
  BaseTree();
  BaseTree(const TString cinput, const TString treename, const TString failedtreename, const TString countersname);
  ~BaseTree();

  template<typename T> bool bookBranch(TString branchname, T valdef);
  template<BranchType T> bool bookBranch(TString branchname);

  bool getSelectedEvent(int ev);
  bool getFailedEvent(int ev);

  int getSelectedNEvents(){ return (tree ? tree->GetEntries() : 0); }
  int getFailedNEvents(){ return (failedtree ? failedtree->GetEntries() : 0); }

  template<typename T> void getVal(TString branchname, T& val);
  template<typename T> void setVal(TString branchname, T const& val);

  void silenceUnused();

  bool isValid(){ return valid; }
  bool branchExists(TString branchname){ return (searchBranchType(branchname)!=BranchType_unknown_t); }

};

BaseTree::BaseTree() :
finput(nullptr),
tree(nullptr),
failedtree(nullptr),
hCounters(nullptr),
valid(false),
receiver(true)
{}
BaseTree::BaseTree(const TString cinput, const TString treename, const TString failedtreename, const TString countersname) :
finput(nullptr),
tree(nullptr),
failedtree(nullptr),
hCounters(nullptr),
valid(false),
receiver(true)
{
  TDirectory* curdir = gDirectory; // Save current directory to return back to it later
  if (!gSystem->AccessPathName(cinput)){
    finput = TFile::Open(cinput, "read");
    if (finput){
      if (finput->IsOpen() && !finput->IsZombie()){
        finput->cd();
        tree = (TTree*)finput->Get(treename);
        valid = (tree!=nullptr);
        if (!valid){ finput->Close(); finput=nullptr; }
        else{
          failedtree = (TTree*)finput->Get(failedtreename);
          hCounters = (TH1F*)finput->Get(countersname);
        }
      }
      else if (finput->IsOpen()){ finput->Close(); finput=nullptr; }
    }
  }
  curdir->cd(); // Return back to the directory before opening the input file
}
BaseTree::~BaseTree(){
  CalcHelpers::cleanUnorderedMap(valshorts);
  CalcHelpers::cleanUnorderedMap(valuints);
  CalcHelpers::cleanUnorderedMap(valints);
  CalcHelpers::cleanUnorderedMap(valfloats);
  CalcHelpers::cleanUnorderedMap(valdoubles);
  if (!receiver){
    CalcHelpers::cleanUnorderedMap(valVshorts);
    CalcHelpers::cleanUnorderedMap(valVuints);
    CalcHelpers::cleanUnorderedMap(valVints);
    CalcHelpers::cleanUnorderedMap(valVfloats);
    CalcHelpers::cleanUnorderedMap(valVdoubles);
  }
}

BaseTree::BranchType BaseTree::searchBranchType(TString branchname){
  if (valshorts.find(branchname)!=valshorts.end()) return BranchType_short_t;
  else if (valuints.find(branchname)!=valuints.end()) return BranchType_uint_t;
  else if (valints.find(branchname)!=valints.end()) return BranchType_int_t;
  else if (valfloats.find(branchname)!=valfloats.end()) return BranchType_float_t;
  else if (valdoubles.find(branchname)!=valdoubles.end()) return BranchType_double_t;

  else if (valVshorts.find(branchname)!=valVshorts.end()) return BranchType_vshort_t;
  else if (valVuints.find(branchname)!=valVuints.end()) return BranchType_vuint_t;
  else if (valVints.find(branchname)!=valVints.end()) return BranchType_vint_t;
  else if (valVfloats.find(branchname)!=valVfloats.end()) return BranchType_vfloat_t;
  else if (valVdoubles.find(branchname)!=valVdoubles.end()) return BranchType_vdouble_t;

  else return BranchType_unknown_t;
}

bool BaseTree::getSelectedEvent(int ev){
  resetBranches();
  if (tree && ev<tree->GetEntries()) return (tree->GetEntry(ev)>0);
  return false;
}
bool BaseTree::getFailedEvent(int ev){
  resetBranches();
  if (failedtree && ev<failedtree->GetEntries()) return (failedtree->GetEntry(ev)>0);
  return false;
}

template<> void BaseTree::resetBranch<BaseTree::BranchType_short_t>(){ for (auto& it:valshorts){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vshort_t>(){ for (auto& it:valVshorts){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_uint_t>(){ for (auto& it:valuints){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vuint_t>(){ for (auto& it:valVuints){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_int_t>(){ for (auto& it:valints){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vint_t>(){ for (auto& it:valVints){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_float_t>(){ for (auto& it:valfloats){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vfloat_t>(){ for (auto& it:valVfloats){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_double_t>(){ for (auto& it:valdoubles){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vdouble_t>(){ for (auto& it:valVdoubles){ if (it.second) it.second->clear(); } }

void BaseTree::resetBranches(){
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


template<> bool BaseTree::bookBranch<short>(TString branchname, short valdef){
  if (valshorts.find(branchname)==valshorts.end()) valshorts[branchname] = new pair<short, short>(valdef, valdef);
  else{ valshorts[branchname]->first=valdef; valshorts[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(valshorts[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valshorts[branchname]->first));
  return true;
}
template<> bool BaseTree::bookBranch<unsigned int>(TString branchname, unsigned int valdef){
  if (valuints.find(branchname)==valuints.end()) valuints[branchname] = new pair<unsigned int, unsigned int>(valdef, valdef);
  else{ valuints[branchname]->first=valdef; valuints[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(valuints[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valuints[branchname]->first));
  return true;
}
template<> bool BaseTree::bookBranch<int>(TString branchname, int valdef){
  if (valints.find(branchname)==valints.end()) valints[branchname] = new pair<int, int>(valdef, valdef);
  else{ valints[branchname]->first=valdef; valints[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(valints[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valints[branchname]->first));
  return true;
}
template<> bool BaseTree::bookBranch<float>(TString branchname, float valdef){
  if (valfloats.find(branchname)==valfloats.end()) valfloats[branchname] = new pair<float, float>(valdef, valdef);
  else{ valfloats[branchname]->first=valdef; valfloats[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(valfloats[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valfloats[branchname]->first));
  return true;
}
template<> bool BaseTree::bookBranch<double>(TString branchname, double valdef){
  if (valdoubles.find(branchname)==valdoubles.end()) valdoubles[branchname] = new pair<double, double>(valdef, valdef);
  else{ valdoubles[branchname]->first=valdef; valdoubles[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(valdoubles[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valdoubles[branchname]->first));
  return true;
}

template<> bool BaseTree::bookBranch<vector<short>*>(TString branchname, vector<short>*/* valdef*/){
  valVshorts[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVshorts[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVshorts[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<vector<unsigned int>*>(TString branchname, vector<unsigned int>*/* valdef*/){
  valVuints[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVuints[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVuints[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<vector<int>*>(TString branchname, vector<int>*/* valdef*/){
  valVints[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVints[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVints[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<vector<float>*>(TString branchname, vector<float>*/* valdef*/){
  valVfloats[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVfloats[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVfloats[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<vector<double>*>(TString branchname, vector<double>*/* valdef*/){
  valVdoubles[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVdoubles[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVdoubles[branchname]));
  return true;
}

template<> bool BaseTree::bookBranch<BaseTree::BranchType_short_t>(TString branchname){ return this->bookBranch<short>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_uint_t>(TString branchname){ return this->bookBranch<unsigned int>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_int_t>(TString branchname){ return this->bookBranch<int>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_float_t>(TString branchname){ return this->bookBranch<float>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_double_t>(TString branchname){ return this->bookBranch<double>(branchname, 0); }

template<> bool BaseTree::bookBranch<BaseTree::BranchType_vshort_t>(TString branchname){ return this->bookBranch<vector<short>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vuint_t>(TString branchname){ return this->bookBranch<vector<unsigned int>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vint_t>(TString branchname){ return this->bookBranch<vector<int>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vfloat_t>(TString branchname){ return this->bookBranch<vector<float>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vdouble_t>(TString branchname){ return this->bookBranch<vector<double>*>(branchname, 0); }


template<> void BaseTree::getVal<short>(TString branchname, short& val){
  if (searchBranchType(branchname)==BranchType_short_t){ auto& tmp = valshorts[branchname]; if (tmp) val=tmp->first; }
}
template<> void BaseTree::getVal<vector<short>*>(TString branchname, vector<short>*& val){
  if (searchBranchType(branchname)==BranchType_vshort_t) val = valVshorts[branchname];
}
template<> void BaseTree::getVal<uint>(TString branchname, unsigned int& val){
  if (searchBranchType(branchname)==BranchType_uint_t){ auto& tmp = valuints[branchname]; if (tmp) val=tmp->first; }
}
template<> void BaseTree::getVal<vector<unsigned int>*>(TString branchname, vector<unsigned int>*& val){
  if (searchBranchType(branchname)==BranchType_vuint_t) val = valVuints[branchname];
}
template<> void BaseTree::getVal<int>(TString branchname, int& val){
  if (searchBranchType(branchname)==BranchType_int_t){ auto& tmp = valints[branchname]; if (tmp) val=tmp->first; }
}
template<> void BaseTree::getVal<vector<int>*>(TString branchname, vector<int>*& val){
  if (searchBranchType(branchname)==BranchType_vint_t) val = valVints[branchname];
}
template<> void BaseTree::getVal<float>(TString branchname, float& val){
  if (searchBranchType(branchname)==BranchType_float_t){ auto& tmp = valfloats[branchname]; if (tmp) val=tmp->first; }
}
template<> void BaseTree::getVal<vector<float>*>(TString branchname, vector<float>*& val){
  if (searchBranchType(branchname)==BranchType_vfloat_t) val = valVfloats[branchname];
}
template<> void BaseTree::getVal<double>(TString branchname, double& val){
  if (searchBranchType(branchname)==BranchType_double_t){ auto& tmp = valdoubles[branchname]; if (tmp) val=tmp->first; }
}
template<> void BaseTree::getVal<vector<double>*>(TString branchname, vector<double>*& val){
  if (searchBranchType(branchname)==BranchType_vdouble_t) val = valVdoubles[branchname];
}

template<> void BaseTree::setVal<short>(TString branchname, short const& val){
  if (searchBranchType(branchname)==BranchType_short_t){ auto& tmp = valshorts[branchname]; if (tmp) tmp->first=val; }
}
template<> void BaseTree::setVal<uint>(TString branchname, unsigned int const& val){
  if (searchBranchType(branchname)==BranchType_uint_t){ auto& tmp = valuints[branchname]; if (tmp) tmp->first=val; }
}
template<> void BaseTree::setVal<int>(TString branchname, int const& val){
  if (searchBranchType(branchname)==BranchType_int_t){ auto& tmp = valints[branchname]; if (tmp) tmp->first=val; }
}
template<> void BaseTree::setVal<float>(TString branchname, float const& val){
  if (searchBranchType(branchname)==BranchType_float_t){ auto& tmp = valfloats[branchname]; if (tmp) tmp->first=val; }
}
template<> void BaseTree::setVal<double>(TString branchname, double const& val){
  if (searchBranchType(branchname)==BranchType_double_t){ auto& tmp = valdoubles[branchname]; if (tmp) tmp->first=val; }
}
template<> void BaseTree::setVal<vector<short>*>(TString branchname, vector<short>* const& val){
  if (searchBranchType(branchname)==BranchType_vshort_t){
    if (valVshorts[branchname] && val) valVshorts[branchname]->assign(val->begin(), val->end());
  }
}
template<> void BaseTree::setVal<vector<unsigned int>*>(TString branchname, vector<unsigned int>* const& val){
  if (searchBranchType(branchname)==BranchType_vuint_t){
    if (valVuints[branchname] && val) valVuints[branchname]->assign(val->begin(), val->end());
  }
}
template<> void BaseTree::setVal<vector<int>*>(TString branchname, vector<int>* const& val){
  if (searchBranchType(branchname)==BranchType_vint_t){
    if (valVints[branchname] && val) valVints[branchname]->assign(val->begin(), val->end());
  }
}
template<> void BaseTree::setVal<vector<float>*>(TString branchname, vector<float>* const& val){
  if (searchBranchType(branchname)==BranchType_vfloat_t){
    if (valVfloats[branchname] && val) valVfloats[branchname]->assign(val->begin(), val->end());
  }
}
template<> void BaseTree::setVal<vector<double>*>(TString branchname, vector<double>* const& val){
  if (searchBranchType(branchname)==BranchType_vdouble_t){
    if (valVdoubles[branchname] && val) valVdoubles[branchname]->assign(val->begin(), val->end());
  }
}


void BaseTree::silenceUnused(){
  const unsigned int ntrees = 2;
  TTree* trees[ntrees]={
    tree,
    failedtree
  };
  for (unsigned int it=0; it<ntrees; it++){
    const TList* blist = (const TList*)trees[it]->GetListOfBranches();
    for (int ib=0; ib<blist->GetSize(); ib++){
      TString bname = blist->At(ib)->GetName();
      if (searchBranchType(bname)==BranchType_unknown_t) trees[it]->SetBranchStatus(bname, 0);
    }
  }

}

#endif
