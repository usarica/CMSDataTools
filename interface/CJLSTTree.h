#ifndef CJLSTTREE_H
#define CJLSTTREE_H

#include "SampleHelpers.h"
#include "Samples.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TList.h"


class CJLSTTree{
protected:
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
  CJLSTTree(TString strsample); // This is the constructor to receive CJLST trees
  ~CJLSTTree();

  template<typename T> bool bookBranch(TString branchname, T valdef);
  template<BranchType T> bool bookBranch(TString branchname);

  bool getSelectedEvent(int ev);
  bool getFailedEvent(int ev);

  template<typename T> void getVal(TString branchname, T& val);
  template<typename T> void setVal(TString branchname, T const& val);

  void silenceUnused();

};


CJLSTTree::CJLSTTree(TString strsample) :
finput(nullptr),
tree(nullptr),
failedtree(nullptr),
hCounters(nullptr),
valid(false),
receiver(true)
{
  const TString TREE_NAME = "ZZTree/candTree";
  const TString TREE_FAILED_NAME = "ZZTree/candTree_failed";
  const TString COUNTERS_NAME = "ZZTree/Counters";

  TDirectory* curdir = gDirectory; // Save current directory to return back to it later
  TString cinput = CJLSTsamplesdir + "/" + strsample + "/ZZ4lAnalysis.root";
  if (!gSystem->AccessPathName(cinput)){
    finput = TFile::Open(cinput, "read");
    if (finput){
      if (finput->IsOpen() && !finput->IsZombie()){
        finput->cd();
        tree = (TTree*)finput->Get(TREE_NAME);
        valid = (tree!=nullptr);
        if (!valid){ finput->Close(); finput=nullptr; }
        else{
          failedtree = (TTree*)finput->Get(TREE_FAILED_NAME);
          hCounters = (TH1F*)finput->Get(COUNTERS_NAME);
        }
      }
      else if (finput->IsOpen()){ finput->Close(); finput=nullptr; }
    }
  }
  curdir->cd(); // Return back to the directory before opening the input file
}
CJLSTTree::~CJLSTTree(){
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

CJLSTTree::BranchType CJLSTTree::searchBranchType(TString branchname){
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

bool CJLSTTree::getSelectedEvent(int ev){
  resetBranches();
  if (tree && ev<tree->GetEntries()) return (tree->GetEntry(ev)>0);
  return false;
}
bool CJLSTTree::getFailedEvent(int ev){
  resetBranches();
  if (failedtree && ev<failedtree->GetEntries()) return (failedtree->GetEntry(ev)>0);
  return false;
}

template<> void CJLSTTree::resetBranch<CJLSTTree::BranchType_short_t>(){ for (auto& it:valshorts){ if (it.second){ it.second->first=it.second->second; } } }
template<> void CJLSTTree::resetBranch<CJLSTTree::BranchType_vshort_t>(){ for (auto& it:valVshorts){ if (it.second) it.second->clear(); } }
template<> void CJLSTTree::resetBranch<CJLSTTree::BranchType_uint_t>(){ for (auto& it:valuints){ if (it.second){ it.second->first=it.second->second; } } }
template<> void CJLSTTree::resetBranch<CJLSTTree::BranchType_vuint_t>(){ for (auto& it:valVuints){ if (it.second) it.second->clear(); } }
template<> void CJLSTTree::resetBranch<CJLSTTree::BranchType_int_t>(){ for (auto& it:valints){ if (it.second){ it.second->first=it.second->second; } } }
template<> void CJLSTTree::resetBranch<CJLSTTree::BranchType_vint_t>(){ for (auto& it:valVints){ if (it.second) it.second->clear(); } }
template<> void CJLSTTree::resetBranch<CJLSTTree::BranchType_float_t>(){ for (auto& it:valfloats){ if (it.second){ it.second->first=it.second->second; } } }
template<> void CJLSTTree::resetBranch<CJLSTTree::BranchType_vfloat_t>(){ for (auto& it:valVfloats){ if (it.second) it.second->clear(); } }
template<> void CJLSTTree::resetBranch<CJLSTTree::BranchType_double_t>(){ for (auto& it:valdoubles){ if (it.second){ it.second->first=it.second->second; } } }
template<> void CJLSTTree::resetBranch<CJLSTTree::BranchType_vdouble_t>(){ for (auto& it:valVdoubles){ if (it.second) it.second->clear(); } }

void CJLSTTree::resetBranches(){
  this->resetBranch<CJLSTTree::BranchType_short_t>();
  this->resetBranch<CJLSTTree::BranchType_uint_t>();
  this->resetBranch<CJLSTTree::BranchType_int_t>();
  this->resetBranch<CJLSTTree::BranchType_float_t>();
  this->resetBranch<CJLSTTree::BranchType_double_t>();

  this->resetBranch<CJLSTTree::BranchType_vshort_t>();
  this->resetBranch<CJLSTTree::BranchType_vuint_t>();
  this->resetBranch<CJLSTTree::BranchType_vint_t>();
  this->resetBranch<CJLSTTree::BranchType_vfloat_t>();
  this->resetBranch<CJLSTTree::BranchType_vdouble_t>();
}


template<> bool CJLSTTree::bookBranch<short>(TString branchname, short valdef){
  valshorts[branchname] = new pair<short,short>(valdef,valdef);
  SampleHelpers::bookBranch(tree, branchname, &(valshorts[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valshorts[branchname]->first));
  return true;
}
template<> bool CJLSTTree::bookBranch<unsigned int>(TString branchname, unsigned int valdef){
  valuints[branchname] = new pair<unsigned int, unsigned int>(valdef, valdef);
  SampleHelpers::bookBranch(tree, branchname, &(valuints[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valuints[branchname]->first));
  return true;
}
template<> bool CJLSTTree::bookBranch<int>(TString branchname, int valdef){
  valints[branchname] = new pair<int, int>(valdef, valdef);
  SampleHelpers::bookBranch(tree, branchname, &(valints[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valints[branchname]->first));
  return true;
}
template<> bool CJLSTTree::bookBranch<float>(TString branchname, float valdef){
  valfloats[branchname] = new pair<float, float>(valdef, valdef);
  SampleHelpers::bookBranch(tree, branchname, &(valfloats[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valfloats[branchname]->first));
  return true;
}
template<> bool CJLSTTree::bookBranch<double>(TString branchname, double valdef){
  valdoubles[branchname] = new pair<double, double>(valdef, valdef);
  SampleHelpers::bookBranch(tree, branchname, &(valdoubles[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valdoubles[branchname]->first));
  return true;
}

template<> bool CJLSTTree::bookBranch<vector<short>*>(TString branchname, vector<short>*/* valdef*/){
  valVshorts[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVshorts[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVshorts[branchname]));
  return true;
}
template<> bool CJLSTTree::bookBranch<vector<unsigned int>*>(TString branchname, vector<unsigned int>*/* valdef*/){
  valVuints[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVuints[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVuints[branchname]));
  return true;
}
template<> bool CJLSTTree::bookBranch<vector<int>*>(TString branchname, vector<int>*/* valdef*/){
  valVints[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVints[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVints[branchname]));
  return true;
}
template<> bool CJLSTTree::bookBranch<vector<float>*>(TString branchname, vector<float>*/* valdef*/){
  valVfloats[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVfloats[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVfloats[branchname]));
  return true;
}
template<> bool CJLSTTree::bookBranch<vector<double>*>(TString branchname, vector<double>*/* valdef*/){
  valVdoubles[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVdoubles[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVdoubles[branchname]));
  return true;
}

template<> bool CJLSTTree::bookBranch<CJLSTTree::BranchType_short_t>(TString branchname){ return this->bookBranch<short>(branchname, 0); }
template<> bool CJLSTTree::bookBranch<CJLSTTree::BranchType_uint_t>(TString branchname){ return this->bookBranch<unsigned int>(branchname, 0); }
template<> bool CJLSTTree::bookBranch<CJLSTTree::BranchType_int_t>(TString branchname){ return this->bookBranch<int>(branchname, 0); }
template<> bool CJLSTTree::bookBranch<CJLSTTree::BranchType_float_t>(TString branchname){ return this->bookBranch<float>(branchname, 0); }
template<> bool CJLSTTree::bookBranch<CJLSTTree::BranchType_double_t>(TString branchname){ return this->bookBranch<double>(branchname, 0); }

template<> bool CJLSTTree::bookBranch<CJLSTTree::BranchType_vshort_t>(TString branchname){ return this->bookBranch<vector<short>*>(branchname, 0); }
template<> bool CJLSTTree::bookBranch<CJLSTTree::BranchType_vuint_t>(TString branchname){ return this->bookBranch<vector<unsigned int>*>(branchname, 0); }
template<> bool CJLSTTree::bookBranch<CJLSTTree::BranchType_vint_t>(TString branchname){ return this->bookBranch<vector<int>*>(branchname, 0); }
template<> bool CJLSTTree::bookBranch<CJLSTTree::BranchType_vfloat_t>(TString branchname){ return this->bookBranch<vector<float>*>(branchname, 0); }
template<> bool CJLSTTree::bookBranch<CJLSTTree::BranchType_vdouble_t>(TString branchname){ return this->bookBranch<vector<double>*>(branchname, 0); }


template<> void CJLSTTree::getVal<short>(TString branchname, short& val){
  if (searchBranchType(branchname)==BranchType_short_t){ auto& tmp = valshorts[branchname]; if (tmp) val=tmp->first; }
}
template<> void CJLSTTree::getVal<vector<short>*>(TString branchname, vector<short>*& val){
  if (searchBranchType(branchname)==BranchType_vshort_t) val = valVshorts[branchname];
}
template<> void CJLSTTree::getVal<uint>(TString branchname, unsigned int& val){
  if (searchBranchType(branchname)==BranchType_uint_t){ auto& tmp = valuints[branchname]; if (tmp) val=tmp->first; }
}
template<> void CJLSTTree::getVal<vector<unsigned int>*>(TString branchname, vector<unsigned int>*& val){
  if (searchBranchType(branchname)==BranchType_vuint_t) val = valVuints[branchname];
}
template<> void CJLSTTree::getVal<int>(TString branchname, int& val){
  if (searchBranchType(branchname)==BranchType_int_t){ auto& tmp = valints[branchname]; if (tmp) val=tmp->first; }
}
template<> void CJLSTTree::getVal<vector<int>*>(TString branchname, vector<int>*& val){
  if (searchBranchType(branchname)==BranchType_vint_t) val = valVints[branchname];
}
template<> void CJLSTTree::getVal<float>(TString branchname, float& val){
  if (searchBranchType(branchname)==BranchType_float_t){ auto& tmp = valfloats[branchname]; if (tmp) val=tmp->first; }
}
template<> void CJLSTTree::getVal<vector<float>*>(TString branchname, vector<float>*& val){
  if (searchBranchType(branchname)==BranchType_vfloat_t) val = valVfloats[branchname];
}
template<> void CJLSTTree::getVal<double>(TString branchname, double& val){
  if (searchBranchType(branchname)==BranchType_double_t){ auto& tmp = valdoubles[branchname]; if (tmp) val=tmp->first; }
}
template<> void CJLSTTree::getVal<vector<double>*>(TString branchname, vector<double>*& val){
  if (searchBranchType(branchname)==BranchType_vdouble_t) val = valVdoubles[branchname];
}

template<> void CJLSTTree::setVal<short>(TString branchname, short const& val){
  if (searchBranchType(branchname)==BranchType_short_t){ auto& tmp = valshorts[branchname]; if (tmp) tmp->first=val; }
}
template<> void CJLSTTree::setVal<uint>(TString branchname, unsigned int const& val){
  if (searchBranchType(branchname)==BranchType_uint_t){ auto& tmp = valuints[branchname]; if (tmp) tmp->first=val; }
}
template<> void CJLSTTree::setVal<int>(TString branchname, int const& val){
  if (searchBranchType(branchname)==BranchType_int_t){ auto& tmp = valints[branchname]; if (tmp) tmp->first=val; }
}
template<> void CJLSTTree::setVal<float>(TString branchname, float const& val){
  if (searchBranchType(branchname)==BranchType_float_t){ auto& tmp = valfloats[branchname]; if (tmp) tmp->first=val; }
}
template<> void CJLSTTree::setVal<double>(TString branchname, double const& val){
  if (searchBranchType(branchname)==BranchType_double_t){ auto& tmp = valdoubles[branchname]; if (tmp) tmp->first=val; }
}
template<> void CJLSTTree::setVal<vector<short>*>(TString branchname, vector<short>* const& val){
  if (searchBranchType(branchname)==BranchType_vshort_t){
    if (valVshorts[branchname] && val) valVshorts[branchname]->assign(val->begin(), val->end());
  }
}
template<> void CJLSTTree::setVal<vector<unsigned int>*>(TString branchname, vector<unsigned int>* const& val){
  if (searchBranchType(branchname)==BranchType_vuint_t){
    if (valVuints[branchname] && val) valVuints[branchname]->assign(val->begin(), val->end());
  }
}
template<> void CJLSTTree::setVal<vector<int>*>(TString branchname, vector<int>* const& val){
  if (searchBranchType(branchname)==BranchType_vint_t){
    if (valVints[branchname] && val) valVints[branchname]->assign(val->begin(), val->end());
  }
}
template<> void CJLSTTree::setVal<vector<float>*>(TString branchname, vector<float>* const& val){
  if (searchBranchType(branchname)==BranchType_vfloat_t){
    if (valVfloats[branchname] && val) valVfloats[branchname]->assign(val->begin(), val->end());
  }
}
template<> void CJLSTTree::setVal<vector<double>*>(TString branchname, vector<double>* const& val){
  if (searchBranchType(branchname)==BranchType_vdouble_t){
    if (valVdoubles[branchname] && val) valVdoubles[branchname]->assign(val->begin(), val->end());
  }
}


void CJLSTTree::silenceUnused(){
  const unsigned int nlists = 2;
  TTree* trees[nlists]={
    tree,
    failedtree
  };
  TList* lists[nlists]={
    (TList*)tree->GetListOfBranches(),
    (TList*)failedtree->GetListOfBranches()
  };
  for (unsigned int it=0; it<nlists; it++){
    for (int ib=0; ib<lists[it]->GetSize(); ib++){
      TString bname = lists[it]->At(ib)->GetName();
      if (searchBranchType(bname)==BranchType_unknown_t) trees[it]->SetBranchStatus(bname, 0);
    }
  }

}

#endif
