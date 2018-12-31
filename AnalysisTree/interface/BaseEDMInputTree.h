#ifndef BASEEDMINPUTTREE_H
#define BASEEDMINPUTTREE_H

#include "BaseTree.h"
#include "CMSEDMWrapper.h"
#include "MELAStreamHelpers.hh"


template<typename T, typename U=T> class CMSEDMWrapperLinker{
protected:
  typedef T Wrapped_t;
  typedef U Val_t;
  typedef edm::Wrapper<Wrapped_t> Wrapper_t;

  Wrapper_t* var;
  Val_t* targetVal;

  void assignProductToTarget(Wrapped_t& product);

public:
  CMSEDMWrapperLinker();
  CMSEDMWrapperLinker(Val_t* targetVal_);
  CMSEDMWrapperLinker(CMSEDMWrapperLinker<T, U> const&) = delete;

  void synchronize();
  void reset();

  void print() const;

  Wrapper_t*& getWrapperRef(){ return var; }
  Val_t*& getTargetRef(){ return targetVal; }

  Wrapper_t const* getWrapperRef() const{ return var; }
  Val_t const* getTargetRef() const{ return targetVal; }

};

template<typename T, typename U> CMSEDMWrapperLinker<T, U>::CMSEDMWrapperLinker() :
  var(nullptr),
  targetVal(nullptr)
{}
template<typename T, typename U> CMSEDMWrapperLinker<T, U>::CMSEDMWrapperLinker(Val_t* targetVal_) :
  var(nullptr),
  targetVal(targetVal_)
{}

template<typename T, typename U> void CMSEDMWrapperLinker<T, U>::synchronize(){
  if (!var) return;
  if (var->isPresent()){
    Wrapped_t& product = var->bareProduct();
    this->assignProductToTarget(product);
  }
}
template<typename T, typename U> void CMSEDMWrapperLinker<T, U>::reset(){ var = nullptr; }


class BaseEDMInputTree : public BaseTree{
protected:
  std::unordered_map<TString, CMSEDMWrapperLinker<bool>*> bridgebools;
  std::unordered_map<TString, CMSEDMWrapperLinker<short>*> bridgeshorts;
  std::unordered_map<TString, CMSEDMWrapperLinker<unsigned int>*> bridgeuints;
  std::unordered_map<TString, CMSEDMWrapperLinker<int>*> bridgeints;
  std::unordered_map<TString, CMSEDMWrapperLinker<unsigned long>*> bridgeulongs;
  std::unordered_map<TString, CMSEDMWrapperLinker<long>*> bridgelongs;
  std::unordered_map<TString, CMSEDMWrapperLinker<unsigned long long>*> bridgeulonglongs;
  std::unordered_map<TString, CMSEDMWrapperLinker<long long>*> bridgelonglongs;
  std::unordered_map<TString, CMSEDMWrapperLinker<float>*> bridgefloats;
  std::unordered_map<TString, CMSEDMWrapperLinker<double>*> bridgedoubles;
  std::unordered_map<TString, CMSEDMWrapperLinker<std::string>*> bridgestrings;
  std::unordered_map<TString, CMSEDMWrapperLinker<CMSLorentzVector>*> bridgeCMSLorentzVectors;

  std::unordered_map<TString, CMSEDMWrapperLinker<std::vector<bool>, std::vector<bool>*>*> bridgeVbools;
  std::unordered_map<TString, CMSEDMWrapperLinker<std::vector<short>, std::vector<short>*>*> bridgeVshorts;
  std::unordered_map<TString, CMSEDMWrapperLinker<std::vector<unsigned int>, std::vector<unsigned int>*>*> bridgeVuints;
  std::unordered_map<TString, CMSEDMWrapperLinker<std::vector<int>, std::vector<int>*>*> bridgeVints;
  std::unordered_map<TString, CMSEDMWrapperLinker<std::vector<unsigned long>, std::vector<unsigned long>*>*> bridgeVulongs;
  std::unordered_map<TString, CMSEDMWrapperLinker<std::vector<long>, std::vector<long>*>*> bridgeVlongs;
  std::unordered_map<TString, CMSEDMWrapperLinker<std::vector<unsigned long long>, std::vector<unsigned long long>*>*> bridgeVulonglongs;
  std::unordered_map<TString, CMSEDMWrapperLinker<std::vector<long long>, std::vector<long long>*>*> bridgeVlonglongs;
  std::unordered_map<TString, CMSEDMWrapperLinker<std::vector<float>, std::vector<float>*>*> bridgeVfloats;
  std::unordered_map<TString, CMSEDMWrapperLinker<std::vector<double>, std::vector<double>*>*> bridgeVdoubles;
  std::unordered_map<TString, CMSEDMWrapperLinker<std::vector<std::string>, std::vector<std::string>*>*> bridgeVstrings;
  std::unordered_map<TString, CMSEDMWrapperLinker<std::vector<CMSLorentzVector>, std::vector<CMSLorentzVector>*>*> bridgeVCMSLorentzVectors;

  template<BranchType T> void resetEDMBranch();
  void resetBranches();
  template<BranchType T> void removeEDMBranch(TString branchname);

  void synchronizeEDMBranches();

public:
  BaseEDMInputTree();
  BaseEDMInputTree(const TString cinput, const TString treename, const TString failedtreename, const TString countersname);
  virtual ~BaseEDMInputTree();

  template<typename T> bool bookEDMBranch(TString branchname, T valdef);
  template<BranchType T> bool bookEDMBranch(TString branchname);

  bool getSelectedEvent(int ev);
  bool getFailedEvent(int ev);

  void print() const;

};


#endif
