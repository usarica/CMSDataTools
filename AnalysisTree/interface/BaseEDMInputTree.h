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
#define SIMPLE_DATA_INPUT_DIRECTIVE(name, type, default_value) std::unordered_map<TString, CMSEDMWrapperLinker<type>*> bridge##name##s;
#define VECTOR_DATA_INPUT_DIRECTIVE(name, type) std::unordered_map<TString, CMSEDMWrapperLinker<type, type*>*> bridgeV##name##s;
#define DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(name, type) std::unordered_map<TString, CMSEDMWrapperLinker<type, type*>*> bridgeVV##name##s;
  SIMPLE_DATA_INPUT_DIRECTIVES
  VECTOR_DATA_INPUT_DIRECTIVES
  DOUBLEVECTOR_DATA_INPUT_DIRECTIVES
#undef SIMPLE_DATA_INPUT_DIRECTIVE
#undef VECTOR_DATA_INPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_INPUT_DIRECTIVE

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
