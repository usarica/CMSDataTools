#ifndef IVYBASE_H
#define IVYBASE_H

#include "TVar.hh"
#include "BaseTree.h"
#include "HelperFunctionsCore.h"
#include "MELAStreamHelpers.hh"


class IvyBase{
protected:
  TVar::VerbosityLevel verbosity;
  BaseTree* currentTree;

  // Consumes
  std::vector<TString> sloppyConsumes; // In case some variables are known to be absent in some trees

#define SIMPLE_DATA_INPUT_DIRECTIVE(name, type, default_value) std::unordered_map<TString, type*> val##name##s;
#define VECTOR_DATA_INPUT_DIRECTIVE(name, type) std::unordered_map<TString, type**> valV##name##s;
#define DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(name, type) std::unordered_map<TString, type**> valVV##name##s;
  SIMPLE_DATA_INPUT_DIRECTIVES
  VECTOR_DATA_INPUT_DIRECTIVES
  DOUBLEVECTOR_DATA_INPUT_DIRECTIVES
#undef SIMPLE_DATA_INPUT_DIRECTIVE
#undef VECTOR_DATA_INPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_INPUT_DIRECTIVE

  template<typename T> bool linkConsumed(BaseTree* tree);
  bool linkConsumes(BaseTree* tree);

  // Get consumed map
  template<typename T> void getConsumedMap(std::unordered_map<TString, T*>*& theMap);
  template<typename T> void getConsumedMap(std::unordered_map<TString, T*> const*& theMap) const;

  // Get consumed reference
  template<typename T> bool getConsumed(TString name, T*& val) const;
  template<typename T> bool getConsumed(TString name, T const*& val) const;

  // Get consumed value
  template<typename T> bool getConsumedValue(TString name, T& val) const;

public:
  // Constructors
  IvyBase();

  // Destructors
  virtual ~IvyBase();

  // Add the necessary objects
  template<typename T> void addConsumed(TString name);
  void defineConsumedSloppy(TString name);

  // Verbosity
  void setVerbosity(TVar::VerbosityLevel flag){ verbosity=flag; }
  TVar::VerbosityLevel getVerbosity() const{ return verbosity; }

  // Tree
  virtual bool wrapTree(BaseTree* tree);
  BaseTree* getWrappedTree(){ return currentTree; }

};

template<typename T> bool IvyBase::getConsumedValue(TString name, T& val) const{
  T const* ref;
  bool res = this->getConsumed<T>(name, ref);
  if (res && ref!=nullptr){
    val=*ref;
    return true;
  }
  else{
    if (!res && verbosity>=TVar::ERROR) MELAStreamHelpers::MELAerr << "IvyBase::getConsumedValue: Cannot consume " << name << std::endl;
    return res;
  }
}


#endif
