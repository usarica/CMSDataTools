#ifndef IVYBASE_HPP
#define IVYBASE_HPP

#include "IvyBase.h"


using namespace std;
using namespace MELAStreamHelpers;


template<typename T> void IvyBase::addConsumed(TString name){
  std::unordered_map<TString, T*>* theMap = nullptr;
  IvyBase::getConsumedMap<T>(theMap);
  if (theMap) (*theMap)[name] = nullptr;
  else if (verbosity>=TVar::ERROR) MELAerr << "IvyBase::addConsumed(" << name << "): Map could not be found." << endl;
}
template<typename T> bool IvyBase::linkConsumed(BaseTree* tree){
  bool result=true;
  std::unordered_map<TString, T*>* theMap = nullptr;
  IvyBase::getConsumedMap<T>(theMap);
  if (theMap){
    for (typename std::unordered_map<TString, T*>::iterator it=theMap->begin(); it!=theMap->end(); it++){
      if (tree->branchExists(it->first)){
        tree->getValRef<T>(it->first, it->second);
        result &= true;
        if (verbosity>=TVar::INFO) MELAout << "IvyBase::linkConsumed(" << tree->sampleIdentifier << "): Linking successful for variable " << it->first << " -> " << it->second << endl;
      }
      else if (std::find(this->sloppyConsumes.begin(), this->sloppyConsumes.end(), it->first)!=this->sloppyConsumes.end()){
        it->second=nullptr;
        result &= true;
        if (verbosity>=TVar::INFO) MELAout << "IvyBase::linkConsumed(" << tree->sampleIdentifier << "): Linking failed for variable " << it->first << ", but the variable is sloppy." << endl;
      }
      else{
        result &= false;
        if (verbosity>=TVar::ERROR) MELAerr << "IvyBase::linkConsumed(" << tree->sampleIdentifier << "): Linking failed for variable " << it->first << endl;
      }
    }
  }
  else{
    result = false;
    if (verbosity>=TVar::ERROR) MELAerr << "IvyBase::linkConsumed(" << tree->sampleIdentifier << "): Map could not be found." << endl;
  }
  return result;
}
template<typename T> bool IvyBase::getConsumed(TString name, T*& val) const{
  val = nullptr; // Reset value
  bool isSloppy = (std::find(this->sloppyConsumes.begin(), this->sloppyConsumes.end(), name)!=this->sloppyConsumes.end());
  std::unordered_map<TString, T*> const* theMap = nullptr;
  IvyBase::getConsumedMap<T>(theMap);
  if (theMap){
    typename unordered_map<TString, T*>::const_iterator it_val;
    if (HelperFunctions::getUnorderedMapIterator(name, *theMap, it_val)){
      val = it_val->second;
      return true;
    }
    else if (isSloppy) return true;
  }
  else if (verbosity>=TVar::ERROR) MELAerr << "IvyBase::getConsumed(" << name << "): Map could not be found." << endl;
  return false;
}
template<typename T> bool IvyBase::getConsumed(TString name, T const*& val) const{
  val = nullptr; // Reset value
  bool isSloppy = (std::find(this->sloppyConsumes.begin(), this->sloppyConsumes.end(), name)!=this->sloppyConsumes.end());
  std::unordered_map<TString, T*> const* theMap = nullptr;
  IvyBase::getConsumedMap<T>(theMap);
  if (theMap){
    typename unordered_map<TString, T*>::const_iterator it_val;
    if (HelperFunctions::getUnorderedMapIterator(name, *theMap, it_val)){
      val = it_val->second;
      return true;
    }
    else if (isSloppy) return true;
  }
  else if (verbosity>=TVar::ERROR) MELAerr << "IvyBase::getConsumed(" << name << "): Map could not be found." << endl;
  return false;
}


#define SIMPLE_DATA_INPUT_DIRECTIVE(name, type, default_value) \
template<> void IvyBase::getConsumedMap<type>(std::unordered_map<TString, type*>*& theMap){ theMap = &val##name##s; } \
template<> void IvyBase::getConsumedMap<type>(std::unordered_map<TString, type*> const*& theMap) const{ theMap = &val##name##s; } \
template void IvyBase::addConsumed<type>(TString); \
template bool IvyBase::linkConsumed<type>(BaseTree* tree); \
template bool IvyBase::getConsumed<type>(TString, type*& val) const; \
template bool IvyBase::getConsumed<type>(TString, type const*& val) const; \


#define VECTOR_DATA_INPUT_DIRECTIVE(name, type) \
template<> void IvyBase::getConsumedMap<type*>(std::unordered_map<TString, type**>*& theMap){ theMap = &valV##name##s; } \
template<> void IvyBase::getConsumedMap<type*>(std::unordered_map<TString, type**> const*& theMap) const{ theMap = &valV##name##s; } \
template void IvyBase::addConsumed<type*>(TString); \
template bool IvyBase::linkConsumed<type*>(BaseTree* tree); \
template bool IvyBase::getConsumed<type*>(TString, type**& val) const; \
template bool IvyBase::getConsumed<type*>(TString, type* const*& val) const; \


#define DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(name, type) \
template<> void IvyBase::getConsumedMap<type*>(std::unordered_map<TString, type**>*& theMap){ theMap = &valVV##name##s; } \
template<> void IvyBase::getConsumedMap<type*>(std::unordered_map<TString, type**> const*& theMap) const{ theMap = &valVV##name##s; } \
template void IvyBase::addConsumed<type*>(TString); \
template bool IvyBase::linkConsumed<type*>(BaseTree* tree); \
template bool IvyBase::getConsumed<type*>(TString, type**& val) const; \
template bool IvyBase::getConsumed<type*>(TString, type* const*& val) const; \


SIMPLE_DATA_INPUT_DIRECTIVES
VECTOR_DATA_INPUT_DIRECTIVES
DOUBLEVECTOR_DATA_INPUT_DIRECTIVES


#undef SIMPLE_DATA_INPUT_DIRECTIVE
#undef VECTOR_DATA_INPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_INPUT_DIRECTIVE


#endif
