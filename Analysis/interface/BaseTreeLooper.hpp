#ifndef BASETREELOOPER_HPP
#define BASETREELOOPER_HPP

#include "BaseTreeLooper.h"
#include "MELAStreamHelpers.hh"


using namespace MELAStreamHelpers;


template<> void BaseTreeLooper::getConsumedMap<short>(std::unordered_map<TString, short*>*& theMap){ theMap = &valshorts; }
template<> void BaseTreeLooper::getConsumedMap<unsigned int>(std::unordered_map<TString, unsigned int*>*& theMap){ theMap = &valuints; }
template<> void BaseTreeLooper::getConsumedMap<int>(std::unordered_map<TString, int*>*& theMap){ theMap = &valints; }
template<> void BaseTreeLooper::getConsumedMap<unsigned long>(std::unordered_map<TString, unsigned long*>*& theMap){ theMap = &valulongs; }
template<> void BaseTreeLooper::getConsumedMap<long>(std::unordered_map<TString, long*>*& theMap){ theMap = &vallongs; }
template<> void BaseTreeLooper::getConsumedMap<long long>(std::unordered_map<TString, long long*>*& theMap){ theMap = &vallonglongs; }
template<> void BaseTreeLooper::getConsumedMap<float>(std::unordered_map<TString, float*>*& theMap){ theMap = &valfloats; }
template<> void BaseTreeLooper::getConsumedMap<double>(std::unordered_map<TString, double*>*& theMap){ theMap = &valdoubles; }
template<> void BaseTreeLooper::getConsumedMap<std::vector<short>>(std::unordered_map<TString, std::vector<short>*>*& theMap){ theMap = &valVshorts; }
template<> void BaseTreeLooper::getConsumedMap<std::vector<unsigned int>>(std::unordered_map<TString, std::vector<unsigned int>*>*& theMap){ theMap = &valVuints; }
template<> void BaseTreeLooper::getConsumedMap<std::vector<int>>(std::unordered_map<TString, std::vector<int>*>*& theMap){ theMap = &valVints; }
template<> void BaseTreeLooper::getConsumedMap<std::vector<unsigned long>>(std::unordered_map<TString, std::vector<unsigned long>*>*& theMap){ theMap = &valVulongs; }
template<> void BaseTreeLooper::getConsumedMap<std::vector<long>>(std::unordered_map<TString, std::vector<long>*>*& theMap){ theMap = &valVlongs; }
template<> void BaseTreeLooper::getConsumedMap<std::vector<long long>>(std::unordered_map<TString, std::vector<long long>*>*& theMap){ theMap = &valVlonglongs; }
template<> void BaseTreeLooper::getConsumedMap<std::vector<float>>(std::unordered_map<TString, std::vector<float>*>*& theMap){ theMap = &valVfloats; }
template<> void BaseTreeLooper::getConsumedMap<std::vector<double>>(std::unordered_map<TString, std::vector<double>*>*& theMap){ theMap = &valVdoubles; }

template<> void BaseTreeLooper::getConsumedMap<short>(std::unordered_map<TString, short*> const*& theMap) const{ theMap = &valshorts; }
template<> void BaseTreeLooper::getConsumedMap<unsigned int>(std::unordered_map<TString, unsigned int*> const*& theMap) const{ theMap = &valuints; }
template<> void BaseTreeLooper::getConsumedMap<int>(std::unordered_map<TString, int*> const*& theMap) const{ theMap = &valints; }
template<> void BaseTreeLooper::getConsumedMap<unsigned long>(std::unordered_map<TString, unsigned long*> const*& theMap) const{ theMap = &valulongs; }
template<> void BaseTreeLooper::getConsumedMap<long>(std::unordered_map<TString, long*> const*& theMap) const{ theMap = &vallongs; }
template<> void BaseTreeLooper::getConsumedMap<long long>(std::unordered_map<TString, long long*> const*& theMap) const{ theMap = &vallonglongs; }
template<> void BaseTreeLooper::getConsumedMap<float>(std::unordered_map<TString, float*> const*& theMap) const{ theMap = &valfloats; }
template<> void BaseTreeLooper::getConsumedMap<double>(std::unordered_map<TString, double*> const*& theMap) const{ theMap = &valdoubles; }
template<> void BaseTreeLooper::getConsumedMap<std::vector<short>>(std::unordered_map<TString, std::vector<short>*> const*& theMap) const{ theMap = &valVshorts; }
template<> void BaseTreeLooper::getConsumedMap<std::vector<unsigned int>>(std::unordered_map<TString, std::vector<unsigned int>*> const*& theMap) const{ theMap = &valVuints; }
template<> void BaseTreeLooper::getConsumedMap<std::vector<int>>(std::unordered_map<TString, std::vector<int>*> const*& theMap) const{ theMap = &valVints; }
template<> void BaseTreeLooper::getConsumedMap<std::vector<unsigned long>>(std::unordered_map<TString, std::vector<unsigned long>*> const*& theMap) const{ theMap = &valVulongs; }
template<> void BaseTreeLooper::getConsumedMap<std::vector<long>>(std::unordered_map<TString, std::vector<long>*> const*& theMap) const{ theMap = &valVlongs; }
template<> void BaseTreeLooper::getConsumedMap<std::vector<long long>>(std::unordered_map<TString, std::vector<long long>*> const*& theMap) const{ theMap = &valVlonglongs; }
template<> void BaseTreeLooper::getConsumedMap<std::vector<float>>(std::unordered_map<TString, std::vector<float>*> const*& theMap) const{ theMap = &valVfloats; }
template<> void BaseTreeLooper::getConsumedMap<std::vector<double>>(std::unordered_map<TString, std::vector<double>*> const*& theMap) const{ theMap = &valVdoubles; }

template<typename T> void BaseTreeLooper::addConsumed(TString name){
  std::unordered_map<TString, T*>* theMap = nullptr;
  BaseTreeLooper::getConsumedMap<T>(theMap);
  if (theMap) (*theMap)[name] = nullptr;
  else MELAerr << "BaseTreeLooper::addConsumed(" << name << "): Map could not be found." << endl;
}
template void BaseTreeLooper::addConsumed<short>(TString name);
template void BaseTreeLooper::addConsumed<unsigned int>(TString name);
template void BaseTreeLooper::addConsumed<int>(TString name);
template void BaseTreeLooper::addConsumed<unsigned long>(TString name);
template void BaseTreeLooper::addConsumed<long>(TString name);
template void BaseTreeLooper::addConsumed<long long>(TString name);
template void BaseTreeLooper::addConsumed<float>(TString name);
template void BaseTreeLooper::addConsumed<double>(TString name);
template void BaseTreeLooper::addConsumed<std::vector<short>>(TString name);
template void BaseTreeLooper::addConsumed<std::vector<unsigned int>>(TString name);
template void BaseTreeLooper::addConsumed<std::vector<int>>(TString name);
template void BaseTreeLooper::addConsumed<std::vector<unsigned long>>(TString name);
template void BaseTreeLooper::addConsumed<std::vector<long>>(TString name);
template void BaseTreeLooper::addConsumed<std::vector<long long>>(TString name);
template void BaseTreeLooper::addConsumed<std::vector<float>>(TString name);
template void BaseTreeLooper::addConsumed<std::vector<double>>(TString name);

template<typename T> bool BaseTreeLooper::linkConsumed(CJLSTTree* tree){
  bool result=true;
  std::unordered_map<TString, T*>* theMap = nullptr;
  BaseTreeLooper::getConsumedMap<T>(theMap);
  if (theMap){
    for (typename std::unordered_map<TString, T*>::iterator it=theMap->begin(); it!=theMap->end(); it++){
      if (tree->branchExists(it->first)){ tree->getValRef<T>(it->first, it->second); result &= true; }
      else{ result &= false; MELAerr << "BaseTreeLooper::linkConsumed(" << tree->sampleIdentifier << "): Linking failed for variable " << it->first << endl; }
    }
  }
  else{
    MELAerr << "BaseTreeLooper::linkConsumed(" << tree->sampleIdentifier << "): Map could not be found." << endl;
    result = false;
  }
  return result;
}
template bool BaseTreeLooper::linkConsumed<short>(CJLSTTree* tree);
template bool BaseTreeLooper::linkConsumed<unsigned int>(CJLSTTree* tree);
template bool BaseTreeLooper::linkConsumed<int>(CJLSTTree* tree);
template bool BaseTreeLooper::linkConsumed<unsigned long>(CJLSTTree* tree);
template bool BaseTreeLooper::linkConsumed<long>(CJLSTTree* tree);
template bool BaseTreeLooper::linkConsumed<long long>(CJLSTTree* tree);
template bool BaseTreeLooper::linkConsumed<float>(CJLSTTree* tree);
template bool BaseTreeLooper::linkConsumed<double>(CJLSTTree* tree);
template bool BaseTreeLooper::linkConsumed<std::vector<short>>(CJLSTTree* tree);
template bool BaseTreeLooper::linkConsumed<std::vector<unsigned int>>(CJLSTTree* tree);
template bool BaseTreeLooper::linkConsumed<std::vector<int>>(CJLSTTree* tree);
template bool BaseTreeLooper::linkConsumed<std::vector<unsigned long>>(CJLSTTree* tree);
template bool BaseTreeLooper::linkConsumed<std::vector<long>>(CJLSTTree* tree);
template bool BaseTreeLooper::linkConsumed<std::vector<long long>>(CJLSTTree* tree);
template bool BaseTreeLooper::linkConsumed<std::vector<float>>(CJLSTTree* tree);
template bool BaseTreeLooper::linkConsumed<std::vector<double>>(CJLSTTree* tree);


template<typename T> bool BaseTreeLooper::getConsumed(TString name, T const*& val) const{
  std::unordered_map<TString, T*> const* theMap = nullptr;
  BaseTreeLooper::getConsumedMap<T>(theMap);
  if (theMap){
    typename unordered_map<TString, T*>::const_iterator it_val;
    if (HelperFunctions::getUnorderedMapIterator(name, *theMap, it_val)){
      val = it_val->second;
      return true;
    }
  }
  else MELAerr << "BaseTreeLooper::getConsumed(" << name << "): Map could not be found." << endl;
  return false;
}
template bool BaseTreeLooper::getConsumed<short>(TString name, short const*& val) const;
template bool BaseTreeLooper::getConsumed<unsigned int>(TString name, unsigned int const*& val) const;
template bool BaseTreeLooper::getConsumed<int>(TString name, int const*& val) const;
template bool BaseTreeLooper::getConsumed<unsigned long>(TString name, unsigned long const*& val) const;
template bool BaseTreeLooper::getConsumed<long>(TString name, long const*& val) const;
template bool BaseTreeLooper::getConsumed<long long>(TString name, long long const*& val) const;
template bool BaseTreeLooper::getConsumed<float>(TString name, float const*& val) const;
template bool BaseTreeLooper::getConsumed<double>(TString name, double const*& val) const;
template bool BaseTreeLooper::getConsumed<std::vector<short>>(TString name, std::vector<short> const*& val) const;
template bool BaseTreeLooper::getConsumed<std::vector<unsigned int>>(TString name, std::vector<unsigned int> const*& val) const;
template bool BaseTreeLooper::getConsumed<std::vector<int>>(TString name, std::vector<int> const*& val) const;
template bool BaseTreeLooper::getConsumed<std::vector<unsigned long>>(TString name, std::vector<unsigned long> const*& val) const;
template bool BaseTreeLooper::getConsumed<std::vector<long>>(TString name, std::vector<long> const*& val) const;
template bool BaseTreeLooper::getConsumed<std::vector<long long>>(TString name, std::vector<long long> const*& val) const;
template bool BaseTreeLooper::getConsumed<std::vector<float>>(TString name, std::vector<float> const*& val) const;
template bool BaseTreeLooper::getConsumed<std::vector<double>>(TString name, std::vector<double> const*& val) const;


#endif
