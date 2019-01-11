#ifndef IVYBASE_HPP
#define IVYBASE_HPP

#include "IvyBase.h"


using namespace std;
using namespace MELAStreamHelpers;


template<> void IvyBase::getConsumedMap<TBits>(std::unordered_map<TString, TBits*>*& theMap){ theMap = &valTBitss; }
template<> void IvyBase::getConsumedMap<bool>(std::unordered_map<TString, bool*>*& theMap){ theMap = &valbools; }
template<> void IvyBase::getConsumedMap<short>(std::unordered_map<TString, short*>*& theMap){ theMap = &valshorts; }
template<> void IvyBase::getConsumedMap<unsigned int>(std::unordered_map<TString, unsigned int*>*& theMap){ theMap = &valuints; }
template<> void IvyBase::getConsumedMap<int>(std::unordered_map<TString, int*>*& theMap){ theMap = &valints; }
template<> void IvyBase::getConsumedMap<unsigned long>(std::unordered_map<TString, unsigned long*>*& theMap){ theMap = &valulongs; }
template<> void IvyBase::getConsumedMap<long>(std::unordered_map<TString, long*>*& theMap){ theMap = &vallongs; }
template<> void IvyBase::getConsumedMap<unsigned long long>(std::unordered_map<TString, unsigned long long*>*& theMap){ theMap = &valulonglongs; }
template<> void IvyBase::getConsumedMap<long long>(std::unordered_map<TString, long long*>*& theMap){ theMap = &vallonglongs; }
template<> void IvyBase::getConsumedMap<float>(std::unordered_map<TString, float*>*& theMap){ theMap = &valfloats; }
template<> void IvyBase::getConsumedMap<double>(std::unordered_map<TString, double*>*& theMap){ theMap = &valdoubles; }
template<> void IvyBase::getConsumedMap<std::string>(std::unordered_map<TString, std::string*>*& theMap){ theMap = &valstrings; }
template<> void IvyBase::getConsumedMap<TString>(std::unordered_map<TString, TString*>*& theMap){ theMap = &valTStrings; }
template<> void IvyBase::getConsumedMap<CMSLorentzVector>(std::unordered_map<TString, CMSLorentzVector*>*& theMap){ theMap = &valCMSLorentzVectors; }
template<> void IvyBase::getConsumedMap<std::vector<bool>*>(std::unordered_map<TString, std::vector<bool>**>*& theMap){ theMap = &valVbools; }
template<> void IvyBase::getConsumedMap<std::vector<short>*>(std::unordered_map<TString, std::vector<short>**>*& theMap){ theMap = &valVshorts; }
template<> void IvyBase::getConsumedMap<std::vector<unsigned int>*>(std::unordered_map<TString, std::vector<unsigned int>**>*& theMap){ theMap = &valVuints; }
template<> void IvyBase::getConsumedMap<std::vector<int>*>(std::unordered_map<TString, std::vector<int>**>*& theMap){ theMap = &valVints; }
template<> void IvyBase::getConsumedMap<std::vector<unsigned long>*>(std::unordered_map<TString, std::vector<unsigned long>**>*& theMap){ theMap = &valVulongs; }
template<> void IvyBase::getConsumedMap<std::vector<long>*>(std::unordered_map<TString, std::vector<long>**>*& theMap){ theMap = &valVlongs; }
template<> void IvyBase::getConsumedMap<std::vector<unsigned long long>*>(std::unordered_map<TString, std::vector<unsigned long long>**>*& theMap){ theMap = &valVulonglongs; }
template<> void IvyBase::getConsumedMap<std::vector<long long>*>(std::unordered_map<TString, std::vector<long long>**>*& theMap){ theMap = &valVlonglongs; }
template<> void IvyBase::getConsumedMap<std::vector<float>*>(std::unordered_map<TString, std::vector<float>**>*& theMap){ theMap = &valVfloats; }
template<> void IvyBase::getConsumedMap<std::vector<double>*>(std::unordered_map<TString, std::vector<double>**>*& theMap){ theMap = &valVdoubles; }
template<> void IvyBase::getConsumedMap<std::vector<std::string>*>(std::unordered_map<TString, std::vector<std::string>**>*& theMap){ theMap = &valVstrings; }
template<> void IvyBase::getConsumedMap<std::vector<TString>*>(std::unordered_map<TString, std::vector<TString>**>*& theMap){ theMap = &valVTStrings; }
template<> void IvyBase::getConsumedMap<std::vector<CMSLorentzVector>*>(std::unordered_map<TString, std::vector<CMSLorentzVector>**>*& theMap){ theMap = &valVCMSLorentzVectors; }
template<> void IvyBase::getConsumedMap<std::vector<std::vector<bool>>*>(std::unordered_map<TString, std::vector<std::vector<bool>>**>*& theMap){ theMap = &valVVbools; }
template<> void IvyBase::getConsumedMap<std::vector<std::vector<short>>*>(std::unordered_map<TString, std::vector<std::vector<short>>**>*& theMap){ theMap = &valVVshorts; }
template<> void IvyBase::getConsumedMap<std::vector<std::vector<unsigned int>>*>(std::unordered_map<TString, std::vector<std::vector<unsigned int>>**>*& theMap){ theMap = &valVVuints; }
template<> void IvyBase::getConsumedMap<std::vector<std::vector<int>>*>(std::unordered_map<TString, std::vector<std::vector<int>>**>*& theMap){ theMap = &valVVints; }
template<> void IvyBase::getConsumedMap<std::vector<std::vector<unsigned long>>*>(std::unordered_map<TString, std::vector<std::vector<unsigned long>>**>*& theMap){ theMap = &valVVulongs; }
template<> void IvyBase::getConsumedMap<std::vector<std::vector<long>>*>(std::unordered_map<TString, std::vector<std::vector<long>>**>*& theMap){ theMap = &valVVlongs; }
template<> void IvyBase::getConsumedMap<std::vector<std::vector<unsigned long long>>*>(std::unordered_map<TString, std::vector<std::vector<unsigned long long>>**>*& theMap){ theMap = &valVVulonglongs; }
template<> void IvyBase::getConsumedMap<std::vector<std::vector<long long>>*>(std::unordered_map<TString, std::vector<std::vector<long long>>**>*& theMap){ theMap = &valVVlonglongs; }
template<> void IvyBase::getConsumedMap<std::vector<std::vector<float>>*>(std::unordered_map<TString, std::vector<std::vector<float>>**>*& theMap){ theMap = &valVVfloats; }
template<> void IvyBase::getConsumedMap<std::vector<std::vector<double>>*>(std::unordered_map<TString, std::vector<std::vector<double>>**>*& theMap){ theMap = &valVVdoubles; }

template<> void IvyBase::getConsumedMap<TBits>(std::unordered_map<TString, TBits*> const*& theMap) const{ theMap = &valTBitss; }
template<> void IvyBase::getConsumedMap<bool>(std::unordered_map<TString, bool*> const*& theMap) const{ theMap = &valbools; }
template<> void IvyBase::getConsumedMap<short>(std::unordered_map<TString, short*> const*& theMap) const{ theMap = &valshorts; }
template<> void IvyBase::getConsumedMap<unsigned int>(std::unordered_map<TString, unsigned int*> const*& theMap) const{ theMap = &valuints; }
template<> void IvyBase::getConsumedMap<int>(std::unordered_map<TString, int*> const*& theMap) const{ theMap = &valints; }
template<> void IvyBase::getConsumedMap<unsigned long>(std::unordered_map<TString, unsigned long*> const*& theMap) const{ theMap = &valulongs; }
template<> void IvyBase::getConsumedMap<long>(std::unordered_map<TString, long*> const*& theMap) const{ theMap = &vallongs; }
template<> void IvyBase::getConsumedMap<unsigned long long>(std::unordered_map<TString, unsigned long long*> const*& theMap) const{ theMap = &valulonglongs; }
template<> void IvyBase::getConsumedMap<long long>(std::unordered_map<TString, long long*> const*& theMap) const{ theMap = &vallonglongs; }
template<> void IvyBase::getConsumedMap<float>(std::unordered_map<TString, float*> const*& theMap) const{ theMap = &valfloats; }
template<> void IvyBase::getConsumedMap<double>(std::unordered_map<TString, double*> const*& theMap) const{ theMap = &valdoubles; }
template<> void IvyBase::getConsumedMap<std::string>(std::unordered_map<TString, std::string*> const*& theMap) const{ theMap = &valstrings; }
template<> void IvyBase::getConsumedMap<TString>(std::unordered_map<TString, TString*> const*& theMap) const{ theMap = &valTStrings; }
template<> void IvyBase::getConsumedMap<CMSLorentzVector>(std::unordered_map<TString, CMSLorentzVector*> const*& theMap) const{ theMap = &valCMSLorentzVectors; }
template<> void IvyBase::getConsumedMap<std::vector<bool>*>(std::unordered_map<TString, std::vector<bool>**> const*& theMap) const{ theMap = &valVbools; }
template<> void IvyBase::getConsumedMap<std::vector<short>*>(std::unordered_map<TString, std::vector<short>**> const*& theMap) const{ theMap = &valVshorts; }
template<> void IvyBase::getConsumedMap<std::vector<unsigned int>*>(std::unordered_map<TString, std::vector<unsigned int>**> const*& theMap) const{ theMap = &valVuints; }
template<> void IvyBase::getConsumedMap<std::vector<int>*>(std::unordered_map<TString, std::vector<int>**> const*& theMap) const{ theMap = &valVints; }
template<> void IvyBase::getConsumedMap<std::vector<unsigned long>*>(std::unordered_map<TString, std::vector<unsigned long>**> const*& theMap) const{ theMap = &valVulongs; }
template<> void IvyBase::getConsumedMap<std::vector<long>*>(std::unordered_map<TString, std::vector<long>**> const*& theMap) const{ theMap = &valVlongs; }
template<> void IvyBase::getConsumedMap<std::vector<unsigned long long>*>(std::unordered_map<TString, std::vector<unsigned long long>**> const*& theMap) const{ theMap = &valVulonglongs; }
template<> void IvyBase::getConsumedMap<std::vector<long long>*>(std::unordered_map<TString, std::vector<long long>**> const*& theMap) const{ theMap = &valVlonglongs; }
template<> void IvyBase::getConsumedMap<std::vector<float>*>(std::unordered_map<TString, std::vector<float>**> const*& theMap) const{ theMap = &valVfloats; }
template<> void IvyBase::getConsumedMap<std::vector<double>*>(std::unordered_map<TString, std::vector<double>**> const*& theMap) const{ theMap = &valVdoubles; }
template<> void IvyBase::getConsumedMap<std::vector<std::string>*>(std::unordered_map<TString, std::vector<std::string>**> const*& theMap) const{ theMap = &valVstrings; }
template<> void IvyBase::getConsumedMap<std::vector<TString>*>(std::unordered_map<TString, std::vector<TString>**> const*& theMap) const{ theMap = &valVTStrings; }
template<> void IvyBase::getConsumedMap<std::vector<CMSLorentzVector>*>(std::unordered_map<TString, std::vector<CMSLorentzVector>**> const*& theMap) const{ theMap = &valVCMSLorentzVectors; }
template<> void IvyBase::getConsumedMap<std::vector<std::vector<bool>>*>(std::unordered_map<TString, std::vector<std::vector<bool>>**> const*& theMap) const{ theMap = &valVVbools; }
template<> void IvyBase::getConsumedMap<std::vector<std::vector<short>>*>(std::unordered_map<TString, std::vector<std::vector<short>>**> const*& theMap) const{ theMap = &valVVshorts; }
template<> void IvyBase::getConsumedMap<std::vector<std::vector<unsigned int>>*>(std::unordered_map<TString, std::vector<std::vector<unsigned int>>**> const*& theMap) const{ theMap = &valVVuints; }
template<> void IvyBase::getConsumedMap<std::vector<std::vector<int>>*>(std::unordered_map<TString, std::vector<std::vector<int>>**> const*& theMap) const{ theMap = &valVVints; }
template<> void IvyBase::getConsumedMap<std::vector<std::vector<unsigned long>>*>(std::unordered_map<TString, std::vector<std::vector<unsigned long>>**> const*& theMap) const{ theMap = &valVVulongs; }
template<> void IvyBase::getConsumedMap<std::vector<std::vector<long>>*>(std::unordered_map<TString, std::vector<std::vector<long>>**> const*& theMap) const{ theMap = &valVVlongs; }
template<> void IvyBase::getConsumedMap<std::vector<std::vector<unsigned long long>>*>(std::unordered_map<TString, std::vector<std::vector<unsigned long long>>**> const*& theMap) const{ theMap = &valVVulonglongs; }
template<> void IvyBase::getConsumedMap<std::vector<std::vector<long long>>*>(std::unordered_map<TString, std::vector<std::vector<long long>>**> const*& theMap) const{ theMap = &valVVlonglongs; }
template<> void IvyBase::getConsumedMap<std::vector<std::vector<float>>*>(std::unordered_map<TString, std::vector<std::vector<float>>**> const*& theMap) const{ theMap = &valVVfloats; }
template<> void IvyBase::getConsumedMap<std::vector<std::vector<double>>*>(std::unordered_map<TString, std::vector<std::vector<double>>**> const*& theMap) const{ theMap = &valVVdoubles; }


template<typename T> void IvyBase::addConsumed(TString name){
  std::unordered_map<TString, T*>* theMap = nullptr;
  IvyBase::getConsumedMap<T>(theMap);
  if (theMap) (*theMap)[name] = nullptr;
  else if (verbosity>=TVar::ERROR) MELAerr << "IvyBase::addConsumed(" << name << "): Map could not be found." << endl;
}
template void IvyBase::addConsumed<TBits>(TString name);
template void IvyBase::addConsumed<bool>(TString name);
template void IvyBase::addConsumed<short>(TString name);
template void IvyBase::addConsumed<unsigned int>(TString name);
template void IvyBase::addConsumed<int>(TString name);
template void IvyBase::addConsumed<unsigned long>(TString name);
template void IvyBase::addConsumed<long>(TString name);
template void IvyBase::addConsumed<unsigned long long>(TString name);
template void IvyBase::addConsumed<long long>(TString name);
template void IvyBase::addConsumed<float>(TString name);
template void IvyBase::addConsumed<double>(TString name);
template void IvyBase::addConsumed<std::string>(TString name);
template void IvyBase::addConsumed<TString>(TString name);
template void IvyBase::addConsumed<CMSLorentzVector>(TString name);
template void IvyBase::addConsumed<std::vector<bool>*>(TString name);
template void IvyBase::addConsumed<std::vector<short>*>(TString name);
template void IvyBase::addConsumed<std::vector<unsigned int>*>(TString name);
template void IvyBase::addConsumed<std::vector<int>*>(TString name);
template void IvyBase::addConsumed<std::vector<unsigned long>*>(TString name);
template void IvyBase::addConsumed<std::vector<long>*>(TString name);
template void IvyBase::addConsumed<std::vector<unsigned long long>*>(TString name);
template void IvyBase::addConsumed<std::vector<long long>*>(TString name);
template void IvyBase::addConsumed<std::vector<float>*>(TString name);
template void IvyBase::addConsumed<std::vector<double>*>(TString name);
template void IvyBase::addConsumed<std::vector<std::string>*>(TString name);
template void IvyBase::addConsumed<std::vector<TString>*>(TString name);
template void IvyBase::addConsumed<std::vector<CMSLorentzVector>*>(TString name);
template void IvyBase::addConsumed<std::vector<std::vector<bool>>*>(TString name);
template void IvyBase::addConsumed<std::vector<std::vector<short>>*>(TString name);
template void IvyBase::addConsumed<std::vector<std::vector<unsigned int>>*>(TString name);
template void IvyBase::addConsumed<std::vector<std::vector<int>>*>(TString name);
template void IvyBase::addConsumed<std::vector<std::vector<unsigned long>>*>(TString name);
template void IvyBase::addConsumed<std::vector<std::vector<long>>*>(TString name);
template void IvyBase::addConsumed<std::vector<std::vector<unsigned long long>>*>(TString name);
template void IvyBase::addConsumed<std::vector<std::vector<long long>>*>(TString name);
template void IvyBase::addConsumed<std::vector<std::vector<float>>*>(TString name);
template void IvyBase::addConsumed<std::vector<std::vector<double>>*>(TString name);

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
template bool IvyBase::linkConsumed<TBits>(BaseTree* tree);
template bool IvyBase::linkConsumed<bool>(BaseTree* tree);
template bool IvyBase::linkConsumed<short>(BaseTree* tree);
template bool IvyBase::linkConsumed<unsigned int>(BaseTree* tree);
template bool IvyBase::linkConsumed<int>(BaseTree* tree);
template bool IvyBase::linkConsumed<unsigned long>(BaseTree* tree);
template bool IvyBase::linkConsumed<long>(BaseTree* tree);
template bool IvyBase::linkConsumed<unsigned long long>(BaseTree* tree);
template bool IvyBase::linkConsumed<long long>(BaseTree* tree);
template bool IvyBase::linkConsumed<float>(BaseTree* tree);
template bool IvyBase::linkConsumed<double>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::string>(BaseTree* tree);
template bool IvyBase::linkConsumed<TString>(BaseTree* tree);
template bool IvyBase::linkConsumed<CMSLorentzVector>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<bool>*>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<short>*>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<unsigned int>*>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<int>*>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<unsigned long>*>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<long>*>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<unsigned long long>*>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<long long>*>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<float>*>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<double>*>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<std::string>*>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<TString>*>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<CMSLorentzVector>*>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<std::vector<bool>>*>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<std::vector<short>>*>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<std::vector<unsigned int>>*>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<std::vector<int>>*>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<std::vector<unsigned long>>*>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<std::vector<long>>*>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<std::vector<unsigned long long>>*>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<std::vector<long long>>*>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<std::vector<float>>*>(BaseTree* tree);
template bool IvyBase::linkConsumed<std::vector<std::vector<double>>*>(BaseTree* tree);


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
  }
  else if (isSloppy) return true;
  else if (verbosity>=TVar::ERROR) MELAerr << "IvyBase::getConsumed(" << name << "): Map could not be found." << endl;
  return false;
}
template bool IvyBase::getConsumed<TBits>(TString name, TBits*& val) const;
template bool IvyBase::getConsumed<bool>(TString name, bool*& val) const;
template bool IvyBase::getConsumed<short>(TString name, short*& val) const;
template bool IvyBase::getConsumed<unsigned int>(TString name, unsigned int*& val) const;
template bool IvyBase::getConsumed<int>(TString name, int*& val) const;
template bool IvyBase::getConsumed<unsigned long>(TString name, unsigned long*& val) const;
template bool IvyBase::getConsumed<long>(TString name, long*& val) const;
template bool IvyBase::getConsumed<unsigned long long>(TString name, unsigned long long*& val) const;
template bool IvyBase::getConsumed<long long>(TString name, long long*& val) const;
template bool IvyBase::getConsumed<float>(TString name, float*& val) const;
template bool IvyBase::getConsumed<double>(TString name, double*& val) const;
template bool IvyBase::getConsumed<std::string>(TString name, std::string*& val) const;
template bool IvyBase::getConsumed<TString>(TString name, TString*& val) const;
template bool IvyBase::getConsumed<CMSLorentzVector>(TString name, CMSLorentzVector*& val) const;
template bool IvyBase::getConsumed<std::vector<bool>*>(TString name, std::vector<bool>**& val) const;
template bool IvyBase::getConsumed<std::vector<short>*>(TString name, std::vector<short>**& val) const;
template bool IvyBase::getConsumed<std::vector<unsigned int>*>(TString name, std::vector<unsigned int>**& val) const;
template bool IvyBase::getConsumed<std::vector<int>*>(TString name, std::vector<int>**& val) const;
template bool IvyBase::getConsumed<std::vector<unsigned long>*>(TString name, std::vector<unsigned long>**& val) const;
template bool IvyBase::getConsumed<std::vector<long>*>(TString name, std::vector<long>**& val) const;
template bool IvyBase::getConsumed<std::vector<unsigned long long>*>(TString name, std::vector<unsigned long long>**& val) const;
template bool IvyBase::getConsumed<std::vector<long long>*>(TString name, std::vector<long long>**& val) const;
template bool IvyBase::getConsumed<std::vector<float>*>(TString name, std::vector<float>**& val) const;
template bool IvyBase::getConsumed<std::vector<double>*>(TString name, std::vector<double>**& val) const;
template bool IvyBase::getConsumed<std::vector<std::string>*>(TString name, std::vector<std::string>**& val) const;
template bool IvyBase::getConsumed<std::vector<TString>*>(TString name, std::vector<TString>**& val) const;
template bool IvyBase::getConsumed<std::vector<CMSLorentzVector>*>(TString name, std::vector<CMSLorentzVector>**& val) const;
template bool IvyBase::getConsumed<std::vector<std::vector<bool>>*>(TString name, std::vector<std::vector<bool>>**& val) const;
template bool IvyBase::getConsumed<std::vector<std::vector<short>>*>(TString name, std::vector<std::vector<short>>**& val) const;
template bool IvyBase::getConsumed<std::vector<std::vector<unsigned int>>*>(TString name, std::vector<std::vector<unsigned int>>**& val) const;
template bool IvyBase::getConsumed<std::vector<std::vector<int>>*>(TString name, std::vector<std::vector<int>>**& val) const;
template bool IvyBase::getConsumed<std::vector<std::vector<unsigned long>>*>(TString name, std::vector<std::vector<unsigned long>>**& val) const;
template bool IvyBase::getConsumed<std::vector<std::vector<long>>*>(TString name, std::vector<std::vector<long>>**& val) const;
template bool IvyBase::getConsumed<std::vector<std::vector<unsigned long long>>*>(TString name, std::vector<std::vector<unsigned long long>>**& val) const;
template bool IvyBase::getConsumed<std::vector<std::vector<long long>>*>(TString name, std::vector<std::vector<long long>>**& val) const;
template bool IvyBase::getConsumed<std::vector<std::vector<float>>*>(TString name, std::vector<std::vector<float>>**& val) const;
template bool IvyBase::getConsumed<std::vector<std::vector<double>>*>(TString name, std::vector<std::vector<double>>**& val) const;

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
  }
  else if (isSloppy) return true;
  else if (verbosity>=TVar::ERROR) MELAerr << "IvyBase::getConsumed(" << name << "): Map could not be found." << endl;
  return false;
}
template bool IvyBase::getConsumed<TBits>(TString name, TBits const*& val) const;
template bool IvyBase::getConsumed<bool>(TString name, bool const*& val) const;
template bool IvyBase::getConsumed<short>(TString name, short const*& val) const;
template bool IvyBase::getConsumed<unsigned int>(TString name, unsigned int const*& val) const;
template bool IvyBase::getConsumed<int>(TString name, int const*& val) const;
template bool IvyBase::getConsumed<unsigned long>(TString name, unsigned long const*& val) const;
template bool IvyBase::getConsumed<long>(TString name, long const*& val) const;
template bool IvyBase::getConsumed<unsigned long long>(TString name, unsigned long long const*& val) const;
template bool IvyBase::getConsumed<long long>(TString name, long long const*& val) const;
template bool IvyBase::getConsumed<float>(TString name, float const*& val) const;
template bool IvyBase::getConsumed<double>(TString name, double const*& val) const;
template bool IvyBase::getConsumed<std::string>(TString name, std::string const*& val) const;
template bool IvyBase::getConsumed<TString>(TString name, TString const*& val) const;
template bool IvyBase::getConsumed<CMSLorentzVector>(TString name, CMSLorentzVector const*& val) const;
template bool IvyBase::getConsumed<std::vector<bool>*>(TString name, std::vector<bool>* const*& val) const;
template bool IvyBase::getConsumed<std::vector<short>*>(TString name, std::vector<short>* const*& val) const;
template bool IvyBase::getConsumed<std::vector<unsigned int>*>(TString name, std::vector<unsigned int>* const*& val) const;
template bool IvyBase::getConsumed<std::vector<int>*>(TString name, std::vector<int>* const*& val) const;
template bool IvyBase::getConsumed<std::vector<unsigned long>*>(TString name, std::vector<unsigned long>* const*& val) const;
template bool IvyBase::getConsumed<std::vector<long>*>(TString name, std::vector<long>* const*& val) const;
template bool IvyBase::getConsumed<std::vector<unsigned long long>*>(TString name, std::vector<unsigned long long>* const*& val) const;
template bool IvyBase::getConsumed<std::vector<long long>*>(TString name, std::vector<long long>* const*& val) const;
template bool IvyBase::getConsumed<std::vector<float>*>(TString name, std::vector<float>* const*& val) const;
template bool IvyBase::getConsumed<std::vector<double>*>(TString name, std::vector<double>* const*& val) const;
template bool IvyBase::getConsumed<std::vector<std::string>*>(TString name, std::vector<std::string>* const*& val) const;
template bool IvyBase::getConsumed<std::vector<TString>*>(TString name, std::vector<TString>* const*& val) const;
template bool IvyBase::getConsumed<std::vector<CMSLorentzVector>*>(TString name, std::vector<CMSLorentzVector>* const*& val) const;
template bool IvyBase::getConsumed<std::vector<std::vector<bool>>*>(TString name, std::vector<std::vector<bool>>* const*& val) const;
template bool IvyBase::getConsumed<std::vector<std::vector<short>>*>(TString name, std::vector<std::vector<short>>* const*& val) const;
template bool IvyBase::getConsumed<std::vector<std::vector<unsigned int>>*>(TString name, std::vector<std::vector<unsigned int>>* const*& val) const;
template bool IvyBase::getConsumed<std::vector<std::vector<int>>*>(TString name, std::vector<std::vector<int>>* const*& val) const;
template bool IvyBase::getConsumed<std::vector<std::vector<unsigned long>>*>(TString name, std::vector<std::vector<unsigned long>>* const*& val) const;
template bool IvyBase::getConsumed<std::vector<std::vector<long>>*>(TString name, std::vector<std::vector<long>>* const*& val) const;
template bool IvyBase::getConsumed<std::vector<std::vector<unsigned long long>>*>(TString name, std::vector<std::vector<unsigned long long>>* const*& val) const;
template bool IvyBase::getConsumed<std::vector<std::vector<long long>>*>(TString name, std::vector<std::vector<long long>>* const*& val) const;
template bool IvyBase::getConsumed<std::vector<std::vector<float>>*>(TString name, std::vector<std::vector<float>>* const*& val) const;
template bool IvyBase::getConsumed<std::vector<std::vector<double>>*>(TString name, std::vector<std::vector<double>>* const*& val) const;

#endif
