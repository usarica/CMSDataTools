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

  std::unordered_map<TString, TBits*> valTBitss;
  std::unordered_map<TString, bool*> valbools;
  std::unordered_map<TString, short*> valshorts;
  std::unordered_map<TString, unsigned int*> valuints;
  std::unordered_map<TString, int*> valints;
  std::unordered_map<TString, unsigned long*> valulongs;
  std::unordered_map<TString, long*> vallongs;
  std::unordered_map<TString, unsigned long long*> valulonglongs;
  std::unordered_map<TString, long long*> vallonglongs;
  std::unordered_map<TString, float*> valfloats;
  std::unordered_map<TString, double*> valdoubles;
  std::unordered_map<TString, std::string*> valstrings;
  std::unordered_map<TString, TString*> valTStrings;
  std::unordered_map<TString, CMSLorentzVector*> valCMSLorentzVectors;

  std::unordered_map<TString, std::vector<bool>**> valVbools;
  std::unordered_map<TString, std::vector<short>**> valVshorts;
  std::unordered_map<TString, std::vector<unsigned int>**> valVuints;
  std::unordered_map<TString, std::vector<int>**> valVints;
  std::unordered_map<TString, std::vector<unsigned long>**> valVulongs;
  std::unordered_map<TString, std::vector<long>**> valVlongs;
  std::unordered_map<TString, std::vector<unsigned long long>**> valVulonglongs;
  std::unordered_map<TString, std::vector<long long>**> valVlonglongs;
  std::unordered_map<TString, std::vector<float>**> valVfloats;
  std::unordered_map<TString, std::vector<double>**> valVdoubles;
  std::unordered_map<TString, std::vector<std::string>**> valVstrings;
  std::unordered_map<TString, std::vector<TString>**> valVTStrings;
  std::unordered_map<TString, std::vector<CMSLorentzVector>**> valVCMSLorentzVectors;

  std::unordered_map<TString, std::vector<std::vector<bool>>**> valVVbools;
  std::unordered_map<TString, std::vector<std::vector<short>>**> valVVshorts;
  std::unordered_map<TString, std::vector<std::vector<unsigned int>>**> valVVuints;
  std::unordered_map<TString, std::vector<std::vector<int>>**> valVVints;
  std::unordered_map<TString, std::vector<std::vector<unsigned long>>**> valVVulongs;
  std::unordered_map<TString, std::vector<std::vector<long>>**> valVVlongs;
  std::unordered_map<TString, std::vector<std::vector<unsigned long long>>**> valVVulonglongs;
  std::unordered_map<TString, std::vector<std::vector<long long>>**> valVVlonglongs;
  std::unordered_map<TString, std::vector<std::vector<float>>**> valVVfloats;
  std::unordered_map<TString, std::vector<std::vector<double>>**> valVVdoubles;


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