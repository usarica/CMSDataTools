#ifndef IVYBASE_H
#define IVYBASE_H

#include "TVar.hh"
#include "BaseTree.h"
#include "HelperFunctionsCore.h"


class IvyBase{
protected:
  TVar::VerbosityLevel verbosity;
  BaseTree* currentTree;

  // Consumes
  std::vector<TString> sloppyConsumes; // In case some variables are known to be absent in some trees

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
  std::unordered_map<TString, CMSLorentzVector*> valCMSLorentzVectors;

  std::unordered_map<TString, std::vector<bool>* const*> valVbools;
  std::unordered_map<TString, std::vector<short>* const*> valVshorts;
  std::unordered_map<TString, std::vector<unsigned int>* const*> valVuints;
  std::unordered_map<TString, std::vector<int>* const*> valVints;
  std::unordered_map<TString, std::vector<unsigned long>* const*> valVulongs;
  std::unordered_map<TString, std::vector<long>* const*> valVlongs;
  std::unordered_map<TString, std::vector<unsigned long long>* const*> valVulonglongs;
  std::unordered_map<TString, std::vector<long long>* const*> valVlonglongs;
  std::unordered_map<TString, std::vector<float>* const*> valVfloats;
  std::unordered_map<TString, std::vector<double>* const*> valVdoubles;
  std::unordered_map<TString, std::vector<std::string>* const*> valVstrings;
  std::unordered_map<TString, std::vector<CMSLorentzVector>* const*> valVCMSLorentzVectors;

  template<typename T> bool linkConsumed(BaseTree* tree);
  bool linkConsumes(BaseTree* tree);

  // Get consumed map
  template<typename T> void getConsumedMap(std::unordered_map<TString, T*>*& theMap);
  template<typename T> void getConsumedMap(std::unordered_map<TString, T*> const*& theMap) const;

  // Get consumed
  template<typename T> bool getConsumed(TString name, T const*& val) const;
  template<typename T> bool getConsumed(TString name, T*& val) const;

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


#endif
