#ifndef BASETREE_H
#define BASETREE_H

#include <vector>
#include <fstream>
#include <cstdlib>
#include <unordered_map>
#include "StdExtensions.h"
#include "SimpleEntry.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "CMSLorentzVector.h"


class BaseTree{
public:
  enum BranchType{
    BranchType_bool_t,
    BranchType_short_t,
    BranchType_uint_t,
    BranchType_int_t,
    BranchType_ulong_t,
    BranchType_long_t,
    BranchType_ulonglong_t,
    BranchType_longlong_t,
    BranchType_float_t,
    BranchType_double_t,
    BranchType_string_t,
    BranchType_TString_t,
    BranchType_CMSLorentzVector_t,

    BranchType_vbool_t,
    BranchType_vshort_t,
    BranchType_vuint_t,
    BranchType_vint_t,
    BranchType_vulong_t,
    BranchType_vlong_t,
    BranchType_vulonglong_t,
    BranchType_vlonglong_t,
    BranchType_vfloat_t,
    BranchType_vdouble_t,
    BranchType_vstring_t,
    BranchType_vTString_t,
    BranchType_vCMSLorentzVector_t,

    BranchType_vvbool_t,
    BranchType_vvshort_t,
    BranchType_vvuint_t,
    BranchType_vvint_t,
    BranchType_vvulong_t,
    BranchType_vvlong_t,
    BranchType_vvulonglong_t,
    BranchType_vvlonglong_t,
    BranchType_vvfloat_t,
    BranchType_vvdouble_t,

    BranchType_unknown_t
  };

  TString sampleIdentifier;

protected:
  TFile* finput;
  TTree* tree;
  TTree* failedtree;
  TH1F* hCounters;
  bool valid;
  const bool receiver;

  int currentEvent;
  TTree* currentTree;

  std::unordered_map<TString, std::pair<bool, bool>*> valbools;
  std::unordered_map<TString, std::pair<short, short>*> valshorts;
  std::unordered_map<TString, std::pair<unsigned int, unsigned int>*> valuints;
  std::unordered_map<TString, std::pair<int, int>*> valints;
  std::unordered_map<TString, std::pair<unsigned long, unsigned long>*> valulongs;
  std::unordered_map<TString, std::pair<long, long>*> vallongs;
  std::unordered_map<TString, std::pair<unsigned long long, unsigned long long>*> valulonglongs;
  std::unordered_map<TString, std::pair<long long, long long>*> vallonglongs;
  std::unordered_map<TString, std::pair<float, float>*> valfloats;
  std::unordered_map<TString, std::pair<double, double>*> valdoubles;
  std::unordered_map<TString, std::pair<std::string, std::string>*> valstrings;
  std::unordered_map<TString, std::pair<TString, TString>*> valTStrings;
  std::unordered_map<TString, std::pair<CMSLorentzVector, CMSLorentzVector>*> valCMSLorentzVectors;

  std::unordered_map<TString, std::vector<bool>*> valVbools;
  std::unordered_map<TString, std::vector<short>*> valVshorts;
  std::unordered_map<TString, std::vector<unsigned int>*> valVuints;
  std::unordered_map<TString, std::vector<int>*> valVints;
  std::unordered_map<TString, std::vector<unsigned long>*> valVulongs;
  std::unordered_map<TString, std::vector<long>*> valVlongs;
  std::unordered_map<TString, std::vector<unsigned long long>*> valVulonglongs;
  std::unordered_map<TString, std::vector<long long>*> valVlonglongs;
  std::unordered_map<TString, std::vector<float>*> valVfloats;
  std::unordered_map<TString, std::vector<double>*> valVdoubles;
  std::unordered_map<TString, std::vector<std::string>*> valVstrings;
  std::unordered_map<TString, std::vector<TString>*> valVTStrings;
  std::unordered_map<TString, std::vector<CMSLorentzVector>*> valVCMSLorentzVectors;

  std::unordered_map<TString, std::vector<std::vector<bool>>*> valVVbools;
  std::unordered_map<TString, std::vector<std::vector<short>>*> valVVshorts;
  std::unordered_map<TString, std::vector<std::vector<unsigned int>>*> valVVuints;
  std::unordered_map<TString, std::vector<std::vector<int>>*> valVVints;
  std::unordered_map<TString, std::vector<std::vector<unsigned long>>*> valVVulongs;
  std::unordered_map<TString, std::vector<std::vector<long>>*> valVVlongs;
  std::unordered_map<TString, std::vector<std::vector<unsigned long long>>*> valVVulonglongs;
  std::unordered_map<TString, std::vector<std::vector<long long>>*> valVVlonglongs;
  std::unordered_map<TString, std::vector<std::vector<float>>*> valVVfloats;
  std::unordered_map<TString, std::vector<std::vector<double>>*> valVVdoubles;

  BranchType searchBranchType(TString branchname) const;

  void getValidBranchNamesWithoutAlias(TTree* t, std::vector<TString>& res) const;

  template<typename T> bool getBranchCIterator(TString branchname, typename std::unordered_map<TString, T>::iterator& it);
  template<typename T> bool getBranchCIterator(TString branchname, typename std::unordered_map<TString, T>::const_iterator& it) const;

  template<BranchType T> void resetBranch();
  virtual void resetBranches();

  template<BranchType T> void removeBranch(TString branchname);

public:
  BaseTree();
  BaseTree(const TString cinput, const TString treename, const TString failedtreename, const TString countersname);
  BaseTree(const TString treename); // Output constructor
  virtual ~BaseTree();

  template<typename T> bool bookBranch(TString branchname, T valdef);
  template<BranchType T> bool bookBranch(TString branchname);

  template<typename T> bool putBranch(TString branchname, T valdef);
  template<BranchType T> bool putBranch(TString branchname);

  TFile* getInputFile();
  TTree* getSelectedTree();
  TTree* getFailedTree();
  TFile const* getInputFile() const;
  TTree const* getSelectedTree() const;
  TTree const* getFailedTree() const;

  bool getSelectedEvent(int ev);
  bool getFailedEvent(int ev);
  bool getEvent(int ev);
  void refreshCurrentEvent();

  int getSelectedNEvents() const;
  int getFailedNEvents() const;
  int getNEvents() const;

  template<typename T> void getVal(TString branchname, T& val) const;
  template<typename T> void setVal(TString branchname, T const& val);
  template<typename T> void getValRef(TString branchname, T*& val) const;
  template<typename T> void getValRef(TString branchname, T*& val);

  void silenceUnused();
  void releaseBranch(TString branchname);

  bool isValid() const;
  bool branchExists(TString branchname, BranchType* type=nullptr);

  virtual bool isValidEvent() const;

  void fill();
  void writeToFile(TFile* file);

  virtual void print() const;

  static void writeSimpleEntries(std::vector<SimpleEntry>::iterator const& vecBegin, std::vector<SimpleEntry>::iterator const& vecEnd, BaseTree* const& tree);

};


#endif
