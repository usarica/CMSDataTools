#ifndef BASETREE_H
#define BASETREE_H

#include "SampleHelpers.h"
#include "Samples.h"
#include "TSystem.h"
#include "TDirectory.h"


class BaseTree{
public:
  enum BranchType{
    BranchType_bool_t,
    BranchType_short_t,
    BranchType_uint_t,
    BranchType_int_t,
    BranchType_float_t,
    BranchType_double_t,

    BranchType_vbool_t,
    BranchType_vshort_t,
    BranchType_vuint_t,
    BranchType_vint_t,
    BranchType_vfloat_t,
    BranchType_vdouble_t,

    BranchType_unknown_t
  };

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
  std::unordered_map<TString, std::pair<float, float>*> valfloats;
  std::unordered_map<TString, std::pair<double, double>*> valdoubles;

  std::unordered_map<TString, std::vector<bool>*> valVbools;
  std::unordered_map<TString, std::vector<short>*> valVshorts;
  std::unordered_map<TString, std::vector<unsigned int>*> valVuints;
  std::unordered_map<TString, std::vector<int>*> valVints;
  std::unordered_map<TString, std::vector<float>*> valVfloats;
  std::unordered_map<TString, std::vector<double>*> valVdoubles;

  BranchType searchBranchType(TString branchname) const;

  template<typename T> bool getBranchCIterator(TString branchname, typename std::unordered_map<TString, T>::iterator& it);
  template<typename T> bool getBranchCIterator(TString branchname, typename std::unordered_map<TString, T>::const_iterator& it) const;

  template<BranchType T> void resetBranch();
  void resetBranches();

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

  bool getSelectedEvent(int ev);
  bool getFailedEvent(int ev);
  bool getEvent(int ev);
  void refreshCurrentEvent();

  int getSelectedNEvents();
  int getFailedNEvents();
  int getNEvents();

  template<typename T> void getVal(TString branchname, T& val) const;
  template<typename T> void setVal(TString branchname, T const& val);
  template<typename T> void getValRef(TString branchname, T*& val) const;

  void silenceUnused();
  void releaseBranch(TString branchname);

  bool isValid() const;
  bool branchExists(TString branchname, BranchType* type=nullptr);

  virtual bool isValidEvent() const;

  void fill();
  void writeToFile(TFile* file);

  static void writeSimpleEntries(std::vector<SimpleEntry>::iterator const& vecBegin, std::vector<SimpleEntry>::iterator const& vecEnd, BaseTree* const& tree);

};


#endif
