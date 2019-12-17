#ifndef BASETREE_H
#define BASETREE_H

#include <vector>
#include <fstream>
#include <cstdlib>
#include <unordered_map>
#include "StdExtensions.h"
#include "SimpleEntry.h"
#include "TFile.h"
#include "TBits.h"
#include "TTree.h"
#include "TH1F.h"
#include "CMSLorentzVector.h"
#include "AnalysisDataTypes.hh"


class BaseTree{
protected:
  static bool robustSaveWrite;

public:

#define SIMPLE_DATA_INPUT_DIRECTIVE(name, type, default_value) BranchType_##name##_t,
#define VECTOR_DATA_INPUT_DIRECTIVE(name, type) BranchType_v##name##_t,
#define DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(name, type) BranchType_vv##name##_t,
  enum BranchType{
    SIMPLE_DATA_INPUT_DIRECTIVES
    VECTOR_DATA_INPUT_DIRECTIVES
    DOUBLEVECTOR_DATA_INPUT_DIRECTIVES
    BranchType_unknown_t
  };
#undef SIMPLE_DATA_INPUT_DIRECTIVE
#undef VECTOR_DATA_INPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_INPUT_DIRECTIVE

  TString sampleIdentifier;

protected:
  TFile* finput;
  TTree* tree;
  TTree* failedtree;
  TH1F* hCounters;
  bool valid;
  const bool receiver;
  bool acquireTreePossession; // If true, deletes trees and histograms upon destruction.

  int currentEvent;
  TTree* currentTree;

#define SIMPLE_DATA_INPUT_DIRECTIVE(name, type, default_value) std::unordered_map<TString, std::pair<type, type>*> val##name##s;
#define VECTOR_DATA_INPUT_DIRECTIVE(name, type) std::unordered_map<TString, type*> valV##name##s;
#define DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(name, type) std::unordered_map<TString, type*> valVV##name##s;

  SIMPLE_DATA_INPUT_DIRECTIVES
  VECTOR_DATA_INPUT_DIRECTIVES
  DOUBLEVECTOR_DATA_INPUT_DIRECTIVES

#undef SIMPLE_DATA_INPUT_DIRECTIVE
#undef VECTOR_DATA_INPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_INPUT_DIRECTIVE

  BranchType searchBranchType(TString branchname) const;

  void getValidBranchNamesWithoutAlias(TTree* t, std::vector<TString>& res, bool check_linked) const;

  template<typename T> bool getBranchCIterator(TString branchname, typename std::unordered_map<TString, T>::iterator& it);
  template<typename T> bool getBranchCIterator(TString branchname, typename std::unordered_map<TString, T>::const_iterator& it) const;

  template<BranchType T> void resetBranch();
  virtual void resetBranches();

  template<BranchType T> void removeBranch(TString branchname);

public:
  BaseTree();
  BaseTree(const TString cinput, const TString treename, const TString failedtreename, const TString countersname);
  BaseTree(const TString treename); // Output constructor
  BaseTree(TFile* finput_, TTree* tree_, TTree* failedtree_, TH1F* hCounters_, bool receiver_override); // Mixed definition constructor
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

  // Overloads of getNEvents
  virtual unsigned int getNGenNoPU();
  virtual float getNGenWithPU();

  template<typename T> void getVal(TString branchname, T& val) const;
  template<typename T> void setVal(TString branchname, T const& val);
  template<typename T> void getValRef(TString branchname, T*& val) const;
  template<typename T> void getValRef(TString branchname, T*& val);

  void silenceUnused();
  void unmuteAllBranches();
  void releaseBranch(TString branchname);

  void setAcquireTreePossession(bool);

  void setAutoSave(Long64_t fsave);
  long long doAutoSave(const char* opts="");

  void setAutoFlush(Long64_t fflush);
  int doFlushBaskets();

  bool isValid() const;
  bool branchExists(TString branchname, BranchType* type=nullptr);
  void getValidBranchNamesWithoutAlias(std::vector<TString>& res, bool check_linked) const;

  virtual bool isValidEvent() const;

  void fill();
  void writeToFile(TFile* file);

  virtual void print() const;

  static void setRobustSaveWrite(bool flag);
  static void writeSimpleEntries(std::vector<SimpleEntry>::iterator const& vecBegin, std::vector<SimpleEntry>::iterator const& vecEnd, BaseTree* const& tree);

};


#endif
