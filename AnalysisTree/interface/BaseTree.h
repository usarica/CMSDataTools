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
  static bool doRobustInputCheck;

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
  std::vector<TTree*> treelist;
  TTree* tree;
  TTree* failedtree;
  TH1F* hCounters;
  bool valid;
  const bool receiver;
  bool acquireTreePossession; // If true, deletes trees and histograms upon destruction.
  bool isTChain;

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

  template<BranchType T> void removeBranch(TString branchname);

  // [tree][file] is the order of output vector indices
  bool getValidFilesForTreeList(TString cinput, std::vector<TString> const& treenames, std::vector< std::vector<TString> >& res) const;

public:
  BaseTree();
  BaseTree(const TString cinput, const TString treename, const TString failedtreename, const TString countersname);
  BaseTree(const TString cinput, std::vector<TString> const& treenames, const TString countersname);
  BaseTree(const TString treename); // Output constructor
  BaseTree(TFile* finput_, TTree* tree_, TTree* failedtree_, TH1F* hCounters_, bool receiver_override); // Mixed definition constructor
  BaseTree(TFile* finput_, std::vector<TTree*> const& treelist_, TH1F* hCounters_, bool receiver_override); // Mixed definition constructor
  virtual ~BaseTree();

  template<typename T> bool bookBranch(TString branchname, T valdef);
  template<BranchType T> bool bookBranch(TString branchname);

  template<typename T> bool putBranch(TString branchname, T valdef);
  template<BranchType T> bool putBranch(TString branchname);

  virtual void resetBranches();

  TFile* getInputFile();
  TTree* getSelectedTree();
  TTree* getFailedTree();
  std::vector<TTree*>& getValidTrees();
  TFile const* getInputFile() const;
  TTree const* getSelectedTree() const;
  TTree const* getFailedTree() const;
  std::vector<TTree*> const& getValidTrees() const;

  bool getSelectedEvent(int ev);
  bool getFailedEvent(int ev);
  bool getEvent(int ev);
  void refreshCurrentEvent();

  // Update single branch for an event
  // Should work similar to getEvent, but for a single branch.
  bool updateBranch(int ev, TString const& bname, bool check_linked = true);

  // For caching purposes
  bool getCurrentEventInfo(TTree*& currentTree_, int& currentEvent_) const;
  bool isSameEvent(TTree* const& currentTree_, int const& currentEvent_) const;

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
  // This function below is for expert-use only!
  // The difference of the function below from the one above is that
  // the one above resets branch addresses
  // whereas the one below can be used to suppress branches temporarily.
  void muteAllBranchesExcept(std::vector<TString> const& bnames_excepted);
  // This function is to recover branches temporarily muted using muteAllBranchesExcept.
  void recoverMutedBranches(){ this->silenceUnused(); }

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
  void writeToDirectory(TDirectory* dir);

  virtual void print() const;

  static void setRobustSaveWrite(bool flag);
  static void setRobustInputCheck(bool flag);
  static void writeSimpleEntries(std::vector<SimpleEntry>::const_iterator const& vecBegin, std::vector<SimpleEntry>::const_iterator const& vecEnd, BaseTree* const& tree_, bool createBranches=true);

};


#endif
