#ifndef SIMPLEENTRY_H
#define SIMPLEENTRY_H

#include <vector>
#include <unordered_map>
#include "AnalysisDataTypes.hh"
#include "TBits.h"
#include "TString.h"
#include "TTree.h"
#include "HelperFunctionsCore.h"


struct SimpleEntry{

#define SIMPLE_DATA_OUTPUT_DIRECTIVE(name_t, type) std::unordered_map<TString, type> named##name_t##s;
#define VECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) std::unordered_map<TString, type> namedV##name_t##s;
#define DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) std::unordered_map<TString, type> namedVV##name_t##s;
  SIMPLE_DATA_OUTPUT_DIRECTIVES
  VECTOR_DATA_OUTPUT_DIRECTIVES
  DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVES
#undef SIMPLE_DATA_OUTPUT_DIRECTIVE
#undef VECTOR_DATA_OUTPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE

  int id;
  float trackingval;
  float weight;
  std::vector<float> recoval;

  SimpleEntry();
  SimpleEntry(int id_, float trackingval_, float weight_=1);
  SimpleEntry(int id_, float trackingval_, std::vector<float> recoval_, float weight_=1);
  SimpleEntry(SimpleEntry const& other);
  SimpleEntry(SimpleEntry&& other);

  void swap(SimpleEntry& other);
  SimpleEntry& operator=(const SimpleEntry& other);

  bool operator != (const SimpleEntry& other)const;
  bool operator == (const SimpleEntry& other)const;
  bool operator > (const SimpleEntry& other)const;
  bool operator >= (const SimpleEntry& other)const;
  bool operator < (const SimpleEntry& other)const;
  bool operator <= (const SimpleEntry& other)const;

  template<typename T> std::unordered_map<TString, T> const& getNamedMap() const;
  template<typename T> std::unordered_map<TString, T>& getNamedMap();

  template<typename T> void setNamedVal(TString strname, T const& val);
  template<typename T> void getNamedVal(TString strname, T& val) const;

  static void cropByTrueVal(std::vector<SimpleEntry>& vec, float minval, float maxval);
  static void writeToTree(std::vector<SimpleEntry>::const_iterator const& vecBegin, std::vector<SimpleEntry>::const_iterator const& vecEnd, TTree* const& tree);
  void print();

};


struct ExtBin{
  double binlow;
  double binhigh;
  std::vector<SimpleEntry> collection;
  void addEvent(SimpleEntry evt);
};



#endif
