#ifndef SIMPLEENTRY_H
#define SIMPLEENTRY_H

#include <vector>
#include <unordered_map>
#include "TString.h"
#include "TTree.h"
#include "HelperFunctionsCore.h"


struct SimpleEntry{
  int id;
  float trackingval;
  float weight;

  std::vector<float> recoval;

  std::unordered_map<TString, bool> namedbools;
  std::unordered_map<TString, short> namedshorts;
  std::unordered_map<TString, unsigned int> nameduints;
  std::unordered_map<TString, int> namedints;
  std::unordered_map<TString, unsigned long> namedulongs;
  std::unordered_map<TString, long> namedlongs;
  std::unordered_map<TString, long long> namedlonglongs;
  std::unordered_map<TString, float> namedfloats;
  std::unordered_map<TString, double> nameddoubles;

  std::unordered_map<TString, std::vector<bool>> namedVbools;
  std::unordered_map<TString, std::vector<short>> namedVshorts;
  std::unordered_map<TString, std::vector<unsigned int>> namedVuints;
  std::unordered_map<TString, std::vector<int>> namedVints;
  std::unordered_map<TString, std::vector<unsigned long>> namedVulongs;
  std::unordered_map<TString, std::vector<long>> namedVlongs;
  std::unordered_map<TString, std::vector<long long>> namedVlonglongs;
  std::unordered_map<TString, std::vector<float>> namedVfloats;
  std::unordered_map<TString, std::vector<double>> namedVdoubles;

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
