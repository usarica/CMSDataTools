#ifndef SIMPLEENTRY_H
#define SIMPLEENTRY_H

#include <vector>
#include <unordered_map>
#include "TString.h"
#include "HelperFunctionsCore.h"


struct SimpleEntry{
  int id;
  float trackingval;
  std::vector<float> recoval;
  float weight;

  std::unordered_map<TString, bool> namedbools;
  std::unordered_map<TString, unsigned int> nameduints;
  std::unordered_map<TString, short> namedshorts;
  std::unordered_map<TString, int> namedints;
  std::unordered_map<TString, float> namedfloats;
  std::unordered_map<TString, double> nameddoubles;

  std::unordered_map<TString, std::vector<bool>> namedVbools;
  std::unordered_map<TString, std::vector<unsigned int>> namedVuints;
  std::unordered_map<TString, std::vector<short>> namedVshorts;
  std::unordered_map<TString, std::vector<int>> namedVints;
  std::unordered_map<TString, std::vector<float>> namedVfloats;
  std::unordered_map<TString, std::vector<double>> namedVdoubles;

  SimpleEntry();
  SimpleEntry(int id_, float trackingval_, float weight_=1);
  SimpleEntry(int id_, float trackingval_, std::vector<float> recoval_, float weight_=1);

  bool operator != (const SimpleEntry& other)const;
  bool operator == (const SimpleEntry& other)const;
  bool operator > (const SimpleEntry& other)const;
  bool operator >= (const SimpleEntry& other)const;
  bool operator < (const SimpleEntry& other)const;
  bool operator <= (const SimpleEntry& other)const;

  template<typename T> void setNamedVal(TString strname, T const& val);
  template<typename T> void getNamedVal(TString strname, T& val) const;

  static void cropByTrueVal(std::vector<SimpleEntry>& vec, float minval, float maxval);
  void print();

};


struct ExtBin{
  double binlow;
  double binhigh;
  std::vector<SimpleEntry> collection;
  void addEvent(SimpleEntry evt);
};



#endif
