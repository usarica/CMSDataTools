#ifndef SIMPLEENTRY_H
#define SIMPLEENTRY_H

#include <vector>
#include <unordered_map>
#include "TString.h"
#include "StdExtensions.h"

struct SimpleEntry{
  int id;
  float trackingval;
  std::vector<float> recoval;
  float weight;

  std::unordered_map<TString, unsigned int> nameduints;
  std::unordered_map<TString, short> namedshorts;
  std::unordered_map<TString, int> namedints;
  std::unordered_map<TString, float> namedfloats;
  std::unordered_map<TString, double> nameddoubles;

  SimpleEntry();
  SimpleEntry(int id_, float trackingval_, std::vector<float> recoval_, float weight_=1);

  bool operator != (const SimpleEntry& other)const;
  bool operator == (const SimpleEntry& other)const;
  bool operator > (const SimpleEntry& other)const;
  bool operator >= (const SimpleEntry& other)const;
  bool operator < (const SimpleEntry& other)const;
  bool operator <= (const SimpleEntry& other)const;

  void setNamedVal(TString strname, unsigned int& val);
  void setNamedVal(TString strname, short& val);
  void setNamedVal(TString strname, int& val);
  void setNamedVal(TString strname, float& val);
  void setNamedVal(TString strname, double& val);

  void getNamedVal(TString strname, unsigned int& val);
  void getNamedVal(TString strname, short& val);
  void getNamedVal(TString strname, int& val);
  void getNamedVal(TString strname, float& val);
  void getNamedVal(TString strname, double& val);

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
