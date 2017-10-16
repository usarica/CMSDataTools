#ifndef SIMPLEENTRY_HPP
#define SIMPLEENTRY_HPP

#include "SimpleEntry.h"


template<> void SimpleEntry::setNamedVal<bool>(TString strname, bool const& val){
  typedef bool V;
  std::unordered_map<TString, V>::iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, namedbools, it)) it->second = val;
}
template<> void SimpleEntry::getNamedVal<bool>(TString strname, bool& val) const{
  typedef bool V;
  std::unordered_map<TString, V>::const_iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, namedbools, it)) val = it->second;
}
template<> void SimpleEntry::setNamedVal<short>(TString strname, short const& val){
  typedef short V;
  std::unordered_map<TString, V>::iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, namedshorts, it)) it->second = val;
}
template<> void SimpleEntry::getNamedVal<short>(TString strname, short& val) const{
  typedef short V;
  std::unordered_map<TString, V>::const_iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, namedshorts, it)) val = it->second;
}
template<> void SimpleEntry::setNamedVal<unsigned int>(TString strname, unsigned int const& val){
  typedef unsigned int V;
  std::unordered_map<TString, V>::iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, nameduints, it)) it->second = val;
}
template<> void SimpleEntry::getNamedVal<unsigned int>(TString strname, unsigned int& val) const{
  typedef unsigned int V;
  std::unordered_map<TString, V>::const_iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, nameduints, it)) val = it->second;
}
template<> void SimpleEntry::setNamedVal<int>(TString strname, int const& val){
  typedef int V;
  std::unordered_map<TString, V>::iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, namedints, it)) it->second = val;
}
template<> void SimpleEntry::getNamedVal<int>(TString strname, int& val) const{
  typedef int V;
  std::unordered_map<TString, V>::const_iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, namedints, it)) val = it->second;
}
template<> void SimpleEntry::setNamedVal<float>(TString strname, float const& val){
  typedef float V;
  std::unordered_map<TString, V>::iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, namedfloats, it)) it->second = val;
}
template<> void SimpleEntry::getNamedVal<float>(TString strname, float& val) const{
  typedef float V;
  std::unordered_map<TString, V>::const_iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, namedfloats, it)) val = it->second;
}
template<> void SimpleEntry::setNamedVal<double>(TString strname, double const& val){
  typedef double V;
  std::unordered_map<TString, V>::iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, nameddoubles, it)) it->second = val;
}
template<> void SimpleEntry::getNamedVal<double>(TString strname, double& val) const{
  typedef double V;
  std::unordered_map<TString, V>::const_iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, nameddoubles, it)) val = it->second;
}

template<> void SimpleEntry::setNamedVal<std::vector<bool>>(TString strname, std::vector<bool> const& val){
  typedef std::vector<bool> V;
  std::unordered_map<TString, V>::iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, namedVbools, it)) it->second = val;
}
template<> void SimpleEntry::getNamedVal<std::vector<bool>>(TString strname, std::vector<bool>& val) const{
  typedef std::vector<bool> V;
  std::unordered_map<TString, V>::const_iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, namedVbools, it)) val = it->second;
}
template<> void SimpleEntry::setNamedVal<std::vector<short>>(TString strname, std::vector<short> const& val){
  typedef std::vector<short> V;
  std::unordered_map<TString, V>::iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, namedVshorts, it)) it->second = val;
}
template<> void SimpleEntry::getNamedVal<std::vector<short>>(TString strname, std::vector<short>& val) const{
  typedef std::vector<short> V;
  std::unordered_map<TString, V>::const_iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, namedVshorts, it)) val = it->second;
}
template<> void SimpleEntry::setNamedVal<std::vector<unsigned int>>(TString strname, std::vector<unsigned int> const& val){
  typedef std::vector<unsigned int> V;
  std::unordered_map<TString, V>::iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, namedVuints, it)) it->second = val;
}
template<> void SimpleEntry::getNamedVal<std::vector<unsigned int>>(TString strname, std::vector<unsigned int>& val) const{
  typedef std::vector<unsigned int> V;
  std::unordered_map<TString, V>::const_iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, namedVuints, it)) val = it->second;
}
template<> void SimpleEntry::setNamedVal<std::vector<int>>(TString strname, std::vector<int> const& val){
  typedef std::vector<int> V;
  std::unordered_map<TString, V>::iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, namedVints, it)) it->second = val;
}
template<> void SimpleEntry::getNamedVal<std::vector<int>>(TString strname, std::vector<int>& val) const{
  typedef std::vector<int> V;
  std::unordered_map<TString, V>::const_iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, namedVints, it)) val = it->second;
}
template<> void SimpleEntry::setNamedVal<std::vector<float>>(TString strname, std::vector<float> const& val){
  typedef std::vector<float> V;
  std::unordered_map<TString, V>::iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, namedVfloats, it)) it->second = val;
}
template<> void SimpleEntry::getNamedVal<std::vector<float>>(TString strname, std::vector<float>& val) const{
  typedef std::vector<float> V;
  std::unordered_map<TString, V>::const_iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, namedVfloats, it)) val = it->second;
}
template<> void SimpleEntry::setNamedVal<std::vector<double>>(TString strname, std::vector<double> const& val){
  typedef std::vector<double> V;
  std::unordered_map<TString, V>::iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, namedVdoubles, it)) it->second = val;
}
template<> void SimpleEntry::getNamedVal<std::vector<double>>(TString strname, std::vector<double>& val) const{
  typedef std::vector<double> V;
  std::unordered_map<TString, V>::const_iterator it;
  if (HelperFunctions::getUnorderedMapIterator<V>(strname, namedVdoubles, it)) val = it->second;
}


#endif
