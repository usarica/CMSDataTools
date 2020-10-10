#ifndef HELPERFUNCTIONSCORE_H
#define HELPERFUNCTIONSCORE_H

#include <string>
#include <sstream>
#include <vector>
#include <utility>
#include <algorithm>
#include <cctype>
#include <unordered_map>
#include "StdExtensions.h"
#include "TUtilHelpers.hh"
#include "CMSLorentzVector.h"


// Behold the power of the member checker!
// From https://stackoverflow.com/questions/1005476/how-to-detect-whether-there-is-a-specific-member-variable-in-class
// You need to use the first one to construct the structs and the second one in the if-statement
#define CONSTRUCT_MEMBER_CHECKER(member_type, member) \
template<typename class_type, typename = void> struct has_##member : std::false_type {}; \
template<typename class_type> struct has_##member<class_type, decltype((void)class_type::member, member_type()) > : std::true_type {};
#define HAS_MEMBER(class_type, member_type, member) has_##member<class_type, member_type>::value
// Usage example:
// [Outside of function definitions]
// CONSTRUCT_MEMBER_CHECKER(int, id)
// [Inside the function]
// HAS_MEMBER(SimpleEntry, int, id) -> returns true
// HAS_MEMBER(int, int, id) -> returns false


namespace HelperFunctions{

  template<typename T> bool getUnorderedMapIterator(TString name, const std::unordered_map<TString, T>& theMap, typename std::unordered_map<TString, T>::const_iterator& it);
  template<typename T> bool getUnorderedMapIterator(TString name, std::unordered_map<TString, T>& theMap, typename std::unordered_map<TString, T>::iterator& it);

  template<typename T> void resetPointer(T*& ptr);

  template<typename T, typename U> bool replaceString(T& strinput, U strTakeOut, U strPutIn);
  template<> bool replaceString<TString, const TString>(TString& strinput, const TString strTakeOut, const TString strPutIn);
  template<> bool replaceString<TString, const char*>(TString& strinput, const char* strTakeOut, const char* strPutIn);
  template<> bool replaceString<std::string, const std::string>(std::string& strinput, const std::string strTakeOut, const std::string strPutIn);
  template<> bool replaceString<std::string, const char*>(std::string& strinput, const char* strTakeOut, const char* strPutIn);

  template<typename T> void lstrip(T& str, const char* chars=nullptr);
  template<typename T> void rstrip(T& str, const char* chars=nullptr);
  template<> void lstrip<std::string>(std::string& str, const char* chars);
  template<> void lstrip<TString>(TString& str, const char* chars);
  template<> void rstrip<std::string>(std::string& str, const char* chars);
  template<> void rstrip<TString>(TString& str, const char* chars);

  template<typename T> void castStringToValue(std::string const& name, T& val);
  template<> void castStringToValue(std::string const& name, bool& val);
  template<typename T> void castStringToValue(TString const& name, T& val);
  template<typename T> void castStringToValue(const char* name, T& val);

  template<typename T> void lowercase(T const& name, T& val);
  template<> void lowercase(std::string const& name, std::string& val);
  template<> void lowercase(TString const& name, TString& val);
  template<> void lowercase(const char* const& name, const char*& val);

  // Value-to-string casts
  template<typename T> std::string castValueToString(T const& val, unsigned short max_decimals=4);

}

template<typename T> void HelperFunctions::resetPointer(T*& ptr){ delete ptr; ptr=nullptr; }

template<typename T> bool HelperFunctions::getUnorderedMapIterator(TString name, const std::unordered_map<TString, T>& theMap, typename std::unordered_map<TString, T>::const_iterator& it){
  it = theMap.find(name);
  return (it!=theMap.cend());
}
template<typename T> bool HelperFunctions::getUnorderedMapIterator(TString name, std::unordered_map<TString, T>& theMap, typename std::unordered_map<TString, T>::iterator& it){
  it = theMap.find(name);
  return (it!=theMap.end());
}

template bool HelperFunctions::getUnorderedMapIterator<bool>(TString name, const std::unordered_map<TString, bool>& theMap, std::unordered_map<TString, bool>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<bool>(TString name, std::unordered_map<TString, bool>& theMap, std::unordered_map<TString, bool>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<short>(TString name, const std::unordered_map<TString, short>& theMap, std::unordered_map<TString, short>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<short>(TString name, std::unordered_map<TString, short>& theMap, std::unordered_map<TString, short>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<unsigned int>(TString name, const std::unordered_map<TString, unsigned int>& theMap, std::unordered_map<TString, unsigned int>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<unsigned int>(TString name, std::unordered_map<TString, unsigned int>& theMap, std::unordered_map<TString, unsigned int>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<int>(TString name, const std::unordered_map<TString, int>& theMap, std::unordered_map<TString, int>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<int>(TString name, std::unordered_map<TString, int>& theMap, std::unordered_map<TString, int>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<unsigned long>(TString name, const std::unordered_map<TString, unsigned long>& theMap, std::unordered_map<TString, unsigned long>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<unsigned long>(TString name, std::unordered_map<TString, unsigned long>& theMap, std::unordered_map<TString, unsigned long>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<long>(TString name, const std::unordered_map<TString, long>& theMap, std::unordered_map<TString, long>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<long>(TString name, std::unordered_map<TString, long>& theMap, std::unordered_map<TString, long>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<unsigned long long>(TString name, const std::unordered_map<TString, unsigned long long>& theMap, std::unordered_map<TString, unsigned long long>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<unsigned long long>(TString name, std::unordered_map<TString, unsigned long long>& theMap, std::unordered_map<TString, unsigned long long>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<long long>(TString name, const std::unordered_map<TString, long long>& theMap, std::unordered_map<TString, long long>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<long long>(TString name, std::unordered_map<TString, long long>& theMap, std::unordered_map<TString, long long>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<float>(TString name, const std::unordered_map<TString, float>& theMap, std::unordered_map<TString, float>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<float>(TString name, std::unordered_map<TString, float>& theMap, std::unordered_map<TString, float>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<double>(TString name, const std::unordered_map<TString, double>& theMap, std::unordered_map<TString, double>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<double>(TString name, std::unordered_map<TString, double>& theMap, std::unordered_map<TString, double>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<CMSLorentzVector>(TString name, const std::unordered_map<TString, CMSLorentzVector>& theMap, std::unordered_map<TString, CMSLorentzVector>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<CMSLorentzVector>(TString name, std::unordered_map<TString, CMSLorentzVector>& theMap, std::unordered_map<TString, CMSLorentzVector>::iterator& it);

template bool HelperFunctions::getUnorderedMapIterator<bool*>(TString name, const std::unordered_map<TString, bool*>& theMap, std::unordered_map<TString, bool*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<bool*>(TString name, std::unordered_map<TString, bool*>& theMap, std::unordered_map<TString, bool*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<short*>(TString name, const std::unordered_map<TString, short*>& theMap, std::unordered_map<TString, short*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<short*>(TString name, std::unordered_map<TString, short*>& theMap, std::unordered_map<TString, short*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<unsigned int*>(TString name, const std::unordered_map<TString, unsigned int*>& theMap, std::unordered_map<TString, unsigned int*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<unsigned int*>(TString name, std::unordered_map<TString, unsigned int*>& theMap, std::unordered_map<TString, unsigned int*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<int*>(TString name, const std::unordered_map<TString, int*>& theMap, std::unordered_map<TString, int*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<int*>(TString name, std::unordered_map<TString, int*>& theMap, std::unordered_map<TString, int*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<unsigned long*>(TString name, const std::unordered_map<TString, unsigned long*>& theMap, std::unordered_map<TString, unsigned long*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<unsigned long*>(TString name, std::unordered_map<TString, unsigned long*>& theMap, std::unordered_map<TString, unsigned long*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<long*>(TString name, const std::unordered_map<TString, long*>& theMap, std::unordered_map<TString, long*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<long*>(TString name, std::unordered_map<TString, long*>& theMap, std::unordered_map<TString, long*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<unsigned long long*>(TString name, const std::unordered_map<TString, unsigned long long*>& theMap, std::unordered_map<TString, unsigned long long*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<unsigned long long*>(TString name, std::unordered_map<TString, unsigned long long*>& theMap, std::unordered_map<TString, unsigned long long*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<long long*>(TString name, const std::unordered_map<TString, long long*>& theMap, std::unordered_map<TString, long long*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<long long*>(TString name, std::unordered_map<TString, long long*>& theMap, std::unordered_map<TString, long long*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<float*>(TString name, const std::unordered_map<TString, float*>& theMap, std::unordered_map<TString, float*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<float*>(TString name, std::unordered_map<TString, float*>& theMap, std::unordered_map<TString, float*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<double*>(TString name, const std::unordered_map<TString, double*>& theMap, std::unordered_map<TString, double*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<double*>(TString name, std::unordered_map<TString, double*>& theMap, std::unordered_map<TString, double*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<CMSLorentzVector*>(TString name, const std::unordered_map<TString, CMSLorentzVector*>& theMap, std::unordered_map<TString, CMSLorentzVector*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<CMSLorentzVector*>(TString name, std::unordered_map<TString, CMSLorentzVector*>& theMap, std::unordered_map<TString, CMSLorentzVector*>::iterator& it);

template bool HelperFunctions::getUnorderedMapIterator<std::vector<bool>>(TString name, const std::unordered_map<TString, std::vector<bool>>& theMap, std::unordered_map<TString, std::vector<bool>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<bool>>(TString name, std::unordered_map<TString, std::vector<bool>>& theMap, std::unordered_map<TString, std::vector<bool>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<short>>(TString name, const std::unordered_map<TString, std::vector<short>>& theMap, std::unordered_map<TString, std::vector<short>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<short>>(TString name, std::unordered_map<TString, std::vector<short>>& theMap, std::unordered_map<TString, std::vector<short>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<unsigned int>>(TString name, const std::unordered_map<TString, std::vector<unsigned int>>& theMap, std::unordered_map<TString, std::vector<unsigned int>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<unsigned int>>(TString name, std::unordered_map<TString, std::vector<unsigned int>>& theMap, std::unordered_map<TString, std::vector<unsigned int>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<int>>(TString name, const std::unordered_map<TString, std::vector<int>>& theMap, std::unordered_map<TString, std::vector<int>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<int>>(TString name, std::unordered_map<TString, std::vector<int>>& theMap, std::unordered_map<TString, std::vector<int>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<unsigned long>>(TString name, const std::unordered_map<TString, std::vector<unsigned long>>& theMap, std::unordered_map<TString, std::vector<unsigned long>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<unsigned long>>(TString name, std::unordered_map<TString, std::vector<unsigned long>>& theMap, std::unordered_map<TString, std::vector<unsigned long>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<long>>(TString name, const std::unordered_map<TString, std::vector<long>>& theMap, std::unordered_map<TString, std::vector<long>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<long>>(TString name, std::unordered_map<TString, std::vector<long>>& theMap, std::unordered_map<TString, std::vector<long>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<unsigned long long>>(TString name, const std::unordered_map<TString, std::vector<unsigned long long>>& theMap, std::unordered_map<TString, std::vector<unsigned long long>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<unsigned long long>>(TString name, std::unordered_map<TString, std::vector<unsigned long long>>& theMap, std::unordered_map<TString, std::vector<unsigned long long>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<long long>>(TString name, const std::unordered_map<TString, std::vector<long long>>& theMap, std::unordered_map<TString, std::vector<long long>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<long long>>(TString name, std::unordered_map<TString, std::vector<long long>>& theMap, std::unordered_map<TString, std::vector<long long>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<float>>(TString name, const std::unordered_map<TString, std::vector<float>>& theMap, std::unordered_map<TString, std::vector<float>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<float>>(TString name, std::unordered_map<TString, std::vector<float>>& theMap, std::unordered_map<TString, std::vector<float>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<double>>(TString name, const std::unordered_map<TString, std::vector<double>>& theMap, std::unordered_map<TString, std::vector<double>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<double>>(TString name, std::unordered_map<TString, std::vector<double>>& theMap, std::unordered_map<TString, std::vector<double>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<CMSLorentzVector>>(TString name, const std::unordered_map<TString, std::vector<CMSLorentzVector>>& theMap, std::unordered_map<TString, std::vector<CMSLorentzVector>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<CMSLorentzVector>>(TString name, std::unordered_map<TString, std::vector<CMSLorentzVector>>& theMap, std::unordered_map<TString, std::vector<CMSLorentzVector>>::iterator& it);

template bool HelperFunctions::getUnorderedMapIterator<std::vector<bool>*>(TString name, const std::unordered_map<TString, std::vector<bool>*>& theMap, std::unordered_map<TString, std::vector<bool>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<bool>*>(TString name, std::unordered_map<TString, std::vector<bool>*>& theMap, std::unordered_map<TString, std::vector<bool>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<short>*>(TString name, const std::unordered_map<TString, std::vector<short>*>& theMap, std::unordered_map<TString, std::vector<short>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<short>*>(TString name, std::unordered_map<TString, std::vector<short>*>& theMap, std::unordered_map<TString, std::vector<short>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<unsigned int>*>(TString name, const std::unordered_map<TString, std::vector<unsigned int>*>& theMap, std::unordered_map<TString, std::vector<unsigned int>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<unsigned int>*>(TString name, std::unordered_map<TString, std::vector<unsigned int>*>& theMap, std::unordered_map<TString, std::vector<unsigned int>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<int>*>(TString name, const std::unordered_map<TString, std::vector<int>*>& theMap, std::unordered_map<TString, std::vector<int>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<int>*>(TString name, std::unordered_map<TString, std::vector<int>*>& theMap, std::unordered_map<TString, std::vector<int>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<unsigned long>*>(TString name, const std::unordered_map<TString, std::vector<unsigned long>*>& theMap, std::unordered_map<TString, std::vector<unsigned long>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<unsigned long>*>(TString name, std::unordered_map<TString, std::vector<unsigned long>*>& theMap, std::unordered_map<TString, std::vector<unsigned long>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<long>*>(TString name, const std::unordered_map<TString, std::vector<long>*>& theMap, std::unordered_map<TString, std::vector<long>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<long>*>(TString name, std::unordered_map<TString, std::vector<long>*>& theMap, std::unordered_map<TString, std::vector<long>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<unsigned long long>*>(TString name, const std::unordered_map<TString, std::vector<unsigned long long>*>& theMap, std::unordered_map<TString, std::vector<unsigned long long>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<unsigned long long>*>(TString name, std::unordered_map<TString, std::vector<unsigned long long>*>& theMap, std::unordered_map<TString, std::vector<unsigned long long>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<long long>*>(TString name, const std::unordered_map<TString, std::vector<long long>*>& theMap, std::unordered_map<TString, std::vector<long long>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<long long>*>(TString name, std::unordered_map<TString, std::vector<long long>*>& theMap, std::unordered_map<TString, std::vector<long long>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<float>*>(TString name, const std::unordered_map<TString, std::vector<float>*>& theMap, std::unordered_map<TString, std::vector<float>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<float>*>(TString name, std::unordered_map<TString, std::vector<float>*>& theMap, std::unordered_map<TString, std::vector<float>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<double>*>(TString name, const std::unordered_map<TString, std::vector<double>*>& theMap, std::unordered_map<TString, std::vector<double>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<double>*>(TString name, std::unordered_map<TString, std::vector<double>*>& theMap, std::unordered_map<TString, std::vector<double>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<CMSLorentzVector>*>(TString name, const std::unordered_map<TString, std::vector<CMSLorentzVector>*>& theMap, std::unordered_map<TString, std::vector<CMSLorentzVector>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<CMSLorentzVector>*>(TString name, std::unordered_map<TString, std::vector<CMSLorentzVector>*>& theMap, std::unordered_map<TString, std::vector<CMSLorentzVector>*>::iterator& it);

template bool HelperFunctions::getUnorderedMapIterator<std::pair<bool, bool>>(TString name, const std::unordered_map<TString, std::pair<bool, bool>>& theMap, std::unordered_map<TString, std::pair<bool, bool>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<bool, bool>>(TString name, std::unordered_map<TString, std::pair<bool, bool>>& theMap, std::unordered_map<TString, std::pair<bool, bool>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<short, short>>(TString name, const std::unordered_map<TString, std::pair<short, short>>& theMap, std::unordered_map<TString, std::pair<short, short>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<short, short>>(TString name, std::unordered_map<TString, std::pair<short, short>>& theMap, std::unordered_map<TString, std::pair<short, short>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<unsigned int, unsigned int>>(TString name, const std::unordered_map<TString, std::pair<unsigned int, unsigned int>>& theMap, std::unordered_map<TString, std::pair<unsigned int, unsigned int>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<unsigned int, unsigned int>>(TString name, std::unordered_map<TString, std::pair<unsigned int, unsigned int>>& theMap, std::unordered_map<TString, std::pair<unsigned int, unsigned int>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<int, int>>(TString name, const std::unordered_map<TString, std::pair<int, int>>& theMap, std::unordered_map<TString, std::pair<int, int>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<int, int>>(TString name, std::unordered_map<TString, std::pair<int, int>>& theMap, std::unordered_map<TString, std::pair<int, int>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<unsigned long, unsigned long>>(TString name, const std::unordered_map<TString, std::pair<unsigned long, unsigned long>>& theMap, std::unordered_map<TString, std::pair<unsigned long, unsigned long>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<unsigned long, unsigned long>>(TString name, std::unordered_map<TString, std::pair<unsigned long, unsigned long>>& theMap, std::unordered_map<TString, std::pair<unsigned long, unsigned long>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<long, long>>(TString name, const std::unordered_map<TString, std::pair<long, long>>& theMap, std::unordered_map<TString, std::pair<long, long>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<long, long>>(TString name, std::unordered_map<TString, std::pair<long, long>>& theMap, std::unordered_map<TString, std::pair<long, long>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<unsigned long long, unsigned long long>>(TString name, const std::unordered_map<TString, std::pair<unsigned long long, unsigned long long>>& theMap, std::unordered_map<TString, std::pair<unsigned long long, unsigned long long>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<unsigned long long, unsigned long long>>(TString name, std::unordered_map<TString, std::pair<unsigned long long, unsigned long long>>& theMap, std::unordered_map<TString, std::pair<unsigned long long, unsigned long long>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<long long, long long>>(TString name, const std::unordered_map<TString, std::pair<long long, long long>>& theMap, std::unordered_map<TString, std::pair<long long, long long>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<long long, long long>>(TString name, std::unordered_map<TString, std::pair<long long, long long>>& theMap, std::unordered_map<TString, std::pair<long long, long long>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<float, float>>(TString name, const std::unordered_map<TString, std::pair<float, float>>& theMap, std::unordered_map<TString, std::pair<float, float>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<float, float>>(TString name, std::unordered_map<TString, std::pair<float, float>>& theMap, std::unordered_map<TString, std::pair<float, float>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<double, double>>(TString name, const std::unordered_map<TString, std::pair<double, double>>& theMap, std::unordered_map<TString, std::pair<double, double>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<double, double>>(TString name, std::unordered_map<TString, std::pair<double, double>>& theMap, std::unordered_map<TString, std::pair<double, double>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<CMSLorentzVector, CMSLorentzVector>>(TString name, const std::unordered_map<TString, std::pair<CMSLorentzVector, CMSLorentzVector>>& theMap, std::unordered_map<TString, std::pair<CMSLorentzVector, CMSLorentzVector>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<CMSLorentzVector, CMSLorentzVector>>(TString name, std::unordered_map<TString, std::pair<CMSLorentzVector, CMSLorentzVector>>& theMap, std::unordered_map<TString, std::pair<CMSLorentzVector, CMSLorentzVector>>::iterator& it);

template bool HelperFunctions::getUnorderedMapIterator<std::pair<bool, bool>*>(TString name, const std::unordered_map<TString, std::pair<bool, bool>*>& theMap, std::unordered_map<TString, std::pair<bool, bool>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<bool, bool>*>(TString name, std::unordered_map<TString, std::pair<bool, bool>*>& theMap, std::unordered_map<TString, std::pair<bool, bool>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<short, short>*>(TString name, const std::unordered_map<TString, std::pair<short, short>*>& theMap, std::unordered_map<TString, std::pair<short, short>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<short, short>*>(TString name, std::unordered_map<TString, std::pair<short, short>*>& theMap, std::unordered_map<TString, std::pair<short, short>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<unsigned int, unsigned int>*>(TString name, const std::unordered_map<TString, std::pair<unsigned int, unsigned int>*>& theMap, std::unordered_map<TString, std::pair<unsigned int, unsigned int>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<unsigned int, unsigned int>*>(TString name, std::unordered_map<TString, std::pair<unsigned int, unsigned int>*>& theMap, std::unordered_map<TString, std::pair<unsigned int, unsigned int>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<int, int>*>(TString name, const std::unordered_map<TString, std::pair<int, int>*>& theMap, std::unordered_map<TString, std::pair<int, int>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<int, int>*>(TString name, std::unordered_map<TString, std::pair<int, int>*>& theMap, std::unordered_map<TString, std::pair<int, int>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<unsigned long, unsigned long>*>(TString name, const std::unordered_map<TString, std::pair<unsigned long, unsigned long>*>& theMap, std::unordered_map<TString, std::pair<unsigned long, unsigned long>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<unsigned long, unsigned long>*>(TString name, std::unordered_map<TString, std::pair<unsigned long, unsigned long>*>& theMap, std::unordered_map<TString, std::pair<unsigned long, unsigned long>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<long, long>*>(TString name, const std::unordered_map<TString, std::pair<long, long>*>& theMap, std::unordered_map<TString, std::pair<long, long>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<long, long>*>(TString name, std::unordered_map<TString, std::pair<long, long>*>& theMap, std::unordered_map<TString, std::pair<long, long>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<unsigned long long, unsigned long long>*>(TString name, const std::unordered_map<TString, std::pair<unsigned long long, unsigned long long>*>& theMap, std::unordered_map<TString, std::pair<unsigned long long, unsigned long long>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<unsigned long long, unsigned long long>*>(TString name, std::unordered_map<TString, std::pair<unsigned long long, unsigned long long>*>& theMap, std::unordered_map<TString, std::pair<unsigned long long, unsigned long long>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<long long, long long>*>(TString name, const std::unordered_map<TString, std::pair<long long, long long>*>& theMap, std::unordered_map<TString, std::pair<long long, long long>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<long long, long long>*>(TString name, std::unordered_map<TString, std::pair<long long, long long>*>& theMap, std::unordered_map<TString, std::pair<long long, long long>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<float, float>*>(TString name, const std::unordered_map<TString, std::pair<float, float>*>& theMap, std::unordered_map<TString, std::pair<float, float>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<float, float>*>(TString name, std::unordered_map<TString, std::pair<float, float>*>& theMap, std::unordered_map<TString, std::pair<float, float>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<double, double>*>(TString name, const std::unordered_map<TString, std::pair<double, double>*>& theMap, std::unordered_map<TString, std::pair<double, double>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<double, double>*>(TString name, std::unordered_map<TString, std::pair<double, double>*>& theMap, std::unordered_map<TString, std::pair<double, double>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<CMSLorentzVector, CMSLorentzVector>*>(TString name, const std::unordered_map<TString, std::pair<CMSLorentzVector, CMSLorentzVector>*>& theMap, std::unordered_map<TString, std::pair<CMSLorentzVector, CMSLorentzVector>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<CMSLorentzVector, CMSLorentzVector>*>(TString name, std::unordered_map<TString, std::pair<CMSLorentzVector, CMSLorentzVector>*>& theMap, std::unordered_map<TString, std::pair<CMSLorentzVector, CMSLorentzVector>*>::iterator& it);


template<typename T> void HelperFunctions::castStringToValue(std::string const& name, T& val){ std::stringstream ss(name); ss >> val; }
template<typename T> void HelperFunctions::castStringToValue(TString const& name, T& val){ std::string s(name.Data()); HelperFunctions::castStringToValue(s, val); }
template<typename T> void HelperFunctions::castStringToValue(const char* name, T& val){ std::string s(name); HelperFunctions::castStringToValue(s, val); }
template void HelperFunctions::castStringToValue(TString const& name, bool& val);
template void HelperFunctions::castStringToValue(const char* name, bool& val);
template void HelperFunctions::castStringToValue(std::string const& name, short& val);
template void HelperFunctions::castStringToValue(TString const& name, short& val);
template void HelperFunctions::castStringToValue(const char* name, short& val);
template void HelperFunctions::castStringToValue(std::string const& name, unsigned int& val);
template void HelperFunctions::castStringToValue(TString const& name, unsigned int& val);
template void HelperFunctions::castStringToValue(const char* name, unsigned int& val);
template void HelperFunctions::castStringToValue(std::string const& name, int& val);
template void HelperFunctions::castStringToValue(TString const& name, int& val);
template void HelperFunctions::castStringToValue(const char* name, int& val);
template void HelperFunctions::castStringToValue(std::string const& name, float& val);
template void HelperFunctions::castStringToValue(TString const& name, float& val);
template void HelperFunctions::castStringToValue(const char* name, float& val);
template void HelperFunctions::castStringToValue(std::string const& name, double& val);
template void HelperFunctions::castStringToValue(TString const& name, double& val);
template void HelperFunctions::castStringToValue(const char* name, double& val);

template<typename T> std::string HelperFunctions::castValueToString(T const& val, unsigned short max_decimals){
  double decimals = std::abs(val - T((int) val));
  if (decimals == 0.) return Form("%.0f", double(val));
  decimals += 1e-10; // Machine precision correction...
  unsigned short base10exponent = std::ceil(std::abs(std::log10(decimals)));
  unsigned short iprec = std::min(max_decimals, static_cast<unsigned short>(base10exponent+4));
  TString strprintf = Form("%s%u%s", "%.", iprec, "f");
  unsigned int remainder_prevtoLastDigit = static_cast<int>(decimals*std::pow(10., iprec+1)) % 5;
  double addval = (remainder_prevtoLastDigit==0 ? std::pow(10., -(iprec+1)) : 0.); // Form is smart enough to round 0.00006 to 0.0001, but 0.00005 becomes 0.0000.
  std::string res = Form(strprintf.Data(), static_cast<double>(val)+addval);
  while (res.back()=='0') res.pop_back();
  return res;
}

#endif
