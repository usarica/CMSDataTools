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

  template<typename T, typename U> bool getUnorderedMapIterator(T const& theKey, std::unordered_map<T, U> const& theMap, typename std::unordered_map<T, U>::const_iterator& it);
  template<typename T, typename U> bool getUnorderedMapIterator(T const& theKey, std::unordered_map<T, U>& theMap, typename std::unordered_map<T, U>::iterator& it);

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

template<typename T, typename U> bool HelperFunctions::getUnorderedMapIterator(T const& theKey, std::unordered_map<T, U> const& theMap, typename std::unordered_map<T, U>::const_iterator& it){
  it = theMap.find(theKey);
  return (it!=theMap.cend());
}
template<typename T, typename U> bool HelperFunctions::getUnorderedMapIterator(T const& theKey, std::unordered_map<T, U>& theMap, typename std::unordered_map<T, U>::iterator& it){
  it = theMap.find(theKey);
  return (it!=theMap.end());
}

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
