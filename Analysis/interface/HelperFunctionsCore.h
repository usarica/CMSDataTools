#ifndef HELPERFUNCTIONSCORE_H
#define HELPERFUNCTIONSCORE_H

#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include "StdExtensions.h"


namespace HelperFunctions{

  template<typename T> bool getUnorderedMapIterator(TString name, const std::unordered_map<TString, T>& theMap, typename std::unordered_map<TString, T>::const_iterator& it);
  template<typename T> bool getUnorderedMapIterator(TString name, std::unordered_map<TString, T>& theMap, typename std::unordered_map<TString, T>::iterator& it);

  template<typename T, typename U> void replaceString(T& strinput, U strTakeOut, U strPutIn);
  template<> void replaceString<TString, const TString>(TString& strinput, const TString strTakeOut, const TString strPutIn);
  template<> void replaceString<TString, const char*>(TString& strinput, const char* strTakeOut, const char* strPutIn);
  template<> void replaceString<std::string, const std::string>(std::string& strinput, const std::string strTakeOut, const std::string strPutIn);
  template<> void replaceString<std::string, const char*>(std::string& strinput, const char* strTakeOut, const char* strPutIn);

}

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
template bool HelperFunctions::getUnorderedMapIterator<float>(TString name, const std::unordered_map<TString, float>& theMap, std::unordered_map<TString, float>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<float>(TString name, std::unordered_map<TString, float>& theMap, std::unordered_map<TString, float>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<double>(TString name, const std::unordered_map<TString, double>& theMap, std::unordered_map<TString, double>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<double>(TString name, std::unordered_map<TString, double>& theMap, std::unordered_map<TString, double>::iterator& it);

template bool HelperFunctions::getUnorderedMapIterator<std::vector<bool>>(TString name, const std::unordered_map<TString, std::vector<bool>>& theMap, std::unordered_map<TString, std::vector<bool>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<bool>>(TString name, std::unordered_map<TString, std::vector<bool>>& theMap, std::unordered_map<TString, std::vector<bool>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<short>>(TString name, const std::unordered_map<TString, std::vector<short>>& theMap, std::unordered_map<TString, std::vector<short>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<short>>(TString name, std::unordered_map<TString, std::vector<short>>& theMap, std::unordered_map<TString, std::vector<short>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<unsigned int>>(TString name, const std::unordered_map<TString, std::vector<unsigned int>>& theMap, std::unordered_map<TString, std::vector<unsigned int>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<unsigned int>>(TString name, std::unordered_map<TString, std::vector<unsigned int>>& theMap, std::unordered_map<TString, std::vector<unsigned int>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<int>>(TString name, const std::unordered_map<TString, std::vector<int>>& theMap, std::unordered_map<TString, std::vector<int>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<int>>(TString name, std::unordered_map<TString, std::vector<int>>& theMap, std::unordered_map<TString, std::vector<int>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<float>>(TString name, const std::unordered_map<TString, std::vector<float>>& theMap, std::unordered_map<TString, std::vector<float>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<float>>(TString name, std::unordered_map<TString, std::vector<float>>& theMap, std::unordered_map<TString, std::vector<float>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<double>>(TString name, const std::unordered_map<TString, std::vector<double>>& theMap, std::unordered_map<TString, std::vector<double>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<double>>(TString name, std::unordered_map<TString, std::vector<double>>& theMap, std::unordered_map<TString, std::vector<double>>::iterator& it);

template bool HelperFunctions::getUnorderedMapIterator<std::vector<bool>*>(TString name, const std::unordered_map<TString, std::vector<bool>*>& theMap, std::unordered_map<TString, std::vector<bool>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<bool>*>(TString name, std::unordered_map<TString, std::vector<bool>*>& theMap, std::unordered_map<TString, std::vector<bool>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<short>*>(TString name, const std::unordered_map<TString, std::vector<short>*>& theMap, std::unordered_map<TString, std::vector<short>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<short>*>(TString name, std::unordered_map<TString, std::vector<short>*>& theMap, std::unordered_map<TString, std::vector<short>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<unsigned int>*>(TString name, const std::unordered_map<TString, std::vector<unsigned int>*>& theMap, std::unordered_map<TString, std::vector<unsigned int>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<unsigned int>*>(TString name, std::unordered_map<TString, std::vector<unsigned int>*>& theMap, std::unordered_map<TString, std::vector<unsigned int>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<int>*>(TString name, const std::unordered_map<TString, std::vector<int>*>& theMap, std::unordered_map<TString, std::vector<int>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<int>*>(TString name, std::unordered_map<TString, std::vector<int>*>& theMap, std::unordered_map<TString, std::vector<int>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<float>*>(TString name, const std::unordered_map<TString, std::vector<float>*>& theMap, std::unordered_map<TString, std::vector<float>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<float>*>(TString name, std::unordered_map<TString, std::vector<float>*>& theMap, std::unordered_map<TString, std::vector<float>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<double>*>(TString name, const std::unordered_map<TString, std::vector<double>*>& theMap, std::unordered_map<TString, std::vector<double>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::vector<double>*>(TString name, std::unordered_map<TString, std::vector<double>*>& theMap, std::unordered_map<TString, std::vector<double>*>::iterator& it);

template bool HelperFunctions::getUnorderedMapIterator<std::pair<bool, bool>>(TString name, const std::unordered_map<TString, std::pair<bool, bool>>& theMap, std::unordered_map<TString, std::pair<bool, bool>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<bool, bool>>(TString name, std::unordered_map<TString, std::pair<bool, bool>>& theMap, std::unordered_map<TString, std::pair<bool, bool>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<short, short>>(TString name, const std::unordered_map<TString, std::pair<short, short>>& theMap, std::unordered_map<TString, std::pair<short, short>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<short, short>>(TString name, std::unordered_map<TString, std::pair<short, short>>& theMap, std::unordered_map<TString, std::pair<short, short>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<unsigned int, unsigned int>>(TString name, const std::unordered_map<TString, std::pair<unsigned int, unsigned int>>& theMap, std::unordered_map<TString, std::pair<unsigned int, unsigned int>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<unsigned int, unsigned int>>(TString name, std::unordered_map<TString, std::pair<unsigned int, unsigned int>>& theMap, std::unordered_map<TString, std::pair<unsigned int, unsigned int>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<int, int>>(TString name, const std::unordered_map<TString, std::pair<int, int>>& theMap, std::unordered_map<TString, std::pair<int, int>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<int, int>>(TString name, std::unordered_map<TString, std::pair<int, int>>& theMap, std::unordered_map<TString, std::pair<int, int>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<float, float>>(TString name, const std::unordered_map<TString, std::pair<float, float>>& theMap, std::unordered_map<TString, std::pair<float, float>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<float, float>>(TString name, std::unordered_map<TString, std::pair<float, float>>& theMap, std::unordered_map<TString, std::pair<float, float>>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<double, double>>(TString name, const std::unordered_map<TString, std::pair<double, double>>& theMap, std::unordered_map<TString, std::pair<double, double>>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<double, double>>(TString name, std::unordered_map<TString, std::pair<double, double>>& theMap, std::unordered_map<TString, std::pair<double, double>>::iterator& it);

template bool HelperFunctions::getUnorderedMapIterator<std::pair<bool, bool>*>(TString name, const std::unordered_map<TString, std::pair<bool, bool>*>& theMap, std::unordered_map<TString, std::pair<bool, bool>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<bool, bool>*>(TString name, std::unordered_map<TString, std::pair<bool, bool>*>& theMap, std::unordered_map<TString, std::pair<bool, bool>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<short, short>*>(TString name, const std::unordered_map<TString, std::pair<short, short>*>& theMap, std::unordered_map<TString, std::pair<short, short>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<short, short>*>(TString name, std::unordered_map<TString, std::pair<short, short>*>& theMap, std::unordered_map<TString, std::pair<short, short>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<unsigned int, unsigned int>*>(TString name, const std::unordered_map<TString, std::pair<unsigned int, unsigned int>*>& theMap, std::unordered_map<TString, std::pair<unsigned int, unsigned int>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<unsigned int, unsigned int>*>(TString name, std::unordered_map<TString, std::pair<unsigned int, unsigned int>*>& theMap, std::unordered_map<TString, std::pair<unsigned int, unsigned int>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<int, int>*>(TString name, const std::unordered_map<TString, std::pair<int, int>*>& theMap, std::unordered_map<TString, std::pair<int, int>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<int, int>*>(TString name, std::unordered_map<TString, std::pair<int, int>*>& theMap, std::unordered_map<TString, std::pair<int, int>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<float, float>*>(TString name, const std::unordered_map<TString, std::pair<float, float>*>& theMap, std::unordered_map<TString, std::pair<float, float>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<float, float>*>(TString name, std::unordered_map<TString, std::pair<float, float>*>& theMap, std::unordered_map<TString, std::pair<float, float>*>::iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<double, double>*>(TString name, const std::unordered_map<TString, std::pair<double, double>*>& theMap, std::unordered_map<TString, std::pair<double, double>*>::const_iterator& it);
template bool HelperFunctions::getUnorderedMapIterator<std::pair<double, double>*>(TString name, std::unordered_map<TString, std::pair<double, double>*>& theMap, std::unordered_map<TString, std::pair<double, double>*>::iterator& it);

#endif
