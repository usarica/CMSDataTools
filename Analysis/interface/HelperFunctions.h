#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include <cmath>
#include <cassert>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"
#include "TF1.h"
#include "TSpline.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "StdExtensions.h"

namespace HelperFunctions{

  template<typename T> void appendVector(std::vector<T>& a, std::vector<T>& b);

  template<typename T> void addByLowest(std::vector<T>& valArray, T val, bool unique);

  template<typename T, typename U> void addByLowest(std::vector<std::pair<T, U>>& valArray, T val, U index);

  template<typename T, typename U> void addByLowest(std::vector<std::pair<T, U>>& valArray, std::vector<std::pair<T, U>>& inArray, bool consecutive=false, bool inputordered=false);

  template<typename T> bool checkListVariable(const std::vector<T>& list, const T& var);


  void splitOption(const std::string rawoption, std::string& wish, std::string& value, char delimiter);
  void splitOptionRecursive(const std::string rawoption, std::vector<std::string>& splitoptions, char delimiter);

  template<typename T, typename U> void cleanUnorderedMap(std::unordered_map<T, U> um);

  // Non-zero and NaN/Inf checkers
  template<typename T> bool checkNonZero(std::vector<T> const& vars);
  template<typename T> bool checkNanInf(std::vector<T> const& vars);

}

template<typename T> void HelperFunctions::appendVector(std::vector<T>& a, std::vector<T>& b){ a.insert(a.end(), b.begin(), b.end()); }

template<typename T> void HelperFunctions::addByLowest(std::vector<T>& valArray, T val, bool unique){
  bool inserted = false;
  for (typename std::vector<T>::iterator it = valArray.begin(); it<valArray.end(); it++){
    if (*it>val || (!unique && *it==val)){
      inserted=true;
      valArray.insert(it, val);
      break;
    }
  }
  if (!inserted) valArray.push_back(val);
}

template<typename T, typename U> void HelperFunctions::addByLowest(std::vector<std::pair<T, U>>& valArray, T val, U index){
  bool inserted = false;
  for (typename std::vector<std::pair<T, U>>::iterator it = valArray.begin(); it<valArray.end(); it++){
    if ((*it).first>=val){
      inserted=true;
      if ((*it).second!=index) valArray.insert(it, std::pair<T, U>(val, index));
      break;
    }
  }
  if (!inserted) valArray.push_back(std::pair<T, U>(val, index));
}

template<typename T, typename U> void HelperFunctions::addByLowest(std::vector<std::pair<T, U>>& valArray, std::vector<std::pair<T, U>>& inArray, bool consecutive, bool inputordered){
  if (consecutive){
    bool inserted = false;
    typename std::vector<std::pair<T, U>>::iterator inbegin = inArray.begin();
    typename std::vector<std::pair<T, U>>::iterator inend = inArray.end();
    for (typename std::vector<std::pair<T, U>>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if ((*it).first>=(*inbegin).first){
        inserted=true;
        if ((*it).second!=(*inbegin).second) valArray.insert(it, inbegin, inend);
        break;
      }
    }
    if (!inserted) appendVector<std::pair<T, U>>(valArray, inArray);
  }
  else if (!inputordered){
    for (typename std::vector<std::pair<T, U>>::iterator init = inArray.begin(); init<inArray.end(); init++){
      bool inserted = false;
      for (typename std::vector<std::pair<T, U>>::iterator it = valArray.begin(); it<valArray.end(); it++){
        if ((*it).first>=(*init).first){
          inserted=true;
          if ((*it).second!=(*init).second) valArray.insert(it, *init);
          break;
        }
      }
      if (!inserted) valArray.push_back(*init);
    }
  }
  else if (inArray.size()>0){
    typename std::vector<std::pair<T, U>>::iterator infirst = inArray.begin();
    typename std::vector<std::pair<T, U>>::iterator inlast = inArray.end()-1;
    typename std::vector<std::pair<T, U>>::iterator valfirst = valArray.begin();
    typename std::vector<std::pair<T, U>>::iterator vallast = valArray.end()-1;
    while ((*valfirst).first<(*infirst).first) valfirst++;
    while ((*vallast).first>=(*inlast).first) vallast--;
    vallast++;
    inlast++;

    for (typename std::vector<std::pair<T, U>>::iterator init = infirst; init<inlast; init++){
      bool inserted = false;
      for (typename std::vector<std::pair<T, U>>::iterator it = valfirst; it<vallast; it++){
        if ((*it).first>=(*init).first){
          inserted=true;
          if ((*it).second!=(*init).second) valArray.insert(it, *init);
          break;
        }
      }
      if (!inserted) valArray.insert(vallast, *init);
    }
  }
}

template<typename T> bool HelperFunctions::checkListVariable(const std::vector<T>& list, const T& var){
  for (unsigned int v=0; v<list.size(); v++){
    if (list.at(v)==var) return true; // Look for exact match
  }
  return false;
}

template<typename T, typename U> void HelperFunctions::cleanUnorderedMap(std::unordered_map<T, U> um){ for (auto& it:um){ delete it.second; it.second=0; } }

// Non-zero and NaN/Inf checkers
template<typename T> bool HelperFunctions::checkNonZero(std::vector<T> const& vars){
  for (T const& v:vars){
    if (v<0.){
      std::cerr << "checkNonZero found value < 0" << std::endl;
      return false;
    }
  }
  return true;
}
template<typename T> bool HelperFunctions::checkNanInf(std::vector<T> const& vars){
  for (T const& v:vars){
    if (std::isnan(v) || std::isinf(v)) return false;
  }
  return true;
}

/****************************************************************/
// Explicit instantiations
template void HelperFunctions::addByLowest<double, int>(std::vector<std::pair<double, int>>& valArray, double val, int index);
template void HelperFunctions::addByLowest<double, double>(std::vector<std::pair<double, double>>& valArray, double val, double index);
template void HelperFunctions::addByLowest<double, int>(std::vector<std::pair<double, int>>& valArray, std::vector<std::pair<double, int>>& inArray, bool consecutive, bool inputordered);
template void HelperFunctions::addByLowest<double, double>(std::vector<std::pair<double, double>>& valArray, std::vector<std::pair<double, double>>& inArray, bool consecutive, bool inputordered);

template bool HelperFunctions::checkListVariable<std::string>(const std::vector<std::string>& list, const std::string& var);
template bool HelperFunctions::checkListVariable<double>(const std::vector<double>& list, const double& var);

template bool HelperFunctions::checkNonZero<short>(std::vector<short> const& vars);
template bool HelperFunctions::checkNanInf<short>(std::vector<short> const& vars);
template bool HelperFunctions::checkNonZero<unsigned int>(std::vector<unsigned int> const& vars);
template bool HelperFunctions::checkNanInf<unsigned int>(std::vector<unsigned int> const& vars);
template bool HelperFunctions::checkNonZero<int>(std::vector<int> const& vars);
template bool HelperFunctions::checkNanInf<int>(std::vector<int> const& vars);
template bool HelperFunctions::checkNonZero<float>(std::vector<float> const& vars);
template bool HelperFunctions::checkNanInf<float>(std::vector<float> const& vars);
template bool HelperFunctions::checkNonZero<double>(std::vector<double> const& vars);
template bool HelperFunctions::checkNanInf<double>(std::vector<double> const& vars);

template void HelperFunctions::cleanUnorderedMap<TString, short*>(std::unordered_map<TString, short*> um);
template void HelperFunctions::cleanUnorderedMap<TString, unsigned int*>(std::unordered_map<TString, unsigned int*> um);
template void HelperFunctions::cleanUnorderedMap<TString, int*>(std::unordered_map<TString, int*> um);
template void HelperFunctions::cleanUnorderedMap<TString, float*>(std::unordered_map<TString, float*> um);
template void HelperFunctions::cleanUnorderedMap<TString, double*>(std::unordered_map<TString, double*> um);
template void HelperFunctions::cleanUnorderedMap<TString, std::pair<short, short>*>(std::unordered_map<TString, std::pair<short, short>*> um);
template void HelperFunctions::cleanUnorderedMap<TString, std::pair<unsigned int, unsigned int>*>(std::unordered_map<TString, std::pair<unsigned int, unsigned int>*> um);
template void HelperFunctions::cleanUnorderedMap<TString, std::pair<int, int>*>(std::unordered_map<TString, std::pair<int, int>*> um);
template void HelperFunctions::cleanUnorderedMap<TString, std::pair<float, float>*>(std::unordered_map<TString, std::pair<float, float>*> um);
template void HelperFunctions::cleanUnorderedMap<TString, std::pair<double, double>*>(std::unordered_map<TString, std::pair<double, double>*> um);
template void HelperFunctions::cleanUnorderedMap<TString, std::vector<short>*>(std::unordered_map<TString, std::vector<short>*> um);
template void HelperFunctions::cleanUnorderedMap<TString, std::vector<unsigned int>*>(std::unordered_map<TString, std::vector<unsigned int>*> um);
template void HelperFunctions::cleanUnorderedMap<TString, std::vector<int>*>(std::unordered_map<TString, std::vector<int>*> um);
template void HelperFunctions::cleanUnorderedMap<TString, std::vector<float>*>(std::unordered_map<TString, std::vector<float>*> um);
template void HelperFunctions::cleanUnorderedMap<TString, std::vector<double>*>(std::unordered_map<TString, std::vector<double>*> um);

#endif
