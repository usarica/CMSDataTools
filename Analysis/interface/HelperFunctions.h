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
#include "SimpleEntry.h"

namespace HelperFunctions{

  template<typename T> void appendVector(std::vector<T>& a, std::vector<T>& b);

  template<typename T> void addByLowest(std::vector<T>& valArray, T val, bool unique);

  template<typename T, typename U> void addByLowest(std::vector<std::pair<T, U>>& valArray, T val, U index);

  template<typename T, typename U> void addByLowest(std::vector<std::pair<T, U>>& valArray, std::vector<std::pair<T, U>>& inArray, bool consecutive=false, bool inputordered=false);

  template<typename T> bool checkListVariable(const std::vector<T>& list, const T& var);

  template<typename T, typename U> void cleanUnorderedMap(std::unordered_map<T, U> um);

  // Non-zero and NaN/Inf checkers
  template<typename T> bool checkNonZero(std::vector<T> const& vars);
  template<typename T> bool checkNanInf(std::vector<T> const& vars);

  // TGraph functions
  template<typename T> TGraph* makeGraphFromPair(std::vector<std::pair<T, T>> points, TString name);

  template<typename T> void addPointsBetween(T*& tgOriginal, double xmin, double xmax, unsigned int nadd);
  template<> void addPointsBetween(TGraph*& tgOriginal, double xmin, double xmax, unsigned int nadd);
  template<> void addPointsBetween(TGraphErrors*& tgOriginal, double xmin, double xmax, unsigned int nadd);

  // Spline functions
  template<int N> TF1* getFcn_a0plusa1overXN(TSpline3* sp, double xmin, double xmax, bool useLowBound);
  template<int N> TF1* getFcn_a0plusa1timesXN(TSpline3* sp, double xmin, double xmax, bool useLowBound);

  // Non-template functions
  void splitOption(const std::string rawoption, std::string& wish, std::string& value, char delimiter);

  void splitOptionRecursive(const std::string rawoption, std::vector<std::string>& splitoptions, char delimiter);

  TSpline3* convertGraphToSpline3(TGraph* tg, bool faithfulFirst=false, bool faithfulSecond=false, double* dfirst=nullptr, double* dlast=nullptr);

  TGraphErrors* makeGraphFromTH1(TH1* hx, TH1* hy, TString name);

  TGraph* multiplyTGraphs(TGraph* tgfirst, TGraph* tgsecond);

  TGraph* divideTGraphs(TGraph* tgnum, TGraph* tgdenom, double powernum=1, double powerdenom=1);

  TGraphErrors* addPoint(TGraphErrors* tgSlice, double x);

  void addPoint(TGraph*& tg, double x, double y);

  void addPoint(TGraphErrors*& tg, double x, double y, double ex, double ey);

  TF1* getFcn_a0plusa1expX(TSpline3* sp, double xmin, double xmax, bool useLowBound);
  TF1* getFcn_a0timesexpa1X(TSpline3* sp, double xmin, double xmax, bool useLowBound);
}

template<typename T> void HelperFunctions::addPointsBetween(T*& tgOriginal, double xmin, double xmax, unsigned int nadd){} // Dummy definition for generic types

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

// TGraph functions
template<typename T> TGraph* HelperFunctions::makeGraphFromPair(std::vector<std::pair<T, T>> points, TString name){
  if (points.size()==0) return 0;
  unsigned int nbins = points.size();
  double* xy[2];
  for (unsigned int ix=0; ix<2; ix++) xy[ix] = new double[nbins];
  for (unsigned int bin=0; bin<nbins; bin++){
    xy[0][bin] = points[bin].first;
    xy[1][bin] = points[bin].second;
  }
  TGraph* tg = new TGraph(nbins, xy[0], xy[1]);
  tg->SetName(name);
  for (unsigned int ix=0; ix<2; ix++) delete[] xy[ix];
  return tg;
}

// Spline functions
template<int N> TF1* HelperFunctions::getFcn_a0plusa1overXN(TSpline3* sp, double xmin, double xmax, bool useLowBound){
  double x, y, s;
  if (useLowBound) x = sp->GetXmin();
  else x = sp->GetXmax();
  y = sp->Eval(x);
  s = sp->Derivative(x);

  // y=a0+a1/x^N
  // s=-N*a1/x^(N+1)
  double a0, a1;
  a0 = y+s*x/(double(N));
  a1 = -s*pow(x, N+1)/(double(N));

  TString fcnName;
  if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
  else fcnName = Form("highFcn_%s", sp->GetName());
  TF1* fcn = new TF1(fcnName, Form("[0]+[1]/pow(x, %i)", N), xmin, xmax);
  fcn->SetParameter(0, a0);
  fcn->SetParameter(1, a1);

  return fcn;
}
template<int N> TF1* HelperFunctions::getFcn_a0plusa1timesXN(TSpline3* sp, double xmin, double xmax, bool useLowBound){ return getFcn_a0plusa1overXN<-N>(sp, xmin, xmax, useLowBound); }

/****************************************************************/
// Explicit instantiations
template void HelperFunctions::addByLowest<SimpleEntry>(std::vector<SimpleEntry>& valArray, SimpleEntry val, bool unique);
template void HelperFunctions::addByLowest<int>(std::vector<int>& valArray, int val, bool unique);
template void HelperFunctions::addByLowest<float>(std::vector<float>& valArray, float val, bool unique);
template void HelperFunctions::addByLowest<double>(std::vector<double>& valArray, double val, bool unique);
template void HelperFunctions::addByLowest<SimpleEntry, int>(std::vector<std::pair<SimpleEntry, int>>& valArray, SimpleEntry val, int index);
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

template TGraph* HelperFunctions::makeGraphFromPair<float>(std::vector<std::pair<float, float>> points, TString name);
template TGraph* HelperFunctions::makeGraphFromPair<double>(std::vector<std::pair<double, double>> points, TString name);

template TF1* HelperFunctions::getFcn_a0plusa1overXN<-1>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1overXN<-2>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1overXN<-3>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1overXN<-4>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1overXN<-5>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1overXN<-6>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1overXN<-7>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1overXN<-8>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1overXN<-9>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1overXN<-10>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1overXN<1>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1overXN<2>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1overXN<3>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1overXN<4>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1overXN<5>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1overXN<6>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1overXN<7>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1overXN<8>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1overXN<9>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1overXN<10>(TSpline3* sp, double xmin, double xmax, bool useLowBound);

template TF1* HelperFunctions::getFcn_a0plusa1timesXN<-1>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1timesXN<-2>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1timesXN<-3>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1timesXN<-4>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1timesXN<-5>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1timesXN<-6>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1timesXN<-7>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1timesXN<-8>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1timesXN<-9>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1timesXN<-10>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1timesXN<1>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1timesXN<2>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1timesXN<3>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1timesXN<4>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1timesXN<5>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1timesXN<6>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1timesXN<7>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1timesXN<8>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1timesXN<9>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
template TF1* HelperFunctions::getFcn_a0plusa1timesXN<10>(TSpline3* sp, double xmin, double xmax, bool useLowBound);


#endif
