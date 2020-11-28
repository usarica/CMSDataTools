#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include <iterator>
#include <ctime>
#include <cmath>
#include <limits>
#include <cassert>
#include <cstdlib>
#include "TROOT.h"
#include "TSystem.h"
#include "TObject.h"
#include "TKey.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TF1.h"
#include "TSpline.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "HelperFunctionsCore.h"
#include "ExtendedBinning.h"
#include "SimpleEntry.h"
#include "MELAStreamHelpers.hh"
#include "Mela.h"


namespace HelperFunctions{

  template<typename T> void appendVector(std::vector<T>& a, std::vector<T> const& b);

  template<typename T> void addByLowest(std::vector<T>& valArray, T val, bool unique);
  template<typename T, typename U> void addByLowest(std::vector<std::pair<T, U>>& valArray, T val, U index);
  template<typename T, typename U> void addByLowest(std::vector<std::pair<T, U>>& valArray, std::vector<std::pair<T, U>>& inArray, bool consecutive=false, bool inputordered=false);

  template<typename T> void addByHighest(std::vector<T>& valArray, T val, bool unique);
  template<typename T, typename U> void addByHighest(std::vector<std::pair<T, U>>& valArray, T val, U index);
  template<typename T, typename U> void addByHighest(std::vector<std::pair<T, U>>& valArray, std::vector<std::pair<T, U>>& inArray, bool consecutive=false, bool inputordered=false);

  template<typename T> bool checkListVariable(std::vector<T> const& list, T const& var);
  template<typename T> bool hasCommonElements(std::vector<T> const& list1, std::vector<T> const& list2);

  template<typename T, typename U> void cleanUnorderedMap(std::unordered_map<T, U>& um);

  // Non-zero and NaN/Inf checkers
  template<typename T> bool checkVarNonNegative(T const& val);
  template<typename T> bool checkNonNegative(std::vector<T> const& vars, int ibegin=-1, int iend=-1);
  template<typename T> bool checkVarNonZero(T const& val);
  template<typename T> bool checkNonZero(std::vector<T> const& vars, int ibegin=-1, int iend=-1);
  template<typename T> bool checkVarPositiveDef(T const& val);
  template<typename T> bool checkPositiveDef(std::vector<T> const& vars, int ibegin=-1, int iend=-1);
  template<typename T> bool checkVarNanInf(T const& val);
  template<typename T> bool checkNanInf(std::vector<T> const& vars);

  template<> bool checkVarNonNegative<TH1F>(TH1F const& val);
  template<> bool checkVarNonNegative<TH2F>(TH2F const& val);
  template<> bool checkVarNonNegative<TH3F>(TH3F const& val);

  // Bit set and test
  template<typename T> void set_bit(T& mask, unsigned int iBit, bool val=true);
  template<typename T> bool test_bit(T const& mask, unsigned int iBit);

  // Vector helpers
  template<typename T> void deltaR(T const& e1, T const& p1, T const& e2, T const& p2, T& res);
  template<typename T> void deltaEta(T const& v1, T const& v2, T& res);
  template<typename T> void deltaPhi(T const& v1, T const& v2, T& res);

  // TGraph functions
  template<typename T> TGraph* makeGraphFromPair(std::vector<std::pair<T, T>> points, TString name);
  template<typename T> TGraphErrors* makeGraphSymErrFromPair(std::vector<std::pair<T, T>> points, std::vector<std::pair<T, T>> errors, TString name);
  template<typename T> TGraphAsymmErrors* makeGraphAsymErrFromPair(std::vector<std::pair<T, T>> points, std::vector<std::pair<T, T>> errorDns, std::vector<std::pair<T, T>> errorUps, TString name);

  template<typename T> void addPointsBetween(T*& tgOriginal, double xmin, double xmax, unsigned int nadd);
  template<> void addPointsBetween(TGraph*& tgOriginal, double xmin, double xmax, unsigned int nadd);
  template<> void addPointsBetween(TGraphErrors*& tgOriginal, double xmin, double xmax, unsigned int nadd);

  template <typename T> double evaluateTObject(T* obj, float val);
  template<> double evaluateTObject<TH1F>(TH1F* obj, float val);
  template<> double evaluateTObject<TGraph>(TGraph* obj, float val);
  template<> double evaluateTObject<TSpline3>(TSpline3* obj, float val);

  template<typename T> bool checkHistogramIntegrity(T const* histo);
  template<> bool checkHistogramIntegrity<TH1>(TH1 const* histo);
  template<> bool checkHistogramIntegrity<TH2>(TH2 const* histo);
  template<> bool checkHistogramIntegrity<TH3>(TH3 const* histo);
  template<> bool checkHistogramIntegrity<TH1F>(TH1F const* histo);
  template<> bool checkHistogramIntegrity<TH2F>(TH2F const* histo);
  template<> bool checkHistogramIntegrity<TH3F>(TH3F const* histo);
  template<> bool checkHistogramIntegrity<TH1D>(TH1D const* histo);
  template<> bool checkHistogramIntegrity<TH2D>(TH2D const* histo);
  template<> bool checkHistogramIntegrity<TH3D>(TH3D const* histo);

  template <typename T> void regularizeHistogram(T*& histo, int nIter_, double threshold_, double acceleration_);
  template<> void regularizeHistogram<TH1F>(TH1F*& histo, int nIter_, double threshold_, double acceleration_);
  template<> void regularizeHistogram<TH2F>(TH2F*& histo, int nIter_, double threshold_, double acceleration_);
  //template<> void regularizeHistogram<TH3F>(TH3F*& histo, int nIter_, double threshold_, double acceleration_);

  template <typename T> void conditionalizeHistogram(T* histo, unsigned int iaxis, std::vector<std::pair<T*, float>> const* conditionalsReference=nullptr, bool useWidth=true, bool useEffErr=false);
  template<> void conditionalizeHistogram<TH2F>(TH2F* histo, unsigned int iaxis, std::vector<std::pair<TH2F*, float>> const* conditionalsReference, bool useWidth, bool useEffErr);
  template<> void conditionalizeHistogram<TH3F>(TH3F* histo, unsigned int iaxis, std::vector<std::pair<TH3F*, float>> const* conditionalsReference, bool useWidth, bool useEffErr);

  template <typename T> void wipeOverUnderFlows(T* hwipe, bool rescale=false, bool addToLastBin=false);
  template<> void wipeOverUnderFlows<TH1F>(TH1F* hwipe, bool rescale, bool addToLastBin);
  template<> void wipeOverUnderFlows<TH2F>(TH2F* hwipe, bool rescale, bool addToLastBin);
  template<> void wipeOverUnderFlows<TH3F>(TH3F* hwipe, bool rescale, bool addToLastBin);

  template <typename T> void divideBinWidth(T* histo);
  template<> void divideBinWidth<TH1F>(TH1F* histo);
  template<> void divideBinWidth<TH2F>(TH2F* histo);
  template<> void divideBinWidth<TH3F>(TH3F* histo);

  template <typename T> void multiplyBinWidth(T* histo);
  template<> void multiplyBinWidth<TH1F>(TH1F* histo);
  template<> void multiplyBinWidth<TH2F>(TH2F* histo);
  template<> void multiplyBinWidth<TH3F>(TH3F* histo);

  template <typename T> double getHistogramIntegralAndError(T const* histo, int ix, int jx, bool useWidth, double* error=nullptr);
  template <typename T> double getHistogramIntegralAndError(T const* histo, int ix, int jx, int iy, int jy, bool useWidth, double* error=nullptr);
  template <typename T> double getHistogramIntegralAndError(T const* histo, int ix, int jx, int iy, int jy, int iz, int jz, bool useWidth, double* error=nullptr);

  template <typename T> float computeIntegral(T* histo, bool useWidth);
  template<> float computeIntegral<TH1F>(TH1F* histo, bool useWidth);
  template<> float computeIntegral<TH2F>(TH2F* histo, bool useWidth);
  template<> float computeIntegral<TH3F>(TH3F* histo, bool useWidth);

  template <typename T> double computeChiSq(T const* h1, T const* h2);
  template<> double computeChiSq<TH1F>(TH1F const* h1, TH1F const* h2);
  template<> double computeChiSq<TH2F>(TH2F const* h1, TH2F const* h2);
  template<> double computeChiSq<TH3F>(TH3F const* h1, TH3F const* h2);

  template <typename T> void divideHistograms(T const* hnum, T const* hden, T*& hAssign, bool useEffErr);
  template<> void divideHistograms<TH1F>(TH1F const* hnum, TH1F const* hden, TH1F*& hAssign, bool useEffErr);
  template<> void divideHistograms<TH2F>(TH2F const* hnum, TH2F const* hden, TH2F*& hAssign, bool useEffErr);
  template<> void divideHistograms<TH3F>(TH3F const* hnum, TH3F const* hden, TH3F*& hAssign, bool useEffErr);

  template <typename T> void multiplyHistograms(T const* h1, T const* h2, T*& hAssign, bool useEffErr);
  template<> void multiplyHistograms<TH1F>(TH1F const* h1, TH1F const* h2, TH1F*& hAssign, bool useEffErr);
  template<> void multiplyHistograms<TH2F>(TH2F const* h1, TH2F const* h2, TH2F*& hAssign, bool useEffErr);
  template<> void multiplyHistograms<TH3F>(TH3F const* h1, TH3F const* h2, TH3F*& hAssign, bool useEffErr);

  template <typename T> void multiplyHistograms(T const* h1, TH1F const* h2, unsigned int matchDimension, T*& hAssign, bool useEffErr);
  template<> void multiplyHistograms<TH2F>(TH2F const* h1, TH1F const* h2, unsigned int matchDimension, TH2F*& hAssign, bool useEffErr);
  template<> void multiplyHistograms<TH3F>(TH3F const* h1, TH1F const* h2, unsigned int matchDimension, TH3F*& hAssign, bool useEffErr);

  template <typename T> void symmetrizeHistogram(T* histo, unsigned int const axis=0);
  template <> void symmetrizeHistogram<TH1F>(TH1F* histo, unsigned int const axis);
  template <> void symmetrizeHistogram<TH2F>(TH2F* histo, unsigned int const axis);
  template <> void symmetrizeHistogram<TH3F>(TH3F* histo, unsigned int const axis);

  template <typename T> void antisymmetrizeHistogram(T* histo, unsigned int const axis=0);
  template <> void antisymmetrizeHistogram<TH1F>(TH1F* histo, unsigned int const axis);
  template <> void antisymmetrizeHistogram<TH2F>(TH2F* histo, unsigned int const axis);
  template <> void antisymmetrizeHistogram<TH3F>(TH3F* histo, unsigned int const axis);

  template <typename T> void getCumulantHistogram(T const* histo, T*& res, TString newname="", std::vector<unsigned int>* condDims=nullptr);
  template <> void getCumulantHistogram<TH1F>(TH1F const* histo, TH1F*& res, TString newname, std::vector<unsigned int>* condDims);
  template <> void getCumulantHistogram<TH2F>(TH2F const* histo, TH2F*& res, TString newname, std::vector<unsigned int>* condDims);
  template <> void getCumulantHistogram<TH3F>(TH3F const* histo, TH3F*& res, TString newname, std::vector<unsigned int>* condDims);

  template <typename T> void translateCumulantToHistogram(T const* histo, T*& res, TString newname="", std::vector<unsigned int>* condDims=nullptr);
  template <> void translateCumulantToHistogram<TH1F>(TH1F const* histo, TH1F*& res, TString newname, std::vector<unsigned int>* condDims);
  template <> void translateCumulantToHistogram<TH2F>(TH2F const* histo, TH2F*& res, TString newname, std::vector<unsigned int>* condDims);
  template <> void translateCumulantToHistogram<TH3F>(TH3F const* histo, TH3F*& res, TString newname, std::vector<unsigned int>* condDims);

  template <typename T> void combineHistogramListByWeightedAverage(std::vector<T const*> const& hList, T*& hAssign, bool useNeff=false);
  template <> void combineHistogramListByWeightedAverage<TProfile>(std::vector<TProfile const*> const& hList, TProfile*& hAssign, bool useNeff);
  template <> void combineHistogramListByWeightedAverage<TH1F>(std::vector<TH1F const*> const& hList, TH1F*& hAssign, bool useNeff);
  template <> void combineHistogramListByWeightedAverage<TH2F>(std::vector<TH2F const*> const& hList, TH2F*& hAssign, bool useNeff);
  template <> void combineHistogramListByWeightedAverage<TH3F>(std::vector<TH3F const*> const& hList, TH3F*& hAssign, bool useNeff);

  template <typename T> void combineHistogramsByWeightedAverage(T const* h1, T const* h2, T*& hAssign, bool useNeff=false);
  template void combineHistogramsByWeightedAverage<TProfile>(TProfile const* h1, TProfile const* h2, TProfile*& hAssign, bool useNeff);
  template void combineHistogramsByWeightedAverage<TH1F>(TH1F const* h1, TH1F const* h2, TH1F*& hAssign, bool useNeff);
  template void combineHistogramsByWeightedAverage<TH2F>(TH2F const* h1, TH2F const* h2, TH2F*& hAssign, bool useNeff);
  template void combineHistogramsByWeightedAverage<TH3F>(TH3F const* h1, TH3F const* h2, TH3F*& hAssign, bool useNeff);

  template <typename T> void findBinContentRange(T const* h, float& bcmin, float& bcmax, bool includeBinErrors, bool includeOverUnderflows, bool onlyPositiveBins);
  template <> void findBinContentRange<TProfile>(TProfile const* h, float& bcmin, float& bcmax, bool includeBinErrors, bool includeOverUnderflows, bool onlyPositiveBins);
  template <> void findBinContentRange<TH1F>(TH1F const* h, float& bcmin, float& bcmax, bool includeBinErrors, bool includeOverUnderflows, bool onlyPositiveBins);
  template <> void findBinContentRange<TH2F>(TH2F const* h, float& bcmin, float& bcmax, bool includeBinErrors, bool includeOverUnderflows, bool onlyPositiveBins);
  template <> void findBinContentRange<TH3F>(TH3F const* h, float& bcmin, float& bcmax, bool includeBinErrors, bool includeOverUnderflows, bool onlyPositiveBins);

  template <typename T> void findBinContentRange(std::vector<T*> hlist, float& bcmin, float& bcmax, bool includeBinErrors, bool includeOverUnderflows, bool onlyPositiveBins);

  // Spline functions
  template<int N> TF1* getFcn_a0plusa1overXN(TSpline3* sp, double xmin, double xmax, bool useLowBound);
  template<int N> TF1* getFcn_a0plusa1timesXN(TSpline3* sp, double xmin, double xmax, bool useLowBound);

  // Non-template functions
  TString todaysdate();

  void progressbar(unsigned int val, unsigned int tot);

  void splitOption(const std::string& rawoption, std::string& wish, std::string& value, char delimiter);
  void splitOption(const TString& rawoption, TString& wish, TString& value, char delimiter);

  void splitOptionRecursive(const std::string& rawoption, std::vector<std::string>& splitoptions, char delimiter, bool uniqueResults=true);
  void splitOptionRecursive(const TString& rawoption, std::vector<TString>& splitoptions, char delimiter, bool uniqueResults=true);

  void getExtendedBinning(TAxis const* theAxis, ExtendedBinning& res);
  void getExtendedBinning(TH1 const* histo, unsigned int iaxis, ExtendedBinning& res);

  TSpline3* convertGraphToSpline3(TGraph const* tg, bool faithfulFirst=false, bool faithfulSecond=false, double* dfirst=nullptr, double* dlast=nullptr);

  void convertTGraphErrorsToTH1F(TGraphErrors const* tg, TH1F*& histo);

  void convertTGraphAsymmErrorsToTH1F(TGraphAsymmErrors const* tg, TH1F*& histo);

  void convertTH1FToTGraphAsymmErrors(TH1F const* histo, TGraphAsymmErrors*& tg, bool errorsOnZero=false, bool useAsymError=false);

  TGraph* createROCFromDistributions(TH1 const* hA, TH1 const* hB, TString name);

  TGraphErrors* makeGraphFromTH1(TH1 const* hx, TH1 const* hy, TString name);

  TGraphErrors* makeGraphFromCumulantHistogram(TH1 const* histo, TString name);

  void multiplyTGraph(TGraph* tg, const double scale);

  TGraph* addTGraphs(TGraph* tgfirst, TGraph* tgsecond);

  TGraph* multiplyTGraphs(TGraph* tgfirst, TGraph* tgsecond);

  TGraph* divideTGraphs(TGraph* tgnum, TGraph* tgdenom, double powernum=1, double powerdenom=1);

  TGraphErrors* addPoint(TGraphErrors* tgSlice, double x);

  void addPoint(TGraph*& tg, double x, double y);

  void addPoint(TGraphErrors*& tg, double x, double y, double ex, double ey);

  TF1* getFcn_a0plusa1expX(TSpline3* sp, double xmin, double xmax, bool useLowBound);
  TF1* getFcn_a0timesexpa1X(TSpline3* sp, double xmin, double xmax, bool useLowBound);
  TF1* getFcn_a0plusa1timesexpXovera2(TSpline3* sp, double xmin, double xmax, bool useLowBound);
  TF1* getFcn_a0plusa1timesatana2timesXminusa3(TSpline3* sp, double xmin, double xmax, bool useLowBound);
  TF1* getFcn_EfficiencyAtan(TSpline3* sp, double xmin, double xmax, bool useLowBound);

  TGraph* genericPatcher(
    TGraph* tg, const TString newname,
    double const xmin, double const xmax,
    TF1* (*lowf)(TSpline3*, double, double, bool),
    TF1* (*highf)(TSpline3*, double, double, bool),
    bool useFaithfulSlopeFirst, bool useFaithfulSlopeSecond,
    std::vector<std::pair<std::pair<double, double>, unsigned int>>* addpoints=nullptr
  );

  void regularizeSlice(
    TGraph* tgSlice,
    std::vector<double>* fixedX=nullptr, double omitbelow=0., double omitabove=0.,
    int nIter_=-1, double threshold_=-1,
    signed char forceUseFaithfulSlopeFirst=-1, signed char forceUseFaithfulSlopeSecond=-1
  );
  void regularizeSlice(
    TGraphErrors* tgSlice,
    std::vector<double>* fixedX=nullptr, double omitbelow=0., double omitabove=0.,
    int nIter_=-1, double threshold_=-1, double acceleration_=-1,
    signed char forceUseFaithfulSlopeFirst=-1, signed char forceUseFaithfulSlopeSecond=-1
  );

  void rebinProfile(TProfile*& prof, const ExtendedBinning& binningX);

  void rebinCumulant(TH1F*& histo, const ExtendedBinning& binningX);
  void rebinCumulant(TH2F*& histo, const ExtendedBinning& binningX, const ExtendedBinning& binningY, std::vector<std::pair<TProfile const*, unsigned int>>* condProfs=nullptr);
  void rebinCumulant(TH3F*& histo, const ExtendedBinning& binningX, const ExtendedBinning& binningY, const ExtendedBinning& binningZ, std::vector<std::pair<TProfile const*, unsigned int>>* condProfs=nullptr);

  void rebinHistogram(TH1F*& histo, const ExtendedBinning& binningX);
  void rebinHistogram(TH2F*& histo, const ExtendedBinning& binningX, const ExtendedBinning& binningY, std::vector<std::pair<TProfile const*, unsigned int>>* condProfs=nullptr);
  void rebinHistogram(TH3F*& histo, const ExtendedBinning& binningX, const ExtendedBinning& binningY, const ExtendedBinning& binningZ, std::vector<std::pair<TProfile const*, unsigned int>>* condProfs=nullptr);

  void rebinHistogram_NoCumulant(TH1F*& histo, const ExtendedBinning& binningX, const TProfile* prof_x);
  void rebinHistogram_NoCumulant(TH2F*& histo, const ExtendedBinning& binningX, const TProfile* prof_x, const ExtendedBinning& binningY, const TProfile* prof_y);
  void rebinHistogram_NoCumulant(TH3F*& histo, const ExtendedBinning& binningX, const TProfile* prof_x, const ExtendedBinning& binningY, const TProfile* prof_y, const ExtendedBinning& binningZ, const TProfile* prof_z);

  template <typename T, typename U> void getGenericHistogramSlice(T*& res, U const* histo, unsigned char XDirection, int iy, int jy, TString newname="");
  template <typename T, typename U> void getGenericHistogramSlice(T*& res, U const* histo, unsigned char XDirection, int iy, int jy, int iz, int jz, TString newname=""); // "y" and "z" are cylical, so if Xdirection==1 (Y), "y"=Z and "z"=X
  template <typename T, typename U> void getGenericHistogramSlice(T*& res, U const* histo, unsigned char XDirection, unsigned char YDirection, int iz, int jz, TString newname="");

  TH1F* getHistogramSlice(TH2F const* histo, unsigned char XDirection, int iy, int jy, TString newname="");
  TH1F* getHistogramSlice(TH3F const* histo, unsigned char XDirection, int iy, int jy, int iz, int jz, TString newname=""); // "y" and "z" are cylical, so if Xdirection==1 (Y), "y"=Z and "z"=X
  TH2F* getHistogramSlice(TH3F const* histo, unsigned char XDirection, unsigned char YDirection, int iz, int jz, TString newname="");

  // Function to calculate error in product/division
  float calculateSimpleProductError(
    float const v1, float const e1, float const p1,
    float const v2, float const e2, float const p2
  );
  // Function to calculate error in efficiency
  float calculateEfficiencyError(
    float const sumW, float const sumWAll,
    float const sumWsq, float const sumWsqAll
    );
  float translateEfficiencyErrorToNumeratorError(
    float const eff, float const sumWAll,
    float const effErr, float const sumWsqAll
  );

  // Function to copy a file
  void CopyFile(TString fname, TTree*(*fcnTree)(TTree*), TDirectory*(*fcnDirectory)(TDirectory*));
  void CopyDirectory(TDirectory* source, TTree*(*fcnTree)(TTree*), TDirectory*(*fcnDirectory)(TDirectory*));

  // Function to distribute a directory/file into multiple chunks. Splits trees.
  void distributeObjects(TDirectory* inputdir, std::vector<TDirectory*> const& outputdirs, TVar::VerbosityLevel verbosity = TVar::ERROR);

  // Function to extract all trees in a file
  void extractTreesFromDirectory(TDirectory* source, std::vector<TTree*>& res, bool doClone=false);

  // Function to extract all histograms from file
  template<typename T> void extractObjectsFromDirectory(TDirectory* source, std::vector<T*>& objlist, TVar::VerbosityLevel verbosity = TVar::ERROR);
  template<typename T> void extractHistogramsFromDirectory(TDirectory* source, std::vector<T*>& histolist, TVar::VerbosityLevel verbosity = TVar::ERROR){ extractObjectsFromDirectory<T>(source, histolist, verbosity); } // Legacy call

}

template<typename T> void HelperFunctions::addPointsBetween(T*& tgOriginal, double xmin, double xmax, unsigned int nadd){} // Dummy definition for generic types

template<typename T> void HelperFunctions::appendVector(std::vector<T>& a, std::vector<T> const& b){ a.insert(a.end(), b.cbegin(), b.cend()); }

template<typename T> void HelperFunctions::addByLowest(std::vector<T>& valArray, T val, bool unique){
  bool inserted = false;
  if (unique){
    for (typename std::vector<T>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if (*it==val){
        inserted=true;
        break;
      }
    }
  }
  if (!inserted){
    for (typename std::vector<T>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if (*it>=val){
        inserted=true;
        valArray.insert(it, val);
        break;
      }
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

template<typename T> void HelperFunctions::addByHighest(std::vector<T>& valArray, T val, bool unique){
  bool inserted = false;
  if (unique){
    for (typename std::vector<T>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if (*it==val){
        inserted=true;
        break;
      }
    }
  }
  if (!inserted){
    for (typename std::vector<T>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if (*it<=val){
        inserted=true;
        valArray.insert(it, val);
        break;
      }
    }
  }
  if (!inserted) valArray.push_back(val);
}
template<typename T, typename U> void HelperFunctions::addByHighest(std::vector<std::pair<T, U>>& valArray, T val, U index){
  bool inserted = false;
  for (typename std::vector<std::pair<T, U>>::iterator it = valArray.begin(); it<valArray.end(); it++){
    if ((*it).first<=val){
      inserted=true;
      if ((*it).second!=index) valArray.insert(it, std::pair<T, U>(val, index));
      break;
    }
  }
  if (!inserted) valArray.push_back(std::pair<T, U>(val, index));
}
template<typename T, typename U> void HelperFunctions::addByHighest(std::vector<std::pair<T, U>>& valArray, std::vector<std::pair<T, U>>& inArray, bool consecutive, bool inputordered){
  if (consecutive){
    bool inserted = false;
    typename std::vector<std::pair<T, U>>::iterator inbegin = inArray.begin();
    typename std::vector<std::pair<T, U>>::iterator inend = inArray.end();
    for (typename std::vector<std::pair<T, U>>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if ((*it).first<=(*inbegin).first){
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
        if ((*it).first<=(*init).first){
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
    while ((*valfirst).first>(*infirst).first) valfirst++;
    while ((*vallast).first<=(*inlast).first) vallast--;
    vallast++;
    inlast++;

    for (typename std::vector<std::pair<T, U>>::iterator init = infirst; init<inlast; init++){
      bool inserted = false;
      for (typename std::vector<std::pair<T, U>>::iterator it = valfirst; it<vallast; it++){
        if ((*it).first<=(*init).first){
          inserted=true;
          if ((*it).second!=(*init).second) valArray.insert(it, *init);
          break;
        }
      }
      if (!inserted) valArray.insert(vallast, *init);
    }
  }
}

template<typename T> bool HelperFunctions::checkListVariable(std::vector<T> const& list, T const& var){ return (std::find(std::begin(list), std::end(list), var)!=std::end(list)); }
template<typename T> bool HelperFunctions::hasCommonElements(std::vector<T> const& list1, std::vector<T> const& list2){
  for (T const& el1:list1){
    for (T const& el2:list2){
      if (el1==el2) return true;
    }
  }
  return false;
}

template<typename T, typename U> void HelperFunctions::cleanUnorderedMap(std::unordered_map<T, U>& um){ for (auto& it:um){ delete it.second; it.second=nullptr; } }

// Non-negative, non-zero, positive-definite and NaN/Inf checkers
template<typename T> bool HelperFunctions::checkVarNonNegative(T const& val){ return (val>=0.); }
template<typename T> bool HelperFunctions::checkNonNegative(std::vector<T> const& vars, int ibegin, int iend){
  int ipos=0;
  for (T const& v:vars){
    if ((ibegin>=0 && ipos<ibegin) || (iend>=0 && ipos>=iend)) continue;
    ipos++;
    if (!checkVarNonNegative<T>(v)) return false;
  }
  return true;
}
template<typename T> bool HelperFunctions::checkVarNonZero(T const& val){ return (val!=0.); }
template<typename T> bool HelperFunctions::checkNonZero(std::vector<T> const& vars, int ibegin, int iend){
  int ipos=0;
  for (T const& v:vars){
    if ((ibegin>=0 && ipos<ibegin) || (iend>=0 && ipos>=iend)) continue;
    ipos++;
    if (!checkVarNonZero<T>(v)) return false;
  }
  return true;
}
template<typename T> bool HelperFunctions::checkVarPositiveDef(T const& val){ return (val>0.); }
template<typename T> bool HelperFunctions::checkPositiveDef(std::vector<T> const& vars, int ibegin, int iend){
  int ipos=0;
  for (T const& v:vars){
    if ((ibegin>=0 && ipos<ibegin) || (iend>=0 && ipos>=iend)) continue;
    ipos++;
    if (!checkVarPositiveDef<T>(v)) return false;
  }
  return true;
}
// Following https://stackoverflow.com/a/42138465
template<typename T> bool HelperFunctions::checkVarNanInf(T const& val){
  return !(std::isnan(val) || std::isinf(val) || (std::numeric_limits<T>::has_quiet_NaN && val==std::numeric_limits<T>::quiet_NaN()));
}
template<typename T> bool HelperFunctions::checkNanInf(std::vector<T> const& vars){
  for (T const& v:vars){ if (!checkVarNanInf<T>(v)) return false; }
  return true;
}

// Bit set and test
template<typename T> void HelperFunctions::set_bit(T& mask, unsigned int iBit, bool val){
  if (val) mask |= (1<<iBit);
  else if (test_bit(mask, iBit)){
    T tmp_mask = (1<<iBit);
    mask = mask ^ tmp_mask;
  }
}
template<typename T> bool HelperFunctions::test_bit(T const& mask, unsigned int iBit){ return (mask >> iBit) & 1; }

// Vector functions
template<typename T> void HelperFunctions::deltaR(T const& e1, T const& p1, T const& e2, T const& p2, T& res){
  T deta;
  deltaEta(e1, e2, deta);
  T dphi;
  deltaPhi(p1, p2, dphi);
  res = std::sqrt(std::pow(deta, 2) + std::pow(dphi, 2));
}
template<typename T> void HelperFunctions::deltaEta(T const& v1, T const& v2, T& res){
  res = v1 - v2;
}
template<typename T> void HelperFunctions::deltaPhi(T const& v1, T const& v2, T& res){
  res = v1 - v2;
  T const pi_val = TMath::Pi();
  T const twopi_val = TMath::Pi()*T(2);
  if (res>pi_val){ while (res>pi_val) res -= twopi_val; }
  else if (res<-pi_val){ while (res<-pi_val) res += twopi_val; }
}

// Histogram functions
template <typename T> double HelperFunctions::getHistogramIntegralAndError(T const* histo, int ix, int jx, bool useWidth, double* error){
  double res=0;
  double reserror=0;
  if (histo){
    if (!useWidth) res=histo->IntegralAndError(ix, jx, reserror, "");
    else{
      int xb[2]={ std::max(1, std::min(histo->GetNbinsX(), ix)), std::max(1, std::min(histo->GetNbinsX(), jx)) };

      double integralinside=0, integralerrorinside=0;
      double integraloutside=0, integralerroroutside=0;
      if (
        xb[0]>=ix && xb[1]<=jx
        ){
        res=histo->IntegralAndError(xb[0], xb[1], reserror, "width");
        integralinside=histo->IntegralAndError(xb[0], xb[1], integralerrorinside, "");
      }
      integraloutside=histo->IntegralAndError(ix, jx, integralerroroutside, "");

      res = res + integraloutside - integralinside;
      reserror = sqrt(std::max(0., pow(reserror, 2) + pow(integralerroroutside, 2) - pow(integralerrorinside, 2)));
    }
  }
  if (error) *error=reserror;
  return res;
}
template <typename T> double HelperFunctions::getHistogramIntegralAndError(T const* histo, int ix, int jx, int iy, int jy, bool useWidth, double* error){
  double res=0;
  double reserror=0;
  if (histo){
    if (!useWidth) res=histo->IntegralAndError(ix, jx, iy, jy, reserror, "");
    else{
      int xb[2]={ std::max(1, std::min(histo->GetNbinsX(), ix)), std::max(1, std::min(histo->GetNbinsX(), jx)) };
      int yb[2]={ std::max(1, std::min(histo->GetNbinsY(), iy)), std::max(1, std::min(histo->GetNbinsY(), jy)) };

      double integralinside=0, integralerrorinside=0;
      double integraloutside=0, integralerroroutside=0;
      if (
        xb[0]>=ix && xb[1]<=jx
        &&
        yb[0]>=iy && yb[1]<=jy
        ){
        res=histo->IntegralAndError(xb[0], xb[1], yb[0], yb[1], reserror, "width");
        integralinside=histo->IntegralAndError(xb[0], xb[1], yb[0], yb[1], integralerrorinside, "");
      }
      integraloutside=histo->IntegralAndError(ix, jx, iy, jy, integralerroroutside, "");

      res = res + integraloutside - integralinside;
      reserror = sqrt(std::max(0., pow(reserror, 2) + pow(integralerroroutside, 2) - pow(integralerrorinside, 2)));
    }
  }
  if (error) *error=reserror;
  return res;
}
template <typename T> double HelperFunctions::getHistogramIntegralAndError(T const* histo, int ix, int jx, int iy, int jy, int iz, int jz, bool useWidth, double* error){
  double res=0;
  double reserror=0;
  if (histo){
    if (!useWidth) res=histo->IntegralAndError(ix, jx, iy, jy, iz, jz, reserror, "");
    else{
      int xb[2]={ std::max(1, std::min(histo->GetNbinsX(), ix)), std::max(1, std::min(histo->GetNbinsX(), jx)) };
      int yb[2]={ std::max(1, std::min(histo->GetNbinsY(), iy)), std::max(1, std::min(histo->GetNbinsY(), jy)) };
      int zb[2]={ std::max(1, std::min(histo->GetNbinsZ(), iz)), std::max(1, std::min(histo->GetNbinsZ(), jz)) };

      double integralinside=0, integralerrorinside=0;
      double integraloutside=0, integralerroroutside=0;
      if (
        xb[0]>=ix && xb[1]<=jx
        &&
        yb[0]>=iy && yb[1]<=jy
        &&
        zb[0]>=iz && zb[1]<=jz
        ){
        res=histo->IntegralAndError(xb[0], xb[1], yb[0], yb[1], zb[0], zb[1], reserror, "width");
        integralinside=histo->IntegralAndError(xb[0], xb[1], yb[0], yb[1], zb[0], zb[1], integralerrorinside, "");
      }
      integraloutside=histo->IntegralAndError(ix, jx, iy, jy, iz, jz, integralerroroutside, "");

      res = res + integraloutside - integralinside;
      reserror = sqrt(std::max(0., pow(reserror, 2) + pow(integralerroroutside, 2) - pow(integralerrorinside, 2)));
    }
  }
  if (error) *error=reserror;
  return res;
}
template <typename T> void HelperFunctions::combineHistogramsByWeightedAverage(T const* h1, T const* h2, T*& hAssign, bool useNeff){
  std::vector<T const*> hlist;
  hlist.push_back(h1);
  hlist.push_back(h2);
  HelperFunctions::combineHistogramListByWeightedAverage<T>(hlist, hAssign, useNeff);
}
template <typename T> void HelperFunctions::findBinContentRange(std::vector<T*> hlist, float& bcmin, float& bcmax, bool includeBinErrors, bool includeOverUnderflows, bool onlyPositiveBins){
  bool firstHistogram = true;
  for (auto*& h:hlist){
    float tmpbcmin=0, tmpbcmax=0;
    findBinContentRange(h, tmpbcmin, tmpbcmax, includeBinErrors, includeOverUnderflows, onlyPositiveBins);
    if (firstHistogram){
      bcmin=tmpbcmin;
      bcmax=tmpbcmax;
      firstHistogram=false;
    }
    else{
      bcmin=std::min(bcmin, tmpbcmin);
      bcmax=std::max(bcmax, tmpbcmax);
    }
  }
}


// TGraph functions
template<typename T> TGraph* HelperFunctions::makeGraphFromPair(std::vector<std::pair<T, T>> points, TString name){
  if (points.empty()) return nullptr;
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
template<typename T> TGraphErrors* HelperFunctions::makeGraphSymErrFromPair(std::vector<std::pair<T, T>> points, std::vector<std::pair<T, T>> errors, TString name){
  if (points.empty() || errors.size()!=points.size()) return nullptr;
  unsigned int nbins = points.size();
  double* xy[4];
  for (unsigned int ix=0; ix<4; ix++) xy[ix] = new double[nbins];
  for (unsigned int bin=0; bin<nbins; bin++){
    xy[0][bin] = points[bin].first;
    xy[1][bin] = points[bin].second;
    xy[2][bin] = errors[bin].first;
    xy[3][bin] = errors[bin].second;
  }
  TGraphErrors* tg = new TGraphErrors(nbins, xy[0], xy[1], xy[2], xy[3]);
  tg->SetName(name);
  for (unsigned int ix=0; ix<4; ix++) delete[] xy[ix];
  return tg;
}
template<typename T> TGraphAsymmErrors* HelperFunctions::makeGraphAsymErrFromPair(std::vector<std::pair<T, T>> points, std::vector<std::pair<T, T>> errorDns, std::vector<std::pair<T, T>> errorUps, TString name){
  if (points.empty() || errorDns.size()!=points.size() || errorUps.size()!=points.size()) return nullptr;
  unsigned int nbins = points.size();
  double* xy[6];
  for (unsigned int ix=0; ix<6; ix++) xy[ix] = new double[nbins];
  for (unsigned int bin=0; bin<nbins; bin++){
    xy[0][bin] = points[bin].first;
    xy[1][bin] = points[bin].second;
    xy[2][bin] = errorDns[bin].first;
    xy[3][bin] = errorUps[bin].first;
    xy[4][bin] = errorDns[bin].second;
    xy[5][bin] = errorUps[bin].second;
  }
  TGraphAsymmErrors* tg = new TGraphAsymmErrors(nbins, xy[0], xy[1], xy[2], xy[3], xy[4], xy[5]);
  tg->SetName(name);
  for (unsigned int ix=0; ix<6; ix++) delete[] xy[ix];
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

// Function to extract all histograms from file
template<typename T> void HelperFunctions::extractObjectsFromDirectory(TDirectory* source, std::vector<T*>& objlist, TVar::VerbosityLevel verbosity){
  // Copy all objects and subdirs of directory source as a subdir of the current directory
  TDirectory* target = gDirectory;
  if (verbosity>=TVar::INFO) source->ls();
  source->cd();
  // Loop on all entries of this directory
  TKey* key;
  TIter nextkey(source->GetListOfKeys());
  std::vector<TString> copiedKeys;
  while ((key = (TKey*) nextkey())){
    const char* classname = key->GetClassName();
    TClass* cl = gROOT->GetClass(classname);
    if (!cl) continue;
    if (cl->InheritsFrom(TDirectory::Class())){
      source->cd(key->GetName());
      TDirectory* subdir = gDirectory;
      source->cd();
      extractObjectsFromDirectory<T>(subdir, objlist, verbosity);
    }
    else if (cl->InheritsFrom(T::Class())){
      T* obj = (T*) source->Get(key->GetName());
      if (!obj) continue;
      TString objname = obj->GetName();
      TString keyname = key->GetName();
      if ((objname=="Graph" && keyname!="Graph") || objname=="") objname = keyname; // Holy jumping monkeys for fake rates
      obj->SetName(objname);
      if (!HelperFunctions::checkListVariable(copiedKeys, keyname)){
        copiedKeys.push_back(keyname);
        objlist.push_back(obj);
      }
    }
  }
  target->cd();
}
template void HelperFunctions::extractObjectsFromDirectory<TH1F>(TDirectory* source, std::vector<TH1F*>& objlist, TVar::VerbosityLevel verbosity);
template void HelperFunctions::extractObjectsFromDirectory<TH2F>(TDirectory* source, std::vector<TH2F*>& objlist, TVar::VerbosityLevel verbosity);
template void HelperFunctions::extractObjectsFromDirectory<TH3F>(TDirectory* source, std::vector<TH3F*>& objlist, TVar::VerbosityLevel verbosity);
template void HelperFunctions::extractObjectsFromDirectory<TH1D>(TDirectory* source, std::vector<TH1D*>& objlist, TVar::VerbosityLevel verbosity);
template void HelperFunctions::extractObjectsFromDirectory<TH2D>(TDirectory* source, std::vector<TH2D*>& objlist, TVar::VerbosityLevel verbosity);
template void HelperFunctions::extractObjectsFromDirectory<TH3D>(TDirectory* source, std::vector<TH3D*>& objlist, TVar::VerbosityLevel verbosity);
// Overloads for TGraph that can be used just like histograms
template void HelperFunctions::extractObjectsFromDirectory<TGraphAsymmErrors>(TDirectory* source, std::vector<TGraphAsymmErrors*>& objlist, TVar::VerbosityLevel verbosity);
template void HelperFunctions::extractObjectsFromDirectory<TGraphErrors>(TDirectory* source, std::vector<TGraphErrors*>& objlist, TVar::VerbosityLevel verbosity);
template void HelperFunctions::extractObjectsFromDirectory<TGraph>(TDirectory* source, std::vector<TGraph*>& objlist, TVar::VerbosityLevel verbosity);
// Overloads for TSpline that can be used just like histograms
template void HelperFunctions::extractObjectsFromDirectory<TSpline3>(TDirectory* source, std::vector<TSpline3*>& objlist, TVar::VerbosityLevel verbosity);
// Overloads for TCanvas that can be used just like histograms
template void HelperFunctions::extractObjectsFromDirectory<TCanvas>(TDirectory* source, std::vector<TCanvas*>& objlist, TVar::VerbosityLevel verbosity);

template <typename T, typename U> void HelperFunctions::getGenericHistogramSlice(T*& res, U const* histo, unsigned char XDirection, int iy, int jy, TString newname){
  using namespace MELAStreamHelpers;

  if (!histo || XDirection>=2) return;
  if (newname=="") newname=Form("Slice_%s_%i_%i_%s", (XDirection==0 ? "X" : "Y"), iy, jy, histo->GetName());

  const TAxis* xaxis=histo->GetXaxis();
  const TAxis* yaxis=histo->GetYaxis();
  std::vector<float> bins;
  if (XDirection==0){
    for (int i=1; i<=xaxis->GetNbins()+1; i++) bins.push_back(xaxis->GetBinLowEdge(i));
    iy = std::max(0, iy); jy = std::min(yaxis->GetNbins()+1, jy);
  }
  else{
    for (int i=1; i<=yaxis->GetNbins()+1; i++) bins.push_back(yaxis->GetBinLowEdge(i));
    iy = std::max(0, iy); jy = std::min(xaxis->GetNbins()+1, jy);
  }
  if (iy>jy) MELAerr << "HelperFunctions::getGenericHistogramSlice: iy>jy!" << std::endl;

  res = new T(newname, "", bins.size()-1, bins.data());

  if (XDirection==0){
    for (int ii=0; ii<=xaxis->GetNbins()+1; ii++){
      double integral=0, integralerror=0;
      integral = getHistogramIntegralAndError(histo, ii, ii, iy, jy, false, &integralerror);
      res->SetBinContent(ii, integral);
      res->SetBinError(ii, integralerror);
    }
  }
  else{
    for (int ii=0; ii<=yaxis->GetNbins()+1; ii++){
      double integral=0, integralerror=0;
      integral = getHistogramIntegralAndError(histo, iy, jy, ii, ii, false, &integralerror);
      res->SetBinContent(ii, integral);
      res->SetBinError(ii, integralerror);
    }
  }
}
template <typename T, typename U> void HelperFunctions::getGenericHistogramSlice(T*& res, U const* histo, unsigned char XDirection, int iy, int jy, int iz, int jz, TString newname){
  using namespace MELAStreamHelpers;

  if (!histo || XDirection>=3) return;
  if (newname=="") newname=Form("Slice_%s_%i_%i_%i_%i_%s", (XDirection==0 ? "X" : (XDirection==1 ? "Y" : "Z")), iy, jy, iz, jz, histo->GetName());

  const TAxis* xaxis;
  const TAxis* yaxis;
  const TAxis* zaxis;
  std::vector<float> bins;
  if (XDirection==0){
    xaxis=histo->GetXaxis();
    yaxis=histo->GetYaxis();
    zaxis=histo->GetZaxis();
  }
  else if (XDirection==1){
    xaxis=histo->GetYaxis();
    yaxis=histo->GetZaxis();
    zaxis=histo->GetXaxis();
  }
  else{
    xaxis=histo->GetZaxis();
    yaxis=histo->GetXaxis();
    zaxis=histo->GetYaxis();
  }
  for (int i=1; i<=xaxis->GetNbins()+1; i++) bins.push_back(xaxis->GetBinLowEdge(i));
  iy = std::max(0, iy); jy = std::min(yaxis->GetNbins()+1, jy);
  iz = std::max(0, iz); jz = std::min(zaxis->GetNbins()+1, jz);
  if (iy>jy) MELAerr << "HelperFunctions::getHistogramSlice: iy>jy!" << std::endl;
  if (iz>jz) MELAerr << "HelperFunctions::getHistogramSlice: iz>jz!" << std::endl;
  
  res = new T(newname, "", bins.size()-1, bins.data());

  for (int ii=0; ii<=xaxis->GetNbins()+1; ii++){
    double integral=0, integralerror=0;
    int IX, JX, IY, JY, IZ, JZ;
    if (XDirection==0){
      IX=ii; JX=ii;
      IY=iy; JY=jy;
      IZ=iz; JZ=jz;
    }
    else if (XDirection==1){
      IX=iz; JX=jz;
      IY=ii; JY=ii;
      IZ=iy; JZ=jy;
    }
    else{
      IX=iy; JX=jy;
      IY=iz; JY=jz;
      IZ=ii; JZ=ii;
    }
    integral = getHistogramIntegralAndError(histo, IX, JX, IY, JY, IZ, JZ, false, &integralerror);
    res->SetBinContent(ii, integral);
    res->SetBinError(ii, integralerror);
  }
}
template <typename T, typename U> void HelperFunctions::getGenericHistogramSlice(T*& res, U const* histo, unsigned char XDirection, unsigned char YDirection, int iz, int jz, TString newname){
  using namespace MELAStreamHelpers;

  if (!histo || XDirection==YDirection || XDirection>=3 || YDirection>=3) return;
  if (newname=="") newname=Form("Slice_%s%s_%i_%i_%s", (XDirection==0 ? "X" : (XDirection==1 ? "Y" : "Z")), (YDirection==0 ? "X" : (YDirection==1 ? "Y" : "Z")), iz, jz, histo->GetName());

  unsigned char ZDirection=3-XDirection-YDirection; // 0+1+2=3
  const TAxis* xaxis;
  const TAxis* yaxis;
  const TAxis* zaxis;
  std::vector<float> xbins, ybins;
  if (XDirection==0) xaxis=histo->GetXaxis();
  else if (XDirection==1) xaxis=histo->GetYaxis();
  else xaxis=histo->GetZaxis();
  if (YDirection==0) yaxis=histo->GetXaxis();
  else if (YDirection==1) yaxis=histo->GetYaxis();
  else yaxis=histo->GetZaxis();
  if (ZDirection==0) zaxis=histo->GetXaxis();
  else if (ZDirection==1) zaxis=histo->GetYaxis();
  else zaxis=histo->GetZaxis();

  for (int i=1; i<=xaxis->GetNbins()+1; i++) xbins.push_back(xaxis->GetBinLowEdge(i));
  for (int i=1; i<=yaxis->GetNbins()+1; i++) ybins.push_back(yaxis->GetBinLowEdge(i));
  iz = std::max(0, iz); std::min(zaxis->GetNbins()+1, jz);
  
  res = new T(newname, "", xbins.size()-1, xbins.data(), ybins.size()-1, ybins.data());

  for (int ii=0; ii<=xaxis->GetNbins()+1; ii++){
    for (int jj=0; jj<=yaxis->GetNbins()+1; jj++){
      double integral=0, integralerror=0;
      int IX=0, JX=0, IY=0, JY=0, IZ=0, JZ=0;
      if (XDirection==0){ IX=ii; JX=ii; }
      else if (XDirection==1){ IY=ii; JY=ii; }
      else{ IZ=ii; JZ=ii; }
      if (YDirection==0){ IX=jj; JX=jj; }
      else if (YDirection==1){ IY=jj; JY=jj; }
      else{ IZ=jj; JZ=jj; }
      if (ZDirection==0){ IX=iz; JX=jz; }
      else if (ZDirection==1){ IY=iz; JY=jz; }
      else{ IZ=iz; JZ=jz; }
      integral = getHistogramIntegralAndError(histo, IX, JX, IY, JY, IZ, JZ, false, &integralerror);
      res->SetBinContent(ii, jj, integral);
      res->SetBinError(ii, jj, integralerror);
    }
  }
}

/****************************************************************/
// Explicit instantiations
template void HelperFunctions::appendVector<TString>(std::vector<TString>& a, std::vector<TString> const& b);

template void HelperFunctions::addByLowest<SimpleEntry>(std::vector<SimpleEntry>& valArray, SimpleEntry val, bool unique);
template void HelperFunctions::addByLowest<int>(std::vector<int>& valArray, int val, bool unique);
template void HelperFunctions::addByLowest<float>(std::vector<float>& valArray, float val, bool unique);
template void HelperFunctions::addByLowest<double>(std::vector<double>& valArray, double val, bool unique);
template void HelperFunctions::addByLowest<SimpleEntry, int>(std::vector<std::pair<SimpleEntry, int>>& valArray, SimpleEntry val, int index);
template void HelperFunctions::addByLowest<double, int>(std::vector<std::pair<double, int>>& valArray, double val, int index);
template void HelperFunctions::addByLowest<double, double>(std::vector<std::pair<double, double>>& valArray, double val, double index);
template void HelperFunctions::addByLowest<double, int>(std::vector<std::pair<double, int>>& valArray, std::vector<std::pair<double, int>>& inArray, bool consecutive, bool inputordered);
template void HelperFunctions::addByLowest<double, double>(std::vector<std::pair<double, double>>& valArray, std::vector<std::pair<double, double>>& inArray, bool consecutive, bool inputordered);

template void HelperFunctions::addByHighest<SimpleEntry>(std::vector<SimpleEntry>& valArray, SimpleEntry val, bool unique);
template void HelperFunctions::addByHighest<int>(std::vector<int>& valArray, int val, bool unique);
template void HelperFunctions::addByHighest<float>(std::vector<float>& valArray, float val, bool unique);
template void HelperFunctions::addByHighest<double>(std::vector<double>& valArray, double val, bool unique);
template void HelperFunctions::addByHighest<SimpleEntry, int>(std::vector<std::pair<SimpleEntry, int>>& valArray, SimpleEntry val, int index);
template void HelperFunctions::addByHighest<double, int>(std::vector<std::pair<double, int>>& valArray, double val, int index);
template void HelperFunctions::addByHighest<double, double>(std::vector<std::pair<double, double>>& valArray, double val, double index);
template void HelperFunctions::addByHighest<double, int>(std::vector<std::pair<double, int>>& valArray, std::vector<std::pair<double, int>>& inArray, bool consecutive, bool inputordered);
template void HelperFunctions::addByHighest<double, double>(std::vector<std::pair<double, double>>& valArray, std::vector<std::pair<double, double>>& inArray, bool consecutive, bool inputordered);

template bool HelperFunctions::checkListVariable<std::string>(const std::vector<std::string>& list, const std::string& var);
template bool HelperFunctions::checkListVariable<TString>(const std::vector<TString>& list, const TString& var);
template bool HelperFunctions::checkListVariable<double>(const std::vector<double>& list, const double& var);

template bool HelperFunctions::checkNonNegative<short>(std::vector<short> const& vars, int ibegin, int iend);
template bool HelperFunctions::checkNonNegative<unsigned int>(std::vector<unsigned int> const& vars, int ibegin, int iend);
template bool HelperFunctions::checkNonNegative<int>(std::vector<int> const& vars, int ibegin, int iend);
template bool HelperFunctions::checkNonNegative<float>(std::vector<float> const& vars, int ibegin, int iend);
template bool HelperFunctions::checkNonNegative<double>(std::vector<double> const& vars, int ibegin, int iend);

template bool HelperFunctions::checkNonZero<short>(std::vector<short> const& vars, int ibegin, int iend);
template bool HelperFunctions::checkNonZero<unsigned int>(std::vector<unsigned int> const& vars, int ibegin, int iend);
template bool HelperFunctions::checkNonZero<int>(std::vector<int> const& vars, int ibegin, int iend);
template bool HelperFunctions::checkNonZero<float>(std::vector<float> const& vars, int ibegin, int iend);
template bool HelperFunctions::checkNonZero<double>(std::vector<double> const& vars, int ibegin, int iend);

template bool HelperFunctions::checkPositiveDef<short>(std::vector<short> const& vars, int ibegin, int iend);
template bool HelperFunctions::checkPositiveDef<unsigned int>(std::vector<unsigned int> const& vars, int ibegin, int iend);
template bool HelperFunctions::checkPositiveDef<int>(std::vector<int> const& vars, int ibegin, int iend);
template bool HelperFunctions::checkPositiveDef<float>(std::vector<float> const& vars, int ibegin, int iend);
template bool HelperFunctions::checkPositiveDef<double>(std::vector<double> const& vars, int ibegin, int iend);

template bool HelperFunctions::checkNanInf<short>(std::vector<short> const& vars);
template bool HelperFunctions::checkNanInf<unsigned int>(std::vector<unsigned int> const& vars);
template bool HelperFunctions::checkNanInf<int>(std::vector<int> const& vars);
template bool HelperFunctions::checkNanInf<float>(std::vector<float> const& vars);
template bool HelperFunctions::checkNanInf<double>(std::vector<double> const& vars);

template double HelperFunctions::getHistogramIntegralAndError<TH1F>(TH1F const* histo, int ix, int jx, bool useWidth, double* error);
template double HelperFunctions::getHistogramIntegralAndError<TH1D>(TH1D const* histo, int ix, int jx, bool useWidth, double* error);
template double HelperFunctions::getHistogramIntegralAndError<TH2F>(TH2F const* histo, int ix, int jx, int iy, int jy, bool useWidth, double* error);
template double HelperFunctions::getHistogramIntegralAndError<TH2D>(TH2D const* histo, int ix, int jx, int iy, int jy, bool useWidth, double* error);
template double HelperFunctions::getHistogramIntegralAndError<TH3F>(TH3F const* histo, int ix, int jx, int iy, int jy, int iz, int jz, bool useWidth, double* error);
template double HelperFunctions::getHistogramIntegralAndError<TH3D>(TH3D const* histo, int ix, int jx, int iy, int jy, int iz, int jz, bool useWidth, double* error);

template TGraph* HelperFunctions::makeGraphFromPair<float>(std::vector<std::pair<float, float>> points, TString name);
template TGraphErrors* HelperFunctions::makeGraphSymErrFromPair<float>(std::vector<std::pair<float, float>> points, std::vector<std::pair<float, float>> errors, TString name);
template TGraphAsymmErrors* HelperFunctions::makeGraphAsymErrFromPair<float>(std::vector<std::pair<float, float>> points, std::vector<std::pair<float, float>> errorDns, std::vector<std::pair<float, float>> errorUps, TString name);
template TGraph* HelperFunctions::makeGraphFromPair<double>(std::vector<std::pair<double, double>> points, TString name);
template TGraphErrors* HelperFunctions::makeGraphSymErrFromPair<double>(std::vector<std::pair<double, double>> points, std::vector<std::pair<double, double>> errors, TString name);
template TGraphAsymmErrors* HelperFunctions::makeGraphAsymErrFromPair<double>(std::vector<std::pair<double, double>> points, std::vector<std::pair<double, double>> errorDns, std::vector<std::pair<double, double>> errorUps, TString name);

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
