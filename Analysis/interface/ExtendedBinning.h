#ifndef EXTENDEDBINNING_H
#define EXTENDEDBINNING_H

#include <vector>
#include <utility>
#include "TString.h"


class ExtendedBinning{
protected:
  std::vector<double> vbinlow; // Size=Nbins+1
  TString label;

public:
  ExtendedBinning(const TString label_="");
  ExtendedBinning(const unsigned int nbins, const double min, const double max, const TString label_=""); // Uniform constructor
  ExtendedBinning(const double* abinlow, const TString label_="");
  ExtendedBinning(const std::vector<double>& vbinlow_, const TString label_="");

  bool isValid() const;

  void setLabel(const TString label_);
  TString getLabel() const;

  double* getBinning();
  const double* getBinning() const;
  std::vector<double> getBinningVector();
  const std::vector<double>& getBinningVector() const;
  template<typename T> std::vector<std::pair<T,T>> getBoundaryPairsList() const;
  unsigned int getNbins() const;

  int getBin(double val) const; // = [ -1,0,...,vbinlow.size() ]
  double getBinLowEdge(const int bin) const;
  double getBinHighEdge(const int bin) const;

  void addBinBoundary(double boundary);

};

template<typename T> std::vector<std::pair<T, T>> ExtendedBinning::getBoundaryPairsList() const{
  std::vector<std::pair<T, T>> res;
  if (vbinlow.size()>1) res.reserve(vbinlow.size()-1);
  for (auto it=vbinlow.begin(); it!=vbinlow.end()-1; it++){
    auto it_next = it+1;
    res.push_back(std::pair<T, T>(T(*it), T(*it_next)));
  }
  return res;
}
template std::vector<std::pair<float, float>> ExtendedBinning::getBoundaryPairsList<float>() const;
template std::vector<std::pair<double, double>> ExtendedBinning::getBoundaryPairsList<double>() const;


#endif

