#ifndef EXTENDEDBINNING_H
#define EXTENDEDBINNING_H

#include <vector>
#include "TString.h"


class ExtendedBinning{
protected:
  std::vector<double> vbinlow; // Size=Nbins+1
  TString label;
  bool isvalid;

public:
  ExtendedBinning();
  ExtendedBinning(const unsigned int nbins, const double min, const double max, const TString label_=""); // Uniform constructor
  ExtendedBinning(const double* abinlow, const TString label_="");
  ExtendedBinning(const std::vector<double>& vbinlow_, const TString label_="");

  bool isValid()const{ return isvalid; }

  void setLabel(const TString label_);

  double* getBinning();
  const double* getBinning() const;
  const unsigned int getNbins() const;
  const TString getLabel(){ return label; }

  int getBin(double val)const; // = [ -1,0,...,vbinlow.size() ]
  double getBinLowEdge(const int bin)const;
  double getBinHighEdge(const int bin)const;

};


#endif

