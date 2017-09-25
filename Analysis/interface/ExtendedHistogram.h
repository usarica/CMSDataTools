#ifndef EXTENDEDHISTOGRAM_H
#define EXTENDEDHISTOGRAM_H

#include "ExtendedBinning.h"
#include "TString.h"


class ExtendedHistogram{
protected:
  TString name, title;


public:
  ExtendedHistogram();
  ExtendedHistogram(const TString name_, const TString title_);
  virtual ~ExtendedHistogram();

  virtual void setBinning(const ExtendedBinning& binning, const int xyz=0, const TString label="")=0; // xyz=0,1,2 for x, y, z
  virtual void build()=0;

};


#endif
