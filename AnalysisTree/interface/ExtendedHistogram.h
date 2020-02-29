#ifndef EXTENDEDHISTOGRAM_H
#define EXTENDEDHISTOGRAM_H

#include "ExtendedBinning.h"
#include "TString.h"


class ExtendedHistogram{
protected:
  TString name;
  TString title;

public:
  ExtendedHistogram();
  ExtendedHistogram(const TString name_, const TString title_);
  ExtendedHistogram(ExtendedHistogram const& other);
  virtual ~ExtendedHistogram();

  const TString& getName() const;
  TString getName();
  const TString& getTitle() const;
  TString getTitle();

  virtual void setNameTitle(const TString name_, const TString title_="");
  virtual void setBinning(const ExtendedBinning& binning, const int xyz=0, const TString label="")=0; // xyz=0,1,2 for x, y, z
  virtual ExtendedBinning const& getBinning(const int xyz=0)=0;
  virtual void build()=0;
  virtual void resetProfiles()=0;

};


#endif
