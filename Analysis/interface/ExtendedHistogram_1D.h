#ifndef EXTENDEDHISTOGRAM_1D_H
#define EXTENDEDHISTOGRAM_1D_H

#include "ExtendedHistogram.h"
#include "TH1F.h"
#include "TProfile.h"


class ExtendedHistogram_1D : public ExtendedHistogram{
protected:
  ExtendedBinning xbinning;
  TH1F* histo;
  TProfile* prof_x;


public:
  ExtendedHistogram_1D();
  ExtendedHistogram_1D(const TString name_, const TString title_);
  virtual ~ExtendedHistogram_1D();

  virtual void setBinning(const ExtendedBinning& binning, const int xyz=0, const TString label="");
  virtual void build();

  TH1F*& getHistogram(){ return histo; }
  const TH1F* getHistogram()const{ return histo; }

  TProfile*& getProfileX(){ return prof_x; }
  const TProfile* getProfileX()const{ return prof_x; }

  void fill(double x, double wgt=1.);

};


#endif
