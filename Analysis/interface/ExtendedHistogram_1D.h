#ifndef EXTENDEDHISTOGRAM_1D_H
#define EXTENDEDHISTOGRAM_1D_H

#include "ExtendedHistogram.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TGraphErrors.h"


class ExtendedHistogram_1D : public ExtendedHistogram{
protected:
  ExtendedBinning xbinning;
  TH1F* histo;
  TProfile* prof_x;
  
public:
  ExtendedHistogram_1D();
  ExtendedHistogram_1D(const TString name_, const TString title_);
  ExtendedHistogram_1D(const TString name_, const TString title_, const ExtendedBinning& xbinning_);
  ExtendedHistogram_1D(ExtendedHistogram_1D const& other);
  virtual ~ExtendedHistogram_1D();

  void swap(ExtendedHistogram_1D& other);
  ExtendedHistogram_1D& operator=(const ExtendedHistogram_1D& other);

  void setNameTitle(const TString name_, const TString title_="");
  virtual void setBinning(const ExtendedBinning& binning, const int xyz=0, const TString label="");
  virtual void build();

  TH1F*& getHistogram(){ return histo; }
  const TH1F* getHistogram()const{ return histo; }

  TProfile*& getProfileX(){ return prof_x; }
  const TProfile* getProfileX() const{ return prof_x; }

  TGraphErrors* getGraph(TString newname="") const;

  void fill(double x, double wgt=1.);

  static ExtendedHistogram_1D divideHistograms(ExtendedHistogram_1D const& h1, ExtendedHistogram_1D const& h2, bool useEffErr, TString newname="");

};


#endif
