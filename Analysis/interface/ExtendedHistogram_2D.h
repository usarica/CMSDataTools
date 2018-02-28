#ifndef EXTENDEDHISTOGRAM_2D_H
#define EXTENDEDHISTOGRAM_2D_H

#include "ExtendedHistogram.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TGraphErrors.h"


class ExtendedHistogram_2D : public ExtendedHistogram{
protected:
  ExtendedBinning xbinning;
  ExtendedBinning ybinning;
  TH2F* histo;
  TProfile* prof_x;
  TProfile* prof_y;

public:
  ExtendedHistogram_2D();
  ExtendedHistogram_2D(const TString name_, const TString title_);
  ExtendedHistogram_2D(const TString name_, const TString title_, const ExtendedBinning& xbinning_, const ExtendedBinning& ybinning_);
  ExtendedHistogram_2D(ExtendedHistogram_2D const& other);
  virtual ~ExtendedHistogram_2D();

  void swap(ExtendedHistogram_2D& other);
  ExtendedHistogram_2D& operator=(const ExtendedHistogram_2D& other);

  void setNameTitle(const TString name_, const TString title_="");
  virtual void setBinning(const ExtendedBinning& binning, const int xyz=0, const TString label="");
  virtual void build();

  TH2F*& getHistogram(){ return histo; }
  const TH2F* getHistogram() const{ return histo; }
  TH2F* getCumulantHistogram(TString newname="") const;

  TProfile*& getProfileX(){ return prof_x; }
  const TProfile* getProfileX() const{ return prof_x; }

  TProfile*& getProfileY(){ return prof_y; }
  const TProfile* getProfileY() const{ return prof_y; }

  void fill(double x, double y, double wgt=1.);

  static ExtendedHistogram_2D divideHistograms(ExtendedHistogram_2D const& h1, ExtendedHistogram_2D const& h2, bool useEffErr, TString newname="");

};


#endif
