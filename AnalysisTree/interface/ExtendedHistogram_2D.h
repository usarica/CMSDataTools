#ifndef EXTENDEDHISTOGRAM_2D_H
#define EXTENDEDHISTOGRAM_2D_H

#include "ExtendedHistogram.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TGraphErrors.h"


template<typename Hist_t=TH2F> class ExtendedHistogram_2D : public ExtendedHistogram{
protected:
  ExtendedBinning xbinning;
  ExtendedBinning ybinning;
  Hist_t* histo;
  TProfile* prof_x;
  TProfile* prof_y;

public:
  ExtendedHistogram_2D();
  ExtendedHistogram_2D(const TString name_, const TString title_);
  ExtendedHistogram_2D(const TString name_, const TString title_, const ExtendedBinning& xbinning_, const ExtendedBinning& ybinning_);
  ExtendedHistogram_2D(ExtendedHistogram_2D<Hist_t> const& other);
  virtual ~ExtendedHistogram_2D();

  void swap(ExtendedHistogram_2D<Hist_t>& other);
  ExtendedHistogram_2D<Hist_t>& operator=(const ExtendedHistogram_2D<Hist_t>& other);

  void setNameTitle(const TString name_, const TString title_="");
  virtual void setBinning(const ExtendedBinning& binning, const int xyz=0, const TString label="");
  virtual ExtendedBinning const& getBinning(const int xyz=0);
  virtual void build();
  virtual void reset();
  virtual void resetProfiles();

  Hist_t*& getHistogram(){ return histo; }
  const Hist_t* getHistogram() const{ return histo; }
  Hist_t* getCumulantHistogram(TString newname="") const;

  TProfile*& getProfileX(){ return prof_x; }
  const TProfile* getProfileX() const{ return prof_x; }

  TProfile*& getProfileY(){ return prof_y; }
  const TProfile* getProfileY() const{ return prof_y; }

  void fill(double x, double y, double wgt=1.);

  void rebin(ExtendedBinning const* binningX, ExtendedBinning const* binningY, signed char condDim=-1);

  void constructFromTree(TTree* tree, float& xvar, float& yvar, float& weight, bool* flag=nullptr, ExtendedBinning const* binningX=nullptr, ExtendedBinning const* binningY=nullptr);

  static ExtendedHistogram_2D<Hist_t> divideHistograms(ExtendedHistogram_2D<Hist_t> const& h1, ExtendedHistogram_2D<Hist_t> const& h2, bool useEffErr, TString newname="");

  static void averageHistograms(ExtendedHistogram_2D<Hist_t>& hTarget, ExtendedHistogram_2D<Hist_t> const& h2, bool useNeff=false);
  static void averageHistograms(ExtendedHistogram_2D<Hist_t>& hTarget, std::vector<ExtendedHistogram_2D<Hist_t> const*> const& hList, bool useNeff=false);

};

typedef ExtendedHistogram_2D<TH2F> ExtendedHistogram_2D_f;
typedef ExtendedHistogram_2D<TH2D> ExtendedHistogram_2D_d;


#include "ExtendedHistogram_2D.hh"


#endif
