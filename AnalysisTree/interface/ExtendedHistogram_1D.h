#ifndef EXTENDEDHISTOGRAM_1D_H
#define EXTENDEDHISTOGRAM_1D_H

#include "ExtendedHistogram.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TGraphErrors.h"


template<typename Hist_t=TH1F> class ExtendedHistogram_1D : public ExtendedHistogram{
protected:
  ExtendedBinning xbinning;
  Hist_t* histo;
  TProfile* prof_x;
  
public:
  ExtendedHistogram_1D();
  ExtendedHistogram_1D(const TString name_, const TString title_);
  ExtendedHistogram_1D(const TString name_, const TString title_, const ExtendedBinning& xbinning_);
  ExtendedHistogram_1D(ExtendedHistogram_1D<Hist_t> const& other);
  virtual ~ExtendedHistogram_1D();

  void swap(ExtendedHistogram_1D<Hist_t>& other);
  ExtendedHistogram_1D<Hist_t>& operator=(const ExtendedHistogram_1D<Hist_t>& other);

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

  TGraphErrors* getGraph(TString newname="") const;

  void fill(double x, double wgt=1.);

  void rebin(ExtendedBinning const* binningX);

  void constructFromTree(TTree* tree, float& xvar, float& weight, bool* flag=nullptr, ExtendedBinning const* binningX=nullptr);

  static ExtendedHistogram_1D<Hist_t> divideHistograms(ExtendedHistogram_1D<Hist_t> const& h1, ExtendedHistogram_1D<Hist_t> const& h2, bool useEffErr, TString newname="");

  static void averageHistograms(ExtendedHistogram_1D<Hist_t>& hTarget, ExtendedHistogram_1D<Hist_t> const& h2, bool useNeff=false);
  static void averageHistograms(ExtendedHistogram_1D<Hist_t>& hTarget, std::vector<ExtendedHistogram_1D<Hist_t> const*> const& hList, bool useNeff=false);

};

typedef ExtendedHistogram_1D<TH1F> ExtendedHistogram_1D_f;
typedef ExtendedHistogram_1D<TH1D> ExtendedHistogram_1D_d;


#include "ExtendedHistogram_1D.hh"


#endif
