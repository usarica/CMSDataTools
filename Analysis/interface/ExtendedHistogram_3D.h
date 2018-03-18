#ifndef EXTENDEDHISTOGRAM_3D_H
#define EXTENDEDHISTOGRAM_3D_H

#include "ExtendedHistogram.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TGraphErrors.h"


class ExtendedHistogram_3D : public ExtendedHistogram{
protected:
  ExtendedBinning xbinning;
  ExtendedBinning ybinning;
  ExtendedBinning zbinning;
  TH3F* histo;
  TProfile* prof_x;
  TProfile* prof_y;
  TProfile* prof_z;

public:
  ExtendedHistogram_3D();
  ExtendedHistogram_3D(const TString name_, const TString title_);
  ExtendedHistogram_3D(const TString name_, const TString title_, const ExtendedBinning& xbinning_, const ExtendedBinning& ybinning_, const ExtendedBinning& zbinning_);
  ExtendedHistogram_3D(ExtendedHistogram_3D const& other);
  virtual ~ExtendedHistogram_3D();

  void swap(ExtendedHistogram_3D& other);
  ExtendedHistogram_3D& operator=(const ExtendedHistogram_3D& other);

  void setNameTitle(const TString name_, const TString title_="");
  virtual void setBinning(const ExtendedBinning& binning, const int xyz=0, const TString label="");
  virtual void build();
  virtual void reset();

  TH3F*& getHistogram(){ return histo; }
  const TH3F* getHistogram() const{ return histo; }
  TH3F* getCumulantHistogram(TString newname="") const;

  TProfile*& getProfileX(){ return prof_x; }
  const TProfile* getProfileX() const{ return prof_x; }

  TProfile*& getProfileY(){ return prof_y; }
  const TProfile* getProfileY() const{ return prof_y; }

  TProfile*& getProfileZ(){ return prof_z; }
  const TProfile* getProfileZ() const{ return prof_z; }

  void fill(double x, double y, double z, double wgt=1.);

  void rebin(ExtendedBinning const* binningX, ExtendedBinning const* binningY, ExtendedBinning const* binningZ, signed char condDim=-1);

  void constructFromTree(TTree* tree, float& xvar, float& yvar, float& zvar, float& weight, bool* flag=nullptr, ExtendedBinning const* binningX=nullptr, ExtendedBinning const* binningY=nullptr, ExtendedBinning const* binningZ=nullptr);

  static ExtendedHistogram_3D divideHistograms(ExtendedHistogram_3D const& h1, ExtendedHistogram_3D const& h2, bool useEffErr, TString newname="");

  static void averageHistograms(ExtendedHistogram_3D& hTarget, ExtendedHistogram_3D const& h2);

};


#endif
