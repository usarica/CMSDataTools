#include "ExtendedHistogram_1D.h"


ExtendedHistogram_1D::ExtendedHistogram_1D() : ExtendedHistogram(), histo(nullptr), prof_x(nullptr){}
ExtendedHistogram_1D::ExtendedHistogram_1D(const TString name_, const TString title_) : ExtendedHistogram(name_, title_), histo(nullptr), prof_x(nullptr){}
ExtendedHistogram_1D::~ExtendedHistogram_1D(){
  delete histo;
  delete prof_x;
}

void ExtendedHistogram_1D::setBinning(const ExtendedBinning& binning, const int xyz, const TString label){
  xbinning = binning;
  xbinning.setLabel(label);
}
void ExtendedHistogram_1D::build(){
  if (xbinning.isValid()){
    const double* xbins = xbinning.getBinning();
    const int nbins = xbinning.getNbins();
    histo = new TH1F(name, title, nbins, xbins); histo->GetXaxis()->SetTitle(xbinning.getLabel()); histo->Sumw2();
    prof_x = new TProfile(name, title, nbins, xbins); prof_x->GetXaxis()->SetTitle(xbinning.getLabel()); prof_x->Sumw2();
  }
}

void ExtendedHistogram_1D::fill(double x, double wgt){
  if (histo && prof_x){
    histo->Fill(x, wgt);
    prof_x->Fill(x, x, wgt);
  }
}
