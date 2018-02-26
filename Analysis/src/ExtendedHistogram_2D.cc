#include "HelperFunctions.h"
#include "ExtendedHistogram_2D.h"


ExtendedHistogram_2D::ExtendedHistogram_2D() : ExtendedHistogram(), histo(nullptr), prof_x(nullptr), prof_y(nullptr){}
ExtendedHistogram_2D::ExtendedHistogram_2D(const TString name_, const TString title_) : ExtendedHistogram(name_, title_), histo(nullptr), prof_x(nullptr), prof_y(nullptr){}
ExtendedHistogram_2D::ExtendedHistogram_2D(const TString name_, const TString title_, const ExtendedBinning& xbinning_, const ExtendedBinning& ybinning_) :
  ExtendedHistogram(name_, title_),
  xbinning(xbinning_),
  ybinning(ybinning_),
  histo(nullptr), prof_x(nullptr), prof_y(nullptr)
{
  build();
}
ExtendedHistogram_2D::ExtendedHistogram_2D(ExtendedHistogram_2D const& other) : ExtendedHistogram(other), xbinning(other.xbinning), ybinning(other.ybinning), histo(nullptr), prof_x(nullptr), prof_y(nullptr){
  if (other.histo) histo = new TH2F(*(other.histo));
  if (other.prof_x) prof_x = new TProfile(*(other.prof_x));
  if (other.prof_y) prof_y = new TProfile(*(other.prof_y));
}
ExtendedHistogram_2D::~ExtendedHistogram_2D(){
  delete histo;
  delete prof_x;
  delete prof_y;
}

void ExtendedHistogram_2D::swap(ExtendedHistogram_2D& other){
  std::swap(name, other.name);
  std::swap(title, other.title);
  std::swap(xbinning, other.xbinning);
  std::swap(ybinning, other.ybinning);
  std::swap(histo, other.histo);
  std::swap(prof_x, other.prof_x);
  std::swap(prof_y, other.prof_y);
}
ExtendedHistogram_2D& ExtendedHistogram_2D::operator=(const ExtendedHistogram_2D& other){
  ExtendedHistogram_2D tmp(other);
  this->swap(tmp);
  return *this;
}

void ExtendedHistogram_2D::setNameTitle(const TString name_, const TString title_){
  ExtendedHistogram::setNameTitle(name, title);
  if (histo && prof_x && prof_y){
    histo->SetName(name); histo->SetTitle(title);
    prof_x->SetName(Form("%s_prof_%s", name.Data(), xbinning.getLabel().Data())); prof_x->SetTitle(title);
    prof_y->SetName(Form("%s_prof_%s", name.Data(), ybinning.getLabel().Data())); prof_y->SetTitle(title);
  }
}
void ExtendedHistogram_2D::setBinning(const ExtendedBinning& binning, const int xyz, const TString label){
  if (xyz==0){
    xbinning = binning;
    if (label!="") xbinning.setLabel(label);
  }
  else{
    ybinning = binning;
    if (label!="") ybinning.setLabel(label);
  }
}
void ExtendedHistogram_2D::build(){
  if (xbinning.isValid() && ybinning.isValid()){
    const double* xbins = xbinning.getBinning();
    const double* ybins = ybinning.getBinning();
    const int nbinsx = xbinning.getNbins();
    const int nbinsy = ybinning.getNbins();
    histo = new TH2F(name, title, nbinsx, xbins, nbinsy, ybins); histo->GetXaxis()->SetTitle(xbinning.getLabel()); histo->GetYaxis()->SetTitle(ybinning.getLabel()); histo->Sumw2();
    prof_x = new TProfile(Form("%s_prof_%s", name.Data(), xbinning.getLabel().Data()), title, nbinsx, xbins); prof_x->GetXaxis()->SetTitle(xbinning.getLabel()); prof_x->Sumw2();
    prof_y = new TProfile(Form("%s_prof_%s", name.Data(), ybinning.getLabel().Data()), title, nbinsy, ybins); prof_y->GetXaxis()->SetTitle(ybinning.getLabel()); prof_y->Sumw2();
  }
}

void ExtendedHistogram_2D::fill(double x, double y, double wgt){
  if (histo && prof_x && prof_y){
    histo->Fill(x, y, wgt);
    prof_x->Fill(x, x, wgt);
    prof_y->Fill(y, y, wgt);
  }
}

TH2F* ExtendedHistogram_2D::getCumulantHistogram(TString newname) const{
  if (!histo) return nullptr;
  if (newname=="") newname=Form("Cumulant_%s", histo->GetName());
  TH2F* res = nullptr;
  HelperFunctions::getCumulantHistogram(histo, res, newname);
  return res;
}

ExtendedHistogram_2D ExtendedHistogram_2D::divideHistograms(ExtendedHistogram_2D const& h1, ExtendedHistogram_2D const& h2, bool useEffErr, TString newname){
  if (newname=="") newname=Form("h_%s_over_%s", h1.name.Data(), h2.name.Data());
  ExtendedHistogram_2D res(h2); res.setNameTitle(newname); res.histo->Reset("ICES");
  if (!h1.histo || !h1.prof_x || !h1.prof_y || !h2.histo || !h2.prof_x || !h2.prof_y) return res;
  HelperFunctions::divideHistograms(h1.histo, h2.histo, res.histo, useEffErr);
  return res;
}
