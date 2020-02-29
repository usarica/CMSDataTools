#include "HelperFunctions.h"
#include "ExtendedHistogram_1D.h"
#include "MELAStreamHelpers.hh"


using namespace HelperFunctions;
using namespace MELAStreamHelpers;


ExtendedHistogram_1D::ExtendedHistogram_1D() : ExtendedHistogram(), histo(nullptr), prof_x(nullptr){}
ExtendedHistogram_1D::ExtendedHistogram_1D(const TString name_, const TString title_) : ExtendedHistogram(name_, title_), histo(nullptr), prof_x(nullptr){}
ExtendedHistogram_1D::ExtendedHistogram_1D(const TString name_, const TString title_, const ExtendedBinning& xbinning_) :
  ExtendedHistogram(name_, title_),
  xbinning(xbinning_),
  histo(nullptr), prof_x(nullptr)
{
  build();
}
ExtendedHistogram_1D::ExtendedHistogram_1D(ExtendedHistogram_1D const& other) : ExtendedHistogram(other), xbinning(other.xbinning), histo(nullptr), prof_x(nullptr){
  if (other.histo) histo = new TH1F(*(other.histo));
  if (other.prof_x) prof_x = new TProfile(*(other.prof_x));
}
ExtendedHistogram_1D::~ExtendedHistogram_1D(){ reset(); }

void ExtendedHistogram_1D::swap(ExtendedHistogram_1D& other){
  std::swap(name, other.name);
  std::swap(title, other.title);
  std::swap(xbinning, other.xbinning);
  std::swap(histo, other.histo);
  std::swap(prof_x, other.prof_x);
}
ExtendedHistogram_1D& ExtendedHistogram_1D::operator=(const ExtendedHistogram_1D& other){
  ExtendedHistogram_1D tmp(other);
  this->swap(tmp);
  return *this;
}

void ExtendedHistogram_1D::setNameTitle(const TString name_, const TString title_){
  ExtendedHistogram::setNameTitle(name_, title_);
  if (histo){ histo->SetName(name); histo->SetTitle(title); }
  if (prof_x){ prof_x->SetName(Form("%s_prof_%s", name.Data(), xbinning.getLabel().Data())); prof_x->SetTitle(title); }
}
void ExtendedHistogram_1D::setBinning(const ExtendedBinning& binning, const int /*xyz*/, const TString label){
  xbinning = binning;
  if (label!="") xbinning.setLabel(label);
}
ExtendedBinning const& ExtendedHistogram_1D::getBinning(const int /*xyz*/){
  return xbinning;
}
void ExtendedHistogram_1D::build(){
  reset();
  if (xbinning.isValid()){
    const double* xbins = xbinning.getBinning();
    const int nbins = xbinning.getNbins();
    histo = new TH1F(name, title, nbins, xbins); histo->GetXaxis()->SetTitle(xbinning.getLabel()); histo->Sumw2();
    prof_x = new TProfile(Form("%s_prof_%s", name.Data(), xbinning.getLabel().Data()), title, nbins, xbins); prof_x->GetXaxis()->SetTitle(xbinning.getLabel()); prof_x->Sumw2();
  }
}
void ExtendedHistogram_1D::reset(){
  resetPointer(histo);
  resetPointer(prof_x);
}
void ExtendedHistogram_1D::resetProfiles(){
  resetPointer(prof_x);
}

void ExtendedHistogram_1D::fill(double x, double wgt){
  if (histo) histo->Fill(x, wgt);
  if (prof_x) prof_x->Fill(x, x, fabs(wgt));
}

void ExtendedHistogram_1D::rebin(ExtendedBinning const* binningX){
  if (binningX && binningX->isValid()){
    if (histo) rebinHistogram(histo, *binningX);
    if (prof_x) rebinProfile(prof_x, *binningX);
  }
}

TGraphErrors* ExtendedHistogram_1D::getGraph(TString newname) const{
  if (!histo || !prof_x) return nullptr;
  if (newname=="") newname=Form("gr_%s_vs_%s", histo->GetName(), prof_x->GetName());
  return HelperFunctions::makeGraphFromTH1(prof_x, histo, newname);
}

TH1F* ExtendedHistogram_1D::getCumulantHistogram(TString newname) const{
  if (!histo) return nullptr;
  if (newname=="") newname=Form("Cumulant_%s", histo->GetName());
  TH1F* res = nullptr;
  HelperFunctions::getCumulantHistogram(histo, res, newname);
  return res;
}

ExtendedHistogram_1D ExtendedHistogram_1D::divideHistograms(ExtendedHistogram_1D const& h1, ExtendedHistogram_1D const& h2, bool useEffErr, TString newname){
  if (newname=="") newname=Form("h_%s_over_%s", h1.name.Data(), h2.name.Data());
  ExtendedHistogram_1D res(h2); res.setNameTitle(newname); res.histo->Reset("ICES");
  if (!h1.histo || !h2.histo) return res;
  HelperFunctions::divideHistograms(h1.histo, h2.histo, res.histo, useEffErr);
  if (!useEffErr){
    if (h1.prof_x && h2.prof_x) combineHistogramsByWeightedAverage(h1.prof_x, h2.prof_x, res.prof_x);
  }
  return res;
}

void ExtendedHistogram_1D::averageHistograms(ExtendedHistogram_1D& hTarget, ExtendedHistogram_1D const& h2, bool useNeff){
  if (!hTarget.histo || !h2.histo) return;
  combineHistogramsByWeightedAverage(hTarget.histo, h2.histo, hTarget.histo, useNeff);
  if (hTarget.prof_x && h2.prof_x) combineHistogramsByWeightedAverage(hTarget.prof_x, h2.prof_x, hTarget.prof_x, useNeff);
}

void ExtendedHistogram_1D::averageHistograms(ExtendedHistogram_1D& hTarget, std::vector<ExtendedHistogram_1D const*> const& hList, bool useNeff){
  if (!hTarget.histo) return;

  std::vector<TH1F const*> hhlist;
  hhlist.push_back(hTarget.histo);

  std::vector<TProfile const*> hprof_xlist;
  hprof_xlist.push_back(hTarget.prof_x);

  for (ExtendedHistogram_1D const* const& h:hList){
    if (h){
      if (h->histo) hhlist.push_back(h->histo);
      if (h->prof_x) hprof_xlist.push_back(h->prof_x);
    }
  }

  combineHistogramListByWeightedAverage(hhlist, hTarget.histo, useNeff);
  if (hTarget.prof_x) combineHistogramListByWeightedAverage(hprof_xlist, hTarget.prof_x, useNeff);
}

void ExtendedHistogram_1D::constructFromTree(TTree* tree, float& xvar, float& weight, bool* flag, ExtendedBinning const* binningX){
  if (!tree) return;
  if (binningX) setBinning(*binningX, 0, binningX->getLabel());
  build();
  double xlow=xbinning.getMin() - xbinning.getBinWidth(0);
  double xhigh=xbinning.getMax() + xbinning.getBinWidth(xbinning.getNbins()-1);
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    if (xvar<xlow || xvar>=xhigh) continue;
    if (!flag || (flag && *flag)) fill(xvar, weight);
  }
}
