#include "HelperFunctions.h"
#include "ExtendedHistogram_2D.h"
#include "MELAStreamHelpers.hh"


using namespace HelperFunctions;
using namespace MELAStreamHelpers;


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
  ExtendedHistogram::setNameTitle(name_, title_);
  if (histo){ histo->SetName(name); histo->SetTitle(title); }
  if (prof_x){ prof_x->SetName(Form("%s_prof_%s", name.Data(), xbinning.getLabel().Data())); prof_x->SetTitle(title); }
  if (prof_y){ prof_y->SetName(Form("%s_prof_%s", name.Data(), ybinning.getLabel().Data())); prof_y->SetTitle(title); }
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
ExtendedBinning const& ExtendedHistogram_2D::getBinning(const int xyz){
  if (xyz==0) return xbinning;
  else return ybinning;
}
void ExtendedHistogram_2D::build(){
  reset();
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
void ExtendedHistogram_2D::reset(){
  delete histo;
  delete prof_x;
  delete prof_y;
}

void ExtendedHistogram_2D::fill(double x, double y, double wgt){
  if (histo) histo->Fill(x, y, wgt);
  if (prof_x) prof_x->Fill(x, x, wgt);
  if (prof_y) prof_y->Fill(y, y, wgt);
}

void ExtendedHistogram_2D::rebin(ExtendedBinning const* binningX, ExtendedBinning const* binningY, signed char condDim){
  if (!binningX) binningX=&xbinning;
  if (!binningY) binningY=&ybinning;
  bool condX=(condDim==0);
  bool condY=(condDim==1);
  if (binningX->isValid() && binningY->isValid()){
    std::vector<std::pair<TProfile const*, unsigned int>> condProfs;
    if (condX) condProfs.emplace_back(prof_x, 0);
    if (condY) condProfs.emplace_back(prof_y, 1);
    if (histo) rebinHistogram(histo, *binningX, *binningY, (condX || condY ? &condProfs : nullptr));
    if (prof_x) rebinProfile(prof_x, *binningX);
    if (prof_y) rebinProfile(prof_y, *binningY);
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
  if (!h1.histo || !h2.histo) return res;
  HelperFunctions::divideHistograms(h1.histo, h2.histo, res.histo, useEffErr);
  if (!useEffErr){
    if (h1.prof_x && h2.prof_x) combineHistogramsByWeightedAverage(h1.prof_x, h2.prof_x, res.prof_x);
    if (h1.prof_y && h2.prof_y) combineHistogramsByWeightedAverage(h1.prof_y, h2.prof_y, res.prof_y);
  }
  return res;
}

void ExtendedHistogram_2D::averageHistograms(ExtendedHistogram_2D& hTarget, ExtendedHistogram_2D const& h2, bool useNeff){
  if (!hTarget.histo || !h2.histo) return;
  combineHistogramsByWeightedAverage(hTarget.histo, h2.histo, hTarget.histo, useNeff);
  if (hTarget.prof_x && h2.prof_x) combineHistogramsByWeightedAverage(hTarget.prof_x, h2.prof_x, hTarget.prof_x, useNeff);
  if (hTarget.prof_y && h2.prof_y) combineHistogramsByWeightedAverage(hTarget.prof_y, h2.prof_y, hTarget.prof_y, useNeff);
}
void ExtendedHistogram_2D::averageHistograms(ExtendedHistogram_2D& hTarget, std::vector<ExtendedHistogram_2D const*> const& hList, bool useNeff){
  if (!hTarget.histo) return;

  std::vector<TH2F const*> hhlist;
  hhlist.push_back(hTarget.histo);

  std::vector<TProfile const*> hprof_xlist;
  hprof_xlist.push_back(hTarget.prof_x);
  std::vector<TProfile const*> hprof_ylist;
  hprof_ylist.push_back(hTarget.prof_y);

  for (ExtendedHistogram_2D const* const& h:hList){
    if (h){
      if (h->histo) hhlist.push_back(h->histo);
      if (h->prof_x) hprof_xlist.push_back(h->prof_x);
      if (h->prof_y) hprof_ylist.push_back(h->prof_y);
    }
  }

  combineHistogramListByWeightedAverage(hhlist, hTarget.histo, useNeff);
  if (hTarget.prof_x) combineHistogramListByWeightedAverage(hprof_xlist, hTarget.prof_x, useNeff);
  if (hTarget.prof_y) combineHistogramListByWeightedAverage(hprof_ylist, hTarget.prof_y, useNeff);
}

void ExtendedHistogram_2D::constructFromTree(TTree* tree, float& xvar, float& yvar, float& weight, bool* flag, ExtendedBinning const* binningX, ExtendedBinning const* binningY){
  if (!tree) return;
  if (binningX) setBinning(*binningX, 0, binningX->getLabel());
  if (binningY) setBinning(*binningY, 1, binningY->getLabel());
  build();
  double xlow=xbinning.getMin() - xbinning.getBinWidth(0);
  double xhigh=xbinning.getMax() + xbinning.getBinWidth(xbinning.getNbins()-1);
  double ylow=ybinning.getMin() - ybinning.getBinWidth(0);
  double yhigh=ybinning.getMax() + ybinning.getBinWidth(ybinning.getNbins()-1);
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    if (xvar<xlow || xvar>=xhigh) continue;
    if (yvar<ylow || yvar>=yhigh) continue;
    if (!flag || (flag && *flag)) fill(xvar, yvar, weight);
  }
}
