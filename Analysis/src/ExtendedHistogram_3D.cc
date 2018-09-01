#include "HelperFunctions.h"
#include "ExtendedHistogram_3D.h"
#include "MELAStreamHelpers.hh"


using namespace HelperFunctions;
using namespace MELAStreamHelpers;


ExtendedHistogram_3D::ExtendedHistogram_3D() : ExtendedHistogram(), histo(nullptr), prof_x(nullptr), prof_y(nullptr), prof_z(nullptr){}
ExtendedHistogram_3D::ExtendedHistogram_3D(const TString name_, const TString title_) : ExtendedHistogram(name_, title_), histo(nullptr), prof_x(nullptr), prof_y(nullptr), prof_z(nullptr){}
ExtendedHistogram_3D::ExtendedHistogram_3D(const TString name_, const TString title_, const ExtendedBinning& xbinning_, const ExtendedBinning& ybinning_, const ExtendedBinning& zbinning_) :
  ExtendedHistogram(name_, title_),
  xbinning(xbinning_),
  ybinning(ybinning_),
  zbinning(zbinning_),
  histo(nullptr), prof_x(nullptr), prof_y(nullptr), prof_z(nullptr)
{
  build();
}
ExtendedHistogram_3D::ExtendedHistogram_3D(ExtendedHistogram_3D const& other) : ExtendedHistogram(other), xbinning(other.xbinning), ybinning(other.ybinning), histo(nullptr), prof_x(nullptr), prof_y(nullptr), prof_z(nullptr){
  if (other.histo) histo = new TH3F(*(other.histo));
  if (other.prof_x) prof_x = new TProfile(*(other.prof_x));
  if (other.prof_y) prof_y = new TProfile(*(other.prof_y));
  if (other.prof_z) prof_z = new TProfile(*(other.prof_z));
}
ExtendedHistogram_3D::~ExtendedHistogram_3D(){
  delete histo;
  delete prof_x;
  delete prof_y;
  delete prof_z;
}

void ExtendedHistogram_3D::swap(ExtendedHistogram_3D& other){
  std::swap(name, other.name);
  std::swap(title, other.title);
  std::swap(xbinning, other.xbinning);
  std::swap(ybinning, other.ybinning);
  std::swap(zbinning, other.zbinning);
  std::swap(histo, other.histo);
  std::swap(prof_x, other.prof_x);
  std::swap(prof_y, other.prof_y);
  std::swap(prof_z, other.prof_z);
}
ExtendedHistogram_3D& ExtendedHistogram_3D::operator=(const ExtendedHistogram_3D& other){
  ExtendedHistogram_3D tmp(other);
  this->swap(tmp);
  return *this;
}

void ExtendedHistogram_3D::setNameTitle(const TString name_, const TString title_){
  ExtendedHistogram::setNameTitle(name_, title_);
  if (histo){ histo->SetName(name); histo->SetTitle(title); }
  if (prof_x){ prof_x->SetName(Form("%s_prof_%s", name.Data(), xbinning.getLabel().Data())); prof_x->SetTitle(title); }
  if (prof_y){ prof_y->SetName(Form("%s_prof_%s", name.Data(), ybinning.getLabel().Data())); prof_y->SetTitle(title); }
  if (prof_z){ prof_z->SetName(Form("%s_prof_%s", name.Data(), zbinning.getLabel().Data())); prof_z->SetTitle(title); }
}
void ExtendedHistogram_3D::setBinning(const ExtendedBinning& binning, const int xyz, const TString label){
  if (xyz==0){
    xbinning = binning;
    if (label!="") xbinning.setLabel(label);
  }
  else if (xyz==1){
    ybinning = binning;
    if (label!="") ybinning.setLabel(label);
  }
  else{
    zbinning = binning;
    if (label!="") zbinning.setLabel(label);
  }
}
ExtendedBinning const& ExtendedHistogram_3D::getBinning(const int xyz){
  if (xyz==0) return xbinning;
  else if (xyz==1) return ybinning;
  else return zbinning;
}
void ExtendedHistogram_3D::build(){
  reset();
  if (xbinning.isValid() && ybinning.isValid() && zbinning.isValid()){
    const double* xbins = xbinning.getBinning();
    const double* ybins = ybinning.getBinning();
    const double* zbins = zbinning.getBinning();
    const int nbinsx = xbinning.getNbins();
    const int nbinsy = ybinning.getNbins();
    const int nbinsz = zbinning.getNbins();
    histo = new TH3F(name, title, nbinsx, xbins, nbinsy, ybins, nbinsz, zbins); histo->GetXaxis()->SetTitle(xbinning.getLabel()); histo->GetYaxis()->SetTitle(ybinning.getLabel()); histo->GetZaxis()->SetTitle(zbinning.getLabel()); histo->Sumw2();
    prof_x = new TProfile(Form("%s_prof_%s", name.Data(), xbinning.getLabel().Data()), title, nbinsx, xbins); prof_x->GetXaxis()->SetTitle(xbinning.getLabel()); prof_x->Sumw2();
    prof_y = new TProfile(Form("%s_prof_%s", name.Data(), ybinning.getLabel().Data()), title, nbinsy, ybins); prof_y->GetXaxis()->SetTitle(ybinning.getLabel()); prof_y->Sumw2();
    prof_z = new TProfile(Form("%s_prof_%s", name.Data(), zbinning.getLabel().Data()), title, nbinsz, zbins); prof_z->GetXaxis()->SetTitle(zbinning.getLabel()); prof_z->Sumw2();
  }
}
void ExtendedHistogram_3D::reset(){
  delete histo;
  delete prof_x;
  delete prof_y;
  delete prof_z;
}

void ExtendedHistogram_3D::fill(double x, double y, double z, double wgt){
  if (histo) histo->Fill(x, y, z, wgt);
  if (prof_x) prof_x->Fill(x, x, fabs(wgt));
  if (prof_y) prof_y->Fill(y, y, fabs(wgt));
  if (prof_z) prof_z->Fill(z, z, fabs(wgt));
}

void ExtendedHistogram_3D::rebin(ExtendedBinning const* binningX, ExtendedBinning const* binningY, ExtendedBinning const* binningZ, signed char condDim){
  if (!binningX) binningX=&xbinning;
  if (!binningY) binningY=&ybinning;
  if (!binningZ) binningZ=&zbinning;
  bool condX=(condDim==0);
  bool condY=(condDim==1);
  bool condZ=(condDim==2);
  if (binningX->isValid() && binningY->isValid() && binningZ->isValid()){
    std::vector<std::pair<TProfile const*, unsigned int>> condProfs;
    if (condX) condProfs.emplace_back(prof_x, 0);
    if (condY) condProfs.emplace_back(prof_y, 1);
    if (condZ) condProfs.emplace_back(prof_z, 2);
    if (histo) rebinHistogram(histo, *binningX, *binningY, *binningZ, (condX || condY || condZ ? &condProfs : nullptr));
    if (prof_x) rebinProfile(prof_x, *binningX);
    if (prof_y) rebinProfile(prof_y, *binningY);
    if (prof_z) rebinProfile(prof_z, *binningZ);
  }
}

TH3F* ExtendedHistogram_3D::getCumulantHistogram(TString newname) const{
  if (!histo) return nullptr;
  if (newname=="") newname=Form("Cumulant_%s", histo->GetName());
  TH3F* res = nullptr;
  HelperFunctions::getCumulantHistogram(histo, res, newname);
  return res;
}

ExtendedHistogram_3D ExtendedHistogram_3D::divideHistograms(ExtendedHistogram_3D const& h1, ExtendedHistogram_3D const& h2, bool useEffErr, TString newname){
  if (newname=="") newname=Form("h_%s_over_%s", h1.name.Data(), h2.name.Data());
  ExtendedHistogram_3D res(h2); res.setNameTitle(newname); res.histo->Reset("ICES");
  if (!h1.histo || !h2.histo) return res;
  HelperFunctions::divideHistograms(h1.histo, h2.histo, res.histo, useEffErr);
  if (!useEffErr){
    if (h1.prof_x && h2.prof_x) combineHistogramsByWeightedAverage(h1.prof_x, h2.prof_x, res.prof_x);
    if (h1.prof_y && h2.prof_y) combineHistogramsByWeightedAverage(h1.prof_y, h2.prof_y, res.prof_y);
    if (h1.prof_z && h2.prof_z) combineHistogramsByWeightedAverage(h1.prof_z, h2.prof_z, res.prof_z);
  }
  return res;
}

void ExtendedHistogram_3D::averageHistograms(ExtendedHistogram_3D& hTarget, ExtendedHistogram_3D const& h2, bool useNeff){
  if (!hTarget.histo || !h2.histo) return;
  combineHistogramsByWeightedAverage(hTarget.histo, h2.histo, hTarget.histo, useNeff);
  if (hTarget.prof_x && h2.prof_x) combineHistogramsByWeightedAverage(hTarget.prof_x, h2.prof_x, hTarget.prof_x, useNeff);
  if (hTarget.prof_y && h2.prof_y) combineHistogramsByWeightedAverage(hTarget.prof_y, h2.prof_y, hTarget.prof_y, useNeff);
  if (hTarget.prof_z && h2.prof_z) combineHistogramsByWeightedAverage(hTarget.prof_z, h2.prof_z, hTarget.prof_z, useNeff);
}
void ExtendedHistogram_3D::averageHistograms(ExtendedHistogram_3D& hTarget, std::vector<ExtendedHistogram_3D const*> const& hList, bool useNeff){
  if (!hTarget.histo) return;

  std::vector<TH3F const*> hhlist;
  hhlist.push_back(hTarget.histo);

  std::vector<TProfile const*> hprof_xlist;
  hprof_xlist.push_back(hTarget.prof_x);
  std::vector<TProfile const*> hprof_ylist;
  hprof_ylist.push_back(hTarget.prof_y);
  std::vector<TProfile const*> hprof_zlist;
  hprof_zlist.push_back(hTarget.prof_z);

  for (ExtendedHistogram_3D const* const& h:hList){
    if (h){
      if (h->histo) hhlist.push_back(h->histo);
      if (h->prof_x) hprof_xlist.push_back(h->prof_x);
      if (h->prof_y) hprof_ylist.push_back(h->prof_y);
      if (h->prof_z) hprof_zlist.push_back(h->prof_z);
    }
  }

  combineHistogramListByWeightedAverage(hhlist, hTarget.histo, useNeff);
  if (hTarget.prof_x) combineHistogramListByWeightedAverage(hprof_xlist, hTarget.prof_x, useNeff);
  if (hTarget.prof_y) combineHistogramListByWeightedAverage(hprof_ylist, hTarget.prof_y, useNeff);
  if (hTarget.prof_z) combineHistogramListByWeightedAverage(hprof_zlist, hTarget.prof_z, useNeff);
}

void ExtendedHistogram_3D::constructFromTree(TTree* tree, float& xvar, float& yvar, float& zvar, float& weight, bool* flag, ExtendedBinning const* binningX, ExtendedBinning const* binningY, ExtendedBinning const* binningZ){
  if (!tree) return;
  if (binningX) setBinning(*binningX, 0, binningX->getLabel());
  if (binningY) setBinning(*binningY, 1, binningY->getLabel());
  if (binningZ) setBinning(*binningZ, 2, binningZ->getLabel());
  build();
  double xlow=xbinning.getMin() - xbinning.getBinWidth(0);
  double xhigh=xbinning.getMax() + xbinning.getBinWidth(xbinning.getNbins()-1);
  double ylow=ybinning.getMin() - ybinning.getBinWidth(0);
  double yhigh=ybinning.getMax() + ybinning.getBinWidth(ybinning.getNbins()-1);
  double zlow=zbinning.getMin() - zbinning.getBinWidth(0);
  double zhigh=zbinning.getMax() + zbinning.getBinWidth(zbinning.getNbins()-1);
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    if (xvar<xlow || xvar>=xhigh) continue;
    if (yvar<ylow || yvar>=yhigh) continue;
    if (zvar<zlow || zvar>=zhigh) continue;
    if (!flag || (flag && *flag)) fill(xvar, yvar, zvar, weight);
  }
}
