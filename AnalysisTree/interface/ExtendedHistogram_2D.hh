#include "HelperFunctions.h"


template<typename Hist_t> ExtendedHistogram_2D<Hist_t>::ExtendedHistogram_2D() : ExtendedHistogram(), histo(nullptr), prof_x(nullptr), prof_y(nullptr){}
template<typename Hist_t> ExtendedHistogram_2D<Hist_t>::ExtendedHistogram_2D(const TString name_, const TString title_) : ExtendedHistogram(name_, title_), histo(nullptr), prof_x(nullptr), prof_y(nullptr){}
template<typename Hist_t> ExtendedHistogram_2D<Hist_t>::ExtendedHistogram_2D(const TString name_, const TString title_, const ExtendedBinning& xbinning_, const ExtendedBinning& ybinning_) :
  ExtendedHistogram(name_, title_),
  xbinning(xbinning_),
  ybinning(ybinning_),
  histo(nullptr), prof_x(nullptr), prof_y(nullptr)
{
  build();
}
template<typename Hist_t> ExtendedHistogram_2D<Hist_t>::ExtendedHistogram_2D(ExtendedHistogram_2D<Hist_t> const& other) : ExtendedHistogram(other), xbinning(other.xbinning), ybinning(other.ybinning), histo(nullptr), prof_x(nullptr), prof_y(nullptr){
  if (other.histo) histo = new Hist_t(*(other.histo));
  if (other.prof_x) prof_x = new TProfile(*(other.prof_x));
  if (other.prof_y) prof_y = new TProfile(*(other.prof_y));
}
template<typename Hist_t> ExtendedHistogram_2D<Hist_t>::~ExtendedHistogram_2D(){ reset(); }

template<typename Hist_t> void ExtendedHistogram_2D<Hist_t>::swap(ExtendedHistogram_2D<Hist_t>& other){
  std::swap(name, other.name);
  std::swap(title, other.title);
  std::swap(xbinning, other.xbinning);
  std::swap(ybinning, other.ybinning);
  std::swap(histo, other.histo);
  std::swap(prof_x, other.prof_x);
  std::swap(prof_y, other.prof_y);
}
template<typename Hist_t> ExtendedHistogram_2D<Hist_t>& ExtendedHistogram_2D<Hist_t>::operator=(const ExtendedHistogram_2D<Hist_t>& other){
  ExtendedHistogram_2D<Hist_t> tmp(other);
  this->swap(tmp);
  return *this;
}

template<typename Hist_t> void ExtendedHistogram_2D<Hist_t>::setNameTitle(const TString name_, const TString title_){
  ExtendedHistogram::setNameTitle(name_, title_);
  if (histo){ histo->SetName(name); histo->SetTitle(title); }
  if (prof_x){ prof_x->SetName(Form("%s_prof_%s", name.Data(), xbinning.getName().Data())); prof_x->SetTitle(title); }
  if (prof_y){ prof_y->SetName(Form("%s_prof_%s", name.Data(), ybinning.getName().Data())); prof_y->SetTitle(title); }
}
template<typename Hist_t> void ExtendedHistogram_2D<Hist_t>::setBinning(const ExtendedBinning& binning, const int xyz, const TString label){
  if (xyz==0){
    xbinning = binning;
    if (label!="") xbinning.setLabel(label);
  }
  else{
    ybinning = binning;
    if (label!="") ybinning.setLabel(label);
  }
}
template<typename Hist_t> ExtendedBinning const& ExtendedHistogram_2D<Hist_t>::getBinning(const int xyz){
  if (xyz==0) return xbinning;
  else return ybinning;
}
template<typename Hist_t> void ExtendedHistogram_2D<Hist_t>::build(){
  reset();
  if (xbinning.isValid() && ybinning.isValid()){
    const double* xbins = xbinning.getBinning();
    const double* ybins = ybinning.getBinning();
    const int nbinsx = xbinning.getNbins();
    const int nbinsy = ybinning.getNbins();
    histo = new Hist_t(name, title, nbinsx, xbins, nbinsy, ybins); histo->GetXaxis()->SetTitle(xbinning.getLabel()); histo->GetYaxis()->SetTitle(ybinning.getLabel()); histo->Sumw2();
    prof_x = new TProfile(Form("%s_prof_%s", name.Data(), xbinning.getName().Data()), title, nbinsx, xbins); prof_x->GetXaxis()->SetTitle(xbinning.getLabel()); prof_x->Sumw2();
    prof_y = new TProfile(Form("%s_prof_%s", name.Data(), ybinning.getName().Data()), title, nbinsy, ybins); prof_y->GetXaxis()->SetTitle(ybinning.getLabel()); prof_y->Sumw2();
  }
}
template<typename Hist_t> void ExtendedHistogram_2D<Hist_t>::reset(){
  HelperFunctions::resetPointer(histo);
  HelperFunctions::resetPointer(prof_x);
  HelperFunctions::resetPointer(prof_y);
}
template<typename Hist_t> void ExtendedHistogram_2D<Hist_t>::resetProfiles(){
  HelperFunctions::resetPointer(prof_x);
  HelperFunctions::resetPointer(prof_y);
}

template<typename Hist_t> void ExtendedHistogram_2D<Hist_t>::fill(double x, double y, double wgt){
  if (histo) histo->Fill(x, y, wgt);
  if (prof_x) prof_x->Fill(x, x, fabs(wgt));
  if (prof_y) prof_y->Fill(y, y, fabs(wgt));
}

template<typename Hist_t> void ExtendedHistogram_2D<Hist_t>::rebin(ExtendedBinning const* binningX, ExtendedBinning const* binningY, signed char condDim){
  if (!binningX) binningX=&xbinning;
  if (!binningY) binningY=&ybinning;
  bool condX=(condDim==0);
  bool condY=(condDim==1);
  if (binningX->isValid() && binningY->isValid()){
    std::vector<std::pair<TProfile const*, unsigned int>> condProfs;
    if (condX) condProfs.emplace_back(prof_x, 0);
    if (condY) condProfs.emplace_back(prof_y, 1);
    if (histo) HelperFunctions::rebinHistogram(histo, *binningX, *binningY, (condX || condY ? &condProfs : nullptr));
    if (prof_x) HelperFunctions::rebinProfile(prof_x, *binningX);
    if (prof_y) HelperFunctions::rebinProfile(prof_y, *binningY);
  }
}


template<typename Hist_t> Hist_t* ExtendedHistogram_2D<Hist_t>::getCumulantHistogram(TString newname) const{
  if (!histo) return nullptr;
  if (newname=="") newname=Form("Cumulant_%s", histo->GetName());
  Hist_t* res = nullptr;
  HelperFunctions::getCumulantHistogram(histo, res, newname);
  return res;
}

template<typename Hist_t> ExtendedHistogram_2D<Hist_t> ExtendedHistogram_2D<Hist_t>::divideHistograms(ExtendedHistogram_2D<Hist_t> const& h1, ExtendedHistogram_2D<Hist_t> const& h2, bool useEffErr, TString newname){
  if (newname=="") newname=Form("h_%s_over_%s", h1.name.Data(), h2.name.Data());
  ExtendedHistogram_2D<Hist_t> res(h2); res.setNameTitle(newname); res.histo->Reset("ICES");
  if (!h1.histo || !h2.histo) return res;
  HelperFunctions::divideHistograms(h1.histo, h2.histo, res.histo, useEffErr);
  if (!useEffErr){
    if (h1.prof_x && h2.prof_x) HelperFunctions::combineHistogramsByWeightedAverage(h1.prof_x, h2.prof_x, res.prof_x);
    if (h1.prof_y && h2.prof_y) HelperFunctions::combineHistogramsByWeightedAverage(h1.prof_y, h2.prof_y, res.prof_y);
  }
  return res;
}

template<typename Hist_t> void ExtendedHistogram_2D<Hist_t>::averageHistograms(ExtendedHistogram_2D<Hist_t>& hTarget, ExtendedHistogram_2D<Hist_t> const& h2, bool useNeff){
  if (!hTarget.histo || !h2.histo) return;
  HelperFunctions::combineHistogramsByWeightedAverage(hTarget.histo, h2.histo, hTarget.histo, useNeff);
  if (hTarget.prof_x && h2.prof_x) HelperFunctions::combineHistogramsByWeightedAverage(hTarget.prof_x, h2.prof_x, hTarget.prof_x, useNeff);
  if (hTarget.prof_y && h2.prof_y) HelperFunctions::combineHistogramsByWeightedAverage(hTarget.prof_y, h2.prof_y, hTarget.prof_y, useNeff);
}
template<typename Hist_t> void ExtendedHistogram_2D<Hist_t>::averageHistograms(ExtendedHistogram_2D<Hist_t>& hTarget, std::vector<ExtendedHistogram_2D<Hist_t> const*> const& hList, bool useNeff){
  if (!hTarget.histo) return;

  std::vector<Hist_t const*> hhlist;
  hhlist.push_back(hTarget.histo);

  std::vector<TProfile const*> hprof_xlist;
  hprof_xlist.push_back(hTarget.prof_x);
  std::vector<TProfile const*> hprof_ylist;
  hprof_ylist.push_back(hTarget.prof_y);

  for (auto const& h:hList){
    if (h){
      if (h->histo) hhlist.push_back(h->histo);
      if (h->prof_x) hprof_xlist.push_back(h->prof_x);
      if (h->prof_y) hprof_ylist.push_back(h->prof_y);
    }
  }

  HelperFunctions::combineHistogramListByWeightedAverage(hhlist, hTarget.histo, useNeff);
  if (hTarget.prof_x) HelperFunctions::combineHistogramListByWeightedAverage(hprof_xlist, hTarget.prof_x, useNeff);
  if (hTarget.prof_y) HelperFunctions::combineHistogramListByWeightedAverage(hprof_ylist, hTarget.prof_y, useNeff);
}

template<typename Hist_t> void ExtendedHistogram_2D<Hist_t>::constructFromTree(TTree* tree, float& xvar, float& yvar, float& weight, bool* flag, ExtendedBinning const* binningX, ExtendedBinning const* binningY){
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
