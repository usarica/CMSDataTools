#include "HelperFunctions.h"
#include "ExtendedHistogram_1D.h"


using namespace HelperFunctions;


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
ExtendedHistogram_1D::~ExtendedHistogram_1D(){
  delete histo;
  delete prof_x;
}

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
  ExtendedHistogram::setNameTitle(name, title);
  if (histo){ histo->SetName(name); histo->SetTitle(title); }
  if (prof_x){ prof_x->SetName(Form("%s_prof_%s", name.Data(), xbinning.getLabel().Data())); prof_x->SetTitle(title); }
}
void ExtendedHistogram_1D::setBinning(const ExtendedBinning& binning, const int xyz, const TString label){
  xbinning = binning;
  if (label!="") xbinning.setLabel(label);
}
void ExtendedHistogram_1D::build(){
  if (xbinning.isValid()){
    const double* xbins = xbinning.getBinning();
    const int nbins = xbinning.getNbins();
    histo = new TH1F(name, title, nbins, xbins); histo->GetXaxis()->SetTitle(xbinning.getLabel()); histo->Sumw2();
    prof_x = new TProfile(Form("%s_prof_%s", name.Data(), xbinning.getLabel().Data()), title, nbins, xbins); prof_x->GetXaxis()->SetTitle(xbinning.getLabel()); prof_x->Sumw2();
  }
}

void ExtendedHistogram_1D::fill(double x, double wgt){
  if (histo) histo->Fill(x, wgt);
  if (prof_x) prof_x->Fill(x, x, wgt);
}

void ExtendedHistogram_1D::rebin(ExtendedBinning const& binningX){
  if (binningX.isValid()){
    if (histo) rebinHistogram(histo, binningX);
    if (prof_x) rebinProfile(prof_x, binningX);
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
    if (h1.prof_x && h2.prof_x){
      for (unsigned int bin=0; bin<=h2.xbinning.getNbins()+1; bin++){
        double val[2]={ h1.prof_x->GetBinContent(bin), h2.prof_x->GetBinContent(bin) };
        double errsq[2]={ pow(h1.prof_x->GetBinError(bin), 2), pow(h2.prof_x->GetBinError(bin), 2) };
        double valnew=0;
        double wnew=0;
        for (unsigned int ii=0; ii<2; ii++){ if (errsq[ii]!=0.){ valnew += val[ii]/errsq[ii]; wnew += 1./errsq[ii]; } }
        if (wnew==0.) valnew = 0;
        else valnew /= wnew;
        res.prof_x->SetBinContent(bin, valnew);
        res.prof_x->SetBinError(bin, sqrt(1./wnew));
      }
    }
  }
  return res;
}
