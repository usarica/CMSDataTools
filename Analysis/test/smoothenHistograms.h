#ifndef SMOOTHENHISTOGRAMS_H
#define SMOOTHENHISTOGRAMS_H

#include "common_includes.h"


ExtendedBinning getFinestBinning(ExtendedBinning const& finalBinning, unsigned int& finalNbins){
  finalNbins=finalBinning.getNbins();
  float varmin=finalBinning.getBinLowEdge(0);
  float varmax=finalBinning.getBinLowEdge(finalNbins);

  TString varname=finalBinning.getLabel(); TString varnamelower=varname; varnamelower.ToLower();
  unsigned int nBinsIntermediate;
  if (varnamelower.Contains("mass")) nBinsIntermediate = 1000*finalNbins;
  else nBinsIntermediate = 100*finalNbins;
  return ExtendedBinning(nBinsIntermediate, varmin, varmax, varname);
}

ExtendedBinning getIntermediateBinning(TH1F const* fineHisto, int effectiveN){
  ExtendedBinning res;
  double integralerror=0;
  double integral=fineHisto->IntegralAndError(1, fineHisto->GetNbinsX(), integralerror); // Do not include underflow/overflow bins; they are not relevant for the binning of the cumulant
  double effectiveNtotal=(integralerror>0. ? pow(integral/integralerror, 2) : fineHisto->GetEntries());
  double effNThreshold;
  if (effectiveN>0 && effectiveNtotal>2.*double(effectiveN)) effNThreshold = effectiveN;
  else effNThreshold=effectiveNtotal/20.;
  MELAout << "Effective total N = " << effectiveNtotal << ", requested effective N = " << effectiveN << ", threshold = " << effNThreshold << endl;

  res.addBinBoundary(fineHisto->GetXaxis()->GetBinLowEdge(1));
  vector<float> effNs;
  float sum=0, sumerrsq=0;
  for (int bin=1; bin<=fineHisto->GetNbinsX(); bin++){
    sum += fineHisto->GetBinContent(bin);
    sumerrsq += pow(fineHisto->GetBinError(bin), 2);
    double effN=(sumerrsq>0. ? pow(sum, 2)/sumerrsq : sum/integral*effectiveNtotal);
    if (effN>effNThreshold || bin==fineHisto->GetNbinsX()){
      effNs.push_back(effN);
      res.addBinBoundary(fineHisto->GetXaxis()->GetBinLowEdge(bin+1));
      sum=0;
      sumerrsq=0;
    }
  }
  // These lines guarantee sum>sumThreshold in every bin
  if (effNs.size()>1 && effNs.at(effNs.size()-1)<effNThreshold){
    effNs.at(effNs.size()-1) += effNs.at(effNs.size()-2);
    effNs.erase(effNs.begin()+effNs.size()-2);
    res.removeBinLowEdge(res.getNbins()-1);
  }
  // Loop over histogram from the beginning to get first non-zero bin
  for (int bin=1; bin<=fineHisto->GetNbinsX(); bin++){
    double val = fineHisto->GetBinContent(bin);
    double valerr = fineHisto->GetBinError(bin);
    if (val==0. && valerr==0.) continue;
    if (bin!=1) res.addBinBoundary(fineHisto->GetXaxis()->GetBinLowEdge(bin));
    break;
  }
  // Loop over histogram from the end to get last non-zero bin
  for (int bin=fineHisto->GetNbinsX(); bin>=1; bin--){
    double val = fineHisto->GetBinContent(bin);
    double valerr = fineHisto->GetBinError(bin);
    if (val==0. && valerr==0.) continue;
    if (bin!=fineHisto->GetNbinsX()) res.addBinBoundary(fineHisto->GetXaxis()->GetBinLowEdge(bin+1));
    break;
  }
  return res;
}

void getIntermediateBinning(TH2F const* fineHisto, ExtendedBinning& xbinning, ExtendedBinning& ybinning, int effectiveNx, int effectiveNy){
  // Do not include overflow/underflow bins; they should not be used for determining the intermediate binning, which doesn't cover the under/overflows.
  TH1F* fineHistoX = getHistogramSlice(fineHisto, 0, 1, fineHisto->GetNbinsY(), "fineHistoX");
  TH1F* fineHistoY = getHistogramSlice(fineHisto, 1, 1, fineHisto->GetNbinsX(), "fineHistoY");

  xbinning=getIntermediateBinning(fineHistoX, effectiveNx);
  ybinning=getIntermediateBinning(fineHistoY, effectiveNy);

  delete fineHistoY;
  delete fineHistoX;
}


TH1F* getSmoothHistogram(
  TTree* tree, TString const hname, TString const htitle, float& weight,
  float& xvar, ExtendedBinning const& finalXBinning,
  int effectiveN=-1
){
  if (!tree) return nullptr;

  // Get fine binning to determine intermediate binning
  if (!finalXBinning.isValid()) return nullptr;
  unsigned int finalNbinsX;
  ExtendedBinning finestXBinning = getFinestBinning(finalXBinning, finalNbinsX);

  TH1F* res=nullptr;
  // Construct fine histogram to determine intermediate binning
  res = new TH1F("fineHisto", finestXBinning.getLabel(), finestXBinning.getNbins(), finestXBinning.getBinning());
  res->Sumw2();
  // First loop over the tree. Notice that the tree is supposed to only contain events that should finally go into the histogram, so create an intermediate tree if necessary
  // xvar and weight should already be pointed by the tree
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    res->Fill(xvar, weight);
  }
  ExtendedBinning intermediateXBinning = getIntermediateBinning(res, effectiveN);
  delete res;

  res = new TH1F(hname, htitle, intermediateXBinning.getNbins(), intermediateXBinning.getBinning());
  res->Sumw2();
  // Second loop over the tree
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    res->Fill(xvar, weight);
  }
  rebinHistogram(res, finalXBinning);

  return res;
}

TH2F* getSmoothHistogram(
  TTree* tree, TString const hname, TString const htitle, float& weight,
  float& xvar, ExtendedBinning const& finalXBinning,
  float& yvar, ExtendedBinning const& finalYBinning,
  int effectiveNx=-1, int effectiveNy=-1
){
  if (!tree) return nullptr;

  // Get fine binning to determine intermediate binning
  if (!finalXBinning.isValid()) return nullptr;
  unsigned int finalNbinsX;
  ExtendedBinning finestXBinning = getFinestBinning(finalXBinning, finalNbinsX);

  if (!finalYBinning.isValid()) return nullptr;
  unsigned int finalNbinsY;
  ExtendedBinning finestYBinning = getFinestBinning(finalYBinning, finalNbinsY);

  MELAout << "getSmoothHistogram: Final X binning requested: [ " << finalXBinning.getBinningVector() << " ]" << endl;
  //MELAout << "\t- Fine X binning: [ " << finestXBinning.getBinningVector() << " ]" << endl;
  MELAout << "getSmoothHistogram: Final Y binning requested: [ " << finalYBinning.getBinningVector() << " ]" << endl;
  //MELAout << "\t- Fine Y binning: [ " << finestYBinning.getBinningVector() << " ]" << endl;

  TH2F* res=nullptr;
  // Construct fine histogram to determine intermediate binning
  res = new TH2F("fineHisto", finestXBinning.getLabel()+":"+finestYBinning.getLabel(), finestXBinning.getNbins(), finestXBinning.getBinning(), finestYBinning.getNbins(), finestYBinning.getBinning());
  res->Sumw2();
  // First loop over the tree. Notice that the tree is supposed to only contain events that should finally go into the histogram, so create an intermediate tree if necessary
  // xvar and weight should already be pointed by the tree
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    res->Fill(xvar, yvar, weight);
  }
  ExtendedBinning intermediateXBinning(finestXBinning.getLabel()), intermediateYBinning(finestYBinning.getLabel());
  getIntermediateBinning(res, intermediateXBinning, intermediateYBinning, effectiveNx, effectiveNy);
  delete res;

  MELAout << "getSmoothHistogram: Intermediate X binning: [ " << intermediateXBinning.getBinningVector() << " ]" << endl;
  MELAout << "getSmoothHistogram: Intermediate Y binning: [ " << intermediateYBinning.getBinningVector() << " ]" << endl;

  res = new TH2F(hname, htitle, intermediateXBinning.getNbins(), intermediateXBinning.getBinning(), intermediateYBinning.getNbins(), intermediateYBinning.getBinning());
  res->Sumw2();
  // Second loop over the tree
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    res->Fill(xvar, yvar, weight);
  }
  rebinHistogram(res, finalXBinning, finalYBinning);

  return res;
}


#endif
