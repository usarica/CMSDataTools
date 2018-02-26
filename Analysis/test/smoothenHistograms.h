#ifndef SMOOTHENHISTOGRAMS_H
#define SMOOTHENHISTOGRAMS_H

#include "common_includes.h"


ExtendedBinning getFinestBinning(ExtendedBinning const& finalBinning, unsigned int& finalNbins){
  finalNbins=finalBinning.getNbins();
  float varmin=finalBinning.getBinLowEdge(0);
  float varmax=finalBinning.getBinLowEdge(finalNbins);

  TString varname=finalBinning.getLabel(); TString varnamelower=varname; varnamelower.ToLower();
  unsigned int nBinsIntermediate=1000*finalNbins;
  if (varnamelower.Contains("mass")) nBinsIntermediate = (varmax-varmin)*100.; // Bin mass every 1 MeV, notice mass unit is GeV
  return ExtendedBinning(nBinsIntermediate, varmin, varmax, varname);
}

ExtendedBinning getIntermediateBinning(TH1F const* fineHisto, unsigned int finalNbinsRequested){
  ExtendedBinning res;
  double integral=fineHisto->Integral(1, fineHisto->GetNbinsX()); // Do not include underflow/overflow bins; they are not relevant for the binning of the cumulant
  double sumThreshold;
  unsigned int nDivisions = std::min(finalNbinsRequested, (unsigned int) 20);
  sumThreshold=integral/double(nDivisions);

  res.addBinBoundary(fineHisto->GetXaxis()->GetBinLowEdge(0));
  vector<float> sums;
  float sum=0;
  for (int bin=1; bin<=fineHisto->GetNbinsX(); bin++){
    sum += fineHisto->GetBinContent(bin);
    if (sum>sumThreshold || bin==fineHisto->GetNbinsX()){
      sums.push_back(sum);
      res.addBinBoundary(fineHisto->GetXaxis()->GetBinLowEdge(bin+1));
      sum=0;
    }
  }
  // These lines guarantee sum>sumThreshold in every bin
  if (sums.size()>1 && sums.at(sums.size()-1)<sumThreshold){
    sums.at(sums.size()-1) += sums.at(sums.size()-2);
    sums.erase(sums.begin()+sums.size()-2);
    res.removeBinLowEdge(res.getNbins());
  }
  return res;
}

void getIntermediateBinning(TH2F const* fineHisto, ExtendedBinning& xbinning, ExtendedBinning& ybinning, unsigned int finalNbinsRequestedX, unsigned int finalNbinsRequestedY){
  const TAxis* xaxis=fineHisto->GetXaxis();
  const TAxis* yaxis=fineHisto->GetYaxis();
  vector<float> xbins, ybins;
  for (int ix=1; ix<=xaxis->GetNbins()+1; ix++) xbins.push_back(xaxis->GetBinLowEdge(ix));
  for (int iy=1; iy<=yaxis->GetNbins()+1; iy++) ybins.push_back(yaxis->GetBinLowEdge(iy));

  TH1F* fineHistoX = new TH1F("fineHistoX", "", xbins.size()-1, xbins.data());
  for (int ix=0; ix<=xaxis->GetNbins()+1; ix++){
    double integral=0, integralerror=0;
    integral = fineHisto->IntegralAndError(ix, ix, 0, yaxis->GetNbins()+1, integralerror);
    fineHistoX->SetBinContent(ix, integral);
    fineHistoX->SetBinError(ix, integralerror);
  }
  TH1F* fineHistoY = new TH1F("fineHistoY", "", ybins.size()-1, ybins.data());
  for (int iy=0; iy<=yaxis->GetNbins()+1; iy++){
    double integral=0, integralerror=0;
    integral = fineHisto->IntegralAndError(0, xaxis->GetNbins()+1, iy, iy, integralerror);
    fineHistoY->SetBinContent(iy, integral);
    fineHistoY->SetBinError(iy, integralerror);
  }

  xbinning=getIntermediateBinning(fineHistoX, finalNbinsRequestedX);
  ybinning=getIntermediateBinning(fineHistoY, finalNbinsRequestedY);

  delete fineHistoY;
  delete fineHistoX;
}


TH1F* getSmoothHistogram(
  TTree* tree, TString const hname, TString const htitle, float& weight,
  float& xvar, ExtendedBinning const& finalXBinning
){
  if (!tree) return nullptr;

  // Get fine binning to determine intermediate binning
  if (!finalXBinning.isValid()) return nullptr;
  unsigned int finalNbinsX;
  ExtendedBinning finestXBinning = getFinestBinning(finalXBinning, finalNbinsX);

  TH1F* res=nullptr;
  // Construct fine histogram to determine intermediate binning
  res = new TH1F("fineHisto", finestXBinning.getLabel(), finestXBinning.getNbins(), finestXBinning.getBinning());
  // Loop over the tree. Notice that the tree is supposed to only contain events that should finally go into the histogram, so create an intermediate tree if necessary
  // xvar and weight should already be pointed by the tree
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    res->Fill(xvar, weight);
  }
  ExtendedBinning intermediateXBinning = getIntermediateBinning(res, finalNbinsX);
  delete res;

  res = new TH1F(hname, htitle, intermediateXBinning.getNbins(), intermediateXBinning.getBinning());

  rebinHistogram(res, finalXBinning);

  return res;
}

TH2F* getSmoothHistogram(
  TTree* tree, TString const hname, TString const htitle, float& weight,
  float& xvar, ExtendedBinning const& finalXBinning,
  float& yvar, ExtendedBinning const& finalYBinning
){
  if (!tree) return nullptr;

  // Get fine binning to determine intermediate binning
  if (!finalXBinning.isValid()) return nullptr;
  unsigned int finalNbinsX;
  ExtendedBinning finestXBinning = getFinestBinning(finalXBinning, finalNbinsX);

  if (!finalYBinning.isValid()) return nullptr;
  unsigned int finalNbinsY;
  ExtendedBinning finestYBinning = getFinestBinning(finalYBinning, finalNbinsY);

  TH2F* res=nullptr;
  // Construct fine histogram to determine intermediate binning
  res = new TH2F("fineHisto", finestXBinning.getLabel()+":"+finestYBinning.getLabel(), finestXBinning.getNbins(), finestXBinning.getBinning(), finestYBinning.getNbins(), finestYBinning.getBinning());
  // Loop over the tree. Notice that the tree is supposed to only contain events that should finally go into the histogram, so create an intermediate tree if necessary
  // xvar and weight should already be pointed by the tree
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    res->Fill(xvar, yvar, weight);
  }
  ExtendedBinning intermediateXBinning(finestXBinning.getLabel()), intermediateYBinning(finestYBinning.getLabel());
  getIntermediateBinning(res, intermediateXBinning, intermediateYBinning, finalNbinsX, finalNbinsY);
  delete res;

  res = new TH2F(hname, htitle, intermediateXBinning.getNbins(), intermediateXBinning.getBinning(), intermediateYBinning.getNbins(), intermediateYBinning.getBinning());

  rebinHistogram(res, finalXBinning, finalYBinning);

  return res;
}


#endif
