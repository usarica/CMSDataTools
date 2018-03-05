#ifndef SMOOTHENHISTOGRAMS_H
#define SMOOTHENHISTOGRAMS_H

#include "common_includes.h"


ExtendedBinning getFinestBinning(ExtendedBinning const& finalBinning){
  unsigned int finalNbins=finalBinning.getNbins();
  float varmin=finalBinning.getBinLowEdge(0);
  float varmax=finalBinning.getBinLowEdge(finalNbins);

  TString varname=finalBinning.getLabel(); TString varnamelower=varname; varnamelower.ToLower();
  unsigned int nBinsIntermediate;
  if (varnamelower.Contains("mass")) nBinsIntermediate = 1000*finalNbins;
  else nBinsIntermediate = std::max((unsigned int) 100, 10*finalNbins);
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
  /*
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
  */
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

void getIntermediateBinning(TH3F const* fineHisto, ExtendedBinning& xbinning, ExtendedBinning& ybinning, ExtendedBinning& zbinning, int effectiveNx, int effectiveNy, int effectiveNz){
  // Do not include overflow/underflow bins; they should not be used for determining the intermediate binning, which doesn't cover the under/overflows.
  TH1F* fineHistoX = getHistogramSlice(fineHisto, 0, 1, fineHisto->GetNbinsY(), 1, fineHisto->GetNbinsZ(), "fineHistoX");
  TH1F* fineHistoY = getHistogramSlice(fineHisto, 1, 1, fineHisto->GetNbinsZ(), 1, fineHisto->GetNbinsX(), "fineHistoY");
  TH1F* fineHistoZ = getHistogramSlice(fineHisto, 2, 1, fineHisto->GetNbinsX(), 1, fineHisto->GetNbinsY(), "fineHistoZ");

  xbinning=getIntermediateBinning(fineHistoX, effectiveNx);
  ybinning=getIntermediateBinning(fineHistoY, effectiveNy);
  zbinning=getIntermediateBinning(fineHistoZ, effectiveNz);

  delete fineHistoZ;
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
  unsigned int finalNbinsX=finalXBinning.getNbins();
  ExtendedBinning finestXBinning = getFinestBinning(finalXBinning);

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
  int effectiveNx=-1, int effectiveNy=-1,
  bool smoothCondX=false, bool smoothCondY=false
){
  if (!tree) return nullptr;

  // Get fine binning to determine intermediate binning
  if (!finalXBinning.isValid()) return nullptr;
  unsigned int finalNbinsX = finalXBinning.getNbins();
  //ExtendedBinning finestXBinning = (!smoothCondX ? getFinestBinning(finalXBinning) : finalXBinning);
  ExtendedBinning finestXBinning = getFinestBinning(finalXBinning);
  finestXBinning.setLabel(finalXBinning.getLabel()+"_fine");

  if (!finalYBinning.isValid()) return nullptr;
  unsigned int finalNbinsY = finalYBinning.getNbins();
  //ExtendedBinning finestYBinning = (!smoothCondY ? getFinestBinning(finalYBinning) : finalYBinning);
  ExtendedBinning finestYBinning = getFinestBinning(finalYBinning);
  finestYBinning.setLabel(finalYBinning.getLabel()+"_fine");

  MELAout << "getSmoothHistogram: Final X binning requested: [ " << finalXBinning.getBinningVector() << " ]" << endl;
  //MELAout << "\t- Fine X binning: [ " << finestXBinning.getBinningVector() << " ]" << endl;
  MELAout << "getSmoothHistogram: Final Y binning requested: [ " << finalYBinning.getBinningVector() << " ]" << endl;
  //MELAout << "\t- Fine Y binning: [ " << finestYBinning.getBinningVector() << " ]" << endl;

  // Construct fine histogram to determine intermediate binning
  ExtendedHistogram_2D res("res", "", finestXBinning, finestYBinning);
  ExtendedHistogram_2D resUns("resUns", "", finalXBinning, finalYBinning);

  res.constructFromTree(tree, xvar, yvar, weight);

  ExtendedBinning intermediateXBinning(finalXBinning.getLabel()), intermediateYBinning(finalYBinning.getLabel());
  getIntermediateBinning(res.getHistogram(), intermediateXBinning, intermediateYBinning, effectiveNx, effectiveNy);
  intermediateXBinning.setLabel(finalXBinning.getLabel()+"_interm");
  intermediateYBinning.setLabel(finalYBinning.getLabel()+"_interm");

  if (!smoothCondX && !smoothCondY){
    MELAout << "getSmoothHistogram: Intermediate X binning: [ " << intermediateXBinning.getBinningVector() << " ]" << endl;
    MELAout << "getSmoothHistogram: Intermediate Y binning: [ " << intermediateYBinning.getBinningVector() << " ]" << endl;

    res.constructFromTree(tree, xvar, yvar, weight, &intermediateXBinning, &intermediateYBinning);
    res.rebin(&finalXBinning, &finalYBinning);

    res.setNameTitle(hname, htitle);
    TH2F* hres = new TH2F(*(res.getHistogram()));
    return hres;
  }
  else if (smoothCondX){
    MELAout << "getSmoothHistogram: Intermediate X binning: [ " << intermediateXBinning.getBinningVector() << " ]" << endl;

    resUns.constructFromTree(tree, xvar, yvar, weight);
    TH1F* hProjOutX = getHistogramSlice(resUns.getHistogram(), 0, 0, resUns.getHistogram()->GetNbinsY()+1, "hProjOutX");

    MELAout << "\t- Rebinning res" << endl;
    res.constructFromTree(tree, xvar, yvar, weight, &intermediateXBinning, &finestYBinning);
    MELAout << "\t- Conditionalizing res" << endl;
    conditionalizeHistogram<TH2F>(res.getHistogram(), 0, nullptr, false);

    TH2F* hres = new TH2F("hres", "", intermediateXBinning.getNbins(), intermediateXBinning.getBinning(), finalYBinning.getNbins(), finalYBinning.getBinning());
    for (unsigned int ix=0; ix<=intermediateXBinning.getNbins()+1; ix++){
      MELAout << "\t- Constructing slice " << ix << endl;

      TH1F* hProj = getHistogramSlice(res.getHistogram(), 1, ix, ix, "hProj");
      ExtendedBinning intermediateBinning = getIntermediateBinning(hProj, effectiveNy);
      MELAout << "\t- Intermediate Y binning: [ " << intermediateBinning.getBinningVector() << " ]" << endl;
      rebinHistogram(hProj, intermediateBinning);
      rebinHistogram(hProj, finalYBinning);
      for (unsigned int iy=0; iy<=finalYBinning.getNbins()+1; iy++){
        double bincontent = hProj->GetBinContent(iy);
        double binerror = hProj->GetBinError(iy);
        hres->SetBinContent(ix, iy, bincontent);
        hres->SetBinError(ix, iy, binerror);
      }
      delete hProj;
    }
    res.rebin(nullptr, &finalYBinning);
    {
      vector<pair<TProfile const*, unsigned int>> condProfs;
      condProfs.emplace_back(res.getProfileX(), 0);
      rebinHistogram(hres, finalXBinning, finalYBinning, &condProfs);
    }
    //rebinHistogram_NoCumulant(hres, finalXBinning, res.getProfileX(), finalYBinning, res.getProfileY());
    for (unsigned int ix=0; ix<=finalXBinning.getNbins()+1; ix++){
      double sumWAll = hProjOutX->GetBinContent(ix);
      double sumWsqAll = pow(hProjOutX->GetBinError(ix), 2);
      for (unsigned int iy=0; iy<=finalYBinning.getNbins()+1; iy++){
        double bincontent = hres->GetBinContent(ix, iy);
        double binerror = hres->GetBinError(ix, iy);
        hres->SetBinContent(ix, iy, bincontent*sumWAll);
        hres->SetBinError(ix, iy, translateEfficiencyErrorToNumeratorError(bincontent, sumWAll, binerror, sumWsqAll));
      }
    }
    delete hProjOutX;
    hres->SetName(hname);
    hres->SetTitle(htitle);
    return hres;
  }
  else if (smoothCondY){
    MELAout << "getSmoothHistogram: Intermediate Y binning: [ " << intermediateYBinning.getBinningVector() << " ]" << endl;

    resUns.constructFromTree(tree, xvar, yvar, weight);
    TH1F* hProjOutY = getHistogramSlice(resUns.getHistogram(), 1, 0, resUns.getHistogram()->GetNbinsX()+1, "hProjOutY");

    MELAout << "\t- Rebinning res" << endl;
    res.constructFromTree(tree, xvar, yvar, weight, &finestXBinning, &intermediateYBinning);
    MELAout << "\t- Conditionalizing res" << endl;
    conditionalizeHistogram<TH2F>(res.getHistogram(), 1, nullptr, false);

    TH2F* hres = new TH2F("hres", "", finalXBinning.getNbins(), finalXBinning.getBinning(), intermediateYBinning.getNbins(), intermediateYBinning.getBinning());
    for (unsigned int iy=0; iy<=intermediateYBinning.getNbins()+1; iy++){
      MELAout << "\t- Constructing slice " << iy << endl;

      TH1F* hProj = getHistogramSlice(res.getHistogram(), 0, iy, iy, "hProj");
      ExtendedBinning intermediateBinning = getIntermediateBinning(hProj, effectiveNx);
      MELAout << "\t- Intermediate X binning: [ " << intermediateBinning.getBinningVector() << " ]" << endl;
      rebinHistogram(hProj, intermediateBinning);
      rebinHistogram(hProj, finalXBinning);
      for (unsigned int ix=0; ix<=finalXBinning.getNbins()+1; ix++){
        double bincontent = hProj->GetBinContent(ix);
        double binerror = hProj->GetBinError(ix);
        hres->SetBinContent(ix, iy, bincontent);
        hres->SetBinError(ix, iy, binerror);
      }
      delete hProj;
    }
    res.rebin(&finalXBinning, nullptr);
    {
      vector<pair<TProfile const*, unsigned int>> condProfs;
      condProfs.emplace_back(res.getProfileY(), 1);
      rebinHistogram(hres, finalXBinning, finalYBinning, &condProfs);
    }
    //rebinHistogram_NoCumulant(hres, finalXBinning, res.getProfileX(), finalYBinning, res.getProfileY());
    for (unsigned int iy=0; iy<=finalYBinning.getNbins()+1; iy++){
      double sumWAll = hProjOutY->GetBinContent(iy);
      double sumWsqAll = pow(hProjOutY->GetBinError(iy), 2);
      for (unsigned int ix=0; ix<=finalXBinning.getNbins()+1; ix++){
        double bincontent = hres->GetBinContent(ix, iy);
        double binerror = hres->GetBinError(ix, iy);
        hres->SetBinContent(ix, iy, bincontent*sumWAll);
        hres->SetBinError(ix, iy, translateEfficiencyErrorToNumeratorError(bincontent, sumWAll, binerror, sumWsqAll));
      }
    }
    delete hProjOutY;
    hres->SetName(hname);
    hres->SetTitle(htitle);
    return hres;
  }
  else{
    res.rebin(&finalXBinning, &finalYBinning);
    res.setNameTitle(hname, htitle);
    TH2F* hres = new TH2F(*(res.getHistogram()));
    return hres;
  }
}


#endif
