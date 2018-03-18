#ifndef HISTOGRAMSMOOTHER_H
#define HISTOGRAMSMOOTHER_H

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>
#include "TFile.h"
#include "TString.h"
#include "TSpline.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "HelperFunctions.h"
#include "SampleHelpers.h"
#include "DiscriminantClasses.h"
#include "CategorizationHelpers.h"
#include "ExtendedHistogram_1D.h"
#include "ExtendedHistogram_2D.h"
#include "ExtendedHistogram_3D.h"
#include "TemplateHelpers.h"
#include "MELAStreamHelpers.hh"
#include "Mela.h"


using namespace std;
using namespace HelperFunctions;
using namespace SampleHelpers;
using namespace DiscriminantClasses;
using namespace TemplateHelpers;
using namespace MELAStreamHelpers;


ExtendedBinning getFinestBinning(ExtendedBinning const& finalBinning){
  unsigned int finalNbins=finalBinning.getNbins();
  float varmin=finalBinning.getBinLowEdge(0);
  float varmax=finalBinning.getBinLowEdge(finalNbins);

  TString varname=finalBinning.getLabel()+"_fine"; TString varnamelower=varname; varnamelower.ToLower();
  unsigned int nBinsIntermediate;
  if (varnamelower.Contains("mass")){
    if ((finalBinning.getMax()-finalBinning.getMin())>1000.) nBinsIntermediate = 1000*finalNbins;
    else nBinsIntermediate = 100*finalNbins;
  }
  else nBinsIntermediate = std::max((unsigned int) 100, 10*finalNbins);
  return ExtendedBinning(nBinsIntermediate, varmin, varmax, varname);
}

ExtendedBinning getIntermediateBinning(TH1F const* fineHisto, ExtendedBinning const& finalBinning, int effectiveN){
  TString bname = finalBinning.getLabel(); /*TString bnamelow = bname; bnamelow.ToLower();*/
  ExtendedBinning res(bname+"_interm");
  unsigned int finalNbins=finalBinning.getNbins();

  double integralerror=0;
  double integral=fineHisto->IntegralAndError(1, fineHisto->GetNbinsX(), integralerror); // Do not include underflow/overflow bins; they are not relevant for the binning of the cumulant
  double effectiveNtotal=(integralerror>0. ? pow(integral/integralerror, 2) : fineHisto->GetEntries());
  assert(effectiveNtotal!=0.);
  double combinedErrorMargin=1./sqrt(effectiveNtotal);
  double perbinErrorMarginPossible=std::min(1., combinedErrorMargin*sqrt(double(finalNbins)));

  double effNThreshold=0;
  unsigned int appNbins=0;
  if (effectiveN>0){
    double perbinErrorMarginReq=1./sqrt(double(effectiveN));
    if (perbinErrorMarginPossible<=perbinErrorMarginReq) effNThreshold=effectiveN;
    appNbins=effectiveNtotal/effNThreshold;
  }
  if (effNThreshold==0.){
    if (perbinErrorMarginPossible<=0.1){
      effNThreshold=pow(1./perbinErrorMarginPossible, 2);
      appNbins=effectiveNtotal/effNThreshold;
    }
    else{
      effNThreshold=effectiveNtotal/4.;
      appNbins=effectiveNtotal/effNThreshold; // =4
    }
  }
  double sumWThreshold = integral/double(appNbins);
  double sumWsqThreshold = pow(integralerror, 2)/double(appNbins);
  MELAout << "Effective total N = " << effectiveNtotal << ", requested effective N = " << effectiveN << ", threshold = " << effNThreshold << ", approximate Nbins=" << appNbins << endl;

  ExtendedBinning fineBinning=ExtendedBinning::extractBinning(fineHisto, 0);
  int fineNbins=fineBinning.getNbins();
  vector<double> sumWList, sumWsqList, effNList;
  double sumW=0, sumWsq=0;
  bool hasSumW2 = (fineHisto->GetSumw2N()!=0);
  // Loop over the fine bins
  res.addBinBoundary(fineBinning.getBinLowEdge(0));
  for (int bin=1; bin<=fineNbins; bin++){
    double berrsq=pow(fineHisto->GetBinError(bin), 2);
    double bval=fineHisto->GetBinContent(bin);
    if (berrsq>0.){
      sumW += bval;
      sumWsq += berrsq;
    }
    double effN;
    if (hasSumW2) effN=(sumWsq>0. ? pow(sumW, 2)/sumWsq : 0);
    else effN=sumW/integral*effectiveNtotal;
    if (bin==fineNbins || (effN>=effNThreshold && sumW>=sumWThreshold && sumWsq>=sumWsqThreshold)){
      effNList.push_back(effN);
      sumWList.push_back(sumW);
      sumWsqList.push_back(sumWsq);
      res.addBinBoundary(fineBinning.getBinLowEdge(bin));
      sumW=0;
      sumWsq=0;
    }
  }
  // These lines guarantee sum>sumThreshold in every bin
  assert(effNList.size()==sumWList.size() && effNList.size()==sumWsqList.size());
  if (
    effNList.size()>1
    && (
      effNList.at(effNList.size()-1)<effNThreshold
      ||
      sumWList.at(sumWList.size()-1)<sumWThreshold
      ||
      sumWsqList.at(sumWsqList.size()-1)<sumWsqThreshold
      )
    ){
    sumWList.at(sumWList.size()-1) += sumWList.at(sumWList.size()-2);
    sumWList.erase(sumWList.begin()+sumWList.size()-2);
    sumWsqList.at(sumWsqList.size()-1) += sumWsqList.at(sumWsqList.size()-2);
    sumWsqList.erase(sumWsqList.begin()+sumWsqList.size()-2);
    effNList.clear(); for (unsigned int i=0; i<sumWList.size(); i++) effNList.push_back(pow(sumWList.at(i), 2)/sumWList.at(i));
    res.removeBinLowEdge(res.getNbins()-1);
  }
  return res;
}

void getIntermediateBinning(
  TH2F const* fineHisto,
  ExtendedBinning const& finalXBinning, ExtendedBinning const& finalYBinning,
  std::vector<ExtendedBinning>& xbinninglist, std::vector<ExtendedBinning>& ybinninglist,
  int effectiveNx, int effectiveNy,
  int condDim
){
  // Do not include overflow/underflow bins; they should not be used for determining the intermediate binning, which doesn't cover the under/overflows.
  switch (condDim){
  case 0:
  {
    TH1F* fineHistoX = getHistogramSlice(fineHisto, 0, 1, fineHisto->GetNbinsY(), "fineHistoX");
    ExtendedBinning xbinning=getIntermediateBinning(fineHistoX, finalXBinning, effectiveNx);
    xbinninglist.push_back(xbinning);
    delete fineHistoX;

    ExtendedBinning cdbinning_original = ExtendedBinning::extractBinning(fineHisto, condDim);
    for (int ibin=0; ibin<(int) xbinning.getNbins(); ibin++){
      double blow=xbinning.getBinLowEdge(ibin);
      double bhigh=xbinning.getBinHighEdge(ibin);
      int ix=cdbinning_original.getBin(blow);
      int jx=cdbinning_original.getBin(bhigh);
      if (jx>ix) jx--;
      TH1F* fineHistoY = getHistogramSlice(fineHisto, 1, ix, jx, "fineHistoY");
      ExtendedBinning ybinning=getIntermediateBinning(fineHistoY, finalYBinning, effectiveNy);
      ybinning.setLabel(ybinning.getLabel()+Form("_Xslice_%i", ibin));
      ybinninglist.push_back(ybinning);
      delete fineHistoY;
    }
    break;
  }
  case 1:
  {
    TH1F* fineHistoY = getHistogramSlice(fineHisto, 1, 1, fineHisto->GetNbinsX(), "fineHistoY");
    ExtendedBinning ybinning=getIntermediateBinning(fineHistoY, finalYBinning, effectiveNy);
    ybinninglist.push_back(ybinning);
    delete fineHistoY;

    ExtendedBinning cdbinning_original = ExtendedBinning::extractBinning(fineHisto, condDim);
    for (int ibin=0; ibin<(int) ybinning.getNbins(); ibin++){
      double blow=ybinning.getBinLowEdge(ibin);
      double bhigh=ybinning.getBinHighEdge(ibin);
      int iy=cdbinning_original.getBin(blow);
      int jy=cdbinning_original.getBin(bhigh);
      if (jy>iy) jy--;
      TH1F* fineHistoX = getHistogramSlice(fineHisto, 0, iy, jy, "fineHistoX");
      ExtendedBinning xbinning=getIntermediateBinning(fineHistoX, finalXBinning, effectiveNx);
      xbinning.setLabel(xbinning.getLabel()+Form("_Yslice_%i", ibin));
      xbinninglist.push_back(xbinning);
      delete fineHistoX;
    }

    break;
  }
  default:
  {
    TH1F* fineHistoX = getHistogramSlice(fineHisto, 0, 1, fineHisto->GetNbinsY(), "fineHistoX");
    ExtendedBinning xbinning=getIntermediateBinning(fineHistoX, finalXBinning, effectiveNx);
    xbinninglist.push_back(xbinning);
    delete fineHistoX;

    TH1F* fineHistoY = getHistogramSlice(fineHisto, 1, 1, fineHisto->GetNbinsX(), "fineHistoY");
    ExtendedBinning ybinning=getIntermediateBinning(fineHistoY, finalYBinning, effectiveNy);
    ybinninglist.push_back(ybinning);
    delete fineHistoY;

    break;
  }
  }
}

void getIntermediateBinning(
  TH3F const* fineHisto,
  ExtendedBinning const& finalXBinning, ExtendedBinning const& finalYBinning, ExtendedBinning const& finalZBinning,
  std::vector<ExtendedBinning>& xbinninglist, std::vector<ExtendedBinning>& ybinninglist, std::vector<ExtendedBinning>& zbinninglist,
  int effectiveNx, int effectiveNy, int effectiveNz,
  int condDim
){
  // Do not include overflow/underflow bins; they should not be used for determining the intermediate binning, which doesn't cover the under/overflows.
  switch (condDim){
  case 0:
  {
    TH1F* fineHistoX = getHistogramSlice(fineHisto, 0, 1, fineHisto->GetNbinsY(), 1, fineHisto->GetNbinsZ(), "fineHistoX");
    ExtendedBinning xbinning=getIntermediateBinning(fineHistoX, finalXBinning, effectiveNx);
    xbinninglist.push_back(xbinning);
    delete fineHistoX;

    ExtendedBinning cdbinning_original = ExtendedBinning::extractBinning(fineHisto, condDim);
    for (int ibin=0; ibin<(int) xbinning.getNbins(); ibin++){
      double blow=xbinning.getBinLowEdge(ibin);
      double bhigh=xbinning.getBinHighEdge(ibin);
      int ix=cdbinning_original.getBin(blow);
      int jx=cdbinning_original.getBin(bhigh);
      if (jx>ix) jx--;
      TH2F* fineHistoYZ = getHistogramSlice(fineHisto, 1, 2, ix, jx, "fineHistoYZ");
      getIntermediateBinning(
        fineHistoYZ,
        finalYBinning, finalZBinning,
        ybinninglist, zbinninglist,
        effectiveNy, effectiveNz,
        -1
      );
      for (auto& binning:ybinninglist) binning.setLabel(binning.getLabel()+Form("_Xslice_%i", ibin));
      for (auto& binning:zbinninglist) binning.setLabel(binning.getLabel()+Form("_Xslice_%i", ibin));
      delete fineHistoYZ;
    }
    break;
  }
  case 1:
  {
    TH1F* fineHistoY = getHistogramSlice(fineHisto, 1, 1, fineHisto->GetNbinsZ(), 1, fineHisto->GetNbinsX(), "fineHistoY");
    ExtendedBinning ybinning=getIntermediateBinning(fineHistoY, finalYBinning, effectiveNy);
    ybinninglist.push_back(ybinning);
    delete fineHistoY;

    ExtendedBinning cdbinning_original = ExtendedBinning::extractBinning(fineHisto, condDim);
    for (int ibin=0; ibin<(int) ybinning.getNbins(); ibin++){
      double blow=ybinning.getBinLowEdge(ibin);
      double bhigh=ybinning.getBinHighEdge(ibin);
      int iy=cdbinning_original.getBin(blow);
      int jy=cdbinning_original.getBin(bhigh);
      if (jy>iy) jy--;
      TH2F* fineHistoZX = getHistogramSlice(fineHisto, 2, 0, iy, jy, "fineHistoZX");
      getIntermediateBinning(
        fineHistoZX,
        finalZBinning, finalXBinning,
        zbinninglist, xbinninglist,
        effectiveNz, effectiveNx,
        -1
      );
      for (auto& binning:zbinninglist) binning.setLabel(binning.getLabel()+Form("_Yslice_%i", ibin));
      for (auto& binning:xbinninglist) binning.setLabel(binning.getLabel()+Form("_Yslice_%i", ibin));
      delete fineHistoZX;
    }

    break;
  }
  case 2:
  {
    TH1F* fineHistoZ = getHistogramSlice(fineHisto, 2, 1, fineHisto->GetNbinsX(), 1, fineHisto->GetNbinsY(), "fineHistoZ");
    ExtendedBinning zbinning=getIntermediateBinning(fineHistoZ, finalZBinning, effectiveNz);
    zbinninglist.push_back(zbinning);
    delete fineHistoZ;

    ExtendedBinning cdbinning_original = ExtendedBinning::extractBinning(fineHisto, condDim);
    for (int ibin=0; ibin<(int) zbinning.getNbins(); ibin++){
      double blow=zbinning.getBinLowEdge(ibin);
      double bhigh=zbinning.getBinHighEdge(ibin);
      int iz=cdbinning_original.getBin(blow);
      int jz=cdbinning_original.getBin(bhigh);
      if (jz>iz) jz--;
      TH2F* fineHistoXY = getHistogramSlice(fineHisto, 0, 1, iz, jz, "fineHistoXY");
      getIntermediateBinning(
        fineHistoXY,
        finalXBinning, finalYBinning,
        xbinninglist, ybinninglist,
        effectiveNx, effectiveNy,
        -1
      );
      for (auto& binning:xbinninglist) binning.setLabel(binning.getLabel()+Form("_Zslice_%i", ibin));
      for (auto& binning:ybinninglist) binning.setLabel(binning.getLabel()+Form("_Zslice_%i", ibin));
      delete fineHistoXY;
    }

    break;
  }
  default:
  {
    TH1F* fineHistoX = getHistogramSlice(fineHisto, 0, 1, fineHisto->GetNbinsY(), 1, fineHisto->GetNbinsZ(), "fineHistoX");
    ExtendedBinning xbinning=getIntermediateBinning(fineHistoX, finalXBinning, effectiveNx);
    xbinninglist.push_back(xbinning);
    delete fineHistoX;

    TH1F* fineHistoY = getHistogramSlice(fineHisto, 1, 1, fineHisto->GetNbinsZ(), 1, fineHisto->GetNbinsX(), "fineHistoY");
    ExtendedBinning ybinning=getIntermediateBinning(fineHistoY, finalYBinning, effectiveNy);
    ybinninglist.push_back(ybinning);
    delete fineHistoY;

    TH1F* fineHistoZ = getHistogramSlice(fineHisto, 2, 1, fineHisto->GetNbinsX(), 1, fineHisto->GetNbinsY(), "fineHistoZ");
    ExtendedBinning zbinning=getIntermediateBinning(fineHistoZ, finalZBinning, effectiveNz);
    zbinninglist.push_back(zbinning);
    delete fineHistoZ;

    break;
  }
  }
}


TH1F* getSmoothHistogram(
  TString const hname, TString const htitle, ExtendedBinning const& finalXBinning,
  TTree* tree, float& xvar, float& weight, bool& flag,
  int effectiveNx
){
  assert(tree && finalXBinning.isValid());

  // Get fine binning to determine intermediate binning
  unsigned int finalNbinsX=finalXBinning.getNbins();
  ExtendedBinning finestXBinning = getFinestBinning(finalXBinning);

  TH1F* res=nullptr;
  // Construct fine histogram to determine intermediate binning
  res = new TH1F(Form("fineHistoX_%s", finestXBinning.getLabel().Data()), "", finestXBinning.getNbins(), finestXBinning.getBinning());
  res->Sumw2();
  // First loop over the tree. Notice that the tree is supposed to only contain events that should finally go into the histogram, so create an intermediate tree if necessary
  // xvar and weight should already be pointed by the tree
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    if (flag) res->Fill(xvar, weight);
  }
  ExtendedBinning intermediateXBinning = getIntermediateBinning(res, finalXBinning, effectiveNx);
  delete res;

  res = new TH1F(hname, htitle, intermediateXBinning.getNbins(), intermediateXBinning.getBinning());
  res->Sumw2();
  // Second loop over the tree
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    if (flag) res->Fill(xvar, weight);
  }
  rebinHistogram(res, finalXBinning);

  return res;
}

TH2F* getSmoothHistogram(
  TString const hname, TString const htitle, ExtendedBinning const& finalXBinning, ExtendedBinning const& finalYBinning,
  TTree* tree, float& xvar, float& yvar, float& weight, bool& flag,
  int effectiveNx, int effectiveNy,
  int condDim=-1
){
  assert(tree && finalXBinning.isValid() && finalYBinning.isValid());

  // Get fine binning to determine intermediate binning
  unsigned int finalNbinsX = finalXBinning.getNbins();
  ExtendedBinning finestXBinning = getFinestBinning(finalXBinning);
  finestXBinning.setLabel(finalXBinning.getLabel()+"_fine");

  unsigned int finalNbinsY = finalYBinning.getNbins();
  ExtendedBinning finestYBinning = getFinestBinning(finalYBinning);
  finestYBinning.setLabel(finalYBinning.getLabel()+"_fine");

  // Construct fine histogram to determine intermediate binning
  ExtendedHistogram_2D res("res", "", finestXBinning, finestYBinning);
  res.constructFromTree(tree, xvar, yvar, weight, &flag);

  vector<ExtendedBinning> intermediateXBinningList, intermediateYBinningList;
  getIntermediateBinning(
    res.getHistogram(),
    finalXBinning, finalYBinning,
    intermediateXBinningList, intermediateYBinningList,
    effectiveNx, effectiveNy,
    condDim
  );

  switch (condDim){
  case 0:
  {
    ExtendedBinning& intermediateXBinning = intermediateXBinningList.at(0);
    int nIbinsX = intermediateXBinning.getNbins();
    res.constructFromTree(tree, xvar, yvar, weight, &flag, &intermediateXBinning, &finestYBinning);
    TH2F*& hIntermediate=res.getHistogram();
    vector<TH1F*> hRebinList; hRebinList.reserve(nIbinsX);
    vector<TProfile*> hRebinProfileList; hRebinProfileList.reserve(nIbinsX);
    for (int i=1; i<=nIbinsX; i++){
      ExtendedBinning& binning = intermediateYBinningList.at(i-1);
      TH1F* hRebin = getHistogramSlice(hIntermediate, 1, i, i, Form("hRebinIntermediate_SliceX%i", i));
      TProfile* hRebinProf = new TProfile(*(res.getProfileY())); hRebinProf->SetName(Form("%s_Slice%i", hRebinProf->GetName(), i));
      rebinHistogram(hRebin, binning);
      rebinHistogram(hRebin, finalYBinning); rebinProfile(hRebinProf, finalYBinning);
      hRebinList.push_back(hRebin); hRebinProfileList.push_back(hRebinProf);
    }
    res.constructFromTree(tree, xvar, yvar, weight, &flag, &intermediateXBinning, &finalYBinning);
    // Replace prof_y of res first
    for (int iy=0; iy<=(int) finalYBinning.getNbins()+1; iy++){
      double sumY=0, sumWY=0;
      for (int ix=0; ix<nIbinsX; ix++){
        double val=hRebinProfileList.at(ix)->GetBinContent(iy);
        double err=hRebinProfileList.at(ix)->GetBinError(iy);
        if (err>0.){ val/=pow(err, 2); sumY+=val; sumWY += 1./pow(err, 2); }
      }
      if (sumWY>0.){ sumY/=sumWY; sumWY = 1./sqrt(sumWY); }
      res.getProfileY()->SetBinContent(iy, sumY);
      res.getProfileY()->SetBinError(iy, sumWY);
    }
    for (int ix=1; ix<=nIbinsX; ix++){
      for (int iy=0; iy<=(int) finalYBinning.getNbins()+1; iy++){
        res.getHistogram()->SetBinContent(ix, iy, hRebinList.at(ix)->GetBinContent(iy));
        res.getHistogram()->SetBinError(ix, iy, hRebinList.at(ix)->GetBinError(iy));
      }
    }
    res.rebin(&finalXBinning, &finalYBinning);

    // Cleanup
    for (auto*& h:hRebinList) delete h;
    for (auto*& h:hRebinProfileList) delete h;

    res.setNameTitle(hname, htitle);
    TH2F* hres = new TH2F(*(res.getHistogram()));
    return hres;
    break;
  }
  case 1:
  {
    ExtendedBinning& intermediateYBinning = intermediateYBinningList.at(0);
    int nIbinsY = intermediateYBinning.getNbins();
    res.constructFromTree(tree, xvar, yvar, weight, &flag, &finestXBinning, &intermediateYBinning);
    TH2F*& hIntermediate=res.getHistogram();
    vector<TH1F*> hRebinList; hRebinList.reserve(nIbinsY);
    vector<TProfile*> hRebinProfileList; hRebinProfileList.reserve(nIbinsY);
    for (int i=1; i<=nIbinsY; i++){
      ExtendedBinning& binning = intermediateXBinningList.at(i-1);
      TH1F* hRebin = getHistogramSlice(hIntermediate, 0, i, i, Form("hRebinIntermediate_SliceX%i", i));
      TProfile* hRebinProf = new TProfile(*(res.getProfileX())); hRebinProf->SetName(Form("%s_Slice%i", hRebinProf->GetName(), i));
      rebinHistogram(hRebin, binning);
      rebinHistogram(hRebin, finalXBinning); rebinProfile(hRebinProf, finalXBinning);
      hRebinList.push_back(hRebin); hRebinProfileList.push_back(hRebinProf);
    }
    res.constructFromTree(tree, xvar, yvar, weight, &flag, &finalXBinning, &intermediateYBinning);
    // Replace prof_x of res first
    for (int ix=0; ix<=(int) finalXBinning.getNbins()+1; ix++){
      double sumX=0, sumWX=0;
      for (int iy=0; iy<nIbinsY; iy++){
        double val=hRebinProfileList.at(iy)->GetBinContent(ix);
        double err=hRebinProfileList.at(iy)->GetBinError(ix);
        if (err>0.){ val/=pow(err, 2); sumX+=val; sumWX += 1./pow(err, 2); }
      }
      if (sumWX>0.){ sumX/=sumWX; sumWX = 1./sqrt(sumWX); }
      res.getProfileX()->SetBinContent(ix, sumX);
      res.getProfileX()->SetBinError(ix, sumWX);
    }
    for (int iy=1; iy<=nIbinsY; iy++){
      for (int ix=0; ix<=(int) finalXBinning.getNbins()+1; ix++){
        res.getHistogram()->SetBinContent(ix, iy, hRebinList.at(iy)->GetBinContent(ix));
        res.getHistogram()->SetBinError(ix, iy, hRebinList.at(iy)->GetBinError(ix));
      }
    }
    res.rebin(&finalXBinning, &finalYBinning);

    // Cleanup
    for (auto*& h:hRebinList) delete h;
    for (auto*& h:hRebinProfileList) delete h;

    res.setNameTitle(hname, htitle);
    TH2F* hres = new TH2F(*(res.getHistogram()));
    return hres;
    break;
  }
  default:
  {
    res.constructFromTree(tree, xvar, yvar, weight, &flag, &(intermediateXBinningList.at(0)), &(intermediateYBinningList.at(0)));
    res.rebin(&finalXBinning, &finalYBinning);

    res.setNameTitle(hname, htitle);
    TH2F* hres = new TH2F(*(res.getHistogram()));
    return hres;
  }
  }
}

TH3F* getSmoothHistogram(
  TString const hname, TString const htitle, ExtendedBinning const& finalXBinning, ExtendedBinning const& finalYBinning, ExtendedBinning const& finalZBinning,
  TTree* tree, float& xvar, float& yvar, float& zvar, float& weight, bool& flag,
  int effectiveNx, int effectiveNy, int effectiveNz,
  int condDim=-1
){
  assert(tree && finalXBinning.isValid() && finalYBinning.isValid() && finalZBinning.isValid());

  // Get fine binning to determine intermediate binning
  unsigned int finalNbinsX = finalXBinning.getNbins();
  ExtendedBinning finestXBinning = getFinestBinning(finalXBinning);
  finestXBinning.setLabel(finalXBinning.getLabel()+"_fine");

  unsigned int finalNbinsY = finalYBinning.getNbins();
  ExtendedBinning finestYBinning = getFinestBinning(finalYBinning);
  finestYBinning.setLabel(finalYBinning.getLabel()+"_fine");

  unsigned int finalNbinsZ = finalZBinning.getNbins();
  ExtendedBinning finestZBinning = getFinestBinning(finalZBinning);
  finestZBinning.setLabel(finalZBinning.getLabel()+"_fine");

  // Construct fine histogram to determine intermediate binning
  ExtendedHistogram_3D res("res", "", finestXBinning, finestYBinning, finestZBinning);
  res.constructFromTree(tree, xvar, yvar, zvar, weight, &flag);

  vector<ExtendedBinning> intermediateXBinningList, intermediateYBinningList, intermediateZBinningList;
  getIntermediateBinning(
    res.getHistogram(),
    finalXBinning, finalYBinning, finalZBinning,
    intermediateXBinningList, intermediateYBinningList, intermediateZBinningList,
    effectiveNx, effectiveNy, effectiveNz,
    condDim
  );

  switch (condDim){
  case 0:
  {
    ExtendedBinning& intermediateXBinning = intermediateXBinningList.at(0);
    int nIbinsX = intermediateXBinning.getNbins();
    res.constructFromTree(tree, xvar, yvar, zvar, weight, &flag, &intermediateXBinning, &finestYBinning, &finestZBinning);
    TH3F*& hIntermediate=res.getHistogram();
    vector<TH2F*> hRebinList; hRebinList.reserve(nIbinsX);
    vector<TProfile*> hRebinProfileYList, hRebinProfileZList; hRebinProfileYList.reserve(nIbinsX); hRebinProfileZList.reserve(nIbinsX);
    for (int i=1; i<=nIbinsX; i++){
      ExtendedBinning& binningY = intermediateYBinningList.at(i-1);
      ExtendedBinning& binningZ = intermediateZBinningList.at(i-1);
      TH2F* hRebin = getHistogramSlice(hIntermediate, 1, 2, i, i, Form("hRebinIntermediate_SliceX%i", i));
      TProfile* hRebinProfY = new TProfile(*(res.getProfileY())); hRebinProfY->SetName(Form("%s_Slice%i", hRebinProfY->GetName(), i));
      TProfile* hRebinProfZ = new TProfile(*(res.getProfileZ())); hRebinProfZ->SetName(Form("%s_Slice%i", hRebinProfZ->GetName(), i));
      rebinHistogram(hRebin, binningY, binningZ); rebinHistogram(hRebin, finalYBinning, finalZBinning);
      rebinProfile(hRebinProfY, finalYBinning); rebinProfile(hRebinProfZ, finalZBinning);
      hRebinList.push_back(hRebin);
      hRebinProfileYList.push_back(hRebinProfY);
      hRebinProfileZList.push_back(hRebinProfZ);
    }
    res.constructFromTree(tree, xvar, yvar, zvar, weight, &flag, &intermediateXBinning, &finalYBinning, &finalZBinning);
    // Replace prof_y of res
    for (int iy=0; iy<=(int) finalYBinning.getNbins()+1; iy++){
      double sum=0, sumW=0;
      for (int ix=0; ix<nIbinsX; ix++){
        double val=hRebinProfileYList.at(ix)->GetBinContent(iy);
        double err=hRebinProfileYList.at(ix)->GetBinError(iy);
        if (err>0.){ val/=pow(err, 2); sum+=val; sumW += 1./pow(err, 2); }
      }
      if (sumW>0.){ sum/=sumW; sumW = 1./sqrt(sumW); }
      res.getProfileY()->SetBinContent(iy, sum);
      res.getProfileY()->SetBinError(iy, sumW);
    }
    // Replace prof_z of res
    for (int iz=0; iz<=(int) finalZBinning.getNbins()+1; iz++){
      double sum=0, sumW=0;
      for (int ix=0; ix<nIbinsX; ix++){
        double val=hRebinProfileZList.at(ix)->GetBinContent(iz);
        double err=hRebinProfileZList.at(ix)->GetBinError(iz);
        if (err>0.){ val/=pow(err, 2); sum+=val; sumW += 1./pow(err, 2); }
      }
      if (sumW>0.){ sum/=sumW; sumW = 1./sqrt(sumW); }
      res.getProfileZ()->SetBinContent(iz, sum);
      res.getProfileZ()->SetBinError(iz, sumW);
    }
    for (int ix=1; ix<=nIbinsX; ix++){
      for (int iy=0; iy<=(int) finalYBinning.getNbins()+1; iy++){
        for (int iz=0; iz<=(int) finalZBinning.getNbins()+1; iz++){
          res.getHistogram()->SetBinContent(ix, iy, iz, hRebinList.at(ix)->GetBinContent(iy, iz));
          res.getHistogram()->SetBinError(ix, iy, iz, hRebinList.at(ix)->GetBinError(iy, iz));
        }
      }
    }
    res.rebin(&finalXBinning, &finalYBinning, &finalZBinning);

    // Cleanup
    for (auto*& h:hRebinList) delete h;
    for (auto*& h:hRebinProfileYList) delete h;
    for (auto*& h:hRebinProfileZList) delete h;

    res.setNameTitle(hname, htitle);
    TH3F* hres = new TH3F(*(res.getHistogram()));
    return hres;
    break;
  }
  case 1:
  {
    ExtendedBinning& intermediateYBinning = intermediateYBinningList.at(0);
    int nIbinsY = intermediateYBinning.getNbins();
    res.constructFromTree(tree, xvar, yvar, zvar, weight, &flag, &finestXBinning, &intermediateYBinning, &finestZBinning);
    TH3F*& hIntermediate=res.getHistogram();
    vector<TH2F*> hRebinList; hRebinList.reserve(nIbinsY);
    vector<TProfile*> hRebinProfileXList, hRebinProfileZList; hRebinProfileXList.reserve(nIbinsY); hRebinProfileZList.reserve(nIbinsY);
    for (int i=1; i<=nIbinsY; i++){
      ExtendedBinning& binningX = intermediateXBinningList.at(i-1);
      ExtendedBinning& binningZ = intermediateZBinningList.at(i-1);
      TH2F* hRebin = getHistogramSlice(hIntermediate, 2, 0, i, i, Form("hRebinIntermediate_SliceX%i", i));
      TProfile* hRebinProfX = new TProfile(*(res.getProfileX())); hRebinProfX->SetName(Form("%s_Slice%i", hRebinProfX->GetName(), i));
      TProfile* hRebinProfZ = new TProfile(*(res.getProfileZ())); hRebinProfZ->SetName(Form("%s_Slice%i", hRebinProfZ->GetName(), i));
      rebinHistogram(hRebin, binningZ, binningX); rebinHistogram(hRebin, finalZBinning, finalXBinning);
      rebinProfile(hRebinProfX, finalXBinning); rebinProfile(hRebinProfZ, finalZBinning);
      hRebinList.push_back(hRebin);
      hRebinProfileXList.push_back(hRebinProfX);
      hRebinProfileZList.push_back(hRebinProfZ);
    }
    res.constructFromTree(tree, xvar, yvar, zvar, weight, &flag, &finalXBinning, &intermediateYBinning, &finalZBinning);
    // Replace prof_x of res
    for (int ix=0; ix<=(int) finalXBinning.getNbins()+1; ix++){
      double sum=0, sumW=0;
      for (int iy=0; iy<nIbinsY; iy++){
        double val=hRebinProfileXList.at(iy)->GetBinContent(ix);
        double err=hRebinProfileXList.at(iy)->GetBinError(ix);
        if (err>0.){ val/=pow(err, 2); sum+=val; sumW += 1./pow(err, 2); }
      }
      if (sumW>0.){ sum/=sumW; sumW = 1./sqrt(sumW); }
      res.getProfileX()->SetBinContent(ix, sum);
      res.getProfileX()->SetBinError(ix, sumW);
    }
    // Replace prof_z of res
    for (int iz=0; iz<=(int) finalZBinning.getNbins()+1; iz++){
      double sum=0, sumW=0;
      for (int iy=0; iy<nIbinsY; iy++){
        double val=hRebinProfileZList.at(iy)->GetBinContent(iz);
        double err=hRebinProfileZList.at(iy)->GetBinError(iz);
        if (err>0.){ val/=pow(err, 2); sum+=val; sumW += 1./pow(err, 2); }
      }
      if (sumW>0.){ sum/=sumW; sumW = 1./sqrt(sumW); }
      res.getProfileZ()->SetBinContent(iz, sum);
      res.getProfileZ()->SetBinError(iz, sumW);
    }
    for (int iy=1; iy<=nIbinsY; iy++){
      for (int iz=0; iz<=(int) finalZBinning.getNbins()+1; iz++){
        for (int ix=0; ix<=(int) finalXBinning.getNbins()+1; ix++){
          res.getHistogram()->SetBinContent(ix, iy, iz, hRebinList.at(iy)->GetBinContent(iz, ix));
          res.getHistogram()->SetBinError(ix, iy, iz, hRebinList.at(iy)->GetBinError(iz, ix));
        }
      }
    }
    res.rebin(&finalXBinning, &finalYBinning, &finalZBinning);

    // Cleanup
    for (auto*& h:hRebinList) delete h;
    for (auto*& h:hRebinProfileXList) delete h;
    for (auto*& h:hRebinProfileZList) delete h;

    res.setNameTitle(hname, htitle);
    TH3F* hres = new TH3F(*(res.getHistogram()));
    return hres;
    break;
  }
  case 2:
  {
    ExtendedBinning& intermediateZBinning = intermediateZBinningList.at(0);
    int nIbinsZ = intermediateZBinning.getNbins();
    res.constructFromTree(tree, xvar, yvar, zvar, weight, &flag, &finestXBinning, &finestYBinning, &intermediateZBinning);
    TH3F*& hIntermediate=res.getHistogram();
    vector<TH2F*> hRebinList; hRebinList.reserve(nIbinsZ);
    vector<TProfile*> hRebinProfileXList, hRebinProfileYList; hRebinProfileXList.reserve(nIbinsZ); hRebinProfileYList.reserve(nIbinsZ);
    for (int i=1; i<=nIbinsZ; i++){
      ExtendedBinning& binningX = intermediateXBinningList.at(i-1);
      ExtendedBinning& binningY = intermediateYBinningList.at(i-1);
      TH2F* hRebin = getHistogramSlice(hIntermediate, 0, 1, i, i, Form("hRebinIntermediate_SliceX%i", i));
      TProfile* hRebinProfX = new TProfile(*(res.getProfileX())); hRebinProfX->SetName(Form("%s_Slice%i", hRebinProfX->GetName(), i));
      TProfile* hRebinProfY = new TProfile(*(res.getProfileY())); hRebinProfY->SetName(Form("%s_Slice%i", hRebinProfY->GetName(), i));
      rebinHistogram(hRebin, binningX, binningY); rebinHistogram(hRebin, finalXBinning, finalYBinning);
      rebinProfile(hRebinProfX, finalXBinning); rebinProfile(hRebinProfY, finalYBinning);
      hRebinList.push_back(hRebin);
      hRebinProfileXList.push_back(hRebinProfX);
      hRebinProfileYList.push_back(hRebinProfY);
    }
    res.constructFromTree(tree, xvar, yvar, zvar, weight, &flag, &finalXBinning, &finalYBinning, &intermediateZBinning);
    // Replace prof_x of res
    for (int ix=0; ix<=(int) finalXBinning.getNbins()+1; ix++){
      double sum=0, sumW=0;
      for (int iz=0; iz<nIbinsZ; iz++){
        double val=hRebinProfileXList.at(iz)->GetBinContent(ix);
        double err=hRebinProfileXList.at(iz)->GetBinError(ix);
        if (err>0.){ val/=pow(err, 2); sum+=val; sumW += 1./pow(err, 2); }
      }
      if (sumW>0.){ sum/=sumW; sumW = 1./sqrt(sumW); }
      res.getProfileX()->SetBinContent(ix, sum);
      res.getProfileX()->SetBinError(ix, sumW);
    }
    // Replace prof_y of res
    for (int iy=0; iy<=(int) finalYBinning.getNbins()+1; iy++){
      double sum=0, sumW=0;
      for (int iz=0; iz<nIbinsZ; iz++){
        double val=hRebinProfileYList.at(iz)->GetBinContent(iy);
        double err=hRebinProfileYList.at(iz)->GetBinError(iy);
        if (err>0.){ val/=pow(err, 2); sum+=val; sumW += 1./pow(err, 2); }
      }
      if (sumW>0.){ sum/=sumW; sumW = 1./sqrt(sumW); }
      res.getProfileY()->SetBinContent(iy, sum);
      res.getProfileY()->SetBinError(iy, sumW);
    }
    for (int iz=1; iz<=nIbinsZ; iz++){
      for (int iy=0; iy<=(int) finalYBinning.getNbins()+1; iy++){
        for (int ix=0; ix<=(int) finalXBinning.getNbins()+1; ix++){
          res.getHistogram()->SetBinContent(ix, iy, iz, hRebinList.at(iz)->GetBinContent(ix, iy));
          res.getHistogram()->SetBinError(ix, iy, iz, hRebinList.at(iz)->GetBinError(ix, iy));
        }
      }
    }
    res.rebin(&finalXBinning, &finalYBinning, &finalZBinning);

    // Cleanup
    for (auto*& h:hRebinList) delete h;
    for (auto*& h:hRebinProfileXList) delete h;
    for (auto*& h:hRebinProfileYList) delete h;

    res.setNameTitle(hname, htitle);
    TH3F* hres = new TH3F(*(res.getHistogram()));
    return hres;
    break;
  }
  default:
  {
    res.constructFromTree(tree, xvar, yvar, zvar, weight, &flag, &(intermediateXBinningList.at(0)), &(intermediateYBinningList.at(0)), &(intermediateZBinningList.at(0)));
    res.rebin(&finalXBinning, &finalYBinning, &finalZBinning);

    res.setNameTitle(hname, htitle);
    TH3F* hres = new TH3F(*(res.getHistogram()));
    return hres;
  }
  }
}


#endif
