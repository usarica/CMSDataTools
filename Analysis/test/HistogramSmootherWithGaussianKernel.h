#ifndef HISTOGRAMSMOOTHERWITHGAUSSIANKERNEL_H
#define HISTOGRAMSMOOTHERWITHGAUSSIANKERNEL_H

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
#include "FunctionHelpers.h"
#include "SampleHelpers.h"
#include "DiscriminantClasses.h"
#include "CategorizationHelpers.h"
#include "ExtendedHistogram_1D.h"
#include "ExtendedHistogram_2D.h"
#include "ExtendedHistogram_3D.h"
#include "ExtendedProfileHistogram.h"
#include "TemplateHelpers.h"
#include "MELAStreamHelpers.hh"
#include "Mela.h"


#ifndef GAUSSIANWIDTHPRECISION
#define GAUSSIANWIDTHPRECISION 5.
#endif


using namespace std;
using namespace HelperFunctions;
using namespace FunctionHelpers;
using namespace SampleHelpers;
using namespace DiscriminantClasses;
using namespace TemplateHelpers;
using namespace MELAStreamHelpers;


ExtendedBinning getIntermediateBinning(ExtendedBinning const& binning){
  ExtendedBinning res(binning);
  TString namelower=res.getLabel(); namelower.ToLower();
  if (!namelower.Contains("mass") && !namelower.Contains("pt")) return res;
  res.addBinBoundary(binning.getBinLowEdge(0)-binning.getBinWidth(0));
  res.addBinBoundary(binning.getBinHighEdge(binning.getNbins()-1)+binning.getBinWidth(binning.getNbins()-1));
  MELAout
    << "getIntermediateBinning: Extended binning " << res.getLabel()
    << " [ " << res.getMin() << ", " << res.getMax() << " ]"
    << " (nbins = " << res.getNbins() << ")"
    << " is created." << endl;
  return res;
}


TH1F* getSmoothHistogram(
  TString const hname, TString const htitle, ExtendedBinning const& finalXBinning,
  TTree* tree, float& xvar, float& weight, bool& flag,
  double sigmaXmult=1,
  TH1F** hRawPtr=nullptr
){
  assert(tree && finalXBinning.isValid());

  ExtendedBinning bX=getIntermediateBinning(finalXBinning);
  bool sameXbins=(bX.getNbins()==finalXBinning.getNbins());

  // Construct fine histogram to determine intermediate binning
  MELAout << "getSmoothHistogram: Filling the reference ExtendedProfileHistogram" << endl;
  ExtendedProfileHistogram reference(bX);
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    progressbar(ev, tree->GetEntries());
    if (flag) reference.fill(xvar, fabs(weight));
  }
  TH1F* res = new TH1F(
    hname, htitle,
    finalXBinning.getNbins(), finalXBinning.getBinning()
  );
  res->Sumw2();
  res->GetXaxis()->SetTitle(finalXBinning.getLabel());
  if (hRawPtr){
    *hRawPtr = new TH1F(
      hname+"_raw", htitle,
      finalXBinning.getNbins(), finalXBinning.getBinning()
    );
    (*hRawPtr)->Sumw2();
    (*hRawPtr)->GetXaxis()->SetTitle(finalXBinning.getLabel());
  }

  SimpleGaussian gausX(0, 1, SimpleGaussian::kHasLowHighRange, bX.getMin(), bX.getMax());

  MELAout << "getSmoothHistogram: Filling the actual histogram with the help of reference" << endl;
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    progressbar(ev, tree->GetEntries());

    if (!flag) continue;
    if (xvar<bX.getMin() || xvar>=bX.getMax()) continue;
    if (hRawPtr) (*hRawPtr)->Fill(xvar, weight);

    int ix=bX.getBin(xvar);
    double Neff=pow(reference.getBinSumW(ix), 2)/reference.getBinSumWsq(ix);
    double widthGlobalScale = 1./sqrt(Neff);

    double sX = bX.getBinWidth(ix)*widthGlobalScale*sigmaXmult;
    if (std::min(fabs(xvar-bX.getBinLowEdge(ix)), fabs(xvar-bX.getBinHighEdge(ix)))>=sX*GAUSSIANWIDTHPRECISION) sX=0.;
    unsigned int ibegin, iend;
    if (sX!=0. || ix<0){ ibegin=0; iend=bX.getNbins(); }
    else{ ibegin=ix; iend=ix; }
    gausX.setMean(xvar); gausX.setSigma(sX);

    assert(checkVarNanInf(sX));

    for (unsigned int i=ibegin; i<iend; i++){
      double fX = gausX.integralNorm(bX.getBinLowEdge(i), bX.getBinHighEdge(i));
      if (fX>1. || fX<0.){ MELAerr << "fX=" << fX << endl; continue; }
      if (fX==0.) continue;

      double w=fX*weight;
      int ii=i;
      if (sameXbins) ii++;
      {
        double bincontent = res->GetBinContent(ii);
        double binerror = res->GetBinError(ii);
        res->SetBinContent(ii, bincontent+w);
        res->SetBinError(ii, sqrt(pow(bincontent, 2)+pow(w, 2)));
      }
    }

  } // End loop over tree
  return res;
}


TH2F* getSmoothHistogram(
  TString const hname, TString const htitle, ExtendedBinning const& finalXBinning, ExtendedBinning const& finalYBinning,
  TTree* tree, float& xvar, float& yvar, float& weight, bool& flag,
  double sigmaXmult=1, double sigmaYmult=1,
  TH2F** hRawPtr=nullptr
){
  assert(tree && finalXBinning.isValid() && finalYBinning.isValid());

  ExtendedBinning bX=getIntermediateBinning(finalXBinning);
  ExtendedBinning bY=getIntermediateBinning(finalYBinning);
  bool sameXbins=(bX.getNbins()==finalXBinning.getNbins());
  bool sameYbins=(bY.getNbins()==finalYBinning.getNbins());

  // Construct fine histogram to determine intermediate binning
  MELAout << "getSmoothHistogram: Filling the reference ExtendedProfileHistogram" << endl;
  ExtendedProfileHistogram reference(bX, bY);
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    progressbar(ev, tree->GetEntries());
    if (flag) reference.fill(xvar, yvar, fabs(weight));
  }
  TH2F* res = new TH2F(
    hname, htitle,
    finalXBinning.getNbins(), finalXBinning.getBinning(),
    finalYBinning.getNbins(), finalYBinning.getBinning()
  );
  res->Sumw2();
  res->GetXaxis()->SetTitle(finalXBinning.getLabel());
  res->GetYaxis()->SetTitle(finalYBinning.getLabel());
  if (hRawPtr){
    *hRawPtr = new TH2F(
      hname+"_raw", htitle,
      finalXBinning.getNbins(), finalXBinning.getBinning(),
      finalYBinning.getNbins(), finalYBinning.getBinning()
    );
    (*hRawPtr)->Sumw2();
    (*hRawPtr)->GetXaxis()->SetTitle(finalXBinning.getLabel());
    (*hRawPtr)->GetYaxis()->SetTitle(finalYBinning.getLabel());
  }

  SimpleGaussian gausX(0, 1, SimpleGaussian::kHasLowHighRange, bX.getMin(), bX.getMax());
  SimpleGaussian gausY(0, 1, SimpleGaussian::kHasLowHighRange, bY.getMin(), bY.getMax());

  MELAout << "getSmoothHistogram: Filling the actual histogram with the help of reference" << endl;
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    progressbar(ev, tree->GetEntries());

    if (!flag) continue;
    if (xvar<bX.getMin() || xvar>=bX.getMax()) continue;
    if (yvar<bY.getMin() || yvar>=bY.getMax()) continue;
    if (hRawPtr) (*hRawPtr)->Fill(xvar, yvar, weight);

    int ix=bX.getBin(xvar);
    int iy=bY.getBin(yvar);
    double Neff=pow(reference.getBinSumW(ix, iy), 2)/reference.getBinSumWsq(ix, iy);
    double widthGlobalScale = 1./sqrt(Neff);

    double sX = bX.getBinWidth(ix)*widthGlobalScale*sigmaXmult;
    if (std::min(fabs(xvar-bX.getBinLowEdge(ix)), fabs(xvar-bX.getBinHighEdge(ix)))>=sX*GAUSSIANWIDTHPRECISION) sX=0.;
    unsigned int ibegin, iend;
    if (sX!=0. || ix<0){ ibegin=0; iend=bX.getNbins(); }
    else{ ibegin=ix; iend=ix; }
    gausX.setMean(xvar); gausX.setSigma(sX);

    double sY = bY.getBinWidth(iy)*widthGlobalScale*sigmaYmult;
    if (std::min(fabs(yvar-bY.getBinLowEdge(iy)), fabs(yvar-bY.getBinHighEdge(iy)))>=sY*GAUSSIANWIDTHPRECISION) sY=0.;
    unsigned int jbegin, jend;
    if (sY!=0. || iy<0){ jbegin=0; jend=bY.getNbins(); }
    else{ jbegin=iy; jend=iy; }
    gausY.setMean(yvar); gausY.setSigma(sY);

    assert(checkVarNanInf(sX) && checkVarNanInf(sY));

    for (unsigned int i=ibegin; i<iend; i++){
      double fX = gausX.integralNorm(bX.getBinLowEdge(i), bX.getBinHighEdge(i));
      if (fX>1. || fX<0.){ MELAerr << "fX=" << fX << endl; continue; }
      if (fX==0.) continue;
      for (unsigned int j=jbegin; j<jend; j++){
        double fY = gausY.integralNorm(bY.getBinLowEdge(j), bY.getBinHighEdge(j));
        if (fY>1. || fY<0.){ MELAerr << "fY=" << fY << endl; continue; }
        if (fY==0.) continue;

        double fprod=fX*fY;
        double w=fprod*weight;
        int ii=i, jj=j;
        if (sameXbins) ii++;
        if (sameYbins) jj++;
        {
          double bincontent = res->GetBinContent(ii, jj);
          double binerror = res->GetBinError(ii, jj);
          res->SetBinContent(ii, jj, bincontent+w);
          res->SetBinError(ii, jj, sqrt(pow(bincontent, 2)+pow(w, 2)));
        }
      }
    }

  } // End loop over tree
  return res;
}


TH3F* getSmoothHistogram(
  TString const hname, TString const htitle, ExtendedBinning const& finalXBinning, ExtendedBinning const& finalYBinning, ExtendedBinning const& finalZBinning,
  TTree* tree, float& xvar, float& yvar, float& zvar, float& weight, bool& flag,
  double sigmaXmult=1, double sigmaYmult=1, double sigmaZmult=1,
  TH3F** hRawPtr=nullptr
){
  assert(tree && finalXBinning.isValid() && finalYBinning.isValid() && finalZBinning.isValid());

  ExtendedBinning bX=getIntermediateBinning(finalXBinning);
  ExtendedBinning bY=getIntermediateBinning(finalYBinning);
  ExtendedBinning bZ=getIntermediateBinning(finalZBinning);
  bool sameXbins=(bX.getNbins()==finalXBinning.getNbins());
  bool sameYbins=(bY.getNbins()==finalYBinning.getNbins());
  bool sameZbins=(bZ.getNbins()==finalZBinning.getNbins());

  // Construct fine histogram to determine intermediate binning
  MELAout << "getSmoothHistogram: Filling the reference ExtendedProfileHistogram" << endl;
  ExtendedProfileHistogram reference(bX, bY, bZ);
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    progressbar(ev, tree->GetEntries());
    if (flag) reference.fill(xvar, yvar, zvar, fabs(weight));
  }
  TH3F* res = new TH3F(
    hname, htitle,
    finalXBinning.getNbins(), finalXBinning.getBinning(),
    finalYBinning.getNbins(), finalYBinning.getBinning(),
    finalZBinning.getNbins(), finalZBinning.getBinning()
  );
  res->Sumw2();
  res->GetXaxis()->SetTitle(finalXBinning.getLabel());
  res->GetYaxis()->SetTitle(finalYBinning.getLabel());
  res->GetZaxis()->SetTitle(finalZBinning.getLabel());
  if (hRawPtr){
    *hRawPtr = new TH3F(
      hname+"_raw", htitle,
      finalXBinning.getNbins(), finalXBinning.getBinning(),
      finalYBinning.getNbins(), finalYBinning.getBinning(),
      finalZBinning.getNbins(), finalZBinning.getBinning()
    );
    (*hRawPtr)->Sumw2();
    (*hRawPtr)->GetXaxis()->SetTitle(finalXBinning.getLabel());
    (*hRawPtr)->GetYaxis()->SetTitle(finalYBinning.getLabel());
    (*hRawPtr)->GetZaxis()->SetTitle(finalZBinning.getLabel());
  }

  SimpleGaussian gausX(0, 1, SimpleGaussian::kHasLowHighRange, bX.getMin(), bX.getMax());
  SimpleGaussian gausY(0, 1, SimpleGaussian::kHasLowHighRange, bY.getMin(), bY.getMax());
  SimpleGaussian gausZ(0, 1, SimpleGaussian::kHasLowHighRange, bZ.getMin(), bZ.getMax());

  MELAout << "getSmoothHistogram: Filling the actual histogram with the help of reference" << endl;
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    progressbar(ev, tree->GetEntries());

    if (!flag) continue;
    if (xvar<bX.getMin() || xvar>=bX.getMax()) continue;
    if (yvar<bY.getMin() || yvar>=bY.getMax()) continue;
    if (zvar<bZ.getMin() || zvar>=bZ.getMax()) continue;
    if (hRawPtr) (*hRawPtr)->Fill(xvar, yvar, zvar, weight);

    int ix=bX.getBin(xvar);
    int iy=bY.getBin(yvar);
    int iz=bZ.getBin(zvar);
    double Neff=pow(reference.getBinSumW(ix, iy, iz), 2)/reference.getBinSumWsq(ix, iy, iz);
    double widthGlobalScale = 1./sqrt(Neff);

    double sX = bX.getBinWidth(ix)*widthGlobalScale*sigmaXmult;
    if (std::min(fabs(xvar-bX.getBinLowEdge(ix)), fabs(xvar-bX.getBinHighEdge(ix)))>=sX*GAUSSIANWIDTHPRECISION) sX=0.;
    unsigned int ibegin, iend;
    if (sX!=0. || ix<0){ ibegin=0; iend=bX.getNbins(); }
    else{ ibegin=ix; iend=ix; }
    gausX.setMean(xvar); gausX.setSigma(sX);

    double sY = bY.getBinWidth(iy)*widthGlobalScale*sigmaYmult;
    if (std::min(fabs(yvar-bY.getBinLowEdge(iy)), fabs(yvar-bY.getBinHighEdge(iy)))>=sY*GAUSSIANWIDTHPRECISION) sY=0.;
    unsigned int jbegin, jend;
    if (sY!=0. || iy<0){ jbegin=0; jend=bY.getNbins(); }
    else{ jbegin=iy; jend=iy; }
    gausY.setMean(yvar); gausY.setSigma(sY);

    double sZ = bZ.getBinWidth(iz)*widthGlobalScale*sigmaZmult;
    if (std::min(fabs(zvar-bZ.getBinLowEdge(iz)), fabs(zvar-bZ.getBinHighEdge(iz)))>=sZ*GAUSSIANWIDTHPRECISION) sZ=0.;
    unsigned int kbegin, kend;
    if (sZ!=0. || iz<0){ kbegin=0; kend=bZ.getNbins(); }
    else{ kbegin=iz; kend=iz; }
    gausZ.setMean(zvar); gausZ.setSigma(sZ);

    assert(checkVarNanInf(sX) && checkVarNanInf(sY) && checkVarNanInf(sZ));

    for (unsigned int i=ibegin; i<iend; i++){
      double fX = gausX.integralNorm(bX.getBinLowEdge(i), bX.getBinHighEdge(i));
      if (fX>1. || fX<0.){ MELAerr << "fX=" << fX << endl; continue; }
      if (fX==0.) continue;
      for (unsigned int j=jbegin; j<jend; j++){
        double fY = gausY.integralNorm(bY.getBinLowEdge(j), bY.getBinHighEdge(j));
        if (fY>1. || fY<0.){ MELAerr << "fY=" << fY << endl; continue; }
        if (fY==0.) continue;
        for (unsigned int k=kbegin; k<kend; k++){
          double fZ = gausZ.integralNorm(bZ.getBinLowEdge(k), bZ.getBinHighEdge(k));
          if (fZ>1. || fZ<0.){ MELAerr << "fZ=" << fZ << endl; continue; }
          if (fZ==0.) continue;

          double fprod=fX*fY*fZ;
          double w=fprod*weight;
          int ii=i, jj=j, kk=k;
          if (sameXbins) ii++;
          if (sameYbins) jj++;
          if (sameZbins) kk++;
          {
            double bincontent = res->GetBinContent(ii, jj, kk);
            double binerror = res->GetBinError(ii, jj, kk);
            res->SetBinContent(ii, jj, kk, bincontent+w);
            res->SetBinError(ii, jj, kk, sqrt(pow(bincontent, 2)+pow(w, 2)));
          }
        }
      }
    }

  } // End loop over tree
  return res;
}


#endif
