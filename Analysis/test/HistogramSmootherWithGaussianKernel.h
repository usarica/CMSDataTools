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
  double sigmaXmult=1
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
    if (flag) reference.fill(xvar, weight);
  }
  TH1F* res = new TH1F(
    hname, htitle,
    finalXBinning.getNbins(), finalXBinning.getBinning()
  ); res->Sumw2();

  SimpleGaussian gausX(0, 1, SimpleGaussian::kHasLowHighRange, bX.getMin(), bX.getMax());

  MELAout << "getSmoothHistogram: Filling the actual histogram with the help of reference" << endl;
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    progressbar(ev, tree->GetEntries());

    if (!flag) continue;
    if (xvar<bX.getMin() || xvar>=bX.getMax()) continue;

    int ix=bX.getBin(xvar);
    double Neff=pow(reference.getBinSumW(ix), 2)/reference.getBinSumWsq(ix);
    double widthGlobalScale = 1./sqrt(Neff);
    double sX = bX.getBinWidth(ix)*widthGlobalScale*sigmaXmult;
    gausX.setMean(xvar); gausX.setSigma(sX);

    for (unsigned int i=0; i<bX.getNbins(); i++){
      if (sX==0. && (int) i!=ix) continue;
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
  double sigmaXmult=1, double sigmaYmult=1
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
    if (flag) reference.fill(xvar, yvar, weight);
  }
  TH2F* res = new TH2F(
    hname, htitle,
    finalXBinning.getNbins(), finalXBinning.getBinning(),
    finalYBinning.getNbins(), finalYBinning.getBinning()
  ); res->Sumw2();

  SimpleGaussian gausX(0, 1, SimpleGaussian::kHasLowHighRange, bX.getMin(), bX.getMax());
  SimpleGaussian gausY(0, 1, SimpleGaussian::kHasLowHighRange, bY.getMin(), bY.getMax());

  MELAout << "getSmoothHistogram: Filling the actual histogram with the help of reference" << endl;
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    progressbar(ev, tree->GetEntries());

    if (!flag) continue;
    if (xvar<bX.getMin() || xvar>=bX.getMax()) continue;
    if (yvar<bY.getMin() || yvar>=bY.getMax()) continue;

    int ix=bX.getBin(xvar);
    int iy=bY.getBin(yvar);
    double Neff=pow(reference.getBinSumW(ix, iy), 2)/reference.getBinSumWsq(ix, iy);
    double widthGlobalScale = 1./sqrt(Neff);
    double sX = bX.getBinWidth(ix)*widthGlobalScale*sigmaXmult;
    gausX.setMean(xvar); gausX.setSigma(sX);
    double sY = bY.getBinWidth(iy)*widthGlobalScale*sigmaYmult;
    gausY.setMean(yvar); gausY.setSigma(sY);
    assert(checkVarNanInf(sX) && checkVarNanInf(sY));

    for (unsigned int i=0; i<bX.getNbins(); i++){
      if (sX==0. && (int) i!=ix) continue;
      double fX = gausX.integralNorm(bX.getBinLowEdge(i), bX.getBinHighEdge(i));
      if (fX>1. || fX<0.){ MELAerr << "fX=" << fX << endl; continue; }
      if (fX==0.) continue;
      for (unsigned int j=0; j<bY.getNbins(); j++){
        if (sY==0. && (int) j!=iy) continue;
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
  double sigmaXmult=1, double sigmaYmult=1, double sigmaZmult=1
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
    if (flag) reference.fill(xvar, yvar, zvar, weight);
  }
  TH3F* res = new TH3F(
    hname, htitle,
    finalXBinning.getNbins(), finalXBinning.getBinning(),
    finalYBinning.getNbins(), finalYBinning.getBinning(),
    finalZBinning.getNbins(), finalZBinning.getBinning()
  ); res->Sumw2();

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

    int ix=bX.getBin(xvar);
    int iy=bY.getBin(yvar);
    int iz=bZ.getBin(zvar);
    double Neff=pow(reference.getBinSumW(ix, iy, iz), 2)/reference.getBinSumWsq(ix, iy, iz);
    double widthGlobalScale = 1./sqrt(Neff);
    double sX = bX.getBinWidth(ix)*widthGlobalScale*sigmaXmult;
    gausX.setMean(xvar); gausX.setSigma(sX);
    double sY = bY.getBinWidth(iy)*widthGlobalScale*sigmaYmult;
    gausY.setMean(yvar); gausY.setSigma(sY);
    double sZ = bZ.getBinWidth(iz)*widthGlobalScale*sigmaZmult;
    gausZ.setMean(zvar); gausZ.setSigma(sZ);
    assert(checkVarNanInf(sX) && checkVarNanInf(sY) && checkVarNanInf(sZ));

    for (unsigned int i=0; i<bX.getNbins(); i++){
      if (sX==0. && (int) i!=ix) continue;
      double fX = gausX.integralNorm(bX.getBinLowEdge(i), bX.getBinHighEdge(i));
      if (fX>1. || fX<0.){ MELAerr << "fX=" << fX << endl; continue; }
      if (fX==0.) continue;
      for (unsigned int j=0; j<bY.getNbins(); j++){
        if (sY==0. && (int) j!=iy) continue;
        double fY = gausY.integralNorm(bY.getBinLowEdge(j), bY.getBinHighEdge(j));
        if (fY>1. || fY<0.){ MELAerr << "fY=" << fY << endl; continue; }
        if (fY==0.) continue;
        for (unsigned int k=0; k<bZ.getNbins(); k++){
          if (sZ==0. && (int) k!=iz) continue;
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
