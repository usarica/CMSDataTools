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


class TreeHistogramAssociation_1D{
public:
  TString const hname;
  TString const htitle;

  TTree* const tree;
  float& xvar;
  float& weight;
  bool& flag;

  TreeHistogramAssociation_1D(TString const hname_, TString const htitle_, TTree* tree_, float& xvar_, float& weight_, bool& flag_);
};
TreeHistogramAssociation_1D::TreeHistogramAssociation_1D(TString const hname_, TString const htitle_, TTree* tree_, float& xvar_, float& weight_, bool& flag_) :
  hname(hname_), htitle(htitle_),
  tree(tree_), xvar(xvar_),
  weight(weight_), flag(flag_)
{
  assert(tree);
}
class TreeHistogramAssociation_2D : public TreeHistogramAssociation_1D{
public:
  float& yvar;

  TreeHistogramAssociation_2D(TString const hname_, TString const htitle_, TTree* tree_, float& xvar_, float& yvar_, float& weight_, bool& flag_);
};
TreeHistogramAssociation_2D::TreeHistogramAssociation_2D(TString const hname_, TString const htitle_, TTree* tree_, float& xvar_, float& yvar_, float& weight_, bool& flag_) :
  TreeHistogramAssociation_1D(hname_, htitle_, tree_, xvar_, weight_, flag_),
  yvar(yvar_)
{}
class TreeHistogramAssociation_3D : public TreeHistogramAssociation_2D{
public:
  float& zvar;

  TreeHistogramAssociation_3D(TString const hname_, TString const htitle_, TTree* tree_, float& xvar_, float& yvar_, float& zvar_, float& weight_, bool& flag_);
};
TreeHistogramAssociation_3D::TreeHistogramAssociation_3D(TString const hname_, TString const htitle_, TTree* tree_, float& xvar_, float& yvar_, float& zvar_, float& weight_, bool& flag_) :
  TreeHistogramAssociation_2D(hname_, htitle_, tree_, xvar_, yvar_, weight_, flag_),
  zvar(zvar_)
{}


ExtendedBinning getIntermediateBinning(ExtendedBinning const& binning){
  ExtendedBinning res(binning);
  TString namelower=res.getLabel(); namelower.ToLower();
  if ((!namelower.Contains("mass") && !namelower.Contains("pt")) || binning.getNbins()<4) return res;
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
  ExtendedProfileHistogram reference(bX, false);
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
    else{ ibegin=ix; iend=ix+1; }
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
  ExtendedProfileHistogram reference(bX, bY, false);
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
    else{ ibegin=ix; iend=ix+1; }
    gausX.setMean(xvar); gausX.setSigma(sX);

    double sY = bY.getBinWidth(iy)*widthGlobalScale*sigmaYmult;
    if (std::min(fabs(yvar-bY.getBinLowEdge(iy)), fabs(yvar-bY.getBinHighEdge(iy)))>=sY*GAUSSIANWIDTHPRECISION) sY=0.;
    unsigned int jbegin, jend;
    if (sY!=0. || iy<0){ jbegin=0; jend=bY.getNbins(); }
    else{ jbegin=iy; jend=iy+1; }
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
  ExtendedProfileHistogram reference(bX, bY, bZ, false);
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    progressbar(ev, tree->GetEntries());
    if (flag) reference.fill(xvar, yvar, zvar, fabs(weight));
  }

  ExtendedProfileHistogram extres(bX, bY, bZ, false); // For the res histogram
  std::vector<std::vector<std::vector<double>>>& extres_sumW=extres.getSumWContainer();
  std::vector<std::vector<std::vector<double>>>& extres_sumWsq=extres.getSumWsqContainer();
  ExtendedProfileHistogram extresraw(bX, bY, bZ, false); // For the resraw histogram, if exists

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
    if (hRawPtr) extresraw.fill(xvar, yvar, zvar, weight);

    int ix=bX.getBin(xvar);
    int iy=bY.getBin(yvar);
    int iz=bZ.getBin(zvar);
    double Neff=pow(reference.getBinSumW(ix, iy, iz), 2)/reference.getBinSumWsq(ix, iy, iz);
    double widthGlobalScale = 1./sqrt(Neff);

    double sX = bX.getBinWidth(ix)*widthGlobalScale*sigmaXmult;
    if (std::min(fabs(xvar-bX.getBinLowEdge(ix)), fabs(xvar-bX.getBinHighEdge(ix)))>=sX*GAUSSIANWIDTHPRECISION) sX=0.;
    unsigned int ibegin, iend;
    if (sX!=0. || ix<0){ ibegin=0; iend=bX.getNbins(); }
    else{ ibegin=ix; iend=ix+1; }
    gausX.setMean(xvar); gausX.setSigma(sX);

    double sY = bY.getBinWidth(iy)*widthGlobalScale*sigmaYmult;
    if (std::min(fabs(yvar-bY.getBinLowEdge(iy)), fabs(yvar-bY.getBinHighEdge(iy)))>=sY*GAUSSIANWIDTHPRECISION) sY=0.;
    unsigned int jbegin, jend;
    if (sY!=0. || iy<0){ jbegin=0; jend=bY.getNbins(); }
    else{ jbegin=iy; jend=iy+1; }
    gausY.setMean(yvar); gausY.setSigma(sY);

    double sZ = bZ.getBinWidth(iz)*widthGlobalScale*sigmaZmult;
    if (std::min(fabs(zvar-bZ.getBinLowEdge(iz)), fabs(zvar-bZ.getBinHighEdge(iz)))>=sZ*GAUSSIANWIDTHPRECISION) sZ=0.;
    unsigned int kbegin, kend;
    if (sZ!=0. || iz<0){ kbegin=0; kend=bZ.getNbins(); }
    else{ kbegin=iz; kend=iz+1; }
    gausZ.setMean(zvar); gausZ.setSigma(sZ);

    assert(checkVarNanInf(sX) && checkVarNanInf(sY) && checkVarNanInf(sZ));

    { // 3D is slower than 1D and 2D, so fill manually
      unsigned int i=ibegin;
      auto it_i = extres_sumW.begin()+ibegin;
      auto itsq_i = extres_sumWsq.begin()+ibegin;
      while (i<iend){
        double fX = gausX.integralNorm(bX.getBinLowEdge(i), bX.getBinHighEdge(i));
        bool doProceedX=true;
        if (fX>1. || fX<0.){ MELAerr << "fX=" << fX << endl; doProceedX=false; }
        doProceedX &= (fX!=0.);

        if (doProceedX){
          unsigned int j=jbegin;
          auto it_j = it_i->begin()+jbegin;
          auto itsq_j = itsq_i->begin()+jbegin;
          while (j<jend){
            double fY = gausY.integralNorm(bY.getBinLowEdge(j), bY.getBinHighEdge(j));
            bool doProceedY=true;
            if (fY>1. || fY<0.){ MELAerr << "fY=" << fY << endl; doProceedY=false; }
            doProceedY &= (fY!=0.);

            if (doProceedY){
              unsigned int k=kbegin;
              auto it_k = it_j->begin()+kbegin;
              auto itsq_k = itsq_j->begin()+kbegin;
              while (k<kend){
                double fZ = gausZ.integralNorm(bZ.getBinLowEdge(k), bZ.getBinHighEdge(k));
                bool doProceedZ=true;
                if (fZ>1. || fZ<0.){ MELAerr << "fZ=" << fZ << endl; doProceedZ=false; }
                doProceedZ &= (fZ!=0.);

                if (doProceedZ){
                  double fprod=fX*fY*fZ;
                  double w=fprod*weight;
                  *(it_k) += w;
                  *(itsq_k) += pow(w, 2);
                }
                k++; it_k++; itsq_k++;
              }
            }
            j++; it_j++; itsq_j++;
          }
        }
        i++; it_i++; itsq_i++;
      }
    } // End scope of i and iterators

  } // End loop over tree

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
  for (unsigned int i=0; i<bX.getNbins(); i++){
    unsigned int ii=i;
    if (sameXbins) ii++;
    for (unsigned int j=0; j<bY.getNbins(); j++){
      unsigned int jj=j;
      if (sameYbins) jj++;
      for (unsigned int k=0; k<bZ.getNbins(); k++){
        unsigned int kk=k;
        if (sameZbins) kk++;
        res->SetBinContent(ii, jj, kk, extres.getBinSumW(i, j, k));
        res->SetBinError(ii, jj, kk, sqrt(extres.getBinSumWsq(i, j, k)));
        if (hRawPtr){
          (*hRawPtr)->SetBinContent(ii, jj, kk, extresraw.getBinSumW(i, j, k));
          (*hRawPtr)->SetBinError(ii, jj, kk, sqrt(extresraw.getBinSumWsq(i, j, k)));
        }
      }
    }
  }

  return res;
}


#endif
