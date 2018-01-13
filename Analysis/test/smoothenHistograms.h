#ifndef SMOOTHENHISTOGRAMS_H
#define SMOOTHENHISTOGRAMS_H

#include "common_includes.h"


void smoothenHistograms(TH1F* hist){
  float integral_old = hist->Integral(1, hist->GetNbinsX(), "width");
  HelperFunctions::regularizeHistogram(hist, 1000, 0, 0.5, 0);
  float integral_new = hist->Integral(1, hist->GetNbinsX(), "width");
  double scale = 1.;
  if (integral_old!=0.) scale = integral_new/integral_old;
  for (int ix=1; ix<=histo->GetNbinsX(); ix++) hist->SetBinContent(ix, scale*hist->GetBinContent(ix));
}

void smoothenHistograms(TH2F* hist, bool normX, bool normY, std::vector<unsigned int>* iorder){
  if (normX && normY) return;
  else if (!normX && !normY){
    std::vector<unsigned int> axis_order;
    if (iorder && iorder->size()==2) axis_order=*iorder;
    else{
      axis_order.reserve(2);
      axis_order.push_back(0);
      axis_order.push_back(1);
    }
    float integral_old = hist->Integral(1, hist->GetNbinsX(), 1, hist->GetNbinsY(), "width");
    HelperFunctions::regularizeHistogram(hist, 1000, 0, 0.5, axis_order.at(0));
    HelperFunctions::regularizeHistogram(hist, 1000, 0, 0.5, axis_order.at(1));
    float integral_new = hist->Integral(1, hist->GetNbinsX(), 1, hist->GetNbinsY(), "width");
    double scale = 1.;
    if (integral_old!=0.) scale = integral_new/integral_old;
    for (int ix=1; ix<=histo->GetNbinsX(); ix++){ for (int iy=1; iy<=histo->GetNbinsY(); iy++) hist->SetBinContent(ix, iy, scale*hist->GetBinContent(ix, iy)); }
  }
  else if (normX){
    // Smoothen
    HelperFunctions::regularizeHistogram(hist, 1000, 0, 0.5, 1);
    conditionalizeHistogram(hist, 0);
  }
  else{
    // Smoothen
    HelperFunctions::regularizeHistogram(hist, 1000, 0, 0.5, 0);
    conditionalizeHistogram(hist, 1);
  }
}

void smoothenHistograms(TH2F* hist, bool normX, bool normY, bool normZ, std::vector<unsigned int>* iorder){
  if (
    (normX && normY)
    || (normX && normZ)
    || (normY && normZ)
    ) return;
  else if (!normX && !normY && !normZ){
    std::vector<unsigned int> axis_order;
    if (iorder && iorder->size()==3) axis_order=*iorder;
    else{
      axis_order.reserve(3);
      axis_order.push_back(0);
      axis_order.push_back(1);
      axis_order.push_back(2);
    }
    float integral_old = hist->Integral(1, hist->GetNbinsX(), 1, hist->GetNbinsY(), 1, hist->GetNbinsZ(), "width");
    HelperFunctions::regularizeHistogram(hist, 1000, 0, 0.5, axis_order.at(0));
    HelperFunctions::regularizeHistogram(hist, 1000, 0, 0.5, axis_order.at(1));
    HelperFunctions::regularizeHistogram(hist, 1000, 0, 0.5, axis_order.at(2));
    float integral_new = hist->Integral(1, hist->GetNbinsX(), 1, hist->GetNbinsY(), 1, hist->GetNbinsZ(), "width");
    double scale = 1.;
    if (integral_old!=0.) scale = integral_new/integral_old;
    for (int ix=1; ix<=histo->GetNbinsX(); ix++){ for (int iy=1; iy<=histo->GetNbinsY(); iy++){ for (int iz=1; iz<=histo->GetNbinsZ(); iz++) hist->SetBinContent(ix, iy, iz, scale*hist->GetBinContent(ix, iy, iz)); } }
  }
  else if (normX){
    std::vector<unsigned int> axis_order;
    if (iorder && iorder->size()==2 && iorder->at(0)!=0 && iorder->at(1)!=0) axis_order=*iorder;
    else{
      axis_order.reserve(2);
      axis_order.push_back(1);
      axis_order.push_back(2);
    }
    // Smoothen
    HelperFunctions::regularizeHistogram(hist, 1000, 0, 0.5, axis_order.at(0));
    HelperFunctions::regularizeHistogram(hist, 1000, 0, 0.5, axis_order.at(1));
    conditionalizeHistogram(hist, 0);
  }
  else if (normY){
    std::vector<unsigned int> axis_order;
    if (iorder && iorder->size()==2 && iorder->at(0)!=1 && iorder->at(1)!=1) axis_order=*iorder;
    else{
      axis_order.reserve(2);
      axis_order.push_back(0);
      axis_order.push_back(2);
    }
    // Smoothen
    HelperFunctions::regularizeHistogram(hist, 1000, 0, 0.5, axis_order.at(0));
    HelperFunctions::regularizeHistogram(hist, 1000, 0, 0.5, axis_order.at(1));
    conditionalizeHistogram(hist, 1);
  }
  else{
    std::vector<unsigned int> axis_order;
    if (iorder && iorder->size()==2 && iorder->at(0)!=2 && iorder->at(1)!=2) axis_order=*iorder;
    else{
      axis_order.reserve(2);
      axis_order.push_back(0);
      axis_order.push_back(1);
    }
    // Smoothen
    HelperFunctions::regularizeHistogram(hist, 1000, 0, 0.5, axis_order.at(0));
    HelperFunctions::regularizeHistogram(hist, 1000, 0, 0.5, axis_order.at(1));
    conditionalizeHistogram(hist, 2);
  }
}

#endif
