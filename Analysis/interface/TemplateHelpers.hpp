#ifndef TEMPLATEHELPERS_HPP
#define TEMPLATEHELPERS_HPP

#include "TemplateHelpers.h"
#include "HelperFunctions.h"
#include "Samples.h"
#include <cassert>
#include "TAxis.h"


/*********************/
/* General functions */
/*********************/
template<typename T> void TemplateHelpers::setTemplateAxisLabels(T* histo){
  if (!histo) return;
  bool is1D = dynamic_cast<TH1F*>(histo)!=nullptr;
  bool is2D = dynamic_cast<TH2F*>(histo)!=nullptr;
  bool is3D = dynamic_cast<TH3F*>(histo)!=nullptr;
  std::vector<TAxis*> axes;
  if (is1D || is2D || is3D) axes.push_back(histo->GetXaxis());
  if (is2D || is3D) axes.push_back(histo->GetYaxis());
  if (is3D) axes.push_back(histo->GetZaxis());
  for (TAxis* axis:axes){
    TString oldlabel=axis->GetTitle();
    TString newlabel=DiscriminantClasses::getKDLabel(oldlabel);
    if (newlabel=="" && oldlabel=="ZZMass") newlabel = "m_{4l} (GeV)";
    else if (newlabel=="") newlabel=oldlabel;
    axis->SetTitle(newlabel);
  }
}
template void TemplateHelpers::setTemplateAxisLabels<TH1F>(TH1F* histo);
template void TemplateHelpers::setTemplateAxisLabels<TH2F>(TH2F* histo);
template void TemplateHelpers::setTemplateAxisLabels<TH3F>(TH3F* histo);


/*****************/
/* QQ background */
/*****************/
template<typename T> void TemplateHelpers::recombineQQBkgHistogramsToTemplates(std::vector<T*>& vals){
  assert(castQQBkgHypothesisTypeToInt(nQQBkgTypes)==castQQBkgTemplateTypeToInt(nQQBkgTplTypes));
  if ((int) vals.size()!=castQQBkgHypothesisTypeToInt(nQQBkgTypes)) return;
  for (T*& hh:vals){
    HelperFunctions::wipeOverUnderFlows<T>(hh);
    HelperFunctions::divideBinWidth<T>(hh);
    hh->Scale(xsecScale);
    hh->SetTitle(getQQBkgProcessLabel());
    setTemplateAxisLabels<T>(hh);
  }
}
template void TemplateHelpers::recombineQQBkgHistogramsToTemplates<TH1F>(std::vector<TH1F*>& vals);
template void TemplateHelpers::recombineQQBkgHistogramsToTemplates<TH2F>(std::vector<TH2F*>& vals);
template void TemplateHelpers::recombineQQBkgHistogramsToTemplates<TH3F>(std::vector<TH3F*>& vals);


#endif
