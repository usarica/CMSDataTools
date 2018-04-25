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
template<typename T> void TemplateHelpers::setTemplateAxisLabels(T* histo, const CategorizationHelpers::Category category){
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
    if (newlabel==""){
      if (oldlabel=="ZZMass") newlabel = "m_{4l} (GeV)";
      else newlabel=oldlabel;
    }
    else if (newlabel.Contains("[category]")){
      TString strReplacement = CategorizationHelpers::getCategoryLabelForKDs(category);
      HelperFunctions::replaceString<TString, const TString>(newlabel, "[category]", strReplacement);
    }
    axis->SetTitle(newlabel);
  }
}
template void TemplateHelpers::setTemplateAxisLabels<TH1F>(TH1F* histo, const CategorizationHelpers::Category category);
template void TemplateHelpers::setTemplateAxisLabels<TH2F>(TH2F* histo, const CategorizationHelpers::Category category);
template void TemplateHelpers::setTemplateAxisLabels<TH3F>(TH3F* histo, const CategorizationHelpers::Category category);

template<typename T> void TemplateHelpers::doTemplatePostprocessing(T* tpl, const CategorizationHelpers::Category category, bool isMC){
  HelperFunctions::wipeOverUnderFlows(tpl);
  HelperFunctions::divideBinWidth(tpl);
  TemplateHelpers::setTemplateAxisLabels(tpl, category);
  if (isMC){
    tpl->Scale(xsecScale);
  }
}
template void TemplateHelpers::doTemplatePostprocessing<TH1F>(TH1F* tpl, const CategorizationHelpers::Category category, bool isMC);
template void TemplateHelpers::doTemplatePostprocessing<TH2F>(TH2F* tpl, const CategorizationHelpers::Category category, bool isMC);
template void TemplateHelpers::doTemplatePostprocessing<TH3F>(TH3F* tpl, const CategorizationHelpers::Category category, bool isMC);


#endif
