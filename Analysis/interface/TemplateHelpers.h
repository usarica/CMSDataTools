#ifndef TEMPLATEHELPERS_H
#define TEMPLATEHELPERS_H

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "DiscriminantClasses.h"
#include "SampleHelpers.h"
#include "CategorizationHelpers.h"
#include "ACHypothesisHelpers.h"
#include "ProcessHandler.h"


namespace TemplateHelpers{

  /*********************/
  /* General functions */
  /*********************/
  template<typename T> void setTemplateAxisLabels(T* histo);
  void getLikelihoodDiscriminants(const SampleHelpers::Channel channel, const CategorizationHelpers::Category category, const TString strSystematics, std::vector<DiscriminantClasses::KDspecs>& KDlist);
  void getCategorizationDiscriminants(const TString strSystematics, std::vector<DiscriminantClasses::KDspecs>& KDlist);
  void adjustDiscriminantJECVariables(const TString strSystematics, std::vector<DiscriminantClasses::KDspecs>& KDlist);


  /****************/
  /* Gluon fusion */
  /****************/
  const GGProcessHandler OnshellGGProcessHandle(true);
  const GGProcessHandler OffshellGGProcessHandle(false);


  /****************/
  /* EW VV fusion */
  /****************/
  const VVProcessHandler OnshellVVProcessHandle(true);
  const VVProcessHandler OffshellVVProcessHandle(false);


  /*****************/
  /* QQ background */
  /*****************/
  const QQBkgProcessHandler OnshellQQBkgProcessHandle(true);
  const QQBkgProcessHandler OffshellQQBkgProcessHandle(false);

}

#endif
