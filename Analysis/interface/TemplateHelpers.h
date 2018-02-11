#ifndef TEMPLATEHELPERS_H
#define TEMPLATEHELPERS_H

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "ExtendedBinning.h"
#include "DiscriminantClasses.h"
#include "SampleHelpers.h"
#include "CategorizationHelpers.h"
#include "SystematicsHelpers.h"
#include "ACHypothesisHelpers.h"
#include "ProcessHandler.h"


namespace TemplateHelpers{

  /*********************/
  /* General functions */
  /*********************/
  template<typename T> void setTemplateAxisLabels(T* histo);
  void getLikelihoodDiscriminants(const SampleHelpers::Channel channel, const CategorizationHelpers::Category category, const SystematicsHelpers::SystematicVariationTypes syst, std::vector<DiscriminantClasses::KDspecs>& KDlist);
  void getCategorizationDiscriminants(const SystematicsHelpers::SystematicVariationTypes syst, std::vector<DiscriminantClasses::KDspecs>& KDlist);
  ExtendedBinning getDiscriminantBinning(const SampleHelpers::Channel channel, const CategorizationHelpers::Category category, TString const strKD, bool const useOffshell);

  /****************/
  /* Gluon fusion */
  /****************/
  const GGProcessHandler OnshellGGProcessHandle(false);
  const GGProcessHandler OffshellGGProcessHandle(true);


  /****************/
  /* EW VV fusion */
  /****************/
  const VVProcessHandler OnshellVVProcessHandle(false);
  const VVProcessHandler OffshellVVProcessHandle(true);


  /*****************/
  /* QQ background */
  /*****************/
  const QQBkgProcessHandler OnshellQQBkgProcessHandle(false);
  const QQBkgProcessHandler OffshellQQBkgProcessHandle(true);

}

#endif
