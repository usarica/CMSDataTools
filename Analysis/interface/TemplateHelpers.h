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


  /*******/
  /* VBF */
  /*******/
  const VVProcessHandler OnshellVBFProcessHandle(false, ProcessHandler::kVBF);
  const VVProcessHandler OffshellVBFProcessHandle(true, ProcessHandler::kVBF);


  /******/
  /* ZH */
  /******/
  const VVProcessHandler OnshellZHProcessHandle(false, ProcessHandler::kZH);
  const VVProcessHandler OffshellZHProcessHandle(true, ProcessHandler::kZH);


  /******/
  /* WH */
  /******/
  const VVProcessHandler OnshellWHProcessHandle(false, ProcessHandler::kWH);
  const VVProcessHandler OffshellWHProcessHandle(true, ProcessHandler::kWH);


  /*****************/
  /* QQ background */
  /*****************/
  const QQBkgProcessHandler OnshellQQBkgProcessHandle(false);
  const QQBkgProcessHandler OffshellQQBkgProcessHandle(true);


  /******************/
  /* Z+X background */
  /******************/
  const ZXProcessHandler OnshellZXProcessHandle(false);
  const ZXProcessHandler OffshellZXProcessHandle(true);
  ZXFakeRateHandler* getFakeRateHandler(ZXFakeRateHandler::FakeRateMethod FRMethod, SystematicsHelpers::SystematicVariationTypes syst);

}

#endif
