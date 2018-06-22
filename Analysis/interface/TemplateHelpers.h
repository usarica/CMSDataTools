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
  template<typename T> void setTemplateAxisLabels(T* histo, const CategorizationHelpers::Category category);
  template<typename T> void doTemplatePostprocessing(T* tpl, const CategorizationHelpers::Category category, bool isMC);

  void getLikelihoodDiscriminants(const SampleHelpers::Channel channel, const CategorizationHelpers::Category category, const SystematicsHelpers::SystematicVariationTypes syst, std::vector<DiscriminantClasses::KDspecs>& KDlist);
  void getCategorizationDiscriminants(const SystematicsHelpers::SystematicVariationTypes syst, std::vector<DiscriminantClasses::KDspecs>& KDlist);
  ExtendedBinning getDiscriminantFineBinning(const SampleHelpers::Channel channel, const CategorizationHelpers::Category category, ACHypothesisHelpers::ACHypothesis hypo, TString const strKD, CategorizationHelpers::MassRegion const massregion);
  ExtendedBinning getDiscriminantCoarseBinning(const SampleHelpers::Channel channel, const CategorizationHelpers::Category category, ACHypothesisHelpers::ACHypothesis hypo, TString const strKD, CategorizationHelpers::MassRegion const massregion);
  float getDiscriminantSmearingStrengthCoefficient(TString KDname, ProcessHandler::ProcessType proctype, CategorizationHelpers::MassRegion massregion);

  ProcessHandler const* getOnshellProcessHandler(ProcessHandler::ProcessType type);
  ProcessHandler const* getOffshellProcessHandler(ProcessHandler::ProcessType type);
  ProcessHandler const* getProcessHandlerPerMassRegion(ProcessHandler::ProcessType type, CategorizationHelpers::MassRegion massregion);

  /****************/
  /* Gluon fusion */
  /****************/
  const GGProcessHandler OnshellGGProcessHandle(CategorizationHelpers::kOnshell);
  const GGProcessHandler OffshellGGProcessHandle(CategorizationHelpers::kOffshell);


  /****************/
  /* EW VV fusion */
  /****************/
  const VVProcessHandler OnshellVVProcessHandle(CategorizationHelpers::kOnshell);
  const VVProcessHandler OffshellVVProcessHandle(CategorizationHelpers::kOffshell);


  /*******/
  /* VBF */
  /*******/
  const VVProcessHandler OnshellVBFProcessHandle(CategorizationHelpers::kOnshell, ProcessHandler::kVBF);
  const VVProcessHandler OffshellVBFProcessHandle(CategorizationHelpers::kOffshell, ProcessHandler::kVBF);


  /******/
  /* ZH */
  /******/
  const VVProcessHandler OnshellZHProcessHandle(CategorizationHelpers::kOnshell, ProcessHandler::kZH);
  const VVProcessHandler OffshellZHProcessHandle(CategorizationHelpers::kOffshell, ProcessHandler::kZH);


  /******/
  /* WH */
  /******/
  const VVProcessHandler OnshellWHProcessHandle(CategorizationHelpers::kOnshell, ProcessHandler::kWH);
  const VVProcessHandler OffshellWHProcessHandle(CategorizationHelpers::kOffshell, ProcessHandler::kWH);


  /*******/
  /* ttH */
  /*******/
  const TTProcessHandler OnshellTTProcessHandle(CategorizationHelpers::kOnshell);
  const TTProcessHandler OffshellTTProcessHandle(CategorizationHelpers::kOffshell);


  /*******/
  /* bbH */
  /*******/
  const BBProcessHandler OnshellBBProcessHandle(CategorizationHelpers::kOnshell);
  const BBProcessHandler OffshellBBProcessHandle(CategorizationHelpers::kOffshell);


  /*****************/
  /* QQ background */
  /*****************/
  const QQBkgProcessHandler OnshellQQBkgProcessHandle(CategorizationHelpers::kOnshell);
  const QQBkgProcessHandler OffshellQQBkgProcessHandle(CategorizationHelpers::kOffshell);


  /******************/
  /* Z+X background */
  /******************/
  const ZXProcessHandler OnshellZXProcessHandle(CategorizationHelpers::kOnshell);
  const ZXProcessHandler OffshellZXProcessHandle(CategorizationHelpers::kOffshell);
  ZXFakeRateHandler* getFakeRateHandler(ZXFakeRateHandler::FakeRateMethod FRMethod, SystematicsHelpers::SystematicVariationTypes syst);

}

#endif
