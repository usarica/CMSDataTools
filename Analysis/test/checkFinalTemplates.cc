#ifndef CHECKFINALTEMPLATES_VV_H
#define CHECKFINALTEMPLATES_VV_H

#include "common_includes.h"
#include "CheckSetTemplatesCategoryScheme.h"


// Constants to affect the template code
#ifndef outputdir_def
#define outputdir_def
const TString user_output_dir = "output/";
#endif

bool getFile(
  const Channel channel, const Category category, const ACHypothesis hypo, const SystematicVariationTypes syst,
  const unsigned int istage, const TString fixedDate,
  ProcessHandler const* thePerProcessHandle,
  std::vector<TFile*>& finputList
){
  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);
  const TString strACHypo = getACHypothesisName(hypo);
  const TString strStage = Form("Stage%i", istage);
  const TString strSystematics = getSystematicsName(syst);
  const TString strSystematicsOutput = getSystematicsCombineName(category, channel, thePerProcessHandle->getProcessType(), syst);

  // Setup the input directories
  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  TString cinput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/FinalTemplates/" + strStage + "/" + strACHypo;

  TString cinput = Form(
    "%s/HtoZZ%s_%s_FinalTemplates_%s_%s%s",
    cinput_common.Data(),
    strChannel.Data(), strCategory.Data(),
    thePerProcessHandle->getProcessName().Data(),
    strSystematicsOutput.Data(),
    ".root"
  );

  if (gSystem->AccessPathName(cinput)){
    MELAout << "getFilesAndTrees::File " << cinput << " is not found! Run " << strStage << " functions first." << endl;
    return false;
  }
  if (cinput!=""){
    TString timestamp = HostHelpers::GetTimestampConverted(cinput.Data());
    TFile* finput = TFile::Open(cinput, "read");
    if (finput){
      if (!finput->IsZombie()){
        MELAout << "getFile: Opening file " << cinput << " with timestamp " << timestamp << endl;
        finputList.push_back(finput);
      }
      else if (finput->IsOpen()){
        MELAout << "getFilesAndTrees::File " << cinput << " with timestamp " << timestamp << " is zombie! Re-run " << strStage << " functions first." << endl;
        finput->Close();
        return false;
      }
      else{
        MELAout << "getFilesAndTrees::File " << cinput << " with timestamp " << timestamp << " could not be opened! Re-run " << strStage << " functions first." << endl;
        return false;
      }
    }
  }
  return true;
}

void checkFinalTemplates_one(const Channel channel, const Category category, const ACHypothesis hypo, CategorizationHelpers::MassRegion massregion, const unsigned int istage=1, const TString fixedDate=""){
  if (channel==NChannels) return;
  if (!CheckSetTemplatesCategoryScheme(category)) return;

  vector<ProcessHandler::ProcessType> proctypes = { ProcessHandler::kGG, ProcessHandler::kVBF, ProcessHandler::kZH, ProcessHandler::kWH, ProcessHandler::kVV, ProcessHandler::kQQBkg, ProcessHandler::kZX };
  const unsigned int nproctypes=proctypes.size();

  TDirectory* rootdir=gDirectory;
  for (auto& proctype:proctypes){
    if (massregion==kOnshell && proctype==ProcessHandler::kVV) continue;

    vector<SystematicsHelpers::SystematicVariationTypes> allowedSysts = getProcessSystematicVariations(category, channel, proctype, "");
    for (auto& syst:allowedSysts){
      rootdir->cd();

      const TString strSystematics = getSystematicsName(syst);
      const TString strSystematicsOutput = getSystematicsCombineName(category, channel, proctype, syst);

      vector<TFile*> finputList;
      ProcessHandler const* inputProcessHandle=getProcessHandlerPerMassRegion(proctype, massregion);
      vector<TString> tplnamelist = inputProcessHandle->getTemplateNames(hypo, true);
      if (
        !getFile(
          channel, category, hypo, syst,
          istage, fixedDate,
          inputProcessHandle,
          finputList
        )
        ) continue;

      for (TFile*& finput:finputList){
        finput->cd();
        for (auto& tplname:tplnamelist){
          TH2F* htmp_2D;
          TH3F* htmp_3D;
          finput->GetObject(tplname, htmp_2D);
          finput->GetObject(tplname, htmp_3D);
          double integral=0, integralerror=0;
          if (!htmp_2D && !htmp_3D) MELAout << "Template " << tplname << " could not be found! BAD FILE...";
          else if (htmp_2D && htmp_3D) MELAout << "Template " << tplname << " is both 2 and 3D!? BAD FILE...";
          else if (htmp_2D){
            MELAout << "2D template " << tplname << " is GOOD. ";
            integral = getHistogramIntegralAndError(htmp_2D, 1, htmp_2D->GetNbinsX(), 1, htmp_2D->GetNbinsY(), true, &integralerror);
            MELAout << "Integral: " << integral << " +- " << integralerror << ".";
          }
          else /*if (htmp_3D)*/{
            MELAout << "3D template " << tplname << " is GOOD. ";
            integral = getHistogramIntegralAndError(htmp_3D, 1, htmp_3D->GetNbinsX(), 1, htmp_3D->GetNbinsY(), 1, htmp_3D->GetNbinsZ(), true, &integralerror);
            MELAout << "Integral: " << integral << " +- " << integralerror << ".";
          }
          MELAout << endl;
        }
        finput->Close();
      }
    }
  }
  rootdir->cd();
}


void checkFinalTemplates(CategorizationHelpers::MassRegion massregion, const unsigned int istage=1, const TString fixedDate=""){
  vector<Category> allowedCats = getAllowedCategories(globalCategorizationScheme);
  for (int ih=0; ih<nACHypotheses; ih++){
    ACHypothesis hypo = (ACHypothesis) ih;
    for (auto& category:allowedCats){
      for (int ch=0; ch<(int) NChannels; ch++){
        Channel channel = (Channel) ch;
        if (channel==k4l || channel==k2l2l) continue;

        checkFinalTemplates_one(channel, category, hypo, massregion, istage, fixedDate);
      }
    }
  }
}



#endif
