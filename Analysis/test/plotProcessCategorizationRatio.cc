#ifndef PLOTPROCESSCHECKSTAGE_H
#define PLOTPROCESSCHECKSTAGE_H

#include "common_includes.h"
#include "CheckSetTemplatesCategoryScheme.h"
#include "smoothenHistograms.h"
#include "TText.h"
#include "TPaveText.h"
#include "TColor.h"
#include "TArrayI.h"


// Constants to affect the template code
#ifndef outputdir_def
#define outputdir_def
const TString user_output_dir = "output/";
#endif


std::vector<SimpleEntry> orderTreeByMass(TTree* tree);
 
void getProcessCategorizationRatio(
  const Channel channel, const Category category, const ACHypothesis hypo,
  const TString fixedDate,
  ProcessHandler::ProcessType proctype,
  const TString strGenerator
){
  const unsigned int istage=2;
  const SystematicVariationTypes syst=sNominal;

  if (channel==NChannels) return;
  if (category==Inclusive) return;
  if (!CheckSetTemplatesCategoryScheme(category)) return;
  ProcessHandler const* thePerProcessHandle=nullptr;
  switch (proctype){
  case ProcessHandler::kGG:
    thePerProcessHandle = &TemplateHelpers::OffshellGGProcessHandle;
    break;
  case ProcessHandler::kVV:
    thePerProcessHandle = &TemplateHelpers::OffshellVVProcessHandle;
    break;
  case ProcessHandler::kQQBkg:
    thePerProcessHandle = &TemplateHelpers::OffshellQQBkgProcessHandle;
    break;
  default:
    break;
  };
  if (!thePerProcessHandle) return;
  if (!systematicAllowed(category, channel, thePerProcessHandle->getProcessType(), syst)) return;

  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);
  const TString strCategory_Inclusive = getCategoryName(Inclusive);
  const TString strACHypo = getACHypothesisName(hypo);
  TString strSystematics = getSystematicsName(syst);

  // Setup the output directories
  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;
  TString cinput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/";
  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/MassRatios/";

  gSystem->Exec("mkdir -p " + coutput_common);
  TString INPUT_INCLUSIVE_NAME = Form(
    "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_%s_Stage%i%s",
    strChannel.Data(), strCategory_Nominal.Data(),
    getACHypothesisName(fhypo).Data(),
    thePerProcessHandle->getProcessName().Data(),
    strSystematics.Data(),
    strGenerator.Data(),
    istage,
    ".root"
  );
  TString INPUT_NAME = Form(
    "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_%s_Stage%i%s",
    strChannel.Data(), strCategory.Data(),
    getACHypothesisName(fhypo).Data(),
    thePerProcessHandle->getProcessName().Data(),
    strSystematics.Data(),
    strGenerator.Data(),
    istage,
    ".root"
  );
  TString OUTPUT_NAME=INPUT_NAME; replaceString(OUTPUT_NAME, "FinalTemplates", "MassRatiosToInclusive");
  TString OUTPUT_LOG_NAME = OUTPUT_NAME;
  OUTPUT_NAME += "_plot.root";
  OUTPUT_LOG_NAME += "_plot.log";

  // Open the input files
  TString cinput[3] ={
    cinput_common + INPUT_INCLUSIVE_NAME,
    cinput_common + INPUT_NAME
  };
  TFile* finput[3] ={
    TFile::Open(cinput[0], "read"),
    TFile::Open(cinput[1], "read")
  };
  if (!finput[0] || !finput[1]) return;



  TString coutput = coutput_common + OUTPUT_NAME;
  TString coutput_log = coutput_common + OUTPUT_LOG_NAME;
  MELAout.open(coutput_log.Data());
  MELAout << "Opened log file " << coutput_log << endl;

  TFile* foutput = TFile::Open(coutput, "recreate");
  MELAout << "Opened output file " << coutput << endl;
  MELAout << "===============================" << endl;
  MELAout << "CoM Energy: " << theSqrts << " TeV" << endl;
  MELAout << "Decay Channel: " << strChannel << endl;
  MELAout << "===============================" << endl;
  MELAout << endl;

  gStyle->SetOptStat(0);
  {
    int colors[100];
    Double_t Red[]    ={ 0.3, 0.4, 1.0 };
    Double_t Green[]  ={ 0.0, 1.0, 0.8 };
    Double_t Blue[]   ={ 1.0, 0.0, 0.3 };
    Double_t Length[] ={ 0.00, 0.50, 1.00 };
    int FI = TColor::CreateGradientColorTable(3, Length, Red, Green, Blue, 100);
    const unsigned int ncolors = gStyle->GetNumberOfColors();
    if (FI<0) MELAout << "Failed to set color palette" << endl;
    else{
      for (unsigned int ic=0; ic<100; ic++) colors[ic] = FI+ic;
      gStyle->SetPalette(100, colors);
    }
    MELAout << "Ncolors: " << ncolors << endl;
  }

  std::vector<TH1F*>* historatios[2]={
    &(h1dlist[1]),
    &(h1dlist[2])
  };
  plotSystRatioTH1Fs(foutput, coutput_common, Form("c_M4lDistribution_%s_RatioToNominal", thePerProcessHandle->getProcessName().Data()), strSystematics, historatios);

  for (unsigned int is=0; is<3; is++) finput[is]->Close();
  foutput->Close();
  MELAout.close();
}

std::vector<SimpleEntry> orderTreeByMass(TTree* tree){

}


#endif
