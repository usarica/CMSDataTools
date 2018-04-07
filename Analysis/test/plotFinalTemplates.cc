#ifndef PLOTFINALTEMPLATES_VV_H
#define PLOTFINALTEMPLATES_VV_H

#include "common_includes.h"
#include "CheckSetTemplatesCategoryScheme.h"
#include "TText.h"
#include "TPaveText.h"
#include "TColor.h"
#include "TArrayI.h"


// Constants to affect the template code
#ifndef outputdir_def
#define outputdir_def
const TString user_output_dir = "output/";
#endif


struct TemplateCollection{
  ACHypothesis hypo;
  Category category;
  Channel channel;
  ProcessHandler::ProcessType proctype;
  SystematicVariationTypes syst;

  vector<TH2F const*> templateList_2D;
  vector<TH3F const*> templateList_3D;
  vector<TFile*> finputList;

  TemplateCollection(){}
  TemplateCollection(ACHypothesis hypo_, Category category_, Channel channel_, ProcessHandler::ProcessType proctype_, SystematicVariationTypes syst_) :
    hypo(hypo_), category(category_), channel(channel_), proctype(proctype_), syst(syst_)
  {}
  TemplateCollection(TemplateCollection const& other) :
    hypo(other.hypo), category(other.category), channel(other.channel), proctype(other.proctype), syst(other.syst),
    templateList_2D(other.templateList_2D), templateList_3D(other.templateList_3D),
    finputList(other.finputList)
  {}

  void addTemplate(TH2F const* tpl){ templateList_2D.push_back(tpl); }
  void addTemplate(vector<TH2F const*> tpls){ for (TH2F const*& tpl:tpls) templateList_2D.push_back(tpl); }

  void addTemplate(TH3F const* tpl){ templateList_3D.push_back(tpl); }
  void addTemplate(vector<TH3F const*> tpls){ for (TH3F const*& tpl:tpls) templateList_3D.push_back(tpl); }

  void getTemplate(unsigned int i, TH2F const*& tpl){ tpl=nullptr; if (!templateList_2D.empty()) tpl=templateList_2D.at(i); }
  void getTemplate(unsigned int i, TH3F const*& tpl){ tpl=nullptr; if (!templateList_3D.empty()) tpl=templateList_3D.at(i); }

  void addInputFile(TFile* f){ finputList.push_back(f); }
  void addInputFile(vector<TFile*> fs){ for (TFile*& f:fs) finputList.push_back(f); }

  unsigned int getNTemplates(){ return templateList_2D.size()+templateList_3D.size(); }
  unsigned int getDimension(){ return 2*(!templateList_2D.empty())+3*(!templateList_3D.empty()); }

  void close(){ for (TFile*& f:finputList){ if (f && f->IsOpen()) f->Close(); } }
};
typedef unordered_map<SystematicsHelpers::SystematicVariationTypes, TemplateCollection, hash<int>> SystTplColl_t;
typedef unordered_map<Channel, TemplateCollection, hash<int>> ChanTplColl_t;
typedef unordered_map<ProcessHandler::ProcessType, SystTplColl_t, hash<int>> ProcSystTplCol_t;
typedef unordered_map<ProcessHandler::ProcessType, ChanTplColl_t, hash<int>> ProcChanTplCol_t;

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
    MELAerr << "getFilesAndTrees::File " << cinput << " is not found! Run " << strStage << " functions first." << endl;
    return false;
  }
  if (cinput!=""){
    TFile* finput = TFile::Open(cinput, "read");
    if (finput){
      if (!finput->IsZombie()){
        MELAout << "getFile: Opening file " << cinput << endl;
        finputList.push_back(finput);
      }
      else if (finput->IsOpen()){
        MELAerr << "getFilesAndTrees::File " << cinput << " is zombie! Re-run " << strStage << " functions first." << endl;
        finput->Close();
        return false;
      }
      else{
        MELAerr << "getFilesAndTrees::File " << cinput << " could not be opened! Re-run " << strStage << " functions first." << endl;
        return false;
      }
    }
  }
  return true;
}

void plotSystToNominal1D(TString coutdir, SystTplColl_t& tplCollList, ProcessHandler const* thePerProcessHandle, vector<SystematicsHelpers::SystematicVariationTypes> const& allowedSysts, bool onlyNominal);

void plotSumChannels(const Category category, const ACHypothesis hypo, CategorizationHelpers::MassRegion massregion, ProcChanTplCol_t& tplCollList, std::vector<ProcessHandler const*> const& procHandlers);

void plotFinalTemplates_one(const Channel channel, const Category category, const ACHypothesis hypo, CategorizationHelpers::MassRegion massregion, const unsigned int istage=1, const TString fixedDate=""){
  if (channel==NChannels) return;
  if (!CheckSetTemplatesCategoryScheme(category)) return;

  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);
  const TString strACHypo = getACHypothesisName(hypo);
  const TString strStage = Form("Stage%i", istage);

  // Setup the input directories
  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/FinalTemplates/" + strStage + "/Validation/" + strACHypo + "/";
  switch (massregion){
  case kOffshell:
    coutput_common += "Offshell"; break;
  case kOnshell:
    coutput_common += "Onshell"; break;
  default:
    assert(0);
  }
  gSystem->Exec("mkdir -p " + coutput_common);

  TString OUTPUT_NAME = Form(
    "%s/HtoZZ%s_%s_Plots",
    coutput_common.Data(),
    strChannel.Data(), strCategory.Data()
  );
  TString OUTPUT_LOG_NAME = OUTPUT_NAME;
  OUTPUT_NAME += ".root";
  OUTPUT_LOG_NAME += ".log";

  MELAout.open(OUTPUT_LOG_NAME);
  TFile* foutput=TFile::Open(OUTPUT_NAME, "recreate");

  vector<ProcessHandler::ProcessType> proctypes = { ProcessHandler::kGG, ProcessHandler::kVBF, ProcessHandler::kZH, ProcessHandler::kWH, ProcessHandler::kVV, ProcessHandler::kQQBkg, ProcessHandler::kZX };
  const unsigned int nproctypes=proctypes.size();

  ProcSystTplCol_t proctpls;
  for (auto& proctype:proctypes) proctpls[proctype]=SystTplColl_t();

  TDirectory* rootdir=gDirectory;
  for (auto& proctype:proctypes){
    ProcessHandler const* inputProcessHandle=getProcessHandlerPerMassRegion(proctype, massregion);
    vector<TString> tplnamelist = inputProcessHandle->getTemplateNames(hypo, true);

    vector<SystematicsHelpers::SystematicVariationTypes> allowedSysts = getProcessSystematicVariations(category, channel, proctype, "");
    for (auto& syst:allowedSysts){
      rootdir->cd();

      const TString strSystematics = getSystematicsName(syst);
      const TString strSystematicsOutput = getSystematicsCombineName(category, channel, proctype, syst);

      vector<TFile*> finputList;
      if (
        !getFile(
          channel, category, hypo, syst,
          istage, fixedDate,
          inputProcessHandle,
          finputList
        )
        ) continue;
      rootdir->cd();

      proctpls[proctype][syst]=TemplateCollection(hypo, category, channel, proctype, syst);
      assert(finputList.size()<=1);
      for (TFile*& finput:finputList){
        proctpls[proctype][syst].addInputFile(finput);

        finput->cd();
        for (auto& tplname:tplnamelist){
          TH2F* htmp_2D;
          TH3F* htmp_3D;
          finput->GetObject(tplname, htmp_2D);
          finput->GetObject(tplname, htmp_3D);
          if (!htmp_2D && !htmp_3D) MELAerr << "Template " << tplname << " could not be found! BAD FILE..." << endl;
          else if (htmp_2D && htmp_3D) MELAerr << "Template " << tplname << " is both 2 and 3D!? BAD FILE..." << endl;
          else if (htmp_2D) proctpls[proctype][syst].addTemplate(htmp_2D);
          else proctpls[proctype][syst].addTemplate(htmp_3D);
        }
        rootdir->cd();
      }

    }

    foutput->cd();
    plotSystToNominal1D(coutput_common, proctpls[proctype], inputProcessHandle, allowedSysts, false);
    plotSystToNominal1D(coutput_common, proctpls[proctype], inputProcessHandle, allowedSysts, true);

  }
  rootdir->cd();

  // Close all files
  for (auto it_p:proctpls){
    for (auto it_s:it_p.second){
      it_s.second.close();
    }
  }

  foutput->Close();
  MELAout.close();
}

void plotFinalTemplatesStacked_one(const Category category, const ACHypothesis hypo, CategorizationHelpers::MassRegion massregion, const unsigned int istage=1, const TString fixedDate=""){
  if (!CheckSetTemplatesCategoryScheme(category)) return;

  const TString strCategory = getCategoryName(category);
  const TString strACHypo = getACHypothesisName(hypo);
  const TString strStage = Form("Stage%i", istage);

  // Setup the input directories
  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/FinalTemplates/" + strStage + "/Validation/" + strACHypo + "/";
  switch (massregion){
  case kOffshell:
    coutput_common += "Offshell"; break;
  case kOnshell:
    coutput_common += "Onshell"; break;
  default:
    assert(0);
  }
  gSystem->Exec("mkdir -p " + coutput_common);

  TString OUTPUT_NAME = Form(
    "%s/HtoZZ_AllChannels_%s_Plots",
    coutput_common.Data(), strCategory.Data()
  );
  TString OUTPUT_LOG_NAME = OUTPUT_NAME;
  OUTPUT_NAME += ".root";
  OUTPUT_LOG_NAME += ".log";

  MELAout.open(OUTPUT_LOG_NAME);
  TFile* foutput=TFile::Open(OUTPUT_NAME, "recreate");

  vector<ProcessHandler::ProcessType> proctypes;
  proctypes.push_back(ProcessHandler::kGG);
  proctypes.push_back(ProcessHandler::kVV);
  proctypes.push_back(ProcessHandler::kQQBkg);
  proctypes.push_back(ProcessHandler::kZX);

  const unsigned int nproctypes=proctypes.size();

  ProcChanTplCol_t proctpls;
  for (auto& proctype:proctypes) proctpls[proctype]=ChanTplColl_t();

  TDirectory* rootdir=gDirectory;
  std::vector<ProcessHandler const*> procHandlers; procHandlers.reserve(nproctypes);
  for (auto& proctype:proctypes){
    ProcessHandler const* inputProcessHandle=getProcessHandlerPerMassRegion(proctype, massregion);
    procHandlers.push_back(inputProcessHandle);

    vector<TString> tplnamelist = inputProcessHandle->getTemplateNames(hypo, true);

    for (int ch=0; ch<(int) NChannels; ch++){
      Channel channel = (Channel) ch;
      if (channel==k4l || channel==k2l2l) continue;

      rootdir->cd();

      vector<TFile*> finputList;
      if (
        !getFile(
          channel, category, hypo, sNominal,
          istage, fixedDate,
          inputProcessHandle,
          finputList
        )
        ) continue;
      rootdir->cd();

      proctpls[proctype][channel]=TemplateCollection(hypo, category, channel, proctype, sNominal);
      assert(finputList.size()<=1);
      for (TFile*& finput:finputList){
        proctpls[proctype][channel].addInputFile(finput);

        finput->cd();
        for (auto& tplname:tplnamelist){
          TH2F* htmp_2D;
          TH3F* htmp_3D;
          finput->GetObject(tplname, htmp_2D);
          finput->GetObject(tplname, htmp_3D);
          if (!htmp_2D && !htmp_3D) MELAerr << "Template " << tplname << " could not be found! BAD FILE..." << endl;
          else if (htmp_2D && htmp_3D) MELAerr << "Template " << tplname << " is both 2 and 3D!? BAD FILE..." << endl;
          else if (htmp_2D) proctpls[proctype][channel].addTemplate(htmp_2D);
          else proctpls[proctype][channel].addTemplate(htmp_3D);
        }
        rootdir->cd();
      }

    }
  }

  foutput->cd();
  plotSumChannels(category, hypo, massregion, proctpls, procHandlers);

  rootdir->cd();

  // Close all files
  for (auto it_p:proctpls){
    for (auto it_s:it_p.second){
      it_s.second.close();
    }
  }

  vector<TCanvas*> canvasList;
  extractHistogramsFromDirectory(foutput, canvasList);
  for (TCanvas*& canvas:canvasList){
    canvas->SaveAs(coutput_common + "/" + canvas->GetName() + ".png");
    canvas->SaveAs(coutput_common + "/" + canvas->GetName() + ".pdf");
    canvas->Close();
  }

  foutput->Close();
  MELAout.close();
}


void plotFinalTemplates(CategorizationHelpers::MassRegion massregion, const unsigned int istage=1, const TString fixedDate=""){
  vector<Category> allowedCats = getAllowedCategories(globalCategorizationScheme);
  for (int ch=0; ch<(int) NChannels; ch++){
    Channel channel = (Channel) ch;
    if (channel==k4l || channel==k2l2l) continue;
    for (auto& category:allowedCats){
      for (int ih=0; ih<nACHypotheses; ih++){
        ACHypothesis hypo = (ACHypothesis) ih;

        plotFinalTemplates_one(channel, category, hypo, massregion, istage, fixedDate);
      }
    }
  }
}

void plotSystToNominal1D(TString coutdir, SystTplColl_t& tplCollList, ProcessHandler const* thePerProcessHandle, vector<SystematicsHelpers::SystematicVariationTypes> const& allowedSysts, bool onlyNominal){
  MELAout << "Begin plotSystToNominal1D" << endl;

  TDirectory* curdir=gDirectory;
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

  const TString strProcName = thePerProcessHandle->getProcessName();
  const TString ytitle = "N (a.u.)";

  TemplateCollection& tpls_Nominal = tplCollList[sNominal];
  const unsigned int nKDs=tpls_Nominal.getDimension();
  vector<TString> KDset; KDset.push_back("ZZMass"); { vector<TString> KDset2=getACHypothesisKDNameSet(tpls_Nominal.hypo, tpls_Nominal.category, thePerProcessHandle->getProcessMassRegion()); appendVector(KDset, KDset2); }
  if (nKDs==0) return;
  for (unsigned isyst=0; isyst<(allowedSysts.size()-1)/2; isyst++){
    if (onlyNominal && isyst!=0) continue;
    auto const& systDn = allowedSysts.at(isyst*2+1);
    auto const& systUp = allowedSysts.at(isyst*2+2);
    TemplateCollection& tpls_Dn = tplCollList[systDn];
    TemplateCollection& tpls_Up = tplCollList[systUp];

    if (!onlyNominal) MELAout << "\t- Checking systematic " << isyst << endl;
    else MELAout << "\t- Checking systematic nominal" << endl;

    for (unsigned int idim=0; idim<nKDs; idim++){
      MELAout << "\t- Checking dimension " << idim << endl;
      MELAout << "\t- Ntpls: " << tpls_Nominal.getNTemplates() << endl;
      vector<TH1F*> hProj[3];
      double scale=0;
      for (unsigned int itpl=0; itpl<tpls_Nominal.getNTemplates(); itpl++){
        MELAout << "\t- nKDs = " << nKDs << endl;
        if (nKDs==2){
          TH2F const* hOriginal[3];
          tpls_Nominal.getTemplate(itpl, hOriginal[0]);
          tpls_Dn.getTemplate(itpl, hOriginal[1]);
          tpls_Up.getTemplate(itpl, hOriginal[2]);

          for (unsigned int ii=0; ii<3; ii++){
            TString newname = hOriginal[ii]->GetName(); newname += Form("_%i", ii);
            TH2F* hOriginal_copy = new TH2F(*(hOriginal[ii])); multiplyBinWidth(hOriginal_copy);
            hProj[ii].push_back(getHistogramSlice(
              hOriginal_copy,
              idim, 1, (hOriginal[ii]->GetNbinsY())*(idim==0)+(hOriginal[ii]->GetNbinsX())*(idim==1),
              newname
            ));
            delete hOriginal_copy;

            divideBinWidth(hProj[ii].back());
            hProj[ii].back()->SetTitle(hOriginal[0]->GetTitle());
            TString xtitle;
            switch (idim){
            case 0:
              xtitle=hOriginal[0]->GetXaxis()->GetTitle();
              break;
            case 1:
              xtitle=hOriginal[0]->GetYaxis()->GetTitle();
              break;
            }
            hProj[ii].back()->GetXaxis()->SetTitle(xtitle);
            hProj[ii].back()->GetYaxis()->SetTitle(ytitle);
            if (ii==0 && itpl==0) scale = hProj[ii].back()->Integral("width");
            if (scale>0.) hProj[ii].back()->Scale(1./scale);
          }
        }
        else if (nKDs==3){
          TH3F const* hOriginal[3];
          tpls_Nominal.getTemplate(itpl, hOriginal[0]);
          tpls_Dn.getTemplate(itpl, hOriginal[1]);
          tpls_Up.getTemplate(itpl, hOriginal[2]);

          for (unsigned int ii=0; ii<3; ii++){
            TString newname = hOriginal[ii]->GetName(); newname += Form("_%i", ii);
            TH3F* hOriginal_copy = new TH3F(*(hOriginal[ii])); multiplyBinWidth(hOriginal_copy);
            hProj[ii].push_back(getHistogramSlice(
              hOriginal_copy,
              idim,
              1, (hOriginal[ii]->GetNbinsY())*(idim==0)+(hOriginal[ii]->GetNbinsZ())*(idim==1)+(hOriginal[ii]->GetNbinsX())*(idim==2),
              1, (hOriginal[ii]->GetNbinsZ())*(idim==0)+(hOriginal[ii]->GetNbinsX())*(idim==1)+(hOriginal[ii]->GetNbinsY())*(idim==2),
              newname
            ));
            delete hOriginal_copy;

            divideBinWidth(hProj[ii].back());
            hProj[ii].back()->SetTitle(hOriginal[0]->GetTitle());
            TString xtitle;
            switch (idim){
            case 0:
              xtitle=hOriginal[0]->GetXaxis()->GetTitle();
              break;
            case 1:
              xtitle=hOriginal[0]->GetYaxis()->GetTitle();
              break;
            case 2:
              xtitle=hOriginal[0]->GetZaxis()->GetTitle();
              break;
            }
            hProj[ii].back()->GetXaxis()->SetTitle(xtitle);
            hProj[ii].back()->GetYaxis()->SetTitle(ytitle);
            if (ii==0 && itpl==0) scale = hProj[ii].back()->Integral("width");
            if (scale>0.) hProj[ii].back()->Scale(1./scale);
          }
        }
      }

      if (tpls_Nominal.hypo!=kSM){
        if (thePerProcessHandle->getProcessType()==ProcessHandler::kGG){
          double integralSM=hProj[0].at((int) GGProcessHandler::GGTplSig)->Integral("width");
          double integralBSM=hProj[0].at((int) GGProcessHandler::GGTplSigBSM)->Integral("width");
          double scale=integralSM/integralBSM;
          for (int ii=0; ii<3; ii++){
            hProj[ii].at((int) GGProcessHandler::GGTplSigBSM)->Scale(scale);
            hProj[ii].at((int) GGProcessHandler::GGTplSigBSMSMInt_Re)->Scale(sqrt(scale));
            hProj[ii].at((int) GGProcessHandler::GGTplIntBSM_Re)->Scale(sqrt(scale));
          }
        }
        else if (thePerProcessHandle->getProcessType()==ProcessHandler::kVV){
          double integralSM=hProj[0].at((int) VVProcessHandler::VVTplSig)->Integral("width");
          double integralBSM=hProj[0].at((int) VVProcessHandler::VVTplSigBSM)->Integral("width");
          double scale=integralSM/integralBSM;
          for (int ii=0; ii<3; ii++){
            hProj[ii].at((int) VVProcessHandler::VVTplSigBSM)->Scale(scale);
            hProj[ii].at((int) VVProcessHandler::VVTplSigBSMSMInt_ai1_1_Re)->Scale(pow(scale, 0.25));
            hProj[ii].at((int) VVProcessHandler::VVTplSigBSMSMInt_ai1_2_PosDef)->Scale(pow(scale, 0.5));
            hProj[ii].at((int) VVProcessHandler::VVTplSigBSMSMInt_ai1_3_Re)->Scale(pow(scale, 0.75));
            hProj[ii].at((int) VVProcessHandler::VVTplIntBSM_ai1_1_Re)->Scale(pow(scale, 0.25));
            hProj[ii].at((int) VVProcessHandler::VVTplIntBSM_ai1_2_Re)->Scale(pow(scale, 0.5));
          }
        }
      }

      // Plot!
      // Adjust colors
      {
        const unsigned int ncolors = gStyle->GetNumberOfColors();
        const unsigned int stepsize = ncolors/hProj[0].size()/3;
        for (int ii=0; ii<3; ii++){
          for (unsigned int ih=0; ih<hProj[ii].size(); ih++){
            int colorToUse = gStyle->GetColorPalette((3*ih+ii)*stepsize);
            hProj[ii].at(ih)->SetMarkerColor(colorToUse);
            hProj[ii].at(ih)->SetLineColor(colorToUse);
            hProj[ii].at(ih)->SetLineWidth(2);
          }
        }
      }

      const TString strChannel = getChannelName(tpls_Nominal.channel);
      const TString strChannelLabel = getChannelLabel(tpls_Nominal.channel);
      const TString strCategory = getCategoryName(tpls_Nominal.category);
      const TString strCategoryLabel = getCategoryLabel(tpls_Nominal.category);
      const TString strACHypo = getACHypothesisName(tpls_Nominal.hypo);
      TString strSystematics = getSystematicsName(systUp); replaceString(strSystematics, "Up", "");
      if (onlyNominal) strSystematics = getSystematicsName(sNominal);
      TString strSystematicsLabel = getSystematicsLabel(systUp); replaceString(strSystematicsLabel, " up", "");
      TString canvasname = Form(
        "c_%s_%s_%s_%s_%s_%s",
        strACHypo.Data(), strChannel.Data(), strCategory.Data(),
        strProcName.Data(),
        strSystematics.Data(),
        KDset.at(idim).Data()
      );

      MELAout << "\t- Constructing canvas " << canvasname << endl;

      TCanvas* canvas = new TCanvas(canvasname, "", 8, 30, 800, 800);
      canvas->cd();
      gStyle->SetOptStat(0);
      canvas->SetFillColor(0);
      canvas->SetBorderMode(0);
      canvas->SetBorderSize(2);
      canvas->SetTickx(1);
      canvas->SetTicky(1);
      canvas->SetLeftMargin(0.17);
      canvas->SetRightMargin(0.05);
      canvas->SetTopMargin(0.07);
      canvas->SetBottomMargin(0.13);
      canvas->SetFrameFillStyle(0);
      canvas->SetFrameBorderMode(0);
      canvas->SetFrameFillStyle(0);
      canvas->SetFrameBorderMode(0);

      TLegend* legend = new TLegend(
        (thePerProcessHandle->getProcessType()==ProcessHandler::kQQBkg || thePerProcessHandle->getProcessType()==ProcessHandler::kZX ? 0.80 : 0.44),
        0.90-0.10/3.*2.*float(hProj[0].size()),
        (thePerProcessHandle->getProcessType()==ProcessHandler::kQQBkg || thePerProcessHandle->getProcessType()==ProcessHandler::kZX ? 0.90 : 0.75),
        0.90
      );
      legend->SetBorderSize(0);
      legend->SetTextFont(42);
      legend->SetTextSize(0.03);
      legend->SetLineColor(1);
      legend->SetLineStyle(1);
      legend->SetLineWidth(1);
      legend->SetFillColor(0);
      legend->SetFillStyle(0);

      TText* text;
      TPaveText* pt = new TPaveText(0.15, 0.93, 0.85, 1, "brNDC");
      pt->SetBorderSize(0);
      pt->SetFillStyle(0);
      pt->SetTextAlign(12);
      pt->SetTextFont(42);
      pt->SetTextSize(0.045);
      text = pt->AddText(0.025, 0.45, "#font[61]{CMS}");
      text->SetTextSize(0.044);
      text = pt->AddText(0.165, 0.42, "#font[52]{Simulation}");
      text->SetTextSize(0.0315);
      TString cErgTev = Form("#font[42]{1 fb^{-1} %i TeV}", theSqrts);
      text = pt->AddText(0.9, 0.45, cErgTev);
      text->SetTextSize(0.0315);

      TPaveText* ptCatChannel = new TPaveText(0.20, 0.85, 0.50, 0.90, "brNDC");
      ptCatChannel->SetBorderSize(0);
      ptCatChannel->SetFillStyle(0);
      ptCatChannel->SetTextAlign(12);
      ptCatChannel->SetTextFont(42);
      ptCatChannel->SetTextSize(0.045);
      text = ptCatChannel->AddText(0.025, 0.45, Form("%s, %s", strCategoryLabel.Data(), strChannelLabel.Data()));
      text->SetTextSize(0.0315);

      TPaveText* ptSyst=nullptr;
      if (!onlyNominal){
        ptSyst = new TPaveText(0.20, 0.80, 0.50, 0.85, "brNDC");
        ptSyst->SetBorderSize(0);
        ptSyst->SetFillStyle(0);
        ptSyst->SetTextAlign(12);
        ptSyst->SetTextFont(42);
        ptSyst->SetTextSize(0.045);
        text = ptSyst->AddText(0.025, 0.45, strSystematicsLabel+" syst.");
        text->SetTextSize(0.0315);
      }

      bool first=true;
      float minY=9e9, maxY=-9e9;
      for (unsigned int ih=0; ih<hProj[0].size(); ih++){
        for (unsigned int ii=0; ii<3; ii++){
          TH1F*& hh=hProj[ii].at(ih);
          float minX=hh->GetXaxis()->GetBinLowEdge(1);
          float maxX=hh->GetXaxis()->GetBinUpEdge(hh->GetNbinsX());
          if (thePerProcessHandle->getProcessMassRegion()==kOffshell && KDset.at(idim)=="ZZMass") maxX=1500;
          int binXlow = hh->GetXaxis()->FindBin(minX);
          int binXhigh = hh->GetXaxis()->FindBin(maxX);

          hh->GetXaxis()->SetRangeUser(minX, maxX);
          hh->GetXaxis()->SetNdivisions(505);
          hh->GetXaxis()->SetLabelFont(42);
          hh->GetXaxis()->SetLabelOffset(0.007);
          hh->GetXaxis()->SetLabelSize(0.04);
          hh->GetXaxis()->SetTitleSize(0.06);
          hh->GetXaxis()->SetTitleOffset(0.9);
          hh->GetXaxis()->SetTitleFont(42);
          hh->GetYaxis()->SetNdivisions(505);
          hh->GetYaxis()->SetLabelFont(42);
          hh->GetYaxis()->SetLabelOffset(0.007);
          hh->GetYaxis()->SetLabelSize(0.04);
          hh->GetYaxis()->SetTitleSize(0.06);
          hh->GetYaxis()->SetTitleOffset(1.1);
          hh->GetYaxis()->SetTitleFont(42);
          hh->GetYaxis()->SetTitle(ytitle);
          for (int ix=binXlow; ix<binXhigh; ix++){
            float bc = hh->GetBinContent(ix);
            if (bc!=0.){
              minY = std::min(bc, minY);
              maxY = std::max(bc, maxY);
            }
          }
          MELAout << "Min = " << minY << ", max = " << maxY << " after " << hh->GetName() << " (" << hh->GetTitle() << ")" << endl;
          if (ii==0) legend->AddEntry(hh, hh->GetTitle(), "l");
          else if (!onlyNominal && hProj[0].size()>1) legend->AddEntry(hh, " ", "l");
          hh->SetTitle("");
        }
      }
      for (unsigned int ii=0; ii<3; ii++){
        if (ii>0 && onlyNominal) continue;
        for (TH1F*& hh:hProj[ii]){
          hh->GetYaxis()->SetRangeUser(minY*(minY<0. ? 1.1 : 0.9), maxY*1.3);
          if (first){
            hh->Draw("hist");
            first=false;
          }
          else hh->Draw("histsame");
        }
      }
      legend->Draw("same");
      pt->Draw();
      ptCatChannel->Draw();
      if (ptSyst) ptSyst->Draw();
      canvas->RedrawAxis();
      canvas->Modified();
      canvas->Update();
      curdir->WriteTObject(canvas);
      canvas->SaveAs(coutdir + "/" + canvas->GetName() + ".png");
      canvas->SaveAs(coutdir + "/" + canvas->GetName() + ".pdf");
      delete ptSyst;
      delete ptCatChannel;
      delete pt;
      delete legend;
      canvas->Close();

      for (unsigned int ii=0; ii<3; ii++){ for (TH1F*& hh:hProj[ii]) delete hh; }
      MELAout << "\t- Done with dimension " << idim << endl;
    }
  }
}

void plotSumChannels(const Category category, const ACHypothesis hypo, CategorizationHelpers::MassRegion massregion, ProcChanTplCol_t& tplCollList, std::vector<ProcessHandler const*> const& procHandlers){
  MELAout << "Begin plotSumChannels" << endl;

  const TString strCategory = getCategoryName(category);
  const TString strCategoryLabel = getCategoryLabel(category);
  const TString strACHypo = getACHypothesisName(hypo);

  float lumi=1;
  if (theDataPeriod=="2016") lumi=35.8671;
  else assert(0);
  float achypoF=0.15;
  float signal_strength_GG=1;
  float signal_strength_VV=1;
  TString acname="SM";
  TString customHypoTitle_GG, customHypoTitle_VV;
  switch (hypo){
  case ACHypothesisHelpers::kL1:
    acname=Form("f_{#Lambda1}=%.2f", achypoF);
    signal_strength_GG=1.2;
    signal_strength_VV=0.01;
    break;
  case ACHypothesisHelpers::kA2:
    acname=Form("f_{a2}=%.2f", achypoF);
    signal_strength_GG=1.2;
    signal_strength_VV=0.01;
    break;
  case ACHypothesisHelpers::kA3:
    acname=Form("f_{a3}=%.2f", achypoF);
    signal_strength_GG=1.2;
    signal_strength_VV=0.01;
    break;
  default:
    achypoF=0;
    signal_strength_GG=10;
    signal_strength_VV=signal_strength_GG;
    break;
  };
  customHypoTitle_GG = Form((hypo==kSM ? "%s, %s%.0f" : "%s, %s%.2f"), acname.Data(), (hypo==kSM ? "#Gamma_{H}=" : "#mu_{F}="), signal_strength_GG);
  customHypoTitle_VV = Form((hypo==kSM ? "%s, %s%.0f" : "%s, %s%.2f"), acname.Data(), (hypo==kSM ? "#Gamma_{H}=" : "#mu_{V}="), signal_strength_VV);

  vector<TString> KDset; KDset.push_back("ZZMass"); { vector<TString> KDset2=getACHypothesisKDNameSet(hypo, category, massregion); appendVector(KDset, KDset2); }
  const unsigned int nKDs=KDset.size();
  if (nKDs==0) return;

  TDirectory* curdir=gDirectory;
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

  const unsigned int nproctypes=procHandlers.size();

  for (unsigned int idim=0; idim<nKDs; idim++){
    const TString ytitle = "#events / lumi.";
    unordered_map<ProcessHandler::ProcessType, TH1F*, std::hash<int>> proc_hist_SM;
    unordered_map<ProcessHandler::ProcessType, TH1F*, std::hash<int>> proc_hist_BSM;
    for (unsigned int ip=0; ip<procHandlers.size(); ip++){
      auto& thePerProcessHandle = procHandlers.at(ip);
      auto proctype=thePerProcessHandle->getProcessType();
      ChanTplColl_t& chantpls = tplCollList[proctype];
      for (int ch=0; ch<(int) NChannels; ch++){
        Channel channel = (Channel) ch;
        if (channel==k4l || channel==k2l2l) continue;

        TemplateCollection& tpls = chantpls[channel];
        vector<TH1F*> tmptpls; // Template projections per channel in each process
        for (unsigned int itpl=0; itpl<tpls.getNTemplates(); itpl++){
          TH1F* htmp;
          if (nKDs==2){
            TH2F const* hOriginal;
            tpls.getTemplate(itpl, hOriginal);

            TString newname = hOriginal->GetName(); newname += Form("_%i", idim);
            TH2F* hOriginal_copy = new TH2F(*hOriginal); multiplyBinWidth(hOriginal_copy);
            htmp=getHistogramSlice(
              hOriginal_copy,
              idim, 1, (hOriginal->GetNbinsY())*(idim==0)+(hOriginal->GetNbinsX())*(idim==1),
              newname
            );
            delete hOriginal_copy;

            htmp->SetTitle(hOriginal->GetTitle());
            TString xtitle;
            switch (idim){
            case 0:
              xtitle=hOriginal->GetXaxis()->GetTitle();
              break;
            case 1:
              xtitle=hOriginal->GetYaxis()->GetTitle();
              break;
            }
            htmp->GetXaxis()->SetTitle(xtitle);
            htmp->GetYaxis()->SetTitle(ytitle);
          }
          else if (nKDs==3){
            TH3F const* hOriginal;
            tpls.getTemplate(itpl, hOriginal);

            TString newname = hOriginal->GetName(); newname += Form("_%i", idim);
            TH3F* hOriginal_copy = new TH3F(*hOriginal); multiplyBinWidth(hOriginal_copy);
            htmp=getHistogramSlice(
              hOriginal_copy,
              idim,
              1, (hOriginal->GetNbinsY())*(idim==0)+(hOriginal->GetNbinsZ())*(idim==1)+(hOriginal->GetNbinsX())*(idim==2),
              1, (hOriginal->GetNbinsZ())*(idim==0)+(hOriginal->GetNbinsX())*(idim==1)+(hOriginal->GetNbinsY())*(idim==2),
              newname
            );
            delete hOriginal_copy;

            htmp->SetTitle(hOriginal->GetTitle());
            TString xtitle;
            switch (idim){
            case 0:
              xtitle=hOriginal->GetXaxis()->GetTitle();
              break;
            case 1:
              xtitle=hOriginal->GetYaxis()->GetTitle();
              break;
            case 2:
              xtitle=hOriginal->GetZaxis()->GetTitle();
              break;
            }
            htmp->GetXaxis()->SetTitle(xtitle);
            htmp->GetYaxis()->SetTitle(ytitle);
          }
          else assert(0);
          divideBinWidth(htmp);

          // Yield hacks
          if (proctype==ProcessHandler::kQQBkg) htmp->Scale(1000); // FIXME
          if (proctype!=ProcessHandler::kZX) htmp->Scale(lumi);

          tmptpls.push_back(htmp);
        }
        TH1F* hRes_SM = new TH1F(*(tmptpls.front()));
        TH1F* hRes_BSM = new TH1F(*(tmptpls.front()));
        if (proctype==ProcessHandler::kGG){
          ((GGProcessHandler const*) thePerProcessHandle)->getHypothesisHistogramFromTemplates(hRes_SM, tmptpls, hypo, 1, 1, 0);
          hRes_SM->SetTitle("gg #rightarrow 4l total, SM");
          ((GGProcessHandler const*) thePerProcessHandle)->getHypothesisHistogramFromTemplates(hRes_BSM, tmptpls, hypo, 1, signal_strength_GG*(1.-achypoF), signal_strength_GG*achypoF);
          hRes_BSM->SetTitle(Form("gg #rightarrow 4l total, %s", customHypoTitle_GG.Data()));
        }
        else if (proctype==ProcessHandler::kVV){
          ((VVProcessHandler const*) thePerProcessHandle)->getHypothesisHistogramFromTemplates(hRes_SM, tmptpls, hypo, 1, 1, 0);
          hRes_SM->SetTitle("VV #rightarrow 4l total, SM");
          ((VVProcessHandler const*) thePerProcessHandle)->getHypothesisHistogramFromTemplates(hRes_BSM, tmptpls, hypo, 1, signal_strength_VV*(1.-achypoF), signal_strength_VV*achypoF);
          hRes_BSM->SetTitle(Form("VV #rightarrow 4l total, %s", customHypoTitle_VV.Data()));
        }
        for (auto& tmp:tmptpls) delete tmp;

        auto it = proc_hist_SM.find(proctype);
        if (it==proc_hist_SM.end()) proc_hist_SM[proctype]=hRes_SM;
        else{
          it->second->Add(hRes_SM);
          delete hRes_SM;
        }
        it = proc_hist_BSM.find(proctype);
        if (it==proc_hist_BSM.end()) proc_hist_BSM[proctype]=hRes_BSM;
        else{
          it->second->Add(hRes_BSM);
          delete hRes_BSM;
        }
      } // End loop over channels
    } // End loop over processes


    std::vector<TH1F*> histToPlot; histToPlot.reserve(nproctypes+2);
    vector<ProcessHandler::ProcessType> procorder{ ProcessHandler::kQQBkg, ProcessHandler::kZX, ProcessHandler::kGG, ProcessHandler::kVV };
    for (unsigned int ip=0; ip<procorder.size(); ip++){
      if (ip>0) proc_hist_SM[procorder.at(ip)]->Add(proc_hist_SM[procorder.at(ip-1)]);
      histToPlot.push_back(proc_hist_SM[procorder.at(ip)]);
    }
    for (unsigned int ip=1; ip<procorder.size(); ip++){
      proc_hist_BSM[procorder.at(ip)]->Add(proc_hist_BSM[procorder.at(ip-1)]);
      if (ip>1) histToPlot.push_back(proc_hist_BSM[procorder.at(ip)]);
    }
    const unsigned int nplotted=histToPlot.size();
    MELAout << "N plotted: " << nplotted << endl;

    // Plot!
    // Adjust colors
    {
      const unsigned int ncolors = gStyle->GetNumberOfColors();
      const unsigned int stepsize = ncolors/nplotted;
      for (unsigned int ip=0; ip<nplotted; ip++){
        int colorToUse = gStyle->GetColorPalette(ip*stepsize);
        histToPlot.at(ip)->SetMarkerColor(colorToUse);
        histToPlot.at(ip)->SetLineColor(colorToUse);
        histToPlot.at(ip)->SetLineWidth(2);
      }
    }

    TString canvasname = Form(
      "c_%s_%s_%s",
      strACHypo.Data(), strCategory.Data(),
      KDset.at(idim).Data()
    );

    MELAout << "\t- Constructing canvas " << canvasname << endl;

    TCanvas* canvas = new TCanvas(canvasname, "", 8, 30, 800, 800);
    canvas->cd();
    gStyle->SetOptStat(0);
    canvas->SetFillColor(0);
    canvas->SetBorderMode(0);
    canvas->SetBorderSize(2);
    canvas->SetTickx(1);
    canvas->SetTicky(1);
    canvas->SetLeftMargin(0.17);
    canvas->SetRightMargin(0.05);
    canvas->SetTopMargin(0.07);
    canvas->SetBottomMargin(0.13);
    canvas->SetFrameFillStyle(0);
    canvas->SetFrameBorderMode(0);
    canvas->SetFrameFillStyle(0);
    canvas->SetFrameBorderMode(0);

    TLegend* legend = new TLegend(
      0.44,
      0.90-0.10/3.*2.*float(nproctypes),
      0.75,
      0.90
    );
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->SetLineColor(1);
    legend->SetLineStyle(1);
    legend->SetLineWidth(1);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);

    TText* text;
    TPaveText* pt = new TPaveText(0.15, 0.93, 0.85, 1, "brNDC");
    pt->SetBorderSize(0);
    pt->SetFillStyle(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.045);
    text = pt->AddText(0.025, 0.45, "#font[61]{CMS}");
    text->SetTextSize(0.044);
    text = pt->AddText(0.165, 0.42, "#font[52]{Simulation}");
    text->SetTextSize(0.0315);
    TString cErgTev = Form("#font[42]{%.1f fb^{-1} %i TeV}", lumi+0.5, theSqrts);
    text = pt->AddText(0.9, 0.45, cErgTev);
    text->SetTextSize(0.0315);

    bool first=true;
    float minY=9e9, maxY=-9e9;
    for (unsigned int ih=0; ih<nplotted; ih++){
      TH1F*& hh=histToPlot.at(ih);
      float minX=hh->GetXaxis()->GetBinLowEdge(1);
      float maxX=hh->GetXaxis()->GetBinUpEdge(hh->GetNbinsX());
      if (massregion==kOffshell && KDset.at(idim)=="ZZMass") maxX=1500;
      int binXlow = hh->GetXaxis()->FindBin(minX);
      int binXhigh = hh->GetXaxis()->FindBin(maxX);

      hh->GetXaxis()->SetRangeUser(minX, maxX);
      hh->GetXaxis()->SetNdivisions(505);
      hh->GetXaxis()->SetLabelFont(42);
      hh->GetXaxis()->SetLabelOffset(0.007);
      hh->GetXaxis()->SetLabelSize(0.04);
      hh->GetXaxis()->SetTitleSize(0.06);
      hh->GetXaxis()->SetTitleOffset(0.9);
      hh->GetXaxis()->SetTitleFont(42);
      hh->GetYaxis()->SetNdivisions(505);
      hh->GetYaxis()->SetLabelFont(42);
      hh->GetYaxis()->SetLabelOffset(0.007);
      hh->GetYaxis()->SetLabelSize(0.04);
      hh->GetYaxis()->SetTitleSize(0.06);
      hh->GetYaxis()->SetTitleOffset(1.1);
      hh->GetYaxis()->SetTitleFont(42);
      hh->GetYaxis()->SetTitle(ytitle);
      for (int ix=binXlow; ix<binXhigh; ix++){
        float bc = hh->GetBinContent(ix);
        if (bc!=0.){
          minY = std::min(bc, minY);
          maxY = std::max(bc, maxY);
        }
      }
      MELAout << "Min = " << minY << ", max = " << maxY << " after " << hh->GetName() << " (" << hh->GetTitle() << ")" << endl;
      legend->AddEntry(hh, hh->GetTitle(), "l");
      hh->SetTitle("");
    }
    for (unsigned int ih=0; ih<nplotted; ih++){
      TH1F*& hh=histToPlot.at(ih);
      hh->GetYaxis()->SetRangeUser(minY*(minY<0. ? 1.1 : 0.9), maxY*1.3);
      if (first){
        hh->Draw("hist");
        first=false;
      }
      else hh->Draw("histsame");
    }
    legend->Draw("same");
    pt->Draw();
    canvas->RedrawAxis();
    canvas->Modified();
    canvas->Update();
    curdir->WriteTObject(canvas);
    delete pt;
    delete legend;
    canvas->Close();

    for (auto it:proc_hist_SM) delete it.second;
    for (auto it:proc_hist_BSM) delete it.second;

    MELAout << "\t- Done with dimension " << idim << endl;
  } // End loop over dimensions
}

#endif
