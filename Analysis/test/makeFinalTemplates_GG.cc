#ifndef MAKEFINALTEMPLATES_GG_H
#define MAKEFINALTEMPLATES_GG_H

#include "common_includes.h"
#include "CheckSetTemplatesCategoryScheme.h"
#include "acquireProcessMassRatios.cc"


// Constants to affect the template code
#ifndef outputdir_def
#define outputdir_def
const TString user_output_dir = "output/";
#endif


struct MassRatioObject{
  Category const category;
  SystematicVariationTypes const systematic;
  TFile* finput;
  unordered_map<TString, TSpline3*> interpolators;
  MassRatioObject(TString cinput, Category cat, SystematicVariationTypes syst) : category(cat), systematic(syst){
    TDirectory* tmpdir=gDirectory;
    finput = TFile::Open(cinput);
    vector<TSpline3*> tmplist;
    HelperFunctions::extractHistogramsFromDirectory<TSpline3>(finput, tmplist);
    for (TSpline3*& s:tmplist){
      TString sname=s->GetName();
      replaceString(sname, "MassRatio_", ""); // sname becomes name of the tree
      replaceString(sname, "_Smooth_Patched", ""); // sname becomes name of the tree
      interpolators[sname]=s;
    }
    tmpdir->cd(); // Go back to the previous directory
  }
  ~MassRatioObject(){ if (finput) finput->Close(); }
};

void makeFinalTemplates_GG(const Channel channel, const ACHypothesis hypo, const SystematicVariationTypes syst, const unsigned int istage=1, const TString fixedDate=""){
  const ProcessHandler::ProcessType proctype=ProcessHandler::kGG;
  if (channel==NChannels) return;
  GGProcessHandler const* thePerProcessHandle=(GGProcessHandler const*) getOffshellProcessHandler(proctype);
  if (!thePerProcessHandle) return;

  const TString strChannel = getChannelName(channel);
  const TString strACHypo = getACHypothesisName(hypo);
  const TString strStage = Form("Stage%i", istage);
  const TString strSystematics = getSystematicsName(syst);

  // Special variables
  const TString strSMHypo = getACHypothesisName(kSM);
  const TString strSystematics_Nominal = getSystematicsName(sNominal);
  const TString strCategory_Inclusive = getCategoryName(Inclusive);

  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;
  TString cinput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/" + strStage + "/";
  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/FinalTemplates/" + strStage + "/" + strACHypo;
  gSystem->Exec("mkdir -p " + coutput_common);

  TDirectory* rootdir=gDirectory;

  vector<Category> catList = getAllowedCategories(globalCategorizationScheme);
  // Nominal categorization efficiencies
  std::vector<MassRatioObject> CategorizationEfficiencies;
  // Systematic/nominal mass ratios
  std::vector<MassRatioObject> CategorizationSystRatios_POWHEG;
  std::vector<MassRatioObject> CategorizationSystRatios_MCFM;
  for (Category& cat:catList){
    if (!(cat==Inclusive || cat==Untagged)){
      const TString strCategory = getCategoryName(cat);
      TString INPUT_NAME = Form(
        "HtoZZ%s_%s_%s_MassRatiosToNominalInclusive_%s_%s_%s%s",
        strChannel.Data(), strCategory.Data(),
        strACHypo.Data(),
        thePerProcessHandle->getProcessName().Data(),
        strSystematics_Nominal.Data(),
        "POWHEG",
        ".root"
      );
      TString cinput = user_output_dir + sqrtsDir + "Templates/" + strdate + "/MassRatios/" + strStage + "/" + INPUT_NAME;
      if (gSystem->AccessPathName(cinput)){
        acquireMassRatio_ProcessNominalToNominalInclusive(channel, cat, hypo, istage, fixedDate, proctype, "POWHEG");
        if (gSystem->AccessPathName(cinput)){
          MELAerr << "Efficiency file " << cinput << " is not found! Run " << strStage << " functions first." << endl;
          return;
        }
      }
      CategorizationEfficiencies.emplace_back(cinput, cat, sNominal);
    }
    if (syst!=sNominal && systematicAllowed(cat, channel, thePerProcessHandle->getProcessType(), syst, "POWHEG")){
      const TString strCategory = getCategoryName(cat);
      TString INPUT_NAME = Form(
        "HtoZZ%s_%s_%s_MassRatios_SystToNominal_%s_%s_%s%s",
        strChannel.Data(), strCategory.Data(),
        strACHypo.Data(),
        thePerProcessHandle->getProcessName().Data(),
        strSystematics.Data(),
        "POWHEG",
        ".root"
      );
      TString cinput = user_output_dir + sqrtsDir + "Templates/" + strdate + "/MassRatios/" + strStage + "/" + INPUT_NAME;
      if (gSystem->AccessPathName(cinput)){
        acquireMassRatio_ProcessSystToNominal(channel, cat, hypo, syst, istage, fixedDate, proctype, "POWHEG");
        if (gSystem->AccessPathName(cinput)){
          MELAerr << "Systematic ratio file " << cinput << " is not found! Run " << strStage << " functions first." << endl;
          return;
        }
      }
      CategorizationSystRatios_POWHEG.emplace_back(cinput, cat, syst);
    }
    if (syst!=sNominal && systematicAllowed(cat, channel, thePerProcessHandle->getProcessType(), syst, "MCFM")){
      const TString strCategory = getCategoryName(cat);
      TString INPUT_NAME = Form(
        "HtoZZ%s_%s_%s_MassRatios_SystToNominal_%s_%s_%s%s",
        strChannel.Data(), strCategory.Data(),
        strACHypo.Data(),
        thePerProcessHandle->getProcessName().Data(),
        strSystematics.Data(),
        "MCFM",
        ".root"
      );
      TString cinput = user_output_dir + sqrtsDir + "Templates/" + strdate + "/MassRatios/" + strStage + "/" + INPUT_NAME;
      if (gSystem->AccessPathName(cinput)){
        acquireMassRatio_ProcessSystToNominal(channel, cat, hypo, syst, istage, fixedDate, proctype, "MCFM");
        if (gSystem->AccessPathName(cinput)){
          MELAerr << "Systematic ratio file " << cinput << " is not found! Run " << strStage << " functions first." << endl;
          return;
        }
      }
      CategorizationSystRatios_MCFM.emplace_back(cinput, cat, syst);
    }
  }

  // Open nominal inclusive file and get mass shapes
  vector<TString> INPUT_NOMINAL_INCLUSIVE;
  INPUT_NOMINAL_INCLUSIVE.push_back(Form(
    "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_%s%s",
    strChannel.Data(), strCategory_Inclusive.Data(),
    strSMHypo.Data(),
    thePerProcessHandle->getProcessName().Data(),
    strSystematics_Nominal.Data(),
    "MCFM",
    ".root"
  )
  );
  if (hypo!=kSM) INPUT_NOMINAL_INCLUSIVE.push_back(Form(
    "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_%s%s",
    strChannel.Data(), strCategory_Inclusive.Data(),
    strACHypo.Data(),
    thePerProcessHandle->getProcessName().Data(),
    strSystematics_Nominal.Data(),
    "MCFM",
    ".root"
  )
  );

  // Output files
  unordered_map<Category, TFile*, std::hash<int>> foutput;
  for (Category& cat:catList){
    const TString strCategory = getCategoryName(cat);
    TString OUTPUT_NAME = Form(
      "%s/HtoZZ%s_%s_FinalTemplates_%s_%s_%s",
      coutput_common.Data(),
      strChannel.Data(), strCategory.Data(),
      thePerProcessHandle->getProcessName().Data(),
      strSystematics.Data(),
      ".root"
    );
    foutput[cat]=TFile::Open(OUTPUT_NAME, "recreate");
  }
  rootdir->cd(); // Go back to thte root directory

  // Get actual mass distribution from nominal inclusive MCFM sample
  vector<vector<ExtendedHistogram_1D>> hMass_FromNominalInclusive; // [tree][category]
  {
    float weight=0;
    vector<TString> KDset; KDset.push_back("ZZMass");
    unordered_map<TString, float> KDvars;
    for (auto& KDname:KDset) KDvars[KDname]=0;
    vector<ExtendedBinning> KDbinning;
    for (auto& KDname:KDset) KDbinning.push_back(getDiscriminantFineBinning(channel, Inclusive, KDname, true));

    vector<TString> cinputList;
    vector<TFile*> finputList;
    vector<TTree*> treeList;
    for (auto& s:INPUT_NOMINAL_INCLUSIVE) cinputList.push_back(cinput_common + s);
    bool inputValid=true;
    for (auto const& cinput:cinputList){
      if (gSystem->AccessPathName(cinput)){
        MELAerr << "File " << cinput << " is not found! Run " << strStage << " functions first." << endl;
        for (auto& finput:finputList) finput->Close();
        inputValid=false;
        break;
      }
      if (inputValid && cinput!="") finputList.push_back(TFile::Open(cinput, "read"));
      rootdir->cd();
    }
    {
      for (auto& finput:finputList) extractTreesFromDirectory(finput, treeList, false);
      rootdir->cd();
    }

    // Book nominal inclusive histograms
    hMass_FromNominalInclusive.assign(treeList.size(), vector<ExtendedHistogram_1D>());
    for (unsigned int it=0; it<treeList.size(); it++){
      TTree*& tree=treeList.at(it);
      TString treename=tree->GetName();
      MELAout << "Looping over " << treename << " to get mass distributions" << endl;

      for (Category& cat:catList){
        const TString strCategory = getCategoryName(cat);
        TString hname = treename; hname = hname + "_" + strCategory + "_" + strSystematics + "_" + KDset.at(0);
        MELAout << "Setting up mass histogram " << hname << endl;
        hMass_FromNominalInclusive.at(it).emplace_back(hname, hname, KDbinning.at(0));
      }

      tree->SetBranchAddress("weight", &weight);
      for (auto const& KDname:KDset) tree->SetBranchAddress(KDname, &(KDvars[KDname]));
      for (int ev=0; ev<tree->GetEntries(); ev++){
        tree->GetEntry(ev);

        float inclusivenorm=1;
        float one=1;

        for (MassRatioObject& systratio:CategorizationSystRatios_POWHEG){
          if (systratio.category==Inclusive){
            one=systratio.interpolators[treename]->Eval(KDvars[KDset.at(0)]);
            inclusivenorm/=one;
            break;
          }
        }
        for (MassRatioObject& systratio:CategorizationSystRatios_MCFM){
          if (systratio.category==Inclusive){
            inclusivenorm *= systratio.interpolators[treename]->Eval(KDvars[KDset.at(0)]);
            break;
          }
        }

        { // Just to localize icat
          unsigned int icat=0;
          for (Category& cat:catList){
            float extraweight=one;
            if (cat==Untagged){
              for (MassRatioObject& cateff:CategorizationEfficiencies){
                float systadj=1;
                for (MassRatioObject& systratio:CategorizationSystRatios_POWHEG){
                  if (systratio.category==cateff.category){
                    systadj=systratio.interpolators[treename]->Eval(KDvars[KDset.at(0)]);;
                    break;
                  }
                }
                extraweight -= systadj * cateff.interpolators[treename]->Eval(KDvars[KDset.at(0)]);
              }
            }
            else if (cat!=Inclusive){
              extraweight=0;
              for (MassRatioObject& cateff:CategorizationEfficiencies){
                if (cateff.category==cat){
                  extraweight = cateff.interpolators[treename]->Eval(KDvars[KDset.at(0)]);
                  for (MassRatioObject& systratio:CategorizationSystRatios_POWHEG){
                    if (systratio.category==cateff.category){
                      extraweight *= systratio.interpolators[treename]->Eval(KDvars[KDset.at(0)]);;
                      break;
                    }
                  }
                  break;
                }
              }
            }
            if (ev==0) MELAout << "Filling category hist with weight/extraweight/inclusivenorm = " << weight << " / " << extraweight << " / " << inclusivenorm << endl;
            hMass_FromNominalInclusive.at(it).at(icat).fill(KDvars[KDset.at(0)], weight*extraweight*inclusivenorm);
            icat++;
          }
        }
      }
    }

    for (auto& finput:finputList) finput->Close();
  }

  MELAout << "Writing the mass distributions" << endl;
  for (unsigned int icat=0; icat<catList.size(); icat++){
    Category& cat=catList.at(icat);
    foutput[cat]->cd();
    TDirectory* massdistrodir = foutput[cat]->mkdir("MassDistributions"); massdistrodir->cd();
    for (unsigned int it=0; it<hMass_FromNominalInclusive.size(); it++) massdistrodir->WriteTObject(hMass_FromNominalInclusive.at(it).at(icat).getHistogram());
  }
  rootdir->cd();

  MELAout << "Closing the output files" << endl;
  for (Category& cat:catList){ if (foutput[cat]) foutput[cat]->Close(); }
}


#endif
