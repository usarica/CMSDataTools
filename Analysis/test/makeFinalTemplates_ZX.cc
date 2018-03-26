#ifndef MAKEFINALTEMPLATES_ZX_H
#define MAKEFINALTEMPLATES_ZX_H

#include "common_includes.h"
#include "CheckSetTemplatesCategoryScheme.h"
#include "fixTreeWeights.h"
#include "acquireProcessMassRatios.cc"
#include "HistogramSmootherWithGaussianKernel.h"


// Constants to affect the template code
#ifndef outputdir_def
#define outputdir_def
const TString user_output_dir = "output/";
#endif

#ifndef CHISQCUT
#define CHISQCUT 3.841
#endif

// If this is false, errors are over-estimated,
// but this is useful when calculating chi-square over low-statistics samples.
#ifndef USEEFFERRINCOND
#define USEEFFERRINCOND false
#endif

typedef ZXProcessHandler ProcessHandleType;


struct MassRatioObject{
  Category const category;
  SystematicVariationTypes const systematic;
  TFile* finput;
  unordered_map<TString, TSpline3*> interpolators;
  TH2F* hKDRatio2D;
  TH3F* hKDRatio3D;
  MassRatioObject(TString cinput, Category cat, SystematicVariationTypes syst, MassRegion massregion=NMassRegions) : category(cat), systematic(syst), finput(nullptr), hKDRatio2D(nullptr), hKDRatio3D(nullptr){
    MELAout << "Begin MassRatioObject::MassRatioObject" << endl;
    TDirectory* tmpdir=gDirectory;
    finput = TFile::Open(cinput);
    finput->cd();
    vector<TSpline3*> tmplist;
    HelperFunctions::extractHistogramsFromDirectory<TSpline3>(finput, tmplist);
    for (TSpline3*& s:tmplist){
      TString sname=s->GetName();
      replaceString(sname, "MassRatio_", ""); // sname becomes name of the tree
      replaceString(sname, "_Smooth_Patched", ""); // sname becomes name of the tree
      interpolators[sname]=s;
    }
    if (massregion!=NMassRegions){ // Special reweighting for KD shapes, conditional in mass
      TString massregionname = getMassRegionName((CategorizationHelpers::MassRegion) massregion);
      TH2F* hKDRatio2D_tmp;
      TH3F* hKDRatio3D_tmp;
      finput->GetObject(massregionname + "/RatioWithKD", hKDRatio2D_tmp);
      finput->GetObject(massregionname + "/RatioWithKD", hKDRatio3D_tmp);

      if (!hKDRatio2D_tmp && !hKDRatio3D_tmp) MELAerr << "MassRatioObject::MassRatioObject: Acquiring the 'Ratio' histogram did not succeed!" << endl;
      else if (hKDRatio2D_tmp && hKDRatio3D_tmp) MELAerr << "MassRatioObject::MassRatioObject: Acquiring two 'Ratio' histograms! Somthing is wrong..." << endl;
      else if (hKDRatio2D_tmp) MELAout << "MassRatioObject::MassRatioObject: Ratio histogram [" << hKDRatio2D_tmp->GetName() << "] is 2D." << endl;
      else MELAout << "MassRatioObject::MassRatioObject: Ratio histogram [" << hKDRatio3D_tmp->GetName() << "] is 3D." << endl;

      tmpdir->cd(); // Go back to the previous directory
      if (hKDRatio2D_tmp){ hKDRatio2D = new TH2F(*hKDRatio2D_tmp); hKDRatio2D->SetName("RatioWithKD_2D"); }
      if (hKDRatio3D_tmp){ hKDRatio3D = new TH3F(*hKDRatio3D_tmp); hKDRatio3D->SetName("RatioWithKD_3D"); }
    }
    tmpdir->cd(); // Go back to the previous directory
    MELAout << "End MassRatioObject::MassRatioObject" << endl;
  }
  MassRatioObject(MassRatioObject const& other) :
    category(other.category), systematic(other.systematic),
    finput(nullptr),
    interpolators(other.interpolators),
    hKDRatio2D(nullptr), hKDRatio3D(nullptr)
  {
    MELAout << "Begin MassRatioObject::MassRatioObject copy constructor" << endl;
    if (other.hKDRatio2D){ hKDRatio2D=new TH2F(*(other.hKDRatio2D)); hKDRatio2D->SetName(TString(other.hKDRatio2D->GetName()) + "_copy"); }
    if (other.hKDRatio3D){ hKDRatio3D=new TH3F(*(other.hKDRatio3D)); hKDRatio3D->SetName(TString(other.hKDRatio3D->GetName()) + "_copy"); }
    MELAout << "End MassRatioObject::MassRatioObject copy constructor" << endl;
  }
  ~MassRatioObject(){
    MELAout << "Begin ~MassRatioObject::MassRatioObject" << endl;
    if (hKDRatio2D) delete hKDRatio2D;
    if (hKDRatio3D) delete hKDRatio3D;
    if (finput) finput->Close();
    MELAout << "End ~MassRatioObject::MassRatioObject" << endl;
  }
};


template <unsigned char N> void getTemplatesPerCategory(
  TDirectory* rootdir, TFile* foutput,
  ProcessHandleType const*& thePerProcessHandle, Category const& category, ACHypothesis const& hypo,
  std::vector<ProcessHandleType::HypothesisType> const& tplset,
  std::vector<TTree*>& fixedTrees,
  std::vector<TString> const& KDset,
  std::vector<ExtendedBinning> const& KDbinning,
  std::vector<ExtendedHistogram_1D> const& hMass_FromNominalInclusive,
  std::unordered_map<TString, float>& KDvars, float& weight, bool& isCategory
);
template<> void getTemplatesPerCategory<2>(
  TDirectory* rootdir, TFile* foutput,
  ProcessHandleType const*& thePerProcessHandle, Category const& category, ACHypothesis const& hypo,
  std::vector<ProcessHandleType::HypothesisType> const& tplset,
  std::vector<TTree*>& fixedTrees,
  std::vector<TString> const& KDset,
  std::vector<ExtendedBinning> const& KDbinning,
  std::vector<ExtendedHistogram_1D> const& hMass_FromNominalInclusive,
  std::unordered_map<TString, float>& KDvars, float& weight, bool& isCategory
  );
template<> void getTemplatesPerCategory<3>(
  TDirectory* rootdir, TFile* foutput,
  ProcessHandleType const*& thePerProcessHandle, Category const& category, ACHypothesis const& hypo,
  std::vector<ProcessHandleType::HypothesisType> const& tplset,
  std::vector<TTree*>& fixedTrees,
  std::vector<TString> const& KDset,
  std::vector<ExtendedBinning> const& KDbinning,
  std::vector<ExtendedHistogram_1D> const& hMass_FromNominalInclusive,
  std::unordered_map<TString, float>& KDvars, float& weight, bool& isCategory
  );

template <typename T> void PostProcessTemplatesWithPhase(
  TDirectory* rootdir,
  ProcessHandleType const*& thePerProcessHandle, ACHypothesis const& hypo,
  std::vector<ProcessHandleType::HypothesisType> const& tplset,
  std::vector<ExtendedBinning> const& KDbinning,
  std::vector<ExtendedHistogram_1D> const& hMass_FromNominalInclusive,
  std::vector<T>& hTemplates
);
template <> void PostProcessTemplatesWithPhase<ExtendedHistogram_2D>(
  TDirectory* rootdir,
  ProcessHandleType const*& thePerProcessHandle, ACHypothesis const& hypo,
  std::vector<ProcessHandleType::HypothesisType> const& tplset,
  std::vector<ExtendedBinning> const& KDbinning,
  std::vector<ExtendedHistogram_1D> const& hMass_FromNominalInclusive,
  std::vector<ExtendedHistogram_2D>& hTemplates
  );
template <> void PostProcessTemplatesWithPhase<ExtendedHistogram_3D>(
  TDirectory* rootdir,
  ProcessHandleType const*& thePerProcessHandle, ACHypothesis const& hypo,
  std::vector<ProcessHandleType::HypothesisType> const& tplset,
  std::vector<ExtendedBinning> const& KDbinning,
  std::vector<ExtendedHistogram_1D> const& hMass_FromNominalInclusive,
  std::vector<ExtendedHistogram_3D>& hTemplates
  );

void getControl2DXSlices(
  TDirectory* rootdir, TFile* foutput,
  ProcessHandleType const*& thePerProcessHandle,
  std::vector<ExtendedBinning> const& KDbinning,
  std::vector<ExtendedHistogram_3D> const& hTemplates
);


bool getFilesAndTrees(
  const Channel channel, const Category category, const SystematicVariationTypes syst,
  const unsigned int istage, const TString fixedDate,
  ProcessHandler const* inputProcessHandle,
  const TString strGenerator,
  std::vector<TFile*>& finputList, std::vector<TTree*>& treeList,
  bool forceUseNominal=false
){
  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);
  const TString strStage = Form("Stage%i", istage);
  const TString strSystematics = getSystematicsName(syst);
  const TString strSystematics_Nominal = getSystematicsName(sNominal);
  TString strSystematics_effective=strSystematics;
  if (forceUseNominal) strSystematics_effective=strSystematics_Nominal;

  vector<ZXFakeRateHandler::FakeRateMethod> FRMethods;
  FRMethods.push_back(ZXFakeRateHandler::mSS);

  // Setup the output directories
  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  TString cinput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/" + strStage + "/";

  vector<TString> INPUT_NAME;
  for (ZXFakeRateHandler::FakeRateMethod& FRMethod:FRMethods){
    const TString FRMethodName = ZXFakeRateHandler::TranslateFakeRateMethodToString(FRMethod);
    INPUT_NAME.push_back(Form(
      "HtoZZ%s_%s_FinalTemplates_%s_%s_%s_%s%s",
      strChannel.Data(), strCategory.Data(),
      inputProcessHandle->getProcessName().Data(), FRMethodName.Data(),
      strSystematics_effective.Data(),
      strGenerator.Data(),
      ".root"
    )
    );
  }

  for (auto& s:INPUT_NAME){
    TString cinput=cinput_common + s;
    if (gSystem->AccessPathName(cinput)){
      if (category==Untagged){
        const Category category_inc = Inclusive;
        const TString strCategory_inc = getCategoryName(category_inc);
        MELAout << "getFilesAndTrees::File " << cinput << " is not found for untagged category! Attempting to substitute inclusive category..." << endl;
        HelperFunctions::replaceString<TString, const TString>(cinput, strCategory, strCategory_inc);
      }
    }
    if (gSystem->AccessPathName(cinput)){
      MELAerr << "getFilesAndTrees::File " << cinput << " is not found! Run " << strStage << " functions first." << endl;
      treeList.clear();
      return false;
    }
    if (cinput!=""){
      TFile* finput = TFile::Open(cinput, "read");
      if (finput){
        if (!finput->IsZombie()){
          extractTreesFromDirectory(finput, treeList);
          finputList.push_back(finput);
        }
        else if (finput->IsOpen()){
          finput->Close();
          treeList.clear();
          return false;
        }
        else{
          treeList.clear();
          return false;
        }
      }
    }
  }
  return true;
}


void makeFinalTemplates_ZX(const Channel channel, const ACHypothesis hypo, const SystematicVariationTypes syst, CategorizationHelpers::MassRegion massregion, const unsigned int istage=1, const TString fixedDate=""){
  const ProcessHandler::ProcessType proctype=ProcessHandler::kZX;
  if (channel==NChannels) return;
  ProcessHandleType const* inputProcessHandle=(ProcessHandleType const*) getProcessHandlerPerMassRegion(proctype, CategorizationHelpers::kOffshell); // Input is always organized in offshell conventions
  ProcessHandleType const* outputProcessHandle=(ProcessHandleType const*) getProcessHandlerPerMassRegion(proctype, massregion);
  if (!inputProcessHandle || !outputProcessHandle) return;

  vector<Category> catList = getAllowedCategories(globalCategorizationScheme);
  {
    bool doProceed=false;
    for (auto& cat:catList) doProceed |= (systematicAllowed(cat, channel, proctype, syst, "Data"));
    if (!doProceed) return;
  }
  bool needsKDreweighting = false;

  const TString strChannel = getChannelName(channel);
  const TString strACHypo = getACHypothesisName(hypo);
  const TString strStage = Form("Stage%i", istage);
  const TString strSystematics = getSystematicsName(syst);

  // Special variables
  const TString strSMHypo = getACHypothesisName(kSM);
  const TString strSystematics_Nominal = getSystematicsName(sNominal);
  const TString strCategory_Inclusive = getCategoryName(Inclusive);

  vector<ProcessHandleType::HypothesisType> tplset = inputProcessHandle->getHypotheses();

  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;
  TString cinput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/" + strStage + "/";
  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/FinalTemplates/" + strStage + "/" + strACHypo;
  gSystem->Exec("mkdir -p " + coutput_common);

  TDirectory* rootdir=gDirectory;

  TString OUTPUT_LOG_NAME = Form(
    "%s/HtoZZ%s_%s_FinalTemplates_%s_%s_%s",
    coutput_common.Data(),
    strChannel.Data(), "AllCategories",
    outputProcessHandle->getProcessName().Data(),
    strSystematics.Data(),
    ".log"
  );
  MELAout.open(OUTPUT_LOG_NAME.Data());

  // Nominal categorization efficiencies
  std::vector<MassRatioObject> CategorizationEfficiencies;
  // Systematic/nominal mass ratios
  std::vector<MassRatioObject> CategorizationSystRatios;
  for (Category& cat:catList){
    if (!(cat==Inclusive || cat==Untagged)){
      const TString strCategory = getCategoryName(cat);
      TString INPUT_NAME = Form(
        "HtoZZ%s_%s_%s_MassRatiosToNominalInclusive_%s_%s_%s%s",
        strChannel.Data(), strCategory.Data(),
        strACHypo.Data(),
        inputProcessHandle->getProcessName().Data(),
        strSystematics_Nominal.Data(),
        "Data",
        ".root"
      );
      TString cinput = user_output_dir + sqrtsDir + "Templates/" + strdate + "/MassRatios/" + strStage + "/" + INPUT_NAME;
      if (gSystem->AccessPathName(cinput)){
        acquireMassRatio_ProcessNominalToNominalInclusive_one(channel, cat, hypo, istage, fixedDate, proctype, "Data");
        if (gSystem->AccessPathName(cinput)){
          MELAerr << "Efficiency file " << cinput << " is not found! Run " << strStage << " functions first." << endl;
          return;
        }
      }
      CategorizationEfficiencies.emplace_back(cinput, cat, sNominal);
    }
    if (syst!=sNominal && cat!=Untagged && systematicAllowed(cat, channel, inputProcessHandle->getProcessType(), syst, "Data")){
      const TString strCategory = getCategoryName(cat);
      TString INPUT_NAME = Form(
        "HtoZZ%s_%s_%s_MassRatios_SystToNominal_%s_%s_%s%s",
        strChannel.Data(), strCategory.Data(),
        strACHypo.Data(),
        inputProcessHandle->getProcessName().Data(),
        strSystematics.Data(),
        "Data",
        ".root"
      );
      TString cinput = user_output_dir + sqrtsDir + "Templates/" + strdate + "/MassRatios/" + strStage + "/" + INPUT_NAME;
      if (gSystem->AccessPathName(cinput)){
        acquireMassRatio_ProcessSystToNominal_one(channel, cat, hypo, syst, istage, fixedDate, proctype, "Data");
        if (gSystem->AccessPathName(cinput)){
          MELAerr << "Systematic ratio file " << cinput << " is not found! Run " << strStage << " functions first." << endl;
          return;
        }
      }
      CategorizationSystRatios.emplace_back(cinput, cat, syst);
    }
  }

  // Open nominal inclusive file and get mass shapes
  vector<TString> INPUT_NOMINAL_INCLUSIVE;
  vector<ZXFakeRateHandler::FakeRateMethod> FRMethods;
  FRMethods.push_back(ZXFakeRateHandler::mSS);
  for (ZXFakeRateHandler::FakeRateMethod& FRMethod:FRMethods){
    const TString FRMethodName = ZXFakeRateHandler::TranslateFakeRateMethodToString(FRMethod);
    INPUT_NOMINAL_INCLUSIVE.push_back(Form(
      "HtoZZ%s_%s_FinalTemplates_%s_%s_%s_%s%s",
      strChannel.Data(), strCategory_Inclusive.Data(),
      inputProcessHandle->getProcessName().Data(), FRMethodName.Data(),
      strSystematics_Nominal.Data(),
      "Data",
      ".root"
    )
    );
  }

  // Setup colors
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

  // Output files
  unordered_map<Category, TFile*, std::hash<int>> foutput;
  for (Category& cat:catList){
    const TString strCategory = getCategoryName(cat);
    const TString strSystematicsOutput = getSystematicsCombineName(cat, channel, proctype, syst);
    TString OUTPUT_NAME = Form(
      "%s/HtoZZ%s_%s_FinalTemplates_%s_%s_%s",
      coutput_common.Data(),
      strChannel.Data(), strCategory.Data(),
      outputProcessHandle->getProcessName().Data(),
      strSystematicsOutput.Data(),
      ".root"
    );
    foutput[cat]=TFile::Open(OUTPUT_NAME, "recreate");
  }
  rootdir->cd(); // Go back to thte root directory

  // Get actual mass distribution from nominal inclusive sample
  unordered_map<Category, vector<ExtendedHistogram_1D>, std::hash<int>> hMass_FromNominalInclusive; // [tree]
  {
    float weight=0;
    vector<TString> KDset; KDset.push_back("ZZMass");
    unordered_map<TString, float> KDvars;
    for (auto& KDname:KDset) KDvars[KDname]=0;
    vector<ExtendedBinning> KDbinning;
    for (auto& KDname:KDset) KDbinning.push_back(getDiscriminantFineBinning(channel, Inclusive, KDname, massregion));

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
    for (Category& cat:catList) hMass_FromNominalInclusive[cat]=vector<ExtendedHistogram_1D>();
    for (unsigned int it=0; it<treeList.size(); it++){
      TTree*& tree=treeList.at(it);
      TString treename=tree->GetName();
      MELAout << "Looping over " << treename << " to get mass distributions" << endl;

      for (Category& cat:catList){
        const TString strCategory = getCategoryName(cat);
        TString hname = outputProcessHandle->getOutputTreeName();
        hname = hname + "_" + strCategory + "_" + strSystematics + "_" + KDset.at(0);
        MELAout << "Setting up mass histogram " << hname << endl;
        hMass_FromNominalInclusive[cat].emplace_back(hname, hname, KDbinning.at(0));
      }

      tree->SetBranchStatus("*", 0);
      bookBranch(tree, "weight", &weight);
      for (auto const& KDname:KDset) bookBranch(tree, KDname, &(KDvars[KDname]));
      int nEntries=tree->GetEntries();
      for (int ev=0; ev<tree->GetEntries(); ev++){
        tree->GetEntry(ev);
        progressbar(ev, nEntries);

        float inclusivenorm=1;
        float one=1;

        for (MassRatioObject& systratio:CategorizationSystRatios){
          if (systratio.category==Inclusive){
            one=systratio.interpolators[treename]->Eval(KDvars[KDset.at(0)]);
            break;
          }
        }

        for (Category& cat:catList){
          float extraweight=one;
          if (cat==Untagged){
            for (MassRatioObject& cateff:CategorizationEfficiencies){
              float systadj=1;
              for (MassRatioObject& systratio:CategorizationSystRatios){
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
                for (MassRatioObject& systratio:CategorizationSystRatios){
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
          hMass_FromNominalInclusive[cat].at(it).fill(KDvars[KDset.at(0)], weight*extraweight*inclusivenorm);
        }
      } // End loop over tree events
    } // End loop over trees

    for (Category& cat:catList){ // Check integrity of mass histograms
      MELAout << "Checking integrity of mass histograms for category " << getCategoryName(cat) << endl;
      for (ExtendedHistogram_1D const& ehmass:hMass_FromNominalInclusive[cat]){
        if (checkHistogramIntegrity(ehmass.getHistogram())) MELAout << "Integrity of " << ehmass.getName() << " is GOOD." << endl;
        else MELAout << "WARNING: Integrity of " << ehmass.getName() << " is BAD." << endl;
      }
    }

    for (auto& finput:finputList) finput->Close();
  }

  gStyle->SetOptStat(0);

  MELAout << "Writing the mass distributions" << endl;
  for (Category& cat:catList){
    foutput[cat]->cd();
    TDirectory* massdistrodir = foutput[cat]->mkdir("MassDistributions"); massdistrodir->cd();
    for (unsigned int it=0; it<hMass_FromNominalInclusive[cat].size(); it++) massdistrodir->WriteTObject(hMass_FromNominalInclusive[cat].at(it).getHistogram());
  }
  rootdir->cd();

  // Move on to create the conditional templates
  MELAout << "Progressing to create the conditional templates" << endl;
  for (Category& cat:catList){
    const TString strCategory = getCategoryName(cat);
    MELAout << "\t- Begin category " << strCategory << endl;

    float weight=0;
    bool isCategory=(cat==Inclusive);
    TString catFlagName="";
    if (!isCategory) catFlagName = TString("is_") + strCategory + TString("_") + strACHypo;
    vector<TString> KDset; KDset.push_back("ZZMass"); { vector<TString> KDset2=getACHypothesisKDNameSet(hypo, cat, massregion); appendVector(KDset, KDset2); }
    unordered_map<TString, float> KDvars;
    for (auto& KDname:KDset) KDvars[KDname]=0;
    vector<ExtendedBinning> KDbinning;
    for (auto& KDname:KDset) KDbinning.push_back(getDiscriminantFineBinning(channel, cat, KDname, massregion));
    unsigned int nKDs = KDset.size();
    MELAout << "\t- Number of template dimensions = " << nKDs << endl;

    vector<TFile*> finputs;
    vector<TTree*> treeList;
    vector<TTree*> fixedTrees;
    MELAout << "\t- Obtaining samples..." << endl;
    bool success = getFilesAndTrees(
      channel, cat, syst,
      istage, fixedDate, inputProcessHandle, "Data",
      finputs, treeList,
      false
    );
    MELAout << "\t-- " << (success ? "Success!" : "failure!") << endl;

    rootdir->cd();

    MELAout << "\t- Processing templates for category " << strCategory << endl;
    MassRatioObject* systratio_KDfix=nullptr;
    if (needsKDreweighting){
      for (MassRatioObject& systratio:CategorizationSystRatios){
        if (systratio.category==cat){
          if (systratio.hKDRatio2D || systratio.hKDRatio3D){
            systratio_KDfix=&systratio;
            if (systratio_KDfix->hKDRatio2D) MELAout << "systratio_KDfix has 2D histogram " << systratio_KDfix->hKDRatio2D->GetName() << endl;
            if (systratio_KDfix->hKDRatio3D) MELAout << "systratio_KDfix has 3D histogram " << systratio_KDfix->hKDRatio3D->GetName() << endl;
          }
        }
      }
    }

    if (success){
      MELAout << "\t- Fixing Zjets tree weights" << endl;
      for (TTree*& tree:treeList){
        tree->SetBranchStatus("*", 0);
        bookBranch(tree, "weight", &weight);
        for (auto& KDname:KDset) bookBranch(tree, KDname, &(KDvars[KDname]));
        if (catFlagName!=""){
          bookBranch(tree, catFlagName, &isCategory);
          if (!branchExists(tree, catFlagName)) isCategory=true;
        }

        // Fix KD shapes for Pythia or MINLO systematics
        if (systratio_KDfix){
          TTree* tmptree=nullptr;
          switch (nKDs){
          case 2:
            tmptree=fixTreeWeights(systratio_KDfix->hKDRatio2D, tree, KDvars.find(KDset.at(0))->second, KDvars.find(KDset.at(1))->second, weight, isCategory);
            break;
          case 3:
            tmptree=fixTreeWeights(systratio_KDfix->hKDRatio3D, tree, KDvars.find(KDset.at(0))->second, KDvars.find(KDset.at(1))->second, KDvars.find(KDset.at(2))->second, weight, isCategory);
            break;
          }
          if (tmptree){
            tree=tmptree;

            // Need to re-book branches
            bookBranch(tree, "weight", &weight);
            for (auto& KDname:KDset) bookBranch(tree, KDname, &(KDvars.find(KDname)->second));
            if (catFlagName!="") bookBranch(tree, catFlagName, &isCategory);
          }
        }

        fixedTrees.push_back(tree);
      }
      MELAout << "\t- Tree weights fixed" << endl;
    }
    rootdir->cd();

    MELAout << "\t- Attempting templates..." << endl;
    if (nKDs==2) getTemplatesPerCategory<2>(
      rootdir, foutput[cat], outputProcessHandle, cat, hypo,
      tplset,
      fixedTrees,
      KDset, KDbinning,
      hMass_FromNominalInclusive.find(cat)->second,
      KDvars, weight, isCategory
      );
    else if (nKDs==3) getTemplatesPerCategory<3>(
      rootdir, foutput[cat], outputProcessHandle, cat, hypo,
      tplset,
      fixedTrees,
      KDset, KDbinning,
      hMass_FromNominalInclusive.find(cat)->second,
      KDvars, weight, isCategory
      );

    rootdir->cd();
    MELAout << "\t- Templates obtained successfully. Cleaning up..." << endl;
    for (TTree*& tree:fixedTrees) delete tree;
    for (TFile*& finput:finputs) finput->Close();
    MELAout << "\t- Cleanup done." << endl;
    rootdir->cd();
  }

  MELAout << "Closing the output files" << endl;
  for (Category& cat:catList){
    if (foutput[cat]) foutput[cat]->Close();
    else MELAerr << "Something went wrong. Category " << cat << " file is invalid!" << endl;
  }
  MELAout.close();
}

template<> void getTemplatesPerCategory<2>(
  TDirectory* rootdir, TFile* foutput,
  ProcessHandleType const*& thePerProcessHandle, Category const& category, ACHypothesis const& hypo,
  std::vector<ProcessHandleType::HypothesisType> const& tplset,
  std::vector<TTree*>& fixedTrees,
  std::vector<TString> const& KDset,
  std::vector<ExtendedBinning> const& KDbinning,
  std::vector<ExtendedHistogram_1D> const& hMass_FromNominalInclusive,
  std::unordered_map<TString, float>& KDvars, float& weight, bool& isCategory
  ){
  typedef ExtendedHistogram_2D ExtHist_t;
  typedef TH2F TH_t;

  if (fixedTrees.empty()) return;
  const unsigned int ntpls=tplset.size();
  const unsigned int nKDs=KDbinning.size();
  assert(fixedTrees.size()==ntpls);
  assert(nKDs==2);

  // Fill templates
  vector<ExtHist_t> hTemplates;
  hTemplates.reserve(ntpls);
  for (unsigned int t=0; t<ntpls; t++){
    TTree*& tree=fixedTrees.at(t);
    ProcessHandleType::HypothesisType const& treetype = tplset.at(t);
    ProcessHandleType::TemplateType tpltype = ProcessHandleType::castIntToTemplateType(ProcessHandleType::castHypothesisTypeToInt(treetype));
    TString templatename = thePerProcessHandle->getTemplateName() + "";
    TString templatetitle = thePerProcessHandle->getProcessLabel();

    hTemplates.emplace_back(templatename, templatetitle, KDbinning.at(0), KDbinning.at(1));
    TH_t* hSmooth=getSmoothHistogram(
      templatename+"_Smooth", "", KDbinning.at(0), KDbinning.at(1),
      tree, KDvars.find(KDset.at(0))->second, KDvars.find(KDset.at(1))->second, weight, isCategory,
      2, 2
    );
    *(hTemplates.back().getHistogram()) = *hSmooth;
    delete hSmooth;
    hTemplates.back().getHistogram()->SetNameTitle(templatename, templatetitle);
  }
  // Do post-processing
  PostProcessTemplatesWithPhase(
    rootdir,
    thePerProcessHandle, hypo,
    tplset,
    KDbinning,
    hMass_FromNominalInclusive,
    hTemplates
  );

  if (hTemplates.empty()) assert(0);

  // Multiply with mass histogram only after combination
  for (unsigned int t=0; t<ntpls; t++){
    auto& tpl = hTemplates.at(t);
    TH_t*& htpl = tpl.getHistogram();

    ProcessHandleType::HypothesisType const& treetype = tplset.at(t);
    ProcessHandleType::TemplateType tpltype = ProcessHandleType::castIntToTemplateType(ProcessHandleType::castHypothesisTypeToInt(treetype));

    TH1F const* hMass=hMass_FromNominalInclusive.at(t).getHistogram();
    assert(hMass->GetNbinsX()==htpl->GetNbinsX());
    multiplyHistograms(htpl, hMass, 0, htpl, USEEFFERRINCOND);

    double integralerror=0;
    double integral = getHistogramIntegralAndError(htpl, 1, htpl->GetNbinsX(), 1, htpl->GetNbinsY(), false, &integralerror);
    MELAout << "Integral [ " << tpl.getName() << " ] after renormalizing mass: " << integral << " +- " << integralerror << endl;
    MELAout << "Checking integrity of [ " << tpl.getName() << " ]" << endl;
    if (checkHistogramIntegrity(htpl)) MELAout << "Integrity of [ " << tpl.getName() << " ] is GOOD." << endl;
    else MELAout << "WARNING: Integrity of [ " << tpl.getName() << " ] is BAD." << endl;
  }

  {
    foutput->cd();
    for (unsigned int t=0; t<hTemplates.size(); t++){
      auto& tpl = hTemplates.at(t);
      TH_t*& htpl = tpl.getHistogram();

      TString tplname=tpl.getName(); TString tpltitle=tpl.getTitle(); replaceString(tplname, "_POWHEG", ""); replaceString(tplname, "_Data", ""); tpl.setNameTitle(tplname, tpltitle);
      ProcessHandleType::TemplateType tpltype = ProcessHandleType::castIntToTemplateType(t);

      MELAout << "Checking integrity of [ " << htpl->GetName() << " ]" << endl;
      if (checkHistogramIntegrity(htpl)) MELAout << "Integrity of [ " << htpl->GetName() << " ] is GOOD." << endl;
      else MELAout << "WARNING: Integrity of [ " << htpl->GetName() << " ] is BAD." << endl;

      doTemplatePostprocessing(htpl);
      double integralerror=0;
      double integral = getHistogramIntegralAndError(htpl, 1, htpl->GetNbinsX(), 1, htpl->GetNbinsY(), true, &integralerror);
      MELAout << "Integral [ " << htpl->GetName() << " ] before writing: " << integral << " +- " << integralerror << endl;
      // SM on-shell analysis uses conditional templates
      if (hypo==kSM && thePerProcessHandle->getProcessMassRegion()==kOnshell) conditionalizeHistogram<TH_t>(htpl, 0, nullptr, true, USEEFFERRINCOND);

      MELAout << "Final integrity check on [ " << htpl->GetName() << " ]" << endl;
      if (checkHistogramIntegrity(htpl)) MELAout << "Integrity of [ " << htpl->GetName() << " ] is GOOD." << endl;
      else MELAout << "WARNING: Integrity of [ " << htpl->GetName() << " ] is BAD." << endl;
      foutput->WriteTObject(htpl);
    }
    rootdir->cd();
  }
}
template<> void getTemplatesPerCategory<3>(
  TDirectory* rootdir, TFile* foutput,
  ProcessHandleType const*& thePerProcessHandle, Category const& category, ACHypothesis const& hypo,
  std::vector<ProcessHandleType::HypothesisType> const& tplset,
  std::vector<TTree*>& fixedTrees,
  std::vector<TString> const& KDset,
  std::vector<ExtendedBinning> const& KDbinning,
  std::vector<ExtendedHistogram_1D> const& hMass_FromNominalInclusive,
  std::unordered_map<TString, float>& KDvars, float& weight, bool& isCategory
  ){
  typedef ExtendedHistogram_3D ExtHist_t;
  typedef TH3F TH_t;

  if (fixedTrees.empty()) return;
  const unsigned int ntpls=tplset.size();
  const unsigned int nKDs=KDbinning.size();
  assert(fixedTrees.size()==ntpls);
  assert(nKDs==3);

  // Fill templates
  vector<ExtHist_t> hTemplates;
  hTemplates.reserve(ntpls);
  for (unsigned int t=0; t<ntpls; t++){
    TTree*& tree=fixedTrees.at(t);
    ProcessHandleType::HypothesisType const& treetype = tplset.at(t);
    ProcessHandleType::TemplateType tpltype = ProcessHandleType::castIntToTemplateType(ProcessHandleType::castHypothesisTypeToInt(treetype));
    TString templatename = thePerProcessHandle->getTemplateName() + "";
    TString templatetitle = thePerProcessHandle->getProcessLabel();

    hTemplates.emplace_back(templatename, templatetitle, KDbinning.at(0), KDbinning.at(1), KDbinning.at(2));
    TH_t* hSmooth=getSmoothHistogram(
      templatename+"_Smooth", "", KDbinning.at(0), KDbinning.at(1), KDbinning.at(2),
      tree, KDvars.find(KDset.at(0))->second, KDvars.find(KDset.at(1))->second, KDvars.find(KDset.at(2))->second, weight, isCategory,
      2, 2, 2
    );
    *(hTemplates.back().getHistogram()) = *hSmooth;
    delete hSmooth;
    hTemplates.back().getHistogram()->SetNameTitle(templatename, templatetitle);
  }
  // Do post-processing
  PostProcessTemplatesWithPhase(
    rootdir,
    thePerProcessHandle, hypo,
    tplset,
    KDbinning,
    hMass_FromNominalInclusive,
    hTemplates
  );

  assert(!hTemplates.empty());

  // Multiply with mass histogram only after combination
  for (unsigned int t=0; t<ntpls; t++){
    auto& tpl = hTemplates.at(t);
    TH_t*& htpl = tpl.getHistogram();

    ProcessHandleType::HypothesisType const& treetype = tplset.at(t);
    ProcessHandleType::TemplateType tpltype = ProcessHandleType::castIntToTemplateType(ProcessHandleType::castHypothesisTypeToInt(treetype));

    TH1F const* hMass=hMass_FromNominalInclusive.at(t).getHistogram();
    assert(hMass->GetNbinsX()==htpl->GetNbinsX());
    multiplyHistograms(htpl, hMass, 0, htpl, USEEFFERRINCOND);

    double integralerror=0;
    double integral = getHistogramIntegralAndError(htpl, 1, htpl->GetNbinsX(), 1, htpl->GetNbinsY(), 1, htpl->GetNbinsZ(), false, &integralerror);
    MELAout << "Integral [ " << tpl.getName() << " ] after renormalizing mass: " << integral << " +- " << integralerror << endl;
    MELAout << "Checking integrity of [ " << tpl.getName() << " ]" << endl;
    if (checkHistogramIntegrity(htpl)) MELAout << "Integrity of [ " << tpl.getName() << " ] is GOOD." << endl;
    else MELAout << "WARNING: Integrity of [ " << tpl.getName() << " ] is BAD." << endl;
  }

  // Recombine all templates
  {
    foutput->cd();
    for (unsigned int t=0; t<hTemplates.size(); t++){
      auto& tpl = hTemplates.at(t);
      TH_t*& htpl = tpl.getHistogram();

      TString tplname=tpl.getName(); TString tpltitle=tpl.getTitle(); replaceString(tplname, "_POWHEG", ""); replaceString(tplname, "_Data", ""); tpl.setNameTitle(tplname, tpltitle);
      ProcessHandleType::TemplateType tpltype = ProcessHandleType::castIntToTemplateType(t);

      MELAout << "Checking integrity of [ " << htpl->GetName() << " ]" << endl;
      if (checkHistogramIntegrity(htpl)) MELAout << "Integrity of [ " << htpl->GetName() << " ] is GOOD." << endl;
      else MELAout << "WARNING: Integrity of [ " << htpl->GetName() << " ] is BAD." << endl;

      doTemplatePostprocessing(htpl);
      double integralerror=0;
      double integral = getHistogramIntegralAndError(htpl, 1, htpl->GetNbinsX(), 1, htpl->GetNbinsY(), 1, htpl->GetNbinsZ(), true, &integralerror);
      MELAout << "Integral [ " << htpl->GetName() << " ] before writing: " << integral << " +- " << integralerror << endl;
      // SM on-shell analysis uses conditional templates
      if (hypo==kSM && thePerProcessHandle->getProcessMassRegion()==kOnshell) conditionalizeHistogram<TH_t>(htpl, 0, nullptr, true, USEEFFERRINCOND);

      MELAout << "Final integrity check on [ " << htpl->GetName() << " ]" << endl;
      if (checkHistogramIntegrity(htpl)) MELAout << "Integrity of [ " << htpl->GetName() << " ] is GOOD." << endl;
      else MELAout << "WARNING: Integrity of [ " << htpl->GetName() << " ] is BAD." << endl;
      foutput->WriteTObject(htpl);
    }
    rootdir->cd();
    getControl2DXSlices(rootdir, foutput, thePerProcessHandle, KDbinning, hTemplates);
  }
}


template <> void PostProcessTemplatesWithPhase<ExtendedHistogram_2D>(
  TDirectory* rootdir,
  ProcessHandleType const*& thePerProcessHandle, ACHypothesis const& hypo,
  std::vector<ProcessHandleType::HypothesisType> const& tplset,
  std::vector<ExtendedBinning> const& KDbinning,
  std::vector<ExtendedHistogram_1D> const& hMass_FromNominalInclusive,
  std::vector<ExtendedHistogram_2D>& hTemplates
  ){
  typedef TH2F TH_t;

  if (!hTemplates.empty()){
    rootdir->cd();
    unsigned int ntpls=hTemplates.size();
    {
      vector<TH_t*> hTemplateObjects;
      for (auto& tpl:hTemplates){
        TH_t*& htpl = tpl.getHistogram();
        MELAout << "Integral [ " << tpl.getName() << " ] before recombineHistogramsToTemplates: " << getHistogramIntegralAndError(htpl, 1, htpl->GetNbinsX(), 1, htpl->GetNbinsY(), false, nullptr) << endl;
        hTemplateObjects.push_back(htpl);

        MELAout << "Checking integrity of [ " << htpl->GetName() << " ]" << endl;
        if (checkHistogramIntegrity(htpl)) MELAout << "Integrity of [ " << htpl->GetName() << " ] is GOOD." << endl;
        else MELAout << "WARNING: Integrity of [ " << htpl->GetName() << " ] is BAD." << endl;
      }
      thePerProcessHandle->recombineHistogramsToTemplates(hTemplateObjects);
    }
    for (unsigned int t=0; t<ntpls; t++){
      auto& tpl = hTemplates.at(t);
      TProfile* prof_x = tpl.getProfileX();
      TProfile* prof_y = tpl.getProfileY();
      TH_t*& htpl = tpl.getHistogram();

      ProcessHandleType::HypothesisType const& treetype = tplset.at(t);
      ProcessHandleType::TemplateType tpltype = ProcessHandleType::castIntToTemplateType(ProcessHandleType::castHypothesisTypeToInt(treetype));

      vector<pair<TProfile const*, unsigned int>> condProfs; condProfs.push_back(pair<TProfile const*, unsigned int>(prof_x, 0));
      conditionalizeHistogram<TH_t>(htpl, 0, nullptr, false, USEEFFERRINCOND);

      MELAout << "Integral [ " << tpl.getName() << " ] after post-processing function: " << getHistogramIntegralAndError(htpl, 1, htpl->GetNbinsX(), 1, htpl->GetNbinsY(), false, nullptr) << endl;
      MELAout << "Checking integrity of [ " << htpl->GetName() << " ]" << endl;
      if (checkHistogramIntegrity(htpl)) MELAout << "Integrity of [ " << htpl->GetName() << " ] is GOOD." << endl;
      else MELAout << "WARNING: Integrity of [ " << htpl->GetName() << " ] is BAD." << endl;
    }
  }
}
template <> void PostProcessTemplatesWithPhase<ExtendedHistogram_3D>(
  TDirectory* rootdir,
  ProcessHandleType const*& thePerProcessHandle, ACHypothesis const& hypo,
  std::vector<ProcessHandleType::HypothesisType> const& tplset,
  std::vector<ExtendedBinning> const& KDbinning,
  std::vector<ExtendedHistogram_1D> const& hMass_FromNominalInclusive,
  std::vector<ExtendedHistogram_3D>& hTemplates
  ){
  typedef TH3F TH_t;

  if (!hTemplates.empty()){
    rootdir->cd();
    unsigned int ntpls=hTemplates.size();
    {
      vector<TH_t*> hTemplateObjects;
      for (auto& tpl:hTemplates){
        TH_t*& htpl = tpl.getHistogram();
        MELAout << "Integral [ " << tpl.getName() << " ] before recombineHistogramsToTemplates: " << getHistogramIntegralAndError(htpl, 1, htpl->GetNbinsX(), 1, htpl->GetNbinsY(), 1, htpl->GetNbinsZ(), false, nullptr) << endl;
        hTemplateObjects.push_back(htpl);

        MELAout << "Checking integrity of [ " << htpl->GetName() << " ]" << endl;
        if (checkHistogramIntegrity(htpl)) MELAout << "Integrity of [ " << htpl->GetName() << " ] is GOOD." << endl;
        else MELAout << "WARNING: Integrity of [ " << htpl->GetName() << " ] is BAD." << endl;
      }
      thePerProcessHandle->recombineHistogramsToTemplates(hTemplateObjects);
    }
    for (unsigned int t=0; t<ntpls; t++){
      auto& tpl = hTemplates.at(t);
      TProfile* prof_x = tpl.getProfileX();
      TProfile* prof_y = tpl.getProfileY();
      TProfile* prof_z = tpl.getProfileZ();
      TH_t*& htpl = tpl.getHistogram();

      ProcessHandleType::HypothesisType const& treetype = tplset.at(t);
      ProcessHandleType::TemplateType tpltype = ProcessHandleType::castIntToTemplateType(ProcessHandleType::castHypothesisTypeToInt(treetype));

      vector<pair<TProfile const*, unsigned int>> condProfs; condProfs.push_back(pair<TProfile const*, unsigned int>(prof_x, 0));
      conditionalizeHistogram<TH_t>(htpl, 0, nullptr, false, USEEFFERRINCOND);

      MELAout << "Integral [ " << tpl.getName() << " ] after post-processing function: " << getHistogramIntegralAndError(htpl, 1, htpl->GetNbinsX(), 1, htpl->GetNbinsY(), 1, htpl->GetNbinsZ(), false, nullptr) << endl;
      MELAout << "Checking integrity of [ " << htpl->GetName() << " ]" << endl;
      if (checkHistogramIntegrity(htpl)) MELAout << "Integrity of [ " << htpl->GetName() << " ] is GOOD." << endl;
      else MELAout << "WARNING: Integrity of [ " << htpl->GetName() << " ] is BAD." << endl;
    }
  }
}

void getControl2DXSlices(
  TDirectory* rootdir, TFile* foutput,
  ProcessHandleType const*& thePerProcessHandle,
  std::vector<ExtendedBinning> const& KDbinning,
  std::vector<ExtendedHistogram_3D> const& hTemplates
){
  foutput->cd();
  vector<TH3F*> hList;
  for (auto& htpl:hTemplates){
    TString tplname=htpl.getName();
    TH3F* htmp = new TH3F(*(htpl.getHistogram())); htmp->SetName(tplname+"_tmp");
    multiplyBinWidth(htmp);
    hList.push_back(htmp);
  }
  thePerProcessHandle->recombineHistogramsToTemplates(hList);
  for(unsigned int t=0;t<hTemplates.size();t++){
    auto*& htpl=hList.at(t);
    TString tplname=hTemplates.at(t).getName();
    TString tpltitle=hTemplates.at(t).getTitle();
    TString tplxtitle=hTemplates.at(t).getHistogram()->GetXaxis()->GetTitle();
    TString tplytitle=hTemplates.at(t).getHistogram()->GetYaxis()->GetTitle();
    TString tplztitle=hTemplates.at(t).getHistogram()->GetZaxis()->GetTitle();
    TDirectory* savedir=foutput->mkdir(tplname+"_control_" + tplxtitle + "_slices");
    savedir->cd();
    for (int ix=1; ix<=htpl->GetNbinsX(); ix++){
      TString slicename=tplname+"_"+Form("%s_Slice%i", KDbinning.at(0).getLabel().Data(), ix);
      TString slicetitle = tplxtitle + Form(": [%.1f, %.1f]", KDbinning.at(0).getBinLowEdge(ix-1), KDbinning.at(0).getBinLowEdge(ix));
      TH2F* hSlice = getHistogramSlice(htpl, 1, 2, ix, ix, slicename);
      double integral = getHistogramIntegralAndError(hSlice, 1, hSlice->GetNbinsX(), 1, hSlice->GetNbinsY(), false, nullptr);
      if (integral!=0.) hSlice->Scale(1./integral);
      divideBinWidth(hSlice);
      hSlice->SetOption("colz");
      hSlice->SetTitle(slicetitle);
      hSlice->GetXaxis()->SetTitle(tplytitle);
      hSlice->GetYaxis()->SetTitle(tplztitle);
      savedir->WriteTObject(hSlice);
      delete hSlice;
    }
    foutput->cd();
  }
  for (auto& h:hList) delete h;
  rootdir->cd();
}


#endif
