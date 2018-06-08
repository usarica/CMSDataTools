#include "common_includes.h"
#include "TemplatesEventAnalyzer.h"


// Process handle
typedef QQBkgProcessHandler ProcessHandleType;
const ProcessHandleType& theProcess = TemplateHelpers::OffshellQQBkgProcessHandle;

// Constants to affect the template code
#ifndef outputdir_def
#define outputdir_def
const TString user_output_dir = "output/";
#endif

void testCategorization(const Category category){
  const Channel channel = k4mu;
  const SystematicVariationTypes syst = sNominal;
  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);

  // Get list of samples
  vector<TString> strSamples;
  vector<TString> strSampleIdentifiers;
  strSampleIdentifiers.push_back("qq_Bkg_Combined");
  getSamplesList(theSqrts, strSampleIdentifiers, strSamples);

  // Register the discriminants
  vector<KDspecs> KDlist;
  vector<TString> strExtraCatVars_short;
  if (category!=Inclusive){
    getCategorizationDiscriminants(syst, KDlist);
    getExtraCategorizationVariables<short>(globalCategorizationScheme, syst, strExtraCatVars_short);
  }

  // Get the CJLST set
  CJLSTSet* theSampleSet = new CJLSTSet(strSamples);
  // Book common variables
  theSampleSet->bookXS(); // "xsec"
  theSampleSet->bookOverallEventWgt(); // Gen weights "PUWeight", "genHEPMCweight" and reco weights "dataMCWeight", "trigEffWeight"
  for (auto& tree:theSampleSet->getCJLSTTreeList()){
    // Book common variables needed for analysis
    tree->bookBranch<float>("ZZMass", -1);
    tree->bookBranch<short>("Z1Flav", 0);
    tree->bookBranch<short>("Z2Flav", 0);
    // Variables for KDs
    for (auto& KD:KDlist){ for (auto& v:KD.KDvars) tree->bookBranch<float>(v, 0); }
    // Extra categorization variables
    for (auto& s:strExtraCatVars_short) tree->bookBranch<short>(s, -1);
    tree->silenceUnused(); // Will no longer book another branch
  }

  // Construct reweighting variables vector
  for (int t=ProcessHandleType::QQBkg; t<(int) ProcessHandleType::nQQBkgTypes; t++){
    // Build the analyzer and loop over the events
    TemplatesEventAnalyzer theAnalyzer(theSampleSet, channel, category);
    theAnalyzer.setMaximumEvents(10);
    // Book common variables needed for analysis
    theAnalyzer.addConsumed<float>("PUWeight");
    theAnalyzer.addConsumed<float>("genHEPMCweight");
    theAnalyzer.addConsumed<float>("dataMCWeight");
    theAnalyzer.addConsumed<float>("trigEffWeight");
    theAnalyzer.addConsumed<float>("ZZMass");
    theAnalyzer.addConsumed<short>("Z1Flav");
    theAnalyzer.addConsumed<short>("Z2Flav");
    // Add discriminant builders
    for (auto& KD:KDlist){ theAnalyzer.addDiscriminantBuilder(KD.KDname, KD.KD, KD.KDvars); }
    // Add extra categorization variables
    for (auto& s:strExtraCatVars_short) theAnalyzer.addConsumed<short>(s);
    theAnalyzer.setRecordCategorizationKDs(true);
    theAnalyzer.setRecordKDVariables(true);
    // Loop
    theAnalyzer.loop(true, false, true);

    unsigned int ev=0;
    MELAout << "There are " << theAnalyzer.getProducts().size() << " products" << endl;
    for (auto const& theEvent : theAnalyzer.getProducts()){
      cout << "Event " << ev << " ";
      float ZZMassval;
      theEvent.getNamedVal("ZZMass", ZZMassval);
      cout << "ZZMass: " << ZZMassval << endl;
      for (auto& KD:KDlist){
        vector<float> KDvarvals;
        float KDval;
        theEvent.getNamedVal(KD.KDname, KDval);
        cout << KD.KDname << " = " << KDval << " <<\n";
        for (auto& v:KD.KDvars){
          float KDvarval;
          theEvent.getNamedVal(v, KDvarval);
          cout << "\t- " << v << " = " << KDvarval << '\n';
          KDvarvals.push_back(KDvarval);
        }
        cout << "\t- c = " << KD.KD->getCval(ZZMassval);
        cout << endl;
      }
      ev++;
      cout << endl;
    }
  }

  for (auto& KD:KDlist) delete KD.KD;
  delete theSampleSet;
}

