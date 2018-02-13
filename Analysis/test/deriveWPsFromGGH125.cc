#include "common_includes.h"
#include "TemplatesEventAnalyzer.h"


// Constants to affect the template code
#ifndef outputdir_def
#define outputdir_def
const TString user_output_dir = "output/";
#endif

void deriveWPsFromGGH125(){
  // Fixed efficiency points to match
  const float effjjZH = 1.-1.41275535772775410e-1;
  const float effjjWH = 1.-1.47078390611889760e-1;
  const float effjjVBF = 1.-2.60493000952400970e-1;
  const float effjVBF = 1. - 5.770594e-1;
  const float effToMatch[4]={ effjjVBF, effjjZH, effjjWH, effjVBF };

  const Channel channel = NChannels;
  const Category category = Inclusive;
  const SystematicVariationTypes syst = sNominal;
  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);

  const TString strKD="DjjVBF,DjjZH,DjjWH,DjVBF";
  vector<TString> KDnames;
  splitOptionRecursive(strKD, KDnames, ',');
  MELAout << "Processing KDs " << KDnames << endl;

  // Register the discriminants
  vector<KDspecs> KDlist;
  getCategorizationDiscriminants(syst, KDlist);
  {
    KDspecs KDjVBF("DjVBF");
    KDjVBF.KD = constructKDFromType(kDjVBF, "../data/SmoothKDConstant_m4l_DjVBF_13TeV.root", "sp_gr_varReco_Constant_Smooth");
    KDjVBF.KDvars = getKDVars(kDjVBF);
    KDlist.push_back(KDjVBF);
  }

  CJLSTSet* theSampleSet = new CJLSTSet("ggH125");
  // Book common variables
  theSampleSet->bookXS(); // "xsec"
  theSampleSet->bookOverallEventWgt(); // Gen weights "PUWeight", "genHEPMCweight" and reco weights "dataMCWeight", "trigEffWeight"
  for (auto& tree:theSampleSet->getCJLSTTreeList()){
    // Book common variables needed for analysis
    tree->bookBranch<float>("GenHMass", 0);
    tree->bookBranch<float>("ZZMass", -1);
    tree->bookBranch<short>("Z1Flav", 0);
    tree->bookBranch<short>("Z2Flav", 0);
    // Variables for KDs
    for (auto& KD:KDlist){ for (auto& v:KD.KDvars) tree->bookBranch<float>(v, 0); }
    tree->silenceUnused(); // Will no longer book another branch
  }
  vector<SimpleEntry> events;
  TemplatesEventAnalyzer theAnalyzer(theSampleSet, channel, category);
  theAnalyzer.addMassWindow(std::pair<float, float>(105, 140));
  theAnalyzer.setRecordCategorizationKDs(true);
  // Book common variables needed for analysis
  theAnalyzer.addConsumed<float>("PUWeight");
  theAnalyzer.addConsumed<float>("genHEPMCweight");
  theAnalyzer.addConsumed<float>("dataMCWeight");
  theAnalyzer.addConsumed<float>("trigEffWeight");
  theAnalyzer.addConsumed<float>("GenHMass");
  theAnalyzer.addConsumed<float>("ZZMass");
  theAnalyzer.addConsumed<short>("Z1Flav");
  theAnalyzer.addConsumed<short>("Z2Flav");
  // Add discriminant builders
  for (auto& KD:KDlist){ theAnalyzer.addDiscriminantBuilder(KD.KDname, KD.KD, KD.KDvars); }
  // Loop
  theAnalyzer.setExternalProductList(&events);
  theAnalyzer.loop(true, false, true);

  delete theSampleSet;
  for (auto& KD:KDlist) delete KD.KD;

  MELAout << "There are " << events.size() << " products" << endl;

  TFile* foutput = TFile::Open("tmp.root", "recreate");
  vector<TH1F*> hKD;
  for (auto const& KDname:KDnames){
    TH1F* htmp = new TH1F(KDname, "", 100000, 0, 1);
    hKD.push_back(htmp);
  }

  for (SimpleEntry& ev:events){
    float weight;
    ev.getNamedVal("weight", weight);
    for (unsigned int ikd=0; ikd<KDnames.size(); ikd++){
      float KDval;
      ev.getNamedVal(KDnames.at(ikd), KDval);
      if (KDval<0.) continue;
      else if (KDval>=1.) KDval=1.-1e-6;
      if (KDnames.at(ikd)=="DjVBF"){
        float KDval2;
        ev.getNamedVal("DjjVBF", KDval2);
        if (KDval2>=0.) continue;
      }
      hKD[ikd]->Fill(KDval, weight);
    }
  }

  for (unsigned int ikd=0; ikd<KDnames.size(); ikd++){
    MELAout << "Total count with " << KDnames.at(ikd) << ">=0: " << hKD[ikd]->Integral(1, hKD[ikd]->GetNbinsX()) << endl;
    hKD[ikd]->Scale(1./hKD[ikd]->Integral(1, hKD[ikd]->GetNbinsX()));
    float eff=0;
    float diffeff=-2;
    for (int bin=1; bin<=hKD[ikd]->GetNbinsX(); bin++){
      if (diffeff>=0.){
        MELAout << "New WP for " << KDnames.at(ikd) << " = " << hKD[ikd]->GetXaxis()->GetBinLowEdge(bin) << endl;
        break;
      }
      else{
        eff += hKD[ikd]->GetBinContent(bin);
        diffeff = eff-effToMatch[ikd];
      }
    }
    foutput->WriteTObject(hKD[ikd]);
    delete hKD[ikd];
  }

  foutput->Close();
}
