#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "TChain.h"
#include "TString.h"
#include "TSpline.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "./data/ZZ4l_125_6_Samples.h"
#include "./data/FitDimensionsList.h"
#include "./data/HZZ4l_LeptonInterference.h"
#include "./data/Uncertainty_Tables.h"

using namespace std;

//Initializers
bool enableDebugging = false;
int useAnomalousCouplings=kAddfLQ;
bool useDjettagging=true;
enum histtypes{kSigHist,kBkgHist,kIntHist,kBSIHist,kBSI10Hist,kBSI25Hist};
void makeCombineTemplateswithTrees_Modified_MCFM_GenLevelVBF_one(int folder, int erg_tev, int tFitD, int Systematics, bool isSmooth, int Djettag);
double doLinearCombination(double first, double c_first, double second, double c_second, double bkg, int outputtype);
double doPeakCombination(double first, double c_first, double second, double c_second, double bkg, int outputtype);
void oneDlinearcombination(TH1F* first, int firsttype, TH1F* second, int secondtype, TH1F* input, int inputtype, TH1F* finaloutput, int outputtype, TH1F* finaloutput2=0, int output2type=-99);
void twoDlinearcombination(TH2F* first, int firsttype, TH2F* second, int secondtype, TH2F* input, int inputtype, TH2F* finaloutput, int outputtype, TH2F* finaloutput2=0, int output2type=-99);
void threeDlinearcombination(TH3F* first, int firsttype, TH3F* second, int secondtype, TH3F* input, int inputtype, TH3F* finaloutput, int outputtype, TH3F* finaloutput2=0, int output2type=-99);
void wipeOverUnderFlows(TH1F* hwipe);
void wipeOverUnderFlows(TH2F* hwipe);
void wipeOverUnderFlows(TH3F* hwipe);
void progressbar(int val, int tot);

//Make VBF uncertainty
TGraph* make_VBF_UncertaintyGraph(int EnergyIndex, int tUnc){ // tUnc==1,2: QCD,PDF Up, ==-1-2: QCD,PDF Dn, ==3,-3: Up,Dn in quadrature
  double x[NMASSES_VBF];
  double y[NMASSES_VBF];

  TString tgname = "tg_";
  if (abs(tUnc)==1) tgname.Append("QCD");
  else if (abs(tUnc)==2) tgname.Append("PDF");
  else if (abs(tUnc)==3) tgname.Append("All");
  else return 0;
  if (tUnc>0) tgname.Append("Up");
  else tgname.Append("Dn");

  for (int a=0; a<NMASSES_VBF; a++){
    x[a] = VBF_QCD_PDF[EnergyIndex][a][0];

    if (tUnc==1){
      y[a] = VBF_QCD_PDF[EnergyIndex][a][1];
    }
    else if (tUnc==-1){
      y[a] = VBF_QCD_PDF[EnergyIndex][a][2];
    }
    else if (tUnc==2){
      y[a] = VBF_QCD_PDF[EnergyIndex][a][3];
    }
    else if (tUnc==-2){
      y[a] = VBF_QCD_PDF[EnergyIndex][a][4];
    }
    else if (tUnc==3){
      y[a] = sqrt(pow(VBF_QCD_PDF[EnergyIndex][a][1], 2) + pow(VBF_QCD_PDF[EnergyIndex][a][3], 2));
    }
    else if (tUnc==-3){
      y[a] = -sqrt(pow(VBF_QCD_PDF[EnergyIndex][a][2], 2) + pow(VBF_QCD_PDF[EnergyIndex][a][4], 2));
    }
    y[a] = y[a]*0.01 + 1.;
  }

  TGraph* tg = new TGraph(NMASSES_VBF, x, y);
  tg->SetNameTitle(tgname,tgname);
  return tg;
};


//Make lepton interference graph
TGraph* make_HZZ_LeptonInterferenceGraph(){
	float x[leptonInterf_YR3_Size];
	float y[leptonInterf_YR3_Size];
	for(int a=0;a<leptonInterf_YR3_Size;a++){
		x[a] = leptonInterf_YR3[a][0];
		y[a] = leptonInterf_YR3[a][1];
	}
	TGraph* tg = new TGraph(leptonInterf_YR3_Size,x,y);
	tg->SetName("tgHZZ_LeptonInterference");
	tg->SetTitle("H#rightarrowZZ#rightarrow4l Lepton Interference Weight on 4e, 4#mu wrt. 2e2#mu");

	return tg;
};

//Main Function, runs over all desired iterations
void makeCombineTemplateswithTrees_Modified_MCFM_GenLevelVBF(int djet_type=0){
	bool isSmooth=false;
	const int kNumSyst=5;
	int systematics[kNumSyst]={0,1,-1,2,-2};
	for(int i=0;i<kNumSyst;++i){
		for(int usesmooth=0;usesmooth<2;++usesmooth){
			if(usesmooth==0) isSmooth=false;
			if(usesmooth==1) isSmooth=true;
			for(int CoM=7;CoM<9;++CoM){
				for(int channel=0;channel<3;++channel){
					//if(useDjettagging){for(int Djettag=-1;Djettag<2;++Djettag) makeCombineTemplateswithTrees_Modified_MCFM_GenLevelVBF_one(channel,CoM,6,systematics[i],isSmooth,Djettag);}
					//else makeCombineTemplateswithTrees_Modified_MCFM_GenLevelVBF_one(channel,CoM,6,systematics[i],isSmooth,0);	
					makeCombineTemplateswithTrees_Modified_MCFM_GenLevelVBF_one(channel,CoM,6,systematics[i],isSmooth,djet_type);	
				}	
			}
		}
	}
}

//Function to build one template
// folder = 0,1,2 (final state corresponds to 4mu,4e,2mu2e respectively)
// erg_tev = 7,8 (CoM energy)
// tFitD = [0,16] (choice of Discriminant, see FitDimensionsList.h for list; only tFitd works right now)
// Systematics = [-2,2] (Flag for systematics. 0=Nominal, +/-1=QCD, +/-2=PDF)
// isSmooth = true/false (flag to apply smoothing at this stage, both needed for later stage)
// Djettag = -1,0,1 (events fail Djet > 0.5 tag, events are not tagged, events pass Djet > 0.5 tag respectively)
void makeCombineTemplateswithTrees_Modified_MCFM_GenLevelVBF_one(int folder, int erg_tev, int tFitD, int Systematics, bool isSmooth, int Djettag){
	char TREE_NAME[] = "SelectedTree";
	TString INPUT_NAME = "HZZ4lTree_ggTo";
	TString OUTPUT_NAME = "HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_";
	if (useAnomalousCouplings > 0) OUTPUT_NAME += strAnomalousType[useAnomalousCouplings] + "_";
	if (!isSmooth) OUTPUT_NAME += "Raw_";
	OUTPUT_NAME += "_GenLevelVBF_" + TString(strFitDim[tFitD]) + "_";
	TString comstring;
	comstring.Form("%i", erg_tev);
	TString erg_dir;
	erg_dir.Form("LHC_%iTeV/", erg_tev);
	if (Djettag == 0){
		if (Systematics == 0) OUTPUT_NAME += "Nominal";
		if (Systematics == 1) OUTPUT_NAME += "SysUp_ggQCD";
		if (Systematics == -1) OUTPUT_NAME += "SysDown_ggQCD";
		if (Systematics == 2) OUTPUT_NAME += "SysUp_ggPDF";
		if (Systematics == -2) OUTPUT_NAME += "SysDown_ggPDF";
	}
	if (Djettag == -1){
		if (Systematics == 0) OUTPUT_NAME += "Nominal_nonDjet";
		if (Systematics == 1) OUTPUT_NAME += "SysUp_ggQCD_nonDjet";
		if (Systematics == -1) OUTPUT_NAME += "SysDown_ggQCD_nonDjet";
		if (Systematics == 2) OUTPUT_NAME += "SysUp_ggPDF_nonDjet";
		if (Systematics == -2) OUTPUT_NAME += "SysDown_ggPDF_nonDjet";
	}
	if (Djettag == 1){
		if (Systematics == 0) OUTPUT_NAME += "Nominal_Djet";
		if (Systematics == 1) OUTPUT_NAME += "SysUp_ggQCD_Djet";
		if (Systematics == -1) OUTPUT_NAME += "SysDown_ggQCD_Djet";
		if (Systematics == 2) OUTPUT_NAME += "SysUp_ggPDF_Djet";
		if (Systematics == -2) OUTPUT_NAME += "SysDown_ggPDF_Djet";
	}
  TString OUTPUT_LOG_NAME = OUTPUT_NAME;
  OUTPUT_NAME += ".root";
  OUTPUT_LOG_NAME += ".log";
  TString systname[5] ={ "SysDown_ggPDF", "SysDown_ggQCD", "Nominal", "SysUp_ggQCD", "SysUp_ggPDF" };
	TString djetname[3] = { "Djet < 0.5", "Nominal", "Djet>=0.5" };

	TString INPUT_VBFREF_NAME = "HtoZZ4l_MCFM_125p6_ModifiedTemplatesForCombine_";
	if (useAnomalousCouplings > 0) INPUT_VBFREF_NAME += strAnomalousType[useAnomalousCouplings] + "_";
	if (!isSmooth) INPUT_VBFREF_NAME += "Raw_";
	INPUT_VBFREF_NAME += TString(strFitDim[tFitD]) + "_";
	if (Djettag == 0) INPUT_VBFREF_NAME += "Nominal.root";
	if (Djettag == -1) INPUT_VBFREF_NAME += "Nominal_nonDjet.root";
	if (Djettag == 1) INPUT_VBFREF_NAME += "Nominal_Djet.root";

	TString sample_VBF_suffix[4] = {
    "BSI_VBF_Phantom",
    "Bkg_VBF_Phantom",
		"BSI10_VBF_Phantom",
		"BSI25_VBF_Phantom"
	};

	TString cinput_common = user_gg2VV_location + erg_dir + "/" + user_folder[folder] + "/";
	TString cinput_common_recover;
	if (folder != 2) cinput_common_recover = user_gg2VV_location + erg_dir + "/" + user_folder[2] + "/";
	TString cinput_common_qqZZ = user_gg2VV_location;
	TString cinput_common_ZX = user_gg2VV_location;
	TString cinput_common_VBF = user_gg2VV_location; //CHANGE - made to point to correct directory

	int EnergyIndex = 1;
	if (erg_tev == 7) EnergyIndex = 0;
	float lowside[3] = { 220, 230, 240 };
	double mPOLE = 125.6;
	float ZZMass_PeakCut[2] = { 105.6, 140.6 }; // Spin 0 analysis

	float templateWeight = 1;
  float MC_weight;
  float MC_weight_down;
	float MC_weight_up;
	float MC_weight_Kfactor = 1;
	float MC_weight_ggZZLepInt = 1;
	float GenHMass = 0;
  float GenDiJetMass;
  float MC_weight_Kfactor_Norm_down = 1;
	float MC_weight_Kfactor_Norm_up = 1;
	float MC_weight_PDF_Norm_down = 1;
	float MC_weight_PDF_Norm_up = 1;
	float ZZMass = 0;

  TGraph* tg_VBF_unc=0;
  if (Systematics>0) tg_VBF_unc = make_VBF_UncertaintyGraph(EnergyIndex, 3); // Retrieve sum in quadrature!
  else if (Systematics<0) tg_VBF_unc = make_VBF_UncertaintyGraph(EnergyIndex, -3); // Retrieve sum in quadrature!

	TString cinput_KDFactor = "./data/HZZ4l-KDFactorGraph";
	if (EnergyIndex == 0) cinput_KDFactor = cinput_KDFactor + "_7TeV";
	cinput_KDFactor = cinput_KDFactor + ".root";
	TFile* finput_KDFactor = new TFile(cinput_KDFactor, "read");
	TString tgkfname = "KDScale_";
	tgkfname = tgkfname + "AllFlavors_UnNormalized";
	TGraphAsymmErrors* tgkf = (TGraphAsymmErrors*)finput_KDFactor->Get(tgkfname);

	TString Djetcutfilename = "./data/HZZ4l-DjetCutShapes";
	if (EnergyIndex == 0) Djetcutfilename += "_7TeV";
	Djetcutfilename += ".root";
	TFile* Djetcutfile = new TFile(Djetcutfilename, "read");

  TString INPUT_VBFVHSCALE_NAME = "HtoZZ4l_Phantom_125p6_VHDistributions.root";
  TString cinput_vbfvhscale = user_TemplateswithTrees_dir + "../VHContributions/" + erg_dir;
  cinput_vbfvhscale += user_folder[folder] + "/";
  cinput_vbfvhscale += INPUT_VBFVHSCALE_NAME;
  TFile* finput_vbfvhscale = new TFile(cinput_vbfvhscale, "read");
  TH1F* h_VBFVH_Scale = (TH1F*)finput_vbfvhscale->Get("Phantom_VBFVH_ScalingRatio");

	//double overall_VBF_scale=1;
	double nVBFPeak[4] = { 0 };
//	for (int e = 0; e < 2; e++){ for (int ss = 0; ss < 3; ss++){ nSM_ScaledPeak[e][ss] /= luminosity[e]; VBF_Sig_Datacard[e][ss] /= luminosity[e]; } }

	TGraph* tg_interf = make_HZZ_LeptonInterferenceGraph();

	TString cinput_VBF_Sig = "./data/HZZ4l-125_6-" + comstring + "TeV-Sig_MCFM_PhantomVBF_Comparison.root";
	TFile* finput_VBF = new TFile(cinput_VBF_Sig, "read");
	TSpline3* tsp_VBF_Sig = (TSpline3*)finput_VBF->Get("Spline3");

	TH2F* hStore_ZX_Unconditional;
	TH2F* hStore_qqZZ_Unconditional;

	for (int lo = 0; lo < 1; lo++){
		TString coutput_common = user_TemplateswithTrees_dir + erg_dir;
		coutput_common += user_folder[folder] + "/";
		gSystem->Exec("mkdir -p " + coutput_common);

		TString coutput = coutput_common + OUTPUT_NAME;
    TString coutput_log = coutput_common + OUTPUT_LOG_NAME;
    ofstream tout(coutput_log.Data(), ios::out);
    cout << "Opened file " << coutput_log << endl;
    TFile* foutput = new TFile(coutput, "recreate");

    cout<<endl;
    cout<<"==============================="<<endl;
    cout<<"CoM Energy: "<<erg_tev<<" TeV"<<endl;
    cout<<"Decay Channel: "<<user_folder[folder]<<endl;
    cout<<"Systematic: "<<systname[Systematics+2]<<endl;
    cout<<"Djet cut: "<<djetname[Djettag+1]<<endl;
    if (isSmooth) cout<<"Smooth version"<<endl;
    cout<<"==============================="<<endl;
    cout<<endl;
    tout<<"==============================="<<endl;
    tout<<"CoM Energy: "<<erg_tev<<" TeV"<<endl;
    tout<<"Decay Channel: "<<user_folder[folder]<<endl;
    tout<<"Systematic: "<<systname[Systematics+2]<<endl;
    tout<<"Djet cut: "<<djetname[Djettag+1]<<endl;
    if (isSmooth) tout<<"Smooth version"<<endl;
    tout<<"==============================="<<endl;
    tout<<endl;

		//Grab VBF Sig templates made in makeCombineTemplates_Modified_MCFM.c
		TString cinput_VBFRef = user_dir + erg_dir + user_folder[folder] + "/" + INPUT_VBFREF_NAME;
		TFile* finput_VBFRef = new TFile(cinput_VBFRef, "read");
		cout << cinput_VBFRef << endl;
		double nVBF_Sig_Reweighted = 0;
		TH1F* hVBFRef_1D = (TH1F*)finput_VBFRef->Get("T_1D_VBF_1");
		TH2F* hVBFRef = (TH2F*)finput_VBFRef->Get("T_2D_VBF_1");
		nVBF_Sig_Reweighted = hVBFRef->Integral(0, hVBFRef->GetNbinsX() + 1, 0, hVBFRef->GetNbinsY() + 1, "width")*luminosity[EnergyIndex];
		foutput->cd();

		float ZZMass_cut[2] = { lowside[lo], 1600 };
		float ZZwidth = 20.0;
		const int nbinsx = (ZZMass_cut[1] - ZZMass_cut[0]) / ZZwidth;
		float kDXarray[nbinsx + 1];
		for (int bin = 0; bin < nbinsx + 1; bin++){
			kDXarray[bin] = ZZMass_cut[0] + ZZwidth*bin;
		}

		int nbinsy = 30;
		float kDYarray[nbinsy + 1];
		float kDY_bounds[2] = { 0, 1 };
		if (tFitD == 3){ kDY_bounds[0] = -7.0; kDY_bounds[1] = 3.0; }
		if (tFitD == 4 || tFitD == 5 || tFitD == 9 || tFitD == 10 || tFitD == 11 || tFitD == 15){ kDY_bounds[0] = -1.0; kDY_bounds[1] = 1.0; }
		for (int bin = 0; bin < nbinsy + 1; bin++){
			double binwidth = (kDY_bounds[1] - kDY_bounds[0]) / nbinsy;
			kDYarray[bin] = kDY_bounds[0] + binwidth*bin;
		}

		const int kNumTemplates = 8;
		TH1F** D_temp_1D[kNumTemplates];
		TH2F** D_temp_2D[kNumTemplates];
		double overall_scale[kNumTemplates] = { 1 };
		TString templatenames[kNumTemplates] = { "ggF Sig", "gg Bkg", "ggF Int", "qqZZ", "Z+X", "VBF Sig", "VBF Bkg", "VBF Int" };
		// Build template structure, including anomalous couplings
		for (int t = 0; t < kNumTemplates; t++){
			if (t <= 2){ // gg(H)VV
				D_temp_1D[t] = new TH1F*[(nAnomalousCouplingTemplates[useAnomalousCouplings][0])];
				D_temp_2D[t] = new TH2F*[(nAnomalousCouplingTemplates[useAnomalousCouplings][0])];
			}
			else if (t >= 5 && t <= 7){ // VV(H)VV
				D_temp_1D[t] = new TH1F*[(nAnomalousCouplingTemplates[useAnomalousCouplings][1])];
				D_temp_2D[t] = new TH2F*[(nAnomalousCouplingTemplates[useAnomalousCouplings][1])];
			}
			else{ // Anything else
				D_temp_1D[t] = new TH1F*[1];
				D_temp_2D[t] = new TH2F*[1];
			}
		}

		int nEntries;
		double nSig_Simulated = 0;
		double nVBF_Sig_Simulated = 0;
		float fitYval = 0;
		MC_weight = 1;
		MC_weight_down = 1;
		MC_weight_up = 1;
		MC_weight_Kfactor = 1;
		MC_weight_Kfactor_Norm_down = 1;
		MC_weight_Kfactor_Norm_up = 1;
		MC_weight_PDF_Norm_down = 1;
		MC_weight_PDF_Norm_up = 1;
		MC_weight_ggZZLepInt = 1;
		GenHMass = 0;
    GenDiJetMass = 0;
		ZZMass = 0;
		float Djet = 0.;

		double nMZZ220[3] = { 0 };

		//Initialize and grab each of the four Phantom trees (change to FullSim)
		TChain* tree_VBF[4];
		TH1F** h1DVBF[4];
		TH2F** h2DVBF[4];
		// Build template structure, including anomalous couplings
		for (int t = 0; t < 4; t++){ // Notice this is purely VBF
			h1DVBF[t] = new TH1F*[(nAnomalousCouplingTemplates[useAnomalousCouplings][1])];
			h2DVBF[t] = new TH2F*[(nAnomalousCouplingTemplates[useAnomalousCouplings][1])];
		}

		TH1F* h1DVBFSigRatio;
		TH2F* h2DVBFSigRatio;
		if (tFitD == 0){
			h1DVBFSigRatio = new TH1F("h1DSigRatio", "h1DSigRatio", nbinsx, kDXarray);
		}
		else{
			h1DVBFSigRatio = new TH1F("h1DSigRatio", "h1DSigRatio", nbinsy, kDYarray);
			h2DVBFSigRatio = new TH2F("h2DSigRatio", "h2DSigRatio", nbinsx, kDXarray, nbinsy, kDYarray);
		}
		for (int tr = 0; tr < 4; tr++){
			tree_VBF[tr] = new TChain(TREE_NAME);
		}

		for (int decay = 0; decay < 3; decay++){
			TString channelname = user_folder[decay];
			if (decay == 2) channelname = "2e2mu";

			/*
			0: BSI
			1: BKG
			2: BSI10
			3: BSI25
			*/
			for (int type_vbf = 0; type_vbf < 4; type_vbf++){
				TString cinput_VBF_Bkg = cinput_common + "HZZ4lTree_ZZTo" + channelname + "JJ_" + sample_suffix_Phantom[type_vbf] + "_Reprocessed.root";
				tree_VBF[type_vbf]->Add(cinput_VBF_Bkg);
			}
		}

		//Fill Template(s) from each Phantom sample with:
		// gen mZ1/mZ2>120
		// Reweighting for lepton interference (necessary)
		// Possible mZZ cuts?
		//Then write smoothed templates to file 
    double nRead_VBF_SM[4]={ 0 };
		for (int tr = 0; tr < 4; tr++){
			int nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][1];

      tree_VBF[tr]->SetBranchAddress("GenDiJetMass", &GenDiJetMass);
      tree_VBF[tr]->SetBranchAddress("GenHMass", &GenHMass);
      tree_VBF[tr]->SetBranchAddress("ZZMass", &ZZMass);
			tree_VBF[tr]->SetBranchAddress("MC_weight", &MC_weight);
			tree_VBF[tr]->SetBranchAddress("Djet_VAJHU", &Djet);
			if (tFitD != 0) tree_VBF[tr]->SetBranchAddress(strFitDim[tFitD], &fitYval);

			TString templatename_1D_core;
			TString templatename_2D_core;
			TString templatename_1D;
			TString templatename_2D;

			templatename_1D_core = "htemp1DVBF_";
			templatename_1D_core += sample_VBF_suffix[tr];
			templatename_2D_core = "htemp2DVBF_";
			templatename_2D_core += sample_VBF_suffix[tr];

			for (int al = 0; al<nAnomalousLoops; al++){
				templatename_1D = templatename_1D_core;
				templatename_2D = templatename_2D_core;
				if (useAnomalousCouplings == kAddfLQ && al>0){
					templatename_1D.Append(Form("_mZZ2_%i", al));
					templatename_2D.Append(Form("_mZZ2_%i", al));
				}
				if (tFitD == 0){
					h1DVBF[tr][al] = new TH1F(templatename_1D, templatename_1D, nbinsx, kDXarray);
				}
				else{
					h1DVBF[tr][al] = new TH1F(templatename_1D, templatename_1D, nbinsy, kDYarray);
					h2DVBF[tr][al] = new TH2F(templatename_2D, templatename_2D, nbinsx, kDXarray, nbinsy, kDYarray);
				}
			}

			int nVBFEntries = tree_VBF[tr]->GetEntries();
			for (int ev = 0; ev < nVBFEntries; ev++){
//				progressbar(ev, tree_VBF[tr]->GetEntries());
				tree_VBF[tr]->GetEntry(ev);
				if (fitYval != fitYval) continue;
        // Protect against any KD exceeding boundaries
        if (tFitD!=0 && fitYval>=kDYarray[nbinsy]) fitYval = kDYarray[nbinsy] - (kDYarray[nbinsy]-kDYarray[nbinsy-1])*0.1*(ev+1.)/(nVBFEntries+1.);
        if (tFitD!=0 && fitYval<kDYarray[0]) fitYval = kDYarray[0] + (kDYarray[1]-kDYarray[0])*0.1*(ev+1.)/(nVBFEntries+1.);
        if (tFitD!=0 && (fitYval>=kDYarray[nbinsy] || fitYval<kDYarray[0])) cout << "Fix has been numerically unsuccessful for " << tree_VBF[tr]->GetName() << endl;

				double weight = MC_weight;
        int VBF_VH_rewgt_bin = h_VBFVH_Scale->GetXaxis()->FindBin(GenDiJetMass);
        if (VBF_VH_rewgt_bin>h_VBFVH_Scale->GetNbinsX()) VBF_VH_rewgt_bin = h_VBFVH_Scale->GetNbinsX();
        weight *= h_VBFVH_Scale->GetBinContent(VBF_VH_rewgt_bin); // VBF-VH scale for Phantom
        if (tg_VBF_unc!=0 && Systematics!=0) weight *= tg_VBF_unc->Eval(GenHMass);

        if (ev == 0 && enableDebugging) cout << "Weight: " << weight << endl;
				if (ZZMass < ZZMass_PeakCut[1] && ZZMass >= ZZMass_PeakCut[0]){
					nVBFPeak[tr] += weight;
				}
				if (ZZMass >= ZZMass_cut[1] || ZZMass < ZZMass_cut[0]) continue;

				if (Djettag == -1 && Djet >= 0.5) continue;
				if (Djettag == 1 && Djet < 0.5) continue;

				// Anomalous couplings loop
				for (int al = 0; al<nAnomalousLoops; al++){
					double fillWeight = weight;
					if (useAnomalousCouplings == kAddfLQ && al>0){
						fillWeight *= pow(GenHMass / mPOLE, 2 * al);
					}
					if (tFitD == 0) h1DVBF[tr][al]->Fill(ZZMass, fillWeight);
					else{
						h1DVBF[tr][al]->Fill(fitYval, fillWeight);
						h2DVBF[tr][al]->Fill(ZZMass, fitYval, fillWeight);
					}
          if (al==0) nRead_VBF_SM[tr] += fillWeight;
				}
			}
			cout << endl;
      cout << "Read " << nRead_VBF_SM[tr]*luminosity[EnergyIndex] << " for VBF tree " << tr << endl;

			ZZMass = 0;
			MC_weight = 1;
			fitYval = 0;

			for (int al = 0; al < nAnomalousLoops; al++){
				double presmoothInt = h1DVBF[tr][al]->Integral();
				if (isSmooth) h1DVBF[tr][al]->Smooth(1, "k3a");
        wipeOverUnderFlows(h1DVBF[tr][al]);
				double postsmoothInt = h1DVBF[tr][al]->Integral();
				h1DVBF[tr][al]->Scale(presmoothInt / postsmoothInt);
        cout << "Scaled VBF " << tr << " 1D nominal template by " << (presmoothInt / postsmoothInt) << endl;

				foutput->WriteTObject(h1DVBF[tr][al]);
				if (tFitD != 0){
					presmoothInt = h2DVBF[tr][al]->Integral();
					if (isSmooth) h2DVBF[tr][al]->Smooth(1, "k3a");
          wipeOverUnderFlows(h2DVBF[tr][al]);
          postsmoothInt = h2DVBF[tr][al]->Integral();
					h2DVBF[tr][al]->Scale(presmoothInt / postsmoothInt);
					foutput->WriteTObject(h2DVBF[tr][al]);
          cout << "Scaled VBF " << tr << " 2D nominal template by " << (presmoothInt / postsmoothInt) << endl;
        }
			}
		}

		//Make VBF Sig/Int from linear combinations of above templates
		//0: VBF Sig
		//2: VBF Int
		//	 For 7 TeV samples, BSI25, Bkg, and BSI10 are used
		//	 For 8 TeV samples, BSI, Bkg, and BSI10 are used
		if (EnergyIndex == 0){
			int nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][1];
			for (int al = 0; al < nAnomalousLoops; al++){
				oneDlinearcombination(h1DVBF[2][al], kBSI10Hist, h1DVBF[1][al], kBkgHist, h1DVBF[3][al], kBSI25Hist, h1DVBF[0][al], kSigHist, h1DVBF[2][al], kIntHist);
				if (tFitD != 0){
					twoDlinearcombination(h2DVBF[2][al], kBSI10Hist, h2DVBF[1][al], kBkgHist, h2DVBF[3][al], kBSI25Hist, h2DVBF[0][al], kSigHist, h2DVBF[2][al], kIntHist);
				}
			}
			double nVBFsigtemp = doPeakCombination(nVBFPeak[2], 10, nVBFPeak[3], 25, nVBFPeak[1], kSigHist);
			double nVBFinterftemp = doPeakCombination(nVBFPeak[2], 10, nVBFPeak[3], 25, nVBFPeak[1], kIntHist);
			for (int tr = 0; tr < 4; tr++) cout << "nVBFPeak[" << tr << "]: " << nVBFPeak[tr] << endl;
			nVBFPeak[0] = nVBFsigtemp;
			nVBFPeak[2] = nVBFinterftemp;
		}
		else if (EnergyIndex == 1){
			int nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][1];
			for (int al = 0; al < nAnomalousLoops; al++){
				oneDlinearcombination(h1DVBF[0][al], kBSIHist, h1DVBF[1][al], kBkgHist, h1DVBF[2][al], kBSI10Hist, h1DVBF[0][al], kSigHist, h1DVBF[2][al], kIntHist);
				if (tFitD != 0){
					twoDlinearcombination(h2DVBF[0][al], kBSIHist, h2DVBF[1][al], kBkgHist, h2DVBF[2][al], kBSI10Hist, h2DVBF[0][al], kSigHist, h2DVBF[2][al], kIntHist);
				}
			}
			double nVBFsigtemp = doPeakCombination(nVBFPeak[2], 10, nVBFPeak[0], 1, nVBFPeak[1], kSigHist);
			double nVBFinterftemp = doPeakCombination(nVBFPeak[2], 10, nVBFPeak[0], 1, nVBFPeak[1], kIntHist);
			for (int tr = 0; tr < 4; tr++) cout << "nVBFPeak[" << tr << "]: " << nVBFPeak[tr] << endl;
			nVBFPeak[0] = nVBFsigtemp;
			nVBFPeak[2] = nVBFinterftemp;
		}

//    double vbfscale = VBF_Sig_Datacard[EnergyIndex][folder] / (nVBFPeak[0]*luminosity[EnergyIndex]);
    double vbfscale = 1;
    for (int tr = 0; tr < 4; tr++){
			for (int al = 0; al < nAnomalousCouplingTemplates[useAnomalousCouplings][1]; al++){
				h1DVBF[tr][al]->Scale(vbfscale);
				if (tFitD != 0) h2DVBF[tr][al]->Scale(vbfscale);
        if (tr == 0 && al == 0){
          cout << "VBF Signal scale for syst=0: " << vbfscale << endl;
          tout << "VBF Signal scale for syst=0: " << vbfscale << endl;
        }
      }
		}
		if (tFitD != 0){
			nVBF_Sig_Simulated = h2DVBF[0][0]->Integral(1, h2DVBF[0][0]->GetNbinsX(), 0, h2DVBF[0][0]->GetNbinsY() + 1)*luminosity[EnergyIndex];
			cout << h2DVBF[0][0]->Integral()*luminosity[EnergyIndex] << endl;
			cout << h2DVBF[1][0]->Integral()*luminosity[EnergyIndex] << endl;
			cout << h2DVBF[2][0]->Integral()*luminosity[EnergyIndex] << endl;
			cout << h2DVBF[3][0]->Integral()*luminosity[EnergyIndex] << endl;
		}

		//Template and tree filler
		// These integrals are important for scaling reweighted VBF Sig, Bkg or Int wrt signal peak
		double intermediateIntegral_VBF_1D[3] = { 0 };
		double intermediateIntegral_VBF_2D[3] = { 0 };

		for (int t = 0; t < kNumTemplates; t++){
			int nAnomalousLoops = 1;
			if (t <= 2) nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][0];
			else if (t >= 5 && t <= 7) nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][1];

			TChain* tree = new TChain(TREE_NAME);
			TTree** templateTree = new TTree*[nAnomalousLoops];

			TF1* Djetcutshape = 0;
			if (Djettag != 0){
				TString cutname;
				//IF USING MCFM FOR DJET CUT SHAPE
				/*if(t==0) cutname="MCFMSig_Djetcutshape";
				if(t==1) cutname="MCFMBkg_Djetcutshape";
				if(t==2) cutname="MCFMBSI_Djetcutshape";*/
				//IF USING MINLO FOR DJET CUT SHAPE
				if (t == 0 || t == 1 || t == 2) cutname = "MINLO_Djetcutshape";
				Djetcutshape = (TF1*)Djetcutfile->Get(cutname);
			}
			//Get file names, templates, and appropriate branches
			if (t == 3){
				for (int b = 0; b < kNumBkg; b++){
					TString cinput = cinput_common_qqZZ;
					cinput += erg_dir + "/" + user_folder[folder] + "/" + hzz4lprefix + sample_BackgroundFile[b] + "_Reprocessed.root";
					tree->Add(cinput);
					cout << cinput << endl;
				}
			}
			else if (t == 4){
				TString cinput = cinput_common_ZX;
				cinput += erg_dir + "/" + user_folder[folder] + "/HZZ4lTree_DoubleOr_CRZLLTree_Reprocessed.root";
				tree->Add(cinput);
				cout << cinput << endl;
			}
			else if (t < 3){
				int tp = t;
				TString cinput_2e2mu = cinput_common + INPUT_NAME + "2e2mu_" + sample_suffix_MCFM[tp] + "_Reprocessed.root";
				TString cinput_4e = cinput_common + INPUT_NAME + "4e_" + sample_suffix_MCFM[tp] + "_Reprocessed.root";
				TString cinput_4mu = cinput_common + INPUT_NAME + "4mu_" + sample_suffix_MCFM[tp] + "_Reprocessed.root";

				tree->Add(cinput_2e2mu);
				cout << cinput_2e2mu << endl;
				tree->Add(cinput_4e);
				cout << cinput_4e << endl;
				tree->Add(cinput_4mu);
				cout << cinput_4mu << endl;
			}

			//Initialize templates
			TString templatename_1D_core;
			TString templatename_2D_core;
			TString templatename_1D;
			TString templatename_2D;
			if (t < 2){
				templatename_1D_core = Form("T_1D_%i", t + 1);
				templatename_2D_core = Form("T_2D_%i", t + 1);
			}
			else if (t == 2){
				templatename_1D_core = Form("T_1D_%i", 4);
				templatename_2D_core = Form("T_2D_%i", 4);
			}
			else if (t == 3){
				templatename_1D_core = Form("T_1D_%s", "qqZZ");
				templatename_2D_core = Form("T_2D_%s", "qqZZ");
			}
			else if (t == 4){
				templatename_1D_core = Form("T_1D_%s", "ZX");
				templatename_2D_core = Form("T_2D_%s", "ZX");
			}
			else if (t == 5 || t == 6){
				templatename_1D_core = Form("T_1D_VBF_%i", t - 4);
				templatename_2D_core = Form("T_2D_VBF_%i", t - 4);
			}
			else if (t == 7){
				templatename_1D_core = Form("T_1D_VBF_%i", 4);
				templatename_2D_core = Form("T_2D_VBF_%i", 4);
			}

			for (int al = 0; al<nAnomalousLoops; al++){
				templatename_1D = templatename_1D_core;
				templatename_2D = templatename_2D_core;
				if (useAnomalousCouplings == kAddfLQ && al>0){
					templatename_1D.Append(Form("_mZZ2_%i", al));
					templatename_2D.Append(Form("_mZZ2_%i", al));
				}
				if (tFitD == 0){
					D_temp_1D[t][al] = new TH1F(templatename_1D, templatename_1D, nbinsx, kDXarray);
					D_temp_1D[t][al]->GetXaxis()->SetTitle(strFitDim_label[tFitD]);
					TString strTreeName = templatename_1D;
					strTreeName += "_Tree";
					templateTree[al] = new TTree(strTreeName, strTreeName);
				}
				else{
					D_temp_1D[t][al] = new TH1F(templatename_1D, templatename_1D, nbinsy, kDYarray);
					D_temp_1D[t][al]->GetXaxis()->SetTitle(strFitDim_label[tFitD]);

					D_temp_2D[t][al] = new TH2F(templatename_2D, templatename_2D, nbinsx, kDXarray, nbinsy, kDYarray);
					D_temp_2D[t][al]->GetXaxis()->SetTitle(strFitDim_label[0]);
					D_temp_2D[t][al]->GetYaxis()->SetTitle(strFitDim_label[tFitD]);

					TString strTreeName = templatename_2D;
					strTreeName += "_Tree";
					templateTree[al] = new TTree(strTreeName, strTreeName);
				}
				templateTree[al]->Branch("templateWeight", &templateWeight);
				templateTree[al]->Branch("ZZMass", &ZZMass);
				if (tFitD != 0) templateTree[al]->Branch(strFitDim[tFitD], &fitYval);
			}

			if (tFitD != 0){
				if (t == 3){
					TString storeName = templatename_2D_core;
					storeName = storeName + "_UnConditional";
					hStore_qqZZ_Unconditional = new TH2F(storeName, storeName, nbinsx, kDXarray, nbinsy, kDYarray);
					hStore_qqZZ_Unconditional->GetXaxis()->SetTitle(strFitDim_label[0]);
					hStore_qqZZ_Unconditional->GetYaxis()->SetTitle(strFitDim_label[tFitD]);
				}
				if (t == 4){
					TString storeName = templatename_2D_core;
					storeName = storeName + "_UnConditional";
					hStore_ZX_Unconditional = new TH2F(storeName, storeName, nbinsx, kDXarray, nbinsy, kDYarray);
					hStore_ZX_Unconditional->GetXaxis()->SetTitle(strFitDim_label[0]);
					hStore_ZX_Unconditional->GetYaxis()->SetTitle(strFitDim_label[tFitD]);
				}
			}

			//Making templates using appropriate weights
			double nTotal = 0;
			if (t < 5){
				if(t!=4) tree->SetBranchAddress("GenHMass", &GenHMass);
				else GenHMass=0;
				tree->SetBranchAddress("ZZMass", &ZZMass);
				if (t != 4) tree->SetBranchAddress("MC_weight", &MC_weight);
				else tree->SetBranchAddress("ZXfake_weightProper", &MC_weight);
				if (t < 3){
					tree->SetBranchAddress("MC_weight_Kfactor", &MC_weight_Kfactor);
					tree->SetBranchAddress("MC_weight_Kfactor_Norm_down", &MC_weight_Kfactor_Norm_down);
					tree->SetBranchAddress("MC_weight_Kfactor_Norm_up", &MC_weight_Kfactor_Norm_up);
					tree->SetBranchAddress("MC_weight_PDF_Norm_down", &MC_weight_PDF_Norm_down);
					tree->SetBranchAddress("MC_weight_PDF_Norm_up", &MC_weight_PDF_Norm_up);
					tree->SetBranchAddress("MC_weight_down", &MC_weight_down);
					tree->SetBranchAddress("MC_weight_up", &MC_weight_up);
				}
				else{
					MC_weight_ggZZLepInt = 1;
					MC_weight_Kfactor = 1;
					MC_weight_down = 1;
					MC_weight_up = 1;
					MC_weight_Kfactor_Norm_down = 1;
					MC_weight_Kfactor_Norm_up = 1;
					MC_weight_PDF_Norm_down = 1;
					MC_weight_PDF_Norm_up = 1;
				}
				if (t<3 || t>4) tree->SetBranchAddress("MC_weight_ggZZLepInt", &MC_weight_ggZZLepInt);
				else MC_weight_ggZZLepInt = 1;
				if (t > 4) tree->SetBranchAddress("Djet_VAJHU", &Djet);
				else Djet = 0.;

				//ONLY CHANGE!!!
				if (tree->GetBranchStatus(strFitDim[tFitD])) tree->SetBranchAddress(strFitDim[tFitD], &fitYval);
				else if (!tree->GetBranchStatus(strFitDim[tFitD])){
					cerr << "Could NOT find branch named " << strFitDim[tFitD] << "!!! Setting strFitDim[" << tFitD << "] = 0." << endl;
					fitYval = 0;
				}

				nEntries = tree->GetEntries();
				for (int ev = 0; ev < nEntries; ev++){
//					progressbar(ev, tree->GetEntries());
					tree->GetEntry(ev);
					if (fitYval != fitYval) continue;
          // Protect against any KD exceeding boundaries
          if (tFitD!=0 && fitYval>=kDYarray[nbinsy]){
            cout << "Found fitYval == " << fitYval;
            fitYval = kDYarray[nbinsy] - (kDYarray[nbinsy]-kDYarray[nbinsy-1])*0.1*(ev+1.)/(nEntries+1.);
            cout << ". Fixed to " << fitYval << endl;
          }
          if (tFitD!=0 && fitYval<kDYarray[0]) fitYval = kDYarray[0] + (kDYarray[1]-kDYarray[0])*0.1*(ev+1.)/(nEntries+1.);
          if (tFitD!=0 && (fitYval>=kDYarray[nbinsy] || fitYval<kDYarray[0])) cout << "Fix has been numerically unsuccessful for " << tree->GetName() << endl;

					double weight = MC_weight;
					if (t < 3) weight *= MC_weight_ggZZLepInt;
					if (abs(Systematics) != 1 && t < 3) weight *= MC_weight_Kfactor;
					if (Systematics == -1 && t < 3) weight *= (1.0 - ggZZ_Syst_AbsNormSyst[EnergyIndex][0])*MC_weight_Kfactor_Norm_down*tgkf->Eval(125.6);
					if (Systematics == 1 && t < 3) weight *= (1.0 + ggZZ_Syst_AbsNormSyst[EnergyIndex][0])*MC_weight_Kfactor_Norm_up*tgkf->Eval(125.6);
					if (Systematics == -2 && t < 3) weight *= (1.0 - ggZZ_Syst_AbsNormSyst[EnergyIndex][1])*MC_weight_PDF_Norm_down;
					if (Systematics == 2 && t<3) weight *= (1.0 + ggZZ_Syst_AbsNormSyst[EnergyIndex][1])*MC_weight_PDF_Norm_up;
					if (t>4) weight *= tsp_VBF_Sig->Eval(GenHMass)*MC_weight_ggZZLepInt; // REMEMBER, THIS IS FROM GGH TEMPLATES,SHOULD NOT BE THE CASE IN 7 TEV VBF PHANTOM SAMPLES

					if (t == 0 && ZZMass >= ZZMass_PeakCut[0] && ZZMass < ZZMass_PeakCut[1]) nSig_Simulated += weight;
					if (t < 3 && ZZMass >= ZZMass_cut[0]) nMZZ220[t] += weight;

					if ((ZZMass >= ZZMass_cut[1] || ZZMass < ZZMass_cut[0])) continue;

					//Code to add Djet reweighting for ggF
					if (t<3){
						if (Djettag == -1){
							weight *= 1 - Djetcutshape->Eval(ZZMass);
						}
						if (Djettag == 1){
							weight *= Djetcutshape->Eval(ZZMass);
						}
					}
					if (t>4){
						if (Djettag == -1 && Djet >= 0.5) continue;
						if (Djettag == 1 && Djet < 0.5) continue;
					}

					nTotal += weight;

					// Anomalous couplings loop
					for (int al = 0; al<nAnomalousLoops; al++){
						double fillWeight = weight;
						if (useAnomalousCouplings == kAddfLQ && al>0){
							fillWeight *= pow(GenHMass / mPOLE, 2 * al);
						}
						if (tFitD == 0) D_temp_1D[t][al]->Fill(ZZMass, fillWeight);
						else D_temp_1D[t][al]->Fill(fitYval, fillWeight);
						if (tFitD > 0) D_temp_2D[t][al]->Fill(ZZMass, fitYval, fillWeight);
					}
				}
				cout << endl;
			}
      if (t!=4){
        cout << templatenames[t] << " Total Simulated: " << nTotal*luminosity[EnergyIndex] << endl;
        tout << templatenames[t] << " Total Simulated: " << nTotal*luminosity[EnergyIndex] << endl;
      }
      else{
        cout << templatenames[t] << " Total Simulated: " << nTotal << endl;
        tout << templatenames[t] << " Total Simulated: " << nTotal << endl;
      }

			//Reweighting normalization
			if (t < 3){
//        cout << nMZZ220[t] * nSM_ScaledPeak[EnergyIndex][folder] / nSig_Simulated*luminosity[EnergyIndex] << endl;
//        cout << nMZZ220[t] * nSM_ScaledPeak[EnergyIndex][folder] / nSig_Simulated << endl;

        double myscale = nSM_ScaledPeak[EnergyIndex][folder] / (nSig_Simulated*luminosity[EnergyIndex]);
				if (Systematics == -1 && t < 3) myscale *= (1.0 - ggZZ_Syst_AbsNormSyst[EnergyIndex][0]);
				if (Systematics == 1 && t < 3) myscale *= (1.0 + ggZZ_Syst_AbsNormSyst[EnergyIndex][0]);
				if (Systematics == -2 && t < 3) myscale *= (1.0 - ggZZ_Syst_AbsNormSyst[EnergyIndex][1]);
				if (Systematics == 2 && t < 3) myscale *= (1.0 + ggZZ_Syst_AbsNormSyst[EnergyIndex][1]);
				overall_scale[t] = myscale;
				for (int al = 0; al<nAnomalousLoops; al++){
					D_temp_1D[t][al]->Scale(myscale);
					if (tFitD>0) D_temp_2D[t][al]->Scale(myscale);
          cout << "Scaling " << templatenames[t] << " (coupling:" << al << ") by " << myscale << endl;
          tout << "Scaling " << templatenames[t] << " (coupling:" << al << ") by " << myscale << endl;
				}
			}
			else if (t >= 5){
				double myscale = (nVBF_Sig_Reweighted / nVBF_Sig_Simulated);

				cout << "Expected signal VBF offshell from reweighting: " << nVBF_Sig_Reweighted << endl;
				cout << "Observed signal VBF offshell from peak normalization: " << nVBF_Sig_Simulated << endl;
				cout << "Re-wgted/full sim. for " << templatenames[t] << ": " << myscale << endl;
        tout << "Expected signal VBF offshell from reweighting: " << nVBF_Sig_Reweighted << endl;
        tout << "Observed signal VBF offshell from peak normalization: " << nVBF_Sig_Simulated << endl;
        tout << "Re-wgted/full sim. for " << templatenames[t] << ": " << myscale << endl;
        overall_scale[t] = vbfscale;

				for (int al = 0; al<nAnomalousLoops; al++){
					for (int binx = 0; binx <= h1DVBF[t - 5][al]->GetNbinsX() + 1; binx++){ // Include under/overflow bins
						double bincontent = h1DVBF[t - 5][al]->GetBinContent(binx);
						double bincontent_SM = h1DVBF[t - 5][0]->GetBinContent(binx);
						if (Systematics != 0){
							double scalingRatio = 1;
							if (t == 5 && al == 0){
								double othercontent = hVBFRef_1D->Integral(binx, binx, "width");
								if (Systematics>0 && bincontent_SM != 0) scalingRatio = othercontent / bincontent_SM;
								else if (Systematics < 0 && othercontent != 0) scalingRatio = bincontent_SM / othercontent;
								else scalingRatio = 1;
								if (scalingRatio<0) scalingRatio = 0;
								if (scalingRatio>4) scalingRatio = 4;
								h1DVBFSigRatio->SetBinContent(binx, scalingRatio);
							}
							else scalingRatio = h1DVBFSigRatio->GetBinContent(binx);
							bincontent *= scalingRatio;
							if (al == 0) intermediateIntegral_VBF_1D[t - 5] += bincontent;
						}

						D_temp_1D[t][al]->SetBinContent(binx, bincontent);
					}
					if (tFitD > 0){
						for (int binx = 0; binx <= h2DVBF[t - 5][al]->GetNbinsX() + 1; binx++){
							for (int biny = 0; biny <= h2DVBF[t - 5][al]->GetNbinsY() + 1; biny++){
								double bincontent = h2DVBF[t - 5][al]->GetBinContent(binx, biny);
								double bincontent_SM = h2DVBF[t - 5][0]->GetBinContent(binx, biny);
								if (Systematics != 0){
									double scalingRatio = 1;
									if (t == 5 && al == 0){
										double othercontent = hVBFRef->Integral(binx, binx, biny, biny, "width");
										if (Systematics > 0 && bincontent_SM != 0) scalingRatio = othercontent / bincontent_SM;
										else if (Systematics < 0 && othercontent != 0) scalingRatio = bincontent_SM / othercontent;
										else scalingRatio = 1;
										if (scalingRatio<0) scalingRatio = 0;
										if (scalingRatio>4) scalingRatio = 4;
										h2DVBFSigRatio->SetBinContent(binx, biny, scalingRatio);
									}
									else scalingRatio = h2DVBFSigRatio->GetBinContent(binx, biny);
									bincontent *= scalingRatio;
									if (al == 0) intermediateIntegral_VBF_2D[t - 5] += bincontent;
								}

								D_temp_2D[t][al]->SetBinContent(binx, biny, bincontent);
							}
						}
					}
				}
				if (Systematics != 0){
					double myscale_1D = 1;
					if (Systematics > 0) myscale_1D = (nVBF_Sig_Reweighted / (intermediateIntegral_VBF_1D[0] * luminosity[EnergyIndex]));
					else  myscale_1D = ((nVBF_Sig_Simulated / myscale) / (intermediateIntegral_VBF_1D[0] * luminosity[EnergyIndex]));
          cout << "VBF SYSTEMATIC 1D SCALE IS " << myscale_1D << endl;
          tout << "VBF SYSTEMATIC 1D SCALE IS " << myscale_1D << endl;
          for (int al = 0; al<nAnomalousLoops; al++) D_temp_1D[t][al]->Scale(myscale_1D);
					if (t == 5) h1DVBFSigRatio->Scale(myscale_1D);
					if (tFitD>0){
						double myscale_2D = 1;
						if (Systematics > 0) myscale_2D = (nVBF_Sig_Reweighted / (intermediateIntegral_VBF_2D[0] * luminosity[EnergyIndex]));
						else  myscale_2D = ((nVBF_Sig_Simulated / myscale) / (intermediateIntegral_VBF_2D[0] * luminosity[EnergyIndex]));
            cout << "VBF SYSTEMATIC 2D SCALE IS " << myscale_2D << endl;
            tout << "VBF SYSTEMATIC 2D SCALE IS " << myscale_2D << endl;
            for (int al = 0; al < nAnomalousLoops; al++) D_temp_2D[t][al]->Scale(myscale_2D);
						if (t == 5) h2DVBFSigRatio->Scale(myscale_2D);
						if (enableDebugging){
							cout << "DEBUG STATEMENTS:\n";
							cout << "RWGT RATIO: " << h2DVBFSigRatio->GetBinContent(h2DVBFSigRatio->FindBin(220, 0.5)) << endl;
							cout << "OVERALL SCALE RATIO: " << overall_scale[t] << endl;
							cout << "myscale2D: " << myscale_2D << endl;
							cout << "Dtemp: " << D_temp_2D[t][0]->GetBinContent(D_temp_2D[t][0]->FindBin(220, 0.5)) << endl;
							int mybinx = hVBFRef->GetXaxis()->FindBin(220);
							int mybiny = hVBFRef->GetYaxis()->FindBin(0.5);
							cout << "REF: " << hVBFRef->Integral(mybinx, mybinx, mybiny, mybiny, "width") << endl;
						}
					}
				}
			}
			else overall_scale[t] = 1.0;

			//Smooths templates if desired
			for (int al = 0; al<nAnomalousLoops; al++){
				double presmoothInt = D_temp_1D[t][al]->Integral("width");
				if (isSmooth) D_temp_1D[t][al]->Smooth(1, "k3a");
        wipeOverUnderFlows(D_temp_1D[t][al]);
        double postsmoothInt = D_temp_1D[t][al]->Integral("width");
				D_temp_1D[t][al]->Scale(presmoothInt / postsmoothInt);
        if (isSmooth){
          cout << "1D SMOOTHING SCALE: " << presmoothInt / postsmoothInt << endl;
          tout << "1D SMOOTHING SCALE: " << presmoothInt / postsmoothInt << endl;
        }
        if (tFitD>0){
					presmoothInt = D_temp_2D[t][al]->Integral("width");
					if (isSmooth) D_temp_2D[t][al]->Smooth(1, "k3a");
          wipeOverUnderFlows(D_temp_2D[t][al]);
          postsmoothInt = D_temp_2D[t][al]->Integral("width");
					D_temp_2D[t][al]->Scale(presmoothInt / postsmoothInt);
          if (isSmooth){
            cout << "2D SMOOTHING SCALE: " << presmoothInt / postsmoothInt << endl;
            tout << "2D SMOOTHING SCALE: " << presmoothInt / postsmoothInt << endl;
          }
        }
			}

			//Makes unconditionally normalized PDF for backgrounds
			if (t >= 3 && tFitD > 0 && t < 5){
				for (int binx = 0; binx <= D_temp_2D[t][0]->GetNbinsX() + 1; binx++){
					double intBinX = D_temp_2D[t][0]->Integral(binx, binx, 0, D_temp_2D[t][0]->GetNbinsY() + 1);
					for (int biny = 0; biny <= D_temp_2D[t][0]->GetNbinsY() + 1; biny++){
						double bincontent = D_temp_2D[t][0]->GetBinContent(binx, biny);
						if (t == 3) hStore_qqZZ_Unconditional->SetBinContent(binx, biny, bincontent);
						if (t == 4) hStore_ZX_Unconditional->SetBinContent(binx, biny, bincontent);
						if (intBinX != 0) D_temp_2D[t][0]->SetBinContent(binx, biny, bincontent / intBinX);
					}
				}
			}

			//Stores total weights in tree for JB's smoother
			double* nTotalRecorded = new double[nAnomalousLoops];
			for (int al = 0; al < nAnomalousLoops; al++) nTotalRecorded[al] = 0;

			if (t < 5){
				nEntries = tree->GetEntries();
				for (int ev = 0; ev < nEntries; ev++){
					tree->GetEntry(ev);
					if (fitYval != fitYval) continue;
          // Protect against any KD exceeding boundaries
          if (tFitD!=0 && fitYval>=kDYarray[nbinsy]){
            cout << "Found fitYval == " << fitYval;
            fitYval = kDYarray[nbinsy] - (kDYarray[nbinsy]-kDYarray[nbinsy-1])*0.1*(ev+1.)/(nEntries+1.);
            cout << ". Fixed to " << fitYval << endl;
          }
          if (tFitD!=0 && fitYval<kDYarray[0]) fitYval = kDYarray[0] + (kDYarray[1]-kDYarray[0])*0.1*(ev+1.)/(nEntries+1.);
          if (tFitD!=0 && (fitYval>=kDYarray[nbinsy] || fitYval<kDYarray[0])) cout << "Fix has been numerically unsuccessful for " << tree->GetName() << endl;

					if ((ZZMass >= ZZMass_cut[1] || ZZMass < ZZMass_cut[0])) continue;

					double weight = MC_weight;
					if (t < 3) weight *= MC_weight_ggZZLepInt;
					if (abs(Systematics) != 1 && t < 3) weight *= MC_weight_Kfactor;
					if (Systematics == -1 && t < 3) weight *= (1.0 - ggZZ_Syst_AbsNormSyst[EnergyIndex][0])*MC_weight_Kfactor_Norm_down*tgkf->Eval(125.6);
					if (Systematics == 1 && t < 3) weight *= (1.0 + ggZZ_Syst_AbsNormSyst[EnergyIndex][0])*MC_weight_Kfactor_Norm_up*tgkf->Eval(125.6);
					if (Systematics == -2 && t < 3) weight *= (1.0 - ggZZ_Syst_AbsNormSyst[EnergyIndex][1])*MC_weight_PDF_Norm_down;
					if (Systematics == 2 && t < 3) weight *= (1.0 + ggZZ_Syst_AbsNormSyst[EnergyIndex][1])*MC_weight_PDF_Norm_up;

					//Code to add Djet reweighting for ggF
					if (t < 3){
						if (Djettag == -1){
							weight *= 1 - Djetcutshape->Eval(ZZMass);
						}
						if (Djettag == 1){
							weight *= Djetcutshape->Eval(ZZMass);
						}
					}

					// Anomalous couplings loop
					for (int al = 0; al<nAnomalousLoops; al++){
						templateWeight = weight * overall_scale[t];
						if (useAnomalousCouplings == kAddfLQ && al>0){
							templateWeight *= pow(GenHMass / mPOLE, 2 * al);
						}
						templateTree[al]->Fill();
						nTotalRecorded[al] += templateWeight;
					}
				}
			}
			else{
				TH2F* debug;
				if (enableDebugging){
					debug = (TH2F*)D_temp_2D[t][0]->Clone("debughisto");
					debug->Reset("ICESM");
				}

				int treeCode = 0;
				if (t == 5){
					//This is to account for the lack of a BSI25 sample for 8T 2e2mu
					if (EnergyIndex == 0) treeCode = 3;
					else if (EnergyIndex == 1) treeCode = t - 5;
				}
				else if (t < 8) treeCode = t - 5;
				nEntries = tree_VBF[treeCode]->GetEntries();
				for (int ev = 0; ev < nEntries; ev++){
					tree_VBF[treeCode]->GetEntry(ev);
					if (fitYval != fitYval) continue;
          // Protect against any KD exceeding boundaries
          if (tFitD!=0 && fitYval>=kDYarray[nbinsy]){
            cout << "Found fitYval == " << fitYval;
            fitYval = kDYarray[nbinsy] - (kDYarray[nbinsy]-kDYarray[nbinsy-1])*0.1*(ev+1.)/(nEntries+1.);
            cout << ". Fixed to " << fitYval << endl;
          }
          if (tFitD!=0 && fitYval<kDYarray[0]) fitYval = kDYarray[0] + (kDYarray[1]-kDYarray[0])*0.1*(ev+1.)/(nEntries+1.);
          if (tFitD!=0 && (fitYval>=kDYarray[nbinsy] || fitYval<kDYarray[0])) cout << "Fix has been numerically unsuccessful for " << tree_VBF[treeCode]->GetName() << endl;

					if (ZZMass >= ZZMass_cut[1] || ZZMass < ZZMass_cut[0]) continue;

					//Djet cut for VBF
					if (Djettag == -1 && Djet >= 0.5) continue;
					if (Djettag == 1 && Djet < 0.5) continue;

					double weight = MC_weight;
          int VBF_VH_rewgt_bin = h_VBFVH_Scale->GetXaxis()->FindBin(GenDiJetMass);
          if (VBF_VH_rewgt_bin>h_VBFVH_Scale->GetNbinsX()) VBF_VH_rewgt_bin = h_VBFVH_Scale->GetNbinsX();
          weight *= h_VBFVH_Scale->GetBinContent(VBF_VH_rewgt_bin); // VBF-VH scale for Phantom
          if (tg_VBF_unc!=0 && Systematics!=0) weight *= tg_VBF_unc->Eval(GenHMass);
/*
					if (Systematics != 0){
						double sysVBFScale = 1;
						if (tFitD == 0) sysVBFScale = h1DVBFSigRatio->GetBinContent(h1DVBFSigRatio->FindBin(ZZMass));
						else sysVBFScale = h2DVBFSigRatio->GetBinContent(h2DVBFSigRatio->FindBin(ZZMass, fitYval));
						weight *= sysVBFScale;
					}
*/
					// Anomalous couplings loop
					for (int al = 0; al<nAnomalousLoops; al++){
						templateWeight = weight * overall_scale[t];
						if (useAnomalousCouplings == kAddfLQ && al>0){
							templateWeight *= pow(GenHMass / mPOLE, 2 * al);
						}
						templateTree[al]->Fill();
						nTotalRecorded[al] += templateWeight;
						if (al == 0 && enableDebugging) debug->Fill(ZZMass, fitYval, templateWeight);
					}
				}
				if (enableDebugging){
					cout << "DEBUG STATEMENTS (TREE FILL):\n";
					cout << "RWGT RATIO: " << h2DVBFSigRatio->GetBinContent(h2DVBFSigRatio->FindBin(220, 0.5)) << endl;
					cout << "OVERALL SCALE RATIO: " << overall_scale[t] << endl;
					cout << "DEBUG HISTO: " << debug->GetBinContent(debug->FindBin(220, 0.5)) << endl;
					delete debug;
				}
			}

			for (int al = 0; al < nAnomalousLoops; al++){
        if (t != 4) cout << "RECORDED YIELD IN " << templateTree[al]->GetName() << ": " << nTotalRecorded[al] * luminosity[EnergyIndex] << endl;
        else cout << "RECORDED YIELD IN " << templateTree[al]->GetName() << ": " << nTotalRecorded[al] << endl;
        if (t != 4) tout << "RECORDED YIELD IN " << templateTree[al]->GetName() << ": " << nTotalRecorded[al] * luminosity[EnergyIndex] << endl;
        else tout << "RECORDED YIELD IN " << templateTree[al]->GetName() << ": " << nTotalRecorded[al] << endl;
        foutput->WriteTObject(templateTree[al]);
				delete templateTree[al];
			}
			delete[] nTotalRecorded;
			delete[] templateTree;
      if (Djetcutshape!=0) delete Djetcutshape;
			delete tree;
		}
		if (Systematics != 0){
			foutput->WriteTObject(h1DVBFSigRatio);
			if (tFitD != 0) foutput->WriteTObject(h2DVBFSigRatio);
		}
		delete h2DVBFSigRatio;
		delete h1DVBFSigRatio;
		for (int tr = 0; tr < 4; tr++){
			delete tree_VBF[tr];
			for (int al = 0; al < nAnomalousCouplingTemplates[useAnomalousCouplings][1]; al++){
				delete h1DVBF[tr][al];
				if (tFitD != 0) delete h2DVBF[tr][al];
			}
			delete[] h1DVBF[tr];
			if (tFitD != 0) delete[] h2DVBF[tr];
		}

		for (int al = 0; al < nAnomalousCouplingTemplates[useAnomalousCouplings][0]; al++){
			D_temp_1D[2][al]->Add(D_temp_1D[1][al], -1.0);
			D_temp_1D[2][al]->Add(D_temp_1D[0][al], -1.0);
			if (tFitD != 0){
				D_temp_2D[2][al]->Add(D_temp_2D[1][al], -1.0);
				D_temp_2D[2][al]->Add(D_temp_2D[0][al], -1.0);
			}
		}

    cout << "Integrals after everything:\nTemplate\t1D\t2D" << endl;
    tout << "Integrals after everything:\nTemplate\t1D\t2D" << endl;
		//Divides bins by Bin Width and ZX mirroring
		for (int t = 0; t < kNumTemplates; t++){
			int nAnomalousLoops = 1;
			if (t <= 2) nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][0];
			else if (t >= 5 && t <= 7) nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][1];

			for (int al = 0; al < nAnomalousLoops; al++){
				for (int binx = 0; binx < D_temp_1D[t][al]->GetNbinsX(); binx++){
          double binwidthx = D_temp_1D[t][al]->GetXaxis()->GetBinWidth(binx+1);
//          if (tFitD == 0) binwidthx = kDXarray[binx + 1] - kDXarray[binx];
//					else binwidthx = kDYarray[binx + 1] - kDYarray[binx];
					double bincontent = D_temp_1D[t][al]->GetBinContent(binx + 1);
					if (t != 3 && t != 4) bincontent /= binwidthx;
					D_temp_1D[t][al]->SetBinContent(binx + 1, bincontent);
				}
				if (tFitD != 0){
					for (int binx = 0; binx < D_temp_2D[t][al]->GetNbinsX(); binx++){
            double binwidthx = D_temp_2D[t][al]->GetXaxis()->GetBinWidth(binx+1);
//            double binwidthx = kDXarray[binx + 1] - kDXarray[binx];
						for (int biny = 0; biny < D_temp_2D[t][al]->GetNbinsY(); biny++){
              double binwidthy = D_temp_2D[t][al]->GetYaxis()->GetBinWidth(biny+1);
//              double binwidthy = kDYarray[biny + 1] - kDYarray[biny];
							double binwidth = binwidthx*binwidthy;
							double bincontent = D_temp_2D[t][al]->GetBinContent(binx + 1, biny + 1);
							if (t != 3 && t != 4) bincontent /= binwidth;
							D_temp_2D[t][al]->SetBinContent(binx + 1, biny + 1, bincontent);
						}
					}
				}
			}
			if (tFitD != 0 && t == 3){
				for (int binx = 0; binx < hStore_qqZZ_Unconditional->GetNbinsX(); binx++){
          double binwidthx = hStore_qqZZ_Unconditional->GetXaxis()->GetBinWidth(binx+1);
          //double binwidthx = kDXarray[binx + 1] - kDXarray[binx];
					for (int biny = 0; biny < hStore_qqZZ_Unconditional->GetNbinsY(); biny++){
            double binwidthy = hStore_qqZZ_Unconditional->GetYaxis()->GetBinWidth(biny+1);
            //double binwidthy = kDYarray[biny + 1] - kDYarray[biny];
						double binwidth = binwidthx*binwidthy;
						double bincontent = hStore_qqZZ_Unconditional->GetBinContent(binx + 1, biny + 1);
						bincontent /= binwidth;
						hStore_qqZZ_Unconditional->SetBinContent(binx + 1, biny + 1, bincontent);
					}
				}
				//				cout << "qqZZ unconditional template recording is complete" << endl;
			}
			if (tFitD != 0 && t == 4){
				for (int binx = 0; binx < hStore_ZX_Unconditional->GetNbinsX(); binx++){
          double binwidthx = hStore_ZX_Unconditional->GetXaxis()->GetBinWidth(binx+1);
          //double binwidthx = kDXarray[binx + 1] - kDXarray[binx];
					for (int biny = 0; biny < hStore_ZX_Unconditional->GetNbinsY(); biny++){
            double binwidthy = hStore_ZX_Unconditional->GetYaxis()->GetBinWidth(biny+1);
            //double binwidthy = kDYarray[biny + 1] - kDYarray[biny];
						double binwidth = binwidthx*binwidthy;
						double bincontent = hStore_ZX_Unconditional->GetBinContent(binx + 1, biny + 1);
						bincontent /= binwidth;
						hStore_ZX_Unconditional->SetBinContent(binx + 1, biny + 1, bincontent);
					}
				}
				//				cout << "ZX unconditional template recording is complete" << endl;
			}
			if (t == 4 && Systematics != 0){
				double intTZX_1D = D_temp_1D[t][0]->Integral(0, D_temp_1D[t][0]->GetNbinsX() + 1);
				double intTZZQQB_1D = D_temp_1D[3][0]->Integral(0, D_temp_1D[3][0]->GetNbinsX() + 1);
				for (int binx = 0; binx <= D_temp_1D[t][0]->GetNbinsX() + 1; binx++){
					double bincontent = D_temp_1D[t][0]->GetBinContent(binx);
					double bincontent_alt = D_temp_1D[3][0]->GetBinContent(binx);
					bincontent_alt *= intTZX_1D / intTZZQQB_1D;
					double difference = bincontent_alt - bincontent;

					if (Systematics > 0) bincontent += difference;
					else bincontent -= difference;
					if (bincontent < 0) bincontent = 0;
					D_temp_1D[t][0]->SetBinContent(binx, bincontent);
				}
				D_temp_1D[t][0]->Scale(intTZX_1D / D_temp_1D[t][0]->Integral(0, D_temp_1D[t][0]->GetNbinsX() + 1));
				//				cout << "ZX conditional template scaling is complete" << endl;

        intTZX_1D = D_temp_2D[t][0]->Integral();
        intTZZQQB_1D = D_temp_2D[3][0]->Integral();
        for (int binx = 0; binx <= D_temp_2D[t][0]->GetNbinsX() + 1; binx++){
					double* storeOriginal = new double[D_temp_2D[t][0]->GetNbinsY() + 2];
					for (int biny = 0; biny <= D_temp_2D[t][0]->GetNbinsY() + 1; biny++){
						double bincontent = D_temp_2D[t][0]->GetBinContent(binx, biny);
						double bincontent_alt = D_temp_2D[3][0]->GetBinContent(binx, biny);
            bincontent_alt *= intTZX_1D / intTZZQQB_1D;
            double difference = bincontent_alt - bincontent;
						storeOriginal[biny] = bincontent;

						if (Systematics > 0) bincontent += difference;
						else bincontent -= difference;
						if (bincontent < 0) bincontent = 0;
						D_temp_2D[t][0]->SetBinContent(binx, biny, bincontent);
					}
					double intBinX = D_temp_2D[t][0]->Integral(binx, binx, 0, D_temp_2D[t][0]->GetNbinsY() + 1);
					for (int biny = 0; biny <= D_temp_2D[t][0]->GetNbinsY() + 1; biny++){
						double bincontent = D_temp_2D[t][0]->GetBinContent(binx, biny);

						if (intBinX != 0){
							D_temp_2D[t][0]->SetBinContent(binx, biny, bincontent / intBinX);
							double sysRatio = 0;
							if (storeOriginal[biny] != 0) sysRatio = (bincontent / intBinX) / (storeOriginal[biny]);
							double unconditionalbincontent = hStore_ZX_Unconditional->GetBinContent(binx, biny);
							hStore_ZX_Unconditional->SetBinContent(binx, biny, unconditionalbincontent*sysRatio);
						}
						else{
							hStore_ZX_Unconditional->SetBinContent(binx, biny, 0);
						}
					}
					delete[] storeOriginal;
				}
				//				cout << "ZX unconditional template mirroring is complete" << endl;

				if (Systematics < 0){
					double presmoothInt = D_temp_1D[t][0]->Integral();
					if (isSmooth) D_temp_1D[t][0]->Smooth(1, "k3a");
					double postsmoothInt = D_temp_1D[t][0]->Integral();
					D_temp_1D[t][0]->Scale(presmoothInt / postsmoothInt);
				}
				//				cout << "ZX conditional template smoothing is complete" << endl;
			}
			for (int al = 0; al < nAnomalousLoops; al++){
				foutput->WriteTObject(D_temp_1D[t][al]);
				if (tFitD != 0) foutput->WriteTObject(D_temp_2D[t][al]);
			}
			if (tFitD != 0 && t == 3) foutput->WriteTObject(hStore_qqZZ_Unconditional);
			if (tFitD != 0 && t == 4) foutput->WriteTObject(hStore_ZX_Unconditional);

			for (int al = 0; al < nAnomalousLoops; al++){
        cout << templatenames[t] << " (anom. coupl.: " << al << ") integrals:\t";
        tout << templatenames[t] << " (anom. coupl.: " << al << ") integrals:\t";
				if (t != 3 && t != 4){
					cout << D_temp_1D[t][al]->Integral("width")*luminosity[EnergyIndex] << '\t';
					if (tFitD != 0) cout << D_temp_2D[t][al]->Integral("width")*luminosity[EnergyIndex] << endl;
					else cout << endl;
          tout << D_temp_1D[t][al]->Integral("width")*luminosity[EnergyIndex] << '\t';
          if (tFitD != 0) tout << D_temp_2D[t][al]->Integral("width")*luminosity[EnergyIndex] << endl;
          else tout << endl;
        }
				else{
          cout << D_temp_1D[t][al]->Integral() << '\t';
          if (tFitD != 0) cout << D_temp_2D[t][al]->Integral() << endl;
          else cout << endl;
          tout << D_temp_1D[t][al]->Integral() << '\t';
          if (tFitD != 0) tout << D_temp_2D[t][al]->Integral() << endl;
          else tout << endl;
        }
			}
		}
		if (tFitD != 0){
			cout << "Unconditional Integrals are: " << endl;
			cout << "qqZZ: " << hStore_qqZZ_Unconditional->Integral("width")*luminosity[EnergyIndex] << endl;
      cout << "Z+X: " << hStore_ZX_Unconditional->Integral("width")*luminosity[EnergyIndex] << endl;
      tout << "Unconditional Integrals are: " << endl;
      tout << "qqZZ: " << hStore_qqZZ_Unconditional->Integral("width")*luminosity[EnergyIndex] << endl;
      tout << "Z+X: " << hStore_ZX_Unconditional->Integral("width")*luminosity[EnergyIndex] << endl;
    }

		delete hStore_qqZZ_Unconditional;
		delete hStore_ZX_Unconditional;
		for (int t = 0; t < kNumTemplates; t++){
			int nAnomalousLoops = 1;
			if (t <= 2) nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][0];
			else if (t >= 5 && t <= 7) nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][1];
			for (int al = 0; al < nAnomalousLoops; al++){
				delete D_temp_1D[t][al];
				if (tFitD != 0) delete D_temp_2D[t][al];
			}
			delete[] D_temp_1D[t];
			if (tFitD != 0) delete[] D_temp_2D[t];
		}

		finput_VBFRef->Close();
		foutput->Close();
    tout.close();
	}
	delete tgkf;
  delete tg_VBF_unc;
	finput_KDFactor->Close();
  delete tsp_VBF_Sig;
	finput_VBF->Close();
  Djetcutfile->Close();
  delete h_VBFVH_Scale;
  finput_vbfvhscale->Close();
	delete tg_interf;
}

double doLinearCombination(double first, double c_first, double second, double c_second, double bkg, int outputtype){
	double aa[3] = { first, c_first, sqrt(c_first) };
	double bb[3] = { second, c_second, sqrt(c_second) };

	double sig = (aa[0] * bb[2] - aa[2] * bb[0] + (aa[2] - bb[2])*bkg) / (aa[1] * bb[2] - aa[2] * bb[1]);
	double interf = (aa[0] * bb[1] - aa[1] * bb[0] + (aa[1] - bb[1])*bkg) / (aa[2] * bb[1] - aa[1] * bb[2]);
	if (sig < 0){ sig = 0; interf = 0; };

	if (outputtype == kIntHist) return interf;
	else if(outputtype==kSigHist) return sig;
	else return 0;
}

double doPeakCombination(double first, double c_first, double second, double c_second, double bkg, int outputtype){
	double aa[3] = { first, 1, sqrt(c_first) };
	double bb[3] = { second, 1, sqrt(c_second) };

	double sig = (aa[0] * bb[2] - aa[2] * bb[0] + (aa[2] - bb[2])*bkg) / (aa[1] * bb[2] - aa[2] * bb[1]);
	double interf = (aa[0] * bb[1] - aa[1] * bb[0] + (aa[1] - bb[1])*bkg) / (aa[2] * bb[1] - aa[1] * bb[2]);
	if (sig < 0){ sig = 0; interf = 0; };

	cout << "VBF peak as calculated\n";
	cout << "Signal: " << sig << "\tInterf: " << interf << endl;
	
	if (outputtype == kIntHist) return interf;
	else if(outputtype==kSigHist) return sig;
	else return 0;
}


void oneDlinearcombination(TH1F* first, int firsttype, TH1F* second, int secondtype, TH1F* input, int inputtype, TH1F* finaloutput, int outputtype, TH1F* finaloutput2, int output2type){
	TH1F* output = (TH1F*) input->Clone();
	TH1F* output2 = (TH1F*) output->Clone();
	if(outputtype==kIntHist || outputtype==kSigHist){
		if(inputtype==kBSI25Hist && firsttype==kSigHist && secondtype==kBkgHist){
			output->Add(first, -25.0);
			output->Add(second, -1.0);
			output->Scale(0.2);
			if (outputtype == kSigHist){ delete output; output = (TH1F*) first->Clone(); }

			output2->Add(first, -25.0);
			output2->Add(second, -1.0);
			output2->Scale(0.2);
			if (output2type == kSigHist){ delete output2; output2 = (TH1F*) first->Clone(); }
		}
		if(inputtype==kBSI25Hist && firsttype==kBSIHist && secondtype==kBkgHist){
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				double bsi = first->GetBinContent(binx);
				double bkg = second->GetBinContent(binx);
				double bsi25 = output->GetBinContent(binx);

				double weight=doLinearCombination(bsi25,25,bsi,1,bkg,outputtype);
				output->SetBinContent(binx,weight);
				if (finaloutput2 != 0){
					double weight2 = doLinearCombination(bsi25, 25, bsi, 1, bkg, output2type);
					output2->SetBinContent(binx, weight2);
				}
			}
		}		
		if(inputtype==kBSI25Hist && firsttype==kBSI10Hist && secondtype==kBkgHist){
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				double bsi10 = first->GetBinContent(binx);
				double bkg = second->GetBinContent(binx);
				double bsi25 = output->GetBinContent(binx);

				double weight=doLinearCombination(bsi25,25,bsi10,10,bkg,outputtype);
				output->SetBinContent(binx,weight);
				if (finaloutput2 != 0){
					double weight2 = doLinearCombination(bsi25, 25, bsi10, 10, bkg, output2type);
					output2->SetBinContent(binx, weight2);
				}
			}
		}
		if(inputtype==kBSI10Hist && firsttype==kBSIHist && secondtype==kBkgHist){
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				double bsi = first->GetBinContent(binx);
				double bkg = second->GetBinContent(binx);
				double bsi10 = output->GetBinContent(binx);

				double weight=doLinearCombination(bsi10,10,bsi,1,bkg,outputtype);
				output->SetBinContent(binx,weight);
				if (finaloutput2 != 0){
					double weight2 = doLinearCombination(bsi10, 10, bsi, 1, bkg, output2type);
					output2->SetBinContent(binx, weight2);
				}
			}
		}
		for (int binx = 1; binx <= output->GetNbinsX(); binx++){
			finaloutput->SetBinContent(binx,output->GetBinContent(binx));
			if (finaloutput2 != 0) finaloutput2->SetBinContent(binx,output2->GetBinContent(binx));
		}
	}
	else{cout<<"Option not yet supported. Exiting..."<<endl; assert(0);};
	delete output;
	delete output2;
}

void twoDlinearcombination(TH2F* first, int firsttype, TH2F* second, int secondtype, TH2F* input, int inputtype, TH2F* finaloutput, int outputtype, TH2F* finaloutput2, int output2type){
	TH2F* output = (TH2F*) input->Clone();
	TH2F* output2 = (TH2F*) output->Clone();
	if(outputtype==kIntHist || outputtype==kSigHist){
		if(inputtype==kBSI25Hist && firsttype==kSigHist && secondtype==kBkgHist){
			output->Add(first, -25.0);
			output->Add(second, -1.0);
			output->Scale(0.2);
			if (outputtype == kSigHist){ delete output; output = (TH2F*) first->Clone(); }

			output2->Add(first, -25.0);
			output2->Add(second, -1.0);
			output2->Scale(0.2);
			if (output2type == kSigHist){ delete output2; output2 = (TH2F*) first->Clone(); }
		}
		if(inputtype==kBSI25Hist && firsttype==kBSIHist && secondtype==kBkgHist){
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				for (int biny = 1; biny <= output->GetNbinsY(); biny++){
					double bsi = first->GetBinContent(binx,biny);
					double bkg = second->GetBinContent(binx,biny);
					double bsi25 = output->GetBinContent(binx,biny);

					double weight=doLinearCombination(bsi25,25,bsi,1,bkg,outputtype);
					output->SetBinContent(binx,biny,weight);
					if (finaloutput2 != 0){
						double weight2 = doLinearCombination(bsi25, 25, bsi, 1, bkg, output2type);
						output2->SetBinContent(binx, biny, weight2);
					}
				}
			}
		}		
		if(inputtype==kBSI25Hist && firsttype==kBSI10Hist && secondtype==kBkgHist){
			//double scaleval = 1./(-50. + 25.*sqrt(10.));
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				for (int biny = 1; biny <= output->GetNbinsY(); biny++){
					double bsi10 = first->GetBinContent(binx,biny);
					double bkg = second->GetBinContent(binx,biny);
					double bsi25 = output->GetBinContent(binx,biny);

					double weight=doLinearCombination(bsi25,25,bsi10,10,bkg,outputtype);
					output->SetBinContent(binx,biny,weight);
					if (finaloutput2 != 0){
						double weight2 = doLinearCombination(bsi25, 25, bsi10, 10, bkg, output2type);
						output2->SetBinContent(binx, biny, weight2);
					}
				}
			}
		}
		if(inputtype==kBSI10Hist && firsttype==kBSIHist && secondtype==kBkgHist){
			//double scaleval = 1./(10 - sqrt(10));
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				for (int biny = 1; biny <= output->GetNbinsY(); biny++){
					double bsi = first->GetBinContent(binx,biny);
					double bkg = second->GetBinContent(binx,biny);
					double bsi10 = output->GetBinContent(binx,biny);

					double weight=doLinearCombination(bsi10,10,bsi,1,bkg,outputtype);
					output->SetBinContent(binx,biny,weight);
					if (finaloutput2 != 0){
						double weight2 = doLinearCombination(bsi10, 10, bsi, 1, bkg, output2type);
						output2->SetBinContent(binx, biny, weight2);
					}
				}
			}
		}
		for (int binx = 1; binx <= output->GetNbinsX(); binx++){
			for (int biny = 1; biny <= output->GetNbinsY(); biny++){
				finaloutput->SetBinContent(binx, biny, output->GetBinContent(binx, biny));
				if (finaloutput2 != 0) finaloutput2->SetBinContent(binx,biny,output2->GetBinContent(binx,biny));
			}
		}
	}
	else{cout<<"Option not yet supported. Exiting..."<<endl; assert(0);};
	delete output;
	delete output2;
}

void threeDlinearcombination(TH3F* first, int firsttype, TH3F* second, int secondtype, TH3F* input, int inputtype, TH3F* finaloutput, int outputtype, TH3F* finaloutput2, int output2type){
	TH3F* output = (TH3F*) input->Clone();
	TH3F* output2 = (TH3F*) output->Clone();
	if(outputtype==kIntHist || outputtype==kSigHist){
		if(inputtype==kBSI25Hist && firsttype==kSigHist && secondtype==kBkgHist){
			output->Add(first, -25.0);
			output->Add(second, -1.0);
			output->Scale(0.2);
			if (outputtype == kSigHist){ delete output; output = (TH3F*) first->Clone(); }

			output2->Add(first, -25.0);
			output2->Add(second, -1.0);
			output2->Scale(0.2);
			if (output2type == kSigHist){ delete output2; output2 = (TH3F*) first->Clone(); }
		}
		if(inputtype==kBSI25Hist && firsttype==kBSIHist && secondtype==kBkgHist){
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				for (int biny = 1; biny <= output->GetNbinsY(); biny++){
					for (int binz = 1; binz <= output->GetNbinsZ(); binz++){
						double bsi = first->GetBinContent(binx, biny, binz);
						double bkg = second->GetBinContent(binx, biny, binz);
						double bsi25 = output->GetBinContent(binx, biny, binz);

						double weight = doLinearCombination(bsi25, 25, bsi, 1, bkg, outputtype);
						output->SetBinContent(binx, biny, binz, weight);
						if (finaloutput2 != 0){
							double weight2 = doLinearCombination(bsi25, 25, bsi, 1, bkg, output2type);
							output2->SetBinContent(binx, biny, binz, weight2);
						}
					}
				}
			}
		}		
		if(inputtype==kBSI25Hist && firsttype==kBSI10Hist && secondtype==kBkgHist){
			//double scaleval = 1./(-50. + 25.*sqrt(10.));
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				for (int biny = 1; biny <= output->GetNbinsY(); biny++){
					for (int binz = 1; binz <= output->GetNbinsZ(); binz++){
						double bsi10 = first->GetBinContent(binx, biny, binz);
						double bkg = second->GetBinContent(binx, biny, binz);
						double bsi25 = output->GetBinContent(binx, biny, binz);

						double weight = doLinearCombination(bsi25, 25, bsi10, 10, bkg, outputtype);
						output->SetBinContent(binx, biny, binz, weight);
						if (finaloutput2 != 0){
							double weight2 = doLinearCombination(bsi25, 25, bsi10, 10, bkg, output2type);
							output2->SetBinContent(binx, biny, binz, weight2);
						}
					}
				}
			}
		}
		if(inputtype==kBSI10Hist && firsttype==kBSIHist && secondtype==kBkgHist){
			//double scaleval = 1./(10 - sqrt(10));
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				for (int biny = 1; biny <= output->GetNbinsY(); biny++){
					for (int binz = 1; binz <= output->GetNbinsZ(); binz++){
						double bsi = first->GetBinContent(binx, biny, binz);
						double bkg = second->GetBinContent(binx, biny, binz);
						double bsi10 = output->GetBinContent(binx, biny, binz);

						double weight = doLinearCombination(bsi10, 10, bsi, 1, bkg, outputtype);
						output->SetBinContent(binx, biny, binz, weight);
						if (finaloutput2 != 0){
							double weight2 = doLinearCombination(bsi10, 10, bsi, 1, bkg, output2type);
							output2->SetBinContent(binx, biny, binz, weight2);
						}
					}
				}
			}
		}
		for (int binx = 1; binx <= output->GetNbinsX(); binx++){
			for (int biny = 1; biny <= output->GetNbinsY(); biny++){
				for (int binz = 1; binz <= output->GetNbinsZ(); binz++){
					finaloutput->SetBinContent(binx, biny, binz, output->GetBinContent(binx, biny, binz));
					if (finaloutput2 != 0) finaloutput2->SetBinContent(binx,biny,binz,output2->GetBinContent(binx,biny,binz));
				}
			}
		}
	}
	else{cout<<"Option not yet supported. Exiting..."<<endl; assert(0);};
	delete output;
	delete output2;
}

void wipeOverUnderFlows(TH1F* hwipe){
  double integral = hwipe->Integral(0, hwipe->GetNbinsX()+1);
  for (int binx=0; binx<=hwipe->GetNbinsX()+1; binx++){
    if (binx>=1 && binx<=hwipe->GetNbinsX()) continue;
    if (hwipe->GetBinContent(binx)!=0){
      hwipe->SetBinContent(binx, 0);
      cout << hwipe->GetName() << " binX = " << binx << " non-zero." << endl;
    }
  }
  double wipeScale = hwipe->Integral();
  wipeScale = integral / wipeScale;
  hwipe->Scale(wipeScale);
}
void wipeOverUnderFlows(TH2F* hwipe){
  double integral = hwipe->Integral(0, hwipe->GetNbinsX()+1, 0, hwipe->GetNbinsY()+1);
  for (int binx=0; binx<=hwipe->GetNbinsX()+1; binx++){
    if (binx>=1 && binx<=hwipe->GetNbinsX()) continue;
    for (int biny=0; biny<=hwipe->GetNbinsY()+1; biny++){
      if (biny>=1 && biny<=hwipe->GetNbinsY()) continue;
      if (hwipe->GetBinContent(binx, biny)!=0){
        hwipe->SetBinContent(binx, biny, 0);
        cout << hwipe->GetName() << " binX = " << binx << " binY = " << biny << " non-zero." << endl;
      }
    }
  }
  double wipeScale = hwipe->Integral();
  wipeScale = integral / wipeScale;
  hwipe->Scale(wipeScale);
}
void wipeOverUnderFlows(TH3F* hwipe){
  double integral = hwipe->Integral(0, hwipe->GetNbinsX()+1, 0, hwipe->GetNbinsY()+1, 0, hwipe->GetNbinsZ()+1);
  for (int binx=0; binx<=hwipe->GetNbinsX()+1; binx++){
    if (binx>=1 && binx<=hwipe->GetNbinsX()) continue;
    for (int biny=0; biny<=hwipe->GetNbinsY()+1; biny++){
      if (biny>=1 && biny<=hwipe->GetNbinsY()) continue;
      for (int binz=0; binz<=hwipe->GetNbinsZ()+1; binz++){
        if (binz>=1 && binz<=hwipe->GetNbinsZ()) continue;
        if (hwipe->GetBinContent(binx, biny, binz)!=0){
          hwipe->SetBinContent(binx, biny, binz, 0);
          cout << hwipe->GetName() << " binX = " << binx << " binY = " << biny << " binZ = " << binz << " non-zero." << endl;
        }
      }
    }
  }
  double wipeScale = hwipe->Integral();
  wipeScale = integral / wipeScale;
  hwipe->Scale(wipeScale);
}


void progressbar(int val, int tot){
  int percent=floor(0.01*tot);
  if(percent==0) percent=1;

  if(val%percent==0 && val!=tot){
    cout<<"[ "<<setw(3)<<val/percent<<"% |";
    for(int k=1;k<val/percent;k++) cout<<"=";
    if(val%percent!=100) cout<<">";
    for(int k=val/percent+1;k<100;k++) cout<<" ";
    cout<<"| ]";
    fflush(stdout);
    putchar('\r');
  }
  else if(val==tot){
    cout<<"[ 100% |";
    for(int k=0;k<100;k++) cout<<"=";
    cout<<"| ]";
    fflush(stdout);
    putchar('\r');        
  }
}
