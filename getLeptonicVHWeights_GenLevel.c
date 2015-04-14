#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TMath.h"
#include "TSpline.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "TDirectory.h"
#include "./data/ZZ4l_125_6_Samples.h"

using namespace std;

//Initializers
bool enableDebugging = false;
bool useDjettagging=true;
enum histtypes{kSigHist,kBkgHist,kIntHist,kBSIHist,kBSI10Hist,kBSI25Hist};
void produce_MJJdistributions_one(int folder, int erg_tev);
void plotVBFVHReweighted();
double doLinearCombination(double first, double c_first, double second, double c_second, double bkg, int outputtype);
double doPeakCombination(double first, double c_first, double second, double c_second, double bkg, int outputtype);
void oneDlinearcombination(TH1F* first, int firsttype, TH1F* second, int secondtype, TH1F* input, int inputtype, TH1F* finaloutput, int outputtype, TH1F* finaloutput2=0, int output2type=-99);
void twoDlinearcombination(TH2F* first, int firsttype, TH2F* second, int secondtype, TH2F* input, int inputtype, TH2F* finaloutput, int outputtype, TH2F* finaloutput2=0, int output2type=-99);
void threeDlinearcombination(TH3F* first, int firsttype, TH3F* second, int secondtype, TH3F* input, int inputtype, TH3F* finaloutput, int outputtype, TH3F* finaloutput2=0, int output2type=-99);
void progressbar(int val, int tot);


//Main Function, runs over all desired iterations
void produce_MJJdistributions(){
  for (int CoM=7; CoM<9; ++CoM){
    for (int channel=0; channel<3; ++channel){
      produce_MJJdistributions_one(channel, CoM);
    }
  }
}
void getLeptonicVHWeights_GenLevel(){
  produce_MJJdistributions();
  plotVBFVHReweighted();
}
//Function to build one template
// folder = 0,1,2 (final state corresponds to 4mu,4e,2mu2e respectively)
// erg_tev = 7,8 (CoM energy)
void produce_MJJdistributions_one(int folder, int erg_tev){
  gStyle->SetOptStat(0);

  char TREE_NAME[] = "SelectedTree";
	TString OUTPUT_NAME = "HtoZZ4l_Phantom_125p6_VHDistributions_GenLevel";
  TString comstring;
  comstring.Form("%iTeV", erg_tev);
  TString gencomstring;
  gencomstring.Form("%iT/", erg_tev);
  TString erg_dir;
	erg_dir.Form("LHC_%iTeV/", erg_tev);
  OUTPUT_NAME += ".root";

  TString sample_VBF_suffix[4] ={
    "BSI_VBF_Phantom",
    "Bkg_VBF_Phantom",
    "BSI10_VBF_Phantom",
    "BSI25_VBF_Phantom"
  };
  TString sample_VBF_gensuffix[4] ={
    "BSI",
    "Bkg",
    "BSI10",
    "BSI25"
  };

  double xsec_VBF_Phantom[2][3][4]={
    {
      { 2.053859553143833e-1, 1.611617374456183e-1, 2.187258636589030e-1, 2.505576749660674e-1 },
      { 2.050597663931831e-1, 1.610671663387286e-1, 2.178041943144518e-1, 2.502592444436491e-1 },
      { 4.512172284721055e-1, 3.718595894117134e-1, 4.512172284721055e-1, 5.417040330546889e-1 }
    },
    {
      { 2.487975936891337e-1, 1.951954743601816e-1, 2.678337815958082e-1, 3.117289755556746e-1 },
      { 2.484688144236787e-1, 1.952187837277965e-1, 2.680082365969152e-1, 3.117402751925978e-1 },
      { 5.434945696873235e-1, 4.447269876857436e-1, 5.814292186114383e-1, 6.707175001136875e-1 }
    }
  };

  TString cinput_common = "/scratch0/hep/ianderso/LHEFiles/VBF/Phantom/" + gencomstring;

	int EnergyIndex = 1;
	if (erg_tev == 7) EnergyIndex = 0;
	float lowside[3] = { 220, 230, 240 };
	double mPOLE = 125.6;
	float ZZMass_PeakCut[2] = { 125.1, 126.1 }; // Spin 0 analysis

	float templateWeight = 1;
  float GenDiJetMass;
  float MC_weight;
	float ZZMass = 0;

	double nVBFPeak[4] = { 0 };

	for (int lo = 0; lo < 1; lo++){
    TString coutput_common = user_TemplateswithTrees_dir + "../VHContributions/" + erg_dir;
		coutput_common += user_folder[folder] + "/";
		gSystem->Exec("mkdir -p " + coutput_common);

		cout << "===============================" << endl;
		cout << "CoM Energy: " << erg_tev << " TeV" << endl;
		cout << "Decay Channel: " << user_folder[folder] << endl;
		cout << endl;

		TString coutput = coutput_common + OUTPUT_NAME;
		TFile* foutput = new TFile(coutput, "recreate");

    foutput->cd();

		float ZZMass_cut[2] = { lowside[lo], 1600 };
		float ZZwidth = 20.0;
		const int nbinsx = (ZZMass_cut[1] - ZZMass_cut[0]) / ZZwidth;
		float kDXarray[nbinsx + 3];
    kDXarray[0] = ZZMass_PeakCut[0];
    kDXarray[1] = ZZMass_PeakCut[1];
		for (int bin = 2; bin < nbinsx + 3; bin++){
			kDXarray[bin] = ZZMass_cut[0] + ZZwidth*(bin-2);
		}
    float beginWMass = floor(MG_W_Phantom[0]-5.*MG_W_Phantom[1]); // = 70 GeV
    float endZMass = ceil(MG_Z_Phantom[0]+5.*MG_Z_Phantom[1]); // = 104 GeV
//    float WZthreshold = 0.5*(floor(MG_W_Phantom[0]+(MG_Z_Phantom[0]-MG_W_Phantom[0])*MG_W_Phantom[1]/(MG_W_Phantom[1]+MG_Z_Phantom[1]))+ceil(MG_W_Phantom[0]+(MG_Z_Phantom[0]-MG_W_Phantom[0])*MG_W_Phantom[1]/(MG_W_Phantom[1]+MG_Z_Phantom[1]))); // = 85.5 GeV
    float WZthreshold = ceil(MG_W_Phantom[0]+(MG_Z_Phantom[0]-MG_W_Phantom[0])*MG_W_Phantom[1]/(MG_W_Phantom[1]+MG_Z_Phantom[1])); // = 86 GeV
    float sidebandBeginMass = beginWMass-(endZMass-beginWMass)*0.5; // = 53 GeV
    float sidebandEndMass = endZMass+(endZMass-beginWMass)*0.5; // = 121 GeV
    float beginLowHHMass = min((float) 50., sidebandBeginMass);
    float beginHHMass = floor(mPOLE-0.5);
    float endHHMass = floor(mPOLE+0.5);
    float sidebandBeginHHMass = floor(beginHHMass - (endHHMass-beginHHMass)*0.5);
    float sidebandEndHHMass = ceil(endHHMass + (endHHMass-beginHHMass)*0.5);
    cout << "Mass thresholds: "
      << beginLowHHMass << '\t'
      << sidebandBeginMass << '\t'
      << beginWMass << '\t'
      << WZthreshold << '\t'
      << endZMass << '\t'
      << sidebandEndMass << '\t'
      << sidebandBeginHHMass << '\t'
      << beginHHMass << '\t'
      << endHHMass << '\t'
      << sidebandEndHHMass << '\t'
      << endl;

//    const int nbinsy = 260;
    const int nbinsy = 130;
    //    const int nbinsy = 261;
    float kDYarray[nbinsy + 2];
    float kDY_bounds[2] ={ 0, 130 };
 //   float kDY_bounds[2] ={ 0, 130.5 };
    for (int bin = 0; bin < nbinsy + 2; bin++){
			double binwidth = (kDY_bounds[1] - kDY_bounds[0]) / nbinsy;
			kDYarray[bin] = kDY_bounds[0] + binwidth*bin;
		}

		const int kNumTemplates = 3;
		double overall_scale[kNumTemplates] = { 1 };
    TString templatetitles[kNumTemplates+1] ={ "VBF Sig", "VBF Bkg", "VBF Int", "VBF extra" };
    TString templatenames[kNumTemplates+1] ={ "VBF_Sig", "VBF_Bkg", "VBF_Int", "VBF_extra" };
    // Build template structure, including anomalous couplings

		int nEntries;
		MC_weight = 1;
		ZZMass = 0;

		//Initialize and grab each of the four Phantom trees
    TChain* tree_VBF[4];
    TH2F* h2DVBF[4];
    TH1F* hVBF_onshell[4];
    TH1F* hVBF_offshell[4];
    TH1F* hVBF_onshell_LC[4];
    TH1F* hVBF_offshell_LC[4];
    TH1F* hVBF_onshell_scaled[4];
    TH1F* hVBF_offshell_scaled[4];
    TH1F* hVBF_onshell_scaled_wVBF[4];
    TH1F* hVBF_offshell_scaled_wVBF[4];

    for (int tr = 0; tr < 4; tr++){
      /*
      0: BSI
      1: BKG
      2: BSI10
      3: BSI25
      */
      tree_VBF[tr] = new TChain(TREE_NAME);
      for (int decay = 0; decay < 3; decay++){
        if (decay!=folder) continue;
        TString channelname = user_folder[decay];
        if (decay == 2) channelname = "2e2mu";

        TString cinput_VBF = cinput_common + sample_VBF_gensuffix[tr] + "/" + channelname + "/phantom_" + comstring + "_" + channelname + "_" + sample_VBF_gensuffix[tr] + "_false.root";
				tree_VBF[tr]->Add(cinput_VBF);
			}

      tree_VBF[tr]->SetBranchAddress("DiJetMass", &GenDiJetMass);
      tree_VBF[tr]->SetBranchAddress("ZZMass", &ZZMass);


      foutput->cd();
      TString templatename_2D_core;
      TString templatename_2D;

      templatename_2D_core = "VBF_";
      templatename_2D_core += templatenames[tr];

      templatename_2D = templatename_2D_core;
      h2DVBF[tr] = new TH2F(templatename_2D, templatetitles[tr], nbinsx+2, kDXarray, nbinsy+1, kDYarray);
      hVBF_onshell[tr] = new TH1F(Form("%s_onshell_raw", h2DVBF[tr]->GetName()), templatetitles[tr], nbinsy+1, kDYarray);
      hVBF_offshell[tr] = new TH1F(Form("%s_offshell_raw", h2DVBF[tr]->GetName()), templatetitles[tr], nbinsy+1, kDYarray);
      hVBF_onshell_LC[tr] = new TH1F(Form("%s_onshell_LC", h2DVBF[tr]->GetName()), templatetitles[tr], nbinsy+1, kDYarray);
      hVBF_offshell_LC[tr] = new TH1F(Form("%s_offshell_LC", h2DVBF[tr]->GetName()), templatetitles[tr], nbinsy+1, kDYarray);
      hVBF_onshell_scaled[tr] = new TH1F(Form("%s_onshell_scaled", h2DVBF[tr]->GetName()), templatetitles[tr], nbinsy+1, kDYarray);
      hVBF_offshell_scaled[tr] = new TH1F(Form("%s_offshell_scaled", h2DVBF[tr]->GetName()), templatetitles[tr], nbinsy+1, kDYarray);
      hVBF_onshell_scaled_wVBF[tr] = new TH1F(Form("%s_onshell_scaled_wVBF", h2DVBF[tr]->GetName()), templatetitles[tr], nbinsy+1, kDYarray);
      hVBF_offshell_scaled_wVBF[tr] = new TH1F(Form("%s_offshell_scaled_wVBF", h2DVBF[tr]->GetName()), templatetitles[tr], nbinsy+1, kDYarray);

      int nVBFEntries = tree_VBF[tr]->GetEntries();
      double sumFilled=0;
      MC_weight = xsec_VBF_Phantom[EnergyIndex][folder][tr];
      for (int ev = 0; ev < nVBFEntries; ev++){
        progressbar(ev, tree_VBF[tr]->GetEntries());
        tree_VBF[tr]->GetEntry(ev);

        if (ZZMass>=kDXarray[nbinsx + 2] || ZZMass<kDXarray[0]) continue;
        if (GenDiJetMass>=kDYarray[nbinsy]) GenDiJetMass = (kDYarray[nbinsy]+kDYarray[nbinsy+1])/2.;

        double weight = MC_weight / nVBFEntries;

        // Anomalous couplings loop
        double fillWeight = weight;
        h2DVBF[tr]->Fill(ZZMass, GenDiJetMass, fillWeight);
        sumFilled += fillWeight;
      }
      cout << endl;
      cout << "Filled " << sumFilled*luminosity[EnergyIndex] << endl;

      h2DVBF[tr]->Scale(luminosity[EnergyIndex]);

      foutput->cd();
      for (int biny=1; biny<=h2DVBF[tr]->GetNbinsY(); biny++){
        double sumOffshell=h2DVBF[tr]->Integral(3, nbinsx + 2, biny, biny);
        double sumOnshell=h2DVBF[tr]->Integral(1, 1, biny, biny);
        hVBF_onshell[tr]->SetBinContent(biny, sumOnshell);
        hVBF_offshell[tr]->SetBinContent(biny, sumOffshell);
        hVBF_onshell_LC[tr]->SetBinContent(biny, sumOnshell);
        hVBF_offshell_LC[tr]->SetBinContent(biny, sumOffshell);
        hVBF_onshell_scaled[tr]->SetBinContent(biny, sumOnshell);
        hVBF_offshell_scaled[tr]->SetBinContent(biny, sumOffshell);
        hVBF_onshell_scaled_wVBF[tr]->SetBinContent(biny, sumOnshell);
        hVBF_offshell_scaled_wVBF[tr]->SetBinContent(biny, sumOffshell);
      }
      cout << "On-shell integral: " << hVBF_onshell[tr]->Integral() << endl;
      cout << "Off-shell integral: " << hVBF_offshell[tr]->Integral() << endl;

      delete tree_VBF[tr];
    }

		//Make VBF Sig/Int from linear combinations of above templates
		//0: VBF Sig
		//2: VBF Int
		//	 For 7 TeV samples, BSI25, Bkg, and BSI10 are used
		//	 For 8 TeV samples, BSI, Bkg, and BSI10 are used
		if (EnergyIndex == 0){
      oneDlinearcombination(hVBF_offshell_LC[2], kBSI10Hist, hVBF_offshell_LC[1], kBkgHist, hVBF_offshell_LC[3], kBSI25Hist, hVBF_offshell_LC[0], kSigHist, hVBF_offshell_LC[2], kIntHist);

      for (int binx=1; binx<=hVBF_onshell_LC[0]->GetNbinsX(); binx++){
        for (int tr = 0; tr < 4; tr++) nVBFPeak[tr] = hVBF_onshell_LC[tr]->GetBinContent(binx);
        double nVBFsigtemp = doPeakCombination(nVBFPeak[2], 10, nVBFPeak[3], 25, nVBFPeak[1], kSigHist);
        double nVBFinterftemp = doPeakCombination(nVBFPeak[2], 10, nVBFPeak[3], 25, nVBFPeak[1], kIntHist);
        nVBFPeak[0] = nVBFsigtemp;
        nVBFPeak[2] = nVBFinterftemp;
        for (int tr = 0; tr < 4; tr++) hVBF_onshell_LC[tr]->SetBinContent(binx, nVBFPeak[tr]);
      }
    }
		else if (EnergyIndex == 1){
      oneDlinearcombination(hVBF_offshell_LC[0], kBSIHist, hVBF_offshell_LC[1], kBkgHist, hVBF_offshell_LC[2], kBSI10Hist, hVBF_offshell_LC[0], kSigHist, hVBF_offshell_LC[2], kIntHist);

      for (int binx=1; binx<=hVBF_onshell_LC[0]->GetNbinsX(); binx++){
        for (int tr = 0; tr < 4; tr++) nVBFPeak[tr] = hVBF_onshell_LC[tr]->GetBinContent(binx);
        double nVBFsigtemp = doPeakCombination(nVBFPeak[2], 10, nVBFPeak[0], 1, nVBFPeak[1], kSigHist);
        double nVBFinterftemp = doPeakCombination(nVBFPeak[2], 10, nVBFPeak[0], 1, nVBFPeak[1], kIntHist);
        nVBFPeak[0] = nVBFsigtemp;
        nVBFPeak[2] = nVBFinterftemp;
        for (int tr = 0; tr < 4; tr++) hVBF_onshell_LC[tr]->SetBinContent(binx, nVBFPeak[tr]);
      }
    }
    cout << "Finished peak combination" << endl;
    int beginSB = hVBF_onshell_LC[0]->GetXaxis()->FindBin(sidebandBeginMass);
    int beginW = hVBF_onshell_LC[0]->GetXaxis()->FindBin(beginWMass);
    int binWZthreshold = hVBF_onshell_LC[0]->GetXaxis()->FindBin(WZthreshold);
    int endZ = hVBF_onshell_LC[0]->GetXaxis()->FindBin(endZMass);
    int endSB = hVBF_onshell_LC[0]->GetXaxis()->FindBin(sidebandEndMass);
    int beginLowHH = hVBF_onshell_LC[0]->GetXaxis()->FindBin(beginLowHHMass);
    int beginHH = hVBF_onshell_LC[0]->GetXaxis()->FindBin(beginHHMass);
    int endHH = hVBF_onshell_LC[0]->GetXaxis()->FindBin(endHHMass);
    int sidebandBeginHH = hVBF_onshell_LC[0]->GetXaxis()->FindBin(sidebandBeginHHMass);
    int sidebandEndHH = hVBF_onshell_LC[0]->GetXaxis()->FindBin(sidebandEndHHMass);

    double simulatedSB = hVBF_onshell_LC[0]->Integral(beginSB, beginW-1) + hVBF_onshell_LC[0]->Integral(endZ, endSB-1);
    cout << "Total signal SB before scaling: " << simulatedSB << endl;
    double simulatedW = hVBF_onshell_LC[0]->Integral(beginW, binWZthreshold-1);
    double simulatedZ = hVBF_onshell_LC[0]->Integral(binWZthreshold,endZ-1);
    double simulatedWbkg = simulatedSB*(WZthreshold-beginWMass)/(endZMass-beginWMass);
    double simulatedZbkg = simulatedSB*(endZMass-WZthreshold)/(endZMass-beginWMass);
    cout << "W Total / Bkg: " << simulatedW << " / " << simulatedWbkg << endl;
    cout << "Z Total / Bkg: " << simulatedZ << " / " << simulatedZbkg << endl;

    double brval = BR_Table[1][2];
    if (folder==2) brval = BR_Table[1][3];
    brval *= 1000.;
    double ratioW = (XSEC_Table_WH[EnergyIndex][1]*brval*luminosity[EnergyIndex]) / (simulatedW-simulatedWbkg);
    double ratioZ = (XSEC_Table_ZH[EnergyIndex][1]*brval*luminosity[EnergyIndex]) / (simulatedZ-simulatedZbkg);
    double scale_SB_VBF = (simulatedW + simulatedZ)/(ratioW*simulatedW + ratioZ*simulatedZ);
    double scale_SB_HH = 1.5;
    double nVBF_Sig_Simulated = 0;
    cout << "VBF SB scale: " << scale_SB_VBF << endl;
    cout << "WH scale: " << ratioW << endl;
    cout << "ZH scale: " << ratioZ << endl;

    for (int tr = 0; tr < 4; tr++){
      for (int binx=1; binx<=hVBF_onshell_LC[tr]->GetNbinsX(); binx++){
        double bincontent = hVBF_onshell_LC[tr]->GetBinContent(binx);
        if (binx < endZ && binx >= binWZthreshold) bincontent *= ratioZ;
        else if (binx < binWZthreshold && binx >= beginW) bincontent *= ratioW;
        else{
          if ((binx < endSB && binx >= endZ) || (binx < beginW && binx >= beginSB)) bincontent *= scale_SB_VBF;
          if (tr==0){
            double bincount = bincontent;
            if (binx<beginLowHH || (binx>=beginHH && binx<endHH)) bincount=0;
            if ((binx<beginHH && binx>=sidebandBeginHH) || (binx<sidebandEndHH && binx>=endHH)) bincount*=scale_SB_HH;
            nVBF_Sig_Simulated += bincount;
          }
        }
        hVBF_onshell_scaled[tr]->SetBinContent(binx, bincontent);
      }
      for (int binx=1; binx<=hVBF_offshell_LC[tr]->GetNbinsX(); binx++){
        double bincontent = hVBF_offshell_LC[tr]->GetBinContent(binx);
        if ((binx < endSB && binx >= endZ) || (binx < beginW && binx >= beginSB)) bincontent *= scale_SB_VBF;
        if (binx < endZ && binx >= binWZthreshold) bincontent *= ratioZ;
        if (binx < binWZthreshold && binx >= beginW) bincontent *= ratioW;
        hVBF_offshell_scaled[tr]->SetBinContent(binx, bincontent);
      }
    }
    double vbfscale = (XSEC_Table_VBF[EnergyIndex][1]*brval*luminosity[EnergyIndex]) / nVBF_Sig_Simulated;
    cout << "VBF scale: " << vbfscale << endl;

    for (int tr = 0; tr < 4; tr++){
      for (int binx=1; binx<=hVBF_onshell_scaled[tr]->GetNbinsX(); binx++){
        double bincontent = hVBF_onshell_scaled[tr]->GetBinContent(binx);
        if (!(binx < endZ && binx >= binWZthreshold) && !(binx < binWZthreshold && binx >= beginW)) bincontent *= vbfscale;
        if (binx<beginLowHH || (binx>=beginHH && binx<endHH)) bincontent=0;
        if ((binx<beginHH && binx>=sidebandBeginHH) || (binx<sidebandEndHH && binx>=endHH)) bincontent*=scale_SB_HH;
        hVBF_onshell_scaled_wVBF[tr]->SetBinContent(binx, bincontent);
      }
      for (int binx=1; binx<=hVBF_offshell_scaled[tr]->GetNbinsX(); binx++){
        double bincontent = hVBF_offshell_scaled[tr]->GetBinContent(binx);
        if (!(binx < endZ && binx >= binWZthreshold) && !(binx < binWZthreshold && binx >= beginW)) bincontent *= vbfscale;
        if (binx<beginLowHH || (binx>=beginHH && binx<endHH)) bincontent=0;
        if ((binx<beginHH && binx>=sidebandBeginHH) || (binx<sidebandEndHH && binx>=endHH)) bincontent*=scale_SB_HH;
        hVBF_offshell_scaled_wVBF[tr]->SetBinContent(binx, bincontent);
      }
      cout << "Initial yield for " << templatetitles[tr] << ": " << hVBF_offshell_LC[tr]->Integral() << endl;
      cout << "Intermediate yield for " << templatetitles[tr] << " without VBF re-scaling: " << hVBF_offshell_scaled[tr]->Integral() << endl;
      cout << "Final yield for " << templatetitles[tr] << ": " << hVBF_offshell_scaled_wVBF[tr]->Integral() << endl;
      cout << endl;
    }

    TH1F* hRatio = (TH1F*)hVBF_offshell_scaled_wVBF[0]->Clone("Phantom_VBFVH_ScalingRatio");
    hRatio->Divide(hVBF_offshell_LC[0]);
    for (int binx=1; binx<=hRatio->GetNbinsX(); binx++){
      double bincontent = hRatio->GetBinContent(binx);
      if (bincontent<=0 && !(binx<beginLowHH || (binx>=beginHH && binx<endHH))) cout << "Bin " << binx << " is unexpectedly non-positive!" << endl;
    }
    foutput->WriteTObject(hRatio);
    delete hRatio;

    for (int tr = 0; tr < 4; tr++){
      foutput->WriteTObject(h2DVBF[tr]);
      foutput->WriteTObject(hVBF_onshell[tr]);
      foutput->WriteTObject(hVBF_offshell[tr]);
      foutput->WriteTObject(hVBF_onshell_LC[tr]);
      foutput->WriteTObject(hVBF_offshell_LC[tr]);
      foutput->WriteTObject(hVBF_onshell_scaled[tr]);
      foutput->WriteTObject(hVBF_offshell_scaled[tr]);
      foutput->WriteTObject(hVBF_onshell_scaled_wVBF[tr]);
      foutput->WriteTObject(hVBF_offshell_scaled_wVBF[tr]);

      delete h2DVBF[tr];
      delete hVBF_onshell[tr];
      delete hVBF_offshell[tr];
      delete hVBF_onshell_LC[tr];
      delete hVBF_offshell_LC[tr];
      delete hVBF_onshell_scaled[tr];
      delete hVBF_offshell_scaled[tr];
      delete hVBF_onshell_scaled_wVBF[tr];
      delete hVBF_offshell_scaled_wVBF[tr];
    }
    
    foutput->Close();
	}
}

void plotVBFVHReweighted(){
  gROOT->ProcessLine(".x tdrstyle.cc");
  gStyle->SetOptStat(0);

  TString INPUT_NAME = "HtoZZ4l_Phantom_125p6_VHDistributions_GenLevel.root";
  TString OUTPUT_NAME = "HtoZZ4l_Phantom_125p6_VHReweightingPlots_GenLevel.root";

  double mPOLE = 125.6;
  float ZZMass_PeakCut[2] ={ 125.1, 126.1 }; // Spin 0 analysis

  float templateWeight = 1;
  float GenDiJetMass;
  float MC_weight;
  float ZZMass = 0;

  TString coutput_common = user_TemplateswithTrees_dir + "../VHContributions/Plots/GenLevel/";
  gSystem->Exec("mkdir -p " + coutput_common);
  TString coutput = coutput_common + OUTPUT_NAME;
  TFile* foutput = new TFile(coutput, "recreate");

  const int kNumTemplates = 3;
  TString templatetitles[kNumTemplates] ={ "VBF Sig", "VBF Bkg", "VBF Int" };
  TString templatenames[kNumTemplates] ={ "VBF_Sig", "VBF_Bkg", "VBF_Int" };
  TH1F* hVBF_onshell_LC[kNumTemplates] ={ 0 };
  TH1F* hVBF_offshell_LC[kNumTemplates] ={ 0 };
  TH1F* hVBF_onshell_scaled[kNumTemplates] ={ 0 };
  TH1F* hVBF_offshell_scaled[kNumTemplates] ={ 0 };
  TH1F* hVBF_onshell_scaled_wVBF[kNumTemplates] ={ 0 };
  TH1F* hVBF_offshell_scaled_wVBF[kNumTemplates] ={ 0 };

  for (int erg_tev=7; erg_tev<=8; erg_tev++){
    for (int folder=0; folder<3; folder++){
      TString comstring;
      comstring.Form("%i", erg_tev);
      TString erg_dir;
      erg_dir.Form("LHC_%iTeV/", erg_tev);
      int EnergyIndex = 1;
      if (erg_tev == 7) EnergyIndex = 0;
      TString cinput_common = user_TemplateswithTrees_dir + "../VHContributions/" + erg_dir;
      cinput_common += user_folder[folder] + "/";
      TString cinput = cinput_common + INPUT_NAME;
      TFile* finput = new TFile(cinput, "read");
      if (finput->IsZombie()){
        delete finput;
        continue;
      }
      else if (finput==0) continue;

      cout << "Opened file " << finput->GetName() << endl;
      for (int tr=0; tr<kNumTemplates; tr++){
        TString templatename = "VBF_";
        templatename += templatenames[tr];

        TH1F* hVBF_onshell_LC_temp = (TH1F*)finput->Get(Form("%s_onshell_LC", templatename.Data()));
        TH1F* hVBF_offshell_LC_temp = (TH1F*)finput->Get(Form("%s_offshell_LC", templatename.Data()));
        TH1F* hVBF_onshell_scaled_temp = (TH1F*)finput->Get(Form("%s_onshell_scaled", templatename.Data()));
        TH1F* hVBF_offshell_scaled_temp = (TH1F*)finput->Get(Form("%s_offshell_scaled", templatename.Data()));
        TH1F* hVBF_onshell_scaled_wVBF_temp = (TH1F*)finput->Get(Form("%s_onshell_scaled_wVBF", templatename.Data()));
        TH1F* hVBF_offshell_scaled_wVBF_temp = (TH1F*)finput->Get(Form("%s_offshell_scaled_wVBF", templatename.Data()));

        foutput->cd();
        gStyle->SetOptStat(0);
        if (hVBF_onshell_LC[tr]==0) hVBF_onshell_LC[tr] = (TH1F*)hVBF_onshell_LC_temp->Clone(Form("%s_clone", hVBF_onshell_LC_temp->GetName()));
        else hVBF_onshell_LC[tr]->Add(hVBF_onshell_LC_temp);
        delete hVBF_onshell_LC_temp;

        if (hVBF_onshell_scaled[tr]==0) hVBF_onshell_scaled[tr] = (TH1F*)hVBF_onshell_scaled_temp->Clone(Form("%s_clone", hVBF_onshell_scaled_temp->GetName()));
        else hVBF_onshell_scaled[tr]->Add(hVBF_onshell_scaled_temp);
        delete hVBF_onshell_scaled_temp;

        if (hVBF_onshell_scaled_wVBF[tr]==0) hVBF_onshell_scaled_wVBF[tr] = (TH1F*)hVBF_onshell_scaled_wVBF_temp->Clone(Form("%s_clone", hVBF_onshell_scaled_wVBF_temp->GetName()));
        else hVBF_onshell_scaled_wVBF[tr]->Add(hVBF_onshell_scaled_wVBF_temp);
        delete hVBF_onshell_scaled_wVBF_temp;

        if (hVBF_offshell_LC[tr]==0) hVBF_offshell_LC[tr] = (TH1F*)hVBF_offshell_LC_temp->Clone(Form("%s_clone", hVBF_offshell_LC_temp->GetName()));
        else hVBF_offshell_LC[tr]->Add(hVBF_offshell_LC_temp);
        delete hVBF_offshell_LC_temp;

        if (hVBF_offshell_scaled[tr]==0) hVBF_offshell_scaled[tr] = (TH1F*)hVBF_offshell_scaled_temp->Clone(Form("%s_clone", hVBF_offshell_scaled_temp->GetName()));
        else hVBF_offshell_scaled[tr]->Add(hVBF_offshell_scaled_temp);
        delete hVBF_offshell_scaled_temp;

        if (hVBF_offshell_scaled_wVBF[tr]==0) hVBF_offshell_scaled_wVBF[tr] = (TH1F*)hVBF_offshell_scaled_wVBF_temp->Clone(Form("%s_clone", hVBF_offshell_scaled_wVBF_temp->GetName()));
        else hVBF_offshell_scaled_wVBF[tr]->Add(hVBF_offshell_scaled_wVBF_temp);
        delete hVBF_offshell_scaled_wVBF_temp;
      }
      finput->Close();
    }
  }

  TH1F* hVBF_sig[2][3]={
    { hVBF_onshell_LC[0], hVBF_onshell_scaled[0], hVBF_onshell_scaled_wVBF[0] },
    { hVBF_offshell_LC[0], hVBF_offshell_scaled[0], hVBF_offshell_scaled_wVBF[0] }
  };
  TH1F* hVBF_bkg[2][3]={
    { hVBF_onshell_LC[1], hVBF_onshell_scaled[1], hVBF_onshell_scaled_wVBF[1] },
    { hVBF_offshell_LC[1], hVBF_offshell_scaled[1], hVBF_offshell_scaled_wVBF[1] }
  };
  TH1F* hVBF_int[2][3]={
    { hVBF_onshell_LC[2], hVBF_onshell_scaled[2], hVBF_onshell_scaled_wVBF[2] },
    { hVBF_offshell_LC[2], hVBF_offshell_scaled[2], hVBF_offshell_scaled_wVBF[2] }
  };
  double max_plot[2][3] ={ { 0 } };
  double min_plot[2][3] ={ { 0 } };
//  TString strmzztitle[2]={ "105.6<m_{4l}<140.6 GeV", "220<m_{4l}<1600 GeV" };
  TString strmzztitle[2]={ "On-shell", "Off-shell" };
  TString strmzzname[2]={ "Onshell", "Offshell" };
  TString strBSItitle[3]={ "Signal", "Background", "Interference" };
  TString strBSIname[3]={ "Signal", "Background", "Interference" };
  TString strScaleSchemeTitle[3]={ " (Default Phantom)", " (VH Rescaling)", " (+VBF Rescaling)" };
  cout << "Set up canvas gadgets" << endl;
  for (int os=0; os<2; os++){
    cout << strmzzname[os] << endl;
    for (int sc=0; sc<3; sc++){
      cout << strBSItitle[0] << strScaleSchemeTitle[sc] << " mJJ>=130 GeV / mJJ<130 GeV: " << hVBF_sig[os][sc]->GetBinContent(hVBF_sig[os][sc]->GetNbinsX()) << " / " << (hVBF_sig[os][sc]->Integral() - hVBF_sig[os][sc]->GetBinContent(hVBF_sig[os][sc]->GetNbinsX())) << " = " << hVBF_sig[os][sc]->GetBinContent(hVBF_sig[os][sc]->GetNbinsX())/(hVBF_sig[os][sc]->Integral() - hVBF_sig[os][sc]->GetBinContent(hVBF_sig[os][sc]->GetNbinsX())) << endl;
      cout << strBSItitle[1] << strScaleSchemeTitle[sc] << " mJJ>=130 GeV / mJJ<130 GeV: " << hVBF_bkg[os][sc]->GetBinContent(hVBF_bkg[os][sc]->GetNbinsX()) << " / " << (hVBF_bkg[os][sc]->Integral() - hVBF_bkg[os][sc]->GetBinContent(hVBF_bkg[os][sc]->GetNbinsX())) << " = " << hVBF_bkg[os][sc]->GetBinContent(hVBF_bkg[os][sc]->GetNbinsX())/(hVBF_bkg[os][sc]->Integral() - hVBF_bkg[os][sc]->GetBinContent(hVBF_bkg[os][sc]->GetNbinsX())) << endl;
      cout << strBSItitle[2] << strScaleSchemeTitle[sc] << " mJJ>=130 GeV / mJJ<130 GeV: " << hVBF_int[os][sc]->GetBinContent(hVBF_int[os][sc]->GetNbinsX()) << " / " << (hVBF_int[os][sc]->Integral() - hVBF_int[os][sc]->GetBinContent(hVBF_int[os][sc]->GetNbinsX())) << " = " << hVBF_int[os][sc]->GetBinContent(hVBF_int[os][sc]->GetNbinsX())/(hVBF_int[os][sc]->Integral() - hVBF_int[os][sc]->GetBinContent(hVBF_int[os][sc]->GetNbinsX())) << endl;

      hVBF_sig[os][sc]->SetTitle("");
      hVBF_bkg[os][sc]->SetTitle("");
      hVBF_int[os][sc]->SetTitle("");
      hVBF_sig[os][sc]->GetXaxis()->SetRangeUser(22, 129.9);
      hVBF_bkg[os][sc]->GetXaxis()->SetRangeUser(22, 129.9);
      hVBF_int[os][sc]->GetXaxis()->SetRangeUser(22, 129.9);
      hVBF_sig[os][sc]->GetXaxis()->SetTitle("m^{true}_{jj} (GeV)");
      hVBF_bkg[os][sc]->GetXaxis()->SetTitle("m^{true}_{jj} (GeV)");
      hVBF_int[os][sc]->GetXaxis()->SetTitle("m^{true}_{jj} (GeV)");
      double binwidth = hVBF_sig[os][sc]->GetBinWidth(1);
      hVBF_sig[os][sc]->GetYaxis()->SetTitleOffset(1.5);
      hVBF_bkg[os][sc]->GetYaxis()->SetTitleOffset(1.5);
      hVBF_int[os][sc]->GetYaxis()->SetTitleOffset(1.5);
      hVBF_sig[os][sc]->GetYaxis()->SetTitle(Form("Events / %.0f GeV", binwidth));
      hVBF_bkg[os][sc]->GetYaxis()->SetTitle(Form("Events / %.0f GeV", binwidth));
      hVBF_int[os][sc]->GetYaxis()->SetTitle(Form("Events / %.0f GeV", binwidth));
      hVBF_sig[os][sc]->SetLineWidth(2);
      hVBF_bkg[os][sc]->SetLineWidth(2);
      hVBF_int[os][sc]->SetLineWidth(2);
      hVBF_sig[os][sc]->SetLineStyle(1);
      hVBF_bkg[os][sc]->SetLineStyle(1);
      hVBF_int[os][sc]->SetLineStyle(1);
      if (sc==0){
        hVBF_sig[os][sc]->SetLineColor(kRed);
        hVBF_bkg[os][sc]->SetLineColor(kRed);
        hVBF_int[os][sc]->SetLineColor(kRed);
      }
      if (sc==1){
        hVBF_sig[os][sc]->SetLineColor(kViolet);
        hVBF_bkg[os][sc]->SetLineColor(kViolet);
        hVBF_int[os][sc]->SetLineColor(kViolet);
      }
      else if (sc==2){
        hVBF_sig[os][sc]->SetLineColor(kBlue);
        hVBF_bkg[os][sc]->SetLineColor(kBlue);
        hVBF_int[os][sc]->SetLineColor(kBlue);
      }
      for (int bin=1; bin<hVBF_sig[os][sc]->GetNbinsX(); bin++){
        double bcsig = hVBF_sig[os][sc]->GetBinContent(bin);
        double bcbkg = hVBF_bkg[os][sc]->GetBinContent(bin);
        double bcint = hVBF_int[os][sc]->GetBinContent(bin);

        max_plot[os][0] = max(max_plot[os][0], bcsig);
        min_plot[os][0] = min(min_plot[os][0], bcsig);
        max_plot[os][1] = max(max_plot[os][1], bcbkg);
        min_plot[os][1] = min(min_plot[os][1], bcbkg);
        max_plot[os][2] = max(max_plot[os][2], bcint);
        min_plot[os][2] = min(min_plot[os][2], bcint);
      }
      cout << "Set up region " << os << " scheme " << sc << " complete." << endl;
    }
    for (int sc=0; sc<3; sc++){
      hVBF_sig[os][sc]->GetYaxis()->SetRangeUser(min_plot[os][0]*1.5, max_plot[os][0]*1.5);
      hVBF_bkg[os][sc]->GetYaxis()->SetRangeUser(min_plot[os][1]*1.5, max_plot[os][1]*1.5);
      hVBF_int[os][sc]->GetYaxis()->SetRangeUser(min_plot[os][2]*1.5, max_plot[os][2]*1.5);
    }
  }

  foutput->cd();

  for (int os=0; os<2; os++){
    for (int bsi=0; bsi<3; bsi++){
      cout << "Begin plot of region " << os << endl;

      TPaveText* pt = new TPaveText(0.15, 0.93, 0.85, 1, "brNDC");
      pt->SetBorderSize(0);
      pt->SetFillStyle(0);
      pt->SetTextAlign(12);
      pt->SetTextFont(42);
      pt->SetTextSize(0.045);
      TText* text = pt->AddText(0.025, 0.45, "#font[61]{CMS}");
      text->SetTextSize(0.044);
      text = pt->AddText(0.165, 0.42, "#font[52]{Simulation}");
      text->SetTextSize(0.0315);
//      TString cErgTev = "#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}";
      TString cErgTev = "#font[42]{              2e+2#mu 19.7 fb^{-1} (8 TeV)}";
      text = pt->AddText(0.537, 0.45, cErgTev);
      text->SetTextSize(0.0315);

      TString appendName;
      appendName = "_";
      appendName += strmzzname[os];
      appendName += strBSIname[bsi];
      TString canvasname = "cCompareMJJ_GenLevel";
      canvasname.Append(appendName);
      TCanvas* cc = new TCanvas(canvasname, "", 8, 30, 800, 800);
      gStyle->SetOptStat(0);
      cc->cd();
      gStyle->SetOptStat(0);
      cc->SetFillColor(0);
      cc->SetBorderMode(0);
      cc->SetBorderSize(2);
      cc->SetTickx(1);
      cc->SetTicky(1);
      cc->SetLeftMargin(0.17);
      cc->SetRightMargin(0.05);
      cc->SetTopMargin(0.07);
      cc->SetBottomMargin(0.13);
      cc->SetFrameFillStyle(0);
      cc->SetFrameBorderMode(0);
      cc->SetFrameFillStyle(0);
      cc->SetFrameBorderMode(0);

      TLegend *ll;
      ll = new TLegend(0.22, 0.70, 0.60, 0.90);
      ll->SetBorderSize(0);
      ll->SetTextFont(42);
      ll->SetTextSize(0.03);
      ll->SetLineColor(1);
      ll->SetLineStyle(1);
      ll->SetLineWidth(1);
      ll->SetFillColor(0);
      ll->SetFillStyle(0);

      for (int sc=0; sc<3; sc++){
        if (bsi==0){
          TString legendLabel = strBSItitle[bsi] + strScaleSchemeTitle[sc];
          ll->AddEntry(hVBF_sig[os][sc], legendLabel, "l");
          if (sc==0) hVBF_sig[os][sc]->Draw("hist");
          else hVBF_sig[os][sc]->Draw("histsame");
        }
        else if (bsi==1){
          TString legendLabel = strBSItitle[bsi] + strScaleSchemeTitle[sc];
          ll->AddEntry(hVBF_bkg[os][sc], legendLabel, "l");
          if (sc==0) hVBF_bkg[os][sc]->Draw("hist");
          else hVBF_bkg[os][sc]->Draw("histsame");
        }
        else if (bsi==2){
          TString legendLabel = strBSItitle[bsi] + strScaleSchemeTitle[sc];
          ll->AddEntry(hVBF_int[os][sc], legendLabel, "l");
          if (sc==0) hVBF_int[os][sc]->Draw("hist");
          else hVBF_int[os][sc]->Draw("histsame");
        }
      }

      ll->Draw("same");
      pt->Draw();

      TPaveText *pt10 = new TPaveText(0.80, 0.86, 0.90, 0.90, "brNDC");
      pt10->SetBorderSize(0);
      pt10->SetTextAlign(12);
      pt10->SetTextSize(0.03);
      pt10->SetFillStyle(0);
      pt10->SetTextFont(42);
      TText* text10;
      text10 = pt10->AddText(0.01, 0.01, strmzztitle[os]);
      pt10->Draw();

      foutput->WriteTObject(cc);

      delete pt10;
      delete ll;
      cc->Close();
      delete pt;

      cout << "End plot of region " << os << endl;
    }
  }

  for (int tr = 0; tr < kNumTemplates; tr++){
    delete hVBF_onshell_LC[tr];
    delete hVBF_offshell_LC[tr];
    delete hVBF_onshell_scaled[tr];
    delete hVBF_offshell_scaled[tr];
    delete hVBF_onshell_scaled_wVBF[tr];
    delete hVBF_offshell_scaled_wVBF[tr];
  }
  foutput->Close();
}

double doLinearCombination(double first, double c_first, double second, double c_second, double bkg, int outputtype){
	double aa[3] = { first, c_first, sqrt(c_first) };
	double bb[3] = { second, c_second, sqrt(c_second) };

	double sig = (aa[0] * bb[2] - aa[2] * bb[0] + (aa[2] - bb[2])*bkg) / (aa[1] * bb[2] - aa[2] * bb[1]);
	double interf = (aa[0] * bb[1] - aa[1] * bb[0] + (aa[1] - bb[1])*bkg) / (aa[2] * bb[1] - aa[1] * bb[2]);
//	if (sig < 0){ sig = 0; interf = 0; };

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

//	cout << "VBF peak as calculated\n";
//	cout << "Signal: " << sig << "\tInterf: " << interf << endl;
	
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
