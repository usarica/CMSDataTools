#include <iostream>
#include <string>
#include "TMath.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TCanvas.h"

const float PI_VAL = TMath::Pi();

enum sample {
	kSig_VBF,
	kBkg_VBF,
	kBSI_VBF,
	kBSI10_VBF,
	kBSI25_VBF,

	kNumSamples
};
string sampleName[kNumSamples] = {
	"Sig_VBF",
	"Bkg_VBF_Phantom",
	"BSI_VBF_Phantom",
	"BSI10_VBF_Phantom",
	"BSI25_VBF_Phantom"
};
char* sample_suffix[kNumSamples]={
	"Sig_VBF",
	"Bkg_VBF_Phantom",
	"BSI_VBF_Phantom",
	"BSI10_VBF_Phantom",
	"BSI25_VBF_Phantom"
};

float luminosity[2]={5.1,19.712};
float xsec_ggHZZ_CERN[2]={	0.0020659600, 0.002562639 };
float xsec_qqZZ_CERN[2][2] = {
	{0.06609,0.152},
	{0.07691,0.1767}
};
float xsec_ggHZZ_BSI_MCFM[2] = {0.009007274418,0.011840779163}; // 4e+4mu+2e2mu, 2e2mu=4e+4mu=2x4e
float xsec_ggHZZ_MCFM[2] = {0.002549983848,0.003255176587}; // 4e+4mu+2e2mu, 2e2mu=4e+4mu=2x4e
float xsec_ggZZ_MCFM[2] = {0.006885033938,0.009188219964}; // 4e+4mu+2e2mu, 2e2mu=4e+4mu=2x4e
float xsec_ggHZZ_gg2VV[2][2] = {
	{0.00046,0.00084},
	{0.00059,0.00107}
}; // 4e,4mu;2e2mu, 2e2mu=4e+4mu=2x4e
float xsec_ggZZ_gg2VV[2][2] = {
	{0.00129,0.00262},
	{0.00172,0.00349}
}; // 4e,4mu;2e2mu, 2e2mu=4e+4mu=2x4e
float xsec_ggHZZ_BSI_gg2VV[2][2] = {
	{0.00162,0.00325},
	{0.00215,0.00426}
}; // 4e,4mu;2e2mu, 2e2mu=4e+4mu=2x4e

float xsec_VBF[2][kNumSamples]={
	{	2.13912e-4,	2.4691731774265e-4,	4.33242746288933e-4,	0,	0	},
	{	2.73069e-4,	3.2951625929086e-4,	5.09057867811004e-4,	5.8008977494252e-4,	7.44289142386252070e-4	}
};
float xsec_VBF_New[2][3][kNumSamples] = {
	{
		{ 2.13912e-4, 1.611800052932084E-004, 2.053859553143833E-004, 2.187258636589030E-004, 2.505576749660674E-004 },
		{ 2.13912e-4, 1.610671663387286E-004, 2.050597663931831E-004, 2.178041943144518E-004, 2.502592444436491E-004 },
		{ 2.73069e-4, 3.718595894117134E-004, 4.512172284721055E-004, 4.794819871366363E-004, 5.417040330546889E-004 }
	},
	{
		{ 2.13912e-4, 1.952788591835105E-004, 2.487975936891337E-004, 2.678337815958082E-004, 3.117289755556746E-004 },
		{ 2.13912e-4, 1.952187837277965E-004, 2.484688144236787E-004, 2.680082365969152E-004, 3.117402751925978E-004 },
		{ 2.73069e-4, 4.447269876857436E-004, 5.434945696873235E-004, 5.814292186114383E-004, 6.707175001136875E-004 }
	}
};
float pTLE10Scale = 1.00105261;
float mJJcutAcceptance = 1.0/(1.0-0.0547638);

char hzz4lprefix[]="HZZ4lTree_";
string user_dir="/afs/cern.ch/work/u/usarica/WidthAnalysis/VBF/";
string user_dir_VBF="/afs/cern.ch/work/u/usarica/WidthAnalysis/VBF/";
char* user_folder[3]={
	"4mu",
	"4e",
	"2mu2e"
};
