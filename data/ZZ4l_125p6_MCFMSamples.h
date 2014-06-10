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
	kSig,
	kBkg,
	kBSI,

	kNumSamples
};
string sampleName[kNumSamples] = {
	"Sig",
	"Bkg",
	"BSI"
};
char* sample_suffix[kNumSamples]={
	"Sig",
	"Bkg",
	"BSI"
};

float luminosity[2]={5.1,19.712};
float xsec_ggHZZ_BSI_MCFM[2] = {0.009007274418,0.011840779163}; // 4e+4mu+2e2mu, 2e2mu=4e+4mu=2x4e
float xsec_ggHZZ_MCFM[2] = {0.002549983848,0.003255176587}; // 4e+4mu+2e2mu, 2e2mu=4e+4mu=2x4e
float xsec_ggZZ_MCFM[2] = {0.006885033938,0.009188219964}; // 4e+4mu+2e2mu, 2e2mu=4e+4mu=2x4e
float xsec_MCFM[2][3]={
	{	0.002549983848,	0.006885033938,	0.009007274418	},
	{	0.003255176587,	0.009188219964,	0.011840779163	}
};
float pTLE10Scale = 1.00105261;

char hzz4lprefix[]="HZZ4lTree_";
string user_dir="/afs/cern.ch/work/u/usarica/WidthAnalysis/MCFM_ggF/";
char* user_folder[3]={
	"4mu",
	"4e",
	"2mu2e"
};
