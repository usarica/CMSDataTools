HiggsWidth_PostICHEP
====================

Run

1) makeCombineTemplates_Modified_MCFM.c

to get VBF templates that are reweighted from ggH. The reweighting functions for 7 and 8 TeV are

	char cinput_VBF_Sig[1000];
	sprintf(cinput_VBF_Sig,"%s%i%s","./data/HZZ4l-125_6-",erg_tev,"TeV-Sig_MCFM_PhantomVBF_Comparison.root");


Changes needed to synchronize:

ggH yield normalizations should change to 105.6 (inclusive) - 140.6:

7 TeV		Yield
2e2mu		1.4452079
4e		0.6087736
4mu		1.0902689
8 TeV		Scaled
2e2mu		6.6562629
4e		2.6944639
4mu		5.1963998

VBF normalization (same m4l range):

Yield
1.2861921E-01
5.1755897E-02
9.2458836E-02

6.1781689E-01
2.4788553E-01
4.6798807E-01

Options are

(int folder, int erg_tev, int tFitD=0, int Systematics=0, bool isSmooth=false)

folder: 0,1,2 for 2mu2e, 4mu, 4e
erg_tev: 7 or 8
tFitD: Only "6" is supported
isSmooth: Run both to get a smooth ggH reweighted curve (with k3a smoothing; this is only for systematics treatment in the fourth stage)

2) makeCombineTemplateswithTrees_Modified_MCFM_GenLevelVBF.c

Options are te same except

bool useResoVBF=true

Run with isSmooth==true; the trees produced are naturally without smoothing, but the template builder will receive this file to refer to the trees.
Run with isSmooth==false; for the last step to provide all possible templates.

the input probably needs to be re-adjusted so that we use full sim samples for VBF

Changes for yield normalization to the peak are still needed



3) Check out Template Builder and put HiggsWidth into TemplateBuilder/run/. Execute buildTemplate.exe on json files in 7/8TeV folders

4) makeCombineTemplatesSmooth_Modified_MCFM_GenLevelVBF.c 

Run to get te final templates. TemplateBuilder conserves the sum of weights, so whatever was passed before should be preserved. Noticve, raw and smooth templates from the previous steps are also needed.




