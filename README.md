# HiggsWidth_PostICHEP


## KD constants
To create pairwise KD constants, run getKDConstant_D* or write your own. Once the output ROOT file is present, pass a spline trough the TGraph by using SmoothKDConstantProducer_D*.


## Template production
To produce templates for https://github.com/HZZ4l/CreateWidthDatacards, in your CMSSW release, follow these directions:

_NOTE_: Use **legacywithfixes** tag to reproduce templates for width results.

0) Check that data/ZZ4l_125_6_Samples.h is pointing to the correct directories.

1) root -q -b -l getLeptonicVHWeights.c+

2) root -q -b -l makeCombineTemplates_Modified_MCFM.c+

This will get the initial D_Gamma_gg_r10 templates. VBF templates are reweighted from ggH. The reweighting functions for 7 and 8 TeV. It will run over all 4l final states (2e2mu, 4e, 4mu), all systematics (Nominal, PDF up/down, QCD up/down), and with or without k3a smoothing.

3)

a) root -q -b -l "makeCombineTemplateswithTrees_Modified_MCFM_GenLevelVBF.c++(0)"
b) root -q -b -l "makeCombineTemplateswithTrees_Modified_MCFM_GenLevelVBF.c++(-1)";root -q -b -l "makeCombineTemplateswithTrees_Modified_MCFM_GenLevelVBF.c++(1)"

Similar to 2, but will add trees with templates. The legacy version has an option to use resolution smeared samples for VBF, this was made obsolete with full-sim Phantom samples. These files will take up a fair amount of space. Option a will produce templates without any Djet splitting. Option b will produce templates with Djet splitting. _NOTE_: All other code will produce all three as is, these take a considerable amount of time to run so this permits running in parallel.

4) Check out https://github.com/jbsauvan/TemplateBuilder. _NOTE_: To reproduce legacy width studies, use tag **070314_add-fa3-configs**.

To produce correct JSON files for TemplateBuilder, run

	./pushjson.sh <Location of LHC_*TeV directories that contain templates with trees> <Location of TemplateBuilder/run directory>

Then, compile TemplateBuilder and submit jobs to queue to smooth each template.

5) root -q -b makeCombineTemplatesSmooth_Modified_MCFM_GenLevelVBF.c+ 

This will produce the final templates needed for combine. Normalizations will be output as the templates are made, you can also compare different templates using included comparesmoothedtemplates.C script.
