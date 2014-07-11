HiggsWidth_PostICHEP
====================

To produce templates for https://github.com/HZZ4l/CreateWidthDatacards, in your CMSSW release, follow these directions:

_NOTE_: Use **legacywithfixes** tag to reproduce templates for width results.

0) Check that data/ZZ4l_125_6_Samples.h is pointing to the correct directories.

1) root -q -b makeCombineTemplates_Modified_MCFM.c+

This will get the initial D_Gamma_gg_r10 templates. VBF templates are reweighted from ggH. The reweighting functions for 7 and 8 TeV. It will run over all 4l final states (2e2mu, 4e, 4mu), all systematics (Nominal, PDF up/down, QCD up/down), and with or without k3a smoothing.

2) root -q -b makeCombineTemplateswithTrees_Modified_MCFM_GenLevelVBF.c+

Similar to 1, but will add trees with templates. Currently has an option to use resolution smeared samples for VBF, this will be removed once Full-Sim is implemented. These files will take up a fair amount of space.

3) Check out https://github.com/jbsauvan/TemplateBuilder. _NOTE_: To reproduce legacy width studies, use tag **070314_add-fa3-configs**.

To produce correct JSON files for TemplateBuilder, run

	./pushjson.sh <Location of LHC_*TeV directories that contain templates with trees> <Location of TemplateBuilder/run directory>

Then, compile TemplateBuilder and submit jobs to queue to smooth each template.

4) root -q -b makeCombineTemplatesSmooth_Modified_MCFM_GenLevelVBF.c+ 

This will produce the final templates needed for combine. Normalizations will be output as the templates are made, you can also compare different templates using included comparesmoothedtemplates.C script.
