{
  gSystem->Load("$CMSSW_BASE/src/ZZMatrixElement/MELA/test/loadMELA.C");
  gSystem->Load("$CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/loadLib.C");

  gSystem->AddIncludePath("-I$CMSSW_BASE/src/CMSDataTools/AnalysisTree/interface/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/CMSDataTools/AnalysisTree/test/");
  gSystem->Load("libCMSDataToolsAnalysisTree.so");
}
