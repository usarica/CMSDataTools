{
  gSystem->Load("$CMSSW_BASE/src/ZZMatrixElement/MELA/test/loadMELA.C");
  gSystem->Load("$CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/loadLib.C");

  gSystem->AddIncludePath("-I$CMSSW_BASE/src/HiggsWidth_PostICHEP/Analysis/interface/");
  gSystem->Load("libHiggsWidth_PostICHEPAnalysis.so");
}
