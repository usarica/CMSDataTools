{
  gSystem->AddIncludePath("-I$ROOFITSYS/include/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/interface/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/HiggsWidth_PostICHEP/Analysis/interface/");
  gSystem->Load("$CMSSW_BASE/src/ZZMatrixElement/MELA/test/loadMELA.C");
  gSystem->Load("libHiggsWidth_PostICHEPAnalysis.so");
}
