# No imports needed


def processDictionary():
   processDict = {
      "GG" : "ProcessHandler::kGG",
      "VBF" : "ProcessHandler::kVBF",
      "ZH" : "ProcessHandler::kZH",
      "WH" : "ProcessHandler::kWH",
      "TT" : "ProcessHandler::kTT",
      "BB" : "ProcessHandler::kBB",
      "QQBkg" : "ProcessHandler::kQQBkg",
      "ZX" : "ProcessHandler::kZX"
   }
   return processDict


def getAnalysisRegions():
   anaregions = [ "NMassRegions", "kOnshell", "kOffshell" ]
   return anaregions


def getChannelList():
   channels = [ "NChannels", "k2e2mu", "k4e", "k4mu" ]
   return channels


def getCategoryList():
   categories = [ "Inclusive", "HadVHTagged", "JJVBFTagged", "Untagged" ]
   return categories


def getFRMethodList():
   frmethods = [ "ZXFakeRateHandler::NFakeRateMethods", "ZXFakeRateHandler::mSS" ]
   return frmethods


def getSystematicsList():
   systematics = [
      "sNominal",
      "tPDFScaleDn", "tPDFScaleUp",
      "tQCDScaleDn", "tQCDScaleUp",
      "tAsMZDn", "tAsMZUp",
      "tPDFReplicaDn", "tPDFReplicaUp",
      "tPythiaScaleDn", "tPythiaScaleUp",
      "tPythiaTuneDn", "tPythiaTuneUp",
      "tMINLODn", "tMINLOUp",
      "tQQBkgEWCorrDn", "tQQBkgEWCorrUp",
      "eLepSFEleDn", "eLepSFEleUp",
      "eLepSFMuDn", "eLepSFMuUp",
      "eLepScaleEleDn", "eLepScaleEleUp",
      "eLepScaleMuDn", "eLepScaleMuUp",
      "eLepResEleDn", "eLepResEleUp",
      "eLepResMuDn", "eLepResMuUp",
      "eJECDn", "eJECUp",
      "eBTagSFDn", "eBTagSFUp",
      "eZXStatsDn", "eZXStatsUp"
   ]
   return systematics


def getACHypothesisList():
   hypos = [ "nACHypotheses", "kSM", "kL1", "kA2", "kA3", "kL1ZGs" ]
   return hypos


def checkValidRun(syst, cat, ch, proc, generator=""):
   testval=True
   if cat == "Inclusive" and ("eJEC" in syst or "eBTag" in syst or "tMINLO" in syst or "tPythia" in syst):
      testval=False
   if not(proc=="TT" or proc=="BB") and "eBTag" in syst:
      testval=False
   if (ch == "k4e" and ("eLepSFMu" in syst or "eLepScaleMu" in syst or "eLepResMu" in syst)) or (ch == "k4mu" and ("eLepSFEle" in syst or "eLepScaleEle" in syst or "eLepResEle" in syst)):
      testval=False
   if generator == "MCFM" and ("tMINLO" in syst or "tPythia" in syst or "LepScale" in syst or "LepRes" in syst):
      testval=False
   if "tMINLO" in syst and proc != "GG":
      testval=False
   if proc == "BB" and ("tMINLO" in syst or "tPythia" in syst):
      testval=False
   if (proc == "QQBkg" or proc == "ZX") and ("LepScale" in syst or "LepRes" in syst):
      testval=False
   if proc == "QQBkg" and ("tMINLO" in syst or "tPythia" in syst):
      testval=False
   if proc == "ZX" and not(syst=="sNominal" or "ZX" in syst):
      testval=False
   if proc != "ZX" and "ZX" in syst:
      testval=False
   if proc != "QQBkg" and "QQBkg" in syst:
      testval=False
   return testval

