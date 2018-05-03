import sys
import imp
import copy
import os
import shutil
import pickle
import math
import pprint
import subprocess




def processDictionary():
   processDict = {
      "GG" : "ProcessHandler::kGG",
      "VBF" : "ProcessHandler::kVBF",
      "ZH" : "ProcessHandler::kZH",
      "WH" : "ProcessHandler::kWH",
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
      "eLepSFEleDn", "eLepSFEleUp",
      "eLepSFMuDn", "eLepSFMuUp",
      "tPDFScaleDn", "tPDFScaleUp",
      "tQCDScaleDn", "tQCDScaleUp",
      "tAsMZDn", "tAsMZUp",
      "tPDFReplicaDn", "tPDFReplicaUp",
      "tPythiaScaleDn", "tPythiaScaleUp",
      "tPythiaTuneDn", "tPythiaTuneUp",
      "tMINLODn", "tMINLOUp",
      "tQQBkgEWCorrDn", "tQQBkgEWCorrUp",
      "eJECDn", "eJECUp",
      "eZXStatsDn", "eZXStatsUp"
   ]
   return systematics


def getACHypothesisList():
   hypos = [ "nACHypotheses", "kSM", "kL1", "kA2", "kA3" ]
   return hypos


def checkValidRun(syst, cat, ch, proc, generator=""):
   testval=True
   if cat == "Inclusive" and ("eJEC" in syst or "tMINLO" in syst or "tPythia" in syst):
      testval=False
   if (ch == "k4e" and "eLepSFMu" in syst) or (ch == "k4mu" and "eLepSFEle" in syst):
      testval=False
   if generator == "MCFM" and ("tMINLO" in syst or "tPythia" in syst):
      testval=False
   if "tMINLO" in syst and proc != "GG":
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
