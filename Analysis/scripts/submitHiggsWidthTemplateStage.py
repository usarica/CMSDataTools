#!/bin/env python

import sys
import imp
import copy
import os
import shutil
import pickle
import math
import pprint
import subprocess
from datetime import date
from optparse import OptionParser


class StageXBatchManager:
   def __init__(self):
      # define options and arguments ====================================
      self.parser = OptionParser()

      self.parser.add_option("--outdir", dest="outdir", type="string", help="Name of the local output directory for your jobs. This directory will be created automatically.", default="./")

      self.parser.add_option("--process", dest="process", type="string", help="Name of the process")
      self.parser.add_option("--generator", dest="generator", type="string", help="Name of the generator")
      self.parser.add_option("--stage", dest="stage", type="int", default=1, help="Stage 1, 2 (default=1)")
      self.parser.add_option("--fixedDate", dest="fixedDate", type="string", help="Fixed output directory", default="")

      self.parser.add_option("--syst", dest="customSysts", type="string", action="append", help="Systematics to run (default=all turned on)")
      self.parser.add_option("--channel", dest="customChannels", type="string", action="append", help="Channels to run (default=all turned on)")
      self.parser.add_option("--category", dest="customCategories", type="string", action="append", help="Categories to run (default=all turned on)")
      self.parser.add_option("--AChypo", dest="customACHypos", type="string", action="append", help="Anomalous couplings hypotheses to run (default=all turned on)")
      self.parser.add_option("--FRMethod", dest="customFRMethods", type="string", action="append", help="ZX fake rate methods (default=all turned on)")

      self.parser.add_option("--dry", dest="dryRun", action="store_true", default=False, help="Do not submit jobs, just set up the files")
      self.parser.add_option("--interactive", dest="interactive", action="store_true", default=False, help="Do not submit jobs; run them interactively")

      self.parser.add_option("--checkstage", dest="checkstage", action="store_true", default=False, help="Submit checkstage functions instead of stage functions themselves")
      self.parser.add_option("--plotcheckstage", dest="plotcheckstage", action="store_true", default=False, help="Plot checkstage")
      self.parser.add_option("--plotcheckstagesystpairs", dest="plotcheckstagesystpairs", action="store_true", default=False, help="Plot checkstage systematics ratiso to nominal")

      (self.opt,self.args) = self.parser.parse_args()

      if not hasattr(self.opt, "process") or self.opt.process is None:
         sys.exit("Need to set --process option")
      elif self.opt.process != "ZX" and (not hasattr(self.opt, "generator") or self.opt.generator is None):
         sys.exit("Need to set --generator option")

      if self.opt.process == "ZX":
         self.generator = "Data"
      else:
         self.generator = self.opt.generator

      strscript="make{}TemplatesFrom{}".format(self.opt.process, self.generator)
      self.scriptname="{}.cc".format(strscript)

      if not os.path.isfile(self.scriptname):
         sys.exit("Script {} does not exist. Exiting...".format(self.scriptname))

      if self.opt.plotcheckstage or self.opt.plotcheckstagesystpairs:
         if self.opt.plotcheckstage and self.opt.plotcheckstagesystpairs:
            sys.exit("Cannot specify both plotcheckstage and plotcheckstagesystpairs")
         self.opt.checkstage=True
         self.opt.interactive=True

      self.fcnname=""
      if self.opt.checkstage:
         if self.opt.plotcheckstage:
            self.fcnname="plotProcessCheckStage"
         elif self.opt.plotcheckstagesystpairs:
            self.fcnname="plotProcessCheckStage_SystPairs"
         else:
            self.fcnname="{}_checkstage".format(strscript)
      elif self.opt.stage==1:
         self.fcnname="{}_one".format(strscript)
      elif self.opt.stage==2:
         self.fcnname="{}_two".format(strscript)
      if not self.fcnname:
         sys.exit("The function name could not be generated. Exiting...")

      self.mkdir(self.opt.outdir)
      self.rm(self.opt.outdir + '/' + self.fcnname + ".c")
      self.rm(self.opt.outdir + '/' + self.fcnname + "_c.d")
      self.rm(self.opt.outdir + '/' + self.fcnname + "_c.so")
      self.rm(self.opt.outdir + '/' + self.fcnname + "_c_ACLiC_dict_rdict.pcm")
      self.rm(self.opt.outdir + '/' + self.fcnname + ".c")
      self.cp(self.scriptname, self.opt.outdir + '/' + self.fcnname + ".c")
      self.submitJobs()


   def getFcnArguments( self, fname, fcnname ):
      fcnargs=[]
      fcnfound=False
      fcnendfound=False
      with open(fname) as testfile:
         for line in testfile:
            if fcnfound and fcnendfound: break
            if fcnname in line:
               fcnfound=True
            if fcnfound:
               linecpy=line
               linecpy=linecpy.rstrip()
               linecpy=linecpy.replace(' ','')
               linecpy=linecpy.replace(fcnname,'')
               linecpy=linecpy.replace('(','')
               linecpy=linecpy.replace(')','')
               linecpy=linecpy.replace('{','')
               tmpargs=linecpy.split(',')
               for tmparg in tmpargs:
                  fcnargs.append(tmparg.lower())
            if fcnfound and ')' in line:
               fcnendfound=True
      if not (fcnfound and fcnendfound):
         sys.exit("Function {} is not found in file {}!".format(fcnname,fname))
      return fcnargs


   def mkdir( self, dirname ):
      strcmd = 'mkdir -p %s' % dirname
      ret = os.system(strcmd)
      if( ret != 0 ):
         print 'Please remove or rename directory: ', dirname
         sys.exit(4)


   def rm( self, target ):
      strcmd = "rm -rf {}".format(target)
      ret = os.system(strcmd)
      if( ret != 0 ):
         print 'Command {} failed!'.format(strcmd)
         sys.exit(4)


   def cp( self, inname, target ):
      strcmd = "cp {} {}".format(inname, target)
      ret = os.system(strcmd)
      if( ret != 0 ):
         print 'Command {} failed!'.format(strcmd)
         sys.exit(4)


   def submitJobs(self):
      channels = [ "NChannels", "k2e2mu", "k4e", "k4mu" ]
      categories = [ "Inclusive", "HadVHTagged", "JJVBFTagged", "Untagged" ]
      hypos = [ "nACHypotheses", "kSM", "kL1", "kA2", "kA3" ]
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
      frmethods = [ "ZXFakeRateHandler::NFakeRateMethods", "ZXFakeRateHandler::mSS" ]

      fcnargs=self.getFcnArguments(self.scriptname, self.fcnname)
      argstr=""
      for fcnarg in fcnargs:
         tmpargstr=""
         if "channel" in fcnarg:
            tmpargstr = "{channel}"
         elif "category" in fcnarg:
            tmpargstr = "{category}"
         elif "syst" in fcnarg:
            tmpargstr = "{systematic}"
         elif "hypo" in fcnarg:
            tmpargstr = "{achypothesis}"
         elif "frmethod" in fcnarg:
            tmpargstr = "{frmethod}"
         # For the rest of if-statements, do not set tmpargstr; append to argstr directly
         elif "istage" in fcnarg:
            strStage = str(self.opt.stage)
            if argstr:
               argstr = "{},{}".format(argstr, strStage)
            else:
               argstr=strStage
         elif "sqrts" in fcnarg:
            strSqrts = str(self.opt.sqrts)
            if argstr:
               argstr = "{},{}".format(argstr, strSqrts)
            else:
               argstr=strSqrts
         elif "fixeddate" in fcnarg:
            strfixedDate=""
            if self.opt.fixedDate:
               if not self.opt.interactive:
                  strfixedDate="\\\"{}\\\"".format(self.opt.fixedDate)
               else:
                  strfixedDate=r"\\\"{}\\\"".format(self.opt.fixedDate)
            else:
               if not self.opt.interactive:
                  strfixedDate="\\\"\\\""
               else:
                  strfixedDate=r"\\\"\\\""
            if argstr:
               argstr = "{},{}".format(argstr, strfixedDate)
            else:
               argstr=strfixedDate
         if tmpargstr:
            if argstr:
               argstr = "{},{}".format(argstr, tmpargstr)
            else:
               argstr=tmpargstr

      print "Argument string: ",argstr

      for ch in channels:
         if (not "channel" in argstr) and ch!="NChannels":
            break
         elif ("channel" in argstr):
            if ch=="NChannels":
               continue
            elif self.opt.customChannels is not None:
               if not ch in self.opt.customChannels:
                  continue

         for cat in categories:
            if (not "category" in argstr) and cat!="Inclusive":
               break
            elif ("category" in argstr):
               if self.opt.customCategories is not None:
                  if not cat in self.opt.customCategories:
                     continue

            for hypo in hypos:
               if (not "achypothesis" in argstr) and hypo!="nACHypotheses":
                  break
               elif ("achypothesis" in argstr):
                  if hypo=="nACHypotheses":
                     continue
                  elif self.opt.customACHypos is not None:
                     if not hypo in self.opt.customACHypos:
                        continue

               for syst in systematics:
                  if (not "systematic" in argstr) and syst!="sNominal":
                     break
                  elif ("systematic" in argstr):
                     if self.opt.customSysts is not None:
                        if not syst in self.opt.customSysts:
                           continue

                  for frm in frmethods:
                     if (not "frmethod" in argstr) and frm!="ZXFakeRateHandler::NFakeRateMethods":
                        break
                     elif ("frmethod" in argstr):
                        if frm=="ZXFakeRateHandler::NFakeRateMethods":
                           continue
                        elif self.opt.customFRMethods is not None:
                           if not frm in self.opt.customFRMethods:
                              continue

                     # Do not submit unnecessary jobs
                     if cat == "Inclusive" and ("eJEC" in syst or "tMINLO" in syst or "tPythia" in syst):
                        continue
                     if self.opt.stage == 1 and cat == "Untagged" and not(self.opt.process == "ZH" or self.opt.process == "WH"):
                        print "{} category distributions in process {} can be obtained from the distributions of inclusive and other categories.".format(cat, self.opt.process)
                        continue
                     if self.opt.process == "QQBkg" and ("tMINLO" in syst or "tPythia" in syst):
                        continue
                     if self.opt.process == "ZX" and not(syst=="sNominal" or "ZX" in syst):
                        continue
                     if self.opt.process != "ZX" and "ZX" in syst:
                        continue
                     if self.opt.process != "QQBkg" and "QQBkg" in syst:
                        continue
                     if "tMINLO" in syst or "tPythia" in syst:
                        print "{} systematic distributions in process {} are obtained in mass ratios step.".format(syst, self.opt.process)
                        continue

                     strscrcmd = argstr.format(channel=ch,category=cat,achypothesis=hypo,systematic=syst,frmethod=frm)
                     strscrcmd = strscrcmd.replace(' ','') # The command passed to bash script should not contain whitespace itself
                     jobcmd = "submitHiggsWidthTemplateStageGeneric.sh {} \({}\)".format(self.fcnname, strscrcmd)
                     if self.opt.interactive:
                        jobcmd = "root -l -b -q -e \"gROOT->ProcessLine(\\\".x loadLib.C\\\");gROOT->ProcessLine(\\\".x {}.c+({})\\\");\"".format(self.fcnname, strscrcmd)
                     if self.opt.dryRun:
                        jobcmd = "echo " + jobcmd
                     ret = os.system( jobcmd )



if __name__ == '__main__':
   batchManager = StageXBatchManager()
