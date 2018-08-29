#!/bin/env python

import sys
import imp
import copy
import os
import filecmp
import shutil
import pickle
import math
import pprint
import subprocess
from datetime import date
from optparse import OptionParser
from HiggsWidth_PostICHEP.Analysis.ProcessOrganization import *
from HiggsWidth_PostICHEP.Analysis.ProcessHelpers import *


class MassRatioStageXBatchManager:
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

      self.parser.add_option("--docatratio", action="store_true", default=False, help="Do categorization ratios for nominal systematic")
      self.parser.add_option("--docatsystratio", action="store_true", default=False, help="Do categorization systematic ratios against nominal systematic")

      self.parser.add_option("--dry", dest="dryRun", action="store_true", default=False, help="Do not submit jobs, just set up the files")
      self.parser.add_option("--interactive", dest="interactive", action="store_true", default=False, help="Do not submit jobs; run them interactively")
      self.parser.add_option("--norecompile", action="store_true", default=False, help="Do not remove executable and shared objects")

      self.parser.add_option("--batchqueue", type="string", default="default", help="Batch queue")

      (self.opt,self.args) = self.parser.parse_args()

      processDict = processDictionary()

      if self.opt.process is None:
         sys.exit("Need to set --process option")
      elif processDict[self.opt.process] is None:
         sys.exit("Process {} is not recognized".format(self.opt.process))
      elif self.opt.process != "ZX" and self.opt.generator is None:
         sys.exit("Need to set --generator option")
      self.process=processDict[self.opt.process]

      if self.opt.process == "ZX":
         self.generator = "Data"
      else:
         self.generator = self.opt.generator

      strscript="acquireProcessMassRatios.cc"
      self.fcnname=""
      if self.opt.docatratio and self.opt.docatsystratio:
         sys.exit("Only specify --docatratio or --docatsystratio")
      elif self.opt.docatratio:
         self.fcnname = "acquireMassRatio_ProcessNominalToNominalInclusive_one"
      elif self.opt.docatsystratio:
         self.fcnname = "acquireMassRatio_ProcessSystToNominal_one"
      if not self.fcnname:
         sys.exit("The function name could not be generated. Exiting...")

      self.scriptname=strscript
      if not os.path.isfile(self.scriptname):
         sys.exit("Script {} does not exist. Exiting...".format(self.scriptname))

      mkdir(self.opt.outdir)
      if not os.path.isfile(self.opt.outdir + '/loadLib.C') or not filecmp.cmp('loadLib.C',self.opt.outdir + '/loadLib.C'):
         print 'Need a new loadLib.C file'
         cp('loadLib.C',self.opt.outdir + '/loadLib.C')
      self.cpscriptnamebare=self.fcnname + ".c"
      cpscriptname=self.opt.outdir + '/' + self.cpscriptnamebare
      if not (os.path.isfile(cpscriptname) and self.opt.norecompile):
         rm(self.opt.outdir + '/' + self.fcnname + "_c.d")
         rm(self.opt.outdir + '/' + self.fcnname + "_c.so")
         rm(self.opt.outdir + '/' + self.fcnname + "_c_ACLiC_dict_rdict.pcm")
         rm(cpscriptname)
         cp(self.scriptname, cpscriptname)
      else:
         print "Copied script {} already exists, will not recompile".format(cpscriptname)

      currentdir = os.getcwd()
      os.chdir(self.opt.outdir)
      print "Job submission directory: {}".format(os.getcwd())
      self.submitJobs()
      os.chdir(currentdir)
      print "Working directory: {}".format(os.getcwd())


   def submitJobs(self):
      channels = getChannelList()
      categories = getCategoryList()
      hypos = getACHypothesisList()
      systematics = getSystematicsList()

      fcnargs=getFcnArguments(self.cpscriptnamebare, self.fcnname)
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
         elif "proctype" in fcnarg:
            tmpargstr = "{processtype}"
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
            if self.opt.interactive:
               if self.opt.fixedDate:
                  strfixedDate=r"\\\"{}\\\"".format(self.opt.fixedDate)
               else:
                  strfixedDate=r"\\\"\\\""
            else:
               if self.opt.fixedDate:
                  strfixedDate="\\\"{}\\\"".format(self.opt.fixedDate)
               else:
                  strfixedDate="\\\"\\\""
            if argstr:
               argstr = "{},{}".format(argstr, strfixedDate)
            else:
               argstr=strfixedDate
         elif "strgenerator" in fcnarg:
            if self.opt.interactive:
               strGenerator=r"\\\"{}\\\"".format(self.generator)
            else:
               strGenerator="\\\"{}\\\"".format(self.generator)
            if argstr:
               argstr = "{},{}".format(argstr, strGenerator)
            else:
               argstr=strGenerator
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
            if self.opt.docatratio and cat=="Inclusive": continue

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
                  if self.opt.docatratio and syst!="sNominal": continue
                  if self.opt.docatsystratio and syst=="sNominal": continue

                  # Do not submit unnecessary jobs
                  if not checkValidRun(syst, cat, ch, self.opt.process, self.generator): continue
                  if self.opt.stage == 1 and cat == "Untagged" and (not(self.opt.process == "ZH" or self.opt.process == "WH") or "tPythia" in syst or "tMINLO" in syst):
                     print "{} category distributions in process {} can be obtained from the distributions of inclusive and other categories.".format(cat, self.opt.process)
                     continue
                  if "LepScale" in syst or "LepRes" in syst:
                     print "{} systematic distributions in process {} are obtained in resolutions step.".format(syst, self.opt.process)
                     continue

                  strscrcmd = argstr.format(channel=ch,category=cat,achypothesis=hypo,systematic=syst,processtype=self.process)
                  strscrcmd = strscrcmd.replace(' ','') # The command passed to bash script should not contain whitespace itself
                  jobcmd = "submitHiggsWidthTemplateStageGeneric.sh {} \({}\) {}".format(self.fcnname, strscrcmd, self.opt.batchqueue)
                  if self.opt.interactive:
                     jobcmd = "root -l -b -q -e \"gROOT->ProcessLine(\\\".x loadLib.C\\\");gROOT->ProcessLine(\\\".x {}.c+({})\\\");\"".format(self.fcnname, strscrcmd)
                  if self.opt.dryRun:
                     jobcmd = "echo " + jobcmd
                  ret = os.system( jobcmd )



if __name__ == '__main__':
   batchManager = MassRatioStageXBatchManager()
