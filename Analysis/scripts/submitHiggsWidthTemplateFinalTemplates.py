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
from HiggsWidth_PostICHEP.Analysis.ProcessOrganization import *


class FinalTemplatesStageXBatchManager:
   def __init__(self):
      # define options and arguments ====================================
      self.parser = OptionParser()

      self.parser.add_option("--outdir", dest="outdir", type="string", help="Name of the local output directory for your jobs. This directory will be created automatically.", default="./")

      self.parser.add_option("--process", dest="process", type="string", help="Name of the process")
      self.parser.add_option("--stage", dest="stage", type="int", default=1, help="Stage 1, 2 (default=1)")
      self.parser.add_option("--fixedDate", dest="fixedDate", type="string", help="Fixed output directory", default="")

      self.parser.add_option("--syst", dest="customSysts", type="string", action="append", help="Systematics to run (default=all turned on)")
      self.parser.add_option("--channel", dest="customChannels", type="string", action="append", help="Channels to run (default=all turned on)")
      self.parser.add_option("--AChypo", dest="customACHypos", type="string", action="append", help="Anomalous couplings hypotheses to run (default=all turned on)")
      self.parser.add_option("--anaregion", dest="customMassRegions", type="string", action="append", help="Analysis region to produce templates")

      self.parser.add_option("--dry", dest="dryRun", action="store_true", default=False, help="Do not submit jobs, just set up the files")
      self.parser.add_option("--interactive", dest="interactive", action="store_true", default=False, help="Do not submit jobs; run them interactively")
      self.parser.add_option("--norecompile", action="store_true", default=False, help="Do not remove executable and shared objects")

      (self.opt,self.args) = self.parser.parse_args()

      if self.opt.process is None:
         sys.exit("Need to set --process option")

      strscript="makeFinalTemplates_{}".format(self.opt.process)
      self.scriptname="{}.cc".format(strscript)
      self.fcnname=strscript
      if not os.path.isfile(self.scriptname):
         sys.exit("Script {} does not exist. Exiting...".format(self.scriptname))
      if not self.fcnname:
         sys.exit("The function name could not be generated. Exiting...")

      self.mkdir(self.opt.outdir)
      self.cpscriptnamebare=self.fcnname + "_tmp.c"
      self.cpscriptname=self.opt.outdir + '/' + self.cpscriptnamebare
      if not (os.path.isfile(self.cpscriptname) and self.opt.norecompile):
         self.rm(self.opt.outdir + '/' + self.fcnname + "_tmp_c.d")
         self.rm(self.opt.outdir + '/' + self.fcnname + "_tmp_c.so")
         self.rm(self.opt.outdir + '/' + self.fcnname + "_tmp_c_ACLiC_dict_rdict.pcm")
         self.rm(self.cpscriptname)
         self.cp(self.scriptname, self.cpscriptname)
      else:
         print "Copied script {} already exists, will not recompile".format(self.cpscriptname)
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
      anaregions = getAnalysisRegions()
      channels = getChannelList()
      hypos = getACHypothesisList()
      systematics = getSystematicsList()

      fcnargs=self.getFcnArguments(self.scriptname, self.fcnname)
      argstr=""
      for fcnarg in fcnargs:
         tmpargstr=""
         if "channel" in fcnarg:
            tmpargstr = "{channel}"
         elif "syst" in fcnarg:
            tmpargstr = "{systematic}"
         elif "hypo" in fcnarg:
            tmpargstr = "{achypothesis}"
         elif "massregion" in fcnarg:
            tmpargstr = "{anaregion}"
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
                  strfixedDate="\\\\\\\"{}\\\\\\\"".format(self.opt.fixedDate)
               else:
                  strfixedDate="\\\\\\\"\\\\\\\""
               if argstr:
                  argstr = "{},{}".format(argstr, strfixedDate)
               else:
                  argstr=strfixedDate
            else:
               if self.opt.fixedDate:
                  strfixedDate=r"\\\\\\\"{}\\\\\\\"".format(self.opt.fixedDate)
               else:
                  strfixedDate=r"\\\\\\\"\\\\\\\""
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

      for anreg in anaregions:
         if (not "anaregion" in argstr) and anreg!="NMassRegions":
            break
         elif ("anaregion" in argstr):
            if anreg=="NMassRegions":
               continue
            elif self.opt.customMassRegions is not None:
               if not anreg in self.opt.customMassRegions:
                  continue

         for ch in channels:
            if (not "channel" in argstr) and ch!="NChannels":
               break
            elif ("channel" in argstr):
               if ch=="NChannels":
                  continue
               elif self.opt.customChannels is not None:
                  if not ch in self.opt.customChannels:
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

                  if not checkValidRun(syst, "", ch, self.opt.process): continue
                  if (anreg=="kOffshell" or hypo=="kSM") and ("LepScale" in syst or "LepRes" in syst):
                     print "{} systematic distributions in process {} are not handled through templates in hypothesis {} of analysis region {}.".format(syst, self.opt.process, hypo, anreg)
                     continue

                  strscrcmd = argstr.format(channel=ch,achypothesis=hypo,systematic=syst,anaregion=anreg)
                  strscrcmd = strscrcmd.replace(' ','') # The command passed to bash script should not contain whitespace itself
                  jobcmd = "submitHiggsWidthROOTCommand.sh {} {} {}".format(self.cpscriptnamebare, self.fcnname, strscrcmd)
                  if self.opt.interactive:
                     jobcmd = "root -l -b -q -e \"gROOT->ProcessLine(\\\".x loadLib.C\\\");gROOT->ProcessLine(\\\".L {}+\\\");gROOT->ProcessLine(\\\"{}({})\\\");\"".format(self.cpscriptname, self.fcnname, strscrcmd)
                  if self.opt.dryRun:
                     jobcmd = "echo " + jobcmd
                  ret = os.system( jobcmd )



if __name__ == '__main__':
   batchManager = FinalTemplatesStageXBatchManager()
