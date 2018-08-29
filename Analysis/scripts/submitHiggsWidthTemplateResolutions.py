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


class ResolutionHarvesterBatchManager:
   def __init__(self):
      # define options and arguments ====================================
      self.parser = OptionParser()

      self.parser.add_option("--outdir", dest="outdir", type="string", help="Name of the local output directory for your jobs. This directory will be created automatically.", default="./")

      self.parser.add_option("--process", dest="process", type="string", help="Name of the process")
      self.parser.add_option("--generator", dest="generator", type="string", help="Name of the generator")
      self.parser.add_option("--fixedDate", dest="fixedDate", type="string", help="Fixed output directory", default="")

      self.parser.add_option("--channel", dest="customChannels", type="string", action="append", help="Channels to run (default=all turned on)")
      self.parser.add_option("--category", dest="customCategories", type="string", action="append", help="Categories to run (default=all turned on)")

      self.parser.add_option("--dry", dest="dryRun", action="store_true", default=False, help="Do not submit jobs, just set up the files")
      self.parser.add_option("--interactive", dest="interactive", action="store_true", default=False, help="Do not submit jobs; run them interactively")
      self.parser.add_option("--H125only", action="store_true", default=False, help="Parameterize H125 shape only")
      self.parser.add_option("--norecompile", action="store_true", default=False, help="Do not remove executable and shared objects")

      self.parser.add_option("--batchqueue", type="string", default="default", help="Batch queue")

      (self.opt,self.args) = self.parser.parse_args()

      processDict = processDictionary()
      del processDict["ZX"]

      if self.opt.process is None:
         sys.exit("Need to set --process option")
      elif processDict[self.opt.process] is None:
         sys.exit("Process {} is not recognized".format(self.opt.process))
      elif self.opt.process != "ZX" and self.opt.generator is None:
         sys.exit("Need to set --generator option")
      self.process=processDict[self.opt.process]

      self.generator = self.opt.generator

      strscript="acquireResolution.cc"
      self.fcnname=""
      if self.opt.H125only:
         self.fcnname="acquireH125OnshellMassShape_one"
      else:
         self.fcnname="acquireResolution_one"
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

      fcnargs=getFcnArguments(self.scriptname, self.fcnname)
      argstr=""
      for fcnarg in fcnargs:
         tmpargstr=""
         if "channel" in fcnarg:
            tmpargstr = "{channel}"
         elif "category" in fcnarg:
            tmpargstr = "{category}"
         elif "proctype" in fcnarg:
            tmpargstr = "{processtype}"
         # For the rest of if-statements, do not set tmpargstr; append to argstr directly
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

            # Do not submit unnecessary jobs
            strscrcmd = argstr.format(channel=ch,category=cat,processtype=self.process)
            strscrcmd = strscrcmd.replace(' ','') # The command passed to bash script should not contain whitespace itself
            jobcmd = "submitHiggsWidthTemplateStageGeneric.sh {} \({}\) {}".format(self.fcnname, strscrcmd, self.opt.batchqueue)
            if self.opt.interactive:
               jobcmd = "root -l -b -q -e \"gROOT->ProcessLine(\\\".x loadLib.C\\\");gROOT->ProcessLine(\\\".x {}.c+({})\\\");\"".format(self.fcnname, strscrcmd)
            if self.opt.dryRun:
               jobcmd = "echo " + jobcmd
            ret = os.system( jobcmd )



if __name__ == '__main__':
   batchManager = ResolutionHarvesterBatchManager()
