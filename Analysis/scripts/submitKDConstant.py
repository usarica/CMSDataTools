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
from HiggsWidth_PostICHEP.Analysis.ProcessHelpers import *


class StageXBatchManager:
   def __init__(self):
      # define options and arguments ====================================
      self.parser = OptionParser()

      self.parser.add_option("--outdir", dest="outdir", type="string", help="Name of the local output directory for your jobs. This directory will be created automatically.", default="./")

      self.parser.add_option("--KD", dest="KD", type="string", help="KD constant to run (mandatory)")

      self.parser.add_option("--stage", dest="stage", type="int", default=1, help="Stage 1, 2 (default=1)")
      self.parser.add_option("--channel", dest="customChannels", type="string", action="append", help="Channels to run (default=all turned on)")
      self.parser.add_option("--category", dest="customCategories", type="string", action="append", help="Categories to run (default=all turned on)")
      self.parser.add_option("--sqrts", type="int", default=13, help="Sqrts (defaul=13)")
      self.parser.add_option("--writeTrees", action="store_true", default=False, help="Write output trees for debugging")

      self.parser.add_option("--dry", dest="dryRun", action="store_true", default=False, help="Do not submit jobs, just set up the files")
      self.parser.add_option("--interactive", dest="interactive", action="store_true", default=False, help="Do not submit jobs; run them interactively")
      self.parser.add_option("--merge4l", action="store_true", default=False, help="Merge 4e and 4mu channels")

      self.parser.add_option("--batchqueue", type="string", default="default", help="Batch queue")

      (self.opt,self.args) = self.parser.parse_args()

      if not hasattr(self.opt, "KD"):
         sys.exit("Need to set --KD option")

      strscript="calculateKDconstants"
      self.scriptname="{}.cc".format(strscript)
      if not os.path.isfile(self.scriptname):
         sys.exit("Script {} does not exist. Exiting...".format(self.scriptname))

      self.fcnname=""
      if self.opt.stage==1:
         self.fcnname="getKDConstant_{}".format(self.opt.KD)
      elif self.opt.stage==2:
         self.fcnname="SmoothKDConstantProducer_{}".format(self.opt.KD)
      if not self.fcnname:
         sys.exit("The function name could not be generated. Exiting...")

      mkdir(self.opt.outdir)
      if not os.path.isfile(self.opt.outdir + '/loadLib.C') or not filecmp.cmp('loadLib.C',self.opt.outdir + '/loadLib.C'):
         print 'Need a new loadLib.C file'
         cp('loadLib.C',self.opt.outdir + '/loadLib.C')
      rm(self.opt.outdir + '/' + self.fcnname + ".c")
      rm(self.opt.outdir + '/' + self.fcnname + "_c.d")
      rm(self.opt.outdir + '/' + self.fcnname + "_c.so")
      rm(self.opt.outdir + '/' + self.fcnname + "_c_ACLiC_dict_rdict.pcm")
      rm(self.opt.outdir + '/' + self.fcnname + ".c")
      cp(self.scriptname, self.opt.outdir + '/' + self.fcnname + ".c")

      currentdir = os.getcwd()
      os.chdir(self.opt.outdir)
      print "Job submission directory: {}".format(os.getcwd())
      self.submitJobs()
      os.chdir(currentdir)
      print "Working directory: {}".format(os.getcwd())


   def submitJobs(self):
      channels = [ "NChannels", "k2e2mu", "k4e", "k4mu" ]
      if self.opt.merge4l:
         channels = [ "k2l2l", "k4l" ]
      categories = [ "Inclusive", "JJVBFTagged", "HadVHTagged" ]

      fcnargs=getFcnArguments(self.scriptname, self.fcnname)
      argstr=""
      for fcnarg in fcnargs:
         tmpargstr=""
         if "channel" in fcnarg:
            tmpargstr = "{channel}"
         elif "category" in fcnarg:
            tmpargstr = "{category}"
         elif "writefinaltree" in fcnarg:
            tmpargstr = "true" if self.opt.writeTrees else "false"
         elif "sqrts" in fcnarg:
            tmpargstr = str(self.opt.sqrts)
         if tmpargstr:
            if argstr:
               argstr = "{},{}".format(argstr, tmpargstr)
            else:
               argstr=tmpargstr

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
               if cat=="Inclusive":
                  continue
               elif self.opt.customCategories is not None:
                  if not cat in self.opt.customCategories:
                     continue

            strscrcmd = argstr.format(channel=ch,category=cat)
            strscrcmd = strscrcmd.replace(' ','') # The command passed to bash script should not contain whitespace itself
            jobcmd = "submitHiggsWidthTemplateStageGeneric.sh {} \({}\) {}".format(self.fcnname, strscrcmd, self.opt.batchqueue)
            if self.opt.interactive:
               jobcmd = "root -l -b -q -e \"gROOT->ProcessLine(\\\".x loadLib.C\\\");gROOT->ProcessLine(\\\".x {}.c+({})\\\");\"".format(self.fcnname, strscrcmd)
            if self.opt.dryRun:
               jobcmd = "echo " + jobcmd
            ret = os.system( jobcmd )



if __name__ == '__main__':
   batchManager = StageXBatchManager()
