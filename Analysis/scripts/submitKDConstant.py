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


      self.parser.add_option("--KD", dest="KD", type="string", help="KD constant to run (mandatory)")

      self.parser.add_option("--stage", dest="stage", type="int", default=1, help="Stage 1, 2 (default=1)")
      self.parser.add_option("--channel", dest="customChannels", type="string", action="append", help="Channels to run (default=all turned on)")
      self.parser.add_option("--category", dest="customCategories", type="string", action="append", help="Categories to run (default=all turned on)")
      self.parser.add_option("--sqrts", type="int", default=13, help="Sqrts (defaul=13)")
      self.parser.add_option("--writeTrees", action="store_true", default=False, help="Write output trees for debugging")

      self.parser.add_option("--dry", dest="dryRun", action="store_true", default=False, help="Do not submit jobs, just set up the files")
      self.parser.add_option("--interactive", dest="interactive", action="store_true", default=False, help="Do not submit jobs; run them interactively")

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
      with open(fname) as testfile:
         for line in testfile:
            if fcnname in line:
               fcnfound=True
               line=line.rstrip()
               line=line.replace(' ','')
               line=line.replace(fcnname,'')
               line=line.replace('(','')
               line=line.replace(')','')
               line=line.replace('{','')
               tmpargs=line.split(',')
               for tmparg in tmpargs:
                  fcnargs.append(tmparg.lower())
      if not fcnfound:
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
      categories = [ "Inclusive", "JJVBFTagged", "HadVHTagged" ]

      fcnargs=self.getFcnArguments(self.scriptname, self.fcnname)
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
            jobcmd = "submitHiggsWidthTemplateStageGeneric.sh {} \({}\)".format(self.fcnname, strscrcmd)
            if self.opt.interactive:
               jobcmd = "root -l -b -q -e \"gROOT->ProcessLine(\\\".x loadLib.C\\\");gROOT->ProcessLine(\\\".x {}.c+({})\\\");\"".format(self.fcnname, strscrcmd)
            if self.opt.dryRun:
               jobcmd = "echo " + jobcmd
            ret = os.system( jobcmd )



if __name__ == '__main__':
   batchManager = StageXBatchManager()
