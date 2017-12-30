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
      self.parser.add_option("--fixedDate", dest="fixedDate", type="string", help="Fixed output directory", default="")
      self.parser.add_option("--process", dest="process", type="string", help="Name of the process")
      self.parser.add_option("--generator", dest="generator", type="string", help="Name of the generator", default="POWHEG")
      self.parser.add_option("--stage", dest="stage", type="int", default=2, help="Stage 1, 2 or checkstage (x-1) (default=2)")
      self.parser.add_option("--dry", dest="dryRun", type="int", default=0, help="Do not submit jobs, just set up the files")
      self.parser.add_option("--customsyst", dest="customSysts", type="string", action="append", help="Systematics to run (default=all turned on)")

      (self.opt,self.args) = self.parser.parse_args()

      if not hasattr(self.opt, "process"):
         sys.exit("Need to set --process option")

      strscript="make{}TemplatesFrom{}".format(self.opt.process, self.opt.generator)
      self.scriptname="{}.cc".format(strscript)
      if not os.path.isfile(self.scriptname):
         sys.exit("Script {} does not exist. Exiting...".format(self.scriptname))

      self.fcnname=""
      if self.opt.stage==1:
         self.fcnname="{}_one".format(strscript)
      elif self.opt.stage==2:
         self.fcnname="{}_two".format(strscript)
      elif (self.opt.stage==-1 or self.opt.stage==-2):
         self.fcnname="{}_checkstage".format(strscript)
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
      channels = [ "k2e2mu", "k4e", "k4mu" ]
      categories = [ "Inclusive" ] # Not yet ready for tagged categories
      hypos = [ "kSM", "kL!", "kA2", "kA3" ]
      systematics = [ "sNominal", "eLepSFDn", "eLepSFUp", "tPDFScaleDn", "tPDFScaleUp", "tQCDScaleDn", "tQCDScaleUp", "tAsMZDn", "tAsMZUp", "tPDFReplicaDn", "tPDFReplicaUp", "tQQBkgEWCorrDn", "tQQBkgEWCorrUp", "eJECDn", "eJECUp" ]

      for channel in channels:
         for cat in categories:
            for hypo in hypos:
               if self.opt.process == "QQBkg" and hypo != "kSM": # QQbkg is only kSM
                  break
               for syst in systematics:
                  if self.opt.customSysts is not None:
                     if not syst in self.opt.customSysts:
                        continue
                  strscrcmd = "{}, {}".format(channel, cat)
                  if self.opt.process != "QQBkg":
                    strscrcmd = "{}, {}".format(strscrcmd, hypo)
                  strscrcmd = "{}, {}".format(strscrcmd, syst)
                  if "checkstage" in self.fcnname:
                     strscrcmd = "{}, {}".format(strscrcmd, -self.opt.stage)
                  if self.opt.fixedDate:
                     strscrcmd = "{}, \\\"{}\\\"".format(strscrcmd, self.opt.fixedDate)
                  strscrcmd = strscrcmd.replace(' ','') # The command passed to bash script should not contain whitespace itself
                  jobcmd = "submitHiggsWidthTemplateStageGeneric.sh {} \({}\)".format(self.fcnname, strscrcmd)
                  if self.opt.dryRun>0:
                     jobcmd = "echo " + jobcmd
                  ret = os.system( jobcmd )



if __name__ == '__main__':
   batchManager = StageXBatchManager()
