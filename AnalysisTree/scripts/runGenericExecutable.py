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
from CMSDataTools.AnalysisTree.TranslateStringBetweenPythonAndShell import *


class GenericExecutor:
   def __init__(self):
      # define options and arguments ====================================
      self.parser = OptionParser()

      self.parser.add_option("--executable", "--exe", dest="executable", type="string", help="Name of the executable")
      self.parser.add_option("--command", dest="fcncmd", type="string", help="Function arguments", default="")
      self.parser.add_option("--dry", dest="dryRun", action="store_true", default=False, help="Do not submit jobs, just set up the files")

      (self.opt,self.args) = self.parser.parse_args()

      if self.opt.executable is None:
         sys.exit("Need to set --executable option")

      self.run()


   def run(self):
      if self.opt.fcncmd is not None:
         self.opt.fcncmd = translateFromShellToPython(self.opt.fcncmd)
      jobcmd = r""
      jobcmd += r'{} {}'.format(self.opt.executable, self.opt.fcncmd)
      print("Running {}".format(jobcmd))
      if not self.opt.dryRun:
         ret = os.system( jobcmd )



if __name__ == '__main__':
   batchManager = GenericExecutor()
