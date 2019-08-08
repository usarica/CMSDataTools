#!/usr/bin/env python
from CMSDataTools.AnalysisTree.eostools import listFiles
import os
import sys
from optparse import OptionParser


class ListerClass:

   def __init__(self):
      # define options and arguments ====================================
      self.parser = OptionParser()

      self.parser.add_option("--dataset", type="string", help="Data set name")
      self.parser.add_option("--method", type="string", default="dbs", help="Method to list the data files")
      self.parser.add_option("--options", type="string", default=None, help="Other options specific to each method")
      self.parser.add_option("--output", type="string", default=None, help="Output file to record the list")

      (self.opt,self.args) = self.parser.parse_args()

      optchecks=[
         "dataset",
      ]
      for theOpt in optchecks:
         if not hasattr(self.opt, theOpt) or getattr(self.opt, theOpt) is None:
            sys.exit("Need to set --{}".format(theOpt))

      filelist = listFiles(
         sample = self.opt.dataset,
         path = self.opt.method,
         rec = True,
         other_options = self.opt.options
         )

      if self.opt.output is not None:
         foutput = open(self.opt.output, 'w')
         for line in filelist:
            foutput.write(line+'\n')
         foutput.close()
      else:
         #print "File list for {}:".format(self.opt.dataset)
         for line in filelist:
            print line



if __name__ == '__main__':
   theObject = ListerClass()
