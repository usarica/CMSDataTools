#!/bin/env python

import sys
import subprocess
import glob
from optparse import OptionParser


class HadoopChecker:
   def __init__(self):
      # define options and arguments ====================================
      self.parser = OptionParser()

      self.parser.add_option("--maindir", type="string", help="Main directory to check. Uses glob, so need to use regular expressions as necessary")
      self.parser.add_option("--output", type="string", default="", help="Name of the output file to record the list of bad files")
      self.parser.add_option("--remove", action="store_true", default=False, help="Remove the bad files")

      (self.opt,self.args) = self.parser.parse_args()

      if not hasattr(self.opt, "maindir"):
         sys.exit("You must use the --maindir option!")

      if self.opt.remove and not self.opt.output:
         sys.exit("The --remove option can only be specified when an output file is also specified. Exiting...")

      outfile = None
      if self.opt.output != "":
         outfile = open(self.opt.output,'w')

      # Start iterating over the files
      maindir=self.opt.maindir
      print "Checking the expression {} for corruptions".format(maindir)
      indirs = glob.glob(maindir)

      for indir in indirs:
         print " - Checking {}".format(indir)
         p = subprocess.Popen("hdfs fsck {}".format(indir.replace("/hadoop","")).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
         out, err = p.communicate()
         lines = out.split("\n")

         for line in lines:
            if "CORRUPT blockpool" in line:
               theFile = line.split(":")[0]
               badfile = "/hadoop" + theFile
               print "    BAD FILE:",badfile
               if self.opt.remove:
                  print "    => REMOVING THE FILE!"
                  subprocess.call("hdfs dfs -rm " + theFile, shell=True)
               if outfile is not None:
                  outfile.write(badfile)

      if outfile is not None:
         outfile.close()



if __name__ == '__main__':
   theObject = HadoopChecker()
