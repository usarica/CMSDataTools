import sys
import os


def getFcnArguments( fname, fcnname ):
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


def mkdir( dirname ):
   strcmd = 'mkdir -p %s' % dirname
   ret = os.system(strcmd)
   if( ret != 0 ):
      print 'Please remove or rename directory: ', dirname
      sys.exit(4)


def rm( target ):
   strcmd = "rm -rf {}".format(target)
   ret = os.system(strcmd)
   if( ret != 0 ):
      print 'Command {} failed!'.format(strcmd)
      sys.exit(4)


def cp( inname, target ):
   strcmd = "cp {} {}".format(inname, target)
   ret = os.system(strcmd)
   if( ret != 0 ):
      print 'Command {} failed!'.format(strcmd)
      sys.exit(4)

