# No imports needed


def translateFromPythonToShell( strarg ):
   if strarg is None:
      return None
   strarg = strarg.replace("\!",".oO[EXCLAMATION].Oo")
   strarg = strarg.replace("!",".oO[EXCLAMATION].Oo")
   strarg = strarg.replace(" ",".oO[SPACE].Oo")
   strarg = strarg.replace("\"",".oO[DOUBLEQUOTE].Oo")
   strarg = strarg.replace("\\",".oO[BACKLASH].Oo")
   return strarg

def translateFromShellToPython( strarg ):
   if strarg is None:
      return None
   strarg = strarg.replace(".oO[EXCLAMATION].Oo","\!")
   strarg = strarg.replace(".oO[SPACE].Oo"," ")
   strarg = strarg.replace(".oO[DOUBLEQUOTE].Oo","\"")
   strarg = strarg.replace(".oO[BACKLASH].Oo","\\")
   return strarg

def translateROOTArgumentFromShellToPython( strarg ):
   if strarg is None:
      return None
   strarg = strarg.replace(".oO[EXCLAMATION].Oo","!")
   strarg = strarg.replace(".oO[SPACE].Oo"," ")
   strarg = strarg.replace(".oO[DOUBLEQUOTE].Oo",'\\\"')
   strarg = strarg.replace(".oO[BACKLASH].Oo",'\\\\')
   return strarg

