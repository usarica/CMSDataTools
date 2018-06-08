#include <cassert>
#include "ACHypothesisHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


TString ACHypothesisHelpers::getACHypothesisName(ACHypothesisHelpers::ACHypothesis hypo){
  switch (hypo){
  case kSM:
    return "SM";
  case kL1:
    return "L1";
  case kA2:
    return "a2";
  case kA3:
    return "a3";
  case kL1ZGs:
    return "L1ZGs";
  default:
    return "";
  };
}

std::vector<DiscriminantClasses::Type> ACHypothesisHelpers::getACHypothesisKDSet(ACHypothesisHelpers::ACHypothesis hypo, CategorizationHelpers::Category category, CategorizationHelpers::MassRegion massregion){
  std::vector<DiscriminantClasses::Type> res;
  if (massregion==CategorizationHelpers::kOffshell){ // First dimension, "mass", is added separately
    if (category==CategorizationHelpers::Inclusive || category==CategorizationHelpers::Untagged){
      res.push_back(DiscriminantClasses::kDbkgkin);
      switch (hypo){
      case kSM:
        res.push_back(DiscriminantClasses::kCggint);
        break;
      case kL1:
        res.push_back(DiscriminantClasses::kDL1dec);
        break;
      case kA2:
        res.push_back(DiscriminantClasses::kDa2dec);
        break;
      case kA3:
        res.push_back(DiscriminantClasses::kDa3dec);
        break;
      case kL1ZGs:
        res.push_back(DiscriminantClasses::kDL1ZGsdec);
        break;
      default:
        break;
      };
    }
    else if (category==CategorizationHelpers::JJVBFTagged){
      res.push_back(DiscriminantClasses::kDbkgjjEWQCD);
      switch (hypo){
      case kSM:
        res.push_back(DiscriminantClasses::kCjjVBFint);
        break;
      case kL1:
        res.push_back(DiscriminantClasses::kDL1jjVBFdec);
        break;
      case kA2:
        res.push_back(DiscriminantClasses::kDa2jjVBFdec);
        break;
      case kA3:
        res.push_back(DiscriminantClasses::kDa3jjVBFdec);
        break;
      case kL1ZGs:
        res.push_back(DiscriminantClasses::kDL1ZGsjjVBFdec);
        break;
      default:
        break;
      };
    }
    else if (category==CategorizationHelpers::HadVHTagged){
      res.push_back(DiscriminantClasses::kDbkgjjEWQCD);
      switch (hypo){
      case kSM:
        res.push_back(DiscriminantClasses::kCjjVHint);
        break;
      case kL1:
        res.push_back(DiscriminantClasses::kDL1jjVHdec);
        break;
      case kA2:
        res.push_back(DiscriminantClasses::kDa2jjVHdec);
        break;
      case kA3:
        res.push_back(DiscriminantClasses::kDa3jjVHdec);
        break;
      case kL1ZGs:
        res.push_back(DiscriminantClasses::kDL1ZGsjjVHdec);
        break;
      default:
        break;
      };
    }
  }
  else if (massregion==CategorizationHelpers::kOnshell){
    if (category==CategorizationHelpers::Inclusive || category==CategorizationHelpers::Untagged){
      switch (hypo){
      case kSM: // First dimension, "mass", is added separately
        res.push_back(DiscriminantClasses::kDbkgkin);
        break;
      case kL1:
        res.push_back(DiscriminantClasses::kDbkgm4l);
        res.push_back(DiscriminantClasses::kDL1dec);
        res.push_back(DiscriminantClasses::kDa2dec);
        break;
      case kA2:
        res.push_back(DiscriminantClasses::kDbkgm4l);
        res.push_back(DiscriminantClasses::kDa2dec);
        res.push_back(DiscriminantClasses::kCa2decint);
        break;
      case kA3:
        res.push_back(DiscriminantClasses::kDbkgm4l);
        res.push_back(DiscriminantClasses::kDa3dec);
        res.push_back(DiscriminantClasses::kCa3decint);
        break;
      case kL1ZGs:
        res.push_back(DiscriminantClasses::kDbkgm4l);
        res.push_back(DiscriminantClasses::kDL1ZGsdec);
        res.push_back(DiscriminantClasses::kDa2dec);
        //res.push_back(DiscriminantClasses::kCL1ZGsdecint);
        break;
      default:
        break;
      };
    }
    else if (category==CategorizationHelpers::JJVBFTagged){
      switch (hypo){
      case kSM: // First dimension, "mass", is added separately
        res.push_back(DiscriminantClasses::kDbkgjjEWQCD);
        break;
      case kL1:
        res.push_back(DiscriminantClasses::kDbkgm4ljjEWQCD);
        res.push_back(DiscriminantClasses::kDL1jjVBFdec);
        res.push_back(DiscriminantClasses::kDa2jjVBFdec);
        break;
      case kA2:
        res.push_back(DiscriminantClasses::kDbkgm4ljjEWQCD);
        res.push_back(DiscriminantClasses::kDa2jjVBFdec);
        res.push_back(DiscriminantClasses::kCa2jjVBFint);
        break;
      case kA3:
        res.push_back(DiscriminantClasses::kDbkgm4ljjEWQCD);
        res.push_back(DiscriminantClasses::kDa3jjVBFdec);
        res.push_back(DiscriminantClasses::kCa3jjVBFint);
        break;
      case kL1ZGs:
        res.push_back(DiscriminantClasses::kDbkgm4ljjEWQCD);
        res.push_back(DiscriminantClasses::kDL1ZGsjjVBFdec);
        res.push_back(DiscriminantClasses::kDa2jjVBFdec);
        //res.push_back(DiscriminantClasses::kCL1ZGsjjVBFint);
        break;
      default:
        break;
      };
    }
    else if (category==CategorizationHelpers::HadVHTagged){
      switch (hypo){
      case kSM: // First dimension, "mass", is added separately
        res.push_back(DiscriminantClasses::kDbkgjjEWQCD);
        break;
      case kL1:
        res.push_back(DiscriminantClasses::kDbkgm4ljjEWQCD);
        res.push_back(DiscriminantClasses::kDL1jjVHdec);
        res.push_back(DiscriminantClasses::kDa2jjVHdec);
        break;
      case kA2:
        res.push_back(DiscriminantClasses::kDbkgm4ljjEWQCD);
        res.push_back(DiscriminantClasses::kDa2jjVHdec);
        res.push_back(DiscriminantClasses::kCa2jjVHint);
        break;
      case kA3:
        res.push_back(DiscriminantClasses::kDbkgm4ljjEWQCD);
        res.push_back(DiscriminantClasses::kDa3jjVHdec);
        res.push_back(DiscriminantClasses::kCa3jjVHint);
        break;
      case kL1ZGs:
        res.push_back(DiscriminantClasses::kDbkgm4ljjEWQCD);
        res.push_back(DiscriminantClasses::kDL1ZGsjjVHdec);
        res.push_back(DiscriminantClasses::kDa2jjVHdec);
        //res.push_back(DiscriminantClasses::kCL1ZGsjjVHint);
        break;
      default:
        break;
      };
    }
  }
  else{
    MELAerr << "ACHypothesisHelpers::getACHypothesisKDSet: Mass region " << massregion << " is not implemented." << endl;
    assert(0);
  }

  return res;
}

std::vector<TString> ACHypothesisHelpers::getACHypothesisKDNameSet(ACHypothesisHelpers::ACHypothesis hypo, CategorizationHelpers::Category category, CategorizationHelpers::MassRegion massregion){
  std::vector<DiscriminantClasses::Type> KDset = ACHypothesisHelpers::getACHypothesisKDSet(hypo, category, massregion);
  vector<TString> res;
  for (auto& type:KDset) res.push_back(DiscriminantClasses::getKDName(type));
  return res;
}

float ACHypothesisHelpers::getACHypothesisMEHZZGVal(ACHypothesisHelpers::ACHypothesis hypo){
  if (hypo==kL1 || hypo==kL1ZGs) return 1e4;
  else return 1;
}
float ACHypothesisHelpers::getACHypothesisHZZGVal(ACHypothesisHelpers::ACHypothesis hypo){
  switch (hypo){
  case kL1:
    return -1.211020e4;
  case kA2:
    return 1.663195;
  case kA3:
    return 2.55502;
  case kL1ZGs:
    return -7613.351;
  default:
    return 1;
  };
}
