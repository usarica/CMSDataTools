#include "FormFactorHelpers.h"
#include "MELAStreamHelpers.hh"
#include <cassert>


using namespace ACHypothesisHelpers;
using namespace std;
using namespace MELAStreamHelpers;


FormFactorHelpers::FormFactorHandle::FormFactorHandle(FormFactorHelpers::FormFactorType fftype_) : fftype(fftype_), rule(nullptr){
  switch (fftype){
    case ffLQ:
      rule=FormFactorHelpers::getLQFFWeight;
      break;
    default:
      MELAerr << "FormFactorHandle::FormFactorHandle: Form factor " << FormFactorHelpers::getFormFactorName(fftype) << " has no matching function! Aborting..." << endl;
      assert(0);
      break;
  }
}

float FormFactorHelpers::FormFactorHandle::eval(ProcessHandler const* proc, int tpltype, ACHypothesisHelpers::ACHypothesis hypo, std::vector<float*> const& vars) const{ return rule(proc, tpltype, hypo, vars); }

TString FormFactorHelpers::getFormFactorName(FormFactorHelpers::FormFactorType fftype){
  switch (fftype){
  case ffLQ:
    return "LQ";
  default:
    assert(0);
    return "";
  }
}
TString FormFactorHelpers::getFormFactorLabel(FormFactorHelpers::FormFactorType fftype){
  switch (fftype){
  case ffLQ:
    return "f_{#LambdaQ}";
  default:
    assert(0);
    return "";
  }
}


float FormFactorHelpers::getLQFFWeight(ProcessHandler const* proc, int tpltype, ACHypothesisHelpers::ACHypothesis hypo, std::vector<float*> const& vals){
  if (vals.size()!=1){
    MELAerr << "FormFactorHelpers::getLQFFWeight: Variables can only contain GenHMass" << endl;
    assert(0);
  }
  const float& GenHMass=*(vals.front());

  float res=1;
  const GGProcessHandler* proc_GG = dynamic_cast<const GGProcessHandler*>(proc);
  const VVProcessHandler* proc_VV = dynamic_cast<const VVProcessHandler*>(proc);
  const TTProcessHandler* proc_TT = dynamic_cast<const TTProcessHandler*>(proc);
  const BBProcessHandler* proc_BB = dynamic_cast<const BBProcessHandler*>(proc);
  if (proc_GG){
    GGProcessHandler::TemplateType tplcode = proc_GG->castIntToTemplateType(tpltype);
    if (hypo==kSM){ // Special scaling for SM templates
      switch (tplcode){
      case GGProcessHandler::GGTplBkg:
      case GGProcessHandler::GGTplSig:
      case GGProcessHandler::GGTplInt_Re:
        res=1;
        break;
      case GGProcessHandler::GGTplSigBSM:
        res=pow(GenHMass/MHRefVal, 4);
        break;
      case GGProcessHandler::GGTplSigBSMSMInt_Re:
        res=-2.*pow(GenHMass/MHRefVal, 2);
        break;
      case GGProcessHandler::GGTplIntBSM_Re:
        res=-pow(GenHMass/MHRefVal, 2);
        break;
      default:
        break;
      }
    }
    else{
      switch (tplcode){
      case GGProcessHandler::GGTplBkg:
      case GGProcessHandler::GGTplSig:
      case GGProcessHandler::GGTplInt_Re:
        res=1;
        break;
      case GGProcessHandler::GGTplSigBSM:
        res=pow(GenHMass/MHRefVal, 4);
        break;
      case GGProcessHandler::GGTplSigBSMSMInt_Re:
        res=-pow(GenHMass/MHRefVal, 2);
        break;
      case GGProcessHandler::GGTplIntBSM_Re:
        res=-pow(GenHMass/MHRefVal, 2);
        break;
      default:
        break;
      }
    }
    return res;
  }
  else if (proc_VV){
    VVProcessHandler::TemplateType tplcode = proc_VV->castIntToTemplateType(tpltype);
    if (hypo==kSM){ // Special scaling for SM templates
      switch (tplcode){
      case VVProcessHandler::VVTplBkg:
      case VVProcessHandler::VVTplSig:
      case VVProcessHandler::VVTplInt_Re:
        res=1;
        break;
      case VVProcessHandler::VVTplSigBSM:
        res=pow(GenHMass/MHRefVal, 8);
        break;
      case VVProcessHandler::VVTplSigBSMSMInt_ai1_1_Re:
        res=-4.*pow(GenHMass/MHRefVal, 2);
        break;
      case VVProcessHandler::VVTplSigBSMSMInt_ai1_2_PosDef:
        res=6.*pow(GenHMass/MHRefVal, 4);
        break;
      case VVProcessHandler::VVTplSigBSMSMInt_ai1_3_Re:
        res=-4.*pow(GenHMass/MHRefVal, 6);
        break;
      case VVProcessHandler::VVTplIntBSM_ai1_1_Re:
        res=-2.*pow(GenHMass/MHRefVal, 2);
        break;
      case VVProcessHandler::VVTplIntBSM_ai1_2_Re:
        res=pow(GenHMass/MHRefVal, 4);
        break;
      default:
        break;
      }
    }
    else{
      switch (tplcode){
      case VVProcessHandler::VVTplBkg:
      case VVProcessHandler::VVTplSig:
      case VVProcessHandler::VVTplInt_Re:
        res=1;
        break;
      case VVProcessHandler::VVTplSigBSM:
        res=pow(GenHMass/MHRefVal, 8);
        break;
      case VVProcessHandler::VVTplSigBSMSMInt_ai1_1_Re:
        res=-pow(GenHMass/MHRefVal, 2);
        break;
      case VVProcessHandler::VVTplSigBSMSMInt_ai1_2_PosDef:
        res=pow(GenHMass/MHRefVal, 4);
        break;
      case VVProcessHandler::VVTplSigBSMSMInt_ai1_3_Re:
        res=-pow(GenHMass/MHRefVal, 6);
        break;
      case VVProcessHandler::VVTplIntBSM_ai1_1_Re:
        res=-pow(GenHMass/MHRefVal, 2);
        break;
      case VVProcessHandler::VVTplIntBSM_ai1_2_Re:
        res=pow(GenHMass/MHRefVal, 4);
        break;
      default:
        break;
      }
    }
    return res;
  }
  else if (proc_TT){
    TTProcessHandler::TemplateType tplcode = proc_TT->castIntToTemplateType(tpltype);
    if (hypo==kSM){ // Special scaling for SM templates
      switch (tplcode){
      case TTProcessHandler::TTTplSig:
        res=1;
        break;
      case TTProcessHandler::TTTplSigBSM:
        res=pow(GenHMass/MHRefVal, 4);
        break;
      case TTProcessHandler::TTTplSigBSMSMInt_Re:
        res=-2.*pow(GenHMass/MHRefVal, 2);
        break;
      default:
        break;
      }
    }
    else{
      switch (tplcode){
      case TTProcessHandler::TTTplSig:
      case TTProcessHandler::TTTplSigBSM:
        res=pow(GenHMass/MHRefVal, 4);
        break;
      case TTProcessHandler::TTTplSigBSMSMInt_Re:
        res=-pow(GenHMass/MHRefVal, 2);
        break;
      default:
        break;
      }
    }
    return res;
  }
  else if (proc_BB){
    BBProcessHandler::TemplateType tplcode = proc_BB->castIntToTemplateType(tpltype);
    if (hypo==kSM){ // Special scaling for SM templates
      switch (tplcode){
      case BBProcessHandler::BBTplSig:
        res=1;
        break;
      case BBProcessHandler::BBTplSigBSM:
        res=pow(GenHMass/MHRefVal, 4);
        break;
      case BBProcessHandler::BBTplSigBSMSMInt_Re:
        res=-2.*pow(GenHMass/MHRefVal, 2);
        break;
      default:
        break;
      }
    }
    else{
      switch (tplcode){
      case BBProcessHandler::BBTplSig:
      case BBProcessHandler::BBTplSigBSM:
        res=pow(GenHMass/MHRefVal, 4);
        break;
      case BBProcessHandler::BBTplSigBSMSMInt_Re:
        res=-pow(GenHMass/MHRefVal, 2);
        break;
      default:
        break;
      }
    }
    return res;
  }
  else return 1;
}

