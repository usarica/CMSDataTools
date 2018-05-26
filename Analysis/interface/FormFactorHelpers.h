#ifndef FORMFACTORHELPERS_H
#define FORMFACTORHELPERS_H

#include "TString.h"
#include "ProcessHandler.h"
#include "ACHypothesisHelpers.h"


namespace FormFactorHelpers{
  enum FormFactorType{
    ffLQ,
    nFormFactorTypes
  };

  typedef float(*FormFactorFunction_t)(ProcessHandler const*, int, ACHypothesisHelpers::ACHypothesis, std::vector<float*> const&);

  struct FormFactorHandle{
    FormFactorType fftype;
    FormFactorFunction_t rule;

    FormFactorHandle(FormFactorHelpers::FormFactorType fftype_);

    float eval(ProcessHandler const* proc, int tpltype, ACHypothesisHelpers::ACHypothesis hypo, std::vector<float*> const& vars) const;
  };

  TString getFormFactorName(FormFactorHelpers::FormFactorType fftype);
  TString getFormFactorLabel(FormFactorHelpers::FormFactorType fftype);

  float getLQFFWeight(ProcessHandler const* proc, int tpltype, ACHypothesisHelpers::ACHypothesis hypo, std::vector<float*> const& vals);

}


#endif
