#ifndef SYSTEMATICSHELPERS_H
#define SYSTEMATICSHELPERS_H

#include "ReweightingBuilder.h"


namespace SystematicsHelpers{

  class ProcessSystematic{
  protected:
    float(*rule)(CJLSTTree*, const std::vector<ReweightingBuilder*>&);

    std::vector<ReweightingBuilder*> evaluators;

  public:
    ProcessSystematic(const std::vector<ReweightingBuilder*>& inEvaluators, float(*infcn)(CJLSTTree*, const std::vector<ReweightingBuilder*>&));
    virtual float eval(CJLSTTree* theTree) const;

  };

  float getRawSystematic(CJLSTTree* theTree, const std::vector<ReweightingBuilder*>& builders);
  float getNormalizedSystematic(CJLSTTree* theTree, const std::vector<ReweightingBuilder*>& builders);

}

#endif
