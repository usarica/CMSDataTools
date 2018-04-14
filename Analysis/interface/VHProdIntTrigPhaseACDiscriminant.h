#ifndef VHPRODINTTRIGPHASEACDISCRIMINANT_H
#define VHPRODINTTRIGPHASEACDISCRIMINANT_H

#include "Discriminant.h"


class VHProdIntTrigPhaseACDiscriminant : public Discriminant{
protected:
  void eval(const std::vector<float>& vars, const float& valReco);

public:
  VHProdIntTrigPhaseACDiscriminant();

};

typedef VHProdIntTrigPhaseACDiscriminant CaiVHint_t;

#endif
