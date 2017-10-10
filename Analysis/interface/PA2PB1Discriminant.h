#ifndef PA2PB1DISCRIMINANT_H
#define PA2PB1DISCRIMINANT_H

#include "Discriminant.h"


class PA2PB1Discriminant : public Discriminant{
protected:
  void eval(const std::vector<float>& vars, const float& valReco);

public:
  PA2PB1Discriminant(const TString cfilename="", const TString splinename="sp_gr_varReco_Constant_Smooth");
};

typedef PA2PB1Discriminant DjVBF_t;

#endif
