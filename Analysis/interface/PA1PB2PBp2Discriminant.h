#ifndef PA1PB2PBP2DISCRIMINANT_H
#define PA1PB2PBP2DISCRIMINANT_H

#include "Discriminant.h"


class PA1PB2PBp2Discriminant : public Discriminant{
protected:
  void eval(const std::vector<float>& vars, const float& valReco);

public:
  PA1PB2PBp2Discriminant(const TString cfilename="", const TString splinename="sp_gr_varReco_Constant_Smooth");
};

typedef PA1PB2PBp2Discriminant Dbkgdec_t;

#endif
