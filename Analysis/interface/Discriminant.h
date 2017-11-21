#ifndef DISCRIMINANT_H
#define DISCRIMINANT_H

#include <vector>
#include "TFile.h"
#include "TString.h"
#include "TSpline.h"


class Discriminant{
protected:
  TFile* theCFile;
  TSpline3* theC;
  float val;

  virtual void eval(const std::vector<float>& vars, const float& valReco)=0;
  float getCval(const float valReco) const;

public:
  Discriminant(const TString cfilename="", const TString splinename="sp_gr_varReco_Constant_Smooth");
  virtual ~Discriminant();

  operator float() const;
  operator float&();
  operator float*();

  bool operator<(const float& other) const;
  bool operator>(const float& other) const;
  bool operator<=(const float& other) const;
  bool operator>=(const float& other) const;
  bool operator==(const float& other) const;
  bool operator!=(const float& other) const;

  float update(const std::vector<float>& vars, const float valReco);
  float applyAdditionalC(const float cval);

};


#endif
