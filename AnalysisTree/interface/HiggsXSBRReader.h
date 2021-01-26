#ifndef HIGGSXSBRREADER_H
#define HIGGSXSBRREADER_H

#include <string>
#include <vector>
#include "TString.h"
#include "TSpline.h"


class HiggsXSBRReader{
protected:
  std::vector<double> masses;
  std::vector<double> total_widths;
  std::vector<double> partial_widths;

  TSpline3 sp_partial_width;
  TSpline3 sp_total_width;

public:
  HiggsXSBRReader(TString fname, TString partial_width_type);
  float eval_partial_width(float const& mass) const;
  float eval_total_width(float const& mass) const;
  float eval_br(float const& mass) const{ return eval_partial_width(mass) / eval_total_width(mass); }

};

#endif
