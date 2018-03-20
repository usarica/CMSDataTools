#ifndef EXTENDEDPROFILEHISTOGRAM_H
#define EXTENDEDPROFILEHISTOGRAM_H

#include "ExtendedBinning.h"
#include "TString.h"


class ExtendedProfileHistogram{
  ExtendedBinning xbinning;
  ExtendedBinning ybinning;
  ExtendedBinning zbinning;

  std::vector<std::vector<std::vector<double>>> sumW;
  std::vector<std::vector<std::vector<double>>> sumWsq;
  std::vector<std::vector<std::vector<double>>> sumWX;
  std::vector<std::vector<std::vector<double>>> sumWXsq;
  std::vector<std::vector<std::vector<double>>> sumWY;
  std::vector<std::vector<std::vector<double>>> sumWYsq;
  std::vector<std::vector<std::vector<double>>> sumWZ;
  std::vector<std::vector<std::vector<double>>> sumWZsq;

public:
  ExtendedProfileHistogram(ExtendedBinning bX);
  ExtendedProfileHistogram(ExtendedBinning bX, ExtendedBinning bY);
  ExtendedProfileHistogram(ExtendedBinning bX, ExtendedBinning bY, ExtendedBinning bZ);
  ExtendedProfileHistogram(ExtendedProfileHistogram const& other);

  void fill(double x, double w);
  void fill(double x, double y, double w);
  void fill(double x, double y, double z, double w);

  double getBinSumW(int ix, int iy=0, int iz=0);
  double getBinSumWsq(int ix, int iy=0, int iz=0);
  double getBinAvgX(int ix, int iy=0, int iz=0);
  double getBinSigmaX(int ix, int iy=0, int iz=0);
  double getBinAvgY(int ix, int iy=0, int iz=0);
  double getBinSigmaY(int ix, int iy=0, int iz=0);
  double getBinAvgZ(int ix, int iy=0, int iz=0);
  double getBinSigmaZ(int ix, int iy=0, int iz=0);

};


#endif
