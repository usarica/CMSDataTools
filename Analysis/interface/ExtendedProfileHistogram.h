#ifndef EXTENDEDPROFILEHISTOGRAM_H
#define EXTENDEDPROFILEHISTOGRAM_H

#include "ExtendedBinning.h"
#include "TString.h"


class ExtendedProfileHistogram{
  bool fillAverageXYZ;

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
  ExtendedProfileHistogram(ExtendedBinning bX, bool fillAverageXYZ_);
  ExtendedProfileHistogram(ExtendedBinning bX, ExtendedBinning bY, bool fillAverageXYZ_);
  ExtendedProfileHistogram(ExtendedBinning bX, ExtendedBinning bY, ExtendedBinning bZ, bool fillAverageXYZ_);
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

  std::vector<std::vector<std::vector<double>>>& getSumWContainer(){ return sumW; }
  std::vector<std::vector<std::vector<double>>>& getSumWsqContainer(){ return sumWsq; }
  std::vector<std::vector<std::vector<double>>>& getSumWXContainer(){ return sumWX; }
  std::vector<std::vector<std::vector<double>>>& getSumWXsqContainer(){ return sumWXsq; }
  std::vector<std::vector<std::vector<double>>>& getSumWYContainer(){ return sumWY; }
  std::vector<std::vector<std::vector<double>>>& getSumWYsqContainer(){ return sumWYsq; }
  std::vector<std::vector<std::vector<double>>>& getSumWZContainer(){ return sumWZ; }
  std::vector<std::vector<std::vector<double>>>& getSumWZsqContainer(){ return sumWZsq; }

  std::vector<std::vector<std::vector<double>>> const& getSumWContainer() const{ return sumW; }
  std::vector<std::vector<std::vector<double>>> const& getSumWsqContainer() const{ return sumWsq; }
  std::vector<std::vector<std::vector<double>>> const& getSumWXContainer() const{ return sumWX; }
  std::vector<std::vector<std::vector<double>>> const& getSumWXsqContainer() const{ return sumWXsq; }
  std::vector<std::vector<std::vector<double>>> const& getSumWYContainer() const{ return sumWY; }
  std::vector<std::vector<std::vector<double>>> const& getSumWYsqContainer() const{ return sumWYsq; }
  std::vector<std::vector<std::vector<double>>> const& getSumWZContainer() const{ return sumWZ; }
  std::vector<std::vector<std::vector<double>>> const& getSumWZsqContainer() const{ return sumWZsq; }


};


#endif
