#ifndef EXTENDEDPROFILEHISTOGRAM_H
#define EXTENDEDPROFILEHISTOGRAM_H

#include "ExtendedBinning.h"
#include "TString.h"


class ExtendedProfileHistogram{
  // Can store up to 3 average quantities 
  unsigned int fillAverageQ;

  ExtendedBinning xbinning;
  ExtendedBinning ybinning;
  ExtendedBinning zbinning;

  std::vector<std::vector<std::vector<double>>> sumW;
  std::vector<std::vector<std::vector<double>>> sumWsq;
  std::vector<std::vector<std::vector<double>>>* sumWQ;
  std::vector<std::vector<std::vector<double>>>* sumWQsq;

public:
  ExtendedProfileHistogram(ExtendedBinning bX, unsigned int fillAverageQ_);
  ExtendedProfileHistogram(ExtendedBinning bX, ExtendedBinning bY, unsigned int fillAverageQ_);
  ExtendedProfileHistogram(ExtendedBinning bX, ExtendedBinning bY, ExtendedBinning bZ, unsigned int fillAverageQ_);
  ExtendedProfileHistogram(ExtendedProfileHistogram const& other);
  ~ExtendedProfileHistogram();

  void fill(double x, double w);
  void fill(double x, double y, double w);
  void fill(double x, double y, double z, double w);

  void fillQ(unsigned int iQ, double x, double Q, double w);
  void fillQ(unsigned int iQ, double x, double y, double Q, double w);
  void fillQ(unsigned int iQ, double x, double y, double z, double Q, double w);

  double getBinSumW(int ix, int iy=0, int iz=0) const;
  double getBinSumWsq(int ix, int iy=0, int iz=0) const;
  double getBinAvgQ(unsigned int iQ, int ix, int iy=0, int iz=0) const;
  double getBinSigmaQ(unsigned int iQ, int ix, int iy=0, int iz=0) const;

  std::vector<std::vector<std::vector<double>>>& getSumWContainer(){ return sumW; }
  std::vector<std::vector<std::vector<double>>>& getSumWsqContainer(){ return sumWsq; }
  std::vector<std::vector<std::vector<double>>>& getSumWQContainer(unsigned int iQ);
  std::vector<std::vector<std::vector<double>>>& getSumWQsqContainer(unsigned int iQ);

  std::vector<std::vector<std::vector<double>>> const& getSumWContainer() const{ return sumW; }
  std::vector<std::vector<std::vector<double>>> const& getSumWsqContainer() const{ return sumWsq; }
  std::vector<std::vector<std::vector<double>>> const& getSumWQContainer(unsigned int iQ) const;
  std::vector<std::vector<std::vector<double>>> const& getSumWQsqContainer(unsigned int iQ) const;


};


#endif
