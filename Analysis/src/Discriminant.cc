#include "Discriminant.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


Discriminant::Discriminant(
  const TString cfilename, const TString splinename,
  const TString gfilename, const TString gsplinename,
  const float gscale_
) :
  theCFile(nullptr), theC(nullptr), WPCshift(1),
  theGFile(nullptr), theG(nullptr), gscale(gscale_), invertG(false)
{
  if (cfilename!=""){
    MELAout << "Discriminant::Discriminant: Opening " << cfilename << endl;
    theCFile = TFile::Open(cfilename);
    if (theCFile){
      if (theCFile->IsOpen() && !theCFile->IsZombie()){
        theC = (TSpline3*) theCFile->Get(splinename);
        if (!theC){
          MELAerr << "Discriminant::Discriminant: Spline " << splinename << " does not exist!" << endl;
          theC=nullptr;
          theCFile->Close();
          theCFile=nullptr;
        }
        else MELAout << "Discriminant::Discriminant: Acquired " << splinename << endl;
      }
      else if (theCFile->IsOpen()){
        MELAerr << "Discriminant::Discriminant: File " << splinename << " is zombie!" << endl;
        theCFile->Close();
        theCFile=nullptr;
      }
    }
    else MELAerr << "Discriminant::Discriminant: File " << splinename << " could not be opened!" << endl;
  }
  else MELAout << "Discriminant::Discriminant: No c-constants file is specified, defaulting to c=1." << endl;
  if (gfilename!=""){
    MELAout << "Discriminant::Discriminant: Opening " << gfilename << endl;
    theGFile = TFile::Open(gfilename);
    if (theGFile){
      if (theGFile->IsOpen() && !theGFile->IsZombie()){
        theG = (TSpline3*) theGFile->Get(gsplinename);
        if (!theG){
          MELAerr << "Discriminant::Discriminant: Spline " << gsplinename << " does not exist!" << endl;
          theG=nullptr;
          theGFile->Close();
          theGFile=nullptr;
        }
        else MELAout << "Discriminant::Discriminant: Acquired " << gsplinename << endl;
      }
      else if (theGFile->IsOpen()){
        MELAerr << "Discriminant::Discriminant: File " << gsplinename << " is zombie!" << endl;
        theGFile->Close();
        theGFile=nullptr;
      }
    }
    else MELAerr << "Discriminant::Discriminant: File " << gsplinename << " could not be opened!" << endl;
  }
  else MELAout << "Discriminant::Discriminant: No g-constants file is specified, defaulting to g=1." << endl;
  val=-1;
}
Discriminant::~Discriminant(){
  if (theCFile) theCFile->Close();
  if (theGFile) theGFile->Close();
}

Discriminant::operator float() const{ return val; }
Discriminant::operator float&(){ return val; }
Discriminant::operator float*(){ return &val; }

bool Discriminant::operator<(const float& other) const{ return val<other; }
bool Discriminant::operator>(const float& other) const{ return val>other; }
bool Discriminant::operator<=(const float& other) const{ return val<=other; }
bool Discriminant::operator>=(const float& other) const{ return val>=other; }
bool Discriminant::operator==(const float& other) const{ return val==other; }
bool Discriminant::operator!=(const float& other) const{ return val!=other; }

void Discriminant::resetVal(){ val=-999; }
float Discriminant::update(const std::vector<float>& vars, const float valReco){
  this->eval(vars, valReco);
  return val;
}
float Discriminant::getCval(const float valReco) const{
  float res=WPCshift;
  if (theC) res*=theC->Eval(valReco);
  if (theG) res*=pow(theG->Eval(valReco)*gscale, (!invertG ? 1 : -1)*2);
  return res;
}
float Discriminant::applyAdditionalC(const float cval){ val = val/(val+(1.-val)*cval); return val; }
void Discriminant::setWP(float WPval){
  if (WPval<=0. || WPval>=1.) return;
  WPCshift = WPval/(1.-WPval);
}
void Discriminant::setInvertG(bool flag){ invertG=flag; }
