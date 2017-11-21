#include "Discriminant.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


Discriminant::Discriminant(const TString cfilename, const TString splinename) : theCFile(nullptr), theC(nullptr){
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
        MELAerr << "Discriminant::Discriminant: File is zombie!" << endl;
        theCFile->Close();
        theCFile=nullptr;
      }
    }
    else MELAerr << "Discriminant::Discriminant: File could not be opened!" << endl;
  }
  else MELAout << "Discriminant::Discriminant: No c-constants file is specified, defaulting to c=1." << endl;
  val=-1;
}
Discriminant::~Discriminant(){ if (theCFile) theCFile->Close(); }

Discriminant::operator float() const{ return val; }
Discriminant::operator float&(){ return val; }
Discriminant::operator float*(){ return &val; }

bool Discriminant::operator<(const float& other) const{ return val<other; }
bool Discriminant::operator>(const float& other) const{ return val>other; }
bool Discriminant::operator<=(const float& other) const{ return val<=other; }
bool Discriminant::operator>=(const float& other) const{ return val>=other; }
bool Discriminant::operator==(const float& other) const{ return val==other; }
bool Discriminant::operator!=(const float& other) const{ return val!=other; }

float Discriminant::update(const std::vector<float>& vars, const float valReco){
  this->eval(vars, valReco);
  return val;
}
float Discriminant::getCval(const float valReco) const{
  if (theC) return theC->Eval(valReco);
  else return 1;
}
float Discriminant::applyAdditionalC(const float cval){ val = val/(val+(1.-val)*cval); return val; }

