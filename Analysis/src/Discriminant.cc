#include "Discriminant.h"


using namespace std;


Discriminant::Discriminant(const TString cfilename, const TString splinename) : theCFile(nullptr), theC(nullptr){
  if (cfilename!=""){
    theCFile = TFile::Open(cfilename);
    if (theCFile){
      if (theCFile->IsOpen() && !theCFile->IsZombie()){
        theC = (TSpline3*) theCFile->Get(splinename);
        if (!theC){
          theC=nullptr;
          theCFile->Close();
          theCFile=nullptr;
        }
      }
      else if (theCFile->IsOpen()){
        theCFile->Close();
        theCFile=nullptr;
      }
    }
  }
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

