#include "Discriminant.h"


using namespace std;


Discriminant::Discriminant(const TString cfilename, const TString splinename){
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
      theC=nullptr;
    }
  }
  else{
    theCFile=nullptr;
    theC=nullptr;
  }
  val=-1;
}
Discriminant::~Discriminant(){ if (theCFile) theCFile->Close(); }
Discriminant::operator float() const{ return val; }
float Discriminant::update(const std::vector<float>& vars, const float& valReco){
  this->eval(vars, valReco);
  return val;
}
float Discriminant::getCval(const float valReco) const{
  if (theC) return theC->Eval(valReco);
  else return 1;
}
