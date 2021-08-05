#ifndef CMSLORENTZVECTOR_H
#define CMSLORENTZVECTOR_H

#ifdef _COMPILE_STANDALONE_

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > CMSLorentzVector_f;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > CMSLorentzVector_d;

#else

#include <DataFormats/Math/interface/LorentzVector.h>

typedef math::XYZTLorentzVectorF CMSLorentzVector_f;
typedef math::XYZTLorentzVectorD CMSLorentzVector_d;

#endif

typedef CMSLorentzVector_f CMSLorentzVector;

#endif
