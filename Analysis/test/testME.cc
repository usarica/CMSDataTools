#include "common_includes.h"
#include "TemplatesEventAnalyzer.h"
#include "CheckSetTemplatesCategoryScheme.h"
#include "fixTreeWeights.h"


// Process handle
typedef GGProcessHandler ProcessHandleType;
const ProcessHandleType& theProcess = TemplateHelpers::OffshellGGProcessHandle;

void testME(TString generator="MCFM"){
  const ACHypothesis hypo=kSM;
  const std::vector<ProcessHandleType::HypothesisType> tplset = theProcess.getHypothesesForACHypothesis(hypo);
  std::vector<TString> melawgtvars; for (auto& hypotype:tplset) melawgtvars.push_back(theProcess.getMELAHypothesisWeight(hypotype, hypo));
  for (TString const& wgtvar:melawgtvars) MELAout << "Will test " << wgtvar << endl;

  Mela mela(theSqrts, 125., TVar::DEBUG_VERBOSE);

  vector<CJLSTSet*> theSets;
  if (generator=="MCFM"){
    const unsigned int nMCFMChannels=6;
    TString strMCFMChannels[nMCFMChannels]={
      "4mu","4e","2e2mu",
      "2e2tau","2mu2tau","4tau"
    };
    vector<TString> strSamples[nMCFMChannels];
    for (unsigned int ich=0; ich<nMCFMChannels; ich++){
      vector<TString> strSampleIdentifiers;
      strSampleIdentifiers.push_back(Form("gg_MCFM_%s", strMCFMChannels[ich].Data()));
      getSamplesList(theSqrts, strSampleIdentifiers, strSamples[ich]);
    }

    for (unsigned int ich=0; ich<nMCFMChannels; ich++){
      CJLSTSet* theSampleSet = new CJLSTSet(strSamples[ich]);

      for (auto& tree:theSampleSet->getCJLSTTreeList()){
        // Book common variables needed for analysis
        tree->bookBranch<vector<float>*>("LHEMotherPz", nullptr);
        tree->bookBranch<vector<float>*>("LHEMotherE", nullptr);
        tree->bookBranch<vector<short>*>("LHEMotherId", nullptr);
        tree->bookBranch<vector<float>*>("LHEDaughterPt", nullptr);
        tree->bookBranch<vector<float>*>("LHEDaughterEta", nullptr);
        tree->bookBranch<vector<float>*>("LHEDaughterPhi", nullptr);
        tree->bookBranch<vector<float>*>("LHEDaughterMass", nullptr);
        tree->bookBranch<vector<short>*>("LHEDaughterId", nullptr);
        tree->bookBranch<vector<float>*>("LHEAssociatedParticlePt", nullptr);
        tree->bookBranch<vector<float>*>("LHEAssociatedParticleEta", nullptr);
        tree->bookBranch<vector<float>*>("LHEAssociatedParticlePhi", nullptr);
        tree->bookBranch<vector<float>*>("LHEAssociatedParticleMass", nullptr);
        tree->bookBranch<vector<short>*>("LHEAssociatedParticleId", nullptr);

        // Variables for MELA reweighting
        for (TString const& wgtvar:melawgtvars) tree->bookBranch<float>(wgtvar, 0);

        tree->silenceUnused(); // Will no longer book another branch
      }
      theSets.push_back(theSampleSet);
    }
  }
  else{
    const unsigned int nPOWHEGChannels=1;
    vector<TString> strSamples[nPOWHEGChannels];
    for (unsigned int ich=0; ich<nPOWHEGChannels; ich++){
      vector<TString> strSampleIdentifiers;
      strSampleIdentifiers.push_back("gg_Sig_POWHEG");
      getSamplesList(theSqrts, strSampleIdentifiers, strSamples[ich]);
    }

    for (unsigned int ich=0; ich<nPOWHEGChannels; ich++){
      CJLSTSet* theSampleSet = new CJLSTSet(strSamples[ich]);
      for (auto& tree:theSampleSet->getCJLSTTreeList()){
        // Book common variables needed for analysis
        tree->bookBranch<vector<float>*>("LHEMotherPz", nullptr);
        tree->bookBranch<vector<float>*>("LHEMotherE", nullptr);
        tree->bookBranch<vector<short>*>("LHEMotherId", nullptr);
        tree->bookBranch<vector<float>*>("LHEDaughterPt", nullptr);
        tree->bookBranch<vector<float>*>("LHEDaughterEta", nullptr);
        tree->bookBranch<vector<float>*>("LHEDaughterPhi", nullptr);
        tree->bookBranch<vector<float>*>("LHEDaughterMass", nullptr);
        tree->bookBranch<vector<short>*>("LHEDaughterId", nullptr);
        tree->bookBranch<vector<float>*>("LHEAssociatedParticlePt", nullptr);
        tree->bookBranch<vector<float>*>("LHEAssociatedParticleEta", nullptr);
        tree->bookBranch<vector<float>*>("LHEAssociatedParticlePhi", nullptr);
        tree->bookBranch<vector<float>*>("LHEAssociatedParticleMass", nullptr);
        tree->bookBranch<vector<short>*>("LHEAssociatedParticleId", nullptr);

        // Variables for MELA reweighting
        for (TString const& wgtvar:melawgtvars) tree->bookBranch<float>(wgtvar, 0);

        tree->silenceUnused(); // Will no longer book another branch
      }
      theSets.push_back(theSampleSet);
    }
  }
  bool printMomenta=false;
  for (auto& theSampleSet:theSets){
    for (auto& tree:theSampleSet->getCJLSTTreeList()){
      MELAout << "Looping over " << tree->sampleIdentifier << endl;
      int ev=0;
      while (tree->getEvent(ev)){
        if (ev%100000==0) MELAout << " - Event " << ev << endl;
        for (TString const& wgtvar:melawgtvars){
          float w=0;
          tree->getVal(wgtvar, w);
          if (w<=0.){
            MELAout << "Weight " << wgtvar << " = " << w << " at event " << ev << endl;
            printMomenta=true;
          }
        }
        if (printMomenta){
          std::vector<float> const* LHEMotherPz;
          std::vector<float> const* LHEMotherE;
          std::vector<short> const* LHEMotherId;
          std::vector<float> const* LHEDaughterPt;
          std::vector<float> const* LHEDaughterEta;
          std::vector<float> const* LHEDaughterPhi;
          std::vector<float> const* LHEDaughterMass;
          std::vector<short> const* LHEDaughterId;
          std::vector<float> const* LHEAssociatedParticlePt;
          std::vector<float> const* LHEAssociatedParticleEta;
          std::vector<float> const* LHEAssociatedParticlePhi;
          std::vector<float> const* LHEAssociatedParticleMass;
          std::vector<short> const* LHEAssociatedParticleId;
          tree->getVal("LHEMotherPz", LHEMotherPz);
          tree->getVal("LHEMotherE", LHEMotherE);
          tree->getVal("LHEMotherId", LHEMotherId);
          tree->getVal("LHEDaughterPt", LHEDaughterPt);
          tree->getVal("LHEDaughterEta", LHEDaughterEta);
          tree->getVal("LHEDaughterPhi", LHEDaughterPhi);
          tree->getVal("LHEDaughterMass", LHEDaughterMass);
          tree->getVal("LHEDaughterId", LHEDaughterId);
          tree->getVal("LHEAssociatedParticlePt", LHEAssociatedParticlePt);
          tree->getVal("LHEAssociatedParticleEta", LHEAssociatedParticleEta);
          tree->getVal("LHEAssociatedParticlePhi", LHEAssociatedParticlePhi);
          tree->getVal("LHEAssociatedParticleMass", LHEAssociatedParticleMass);
          tree->getVal("LHEAssociatedParticleId", LHEAssociatedParticleId);

          SimpleParticleCollection_t mothers, daughters, aparticles;
          for (unsigned int ipart=0; ipart<LHEMotherId->size(); ipart++){
            TLorentzVector p4(0, 0, LHEMotherPz->at(ipart), LHEMotherE->at(ipart));
            MELAParticle part(LHEMotherId->at(ipart), p4);
            MELAout << "Mother " << ipart << ": " << part << endl;
            mothers.push_back(SimpleParticle_t((part.id!=21 ? 0 : part.id), part.p4));
          }
          for (unsigned int ipart=0; ipart<LHEDaughterId->size(); ipart++){
            TLorentzVector p4; p4.SetPtEtaPhiM(LHEDaughterPt->at(ipart), LHEDaughterEta->at(ipart), LHEDaughterPhi->at(ipart), LHEDaughterMass->at(ipart));
            MELAParticle part(LHEDaughterId->at(ipart), p4);
            MELAout << "Daughter " << ipart << ": " << part << endl;
            daughters.push_back(SimpleParticle_t(part.id, part.p4));
          }
          for (unsigned int ipart=0; ipart<LHEAssociatedParticleId->size(); ipart++){
            TLorentzVector p4; p4.SetPtEtaPhiM(LHEAssociatedParticlePt->at(ipart), LHEAssociatedParticleEta->at(ipart), LHEAssociatedParticlePhi->at(ipart), LHEAssociatedParticleMass->at(ipart));
            MELAParticle part(LHEAssociatedParticleId->at(ipart), p4);
            MELAout << "AssociatedParticle " << ipart << ": " << part << endl;
            aparticles.push_back(SimpleParticle_t(part.id, part.p4));
          }
          
          float pME;
          mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
          mela.setInputEvent(&daughters, &aparticles, &mothers, true);

          pME=-1;
          mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
          mela.computeP(pME, false);
          MELAout << "p_ggZZ_sig: " << pME << endl;

          pME=-1;
          mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
          mela.computeP(pME, false);
          MELAout << "p_ggZZ_bkg: " << pME << endl;

          mela.resetInputEvent();

          break;
        }
        ev++;
      }
      if (printMomenta) break;
    }
    if (printMomenta) break;
  }

  for (auto& theSampleSet:theSets) delete theSampleSet; theSets.clear();
  MELAout.close();
}

