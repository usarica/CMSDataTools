#ifndef SYSTEMATICVARIATIONS_H
#define SYSTEMATICVARIATIONS_H


namespace SystematicsHelpers{

  enum SystematicVariationTypes{
    sNominal,
    tPDFScaleDn, tPDFScaleUp,
    tQCDScaleDn, tQCDScaleUp,
    tAsMZDn, tAsMZUp,
    tPDFReplicaDn, tPDFReplicaUp,
    tPythiaScaleDn, tPythiaScaleUp,
    tPythiaTuneDn, tPythiaTuneUp,
    tQQBkgEWCorrDn, tQQBkgEWCorrUp,
    eLepSFEleDn, eLepSFEleUp,
    eLepSFMuDn, eLepSFMuUp,
    eJECDn, eJECUp,
    eZXStatsDn, eZXStatsUp,
    nSystematicVariations
  };

}

#endif
