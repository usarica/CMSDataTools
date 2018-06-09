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
    tMINLODn, tMINLOUp,
    tQQBkgEWCorrDn, tQQBkgEWCorrUp,
    eLepSFEleDn, eLepSFEleUp,
    eLepSFMuDn, eLepSFMuUp,
    eLepScaleEleDn, eLepScaleEleUp,
    eLepScaleMuDn, eLepScaleMuUp,
    eLepResEleDn, eLepResEleUp,
    eLepResMuDn, eLepResMuUp,
    eJECDn, eJECUp,
    eBTagSFDn, eBTagSFUp,
    eZXStatsDn, eZXStatsUp,
    nSystematicVariations
  };

}

#endif
