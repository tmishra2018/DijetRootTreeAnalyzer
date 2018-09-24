#ifndef analysisClass_h
#define analysisClass_h

#include "baseClass.h"
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <memory>

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"

#include <boost/shared_ptr.hpp>
// typedefs
typedef boost::shared_ptr<fastjet::ClusterSequence>  ClusterSequencePtr;
typedef boost::shared_ptr<fastjet::JetDefinition>    JetDefPtr;
// For JECs
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

//For JER
#include "JetMETCorrections/Modules/interface/JetResolution.h"

using namespace std;

class analysisClass : public baseClass {
public :
  analysisClass(string * inputList, string * cutFile, string * treeName,  string *outputFileName=0, string * cutEfficFile=0);
  virtual ~analysisClass();
  void Loop();
  void correctMETWithTypeI();
private :
  ClusterSequencePtr  fjClusterSeq, fjClusterSeq_shift;
  JetDefPtr           fjJetDefinition;
  // For JECs
  JetCorrectorParameters *L1Par;
  JetCorrectorParameters *L2Par;
  JetCorrectorParameters *L3Par;
  JetCorrectorParameters *L1DATAPar;
  JetCorrectorParameters *L2DATAPar;
  JetCorrectorParameters *L3DATAPar;
  JetCorrectorParameters *L2L3Residual;
  JetCorrectorParameters *L1DATAHLTPar;
  JetCorrectorParameters *L2DATAHLTPar;
  JetCorrectorParameters *L3DATAHLTPar;
  JetCorrectorParameters *L2L3ResidualHLT;
  FactorizedJetCorrector *JetCorrector;
  FactorizedJetCorrector *JetCorrector_data;
  FactorizedJetCorrector *JetCorrector_dataHLT;
  
  FactorizedJetCorrector *JetCorrectortypI;
  FactorizedJetCorrector *JetCorrectortypIL123;
  FactorizedJetCorrector *JetCorrectortypIMC;
  FactorizedJetCorrector *JetCorrectortypIL123MC;
  JetCorrectionUncertainty *unc;
  JetCorrectionUncertainty *unc_TI;
  JetCorrectionUncertainty *unc_RC;
  JetCorrectionUncertainty *unc_MC;
  JetCorrectionUncertainty *unc_MC_TI;
  JetCorrectionUncertainty *unc_MC_RC;
  JetCorrectorParameters *L1JetParForTypeI;
  JetCorrectorParameters *L1JetParForTypeIMC;
  
  JME::JetResolution *resolution;
  JME::JetResolutionScaleFactor *resolution_sf;
  
  JME::JetResolution *resolution_TI;
  JME::JetResolutionScaleFactor *resolution_sf_TI;
  
  
};

#endif

#ifdef analysisClass_cxx

#endif // #ifdef analysisClass_cxx
