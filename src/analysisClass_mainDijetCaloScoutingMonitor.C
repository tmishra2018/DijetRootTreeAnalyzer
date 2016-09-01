#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>

bool verbose;

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile)
{
  std::cout << "analysisClass::analysisClass(): begins " << std::endl;

  std::string jetAlgo = getPreCutString1("jetAlgo");
  double rParam = getPreCutValue1("DeltaR");
  verbose = getPreCutValue1("verbose");

  if( jetAlgo == "AntiKt" )
    fjJetDefinition = JetDefPtr( new fastjet::JetDefinition(fastjet::antikt_algorithm, rParam) );
  else if( jetAlgo == "Kt" )
    fjJetDefinition = JetDefPtr( new fastjet::JetDefinition(fastjet::kt_algorithm, rParam) );
  else 
    fjJetDefinition = JetDefPtr( new fastjet::JetDefinition(fastjet::cambridge_algorithm, rParam) );

  // For JECs
  if( int(getPreCutValue1("useJECs"))==1 )
  {
    std::cout << "Reapplying JECs on the fly" << std::endl;

    /* JEC v6
    std::string L1Path = "data/Summer15_25nsV6_MC/Summer15_25nsV6_MC_L1FastJet_AK4PFchs.txt";
    std::string L2Path = "data/Summer15_25nsV6_MC/Summer15_25nsV6_MC_L2Relative_AK4PFchs.txt"; 
    std::string L3Path = "data/Summer15_25nsV6_MC/Summer15_25nsV6_MC_L3Absolute_AK4PFchs.txt";

    std::string L1DATAPath = "data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_L1FastJet_AK4PFchs.txt";
    std::string L2DATAPath = "data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_L2Relative_AK4PFchs.txt"; 
    std::string L3DATAPath = "data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_L3Absolute_AK4PFchs.txt";
    std::string L2L3ResidualPath = "data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_L2L3Residual_AK4PFchs.txt" ;

    std::string L1DATAHLTPath = "data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_L1FastJet_AK4PFchs.txt";
    std::string L2DATAHLTPath = "data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_L2Relative_AK4PFchs.txt"; 
    std::string L3DATAHLTPath = "data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_L3Absolute_AK4PFchs.txt";
    std::string L2L3ResidualHLTPath = "data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_L2L3Residual_AK4PFchs.txt" ;
    */

    /*
    // HLT + v6Res
    std::string L1Path = "data/Summer15_25nsV6_MC/Summer15_25nsV6_MC_L1FastJet_AK4PFchs.txt";
    std::string L2Path = "data/Summer15_25nsV6_MC/Summer15_25nsV6_MC_L2Relative_AK4PFchs.txt"; 
    std::string L3Path = "data/Summer15_25nsV6_MC/Summer15_25nsV6_MC_L3Absolute_AK4PFchs.txt";

    std::string L1DATAPath = "data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_L1FastJet_AK4PFchs.txt";
    std::string L2DATAPath = "data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_L2Relative_AK4PFchs.txt"; 
    std::string L3DATAPath = "data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_L3Absolute_AK4PFchs.txt";
    std::string L2L3ResidualPath = "data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_L2L3Residual_AK4PFchs.txt" ;

    std::string L1DATAHLTPath = "data/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_MC/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_L1FastJet_AK4PFchs.txt";
    std::string L2DATAHLTPath = "data/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_MC/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_L2Relative_AK4PFchs.txt"; 
    std::string L3DATAHLTPath = "data/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_MC/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_L3Absolute_AK4PFchs.txt";
    std::string L2L3ResidualHLTPath = "data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_L2L3Residual_AK4PFchs.txt" ;
    */

    // HLT + v7Res
    // MC non used
    std::string L1Path = "data/Summer15_25nsV7_MC/Summer15_25nsV7_MC_L1FastJet_AK4PFchs.txt";
    std::string L2Path = "data/Summer15_25nsV7_MC/Summer15_25nsV7_MC_L2Relative_AK4PFchs.txt"; 
    std::string L3Path = "data/Summer15_25nsV7_MC/Summer15_25nsV7_MC_L3Absolute_AK4PFchs.txt";
    // data: PFReco
    std::string L1DATAPath = "data/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L1FastJet_AK4PFchs.txt";
    std::string L2DATAPath = "data/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L2Relative_AK4PFchs.txt"; 
    std::string L3DATAPath = "data/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L3Absolute_AK4PFchs.txt";
    std::string L2L3ResidualPath = "data/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L2L3Residual_AK4PFchs.txt" ;
    // data: CaloHLT
    std::string L1DATAHLTPath = "data/80X_dataRun2_HLT_frozen_v12/80X_dataRun2_HLT_frozen_v12_L1FastJet_AK4CaloHLT.txt";
    std::string L2DATAHLTPath = "data/80X_dataRun2_HLT_frozen_v12/80X_dataRun2_HLT_frozen_v12_L2Relative_AK4CaloHLT.txt"; 
    std::string L3DATAHLTPath = "data/80X_dataRun2_HLT_frozen_v12/80X_dataRun2_HLT_frozen_v12_L3Absolute_AK4CaloHLT.txt";
    std::string L2L3ResidualHLTPath = "data/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L2L3Residual_AK4PFchs.txt";

    L1Par = new JetCorrectorParameters(L1Path);
    L2Par = new JetCorrectorParameters(L2Path);
    L3Par = new JetCorrectorParameters(L3Path);

    L1DATAPar = new JetCorrectorParameters(L1DATAPath);
    L2DATAPar = new JetCorrectorParameters(L2DATAPath);
    L3DATAPar = new JetCorrectorParameters(L3DATAPath);
    L2L3Residual = new JetCorrectorParameters(L2L3ResidualPath);

    L1DATAHLTPar = new JetCorrectorParameters(L1DATAHLTPath);
    L2DATAHLTPar = new JetCorrectorParameters(L2DATAHLTPath);
    L3DATAHLTPar = new JetCorrectorParameters(L3DATAHLTPath);
    L2L3ResidualHLT = new JetCorrectorParameters(L2L3ResidualHLTPath);

    std::vector<JetCorrectorParameters> vPar;
    std::vector<JetCorrectorParameters> vPar_data;
    std::vector<JetCorrectorParameters> vPar_dataHLT;
    vPar.push_back(*L1Par);
    vPar.push_back(*L2Par);
    vPar.push_back(*L3Par);
    //residuals are applied only to data
    vPar_data.push_back(*L1DATAPar);
    vPar_data.push_back(*L2DATAPar);
    vPar_data.push_back(*L3DATAPar);
    vPar_data.push_back(*L2L3Residual);
    //residuals are applied only to data -- HLT
    vPar_dataHLT.push_back(*L1DATAHLTPar);
    vPar_dataHLT.push_back(*L2DATAHLTPar);
    vPar_dataHLT.push_back(*L3DATAHLTPar);
    vPar_dataHLT.push_back(*L2L3ResidualHLT);

    JetCorrector = new FactorizedJetCorrector(vPar);
    JetCorrector_data = new FactorizedJetCorrector(vPar_data);
    JetCorrector_dataHLT = new FactorizedJetCorrector(vPar_dataHLT);

    //uncertainty
    unc = new JetCorrectionUncertainty("data/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_Uncertainty_AK4PFchs.txt");
  }
  
  std::cout << "analysisClass::analysisClass(): ends " << std::endl;
}

analysisClass::~analysisClass()
{
  std::cout << "analysisClass::~analysisClass(): begins " << std::endl;

  std::cout << "analysisClass::~analysisClass(): ends " << std::endl;
}

void analysisClass::Loop()
{
   std::cout << "analysisClass::Loop() begins" <<std::endl;   
    
   if (fChain == 0) return;
   
   //////////book histos here


   // variable binning for mjj trigger efficiency plots
   const int nMassBins = 103;

   double massBoundaries[nMassBins+1] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325,
     354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687,
     1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509,
     4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 
     10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000};

   // char* HLTname[50] = {"noTrig","PFHT475","PFHT800","PFHT650MJJ900","PFHT800_OR_PFHT650MJJ900","PFHT800_noPFHT475", 
   //                      "Mu45Eta2p1", "PFHT800AndMu45Eta2p1"};
   // TH1F* h_mjj_HLTpass[8];
   // char name_histoHLT[50];
   // for (int i=0; i<8; i++){  
   //   sprintf(name_histoHLT,"h_mjj_HLTpass_%s",HLTname[i]);
   //   h_mjj_HLTpass[i]= new TH1F(name_histoHLT,"",103,massBoundaries);
   // }

   //For trigger efficiency measurements
   //No trigger selection applied (full offline selection applied)
   //   TH1F* h_mjj_NoTrigger = new TH1F("h_mjj_NoTrigger","",103,massBoundaries);
   //   TH1F* h_mjj_HLTpass_CaloJet40 = new TH1F("h_mjj_HLTpass_CaloJet40","",103,massBoundaries);
   //   TH1F* h_mjj_HLTpass_CaloJet40_L1HTT = new TH1F("h_mjj_HLTpass_CaloJet40_L1HTT","",103,massBoundaries);  
   //  TH1F* h_mjj_HLTpass_CaloJet40_L1HTT_CaloHT250 = new TH1F("h_mjj_HLTpass_CaloJet40_L1HTT_CaloHT250","",103,massBoundaries);  


   /////////initialize variables

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) { // running over all events
     //  for (Long64_t jentry=0; jentry<2;jentry++) { //runnig over 2 events for debugging
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     // if (Cut(ientry) < 0) continue;
   
     ////////////////////// User's code starts here ///////////////////////
     
     if(verbose) std::cout<< std::endl;
     if(verbose) std::cout<< "run number: "<< runNo << std::endl;

     ///Stuff to be done for every event

     size_t no_jets_ak4= jetPtAK4calo->size();          // HLT jets
     size_t no_jetsReco_ak4=jetPtAK4reco->size();  // Reco jets

     if(verbose) std::cout << "number of jets Reco = "<< no_jetsReco_ak4 << std::endl;
     if(verbose) std::cout << "number of jets HLT = "<< no_jets_ak4 << std::endl;

     vector<TLorentzVector> widejetsReco;
     TLorentzVector wj1Reco, wj2Reco, wdijetReco; 
     TLorentzVector wj1Reco_shift, wj2Reco_shift, wdijetReco_shift; 

     vector<TLorentzVector> widejets;
     TLorentzVector wj1, wj2, wdijet; 
     TLorentzVector wj1_shift, wj2_shift, wdijet_shift; 

     vector<TLorentzVector> AK4recojets;
     TLorentzVector ak4j1Reco, ak4j2Reco, ak4dijetReco;      
     TLorentzVector ak4j3Reco;
     bool thirdRecoJetPresent = false;

     vector<TLorentzVector> AK4jets;
     TLorentzVector ak4j1, ak4j2, ak4dijet;      
     TLorentzVector ak4j3;
     bool thirdJetPresent = false ;

     resetCuts();

     std::vector<double> jecFactors;
     std::vector<double> jecUncertainty;
     std::vector<bool> idCaloJet;

     std::vector<double> jecFactorsReco;
     std::vector<double> jecUncertaintyReco;
     // new JECs could change the jet pT ordering. the vector below
     // holds sorted jet indices after the new JECs had been applied
     std::vector<unsigned> sortedRecoJetIdx;
     std::vector<unsigned> sortedJetIdx;

     if( int(getPreCutValue1("useJECs"))==1 )
       {
	 if(verbose) std::cout << "USE JECS"<<std::endl;

	 // Reco jets
	 // sort jets by increasing pT
	 std::multimap<double, unsigned> sortedRecoJets;
	 for(size_t j=0; j<no_jetsReco_ak4; ++j)
	   { // Reco Jets
	     //MC
	     JetCorrector->setJetEta(jetEtaAK4reco->at(j));
	     JetCorrector->setJetPt(jetPtAK4reco->at(j)/jetJecAK4reco->at(j)); //pTraw
	     JetCorrector->setJetA(jetAreaAK4reco->at(j));
	     JetCorrector->setRho(rhoreco);
	     //Data
  	     JetCorrector_data->setJetEta(jetEtaAK4reco->at(j));
	     JetCorrector_data->setJetPt(jetPtAK4reco->at(j)/jetJecAK4reco->at(j)); //pTraw
	     JetCorrector_data->setJetA(jetAreaAK4reco->at(j));
	     JetCorrector_data->setRho(rhoreco);

  	     //nominal value of JECs
	     double correction;//, old_correction, nominal_correction;
	     //if( int(getPreCutValue1("shiftJECs"))==0 ){
	     if (isData == 1) correction = JetCorrector_data->getCorrection();
	     else correction = JetCorrector->getCorrection();
	     //nominal_correction=correction;
	     //old_correction = jetJecAK4reco->at(j);
	     //}

	     //JEC uncertainties
	     unc->setJetEta(jetEtaAK4reco->at(j));
	     unc->setJetPt(jetPtAK4reco->at(j)/jetJecAK4reco->at(j)*correction);
	     double uncertainty = unc->getUncertainty(true);
	     jecUncertaintyReco.push_back(uncertainty); 

	     // std::cout << "run:" << runNo << "    lumi:" << lumi << "   event:" << evtNo << "   jet pt:" << jetPtAK4reco->at(j)/jetJecAK4reco->at(j)*correction << "   correction:" << correction <<   "   uncertainty:" <<  uncertainty  << "  nominal correction:" << nominal_correction  << " old correction: " << old_correction << std::endl;
	     //use "shifted" JECs for study of systematic uncertainties 
	     if( int(getPreCutValue1("shiftJECs"))==1 ){
	       //flat shift
	       //if (isData == 1) correction = JetCorrector_data->getCorrection() * getPreCutValue2("shiftJECs");
	       //else correction = JetCorrector->getCorrection() * getPreCutValue2("shiftJECs");
	       //shift of the corresponding unc
	       correction = correction + getPreCutValue2("shiftJECs")*uncertainty*correction;
	       //  std::cout << "run:" << runNo << "    lumi:" << lumi << "   event:" << evtNo << "   jet pt:" << jetPtAK3->at(j)/jetJecAK4reco->at(j)*correction << "   correction:" << correction << "   uncertainty:" <<  uncertainty  << std::endl << std::endl;
	       
	     }

	     jecFactorsReco.push_back(correction);
	     sortedRecoJets.insert(std::make_pair((jetPtAK4reco->at(j)/jetJecAK4reco->at(j))*correction, j));

	   } // end loop reco jets

	 // get jet indices in decreasing pT order
	 for(std::multimap<double, unsigned>::const_reverse_iterator it = sortedRecoJets.rbegin(); it != sortedRecoJets.rend(); ++it)
	   sortedRecoJetIdx.push_back(it->second);
	 
	 //HLT jets
	 // sort jets by increasing pT
	 std::multimap<double, unsigned> sortedJets;
	 for(size_t j=0; j<no_jets_ak4; ++j)
	   { // HLT Jets
	     JetCorrector->setJetEta(jetEtaAK4calo->at(j));
	     JetCorrector->setJetPt(jetPtAK4calo->at(j)/jetJecAK4calo->at(j)); //pTraw
	     JetCorrector->setJetA(jetAreaAK4calo->at(j));
	     JetCorrector->setRho(rhocalo);
	     
	     JetCorrector_dataHLT->setJetEta(jetEtaAK4calo->at(j));
	     JetCorrector_dataHLT->setJetPt(jetPtAK4calo->at(j)/jetJecAK4calo->at(j)); //pTraw
	     JetCorrector_dataHLT->setJetA(jetAreaAK4calo->at(j));
	     JetCorrector_dataHLT->setRho(rhocalo);
	     
	     //nominal value of JECs
	     double correction;//, old_correction, nominal_correction;
	     //if( int(getPreCutValue1("shiftJECs"))==0 ){
	     if (isData == 1) correction = JetCorrector_dataHLT->getCorrection();
	     else correction = JetCorrector->getCorrection();
	     //to work with raw jet
	     //	   if (isData == 1) correction = 1;
	     //          else correction = 1;
	     //nominal_correction=correction;
	     //old_correction = jetJecAK4calo->at(j);
	     //}
	     
	     //JEC uncertainties
	     unc->setJetEta(jetEtaAK4calo->at(j));
	     unc->setJetPt(jetPtAK4calo->at(j)/jetJecAK4calo->at(j)*correction);
	     double uncertainty = unc->getUncertainty(true);
	     jecUncertainty.push_back(uncertainty); 
	     
	     // std::cout << "run:" << runNo << "    lumi:" << lumi << "   event:" << evtNo << "   jet pt:" << jetPtAK4calo->at(j)/jetJecAK4calo->at(j)*correction << "   correction:" << correction <<   "   uncertainty:" <<  uncertainty  << "  nominal correction:" << nominal_correction  << " old correction: " << old_correction << std::endl;
	     //use "shifted" JECs for study of systematic uncertainties 
	     if( int(getPreCutValue1("shiftJECs"))==1 ){
	       //flat shift
	       //if (isData == 1) correction = JetCorrector_data->getCorrection() * getPreCutValue2("shiftJECs");
	       //else correction = JetCorrector->getCorrection() * getPreCutValue2("shiftJECs");
	       //shift of the corresponding unc
	       correction = correction + getPreCutValue2("shiftJECs")*uncertainty*correction;
	       //  std::cout << "run:" << runNo << "    lumi:" << lumi << "   event:" << evtNo << "   jet pt:" << jetPtAK3->at(j)/jetJecAK4calo->at(j)*correction << "   correction:" << correction << "   uncertainty:" <<  uncertainty  << std::endl << std::endl;
	       
	     }
	     
	     jecFactors.push_back(correction);

	     bool idval = false;
	     if ( jetHadfAK4calo->at(j) <  getPreCutValue1("hadFraction") && jetEmfAK4calo->at(j) <  getPreCutValue1("emFraction")  ) 
	       idval = true;
	     
	     idCaloJet.push_back(idval);
	     
	     sortedJets.insert(std::make_pair((jetPtAK4calo->at(j)/jetJecAK4calo->at(j))*correction, j));
	     
	   }// end loop HLT jets
	 
	 // get jet indices in decreasing pT order
	 for(std::multimap<double, unsigned>::const_reverse_iterator it = sortedJets.rbegin(); it != sortedJets.rend(); ++it)
	   sortedJetIdx.push_back(it->second);
	 
       }// end use JEC
     else if( int(getPreCutValue1("noJECs"))==1  )
       {
	 if(verbose) std::cout << "No JECS"<<std::endl; // raw jets
	 
	 // sort jets by increasing pT
	 std::multimap<double, unsigned> sortedJets;
	 std::multimap<double, unsigned> sortedRecoJets;
	 for(size_t j=0; j<no_jetsReco_ak4; ++j) //same ordering of original root trees
	   {
	     jecUncertaintyReco.push_back(0.); 
	     jecFactorsReco.push_back(1.);
	     sortedRecoJets.insert(std::make_pair((jetPtAK4reco->at(j)/jetJecAK4reco->at(j)), j)); //raw
	   }       
	 // get jet indices in decreasing pT order
	 for(std::multimap<double, unsigned>::const_reverse_iterator it = sortedRecoJets.rbegin(); it != sortedRecoJets.rend(); ++it)
	   sortedRecoJetIdx.push_back(it->second);
	 // ----------------- fatto Reco  ---- faccio HLT
	 for(size_t j=0; j<no_jets_ak4; ++j) //same ordering of original root trees
	   {
	     jecUncertainty.push_back(0.); 
	     jecFactors.push_back(1.);

	     bool idval = false;
	     if ( jetHadfAK4calo->at(j) <  getPreCutValue1("hadFraction") && jetEmfAK4calo->at(j) <  getPreCutValue1("emFraction")  ) 
	       idval = true;
	     
	     idCaloJet.push_back(idval);
	     
	     sortedJets.insert(std::make_pair((jetPtAK4calo->at(j)/jetJecAK4calo->at(j)), j)); //raw
	   }       
	 // get jet indices in decreasing pT order
	 for(std::multimap<double, unsigned>::const_reverse_iterator it = sortedJets.rbegin(); it != sortedJets.rend(); ++it)
	   sortedJetIdx.push_back(it->second);
       }
     else
       {
	 if(verbose) std::cout << "Third case"<< std::endl;
	 // don't remove the old corrections	 
	 for(size_t j=0; j<no_jetsReco_ak4; ++j) //same ordering of original root trees
	   {
	     jecFactorsReco.push_back(jetJecAK4reco->at(j));
	     jecUncertaintyReco.push_back(0.); 
	     sortedRecoJetIdx.push_back(j);
	   }
	 // ----------------- fatto Reco  ---- faccio HLT
	 for(size_t j=0; j<no_jets_ak4; ++j) //same ordering of original root trees
	   {
	     jecFactors.push_back(jetJecAK4calo->at(j));
	     jecUncertainty.push_back(0.); 

	     bool idval = false;
	     if ( jetHadfAK4calo->at(j) <  getPreCutValue1("hadFraction") && jetEmfAK4calo->at(j) <  getPreCutValue1("emFraction")  ) 
	       idval = true;
	     
	     idCaloJet.push_back(idval);
	     
	     sortedJetIdx.push_back(j);
	   }
       }


     //#############################################################
     //########## NOTE: from now on sortedJetIdx[ijet] should be used
     //#############################################################

     // if(no_jets_ak4>=2){
     //  if(!(fabs(jetEtaAK4reco->at(0)) < getPreCutValue1("jetFidRegion") && idTAK4reco->at(0) == getPreCutValue1("tightJetID"))){
     //    std::cout << " JET 0 FAIL " << jetEtaAK4reco->at(0) << " JET 0  ID " << idTAK4reco->at(0) << std::endl;
     //  }
     //  if(!(fabs(jetEtaAK4reco->at(1)) < getPreCutValue1("jetFidRegion") && idTAK4reco->at(1) == getPreCutValue1("tightJetID"))){
     //    std::cout << " JET 1 FAIL " << jetEtaAK4reco->at(1) << " JET 1  ID " << idTAK4reco->at(1) << std::endl;
     //  }  
     // }

     //count (Reco/HLT) ak4 jets passing pt threshold and id criteria
     int Nak4Reco = 0;
     double HTak4Reco = 0;
     
     int Nak4 = 0;
     double HTak4 = 0;
     
     if(verbose) std::cout << "RECO PF JETS"<< std::endl;
     for(size_t ijet=0; ijet<no_jetsReco_ak4; ++ijet)
       {       // Reco Jets
	 if(verbose){
	   cout << "ijet=" << ijet << " , sortedRecoJetIdx[ijet]=" << sortedRecoJetIdx[ijet] 
		<< ", raw pT=" << jetPtAK4reco->at(sortedRecoJetIdx[ijet])/jetJecAK4reco->at(sortedRecoJetIdx[ijet]) 
		<< ", corrected pT-old=" << jetPtAK4reco->at(sortedRecoJetIdx[ijet] )
		<< ", jecFactors-old=" << jetJecAK4reco->at(sortedRecoJetIdx[ijet])
		<< ", jecFactors=" << jecFactorsReco[sortedRecoJetIdx[ijet]]
		<< ", final corrected pT =" << (jecFactorsReco[sortedRecoJetIdx[ijet]]/jetJecAK4reco->at(sortedRecoJetIdx[ijet]))*jetPtAK4reco->at(sortedRecoJetIdx[ijet])
		<< endl;
	 }	 
	 //////////////cout << "id Tight jet" << sortedRecoJetIdx[1] << " = " << idTAK4reco->at(sortedRecoJetIdx[1]) << endl;
	 if(fabs(jetEtaAK4reco->at(sortedRecoJetIdx[ijet])) < getPreCutValue1("jetFidRegion")
	    && idTAK4reco->at(sortedRecoJetIdx[ijet]) == getPreCutValue1("tightJetID")
	    && (jecFactorsReco[sortedRecoJetIdx[ijet]]/jetJecAK4reco->at(sortedRecoJetIdx[ijet]))*jetPtAK4reco->at(sortedRecoJetIdx[ijet]) > getPreCutValue1("pt1Cut"))
	   {
	     Nak4Reco += 1;
	     HTak4Reco += (jecFactorsReco[sortedRecoJetIdx[ijet]]/jetJecAK4reco->at(sortedRecoJetIdx[ijet]))*jetPtAK4reco->at(sortedRecoJetIdx[ijet]);
	   }
       }
     
     if(verbose) std::cout << "HLT Calo JETS"<< std::endl;
     for(size_t ijet=0; ijet<no_jets_ak4; ++ijet)
       {	 // HLT Jets
	 if(verbose){
	   cout << "ijet=" << ijet << " , sortedJetIdx[ijet]=" << sortedJetIdx[ijet] 
		<< ", raw pT=" << jetPtAK4calo->at(sortedJetIdx[ijet])/jetJecAK4calo->at(sortedJetIdx[ijet]) 
		<< ", corrected pT-old=" << jetPtAK4calo->at(sortedJetIdx[ijet] ) 
		<< ", jecFactors-old=" << jetJecAK4calo->at(sortedJetIdx[ijet])
		<< ", jecFactors=" << jecFactors[sortedJetIdx[ijet]]
		<< ", final corrected pT=" << (jecFactors[sortedJetIdx[ijet]]/jetJecAK4calo->at(sortedJetIdx[ijet]))*jetPtAK4calo->at(sortedJetIdx[ijet])
		<< endl;
	 }
	 //////////////cout << "id Tight jet" << sortedJetIdx[1] << " = " << idTAK4calo->at(sortedJetIdx[1]) << endl;
	 if(fabs(jetEtaAK4calo->at(sortedJetIdx[ijet])) < getPreCutValue1("jetFidRegion")
	    && idCaloJet[sortedJetIdx[ijet]] == getPreCutValue1("tightJetID")
	    && (jecFactors[sortedJetIdx[ijet]]/jetJecAK4calo->at(sortedJetIdx[ijet]))*jetPtAK4calo->at(sortedJetIdx[ijet]) > getPreCutValue1("pt1Cut"))
	   {
	     Nak4 += 1;
	     HTak4 += (jecFactors[sortedJetIdx[ijet]]/jetJecAK4calo->at(sortedJetIdx[ijet]))*jetPtAK4calo->at(sortedJetIdx[ijet]);
	   }
       }
   
     // +++++++++++++++++ Reco analysis
     
     if( int(getPreCutValue1("useFastJet"))==1 )
       {
	 // vector of ak4 jets used for wide jet clustering
	 std::vector<fastjet::PseudoJet> fjInputs, fjInputs_shift;
	 
	 for(size_t j=0; j<no_jetsReco_ak4; ++j)
	   {
	     if( !(jetEtaAK4reco->at(sortedRecoJetIdx[j]) < getPreCutValue1("jetFidRegion")
		   && idTAK4reco->at(sortedRecoJetIdx[j]) == getPreCutValue1("tightJetID")) ) continue;
	     
	     double rescale = (jecFactorsReco[sortedRecoJetIdx[j]]/jetJecAK4reco->at(sortedRecoJetIdx[j]));
	     
	     if( j==0 && !( rescale*jetPtAK4reco->at(sortedRecoJetIdx[j]) > getPreCutValue1("pt0Cut")) ) continue;
	     else if( j==1 && !( rescale*jetPtAK4reco->at(sortedRecoJetIdx[j]) > getPreCutValue1("pt1Cut")) ) continue;
	     else if( !( rescale*jetPtAK4reco->at(sortedRecoJetIdx[j]) > getPreCutValue1("pt1Cut")) ) continue;
	     
	     TLorentzVector tempJet, tempJet_shift;
	     
	     tempJet.SetPtEtaPhiM( rescale*jetPtAK4reco->at(sortedRecoJetIdx[j]) , jetEtaAK4reco->at(sortedRecoJetIdx[j]) , jetPhiAK4reco->at(sortedRecoJetIdx[j]) , rescale*jetMassAK4reco->at(sortedRecoJetIdx[j]));
	     tempJet_shift.SetPtEtaPhiM( (1+ jecUncertaintyReco[sortedRecoJetIdx[j]])* rescale*jetPtAK4reco->at(sortedRecoJetIdx[j]) , jetEtaAK4reco->at(sortedRecoJetIdx[j]) , jetPhiAK4reco->at(sortedRecoJetIdx[j]) ,  (1+ jecUncertaintyReco[sortedRecoJetIdx[j]])* rescale*jetMassAK4reco->at(sortedRecoJetIdx[j]));
	     
	     fjInputs.push_back(fastjet::PseudoJet(tempJet.Px(),tempJet.Py(),tempJet.Pz(),tempJet.E()));
	     fjInputs_shift.push_back(fastjet::PseudoJet(tempJet_shift.Px(),tempJet_shift.Py(),tempJet_shift.Pz(),tempJet_shift.E()));
	   }
	 
	 fjClusterSeq = ClusterSequencePtr( new fastjet::ClusterSequence( fjInputs, *fjJetDefinition ) );
	 fjClusterSeq_shift = ClusterSequencePtr( new fastjet::ClusterSequence( fjInputs_shift, *fjJetDefinition ) );
	 
	 std::vector<fastjet::PseudoJet> inclusiveWideJets = fastjet::sorted_by_pt( fjClusterSeq->inclusive_jets(0.) );
	 std::vector<fastjet::PseudoJet> inclusiveWideJets_shift = fastjet::sorted_by_pt( fjClusterSeq_shift->inclusive_jets(0.) );
	 
	 if( inclusiveWideJets.size()>1 )
	   {
	     wj1Reco.SetPxPyPzE(inclusiveWideJets.at(0).px(), inclusiveWideJets.at(0).py(), inclusiveWideJets.at(0).pz(), inclusiveWideJets.at(0).e());
	     wj2Reco.SetPxPyPzE(inclusiveWideJets.at(1).px(), inclusiveWideJets.at(1).py(), inclusiveWideJets.at(1).pz(), inclusiveWideJets.at(1).e());
	     wj1Reco_shift.SetPxPyPzE(inclusiveWideJets_shift.at(0).px(), inclusiveWideJets_shift.at(0).py(), inclusiveWideJets_shift.at(0).pz(), inclusiveWideJets_shift.at(0).e());
	     wj2Reco_shift.SetPxPyPzE(inclusiveWideJets_shift.at(1).px(), inclusiveWideJets_shift.at(1).py(), inclusiveWideJets_shift.at(1).pz(), inclusiveWideJets_shift.at(1).e());
	   }
       }
     else // don't use FastJet
       {
	 if(verbose) std::cout << "Reco Widejet Construction " << std::endl;

	 TLorentzVector wj1_tmp, wj2_tmp;
	 TLorentzVector wj1_shift_tmp, wj2_shift_tmp;
	 double wideJetDeltaR_ = getPreCutValue1("DeltaR");

	 // building widejet
	 if(no_jetsReco_ak4>=2)
	   {
	     if(fabs(jetEtaAK4reco->at(sortedRecoJetIdx[0])) < getPreCutValue1("jetFidRegion") 
		&& (jecFactorsReco[sortedRecoJetIdx[0]]/jetJecAK4reco->at(sortedRecoJetIdx[0]))*jetPtAK4reco->at(sortedRecoJetIdx[0]) > getPreCutValue1("pt0Cut"))
	       {
		 if(fabs(jetEtaAK4reco->at(sortedRecoJetIdx[1])) < getPreCutValue1("jetFidRegion") 
		    && (jecFactorsReco[sortedRecoJetIdx[1]]/jetJecAK4reco->at(sortedRecoJetIdx[1]))*jetPtAK4reco->at(sortedRecoJetIdx[1]) > getPreCutValue1("pt1Cut"))
		   {

		     TLorentzVector jet1, jet2, jet1_shift, jet2_shift;
		     jet1.SetPtEtaPhiM( (jecFactorsReco[sortedRecoJetIdx[0]]/jetJecAK4reco->at(sortedRecoJetIdx[0])) *jetPtAK4reco->at(sortedRecoJetIdx[0])
					,jetEtaAK4reco->at(sortedRecoJetIdx[0]),jetPhiAK4reco->at(sortedRecoJetIdx[0])
					, (jecFactorsReco[sortedRecoJetIdx[0]]/jetJecAK4reco->at(sortedRecoJetIdx[0])) * jetMassAK4reco->at(sortedRecoJetIdx[0]));
		     jet2.SetPtEtaPhiM( (jecFactorsReco[sortedRecoJetIdx[1]]/jetJecAK4reco->at(sortedRecoJetIdx[1])) *jetPtAK4reco->at(sortedRecoJetIdx[1])
					,jetEtaAK4reco->at(sortedRecoJetIdx[1]),jetPhiAK4reco->at(sortedRecoJetIdx[1])
					, (jecFactorsReco[sortedRecoJetIdx[1]]/jetJecAK4reco->at(sortedRecoJetIdx[1])) * jetMassAK4reco->at(sortedRecoJetIdx[1]));
		     jet1_shift.SetPtEtaPhiM( (1+ jecUncertaintyReco[sortedRecoJetIdx[0]])*(jecFactorsReco[sortedRecoJetIdx[0]]/jetJecAK4reco->at(sortedRecoJetIdx[0])) *jetPtAK4reco->at(sortedRecoJetIdx[0])
					      ,jetEtaAK4reco->at(sortedRecoJetIdx[0]),jetPhiAK4reco->at(sortedRecoJetIdx[0])
					      , (1+ jecUncertaintyReco[sortedRecoJetIdx[0]])*(jecFactorsReco[sortedRecoJetIdx[0]]/jetJecAK4reco->at(sortedRecoJetIdx[0])) * jetMassAK4reco->at(sortedRecoJetIdx[0]));
		     jet2_shift.SetPtEtaPhiM( (1+ jecUncertaintyReco[sortedRecoJetIdx[1]])* (jecFactorsReco[sortedRecoJetIdx[1]]/jetJecAK4reco->at(sortedRecoJetIdx[1])) *jetPtAK4reco->at(sortedRecoJetIdx[1])
					      ,jetEtaAK4reco->at(sortedRecoJetIdx[1]),jetPhiAK4reco->at(sortedRecoJetIdx[1])
					      , (1+ jecUncertaintyReco[sortedRecoJetIdx[1]])*(jecFactorsReco[sortedRecoJetIdx[1]]/jetJecAK4reco->at(sortedRecoJetIdx[1])) * jetMassAK4reco->at(sortedRecoJetIdx[1]));

		     
		     for(Long64_t ijet=0; ijet<no_jetsReco_ak4; ijet++)
		       { //jet loop for ak4
			 TLorentzVector currentJet;
			 
			 if(fabs(jetEtaAK4reco->at(sortedRecoJetIdx[ijet])) < getPreCutValue1("jetFidRegion") 
			    && idTAK4reco->at(sortedRecoJetIdx[ijet]) == getPreCutValue1("tightJetID") 
			    && (jecFactorsReco[sortedRecoJetIdx[ijet]]/jetJecAK4reco->at(sortedRecoJetIdx[ijet]))*jetPtAK4reco->at(sortedRecoJetIdx[ijet]) > getPreCutValue1("pt1Cut"))
			   {
			     TLorentzVector currentJet, currentJet_shift;
			     currentJet.SetPtEtaPhiM( (jecFactorsReco[sortedRecoJetIdx[ijet]]/jetJecAK4reco->at(sortedRecoJetIdx[ijet])) *jetPtAK4reco->at(sortedRecoJetIdx[ijet])
						      ,jetEtaAK4reco->at(sortedRecoJetIdx[ijet]),jetPhiAK4reco->at(sortedRecoJetIdx[ijet])
						      , (jecFactorsReco[sortedRecoJetIdx[ijet]]/jetJecAK4reco->at(sortedRecoJetIdx[ijet])) *jetMassAK4reco->at(sortedRecoJetIdx[ijet]));   
			     currentJet_shift.SetPtEtaPhiM( (1+ jecUncertaintyReco[sortedRecoJetIdx[ijet]])*(jecFactorsReco[sortedRecoJetIdx[ijet]]/jetJecAK4reco->at(sortedRecoJetIdx[ijet])) *jetPtAK4reco->at(sortedRecoJetIdx[ijet])
							    ,jetEtaAK4reco->at(sortedRecoJetIdx[ijet]),jetPhiAK4reco->at(sortedRecoJetIdx[ijet])
							    , (1+ jecUncertaintyReco[sortedRecoJetIdx[ijet]])*(jecFactorsReco[sortedRecoJetIdx[ijet]]/jetJecAK4reco->at(sortedRecoJetIdx[ijet])) *jetMassAK4reco->at(sortedRecoJetIdx[ijet]));   
			   
			     double DeltaR1 = currentJet.DeltaR(jet1);
			     double DeltaR2 = currentJet.DeltaR(jet2);
			     
			     if(DeltaR1 < DeltaR2 && DeltaR1 < wideJetDeltaR_)
			       {
				 wj1_tmp += currentJet;
				 wj1_shift_tmp += currentJet_shift;
			       }
			     else if(DeltaR2 < wideJetDeltaR_)
			       {
				 wj2_tmp += currentJet;
				 wj2_shift_tmp += currentJet_shift;
			       }			 
			   } // if AK4reco jet passes fid and jetid.
		       } //end of ak4 jet loop		     
		     
		     // if(wj1_tmp.Pt()==0 && wj2_tmp.Pt() ==0) 
		     // std::cout << " wj1_tmp.Pt() IN  " <<wj1_tmp.Pt()  << " wj2_tmp.Pt() " <<  wj2_tmp.Pt()  << std::endl;		     
		     
		   } //fid, jet id, pt cut
	       } //fid, jet id, pt cut
	   } // end of two jets.
	 
	 // Re-order the wide jets in pt
	 if( wj1_tmp.Pt() > wj2_tmp.Pt())
	   {
	     wj1Reco = wj1_tmp;
	     wj2Reco = wj2_tmp;
	     wj1Reco_shift = wj1_shift_tmp;
	     wj2Reco_shift = wj2_shift_tmp;
	   }
	 else
	   {
	     wj1Reco = wj2_tmp;
	     wj2Reco = wj1_tmp;
	     wj1Reco_shift = wj2_shift_tmp;
	     wj2Reco_shift = wj1_shift_tmp;
	   }
       }// end building wide jet reco

     if(verbose) std::cout << "widejet1 Reco Pt = "<<wj1Reco.Pt() << std::endl;     
     if(verbose) std::cout << "widejet2 Reco Pt = "<<wj2Reco.Pt() << std::endl;     
     
     double MJJWideReco = 0; 
     double DeltaEtaJJWideReco = 0;
     double DeltaPhiJJWideReco = 0;
     double MJJWideReco_shift = 0; 
     if( wj1Reco.Pt()>0 && wj2Reco.Pt()>0 )
       {
	 // Create dijet system
	 wdijetReco = wj1Reco + wj2Reco;
	 MJJWideReco = wdijetReco.M();
	 DeltaEtaJJWideReco = fabs(wj1Reco.Eta()-wj2Reco.Eta());
	 DeltaPhiJJWideReco = fabs(wj1Reco.DeltaPhi(wj2Reco));
	 
	 wdijetReco_shift = wj1Reco_shift + wj2Reco_shift;
	 MJJWideReco_shift = wdijetReco_shift.M();
	 
	 // Put widejets in the container
	 widejetsReco.push_back( wj1Reco );
	 widejetsReco.push_back( wj2Reco );
       }
       
     //AK4reco jets
     if(no_jetsReco_ak4>=2)
       {
	 if(fabs(jetEtaAK4reco->at(sortedRecoJetIdx[0])) < getPreCutValue1("jetFidRegion") 
	    && (jecFactorsReco[sortedRecoJetIdx[0]]/jetJecAK4reco->at(sortedRecoJetIdx[0]))*jetPtAK4reco->at(sortedRecoJetIdx[0]) > getPreCutValue1("pt0Cut"))
	   {
	     if(fabs(jetEtaAK4reco->at(sortedRecoJetIdx[1])) < getPreCutValue1("jetFidRegion") 
		&& (jecFactorsReco[sortedRecoJetIdx[1]]/jetJecAK4reco->at(sortedRecoJetIdx[1]))*jetPtAK4reco->at(sortedRecoJetIdx[1]) > getPreCutValue1("pt1Cut"))
	       {

		 ak4j1Reco.SetPtEtaPhiM( (jecFactorsReco[sortedRecoJetIdx[0]]/jetJecAK4reco->at(sortedRecoJetIdx[0])) *jetPtAK4reco->at(sortedRecoJetIdx[0])
					 ,jetEtaAK4reco->at(sortedRecoJetIdx[0])
					 ,jetPhiAK4reco->at(sortedRecoJetIdx[0])
					 , (jecFactorsReco[sortedRecoJetIdx[0]]/jetJecAK4reco->at(sortedRecoJetIdx[0])) *jetMassAK4reco->at(sortedRecoJetIdx[0]));
		 ak4j2Reco.SetPtEtaPhiM( (jecFactorsReco[sortedRecoJetIdx[1]]/jetJecAK4reco->at(sortedRecoJetIdx[1])) *jetPtAK4reco->at(sortedRecoJetIdx[1])
					 ,jetEtaAK4reco->at(sortedRecoJetIdx[1])
					 ,jetPhiAK4reco->at(sortedRecoJetIdx[1])
					 , (jecFactorsReco[sortedRecoJetIdx[1]]/jetJecAK4reco->at(sortedRecoJetIdx[1])) *jetMassAK4reco->at(sortedRecoJetIdx[1]));

		 // saving third jet
		 if(no_jetsReco_ak4>=3){
		   if(fabs(jetEtaAK4reco->at(sortedRecoJetIdx[2])) < getPreCutValue1("jetFidRegion")
		      && idTAK4reco->at(sortedRecoJetIdx[2]) == getPreCutValue1("tightJetID") 
		      && (jecFactorsReco[sortedRecoJetIdx[2]]/jetJecAK4reco->at(sortedRecoJetIdx[2]))*jetPtAK4reco->at(sortedRecoJetIdx[2]) > getPreCutValue1("pt1Cut") ){
		     
		     thirdRecoJetPresent = true;

		     ak4j3Reco.SetPtEtaPhiM( (jecFactorsReco[sortedRecoJetIdx[2]]/jetJecAK4reco->at(sortedRecoJetIdx[2])) *jetPtAK4reco->at(sortedRecoJetIdx[2])
					     ,jetEtaAK4reco->at(sortedRecoJetIdx[2])
					     ,jetPhiAK4reco->at(sortedRecoJetIdx[2])
					     , (jecFactorsReco[sortedRecoJetIdx[2]]/jetJecAK4reco->at(sortedRecoJetIdx[2])) *jetMassAK4reco->at(sortedRecoJetIdx[2]));
		   
		   }
		 }		 
	       }
	   }
       }// end of two reco jets
     
     double MJJAK4reco = 0; 
     double DeltaEtaJJAK4reco = 0;
     double DeltaPhiJJAK4reco = 0;
     
     if( ak4j1Reco.Pt()>0 && ak4j2Reco.Pt()>0 )
       {
	 // Create dijet system
	 ak4dijetReco = ak4j1Reco + ak4j2Reco;
	 MJJAK4reco = ak4dijetReco.M();
	 DeltaEtaJJAK4reco = fabs(ak4j1Reco.Eta()-ak4j2Reco.Eta());
	 DeltaPhiJJAK4reco = fabs(ak4j1Reco.DeltaPhi(ak4j2Reco));
	 
	 // Put widejets in the container
	 AK4recojets.push_back( ak4j1Reco );
	 AK4recojets.push_back( ak4j2Reco );
	 if( ak4j3Reco.Pt()>0 )
	   {
	     AK4recojets.push_back( ak4j3Reco );
	   }
       }

     // +++++++++++++++++++++ HLT Analysis

     if(verbose)     std::cout<<"Start HLT analysis"<<std::endl;  

     if( int(getPreCutValue1("useFastJet"))==1 )
       {
	 // vector of ak4 jets used for wide jet clustering
	 std::vector<fastjet::PseudoJet> fjInputs, fjInputs_shift;
	 
	 for(size_t j=0; j<no_jets_ak4; ++j)
	   {
	     if( !(jetEtaAK4calo->at(sortedJetIdx[j]) < getPreCutValue1("jetFidRegion")
		   && idCaloJet[sortedJetIdx[j]] == getPreCutValue1("tightJetID")) ) continue;
	     
	     double rescale = (jecFactors[sortedJetIdx[j]]/jetJecAK4calo->at(sortedJetIdx[j]));
	     
	     if( j==0 && !( rescale*jetPtAK4calo->at(sortedJetIdx[j]) > getPreCutValue1("pt0Cut")) ) continue;
	     else if( j==1 && !( rescale*jetPtAK4calo->at(sortedJetIdx[j]) > getPreCutValue1("pt1Cut")) ) continue;
	     else if( !( rescale*jetPtAK4calo->at(sortedJetIdx[j]) > getPreCutValue1("pt1Cut")) ) continue;
	     
	     TLorentzVector tempJet, tempJet_shift;
	     
	     tempJet.SetPtEtaPhiM( rescale*jetPtAK4calo->at(sortedJetIdx[j]) , jetEtaAK4calo->at(sortedJetIdx[j]) , jetPhiAK4calo->at(sortedJetIdx[j]) , rescale*jetMassAK4calo->at(sortedJetIdx[j]));
	     tempJet_shift.SetPtEtaPhiM( (1+jecUncertainty[sortedJetIdx[j]])* rescale*jetPtAK4calo->at(sortedJetIdx[j]) , jetEtaAK4calo->at(sortedJetIdx[j]) , jetPhiAK4calo->at(sortedJetIdx[j]) ,  (1+jecUncertainty[sortedJetIdx[j]])* rescale*jetMassAK4calo->at(sortedJetIdx[j]));
	     
	     fjInputs.push_back(fastjet::PseudoJet(tempJet.Px(),tempJet.Py(),tempJet.Pz(),tempJet.E()));
	     fjInputs_shift.push_back(fastjet::PseudoJet(tempJet_shift.Px(),tempJet_shift.Py(),tempJet_shift.Pz(),tempJet_shift.E()));
	   }
	 
	 fjClusterSeq = ClusterSequencePtr( new fastjet::ClusterSequence( fjInputs, *fjJetDefinition ) );
	 fjClusterSeq_shift = ClusterSequencePtr( new fastjet::ClusterSequence( fjInputs_shift, *fjJetDefinition ) );
	 
       std::vector<fastjet::PseudoJet> inclusiveWideJets = fastjet::sorted_by_pt( fjClusterSeq->inclusive_jets(0.) );
       std::vector<fastjet::PseudoJet> inclusiveWideJets_shift = fastjet::sorted_by_pt( fjClusterSeq_shift->inclusive_jets(0.) );
       
       if( inclusiveWideJets.size()>1 )
	 {
	   wj1.SetPxPyPzE(inclusiveWideJets.at(0).px(), inclusiveWideJets.at(0).py(), inclusiveWideJets.at(0).pz(), inclusiveWideJets.at(0).e());
	   wj2.SetPxPyPzE(inclusiveWideJets.at(1).px(), inclusiveWideJets.at(1).py(), inclusiveWideJets.at(1).pz(), inclusiveWideJets.at(1).e());
	   wj1_shift.SetPxPyPzE(inclusiveWideJets_shift.at(0).px(), inclusiveWideJets_shift.at(0).py(), inclusiveWideJets_shift.at(0).pz(), inclusiveWideJets_shift.at(0).e());
	   wj2_shift.SetPxPyPzE(inclusiveWideJets_shift.at(1).px(), inclusiveWideJets_shift.at(1).py(), inclusiveWideJets_shift.at(1).pz(), inclusiveWideJets_shift.at(1).e());
	 }
       }
     else // don't use FastJet
       {
	 if(verbose) std::cout << "Widejet Construction " << std::endl;

	 TLorentzVector wj1_tmp, wj2_tmp;
	 TLorentzVector wj1_shift_tmp, wj2_shift_tmp;
	 double wideJetDeltaR_ = getPreCutValue1("DeltaR");
	 // building widejets
	 if(no_jets_ak4>=2)
	   {
	     if(fabs(jetEtaAK4calo->at(sortedJetIdx[0])) < getPreCutValue1("jetFidRegion") 
		&& (jecFactors[sortedJetIdx[0]]/jetJecAK4calo->at(sortedJetIdx[0]))*jetPtAK4calo->at(sortedJetIdx[0]) > getPreCutValue1("pt0Cut"))
	       {
		 if(fabs(jetEtaAK4calo->at(sortedJetIdx[1])) < getPreCutValue1("jetFidRegion") 
		    && (jecFactors[sortedJetIdx[1]]/jetJecAK4calo->at(sortedJetIdx[1]))*jetPtAK4calo->at(sortedJetIdx[1]) > getPreCutValue1("pt1Cut"))
		   {
		     TLorentzVector jet1, jet2, jet1_shift, jet2_shift;
		     jet1.SetPtEtaPhiM( (jecFactors[sortedJetIdx[0]]/jetJecAK4calo->at(sortedJetIdx[0])) *jetPtAK4calo->at(sortedJetIdx[0])
					,jetEtaAK4calo->at(sortedJetIdx[0]),jetPhiAK4calo->at(sortedJetIdx[0])
					, (jecFactors[sortedJetIdx[0]]/jetJecAK4calo->at(sortedJetIdx[0])) * jetMassAK4calo->at(sortedJetIdx[0]));
		     jet2.SetPtEtaPhiM( (jecFactors[sortedJetIdx[1]]/jetJecAK4calo->at(sortedJetIdx[1])) *jetPtAK4calo->at(sortedJetIdx[1])
					,jetEtaAK4calo->at(sortedJetIdx[1]),jetPhiAK4calo->at(sortedJetIdx[1])
					, (jecFactors[sortedJetIdx[1]]/jetJecAK4calo->at(sortedJetIdx[1])) * jetMassAK4calo->at(sortedJetIdx[1]));
		     jet1_shift.SetPtEtaPhiM( (1+jecUncertainty[sortedJetIdx[0]])*(jecFactors[sortedJetIdx[0]]/jetJecAK4calo->at(sortedJetIdx[0])) *jetPtAK4calo->at(sortedJetIdx[0])
					      ,jetEtaAK4calo->at(sortedJetIdx[0]),jetPhiAK4calo->at(sortedJetIdx[0])
					      , (1+jecUncertainty[sortedJetIdx[0]])*(jecFactors[sortedJetIdx[0]]/jetJecAK4calo->at(sortedJetIdx[0])) * jetMassAK4calo->at(sortedJetIdx[0]));
		     jet2_shift.SetPtEtaPhiM( (1+jecUncertainty[sortedJetIdx[1]])* (jecFactors[sortedJetIdx[1]]/jetJecAK4calo->at(sortedJetIdx[1])) *jetPtAK4calo->at(sortedJetIdx[1])
					      ,jetEtaAK4calo->at(sortedJetIdx[1]),jetPhiAK4calo->at(sortedJetIdx[1])
					      , (1+jecUncertainty[sortedJetIdx[1]])*(jecFactors[sortedJetIdx[1]]/jetJecAK4calo->at(sortedJetIdx[1])) * jetMassAK4calo->at(sortedJetIdx[1]));
		     
		     for(Long64_t ijet=0; ijet<no_jets_ak4; ijet++)
		       { //jet loop for ak4
			 TLorentzVector currentJet;
			 
			 if(fabs(jetEtaAK4calo->at(sortedJetIdx[ijet])) < getPreCutValue1("jetFidRegion") 
			    && idCaloJet[sortedJetIdx[ijet]] == getPreCutValue1("tightJetID") 
			    && (jecFactors[sortedJetIdx[ijet]]/jetJecAK4calo->at(sortedJetIdx[ijet]))*jetPtAK4calo->at(sortedJetIdx[ijet]) > getPreCutValue1("pt1Cut"))
			   {
			     TLorentzVector currentJet, currentJet_shift;
			     currentJet.SetPtEtaPhiM( (jecFactors[sortedJetIdx[ijet]]/jetJecAK4calo->at(sortedJetIdx[ijet])) *jetPtAK4calo->at(sortedJetIdx[ijet])
						      ,jetEtaAK4calo->at(sortedJetIdx[ijet]),jetPhiAK4calo->at(sortedJetIdx[ijet])
						      , (jecFactors[sortedJetIdx[ijet]]/jetJecAK4calo->at(sortedJetIdx[ijet])) *jetMassAK4calo->at(sortedJetIdx[ijet]));   
			     currentJet_shift.SetPtEtaPhiM( (1+jecUncertainty[sortedJetIdx[ijet]])*(jecFactors[sortedJetIdx[ijet]]/jetJecAK4calo->at(sortedJetIdx[ijet])) *jetPtAK4calo->at(sortedJetIdx[ijet])
							    ,jetEtaAK4calo->at(sortedJetIdx[ijet]),jetPhiAK4calo->at(sortedJetIdx[ijet])
							    , (1+jecUncertainty[sortedJetIdx[ijet]])*(jecFactors[sortedJetIdx[ijet]]/jetJecAK4calo->at(sortedJetIdx[ijet])) *jetMassAK4calo->at(sortedJetIdx[ijet]));   
			     
			     double DeltaR1 = currentJet.DeltaR(jet1);
			     double DeltaR2 = currentJet.DeltaR(jet2);
			     
			     if(DeltaR1 < DeltaR2 && DeltaR1 < wideJetDeltaR_)
			       {
				 wj1_tmp += currentJet;
				 wj1_shift_tmp += currentJet_shift;
			       }
			     else if(DeltaR2 < wideJetDeltaR_)
			       {
				 wj2_tmp += currentJet;
				 wj2_shift_tmp += currentJet_shift;
			       }			 
			   } // if AK4 jet passes fid and jetid.
		       } //end of ak4 jet loop		     
		     
		     // if(wj1_tmp.Pt()==0 && wj2_tmp.Pt() ==0) 
		     // std::cout << " wj1_tmp.Pt() IN  " <<wj1_tmp.Pt()  << " wj2_tmp.Pt() " <<  wj2_tmp.Pt()  << std::endl;		     
		     
		   } //fid, jet id, pt cut
	       } //fid, jet id, pt cut
	   } // end of two jets.
	 
	 // Re-order the wide jets in pt
	 if( wj1_tmp.Pt() > wj2_tmp.Pt())
	   {
	     wj1 = wj1_tmp;
	     wj2 = wj2_tmp;
	     wj1_shift = wj1_shift_tmp;
	     wj2_shift = wj2_shift_tmp;
	   }
	 else
	   {
	     wj1 = wj2_tmp;
	     wj2 = wj1_tmp;
	     wj1_shift = wj2_shift_tmp;
	     wj2_shift = wj1_shift_tmp;
	   }
       }// end building wide jet 
     
     if(verbose) std::cout << "widejet1 Pt = "<<wj1.Pt() << std::endl;     
     if(verbose) std::cout << "widejet2 Pt = "<<wj2.Pt() << std::endl;     
   
     double MJJWide = 0; 
     double DeltaEtaJJWide = 0;
     double DeltaPhiJJWide = 0;
     double MJJWide_shift = 0; 
     if( wj1.Pt()>0 && wj2.Pt()>0 )
       {
	 // Create dijet system
	 wdijet = wj1 + wj2;
	 MJJWide = wdijet.M();
	 DeltaEtaJJWide = fabs(wj1.Eta()-wj2.Eta());
	 DeltaPhiJJWide = fabs(wj1.DeltaPhi(wj2));
	 
	 wdijet_shift = wj1_shift + wj2_shift;
	 MJJWide_shift = wdijet_shift.M();
	 
	 // Put widejets in the container
	 widejets.push_back( wj1 );
	 widejets.push_back( wj2 );
       }
     
     //AK4 jets
     if(no_jets_ak4>=2)
       {
	 if(fabs(jetEtaAK4calo->at(sortedJetIdx[0])) < getPreCutValue1("jetFidRegion") 
	    && (jecFactors[sortedJetIdx[0]]/jetJecAK4calo->at(sortedJetIdx[0]))*jetPtAK4calo->at(sortedJetIdx[0]) > getPreCutValue1("pt0Cut"))
	   {
	     if(fabs(jetEtaAK4calo->at(sortedJetIdx[1])) < getPreCutValue1("jetFidRegion") 
		&& (jecFactors[sortedJetIdx[1]]/jetJecAK4calo->at(sortedJetIdx[1]))*jetPtAK4calo->at(sortedJetIdx[1]) > getPreCutValue1("pt1Cut"))
	       {
		 ak4j1.SetPtEtaPhiM( (jecFactors[sortedJetIdx[0]]/jetJecAK4calo->at(sortedJetIdx[0])) *jetPtAK4calo->at(sortedJetIdx[0])
				     ,jetEtaAK4calo->at(sortedJetIdx[0])
				     ,jetPhiAK4calo->at(sortedJetIdx[0])
				     , (jecFactors[sortedJetIdx[0]]/jetJecAK4calo->at(sortedJetIdx[0])) *jetMassAK4calo->at(sortedJetIdx[0]));
		 ak4j2.SetPtEtaPhiM( (jecFactors[sortedJetIdx[1]]/jetJecAK4calo->at(sortedJetIdx[1])) *jetPtAK4calo->at(sortedJetIdx[1])
				     ,jetEtaAK4calo->at(sortedJetIdx[1])
				     ,jetPhiAK4calo->at(sortedJetIdx[1])
				     , (jecFactors[sortedJetIdx[1]]/jetJecAK4calo->at(sortedJetIdx[1])) *jetMassAK4calo->at(sortedJetIdx[1]));
		 
		 // saving third jet
		 if(no_jets_ak4>=3){
		   if(fabs(jetEtaAK4calo->at(sortedJetIdx[2])) < getPreCutValue1("jetFidRegion")
		      && idCaloJet[sortedJetIdx[2]] == getPreCutValue1("tightJetID")
		      && (jecFactors[sortedJetIdx[2]]/jetJecAK4calo->at(sortedJetIdx[2]))*jetPtAK4calo->at(sortedJetIdx[2]) > getPreCutValue1("pt1Cut") ){

		     thirdJetPresent = true;
		     
		     ak4j3.SetPtEtaPhiM( (jecFactors[sortedJetIdx[2]]/jetJecAK4calo->at(sortedJetIdx[2])) *jetPtAK4calo->at(sortedJetIdx[2])
					 ,jetEtaAK4calo->at(sortedJetIdx[2])
					 ,jetPhiAK4calo->at(sortedJetIdx[2])
					 , (jecFactors[sortedJetIdx[2]]/jetJecAK4calo->at(sortedJetIdx[2])) *jetMassAK4calo->at(sortedJetIdx[2]));
		   
		   }
		 }
	       }
	   }
       }
     
     double MJJAK4 = 0; 
     double DeltaEtaJJAK4 = 0;
     double DeltaPhiJJAK4 = 0;
     
     if( ak4j1.Pt()>0 && ak4j2.Pt()>0 )
       {
	 // Create dijet system
	 ak4dijet = ak4j1 + ak4j2;
	 MJJAK4 = ak4dijet.M();
	 DeltaEtaJJAK4 = fabs(ak4j1.Eta()-ak4j2.Eta());
	 DeltaPhiJJAK4 = fabs(ak4j1.DeltaPhi(ak4j2));
	 
	 // Put widejets in the container
	 AK4jets.push_back( ak4j1 );
	 AK4jets.push_back( ak4j2 );
	 if( ak4j3.Pt()>0 )
	   {
	     AK4jets.push_back( ak4j3 );
	   }
       }

     //////////////////////
     if(verbose)     std::cout<<"Start fillVariableWithValue"<<std::endl;

     //== Fill Variables ==
     fillVariableWithValue("isData",isData);     
     fillVariableWithValue("run",runNo);     
     fillVariableWithValue("event",evtNo);     
     fillVariableWithValue("lumi",lumi);       
     fillVariableWithValue("nVtx",nvtx);     
     fillVariableWithValue("nVtxreco",nvtxreco);     
     fillVariableWithValue("nJet",widejets.size());
     fillVariableWithValue("nJetreco",widejetsReco.size());
     fillVariableWithValue("Nak4reco",Nak4Reco);
     fillVariableWithValue("Nak4",Nak4);
     fillVariableWithValue( "PassJSON", passJSON (runNo, lumi, isData));

     //directly taken from big root tree (i.e. jec not reapplied)
     fillVariableWithValue("htAK4",htAK4calo); // summing all jets with minimum pT cut and no jetid cut (jec not reapplied)
     fillVariableWithValue("htAK4reco",htAK4reco); // summing all jets with minimum pT cut and no jetid cut (jec not reapplied)
     fillVariableWithValue("mhtAK4",mhtAK4calo); //summing all jets with minimum pT cut and no jetid cut (jec not reapplied)
     fillVariableWithValue("mhtAK4reco",mhtAK4reco); //summing all jets with minimum pT cut and no jetid cut (jec not reapplied)
     fillVariableWithValue("mhtAK4Sig",mhtAK4caloSig); // mhtAK4/htAK4 summing all jets with minimum pT cut and no jetid cut (jec not reapplied)
     fillVariableWithValue("mhtAK4recoSig",mhtAK4recoSig); // mhtAK4/htAK4 summing all jets with minimum pT cut and no jetid cut (jec not reapplied)
     // fillVariableWithValue("offMet",offMet); //mpted from PF candidates (off=offline)
     // fillVariableWithValue("offMetSig",offMetSig); // offMet/offSumEt mputed from PF candidates (off=offline)
     fillVariableWithValue("met",metcalo); //directly taken from event
     fillVariableWithValue("metreco",metreco); //directly taken from event
     fillVariableWithValue("metSig",metcaloSig); // met/offMetSig (to be substituted with met/SumEt from event in future)
     fillVariableWithValue("metrecoSig",metrecoSig); // met/offMetSig (to be substituted with met/SumEt from event in future)
 
     /////////////////////// RECO Variables

     fillVariableWithValue( "NAK4_recojets",AK4recojets.size());
     fillVariableWithValue( "NAK4_jets",AK4jets.size());

     if( AK4recojets.size() >=1 ){
       fillVariableWithValue( "IdTight_recoj1",idTAK4reco->at(sortedRecoJetIdx[0]));
       fillVariableWithValue( "pTAK4_recoj1", AK4recojets[0].Pt());
       fillVariableWithValue( "etaAK4_recoj1", AK4recojets[0].Eta());
       fillVariableWithValue( "phiAK4_recoj1", AK4recojets[0].Phi());
       fillVariableWithValue( "massAK4_recoj1", AK4recojets[0].M());
       //fillVariableWithValue( "jetPtAK4matchCaloJet_j1", jetPtAK4matchCaloJet->at(sortedRecoJetIdx[0]));
       fillVariableWithValue( "jetJecAK4_recoj1", jecFactorsReco[sortedRecoJetIdx[0]] );
       fillVariableWithValue( "jetJecUncAK4_recoj1", jecUncertaintyReco[sortedRecoJetIdx[0]] );
       fillVariableWithValue( "jetCSVAK4_recoj1", jetCSVAK4reco->at(sortedRecoJetIdx[0]) );
       fillVariableWithValue( "chargedHadEnFrac_recoj1", jetChfAK4reco->at(sortedRecoJetIdx[0]));
       fillVariableWithValue( "neutrHadEnFrac_recoj1", jetNhfAK4reco->at(sortedRecoJetIdx[0]));
       double hadf_recoj1 = jetChfAK4reco->at(sortedRecoJetIdx[0]) +jetNhfAK4reco->at(sortedRecoJetIdx[0]);
       fillVariableWithValue( "HadEnFrac_recoj1",  hadf_recoj1);
       fillVariableWithValue( "photonEnFrac_recoj1", jetPhfAK4reco->at(sortedRecoJetIdx[0]));
       fillVariableWithValue( "eleEnFract_recoj1", jetElfAK4reco->at(sortedRecoJetIdx[0]));
       double emf_recoj1 = jetPhfAK4reco->at(sortedRecoJetIdx[0]) +jetElfAK4reco->at(sortedRecoJetIdx[0]);
       fillVariableWithValue( "ElectromFrac_recoj1",  emf_recoj1);
       fillVariableWithValue( "muEnFract_recoj1", jetMufAK4reco->at(sortedRecoJetIdx[0]));
       fillVariableWithValue( "neutrElectromFrac_recoj1", jetNemfAK4reco->at(sortedRecoJetIdx[0]));
       fillVariableWithValue( "chargedElectromFrac_recoj1", jetCemfAK4reco->at(sortedRecoJetIdx[0]));
       fillVariableWithValue( "chargedMult_recoj1", chMultAK4reco->at(sortedRecoJetIdx[0]));
       fillVariableWithValue( "neutrMult_recoj1", neMultAK4reco->at(sortedRecoJetIdx[0]));
       fillVariableWithValue( "photonMult_recoj1", phoMultAK4reco->at(sortedRecoJetIdx[0]));
     }

     if( AK4recojets.size() >=2 ){
       //cout << "IdTight_recoj2 : " << idTAK4->at(sortedRecoJetIdx[1]) << endl << endl;
       fillVariableWithValue( "IdTight_recoj2",idTAK4reco->at(sortedRecoJetIdx[1]));
       fillVariableWithValue( "pTAK4_recoj2", AK4recojets[1].Pt() );
       fillVariableWithValue( "etaAK4_recoj2", AK4recojets[1].Eta());
       fillVariableWithValue( "phiAK4_recoj2", AK4recojets[1].Phi());
       fillVariableWithValue( "massAK4_recoj2", AK4recojets[1].M());
       //fillVariableWithValue( "jetPtAK4matchCaloJet_recoj2", jetPtAK4matchCaloJet->at(sortedRecoJetIdx[1]));
       fillVariableWithValue( "jetJecAK4_recoj2", jecFactorsReco[sortedRecoJetIdx[1]]); 
       fillVariableWithValue( "jetJecUncAK4_recoj2", jecUncertaintyReco[sortedRecoJetIdx[1]] );
       fillVariableWithValue( "jetCSVAK4_recoj2", jetCSVAK4reco->at(sortedRecoJetIdx[1]) );  
       fillVariableWithValue( "chargedHadEnFrac_recoj2", jetChfAK4reco->at(sortedRecoJetIdx[1]));
       fillVariableWithValue( "neutrHadEnFrac_recoj2", jetNhfAK4reco->at(sortedRecoJetIdx[1]));
       double hadf_recoj2 = jetChfAK4reco->at(sortedRecoJetIdx[1]) +jetNhfAK4reco->at(sortedRecoJetIdx[1]);
       fillVariableWithValue( "HadEnFrac_recoj2",  hadf_recoj2);
       fillVariableWithValue( "photonEnFrac_recoj2", jetPhfAK4reco->at(sortedRecoJetIdx[1]));
       fillVariableWithValue( "eleEnFract_recoj2", jetElfAK4reco->at(sortedRecoJetIdx[1]));
       double emf_recoj2 = jetPhfAK4reco->at(sortedRecoJetIdx[1]) +jetElfAK4reco->at(sortedRecoJetIdx[1]);
       fillVariableWithValue( "ElectromFrac_recoj2",  emf_recoj2);
       fillVariableWithValue( "muEnFract_recoj2", jetMufAK4reco->at(sortedRecoJetIdx[1]));
       fillVariableWithValue( "neutrElectromFrac_recoj2", jetNemfAK4reco->at(sortedRecoJetIdx[1]));
       fillVariableWithValue( "chargedElectromFrac_recoj2", jetCemfAK4reco->at(sortedRecoJetIdx[1]));
       fillVariableWithValue( "chargedMult_recoj2", chMultAK4reco->at(sortedRecoJetIdx[1]));
       fillVariableWithValue( "neutrMult_recoj2", neMultAK4reco->at(sortedRecoJetIdx[1]));
       fillVariableWithValue( "photonMult_recoj2", phoMultAK4reco->at(sortedRecoJetIdx[1]));

       fillVariableWithValue( "Dijet_MassAK4reco", MJJAK4reco) ; 
       fillVariableWithValue( "CosThetaStarAK4reco", TMath::TanH( (AK4recojets[0].Eta()-AK4recojets[1].Eta())/2 )); 
       fillVariableWithValue( "deltaETAjjAK4reco", DeltaEtaJJAK4reco ) ;
       fillVariableWithValue( "deltaPHIjjAK4reco", DeltaPhiJJAK4reco) ;
     }

     if( AK4recojets.size() >=3 ){
       fillVariableWithValue( "IsPresent_recoj3", thirdRecoJetPresent);
       fillVariableWithValue( "IdTight_recoj3",idTAK4reco->at(sortedRecoJetIdx[2]));
       fillVariableWithValue( "pTAK4_recoj3", AK4recojets[2].Pt() );
       fillVariableWithValue( "etaAK4_recoj3", AK4recojets[2].Eta());
       fillVariableWithValue( "phiAK4_recoj3", AK4recojets[2].Phi());
       fillVariableWithValue( "massAK4_recoj3", AK4recojets[2].M());
     }

     if( widejetsReco.size() >= 1 ){

         fillVariableWithValue( "pTWJ_recoj1", widejetsReco[0].Pt() );
         fillVariableWithValue( "etaWJ_recoj1", widejetsReco[0].Eta());
	 //no cuts on these variables, just to store in output
         fillVariableWithValue( "massWJ_recoj1", widejetsReco[0].M());
         fillVariableWithValue( "phiWJ_recoj1", widejetsReco[0].Phi());
     }

     if( widejetsReco.size() >= 2 ){
         fillVariableWithValue( "pTWJ_recoj2", widejetsReco[1].Pt() );
         fillVariableWithValue( "etaWJ_recoj2", widejetsReco[1].Eta());
	 fillVariableWithValue( "deltaETAjjreco", DeltaEtaJJWideReco ) ;
         fillVariableWithValue( "mjjreco", MJJWideReco ) ;
         fillVariableWithValue( "mjjreco_shiftJEC", MJJWideReco_shift ) ;
	 //no cuts on these variables, just to store in output
         fillVariableWithValue( "massWJ_recoj2", widejetsReco[1].M());
         fillVariableWithValue( "phiWJ_recoj2", widejetsReco[1].Phi());	
	 //dijet
         fillVariableWithValue( "CosThetaStarWJreco", TMath::TanH( (widejetsReco[0].Eta()-widejetsReco[1].Eta())/2 )); 
	 fillVariableWithValue( "deltaPHIjjreco", DeltaPhiJJWideReco ) ;
	 //fillVariableWithValue( "Dijet_MassAK8", mjjAK8 ) ;  
	 //fillVariableWithValue( "Dijet_MassC", mjjCA8 ) ;
	 // if(wdijet.M()<1){
	 //    std::cout << " INV MASS IS " << wdijet.M() << std::endl;
	 //    std::cout << " Delta Eta IS " << DeltaEtaJJWide << " n is  " << widejets.size() << std::endl;
	 //    std::cout << " INV MASS FROM NTUPLE AK8 " << mjjAK8 << std::endl;
	 //    //std::cout << " INV MASS FROM NTUPLE CA8 " << mjjCA8 << std::endl;
     }

     ////////////// HLT variables

     if( AK4jets.size() >=1 ){
       fillVariableWithValue( "IdTight_j1", idCaloJet[sortedJetIdx[0]]);
       fillVariableWithValue( "pTAK4_j1", AK4jets[0].Pt());
       fillVariableWithValue( "etaAK4_j1", AK4jets[0].Eta());
       fillVariableWithValue( "phiAK4_j1", AK4jets[0].Phi());
       fillVariableWithValue( "massAK4_j1", AK4jets[0].M());
       fillVariableWithValue( "jetJecAK4_j1", jecFactors[sortedJetIdx[0]] );
       fillVariableWithValue( "jetJecUncAK4_j1", jecUncertainty[sortedJetIdx[0]] );
       fillVariableWithValue( "jetCSVAK4_j1", jetCSVAK4calo->at(sortedJetIdx[0]) );
       fillVariableWithValue( "HadEnFrac_j1",    jetHadfAK4calo->at(sortedJetIdx[0]));
       fillVariableWithValue( "ElectromFrac_j1", jetEmfAK4calo->at(sortedJetIdx[0]));
     }

     if( AK4jets.size() >=2 ){
       fillVariableWithValue( "IdTight_j2", idCaloJet[sortedJetIdx[1]]);
       fillVariableWithValue( "pTAK4_j2", AK4jets[1].Pt() );
       fillVariableWithValue( "etaAK4_j2", AK4jets[1].Eta());
       fillVariableWithValue( "phiAK4_j2", AK4jets[1].Phi());
       fillVariableWithValue( "massAK4_j2", AK4jets[1].M());
       fillVariableWithValue( "jetJecAK4_j2", jecFactors[sortedJetIdx[1]]); 
       fillVariableWithValue( "jetJecUncAK4_j2", jecUncertainty[sortedJetIdx[1]] );
       fillVariableWithValue( "jetCSVAK4_j2", jetCSVAK4calo->at(sortedJetIdx[1]) );
       fillVariableWithValue( "HadEnFrac_j2",    jetHadfAK4calo->at(sortedJetIdx[1]));
       fillVariableWithValue( "ElectromFrac_j2", jetEmfAK4calo->at(sortedJetIdx[1]));

       fillVariableWithValue( "Dijet_MassAK4", MJJAK4) ; 
       fillVariableWithValue( "CosThetaStarAK4", TMath::TanH( (AK4jets[0].Eta()-AK4jets[1].Eta())/2 )); 
       fillVariableWithValue( "deltaETAjjAK4", DeltaEtaJJAK4 ) ;
       fillVariableWithValue( "deltaPHIjjAK4", DeltaPhiJJAK4) ;
     }

     if( AK4jets.size() >=3 ){
       fillVariableWithValue( "IsPresent_j3", thirdJetPresent);
       fillVariableWithValue( "IdTight_j3", idCaloJet[sortedJetIdx[2]]);
       fillVariableWithValue( "pTAK4_j3", AK4jets[2].Pt() );
       fillVariableWithValue( "etaAK4_j3", AK4jets[2].Eta());
       fillVariableWithValue( "phiAK4_j3", AK4jets[2].Phi());
       fillVariableWithValue( "massAK4_j3", AK4jets[2].M());
     }

     if( widejets.size() >= 1 ){
         fillVariableWithValue( "pTWJ_j1", widejets[0].Pt() );
         fillVariableWithValue( "etaWJ_j1", widejets[0].Eta());
	 //no cuts on these variables, just to store in output
         fillVariableWithValue( "massWJ_j1", widejets[0].M());
         fillVariableWithValue( "phiWJ_j1", widejets[0].Phi());
     }

     if( widejets.size() >= 2 ){
         fillVariableWithValue( "pTWJ_j2", widejets[1].Pt() );
         fillVariableWithValue( "etaWJ_j2", widejets[1].Eta());
	 fillVariableWithValue( "deltaETAjj", DeltaEtaJJWide ) ;
         fillVariableWithValue( "mjj", MJJWide ) ;
         fillVariableWithValue( "mjj_shiftJEC", MJJWide_shift ) ;
	 //no cuts on these variables, just to store in output
         fillVariableWithValue( "massWJ_j2", widejets[1].M());
         fillVariableWithValue( "phiWJ_j2", widejets[1].Phi());	
	 //dijet
         fillVariableWithValue( "CosThetaStarWJ", TMath::TanH( (widejets[0].Eta()-widejets[1].Eta())/2 )); 
	 fillVariableWithValue( "deltaPHIjj", DeltaPhiJJWide ) ;
	 //fillVariableWithValue( "Dijet_MassAK8", mjjAK8 ) ;  
	 //fillVariableWithValue( "Dijet_MassC", mjjCA8 ) ;
	 // if(wdijet.M()<1){
	 //    std::cout << " INV MASS IS " << wdijet.M() << std::endl;
	 //    std::cout << " Delta Eta IS " << DeltaEtaJJWide << " n is  " << widejets.size() << std::endl;
	 //    std::cout << " INV MASS FROM NTUPLE AK8 " << mjjAK8 << std::endl;
	 //    //std::cout << " INV MASS FROM NTUPLE CA8 " << mjjCA8 << std::endl;
     }


     ////////////////////////////////////////

     //no cuts on these variables, just to store in output
     // if(!isData)
     //   fillVariableWithValue("trueVtx",PileupInteractions->at(idx_InTimeBX));
     // else if(isData)
     //   fillVariableWithValue("trueVtx",999);     

     // Trigger
     //int NtriggerBits = triggerResult->size();
     // to be changed
     if (isData)
       {
	 // ZeroBias
	 fillVariableWithValue("passHLT_ZeroBias_BTagScouting",triggerResult->at(12));
	 fillVariableWithValue("passHLT_ZeroBias_CaloScouting_PFScouting",triggerResult->at(13));
	 //DoubleMu
	 fillVariableWithValue("passHLT_L1DoubleMu_BTagScouting",triggerResult->at(14));
	 fillVariableWithValue("passHLT_L1DoubleMu_CaloScouting_PFScouting",triggerResult->at(15));
	 //CaloJet40
	 fillVariableWithValue("passHLT_CaloJet40_BTagScouting",triggerResult->at(0));
	 fillVariableWithValue("passHLT_CaloJet40_CaloScouting_PFScouting",triggerResult->at(1));
	 fillVariableWithValue("passHLT_CaloJet40_CaloBTagScouting",triggerResult->at(2));
	 //L1HTT
	 fillVariableWithValue("passHLT_L1HTT_BTagScouting",triggerResult->at(3));
	 fillVariableWithValue("passHLT_L1HTT_CaloScouting_PFScouting",triggerResult->at(4));
	 fillVariableWithValue("passHLT_L1HTT_CaloBTagScouting",triggerResult->at(5));
	 // HT250 
	 fillVariableWithValue("passHLT_CaloScoutingHT250",triggerResult->at(6));                         
	 fillVariableWithValue("passHLT_CaloBTagScoutingHT250",triggerResult->at(7));                         
	 // HT410
	 fillVariableWithValue("passHLT_BTagScoutingHT410",triggerResult->at(8));                         
	 fillVariableWithValue("passHLT_PFScoutingHT410",triggerResult->at(9));                         	 
	 // HT650
	 fillVariableWithValue("passHLT_PFHT650",triggerResult->at(20));                         	 

       }

     // Evaluate cuts (but do not apply them)
     evaluateCuts();
     /*
     cout<< "Cuts:"<< endl;
     cout<< "PassJSON: "<< passJSON(runNo, lumi, isData) << endl;
     cout<< "nVtx: "<< nvtx << endl;
     cout<< "nVtxreco: "<< nvtxreco << endl;
     cout<< "IdTIght_j1: "<< idTAK4->at(sortedJetIdx[0]) << endl;
     cout<< "IdTIght_j2: "<< idTAK4->at(sortedJetIdx[1]) << endl;
     cout<< "IdTIght_recoj1: "<< idTAK4reco->at(sortedRecoJetIdx[0]) << endl;
     cout<< "IdTIght_recoj2: "<< idTAK4reco->at(sortedRecoJetIdx[1]) << endl;
     cout<< "nJet: "<< widejets.size() << endl;
     cout<< "nJetReco: "<< widejetsReco.size() << endl;
     if(widejetsReco.size()>=1){
       cout<< "pTWJ_recoj1: "<< AK4recojets[0].Pt()   << endl;
       cout<< "etaWJ_recoj1: "<< AK4recojets[0].Eta()   << endl;
     }
     if(widejetsReco.size()>=2){
       cout<< "pTWJ_recoj2: "<< AK4recojets[1].Pt()   << endl;
       cout<< "etaWJ_recoj2: "<< AK4recojets[1].Eta()   << endl;
     }
     if(widejets.size()>=1){
       cout<< "pTWJ_j1: "<< AK4jets[0].Pt()   << endl;
       cout<< "etaWJ_j1: "<< AK4jets[0].Eta()   << endl;
     }
     if(widejets.size()>=2){
       cout<< "pTWJ_j2: "<< AK4jets[1].Pt()   << endl;
       cout<< "etaWJ_j2: "<< AK4jets[1].Eta()   << endl;
     }
     */

     // optional call to fill a skim with the full content of the input roottuple
     //if( passedCut("nJetFinal") ) fillSkimTree();
       

     // optional call to fill a skim with a subset of the variables defined in the cutFile (use flag SAVE)
     if( passedAllPreviousCuts("mjj") && passedCut("mjj") ) 
       {
	 fillReducedSkimTree();

	 // ===== Take a look at this =====
	 // //Example on how to investigate quickly the data
	 // if(getVariableValue("mjj")>4000)
	 //   {
	 //     //fast creation and filling of histograms
	 //     CreateAndFillUserTH1D("h_dphijj_mjjgt4000", 100, 0, 3.15, getVariableValue("deltaPHIjj"));
	 //     CreateAndFillUserTH1D("h_htak4_mjjgt4000", 1000, 0, 10000, getVariableValue("HTAK4reco"));
	 //     CreateAndFillUserTH1D("h_nvtx_mjjgt4000", 31, -0.5, 30.5, getVariableValue("nVtx"));
	 //   }

       }

     // ===== Example of mjj spectrum after HLT selection =====
     // if( passedAllPreviousCuts("mjj") )
     //   {
     // 	 if(getVariableValue("passHLT")>0)
     // 	   {
     // 	     //fast creation and filling of histograms
     // 	     CreateAndFillUserTH1D("h_mjj_passHLT", getHistoNBins("mjj"), getHistoMin("mjj"), getHistoMax("mjj"), getVariableValue("mjj"));
     // 	   }
     //   }

     // reject events that did not pass level 0 cuts
     //if( !passedCut("0") ) continue;
     // ......
     
     // reject events that did not pass level 1 cuts
     //if( !passedCut("1") ) continue;
     // ......

     // reject events that did not pass the full cut list
     //if( !passedCut("all") ) continue;
     // ......

     // if( widejets.size() >= 2) {
     //  h_nJetFinal->Fill(widejets.size());
     //  h_DijetMass->Fill(wdijet.M());
     //  h_pT1stJet->Fill(widejets[0].Pt());
     //  h_pT2ndJet->Fill(widejets[1].Pt());
     //  h_eta1stJet->Fill(widejets[0].Eta());
     //  h_eta2ndJet->Fill(widejets[1].Eta());
     // }
     ////////////////////// User's code ends here ///////////////////////

   } // End loop over events

   //////////write histos 
   //   h_mjj_NoTrigger -> Write();
   //   h_mjj_HLTpass_CaloJet40 -> Write();
   //   h_mjj_HLTpass_CaloJet40_L1HTT -> Write();
   //   h_mjj_HLTpass_CaloJet40_L1HTT_CaloHT250 -> Write();
 
   // h_nVtx->Write();
   // h_trueVtx->Write();
   // h_nJetFinal->Write();
   // h_pT1stJet->Write();
   // h_pT2ndJet->Write();
   // h_DijetMass->Write();
   // h_eta1stJet->Write();
   // h_eta2ndJet->Write();

   // //pT of both jets, to be built using the histograms produced automatically by baseClass
   // TH1F * h_pTJets = new TH1F ("h_pTJets","", getHistoNBins("pT1stJet"), getHistoMin("pT1stJet"), getHistoMax("pT1stJet"));
   // h_pTJets->Add( & getHisto_noCuts_or_skim("pT1stJet") ); // all histos can be retrieved, see other getHisto_xxxx methods in baseClass.h
   // h_pTJets->Add( & getHisto_noCuts_or_skim("pT2ndJet") );
   // //one could also do:  *h_pTJets = getHisto_noCuts_or_skim("pT1stJet") + getHisto_noCuts_or_skim("pT2ndJet");
   // h_pTJets->Write();
   // //one could also do:   const TH1F& h = getHisto_noCuts_or_skim// and use h


   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
