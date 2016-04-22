#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TProfile.h>
#include "ptBinning.h"
#include "etaBinning.h"
#include "mJJBinning.h"
#include "vertexBinning.h"

bool verbose;
bool AK4Comparison = true;
bool WJComparison = true; 
double DeltaR_;
double alphacut_;
double DeltaEtaJJ_;
double MJJCut_;
EtaBinning mEtaBinning;
PtBinning mPtBinning; 
mJJBinning mMjjBinning; 

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile)
{
  std::cout << "analysisClass::analysisClass(): begins " << std::endl;

  //  std::string jetAlgo = getPreCutString1("jetAlgo");
  //  double rParam = getPreCutValue1("DeltaR");

  alphacut_ = getPreCutValue1("AlphaCut");
  DeltaR_ = getPreCutValue1("DeltaR");
  DeltaEtaJJ_ = getPreCutValue1("DeltaEtaJJ");
  MJJCut_ = getPreCutValue1("MJJCut");
  verbose = getPreCutValue1("verbose");


  //  if( jetAlgo == "AntiKt" )
  //    fjJetDefinition = JetDefPtr( new fastjet::JetDefinition(fastjet::antikt_algorithm, rParam) );
  //  else if( jetAlgo == "Kt" )
  //    fjJetDefinition = JetDefPtr( new fastjet::JetDefinition(fastjet::kt_algorithm, rParam) );
  //  else 
  //    fjJetDefinition = JetDefPtr( new fastjet::JetDefinition(fastjet::cambridge_algorithm, rParam) );
  
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
   double ptBins[] = {30, 100, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2500};                                                           
   int  binnum = sizeof(ptBins)/sizeof(double) -1;                                                                                                      
   TH1F* h_pTJet_Binned = new TH1F("pTJet_Binned","pTJet", binnum, ptBins);

   //   double etaBins[] = {-2.50, -1.50, -1.25, -1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 2.50};  
   double etaBins[] = {0.00, 0.50, 1.00, 1.50, 2.50};  
   int  binnumEta = sizeof(etaBins)/sizeof(double) -1;                                                                                                      

   /////////initialize variables

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) { // running over all events
     //  for (Long64_t jentry=0; jentry<5000;jentry++) { //runnig over few events
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     //std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     // if (Cut(ientry) < 0) continue;

     //     cout<<"Evento "<< jentry << endl;

     // uso eventi pari per estrarre le calibrazioni
     if( (int)jentry%2 ==! 0 ) continue;
     // uso eventi dispari per il closure test
     //     cout<< "Usato "<< jentry << endl;

     ////////////////////// User's code starts here ///////////////////////

     if( PassJSON){     

       if(verbose){
	 std::cout<< std::endl;
	 std::cout<<"Start"<< std::endl;
	 std::cout<<"PassJSON "<< PassJSON<<std::endl;
	 
	 std::cout<<"AK4"<< std::endl;
	 std::cout<<"Reco1"<< std::endl;
	 std::cout<<"Pt: "<< pTAK4_recoj1<<" eta: "<<etaAK4_recoj1<<" phi: "<<phiAK4_recoj1<< std::endl;
	 std::cout<<"Reco2"<< std::endl;
	 std::cout<<"Pt: "<< pTAK4_recoj2<<" eta: "<<etaAK4_recoj2<<" phi: "<<phiAK4_recoj2<< std::endl;
	 std::cout<<"Reco3"<< std::endl;
	 std::cout<<"Pt: "<< pTAK4_recoj3<<" eta: "<<etaAK4_recoj3<<" phi: "<<phiAK4_recoj3<< std::endl;
	 std::cout<<"HLT 1"<< std::endl;
	 std::cout<<"Pt: "<< pTAK4_j1<<" eta: "<<etaAK4_j1<<" phi: "<<phiAK4_j1<< std::endl;
	 std::cout<<"HLT2"<< std::endl;
	 std::cout<<"Pt: "<< pTAK4_j2<<" eta: "<<etaAK4_j2<<" phi: "<<phiAK4_j2<< std::endl;
	 std::cout<<"HLT3"<< std::endl;
	 std::cout<<"Pt: "<< pTAK4_j3<<" eta: "<<etaAK4_j3<<" phi: "<<phiAK4_j3<< std::endl;
       }

       CreateAndFillUserTH1D("mJJ_NoTrigger", 30, 0, 3000, mjj); 

       if( passHLT_ZeroBias_BtagSeq || passHLT_ZeroBias || passHLT_L1DoubleMu_BtagSeq || passHLT_L1DoubleMu || passHLT_CaloJet40_BtagSeq || passHLT_CaloJet40 || passHLT_L1HTT150_BtagSeq || passHLT_L1HTT150 || passHLT_CaloScoutingHT250 || passHLT_HT450_BtagSeq || passHLT_HT450){
	 CreateAndFillUserTH1D("mJJ_passTriggers", 30, 0, 3000, mjj); 
	 CreateAndFillUserTH1D("deltaEtaJJ_passTriggers", 100, 0, 3, deltaETAjj); 
       }
       // #1
       CreateAndFillUserTH1D("passHLT_ZeroBias_BtagSeq", 2 , -0.5, 1.5, passHLT_ZeroBias_BtagSeq); 
       CreateAndFillUserTH1D("passHLT_ZeroBias", 2 , -0.5, 1.5, passHLT_ZeroBias);
       if( passHLT_ZeroBias_BtagSeq || passHLT_ZeroBias ){
	 CreateAndFillUserTH1D("mJJ_passHLT_ZeroBias", 30, 0, 3000, mjj); 
	 CreateAndFillUserTH1D("deltaEtaJJ_passHLT_ZeroBias", 100, 0, 3, deltaETAjj); 
       }
       // #2       
       CreateAndFillUserTH1D("passHLT_L1DoubleMu_BtagSeq", 2 , -0.5, 1.5, passHLT_L1DoubleMu_BtagSeq);
       CreateAndFillUserTH1D("passHLT_L1DoubleMu", 2 , -0.5, 1.5, passHLT_L1DoubleMu);
       if( passHLT_L1DoubleMu_BtagSeq || passHLT_L1DoubleMu ){
	 CreateAndFillUserTH1D("mJJ_passHLT_L1DoubleMu", 30, 0, 3000, mjj); 
	 CreateAndFillUserTH1D("deltaEtaJJ_passHLT_L1DoubleMu", 100, 0, 3, deltaETAjj); 
       }
       // #3
       CreateAndFillUserTH1D("passHLT_CaloJet40_BtagSeq", 2 , -0.5, 1.5, passHLT_CaloJet40_BtagSeq);
       CreateAndFillUserTH1D("passHLT_CaloJet40", 2 , -0.5, 1.5, passHLT_CaloJet40);
       if( passHLT_CaloJet40_BtagSeq || passHLT_CaloJet40 ){
	 CreateAndFillUserTH1D("mJJ_passHLT_CaloJet40", 30, 0, 3000, mjj); 
	 CreateAndFillUserTH1D("deltaEtaJJ_passHLT_CaloJet40", 100, 0, 3, deltaETAjj); 
       }
       // #4
       CreateAndFillUserTH1D("passHLT_L1HTT150_BtagSeq", 2 , -0.5, 1.5, passHLT_L1HTT150_BtagSeq);
       CreateAndFillUserTH1D("passHLT_L1HTT150", 2 , -0.5, 1.5, passHLT_L1HTT150);
       if( passHLT_L1HTT150_BtagSeq || passHLT_L1HTT150 ){
	 CreateAndFillUserTH1D("mJJ_passHLT_L1HTT150", 30, 0, 3000, mjj); 
	 CreateAndFillUserTH1D("deltaEtaJJ_passHLT_L1HTT150", 100, 0, 3, deltaETAjj); 
       }
       // #5
       CreateAndFillUserTH1D("passHLT_CaloScoutingHT250", 2 , -0.5, 1.5,  passHLT_CaloScoutingHT250);
       if(  passHLT_CaloScoutingHT250 ){
	 CreateAndFillUserTH1D("mJJ_passHLT_CaloScoutingHT250", 30, 0, 3000, mjj); 
	 CreateAndFillUserTH1D("deltaEtaJJ_passHLT_CaloScoutingHT250", 100, 0, 3, deltaETAjj); 
       }
       // #6
       CreateAndFillUserTH1D("passHLT_HT450_BtagSeq", 2 , -0.5, 1.5, passHLT_HT450_BtagSeq);
       CreateAndFillUserTH1D("passHLT_HT450", 2 , -0.5, 1.5, passHLT_HT450);
       if( passHLT_HT450_BtagSeq || passHLT_HT450 ){
	 CreateAndFillUserTH1D("mJJ_passHLT_HT450", 30, 0, 3000, mjj); 
	 CreateAndFillUserTH1D("deltaEtaJJ_passHLT_HT450", 100, 0, 3, deltaETAjj); 
       }
       //////////////////// end of trigger study

       //       if( passHLT_ZeroBias_BtagSeq || passHLT_ZeroBias || passHLT_L1DoubleMu_BtagSeq || passHLT_L1DoubleMu || passHLT_CaloJet40_BtagSeq || passHLT_CaloJet40 || passHLT_L1HTT150_BtagSeq || passHLT_L1HTT150 || passHLT_CaloScoutingHT250 || passHLT_HT450_BtagSeq || passHLT_HT450 ){

       // trigger selected
       //   if( passHLT_ZeroBias_BtagSeq || passHLT_ZeroBias ){ // #1
       //   if( passHLT_L1DoubleMu_BtagSeq || passHLT_L1DoubleMu ){ // #2
       //   if( passHLT_CaloJet40_BtagSeq || passHLT_CaloJet40 ){ // #3
       //   if( passHLT_L1HTT150_BtagSeq || passHLT_L1HTT150 ){ // #4
       //   if( passHLT_CaloScoutingHT250 ){ // #5
       //   if( passHLT_HT450_BtagSeq || passHLT_HT450 ){ // #6

	 CreateAndFillUserTH1D("AK4_HTReco", 20 , 0, 2000, htAK4reco);
	 CreateAndFillUserTH1D("AK4_HTHLT", 20 , 0, 2000, htAK4);
	 CreateAndFillUserTH2D("AK4_HT", 20, 0, 2000, 20, 0, 2000, htAK4reco, htAK4);
	 
	 // qua posso aggiungere io taglio in nvertex per studiare la dipendenza dal pile up
	 // n< 5 , 5< n<15 , n>15
	 // HLT nVertex
	 //     if( nVtx > 5) continue; // Loose
	 //     if( nVtx <= 5 || nVtx >15) continue; //  Medium
	 //	 if( nVtx <= 15) continue; // Tight

	 if( AK4Comparison ){
	   // AK4 Jets RECO
	   vector<TLorentzVector> AK4recojets;
	   TLorentzVector ak4j1Reco, ak4j2Reco, ak4dijetReco;             
	   TLorentzVector ak4j3Reco;
	   vector<double> chargedHadEnFracReco ;
	   vector<double> neutrHadEnFracReco ;
	   vector<double> chargedElectromFracReco ;
	   vector<double> neutrElectromFracReco ;
	   vector<double> muEnFracReco ;
	   vector<double> chargedMultReco ;
	   vector<double> neutrMultReco ;
	   vector<double> photonMultReco ;
	   vector<double> jetIdTightReco ;
	   vector<double> jetCSVAK4Reco ;
	   // AK4 Jets HLT	 
	   vector<TLorentzVector> AK4jets;
	   TLorentzVector ak4j1, ak4j2, ak4dijet;      
	   TLorentzVector ak4j3;
	   vector<double> chargedHadEnFrac ;
	   vector<double> neutrHadEnFrac ;
	   vector<double> chargedElectromFrac ;
	   vector<double> neutrElectromFrac ;
	   vector<double> muEnFrac ;
	   vector<double> chargedMult ;
	   vector<double> neutrMult ;
	   vector<double> photonMult ;
	   vector<double> jetIdTight ;
	   vector<double> jetCSVAK4 ;

	   // Reco
	   ak4j1Reco.SetPtEtaPhiM(pTAK4_recoj1, etaAK4_recoj1, phiAK4_recoj1, massAK4_recoj1);
	   ak4j2Reco.SetPtEtaPhiM(pTAK4_recoj2, etaAK4_recoj2, phiAK4_recoj2, massAK4_recoj2);
	   ak4j3Reco.SetPtEtaPhiM(pTAK4_recoj3, etaAK4_recoj3, phiAK4_recoj3, massAK4_recoj3);
	   // Reco jet 1
	   AK4recojets.push_back( ak4j1Reco );
	   chargedHadEnFracReco.push_back(chargedHadEnFrac_recoj1) ;
	   neutrHadEnFracReco.push_back(neutrHadEnFrac_recoj1) ;
	   chargedElectromFracReco.push_back(chargedElectromFrac_recoj1) ;
	   neutrElectromFracReco.push_back(neutrElectromFrac_recoj1) ;
	   muEnFracReco.push_back(muEnFract_recoj1) ;
	   chargedMultReco.push_back(chargedMult_recoj1) ;
	   neutrMultReco.push_back(neutrMult_recoj1) ;
	   photonMultReco.push_back(photonMult_recoj1) ;
	   jetIdTightReco.push_back(IdTight_recoj1) ;
	   jetCSVAK4Reco.push_back(jetCSVAK4_recoj1) ;
	   // Reco jet 2
	   AK4recojets.push_back( ak4j2Reco );	 
	   chargedHadEnFracReco.push_back(chargedHadEnFrac_recoj2) ;
	   neutrHadEnFracReco.push_back(neutrHadEnFrac_recoj2) ;
	   chargedElectromFracReco.push_back(chargedElectromFrac_recoj2) ;
	   neutrElectromFracReco.push_back(neutrElectromFrac_recoj2) ;
	   muEnFracReco.push_back(muEnFract_recoj2) ;
	   chargedMultReco.push_back(chargedMult_recoj2) ;
	   neutrMultReco.push_back(neutrMult_recoj2) ;
	   photonMultReco.push_back(photonMult_recoj2) ;
	   jetIdTightReco.push_back(IdTight_recoj2) ;
	   jetCSVAK4Reco.push_back(jetCSVAK4_recoj2) ;
	   // Reco jet 3
	   AK4recojets.push_back( ak4j3Reco );

	   // HLT
	   ak4j1.SetPtEtaPhiM(pTAK4_j1, etaAK4_j1, phiAK4_j1, massAK4_j1);
	   ak4j2.SetPtEtaPhiM(pTAK4_j2, etaAK4_j2, phiAK4_j2, massAK4_j2);
	   ak4j3.SetPtEtaPhiM(pTAK4_j3, etaAK4_j3, phiAK4_j3, massAK4_j3);
	   // HLT jet 1
	   AK4jets.push_back( ak4j1 );
	   chargedHadEnFrac.push_back(chargedHadEnFrac_j1) ;
	   neutrHadEnFrac.push_back(neutrHadEnFrac_j1) ;
	   chargedElectromFrac.push_back(chargedElectromFrac_j1) ;
	   neutrElectromFrac.push_back(neutrElectromFrac_j1) ;
	   muEnFrac.push_back(muEnFract_j1) ;
	   chargedMult.push_back(chargedMult_j1) ;
	   neutrMult.push_back(neutrMult_j1) ;
	   photonMult.push_back(photonMult_j1) ;
	   jetIdTight.push_back(IdTight_j1) ;
	   jetCSVAK4.push_back(jetCSVAK4_j1) ;
	   // HLT jet 2	     
	   AK4jets.push_back( ak4j2 );
	   chargedHadEnFrac.push_back(chargedHadEnFrac_j2) ;
	   neutrHadEnFrac.push_back(neutrHadEnFrac_j2) ;
	   chargedElectromFrac.push_back(chargedElectromFrac_j2) ;
	   neutrElectromFrac.push_back(neutrElectromFrac_j2) ;
	   muEnFrac.push_back(muEnFract_j2) ;
	   chargedMult.push_back(chargedMult_j2) ;
	   neutrMult.push_back(neutrMult_j2) ;
	   photonMult.push_back(photonMult_j2) ;
	   jetIdTight.push_back(IdTight_j2) ;
	   jetCSVAK4.push_back(jetCSVAK4_j2) ;
	   // HLT jet 3
	   AK4jets.push_back( ak4j3 );

	   // Compare with AK4 Jets
	   //////////// matching	 
	   double DeltaReco1HLT1 = AK4recojets[0].DeltaR(AK4jets[0]);
	   double DeltaReco1HLT2 = AK4recojets[0].DeltaR(AK4jets[1]);
	   double DeltaReco2HLT1 = AK4recojets[1].DeltaR(AK4jets[0]);
	   double DeltaReco2HLT2 = AK4recojets[1].DeltaR(AK4jets[1]);
	     
	   if(verbose){
	     cout << "AK4: " << endl;
	     cout << "deltaR Reco1-HLT1 " <<  DeltaReco1HLT1 << endl;
	     cout << "deltaR Reco1-HLT2 " <<  DeltaReco1HLT2 << endl;
	     cout << "deltaR Reco2-HLT1 " <<  DeltaReco2HLT1 << endl;
	     cout << "deltaR Reco2-HLT2 " <<  DeltaReco2HLT2 << endl;
	   }       
	     
	   CreateAndFillUserTH1D("AK4_deltaR", 500, 0, 5, DeltaReco1HLT1 );
	   CreateAndFillUserTH1D("AK4_deltaR", 500, 0, 5, DeltaReco2HLT1 );

	   vector<int> IdxMatched;
	   vector<double> deltaR_min;
	   if( DeltaReco1HLT1 <= DeltaReco2HLT1 )
	     {
	       IdxMatched.push_back(0) ; // index Reco jets
	       IdxMatched.push_back(1) ; // index Reco jets
	       deltaR_min.push_back( DeltaReco1HLT1);
	       deltaR_min.push_back( DeltaReco2HLT2);
	     }else
	     {
	       IdxMatched.push_back(1) ; // index Reco jets
	       IdxMatched.push_back(0) ; // index Reco jets
	       deltaR_min.push_back( DeltaReco2HLT1);
	       deltaR_min.push_back( DeltaReco1HLT2);
	     }

	   for(int ii=0; ii<deltaR_min.size(); ii++){
	     if(verbose) cout << "deltaR_min " << deltaR_min.at(ii)  << endl;
	     CreateAndFillUserTH1D("AK4_deltaR_minimum", 500, 0, 5, deltaR_min.at(ii) );
	     
	     if(deltaR_min.at(ii) < DeltaR_ ){   
	       CreateAndFillUserTH1D("AK4_deltaR_minimum_AfterMatch", 500, 0, 5, deltaR_min.at(ii));
	       if(verbose){
		 std::cout<<"AK4"<<std::endl;
		 std::cout<<"First jet Reco"<<std::endl;
		 std::cout<<"Index= " << IdxMatched.at(ii) <<std::endl;
		 std::cout<<"Pt= "<< AK4recojets[IdxMatched.at(ii)].Pt()<<" Eta= "<<AK4recojets[IdxMatched.at(ii)].Eta()<< " Phi= "<<AK4recojets[IdxMatched.at(ii)].Phi()  <<std::endl;
		 std::cout<<"HLT" <<std::endl;
		 std::cout<<"Pt= "<< AK4jets[ii].Pt()<<" Eta= "<<AK4jets[ii].Eta()<< " Phi= "<<AK4jets[ii].Phi()  <<std::endl;
	       }

	       if(verbose)                                                                                                                                               
		 {                                                                                                                                                       
                   std::cout<<"pTJet= " << AK4jets[ii].Pt() << std::endl;                                                                                                               
                   if( passHLT_L1DoubleMu_BtagSeq || passHLT_L1DoubleMu ){                                                                                                   
                     cout<<"Trigger DoubleMu passato"<< endl;                                                                                                               
                   }                                                                                                                                                     
                   if( passHLT_HT450_BtagSeq || passHLT_HT450 ){                                                                                                          
                     cout<<"Trigger HT450 passato"<< endl;                                                                                                               
                   }
		 }
	       
	       if(AK4jets[ii].Pt() <= 300 && !passHLT_L1DoubleMu_BtagSeq && !passHLT_L1DoubleMu) continue; 
	       if(AK4jets[ii].Pt() > 300 && !passHLT_HT450_BtagSeq && !passHLT_HT450 ) continue; 
	       
	       if(verbose) cout<<"Evento Passato"<< endl;
	       
	       h_pTJet_Binned -> Fill( AK4jets[ii].Pt() );	       

	       // Reco variables
	       CreateAndFillUserTH1D("AK4_pT_JetReco", 150, 0, 3000, AK4recojets[IdxMatched.at(ii)].Pt() );
	       CreateAndFillUserTH1D("AK4_Eta_JetReco", 20, -2.5, 2.5, AK4recojets[IdxMatched.at(ii)].Eta() );
	       CreateAndFillUserTH1D("AK4_Phi_JetReco", 100, -3.15, 3.15, AK4recojets[IdxMatched.at(ii)].Phi() );
	       CreateAndFillUserTH1D("AK4_chargedHadEnFracReco", 100, 0, 1, chargedHadEnFracReco[IdxMatched.at(ii)] );
	       CreateAndFillUserTH1D("AK4_neutrHadEnFracReco", 100, 0, 1, neutrHadEnFracReco[IdxMatched.at(ii)] );
	       CreateAndFillUserTH1D("AK4_chargedElectromFracReco", 100, 0, 1, chargedElectromFracReco[IdxMatched.at(ii)] );
	       CreateAndFillUserTH1D("AK4_neutrElectromFracReco", 100, 0, 1, neutrElectromFracReco[IdxMatched.at(ii)] );
	       CreateAndFillUserTH1D("AK4_muEnFracReco", 100, 0, 1, muEnFracReco[IdxMatched.at(ii)] );
	       CreateAndFillUserTH1D("AK4_chargedMultReco", 100, 0, 100, chargedMultReco[IdxMatched.at(ii)] );
	       CreateAndFillUserTH1D("AK4_neutrMultReco", 100, 0, 100, neutrMultReco[IdxMatched.at(ii)] );
	       CreateAndFillUserTH1D("AK4_photonMultReco", 100, 0, 100, photonMultReco[IdxMatched.at(ii)] );
	       CreateAndFillUserTH1D("AK4_jetCSVAK4Reco", 100, 0, 2, jetCSVAK4Reco[IdxMatched.at(ii)] );
	       // HLT variables
	       CreateAndFillUserTH1D("AK4_pT_JetHLT", 150, 0, 3000, AK4jets[ii].Pt() );
	       CreateAndFillUserTH1D("AK4_Eta_JetHLT", 20, -2.5, 2.5, AK4jets[ii].Eta() );
	       CreateAndFillUserTH1D("AK4_Phi_JetHLT", 100, -3.15, 3.15, AK4jets[ii].Phi() );
	       CreateAndFillUserTH1D("AK4_chargedHadEnFracHLT", 100, 0, 1, chargedHadEnFrac[ii] );
	       CreateAndFillUserTH1D("AK4_neutrHadEnFracHLT", 100, 0, 1, neutrHadEnFrac[ii] );
	       CreateAndFillUserTH1D("AK4_chargedElectromFracHLT", 100, 0, 1, chargedElectromFrac[ii] );
	       CreateAndFillUserTH1D("AK4_neutrElectromFracHLT", 100, 0, 1, neutrElectromFrac[ii] );
	       CreateAndFillUserTH1D("AK4_muEnFracHLT", 100, 0, 1, muEnFrac[ii] );
	       CreateAndFillUserTH1D("AK4_chargedMultHLT", 100, 0, 100, chargedMult[ii] );
	       CreateAndFillUserTH1D("AK4_neutrMultHLT", 100, 0, 100, neutrMult[ii] );
	       CreateAndFillUserTH1D("AK4_photonMultHLT", 100, 0, 100, photonMult[ii] );
	       CreateAndFillUserTH1D("AK4_jetCSVAK4HLT", 100, 0, 2, jetCSVAK4[ii] );
	       // Reco vs HLT
	       CreateAndFillUserTH2D("AK4_pT_Jet", 150, 0, 3000, 150, 0, 3000, AK4recojets[IdxMatched.at(ii)].Pt(), AK4jets[ii].Pt() );
	       CreateAndFillUserTH2D("AK4_Eta_Jet", 20, -2.5, 2.5, 20, -2.5, 2.5, AK4recojets[IdxMatched.at(ii)].Eta(), AK4jets[ii].Eta() );
	       CreateAndFillUserTH2D("AK4_Phi_Jet", 100, -3.15, 3.15, 100, -3.15, 3.15, AK4recojets[IdxMatched.at(ii)].Phi(), AK4jets[ii].Phi() );
	       CreateAndFillUserTH2D("AK4_chargedHadEnFrac", 100, 0, 1, 100, 0, 1, chargedHadEnFracReco[IdxMatched.at(ii)], chargedHadEnFrac[ii] );
	       CreateAndFillUserTH2D("AK4_neutrHadEnFrac", 100, 0, 1, 100, 0, 1, neutrHadEnFracReco[IdxMatched.at(ii)], neutrHadEnFrac[ii] );
	       CreateAndFillUserTH2D("AK4_chargedElectromFrac", 100, 0, 1, 100, 0, 1, chargedElectromFracReco[IdxMatched.at(ii)], chargedElectromFrac[ii] );
	       CreateAndFillUserTH2D("AK4_neutrElectromFrac", 100, 0, 1, 100, 0, 1, neutrElectromFracReco[IdxMatched.at(ii)], neutrElectromFrac[ii] );
	       CreateAndFillUserTH2D("AK4_muEnFrac", 100, 0, 1, 100, 0, 1, muEnFracReco[IdxMatched.at(ii)], muEnFrac[ii] );
	       CreateAndFillUserTH2D("AK4_chargedMult", 100, 0, 100, 100, 0, 100, chargedMultReco[IdxMatched.at(ii)], chargedMult[ii] );
	       CreateAndFillUserTH2D("AK4_neutrMult", 100, 0, 100, 100, 0, 100, neutrMultReco[IdxMatched.at(ii)], neutrMult[ii] );
	       CreateAndFillUserTH2D("AK4_photonMult", 100, 0, 100, 100, 0, 100, photonMultReco[IdxMatched.at(ii)], photonMult[ii] );
	       CreateAndFillUserTH2D("AK4_jetCSVAK4", 100, 0, 2, 100, 0, 2, jetCSVAK4Reco[IdxMatched.at(ii)], jetCSVAK4[ii] );
	     
	       // Bias in pT AK4
	       // pT(Reco) - pT(HLT) / pT(HLT)	   
	       double pTBias = (AK4recojets[IdxMatched.at(ii)].Pt() - AK4jets[ii].Pt() ) / AK4jets[ii].Pt();
	       if(verbose) std::cout<<"pT Bias= "<<pTBias << std::endl;

	       double CHFBias = ( chargedHadEnFracReco[IdxMatched.at(ii)] - chargedHadEnFrac[ii] ) / chargedHadEnFrac[ii];
	       double NHFBias = ( neutrHadEnFracReco[IdxMatched.at(ii)] - neutrHadEnFrac[ii] ) / neutrHadEnFrac[ii];
	       double CEFBias = ( chargedElectromFracReco[IdxMatched.at(ii)] - chargedElectromFrac[ii] ) / chargedElectromFrac[ii];	     
	       double NEFBias = ( neutrElectromFracReco[IdxMatched.at(ii)] - neutrElectromFrac[ii] ) / neutrElectromFrac[ii];	     
	       double MuFBias = ( muEnFracReco[IdxMatched.at(ii)] - muEnFrac[ii] ) / muEnFrac[ii];	     

	       int AK4_etaBin = mEtaBinning.getBin(  fabs(AK4jets[ii].Eta())  ); // modulo
	       int AK4Reco_etaBin = mEtaBinning.getBin(  fabs(AK4recojets[IdxMatched.at(ii)].Eta())  ); // modulo
	       // int AK4_etaBin = mEtaBinning.getBin( AK4jets[ii].Eta() ); //separo negativi da positivi
	       int AK4_ptBin   = mPtBinning.getPtBin( AK4jets[ii].Pt() );
	       int AK4Reco_ptBin   = mPtBinning.getPtBin( AK4recojets[IdxMatched.at(ii)].Pt() );
	       std::pair<float, float> AK4_ptBins = mPtBinning.getBinValue( AK4_ptBin );
	       std::pair<float, float> AK4Reco_ptBins = mPtBinning.getBinValue( AK4Reco_ptBin );
	       std::string AK4_etaName = mEtaBinning.getBinName( AK4_etaBin );
	       std::string AK4Reco_etaName = mEtaBinning.getBinName( AK4Reco_etaBin );
	       std::string AK4_EtaBin = TString::Format("%s", AK4_etaName.c_str() ).Data();
	       std::string AK4Reco_EtaBin = TString::Format("%s", AK4Reco_etaName.c_str() ).Data();
	       std::string AK4_PtBin = TString::Format("pT_%i_%i", (int) AK4_ptBins.first, (int) AK4_ptBins.second ).Data();
	       std::string AK4Reco_PtBin = TString::Format("pT_%i_%i", (int) AK4Reco_ptBins.first, (int) AK4Reco_ptBins.second ).Data();
	       std::string AK4_Bin = TString::Format("%s_pT_%i_%i", AK4_etaName.c_str(), (int) AK4_ptBins.first, (int) AK4_ptBins.second ).Data(); 

	       if( verbose ){	       
		 cout<<"AK4: Eta =   "<< AK4jets[ii].Eta()  << " pT = "<< AK4jets[ii].Pt() <<endl;
		 cout<<"EtaBin =   "<< AK4_etaBin <<endl;
		 std::pair<float, float> AK4_etaBins = mEtaBinning.getBinValue( AK4_etaBin );
		 cout<<"EtaBin.first  "<< AK4_etaBins.first << "   EtaBin.second   "<< AK4_etaBins.second<<endl;
		 cout<<"pTBin =   " << AK4_ptBin  <<endl;
		 cout<<"ptBin.first  "  << AK4_ptBins.first  << "    ptBin.second   "<< AK4_ptBins.second<<endl;
		 cout<<"pTBin Reco =   " << AK4Reco_ptBin  <<endl;
		 cout<<"ptBinReco.first  "  << AK4Reco_ptBins.first  << "    ptBinReco.second   "<< AK4Reco_ptBins.second<<endl;
		 cout<<"AK4 Eta Bin  "  << AK4_EtaBin.c_str() <<endl;
		 cout<<"AK4 Pt Bin  "  << AK4_PtBin.c_str() <<endl;
		 cout<<"AK4 Reco Pt Bin  "  << AK4Reco_PtBin.c_str() <<endl;
	       }

	       CreateAndFillUserTH1D("AK4_pTBias", 1200, -2, 2, pTBias );
	       CreateAndFillUserTH1D("AK4_CHFBias", 1200, -2., 2., CHFBias );
	       CreateAndFillUserTH1D("AK4_NHFBias", 1200, -2., 2., NHFBias );
	       CreateAndFillUserTH1D("AK4_CEFBias", 1200, -2., 2.,  CEFBias );
	       CreateAndFillUserTH1D("AK4_NEFBias", 1200, -2., 2.,  NEFBias );
	       CreateAndFillUserTH1D("AK4_MuFBias", 1200, -2., 2., MuFBias );	       
	       // devo guardare questi e... --> fit per prendere il picco
	       CreateAndFillUserTH1D( Form("AK4_pTBias_%s", AK4_Bin.c_str()), 1200, -2., 2., pTBias );
	       CreateAndFillUserTH1D( Form("AK4_CHFBias_%s", AK4_Bin.c_str()), 1200, -2., 2., CHFBias );
	       CreateAndFillUserTH1D( Form("AK4_NHFBias_%s", AK4_Bin.c_str()), 1200, -2., 2., NHFBias );
	       CreateAndFillUserTH1D( Form("AK4_CEFBias_%s", AK4_Bin.c_str()), 1200, -2., 2.,  CEFBias );
	       CreateAndFillUserTH1D( Form("AK4_NEFBias_%s", AK4_Bin.c_str()), 1200, -2., 2.,  NEFBias );
	       CreateAndFillUserTH1D( Form("AK4_MuFBias_%s", AK4_Bin.c_str()), 1200, -2., 2., MuFBias );	       

	       //////////// energy fraction vs Eta
	       // charged hadron fraction
	       CreateAndFillUserTProfile("AK4_pTBias_vs_CHF",  100,  0., 1., -2., 2.,  chargedHadEnFrac[ii], pTBias );
	       CreateAndFillUserTProfile( Form("AK4_pTBias_vs_CHF_%s", AK4_PtBin.c_str()),  100,  0., 1., -2., 2.,  chargedHadEnFrac[ii], pTBias );
	       //vs eta
	       CreateAndFillUserTProfile("AK4_CHF_vs_Eta_HLT",  20, -2.5, 2.5, 0., 1.,  AK4jets[ii].Eta(), chargedHadEnFrac[ii] );  
	       CreateAndFillUserTProfile("AK4_CHF_vs_Eta_Reco",  20, -2.5, 2.5, 0., 1., AK4recojets[IdxMatched.at(ii)].Eta(), chargedHadEnFracReco[IdxMatched.at(ii)] );
	       CreateAndFillUserTProfile( Form("AK4_CHF_vs_Eta_%s_HLT", AK4_PtBin.c_str()),  20, -2.5, 2.5, 0., 1., AK4jets[ii].Eta(), chargedHadEnFrac[ii] );
	       CreateAndFillUserTProfile( Form("AK4_CHF_vs_Eta_%s_Reco", AK4Reco_PtBin.c_str()),  20, -2.5, 2.5, 0., 1., AK4recojets[IdxMatched.at(ii)].Eta(), chargedHadEnFracReco[IdxMatched.at(ii)] );
	       // vs pT
	       CreateAndFillUserTProfile("AK4_CHF_vs_Pt_HLT",  binnum, ptBins, 0., 1., AK4jets[ii].Pt(), chargedHadEnFrac[ii]);
	       CreateAndFillUserTProfile("AK4_CHF_vs_Pt_Reco", binnum, ptBins, 0., 1., AK4recojets[IdxMatched.at(ii)].Pt(), chargedHadEnFracReco[IdxMatched.at(ii)]);
	       CreateAndFillUserTProfile( Form("AK4_CHF_vs_Pt_%s_HLT", AK4_EtaBin.c_str()),  binnum, ptBins, 0., 1., AK4jets[ii].Pt(), chargedHadEnFrac[ii] );	      
	       CreateAndFillUserTProfile( Form("AK4_CHF_vs_Pt_%s_Reco", AK4Reco_EtaBin.c_str()), binnum, ptBins, 0., 1., AK4recojets[IdxMatched.at(ii)].Pt(), chargedHadEnFracReco[IdxMatched.at(ii)]  );

	       // neutral hadron fraction
	       CreateAndFillUserTProfile("AK4_pTBias_vs_NHF", 100,  0., 1., -2., 2.,  neutrHadEnFrac[ii], pTBias );
	       CreateAndFillUserTProfile( Form("AK4_pTBias_vs_NHF_%s", AK4_PtBin.c_str()),  100,  0., 1., -2., 2.,  neutrHadEnFrac[ii], pTBias );
	       //vs eta
	       CreateAndFillUserTProfile("AK4_NHF_vs_Eta_HLT",  20, -2.5, 2.5, 0., 1.,  AK4jets[ii].Eta(), neutrHadEnFrac[ii] );   
	       CreateAndFillUserTProfile("AK4_NHF_vs_Eta_Reco",  20, -2.5, 2.5, 0., 1., AK4recojets[IdxMatched.at(ii)].Eta(), neutrHadEnFracReco[IdxMatched.at(ii)] );
	       CreateAndFillUserTProfile( Form("AK4_NHF_vs_Eta_%s_HLT", AK4_PtBin.c_str()),  20, -2.5, 2.5, 0., 1., AK4jets[ii].Eta(), neutrHadEnFrac[ii] ); 
	       CreateAndFillUserTProfile( Form("AK4_NHF_vs_Eta_%s_Reco", AK4Reco_PtBin.c_str()),  20, -2.5, 2.5, 0., 1., AK4recojets[IdxMatched.at(ii)].Eta(), neutrHadEnFracReco[IdxMatched.at(ii)] );
	       // vs pT
	       CreateAndFillUserTProfile("AK4_NHF_vs_Pt_HLT",  binnum, ptBins, 0., 1., AK4jets[ii].Pt(), neutrHadEnFrac[ii]);
	       CreateAndFillUserTProfile("AK4_NHF_vs_Pt_Reco", binnum, ptBins, 0., 1., AK4recojets[IdxMatched.at(ii)].Pt(), neutrHadEnFracReco[IdxMatched.at(ii)]);
	       CreateAndFillUserTProfile( Form("AK4_NHF_vs_Pt_%s_HLT", AK4_EtaBin.c_str()),  binnum, ptBins, 0., 1., AK4jets[ii].Pt(), neutrHadEnFrac[ii] );	      
	       CreateAndFillUserTProfile( Form("AK4_NHF_vs_Pt_%s_Reco", AK4Reco_EtaBin.c_str()), binnum, ptBins, 0., 1., AK4recojets[IdxMatched.at(ii)].Pt(), neutrHadEnFracReco[IdxMatched.at(ii)]  );

	       // charged electrom fraction
	       CreateAndFillUserTProfile("AK4_pTBias_vs_CEF",  100,  0., 1., -2., 2.,  chargedElectromFrac[ii], pTBias );
	       CreateAndFillUserTProfile( Form("AK4_pTBias_vs_CEF_%s", AK4_PtBin.c_str()),  100,  0., 1., -2., 2.,  chargedElectromFrac[ii], pTBias );
	       //vs eta
	       CreateAndFillUserTProfile("AK4_CEF_vs_Eta_HLT",  20, -2.5, 2.5, 0., 1.,  AK4jets[ii].Eta(), chargedElectromFrac[ii] );  
	       CreateAndFillUserTProfile("AK4_CEF_vs_Eta_Reco",  20, -2.5, 2.5, 0., 1., AK4recojets[IdxMatched.at(ii)].Eta(), chargedElectromFracReco[IdxMatched.at(ii)] );
	       CreateAndFillUserTProfile( Form("AK4_CEF_vs_Eta_%s_HLT", AK4_PtBin.c_str()),  20, -2.5, 2.5, 0., 1., AK4jets[ii].Eta(), chargedElectromFrac[ii] );
	       CreateAndFillUserTProfile( Form("AK4_CEF_vs_Eta_%s_Reco", AK4Reco_PtBin.c_str()),  20, -2.5, 2.5, 0., 1., AK4recojets[IdxMatched.at(ii)].Eta(), chargedElectromFracReco[IdxMatched.at(ii)] );
	       // vs pT
	       CreateAndFillUserTProfile("AK4_CEF_vs_Pt_HLT",  binnum, ptBins, 0., 1., AK4jets[ii].Pt(), chargedElectromFrac[ii]);
	       CreateAndFillUserTProfile("AK4_CEF_vs_Pt_Reco", binnum, ptBins, 0., 1., AK4recojets[IdxMatched.at(ii)].Pt(), chargedElectromFracReco[IdxMatched.at(ii)]);
	       CreateAndFillUserTProfile( Form("AK4_CEF_vs_Pt_%s_HLT", AK4_EtaBin.c_str()),  binnum, ptBins, 0., 1., AK4jets[ii].Pt(), chargedElectromFrac[ii] );	      
	       CreateAndFillUserTProfile( Form("AK4_CEF_vs_Pt_%s_Reco", AK4Reco_EtaBin.c_str()), binnum, ptBins, 0., 1., AK4recojets[IdxMatched.at(ii)].Pt(), chargedElectromFracReco[IdxMatched.at(ii)]  );

	       // neutral electrom fraction
	       CreateAndFillUserTProfile("AK4_pTBias_vs_NEF",  100, 0., 1., -2., 2.,  neutrElectromFrac[ii], pTBias );
	       CreateAndFillUserTProfile( Form("AK4_pTBias_vs_NEF_%s", AK4_PtBin.c_str()),  100,  0., 1., -2., 2.,  neutrElectromFrac[ii], pTBias );
	       //vs eta
	       CreateAndFillUserTProfile("AK4_NEF_vs_Eta_HLT",  20, -2.5, 2.5, 0., 1.,  AK4jets[ii].Eta(), neutrElectromFrac[ii] );  
	       CreateAndFillUserTProfile("AK4_NEF_vs_Eta_Reco",  20, -2.5, 2.5, 0., 1., AK4recojets[IdxMatched.at(ii)].Eta(), neutrElectromFracReco[IdxMatched.at(ii)] );
	       CreateAndFillUserTProfile( Form("AK4_NEF_vs_Eta_%s_HLT", AK4_PtBin.c_str()),  20, -2.5, 2.5, 0., 1., AK4jets[ii].Eta(), neutrElectromFrac[ii] );
	       CreateAndFillUserTProfile( Form("AK4_NEF_vs_Eta_%s_Reco", AK4Reco_PtBin.c_str()),  20, -2.5, 2.5, 0., 1., AK4recojets[IdxMatched.at(ii)].Eta(), neutrElectromFracReco[IdxMatched.at(ii)] );
	       // vs pT
	       CreateAndFillUserTProfile("AK4_NEF_vs_Pt_HLT",  binnum, ptBins, 0., 1., AK4jets[ii].Pt(), neutrElectromFrac[ii]);
	       CreateAndFillUserTProfile("AK4_NEF_vs_Pt_Reco", binnum, ptBins, 0., 1., AK4recojets[IdxMatched.at(ii)].Pt(), neutrElectromFracReco[IdxMatched.at(ii)]);
	       CreateAndFillUserTProfile( Form("AK4_NEF_vs_Pt_%s_HLT", AK4_EtaBin.c_str()),  binnum, ptBins, 0., 1., AK4jets[ii].Pt(), neutrElectromFrac[ii] );	      
	       CreateAndFillUserTProfile( Form("AK4_NEF_vs_Pt_%s_Reco", AK4Reco_EtaBin.c_str()), binnum, ptBins, 0., 1., AK4recojets[IdxMatched.at(ii)].Pt(), neutrElectromFracReco[IdxMatched.at(ii)]  );

	       // mu fraction
	       CreateAndFillUserTProfile("AK4_pTBias_vs_MuF", 100,  0., 1., -2., 2.,  muEnFrac[ii], pTBias );
	       CreateAndFillUserTProfile( Form("AK4_pTBias_vs_MuF_%s", AK4_PtBin.c_str()),  100,  0., 1., -2., 2.,  muEnFrac[ii], pTBias );
	       //vs eta
	       CreateAndFillUserTProfile("AK4_MuF_vs_Eta_HLT",  20, -2.5, 2.5, 0., 1.,  AK4jets[ii].Eta(), muEnFrac[ii] );  
	       CreateAndFillUserTProfile("AK4_MuF_vs_Eta_Reco",  20, -2.5, 2.5, 0., 1., AK4recojets[IdxMatched.at(ii)].Eta(), muEnFracReco[IdxMatched.at(ii)] );
	       CreateAndFillUserTProfile( Form("AK4_MuF_vs_Eta_%s_HLT", AK4_PtBin.c_str()),  20, -2.5, 2.5, 0., 1., AK4jets[ii].Eta(), muEnFrac[ii] ); 
	       CreateAndFillUserTProfile( Form("AK4_MuF_vs_Eta_%s_Reco", AK4Reco_PtBin.c_str()),  20, -2.5, 2.5, 0., 1., AK4recojets[IdxMatched.at(ii)].Eta(), muEnFracReco[IdxMatched.at(ii)] );
	       // vs pT
	       CreateAndFillUserTProfile("AK4_MuF_vs_Pt_HLT",  binnum, ptBins, 0., 1., AK4jets[ii].Pt(), muEnFrac[ii]);
	       CreateAndFillUserTProfile("AK4_MuF_vs_Pt_Reco", binnum, ptBins, 0., 1., AK4recojets[IdxMatched.at(ii)].Pt(), muEnFracReco[IdxMatched.at(ii)]);
	       CreateAndFillUserTProfile( Form("AK4_MuF_vs_Pt_%s_HLT", AK4_EtaBin.c_str()),  binnum, ptBins, 0., 1., AK4jets[ii].Pt(), muEnFrac[ii] );	      
	       CreateAndFillUserTProfile( Form("AK4_MuF_vs_Pt_%s_Reco", AK4Reco_EtaBin.c_str()), binnum, ptBins, 0., 1., AK4recojets[IdxMatched.at(ii)].Pt(), muEnFracReco[IdxMatched.at(ii)]  );

	       // bias vs HLT quantity
	       ////// TProfile: pTBias vs nVtx
	       CreateAndFillUserTProfile("AK4_pTBias_vs_Nvtx",  50, 0., 50., -2., 2.,      nVtx, pTBias );
	       CreateAndFillUserTProfile("AK4_pT_vs_NvtxReco", 50, 0., 50., 0., 5000., nVtxreco, AK4recojets[IdxMatched.at(ii)].Pt() );
	       CreateAndFillUserTProfile("AK4_pT_vs_NvtxHLT",  50, 0., 50., 0., 5000., nVtx, AK4jets[ii].Pt() );
	       CreateAndFillUserTProfile( Form("AK4_pTBias_vs_Nvtx_%s", AK4_Bin.c_str()),  50, 0., 50., -2., 2.,      nVtx, pTBias );
	       CreateAndFillUserTProfile( Form("AK4_pT_vs_NvtxReco_%s", AK4_Bin.c_str()), 50, 0., 50., 0., 5000., nVtxreco, AK4recojets[IdxMatched.at(ii)].Pt() );
	       CreateAndFillUserTProfile( Form("AK4_pT_vs_NvtxHLT_%s", AK4_Bin.c_str()),  50, 0., 50., 0., 5000., nVtx, AK4jets[ii].Pt() );
	       ////// TProfile: pTBias vs Eta
	       CreateAndFillUserTProfile("AK4_pTBias_vs_Eta",  20, -2.5, 2.5, -2., 2.,      AK4jets[ii].Eta(), pTBias );
	       CreateAndFillUserTProfile("AK4_pT_vs_EtaReco", 20, -2.5, 2.5, 0., 5000., AK4recojets[IdxMatched.at(ii)].Eta(), AK4recojets[IdxMatched.at(ii)].Pt() );
	       CreateAndFillUserTProfile("AK4_pT_vs_EtaHLT",  20, -2.5, 2.5, 0., 5000., AK4jets[ii].Eta(), AK4jets[ii].Pt() );
	       CreateAndFillUserTProfile( Form("AK4_pTBias_vs_Eta_%s", AK4_PtBin.c_str()),  20, -2.5, 2.5, -2., 2.,      AK4jets[ii].Eta(), pTBias );
	       CreateAndFillUserTProfile( Form("AK4_pT_vs_Eta_%s_Reco", AK4Reco_PtBin.c_str()), 20, -2.5, 2.5, 0., 5000., AK4recojets[IdxMatched.at(ii)].Eta(), AK4recojets[IdxMatched.at(ii)].Pt() );
	       CreateAndFillUserTProfile( Form("AK4_pT_vs_Eta_%s_HLT", AK4_PtBin.c_str()),  20, -2.5, 2.5, 0., 5000., AK4jets[ii].Eta(), AK4jets[ii].Pt() );
	       // ...e questo che e' quello importante -- adesso sto usando il bias 1D
	       CreateAndFillUserTProfile( Form("AK4_pTBias_vs_Eta_VariableBin_%s", AK4_PtBin.c_str()),  binnumEta, etaBins, -2., 2., AK4jets[ii].Eta(), pTBias );

	       ////// TProfile: pTBias vs Pt
	       CreateAndFillUserTProfile("AK4_pTBias_vs_Pt",  binnum, ptBins, -2., 2.,      AK4jets[ii].Pt(), pTBias );
	       CreateAndFillUserTProfile("AK4_pT_vs_PtReco",  binnum, ptBins, 0., 5000., AK4recojets[IdxMatched.at(ii)].Pt(), AK4recojets[IdxMatched.at(ii)].Pt() );
	       CreateAndFillUserTProfile("AK4_pT_vs_PtHLT",   binnum, ptBins, 0., 5000., AK4jets[ii].Pt(), AK4jets[ii].Pt() );	      
	       CreateAndFillUserTProfile( Form("AK4_pTBias_vs_Pt_%s", AK4_EtaBin.c_str()),   binnum, ptBins, -2., 2.,      AK4jets[ii].Pt(), pTBias );
	       CreateAndFillUserTProfile( Form("AK4_pT_vs_Pt_%s_Reco", AK4Reco_EtaBin.c_str()),  binnum, ptBins, 0., 5000., AK4recojets[IdxMatched.at(ii)].Pt(), AK4recojets[IdxMatched.at(ii)].Pt() );
	       CreateAndFillUserTProfile( Form("AK4_pT_vs_Pt_%s_HLT", AK4_EtaBin.c_str()),   binnum, ptBins, 0., 5000., AK4jets[ii].Pt(), AK4jets[ii].Pt() );	      
	       
	       if( fabs(AK4jets[ii].Eta())<=1.5 )
		 {
		   CreateAndFillUserTProfile("AK4_pTBias_vs_Pt_Barrel",  50, 0., 5000., -2., 2., AK4jets[ii].Pt(), pTBias );		 
		 }else if( fabs(AK4jets[ii].Eta())>1.5 && fabs(AK4jets[ii].Eta())<=2.5 )
		 {
		   CreateAndFillUserTProfile("AK4_pTBias_vs_Pt_Endcap",  50, 0., 5000., -2., 2., AK4jets[ii].Pt(), pTBias );		 		 
		 }
	       
	     }//if matching
	   }// loop over jets
	   ////// End AK4 analysis	     
	 } // AK4 jet comparison
	   
	   /////////////////////////////////////////// WIDEJETS
	  
	 if(WJComparison){

	   if(verbose){
	     std::cout<<"WIDEJET"<< std::endl;
	     std::cout<<"Reco 1"<< std::endl;
	     std::cout<<"Pt: "<<pTWJ_recoj1<<" eta: "<<etaWJ_recoj1<<" phi: "<<phiWJ_recoj1<< std::endl;
	     std::cout<<"Reco 2"<< std::endl;
	     std::cout<<"Pt: "<<pTWJ_recoj2<<" eta: "<<etaWJ_recoj2<<" phi: "<<phiWJ_recoj2<< std::endl;
	     std::cout<<"HLT 1"<< std::endl;
	     std::cout<<"Pt: "<<pTWJ_j1<<" eta: "<<etaWJ_j1<<" phi: "<<phiWJ_j1<< std::endl;
	     std::cout<<"HLT 2"<< std::endl;
	     std::cout<<"Pt: "<<pTWJ_j2<<" eta: "<<etaWJ_j2<<" phi: "<<phiWJ_j2<< std::endl;
	   }
	     
	   vector<TLorentzVector> widejetsReco;
	   TLorentzVector wj1Reco, wj2Reco, wdijetReco; 	     
	   vector<TLorentzVector> widejets;
	   TLorentzVector wj1, wj2, wdijet; 
	     
	   //Reco
	   wj1Reco.SetPtEtaPhiM(pTWJ_recoj1, etaWJ_recoj1, phiWJ_recoj1, massWJ_recoj1);
	   wj2Reco.SetPtEtaPhiM(pTWJ_recoj2, etaWJ_recoj2, phiWJ_recoj2, massWJ_recoj2);	     
	   widejetsReco.push_back( wj1Reco );
	   widejetsReco.push_back( wj2Reco );
	     
	   //HLT
	   wj1.SetPtEtaPhiM(pTWJ_j1, etaWJ_j1, phiWJ_j1, massWJ_j1);
	   wj2.SetPtEtaPhiM(pTWJ_j2, etaWJ_j2, phiWJ_j2, massWJ_j2);
	   widejets.push_back( wj1 );
	   widejets.push_back( wj2 );

	   // Compare with WIDEJETS
	   //////////// matching	 
	   double DeltaReco1HLT1 = widejetsReco[0].DeltaR(widejets[0]);
	   double DeltaReco1HLT2 = widejetsReco[0].DeltaR(widejets[1]);
	   double DeltaReco2HLT1 = widejetsReco[1].DeltaR(widejets[0]);
	   double DeltaReco2HLT2 = widejetsReco[1].DeltaR(widejets[1]);
	     
	   if(verbose){
	     cout << "WIDEJETS " << endl;
	     cout << "deltaR Reco1-HLT1 " <<  DeltaReco1HLT1 << endl;
	     cout << "deltaR Reco1-HLT2 " <<  DeltaReco1HLT2 << endl;
	     cout << "deltaR Reco2-HLT1 " <<  DeltaReco2HLT1 << endl;
	     cout << "deltaR Reco2-HLT2 " <<  DeltaReco2HLT2 << endl;
	   }       

	   CreateAndFillUserTH1D("deltaR", 500, 0, 5, DeltaReco1HLT1 );
	   CreateAndFillUserTH1D("deltaR", 500, 0, 5, DeltaReco2HLT1 );

	   vector<int> IdxMatched;
	   vector<double> deltaR_min;
	   if( DeltaReco1HLT1 <= DeltaReco2HLT1 )
	     {
	       IdxMatched.push_back( 0 ) ; // index Reco
	       IdxMatched.push_back( 1 ) ;
	       deltaR_min.push_back( DeltaReco1HLT1 ) ;
	       deltaR_min.push_back( DeltaReco2HLT2 ) ;
	     }else
	     {
	       IdxMatched.push_back( 1 ) ; // index Reco
	       IdxMatched.push_back( 0 ) ;
	       deltaR_min.push_back( DeltaReco2HLT1 ) ;
	       deltaR_min.push_back( DeltaReco1HLT2 ) ;
	     }

	   for(int ii=0; ii<deltaR_min.size(); ii++){ // loop on jets
	     if(verbose) cout << "deltaR_min " << deltaR_min.at(ii) << endl;
	     CreateAndFillUserTH1D("deltaR_minimum", 500, 0, 5, deltaR_min.at(ii) );
	     
	     if(deltaR_min.at(ii) < DeltaR_ ){ // dR check
	       CreateAndFillUserTH1D("deltaR_minimum_AfterMatch", 500, 0, 5, deltaR_min.at(ii) );
	       
	       if(verbose){
		 std::cout<<"WIDEJETS" <<std::endl;
		 std::cout<<"Reco"<<std::endl;
		 std::cout<<"Index= "<< IdxMatched.at(ii)  <<std::endl;
		 std::cout<<"Pt= "<< widejetsReco[IdxMatched.at(ii)].Pt()<<" Eta= "<<widejetsReco[IdxMatched.at(ii)].Eta()<< " Phi= "<<widejetsReco[IdxMatched.at(ii)].Phi()  <<std::endl;
		 std::cout<<"HLT" <<std::endl;
		 std::cout<<"Pt= "<< widejets[ii].Pt()<<" Eta= "<<widejets[ii].Eta()<< " Phi= "<<widejets[ii].Phi()  <<std::endl;
	       }
	       
	       // Reco variables
	       CreateAndFillUserTH1D("pT_WideJetReco", 150, 0, 3000, widejetsReco[IdxMatched.at(ii)].Pt()  );
	       CreateAndFillUserTH1D("Eta_WideJetReco", 20, -2.5, 2.5, widejetsReco[IdxMatched.at(ii)].Eta() );
	       CreateAndFillUserTH1D("Phi_WideJetReco", 100, -3.5, 3.5, widejetsReco[IdxMatched.at(ii)].Phi()  );
	       // HLT variables	       
	       CreateAndFillUserTH1D("pT_WideJetHLT", 150, 0, 3000, widejets[ii].Pt()  );
	       CreateAndFillUserTH1D("Eta_WideJetHLT", 20, -2.5, 2.5, widejets[ii].Eta()  );
	       CreateAndFillUserTH1D("Phi_WideJetHLT", 100, -3.5, 3.5, widejets[ii].Phi()  );
	       // Reco vs HLT
	       CreateAndFillUserTH2D("pT_WideJet", 150, 0, 3000, 150, 0, 3000, widejetsReco[IdxMatched.at(ii)].Pt(), widejets[ii].Pt() );
	       CreateAndFillUserTH2D("Eta_WideJet", 20, -2.5, 2.5, 20, -2.5, 2.5, widejetsReco[IdxMatched.at(ii)].Eta(), widejets[ii].Eta() );
	       CreateAndFillUserTH2D("Phi_WideJet", 100, -3.5, 3.5, 100, -3.5, 3.5, widejetsReco[IdxMatched.at(ii)].Phi(), widejets[ii].Phi() );

	       double pTBias = ( widejetsReco[IdxMatched.at(ii)].Pt() - widejets[ii].Pt() )/ widejets[ii].Pt() ; 
	       if(verbose) std::cout<<"pT Bias= "<<pTBias << std::endl;

	       int etaBin = mEtaBinning.getBin( fabs(widejets[ii].Eta()) ); // modulo
	       // int etaBin = mEtaBinning.getBin( widejets[ii].Eta() ); //separo negativi da positivi
	       int ptBin   = mPtBinning.getPtBin( widejets[ii].Pt() );
	       std::pair<float, float> ptBins = mPtBinning.getBinValue( ptBin );
	       std::string etaName = mEtaBinning.getBinName( etaBin );
	       std::string EtaBin = TString::Format("%s", etaName.c_str() ).Data();
	       std::string PtBin = TString::Format("pT_%i_%i", (int) ptBins.first, (int) ptBins.second ).Data();
	       std::string Bin = TString::Format("%s_pT_%i_%i", etaName.c_str(), (int) ptBins.first, (int) ptBins.second ).Data();

	       if( verbose ){	       
		 cout<<"WideJet: Eta = "<< widejets[ii].Eta() << " pT = "<< widejets[ii].Pt() <<endl;
		 cout<<"EtaBin =   "<< etaBin <<endl;
		 std::pair<float, float> etaBins = mEtaBinning.getBinValue( etaBin);
		 cout<<"EtaBin.first  "<< etaBins.first << "   EtaBin.second   "<< etaBins.second<<endl;
		 cout<<"pTBin =   " << ptBin  <<endl;
		 cout<<"ptBin.first  "  << ptBins.first  << "    ptBin.second   "<< ptBins.second<<endl;
		 cout<<"Eta Bin  "  << EtaBin.c_str() <<endl;
		 cout<<"Pt Bin  "  << PtBin.c_str() <<endl;
		 cout<<"Bin  "  << Bin.c_str() <<endl;
	       }

	       CreateAndFillUserTH1D("pTBias", 300, -2, 2, pTBias );
	       // devo guardare questi e... --> fit crystal ball per prendere il picco
	       CreateAndFillUserTH1D( Form("pTBias_%s", Bin.c_str()), 300, -2., 2., pTBias );

	       ////// TProfile: pTBias vs nVtx
	       CreateAndFillUserTProfile("pTBias_vs_Nvtx",  50, 0., 50., -2., 2.,      nVtx, pTBias );
	       CreateAndFillUserTProfile("pT_vs_NvtxReco", 50, 0., 50., 0., 5000., nVtxreco, widejetsReco[IdxMatched.at(ii)].Pt() );
	       CreateAndFillUserTProfile("pT_vs_NvtxHLT",  50, 0., 50., 0., 5000., nVtx, widejets[ii].Pt() );
	       CreateAndFillUserTProfile( Form("pTBias_vs_Nvtx_%s", Bin.c_str()),  50, 0., 50., -2., 2.,      nVtx, pTBias );
	       CreateAndFillUserTProfile( Form("pT_vs_NvtxReco_%s", Bin.c_str()), 50, 0., 50., 0., 5000., nVtxreco, widejetsReco[IdxMatched.at(ii)].Pt() );
	       CreateAndFillUserTProfile( Form("pT_vs_NvtxHLT_%s", Bin.c_str()),  50, 0., 50., 0., 5000., nVtx, widejets[ii].Pt() );
	       ////// TProfile: pTBias vs Eta
	       CreateAndFillUserTProfile("pTBias_vs_Eta",  20, -2.5, 2.5, -2., 2.,      widejets[ii].Eta(), pTBias );
	       CreateAndFillUserTProfile("pT_vs_EtaReco", 20, -2.5, 2.5, 0., 5000., widejetsReco[IdxMatched.at(ii)].Eta(), widejetsReco[IdxMatched.at(ii)].Pt() );
	       CreateAndFillUserTProfile("pT_vs_EtaHLT",  20, -2.5, 2.5, 0., 5000., widejets[ii].Eta(), widejets[ii].Pt() );
	       CreateAndFillUserTProfile( Form("pTBias_vs_Eta_%s", PtBin.c_str()),  20, -2.5, 2.5, -2., 2.,      widejets[ii].Eta(), pTBias );
	       CreateAndFillUserTProfile( Form("pT_vs_EtaReco_%s", PtBin.c_str()), 20, -2.5, 2.5, 0., 5000., widejetsReco[IdxMatched.at(ii)].Eta(), widejetsReco[IdxMatched.at(ii)].Pt() );
	       CreateAndFillUserTProfile( Form("pT_vs_EtaHLT_%s", PtBin.c_str()),  20, -2.5, 2.5, 0., 5000., widejets[ii].Eta(), widejets[ii].Pt() );
	       // ...e questo che e' quello importante -- adesso sto usando il bias 1D
	       CreateAndFillUserTProfile( Form("pTBias_vs_Eta_VariableBin_%s", PtBin.c_str()),  binnumEta, etaBins, -2., 2., widejets[ii].Eta(), pTBias );
	       ////// TProfile: pTBias vs Pt
	       CreateAndFillUserTProfile("pTBias_vs_Pt",  50, 0., 5000., -2., 2.,      widejets[ii].Pt(), pTBias );
	       CreateAndFillUserTProfile("pT_vs_PtReco", 50, 0., 5000., 0., 5000., widejetsReco[IdxMatched.at(ii)].Pt(), widejetsReco[IdxMatched.at(ii)].Pt() );
	       CreateAndFillUserTProfile("pT_vs_PtHLT",  50, 0., 5000., 0., 5000., widejets[ii].Pt(), widejets[ii].Pt() );	      
	       CreateAndFillUserTProfile( Form("pTBias_vs_Pt_%s", EtaBin.c_str()),  50, 0., 5000., -2., 2.,      widejets[ii].Pt(), pTBias );
	       CreateAndFillUserTProfile( Form("pT_vs_PtReco_%s", EtaBin.c_str()), 50, 0., 5000., 0., 5000., widejetsReco[IdxMatched.at(ii)].Pt(), widejetsReco[IdxMatched.at(ii)].Pt() );
	       CreateAndFillUserTProfile( Form("pT_vs_PtHLT_%s", EtaBin.c_str()),  50, 0., 5000., 0., 5000., widejets[ii].Pt(), widejets[ii].Pt() );	      
	       /*
	       ////// TProfile: pTBias vs mJJ
	       CreateAndFillUserTProfile("pTBias_vs_Mjj",  500, 0., 5000., -2., 2.,      mjj, pTBias );
	       CreateAndFillUserTProfile("pT_vs_MjjReco", 500, 0., 5000., 0., 5000., mjjreco, widejetsReco[IdxMatched.at(ii)].Pt() );
	       CreateAndFillUserTProfile("pT_vs_MjjHLT",  500, 0., 5000., 0., 5000., mjj, widejets[ii].Pt() );	      
	       CreateAndFillUserTProfile( Form("pTBias_vs_Mjj_%s", Bin.c_str()),  500, 0., 5000., -2., 2.,      mjj, pTBias );
	       CreateAndFillUserTProfile( Form("pT_vs_MjjReco_%s", Bin.c_str()), 500, 0., 5000., 0., 5000., mjjreco, widejetsReco[IdxMatched.at(ii)].Pt() );
	       CreateAndFillUserTProfile( Form("pT_vs_MjjHLT_%s", Bin.c_str()),  500, 0., 5000., 0., 5000., mjj, widejets[ii].Pt() );	      
	       */



	     }// if matched
	   }// loop on jets

	   if(verbose){
	     std::cout<<"DeltaR min 0= "<<deltaR_min.at(0) << std::endl;
	     std::cout<<"DeltaR min 1= "<<deltaR_min.at(1) << std::endl;
	     std::cout<<"deltaEtaCut_ = "<< DeltaEtaJJ_ << std::endl;
	     std::cout<<"MJJCut_ = "<< MJJCut_ << std::endl;
	   }	     

	   //chiedo entrambi i jet matcheati per studiare variabili comuni
	   if(deltaR_min.at(0) < DeltaR_ && deltaR_min.at(1) < DeltaR_){
	     if(verbose) std::cout<<"Both jets matched " << std::endl;
	     if( deltaETAjj < DeltaEtaJJ_){ // cut only on HLT
	       if(verbose) std::cout<<"deltaEta Cut Passed " << std::endl;
	       if( mjj > MJJCut_ ){ //mass HLT
		 if(verbose) std::cout<<"mJJCut Passed " << std::endl;
		 
		 /* // trigger selezionato all'inizio // altrimenti posso fare cosi:
		    if(verbose){
		    std::cout<<"mjj= " << mjj << std::endl;
		    if( passHLT_L1HTT150_BtagSeq || passHLT_L1HTT150 ){
		    cout<<"Trigger L1HTT passato"<< endl;
		    }
		    if( passHLT_HT450_BtagSeq ||passHLT_HT450 ){
		    cout<<"Trigger HT450 passato"<< endl;
		    }
		    }		 
		    if(mjj < 600 &&  !passHLT_L1HTT150_BtagSeq && !passHLT_L1HTT150 ) continue;
		    if(mjj > 600 &&  !passHLT_HT450_BtagSeq && !passHLT_HT450 ) continue;
		 */
	
		 //		 if(verbose) cout<<"Evento Passato"<< endl;
		 
		 CreateAndFillUserTH1D("metHLT", 50, 0, 5000, met );
		 CreateAndFillUserTH1D("metReco", 50, 0, 5000, metreco );
		 CreateAndFillUserTH2D("met", 50, 0, 5000, 50, 0, 5000, metreco, met );
		 
		 CreateAndFillUserTH1D("metSigHLT", 100, 0, 2, metSig );
		 CreateAndFillUserTH1D("metSigReco", 100, 0, 2, metrecoSig );
		 CreateAndFillUserTH2D("metSig", 100, 0, 2, 100, 0, 2, metrecoSig, metSig );
		 
		 CreateAndFillUserTH1D("NVertexHLT", 50, 0, 50, nVtx );
		 CreateAndFillUserTH1D("NVertexReco", 50, 0, 50, nVtxreco );
		 CreateAndFillUserTH2D("NVertex", 50, 0, 50, 50, 0, 50, nVtxreco, nVtx );
		 
		 CreateAndFillUserTH1D("deltaETAjjHLT", 100, 0, 3, deltaETAjj );
		 CreateAndFillUserTH1D("deltaETAjjReco", 100, 0, 3, deltaETAjjreco );
		 CreateAndFillUserTH2D("deltaETAjj", 100, 0, 3, 100, 0, 3, deltaETAjjreco, deltaETAjj );
		 
		 CreateAndFillUserTH1D("MJJ_WideJetHLT", 30, 0, 3000, mjj );
		 CreateAndFillUserTH1D("MJJ_WideJetReco", 30, 0, 3000, mjjreco );
		 CreateAndFillUserTH2D("MJJ_WideJet", 30, 0, 3000, 30, 0, 3000, mjjreco, mjj );
		 
		 // bias per la variabile mJJ
		 // mJJ(Reco) - mJJ(HLT) / mJJ(HLT)
		 double mJJBias =  (mjjreco - mjj) / mjj;
		 if(verbose) std::cout<<"mJJ Bias= "<<mJJBias << std::endl;
		 
		 int MjjBin   = mMjjBinning.getBin( mjj );
		 std::pair<float, float> mJJBins = mMjjBinning.getBinValue( MjjBin );
		 std::string mJJBin = TString::Format("mJJ_%i_%i", (int) mJJBins.first, (int) mJJBins.second ).Data();
		 
		 if( verbose ){	       
		   cout<<"WideJet: mJJ = "<< mjj << endl;
		   cout<<"mJJBin =   "<< MjjBin <<endl;
		   cout<<"mJJBin.first  "  << mJJBins.first  << "    mJJBin.second   "<< mJJBins.second<<endl;
		   cout<<"mJJBin  "  << mJJBin.c_str() <<endl;
		 }
		 
		 CreateAndFillUserTH1D("mJJBias", 300, -2, 2, mJJBias);  
		 // devo guardare questi e..
		 CreateAndFillUserTH1D( Form("mJJBias_%s", mJJBin.c_str()), 300, -2., 2., mJJBias );
		 
		 ////// TProfile: mJJBias vs Nvtx
		 CreateAndFillUserTProfile("mJJBias_vs_Nvtx",  50, 0., 50., -2., 2.,      nVtx, mJJBias );
		 CreateAndFillUserTProfile("mJJ_vs_NvtxReco", 50, 0., 50., 0., 5000., nVtxreco, mjjreco );
		 CreateAndFillUserTProfile("mJJ_vs_NvtxHLT",  50, 0., 50., 0., 5000., nVtx, mjj );
		 CreateAndFillUserTProfile( Form("mJJBias_vs_Nvtx_%s", mJJBin.c_str()),  50, 0., 50., -2., 2.,      nVtx, mJJBias );
		 CreateAndFillUserTProfile( Form("mJJ_vs_NvtxReco_%s", mJJBin.c_str()), 50, 0., 50., 0., 5000., nVtxreco, mjjreco );
		 CreateAndFillUserTProfile( Form("mJJ_vs_NvtxHLT_%s", mJJBin.c_str()),  50, 0., 50., 0., 5000., nVtx, mjj );

		 // for sui jet
		 for(int ii=0; ii<deltaR_min.size(); ii++){
		   // 2 entries per massa --> considero entrambi i jet
		   ////// TProfile: mJJBias vs Eta
		   CreateAndFillUserTProfile("mJJBias_vs_Eta",  20, -2.5, 2.5, -2., 2.,      widejets[ii].Eta(), mJJBias );
		   CreateAndFillUserTProfile("mJJ_vs_EtaReco",  20, -2.5, 2.5, 0., 5000., widejetsReco[IdxMatched.at(ii)].Eta(), mjjreco );
		   CreateAndFillUserTProfile("mJJ_vs_EtaHLT",   20, -2.5, 2.5, 0., 5000., widejets[ii].Eta(), mjj );
		   CreateAndFillUserTProfile( Form("mJJBias_vs_Eta_%s", mJJBin.c_str()),  20, -2.5, 2.5, -2., 2.,      widejets[ii].Eta(), mJJBias );
		   CreateAndFillUserTProfile( Form("mJJ_vs_EtaReco_%s", mJJBin.c_str()), 20, -2.5, 2.5, 0., 5000., widejetsReco[IdxMatched.at(ii)].Eta(), mjjreco );
		   CreateAndFillUserTProfile( Form("mJJ_vs_EtaHLT_%s", mJJBin.c_str()),  20, -2.5, 2.5, 0., 5000., widejets[ii].Eta(), mjj );
		   ////// TProfile: mJJBias vs Pt
		   CreateAndFillUserTProfile("mJJBias_vs_Pt",  30, 0., 3000., -2., 2.,      widejets[ii].Pt(), mJJBias );
		   CreateAndFillUserTProfile("mJJ_vs_PtReco", 30, 0., 3000., 0., 3000., widejetsReco[IdxMatched.at(ii)].Pt(), mjjreco );
		   CreateAndFillUserTProfile("mJJ_vs_PtHLT",  30, 0., 3000., 0., 3000., widejets[ii].Pt(), mjj );	      
		   CreateAndFillUserTProfile( Form("mJJBias_vs_Pt_%s", mJJBin.c_str()),  30, 0., 3000., -2., 2.,      widejets[ii].Pt(), mJJBias );
		   CreateAndFillUserTProfile( Form("mJJ_vs_PtReco_%s", mJJBin.c_str()), 30, 0., 3000., 0., 3000., widejetsReco[IdxMatched.at(ii)].Pt(), mjjreco );
		   CreateAndFillUserTProfile( Form("mJJ_vs_PtHLT_%s", mJJBin.c_str()),  30, 0., 3000., 0., 3000., widejets[ii].Pt(), mjj );	      	       
		 }

		 ////// TProfile: mJJBias vs mJJ
		 // questo e' quello importante!!! basta fittare questo e ho le correzioni
		 CreateAndFillUserTProfile("mJJBias_vs_Mjj",  25, 500., 3000., -2., 2.,  mjj, mJJBias );
		 CreateAndFillUserTProfile("mJJ_vs_MjjReco", 50, 0., 5000., 0., 5000., mjjreco, mjjreco );
		 CreateAndFillUserTProfile("mJJ_vs_MjjHLT",  50, 0., 5000., 0., 5000., mjj, mjj );	      
		 CreateAndFillUserTProfile( Form("mJJBias_vs_Mjj_%s", mJJBin.c_str()),  30, 0., 3000., -2., 2.,  mjj, mJJBias );
		 CreateAndFillUserTProfile( Form("mJJ_vs_MjjReco_%s", mJJBin.c_str()), 30, 0., 3000., 0., 3000., mjjreco, mjjreco );
		 CreateAndFillUserTProfile( Form("mJJ_vs_MjjHLT_%s", mJJBin.c_str()),  30, 0., 3000., 0., 3000., mjj, mjj );	      
		 
	       }// mjj > mCUT	  	  
	     }// deltaEtaJJ < cut	
	   }//both matched
	 }// widejet comparison
	  	   
	 //       } // triggers
     } // passJson
       
       /////////////////////////////////////////////////////////////////
     
       //== Fill Variables ==
     fillVariableWithValue("PassJSON", PassJSON);     
     fillVariableWithValue("MJJ_reco", mjjreco);     
     fillVariableWithValue("MJJ", mjj);     
     fillVariableWithValue("deltaETAjj_reco", deltaETAjjreco);     
     fillVariableWithValue("deltaETAjj", deltaETAjj);     
     //	       fillVariableWithValue("WideJet_pTBias", pTBias);     
     //	       fillVariableWithValue("AK4_pTBias", pTBias);     

     ////////////////////////////////////////

     //no cuts on these variables, just to store in output
     // if(!isData)
     //   fillVariableWithValue("trueVtx",PileupInteractions->at(idx_InTimeBX));
     // else if(isData)
     //   fillVariableWithValue("trueVtx",999);     

     // Trigger
     //int NtriggerBits = triggerResult->size();
     /*   
	  if (isData)
	  {
	  fillVariableWithValue("passHLT_ZeroBias_BtagSeq",triggerResult->at(8));// DST_ZeroBias_BTagScouting_v* (run>=259636)
	  fillVariableWithValue("passHLT_ZeroBias",triggerResult->at(7));// DST_ZeroBias_PFScouting_v* (run>=259636)

	  fillVariableWithValue("passHLT_L1DoubleMu_BtagSeq",triggerResult->at(9));// DST_L1DoubleMu_BTagScouting_v* (run>=259636)
	  fillVariableWithValue("passHLT_L1DoubleMu",triggerResult->at(10));// DST_L1DoubleMu_PFScouting_v* (run>=259636)

	  fillVariableWithValue("passHLT_CaloJet40_BtagSeq",triggerResult->at(0));//  DST_CaloJet40_PFReco_PFBTagCSVReco_PFScouting_v* (257933<=run<259636) 
	  //  OR DST_CaloJet40_BTagScouting_v* (run>=259636)
	  fillVariableWithValue("passHLT_CaloJet40",triggerResult->at(1));// DST_CaloJet40_CaloScouting_PFScouting_v*  (run>=259636)

	  fillVariableWithValue("passHLT_L1HTT150_BtagSeq",triggerResult->at(2));// DST_L1HTT125ORHTT150ORHTT175_PFReco_PFBTagCSVReco_PFScouting_v* (257933<=run<259636) 
	  // OR DST_L1HTT_BTagScouting_v* (run>=259636)    
	  fillVariableWithValue("passHLT_L1HTT150",triggerResult->at(3));// DST_L1HTT_CaloScouting_PFScouting_v* (run>=259636)

	  fillVariableWithValue("passHLT_HT450_BtagSeq",triggerResult->at(5));// DST_HT450_PFReco_PFBTagCSVReco_PFScouting_v* (257933<=run<259636) 
	  // OR DST_HT450_BTagScouting_v* (run>=259636)    
	  fillVariableWithValue("passHLT_HT450",triggerResult->at(6));// DST_HT450_PFScouting_v* (run>=259636)        

	  fillVariableWithValue("passHLT_PFHT800",triggerResult->at(13));// HLT_PFHT800_v* (all runs)   
	  fillVariableWithValue("passHLT_PFHT650MJJ950",triggerResult->at(22));// HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v* (all runs)        
	  fillVariableWithValue("passHLT_PFHT650MJJ900",triggerResult->at(23));// HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v* (all runs)        
	  }
     */
     // Evaluate cuts (but do not apply them)
     evaluateCuts();
     
     // optional call to fill a skim with the full content of the input roottuple
     //if( passedCut("nJetFinal") ) fillSkimTree();
     /*
       if( passedCut("PassJSON")
       && passedCut("nVtx") 
       && passedCut("IdTight_j1")
       && passedCut("IdTight_j2")
       && passedCut("nJet")
       && passedCut("pTWJ_j1")
       && passedCut("etaWJ_j1")
       && passedCut("pTWJ_j2")
       && passedCut("etaWJ_j2")
       && getVariableValue("deltaETAjj") <  getPreCutValue1("DetaJJforTrig") ){

       h_mjj_NoTrigger -> Fill(MJJWideReco); 
       
       if( (getVariableValue("passHLT_ZeroBias_BtagSeq")||getVariableValue("passHLT_ZeroBias")) )
       h_mjj_HLTpass_ZeroBias -> Fill(MJJWideReco);  
       if( (getVariableValue("passHLT_ZeroBias_BtagSeq")||getVariableValue("passHLT_ZeroBias")) 
       && (getVariableValue("passHLT_L1HTT150_BtagSeq")||getVariableValue("passHLT_L1HTT150")) )
       h_mjj_HLTpass_ZeroBias_L1HTT150 -> Fill(MJJWideReco);  

       if( (getVariableValue("passHLT_L1HTT150_BtagSeq")||getVariableValue("passHLT_L1HTT150")) )
       h_mjj_HLTpass_L1HTT150 -> Fill(MJJWideReco);  

       if( (getVariableValue("passHLT_L1HTT150_BtagSeq")||getVariableValue("passHLT_L1HTT150")) 
       && (getVariableValue("passHLT_HT450_BtagSeq")||getVariableValue("passHLT_HT450")) )
       h_mjj_HLTpass_L1HTT150_HT450 -> Fill(MJJWideReco);  
       
       }
     */
     // optional call to fill a skim with a subset of the variables defined in the cutFile (use flag SAVE)

     //     if( passedAllPreviousCuts("mjj") && passedCut("mjj") ) 
     //       {
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

     //       }

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

   h_pTJet_Binned -> Write();

   //  h_mjj_NoTrigger -> Write();
   //  h_mjj_HLTpass_ZeroBias -> Write();
   //  h_mjj_HLTpass_ZeroBias_L1HTT150 -> Write();
   //  h_mjj_HLTpass_L1HTT150 -> Write();
   //  h_mjj_HLTpass_L1HTT150_HT450 -> Write();

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

//Federico

/*
// analysis in bin of number of vertices
VertexBinning mVtxBinning;
	     
int vtxBinReco = mVtxBinning.getVertexBin( nVtxreco );
int vtxBinHLT = mVtxBinning.getVertexBin( nVtx );
	     
std::pair<int, int> vtxBinsReco = mVtxBinning.getBinValue(vtxBinReco);
std::pair<int, int> vtxBinsHLT = mVtxBinning.getBinValue(vtxBinHLT);
	     
if(verbose){
cout<<"vtxReco  "<< nVtxreco <<endl;
cout<<"vtx  "<< nVtx <<endl;
cout<<"vtxBin(Reco)  "<< vtxBinReco<<endl;
cout<<"vtxBin.first  "<< vtxBinsReco.first << "    vtxBin.second   "<<vtxBinsReco.second<<endl;
cout<<"vtxBin(HLT)  "<< vtxBinHLT <<endl;
cout<<"vtxBin.first  "<< vtxBinsHLT.first << "    vtxBin.second   "<<vtxBinsHLT.second<<endl;
}
	     
std::string HistoNameReco = TString::Format("pTWJ_recoj1_nvtx_%i_%i",  vtxBinsReco.first, vtxBinsReco.second ).Data();
std::string HistoNameHLT  = TString::Format("pTWJ_j1_nvtx_%i_%i",  vtxBinsHLT.first, vtxBinsHLT.second ).Data();
	     
std::cout << HistoNameReco.c_str()<< std::endl;
std::cout << HistoNameHLT.c_str()<< std::endl;
	     
CreateAndFillUserTH1D(HistoNameReco.c_str(), 300, -2, 2, pTWJ_recoj1 );
CreateAndFillUserTH1D(HistoNameHLT.c_str(), 300, -2, 2, pTWJ_j1 );
*/



/*	     
//////////////////// analysis in bin di pT e Eta	     
EtaBinning mEtaBinning;
PtBinning mPtBinning;
	     
int etaBin = mEtaBinning.getBin( fabs(widejetsReco[IdxMatched.at(ii)].Eta() ) );
int ptBin = mPtBinning.getPtBin( widejetsReco[IdxMatched.at(ii)].Pt() );
	       
std::pair<float, float> ptBins = mPtBinning.getBinValue(ptBin);	     
	       
if(verbose){
std::pair<float, float> etaBins = mEtaBinning.getBinValue(etaBin);
cout<<"etaReco  "<< fabs(widejetsReco[IdxMatched.at(ii)].Eta() )  << "    pT Reco   "<<  widejetsReco[IdxMatched.at(ii)].Pt() <<endl;
cout<<"etaBin  "<< etaBin << "    pTBin   "<<ptBin<<endl;
cout<<"etaBin.first  "<< etaBins.first << "    etaBin.second   "<<etaBins.second<<endl;
cout<<"ptBin.first  "<< ptBins.first << "    ptBin.second   "<<ptBins.second<<endl;
}
	       
std::string etaName = mEtaBinning.getBinName(etaBin); 
std::string HistoName = TString::Format("pTBias_%s_pT_%i_%i", etaName.c_str(), (int) ptBins.first, (int) ptBins.second ).Data();
std::string HistoName2 = TString::Format("pt_WideJet_%s_pT_%i_%i", etaName.c_str(), (int) ptBins.first, (int) ptBins.second ).Data();
	       
if(verbose) std::cout << HistoName.c_str()<< std::endl;
	       
CreateAndFillUserTH1D(HistoName.c_str(), 300, -2, 2, pTBias.at(ii));
CreateAndFillUserTH2D(HistoName2.c_str(), 300, 0, 3000, 300, 0, 3000, widejetsReco[IdxMatched.at(ii)].Pt(), widejets[ii].Pt() );
*/
//	     }// loop on jets


	   /*
	   // only for AK4: jets Asymmetry --> Resolution
	   // Asymmetry only for events with 2 matched jets 
	   if(deltaR_min.at(0) < DeltaR_ && deltaR_min.at(1) < DeltaR_ ){   

	     int random;
	     random = lrand48()%(2);
	     if(verbose) std::cout<< "random = " <<random << std::endl;
	     // Reco
	     double pTaveReco = ( ak4j1Reco.Pt() + ak4j2Reco.Pt() ) / 2;
	     if(verbose) std::cout<< "Average pT Reco = " << pTaveReco << std::endl;
	     CreateAndFillUserTH1D("pTaveReco", 30 , 0, 3000, pTaveReco);
	     
	     double alphaReco = ak4j3Reco.Pt() / pTaveReco ;
	     if(verbose) std::cout<< "alpha Reco = " << alphaReco << std::endl;
	     CreateAndFillUserTH1D("alphaReco", 100 , 0, 2, alphaReco);
	     double AsymmetryReco = 999;
	     
	     int etaBin_j1Reco = mEtaBinning.getBin( fabs(ak4j1Reco.Eta()) );
	     int etaBin_j2Reco = mEtaBinning.getBin( fabs(ak4j2Reco.Eta()) );
	     int ptBin_Reco = mPtBinning.getPtBin( pTaveReco );
	     std::pair<float, float> ptBins_Reco = mPtBinning.getBinValue(ptBin_Reco);

	     if(verbose){       
	       cout<<"|Eta_ak4j1Reco| =   "<<fabs(ak4j1Reco.Eta())<<endl;
	       cout<<"|Eta_ak4j2Reco| =   "<<fabs(ak4j2Reco.Eta())<<endl;
	       std::pair<float, float> etaBins_j1 = mEtaBinning.getBinValue(etaBin_j1Reco);
	       std::pair<float, float> etaBins_j2 = mEtaBinning.getBinValue(etaBin_j2Reco);
	       cout<<"EtaBin_j1 =   "<<etaBin_j1Reco<<endl;
	       cout<<"EtaBin_j2 =   "<<etaBin_j2Reco<<endl;
	       cout<<"pTBin =   "<<ptBin_Reco<<endl;
	       cout<<"j1: etaBin.first  "<< etaBins_j1.first << "    etaBin.second   "<<etaBins_j1.second<<endl;
	       cout<<"j2: etaBin.first  "<< etaBins_j2.first << "    etaBin.second   "<<etaBins_j2.second<<endl;
	       cout<<"ptBin.first  "<< ptBins_Reco.first << "    ptBin.second   "<<ptBins_Reco.second<<endl;
	     }

	     /// HLT
	     double pTave = ( ak4j1.Pt() + ak4j2.Pt() ) /2;
	     if(verbose) std::cout<< "Average pT = " << pTave << std::endl;
	     CreateAndFillUserTH1D("pTave", 30 , 0, 3000, pTave);
	       
	     double alpha = ak4j3.Pt() / pTave ;
	     if(verbose) std::cout<< "alpha = " << alpha << std::endl;
	     CreateAndFillUserTH1D("alpha", 100 , 0, 2, alpha);
	     double Asymmetry = 999;
	       
	     int etaBin_j1 = mEtaBinning.getBin( fabs(ak4j1.Eta()) );
	     int etaBin_j2 = mEtaBinning.getBin( fabs(ak4j2.Eta()) );
	     int ptBin = mPtBinning.getPtBin( pTave );
	     std::pair<float, float> ptBins = mPtBinning.getBinValue(ptBin);	       
	   
	     if(verbose){
	       cout<<"|Eta_ak4j1| =   "<<fabs(ak4j1.Eta())<<endl;
	       cout<<"|Eta_ak4j2| =   "<<fabs(ak4j2.Eta())<<endl;
	       std::pair<float, float> etaBins_j1 = mEtaBinning.getBinValue(etaBin_j1);
	       std::pair<float, float> etaBins_j2 = mEtaBinning.getBinValue(etaBin_j2);
	       cout<<"EtaBin_j1 =   "<<etaBin_j1<<endl;
	       cout<<"EtaBin_j2 =   "<<etaBin_j2<<endl;
	       cout<<"pTBin =   "<<ptBin <<endl;
	       cout<<"j1: etaBin.first  "<< etaBins_j1.first << "    etaBin.second   "<<etaBins_j1.second<<endl;
	       cout<<"j2: etaBin.first  "<< etaBins_j2.first << "    etaBin.second   "<<etaBins_j2.second<<endl;
	       cout<<"ptBin.first  "<< ptBins.first << "    ptBin.second   "<<ptBins.second<<endl;
	     }

	     // seleziono gli stessi eventi	     
	     if(etaBin_j1Reco == etaBin_j2Reco){ // same eta bin RECO
	       std::string etaNameReco = mEtaBinning.getBinName(etaBin_j1Reco);	     
	       if(etaBin_j1 == etaBin_j2){ // same eta bin HLT
		 std::string etaName = mEtaBinning.getBinName(etaBin_j1);	     
		 if(alphaReco < alphacut_ ){	// alpha cut RECO
		   if(alpha < alphacut_ ){	// aloha cut HLT
		    
		     if(random == 0){
		       AsymmetryReco = ( ak4j1Reco.Pt() - ak4j2Reco.Pt() ) / ( ak4j1Reco.Pt() + ak4j2Reco.Pt() ) ;
		       Asymmetry        = ( ak4j1.Pt() - ak4j2.Pt() ) / ( ak4j1.Pt() + ak4j2.Pt() ) ;
		     }else if(random == 1){
		       AsymmetryReco = ( ak4j2Reco.Pt() - ak4j1Reco.Pt() ) / ( ak4j1Reco.Pt() + ak4j2Reco.Pt() ) ; 
		       Asymmetry        = ( ak4j2.Pt() - ak4j1.Pt() ) / ( ak4j1.Pt() + ak4j2.Pt() ) ; 
		     }
		     std::string HistoNameReco = TString::Format("AsymmetryReco_%s_pT_%i_%i", etaNameReco.c_str(), (int) ptBins_Reco.first, (int) ptBins_Reco.second ).Data();	 	     
		     if(verbose) std::cout << HistoNameReco.c_str()<< std::endl;
		     if(verbose) std::cout<< "Asymmetry Reco = " << AsymmetryReco << std::endl; 
		     CreateAndFillUserTH1D(HistoNameReco.c_str(), 300, -2, 2, AsymmetryReco);
		     
		     std::string HistoName = TString::Format("Asymmetry_%s_pT_%i_%i", etaName.c_str(), (int) ptBins.first, (int) ptBins.second ).Data();	 	     
		     if(verbose) std::cout << HistoName.c_str()<< std::endl;
		     if(verbose) std::cout<< "Asymmetry  = " << Asymmetry << std::endl; 
		     CreateAndFillUserTH1D(HistoName.c_str(), 300, -2, 2, Asymmetry);
		     
		   }else{//alpha cut HLT
		     if(verbose) std::cout<< "HLT: alpha too large = "<<alpha<<std::endl;
		   }
		 }else{
		   if(verbose) std::cout<< "Reco: alpha too large = "<<alphaReco<<std::endl;
		 }
	       }// eta bin HLT
	     } // eta bin Reco
	   }// both matched jets 	     
*/
