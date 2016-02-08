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
#include "vertexBinning.h"

bool verbose = true;
int Counter = 0;   // debugging -- can be erased    
double DeltaR_;
double MJJCut_;

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile)
{
  std::cout << "analysisClass::analysisClass(): begins " << std::endl;

  //  std::string jetAlgo = getPreCutString1("jetAlgo");
  //  double rParam = getPreCutValue1("DeltaR");

  DeltaR_ = getPreCutValue1("DeltaR");
  MJJCut_ = getPreCutValue1("MJJCut");

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

   TProfile *pTResolution_vs_Nvtx = new TProfile("pTResolution_vs_Nvtx","pT Resolution vs Nvtx", 50, 0., 50., -2., 2.);
   TProfile *pT_vs_NvtxReco = new TProfile("pT_vs_NvtxReco","pT vs Nvtx", 50, 0., 50., 0., 900.);
   TProfile *pT_vs_NvtxHLT = new TProfile("pT_vs_NvtxHLT","pT vs Nvtx", 50, 0., 50., 0., 900.);

   TProfile *pTResolution_vs_Nvtx_Barrel = new TProfile("pTResolution_vs_Nvtx_Barrel","pT Resolution vs Nvtx", 50, 0., 50., -2., 2.);
   TProfile *pT_vs_Nvtx_BarrelReco = new TProfile("pT_vs_Nvtx_BarrelReco","pT vs Nvtx", 50, 0., 50., 0., 900.);
   TProfile *pT_vs_Nvtx_BarrelHLT = new TProfile("pT_vs_Nvtx_BarrelHLT","pT vs Nvtx", 50, 0., 50., 0., 900.);

   TProfile *pTResolution_vs_Nvtx_Endcap = new TProfile("pTResolution_vs_Nvtx_Endcap","pT Resolution vs Nvtx", 50, 0., 50., -2., 2.);
   TProfile *pT_vs_Nvtx_EndcapReco = new TProfile("pT_vs_Nvtx_EndcapReco","pT vs Nvtx", 50, 0., 50., 0., 900.);
   TProfile *pT_vs_Nvtx_EndcapHLT = new TProfile("pT_vs_Nvtx_EndcapHLT","pT vs Nvtx", 50, 0., 50., 0., 900.);

   TProfile *AK4_pTResolution_vs_Nvtx = new TProfile("AK4_pTResolution_vs_Nvtx","pT Resolution vs Nvtx", 50, 0., 50., -2., 2.);
   TProfile *AK4_pT_vs_NvtxReco = new TProfile("AK4_pT_vs_NvtxReco","pT vs Nvtx", 50, 0., 50., 0., 900.);
   TProfile *AK4_pT_vs_NvtxHLT = new TProfile("AK4_pT_vs_NvtxHLT","pT vs Nvtx", 50, 0., 50., 0., 900.);


   // TH1F *h_nJetFinal = new TH1F ("h_nJetFinal","",10,0,10);
   // h_nJetFinal->Sumw2();      

   // variable binning for mjj trigger efficiency plots
   //   const int nMassBins = 103;

   //   double massBoundaries[nMassBins+1] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325,
   //  354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687,
   //  1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509,
   //  4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 
   //  10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000};

   // char* HLTname[50] = {"noTrig","PFHT475","PFHT800","PFHT650MJJ900","PFHT800_OR_PFHT650MJJ900","PFHT800_noPFHT475", 
   //                      "Mu45Eta2p1", "PFHT800AndMu45Eta2p1"};
   // TH1F* h_mjj_HLTpass[8];
   // char name_histoHLT[50];
   // for (int i=0; i<8; i++){  
   //   sprintf(name_histoHLT,"h_mjj_HLTpass_%s",HLTname[i]);
   //   h_mjj_HLTpass[i]= new TH1F(name_histoHLT,"",103,massBoundaries);
   // }

   /////////initialize variables

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) { // running over all events
   //   for (Long64_t jentry=0; jentry<10;jentry++) { //runnig over 2 events
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     // if (Cut(ientry) < 0) continue;
   
     ////////////////////// User's code starts here ///////////////////////
     
     //     verbose = true;

     if(verbose){
     std::cout<< std::endl;
     std::cout<<"Start"<< std::endl;
     std::cout<<"PassJSON "<< PassJSON<<std::endl;
     std::cout<<"RECO"<< std::endl;
     std::cout<<"AK4 jet1"<< std::endl;
     std::cout<<"Pt: "<< pTAK4_recoj1<<" eta: "<<etaAK4_recoj1<<" phi: "<<phiAK4_recoj1<< std::endl;
     std::cout<<"AK4 jet2"<< std::endl;
     std::cout<<"Pt: "<< pTAK4_recoj2<<" eta: "<<etaAK4_recoj2<<" phi: "<<phiAK4_recoj2<< std::endl;
     std::cout<<"Widejet 1"<< std::endl;
     std::cout<<"Pt: "<<pTWJ_recoj1<<" eta: "<<etaWJ_recoj1<<" phi: "<<phiWJ_recoj1<< std::endl;
     std::cout<<"Widejet 2"<< std::endl;
     std::cout<<"Pt: "<<pTWJ_recoj2<<" eta: "<<etaWJ_recoj2<<" phi: "<<phiWJ_recoj2<< std::endl;
     std::cout<<"HLT"<< std::endl;
     std::cout<<"AK4 jet1"<< std::endl;
     std::cout<<"Pt: "<< pTAK4_j1<<" eta: "<<etaAK4_j1<<" phi: "<<phiAK4_j1<< std::endl;
     std::cout<<"AK4 jet2"<< std::endl;
     std::cout<<"Pt: "<< pTAK4_j2<<" eta: "<<etaAK4_j2<<" phi: "<<phiAK4_j2<< std::endl;
     std::cout<<"Widejet 1"<< std::endl;
     std::cout<<"Pt: "<<pTWJ_j1<<" eta: "<<etaWJ_j1<<" phi: "<<phiWJ_j1<< std::endl;
     std::cout<<"Widejet 2"<< std::endl;
     std::cout<<"Pt: "<<pTWJ_j2<<" eta: "<<etaWJ_j2<<" phi: "<<phiWJ_j2<< std::endl;
     std::cout<<std::endl;
     }

     //     if( !PassJSON) std::cout << "TEST  CONFIGURATION: Pass JSON FALSE"<<std::endl;

     if( PassJSON){
       if( mjjreco > MJJCut_ && mjj > MJJCut_){
	 
	 // AK4 Jets
	 vector<TLorentzVector> AK4recojets;
	 TLorentzVector ak4j1Reco, ak4j2Reco, ak4dijetReco;             
	 // JetID
	 vector<double> chargedHadEnFracReco ;
	 vector<double> neutrHadEnFracReco ;
	 vector<double> chargedElectromFracReco ;
	 vector<double> neutrElectromFracReco ;
	 vector<double> muEnFracReco ;
	 vector<double> chargedMultReco ;
	 vector<double> neutrMultReco ;
	 vector<double> photonMultReco ;
	 // btag
	 vector<double> jetCSVAK4Reco ;
	 
	 vector<TLorentzVector> AK4jets;
	 TLorentzVector ak4j1, ak4j2, ak4dijet;      
	 // JetID
	 vector<double> chargedHadEnFrac ;
	 vector<double> neutrHadEnFrac ;
	 vector<double> chargedElectromFrac ;
	 vector<double> neutrElectromFrac ;
	 vector<double> muEnFrac ;
	 vector<double> chargedMult ;
	 vector<double> neutrMult ;
	 vector<double> photonMult ;
	 //btag
	 vector<double> jetCSVAK4 ;
	 
	 // Reco
	 ak4j1Reco.SetPtEtaPhiM(pTAK4_recoj1, etaAK4_recoj1, phiAK4_recoj1, massAK4_recoj1);
	 ak4j2Reco.SetPtEtaPhiM(pTAK4_recoj2, etaAK4_recoj2, phiAK4_recoj2, massAK4_recoj2);
	 
	 if( ak4j1Reco.Pt()>0 && ak4j2Reco.Pt()>0 )
	   {
	     // Put widejets in the container
	     AK4recojets.push_back( ak4j1Reco );
	     chargedHadEnFracReco.push_back(chargedHadEnFrac_recoj1) ;
	     neutrHadEnFracReco.push_back(neutrHadEnFrac_recoj1) ;
	     chargedElectromFracReco.push_back(chargedElectromFrac_recoj1) ;
	     neutrElectromFracReco.push_back(neutrElectromFrac_recoj1) ;
	     muEnFracReco.push_back(muEnFract_recoj1) ;
	     chargedMultReco.push_back(chargedMult_recoj1) ;
	     neutrMultReco.push_back(neutrMult_recoj1) ;
	     photonMultReco.push_back(photonMult_recoj1) ;
	     jetCSVAK4Reco.push_back(jetCSVAK4_recoj1) ;
	     
	     AK4recojets.push_back( ak4j2Reco );	 
	     chargedHadEnFracReco.push_back(chargedHadEnFrac_recoj2) ;
	     neutrHadEnFracReco.push_back(neutrHadEnFrac_recoj2) ;
	     chargedElectromFracReco.push_back(chargedElectromFrac_recoj2) ;
	     neutrElectromFracReco.push_back(neutrElectromFrac_recoj2) ;
	     muEnFracReco.push_back(muEnFract_recoj2) ;
	     chargedMultReco.push_back(chargedMult_recoj2) ;
	     neutrMultReco.push_back(neutrMult_recoj2) ;
	     photonMultReco.push_back(photonMult_recoj2) ;
	     jetCSVAK4Reco.push_back(jetCSVAK4_recoj2) ;
	   }
	 
	 // HLT
	 ak4j1.SetPtEtaPhiM(pTAK4_j1, etaAK4_j1, phiAK4_j1, massAK4_j1);
	 ak4j2.SetPtEtaPhiM(pTAK4_j2, etaAK4_j2, phiAK4_j2, massAK4_j2);
	 
	 if( ak4j1.Pt()>0 && ak4j2.Pt()>0 )
	   {
	     // Put widejets in the container
	     AK4jets.push_back( ak4j1 );
	     chargedHadEnFrac.push_back(chargedHadEnFrac_j1) ;
	     neutrHadEnFrac.push_back(neutrHadEnFrac_j1) ;
	     chargedElectromFrac.push_back(chargedElectromFrac_j1) ;
	     neutrElectromFrac.push_back(neutrElectromFrac_j1) ;
	     muEnFrac.push_back(muEnFract_j1) ;
	     chargedMult.push_back(chargedMult_j1) ;
	     neutrMult.push_back(neutrMult_j1) ;
	     photonMult.push_back(photonMult_j1) ;
	     jetCSVAK4.push_back(jetCSVAK4_j1) ;
	     
	     AK4jets.push_back( ak4j2 );
	     chargedHadEnFrac.push_back(chargedHadEnFrac_j2) ;
	     neutrHadEnFrac.push_back(neutrHadEnFrac_j2) ;
	     chargedElectromFrac.push_back(chargedElectromFrac_j2) ;
	     neutrElectromFrac.push_back(neutrElectromFrac_j2) ;
	     muEnFrac.push_back(muEnFract_j2) ;
	     chargedMult.push_back(chargedMult_j2) ;
	     neutrMult.push_back(neutrMult_j2) ;
	     photonMult.push_back(photonMult_j2) ;
	     jetCSVAK4.push_back(jetCSVAK4_j2) ;
	   }
	 
	 // Compare with AK4 Jets
	 if( AK4recojets.size() !=0  && AK4jets.size() !=0 ) {
	   
	   //////////// matching	 
	   double DeltaReco1HLT1 = AK4recojets[0].DeltaR(AK4jets[0]);
	   double DeltaReco1HLT2 = AK4recojets[0].DeltaR(AK4jets[1]);
	   double DeltaReco2HLT1 = AK4recojets[1].DeltaR(AK4jets[0]);
	   double DeltaReco2HLT2 = AK4recojets[1].DeltaR(AK4jets[1]);
	   
	   if(verbose){
	     //	   cout << endl;
	     cout << "AK4: " << endl;
	     cout << "deltaR Reco1-HLT1 " <<  DeltaReco1HLT1 << endl;
	     cout << "deltaR Reco1-HLT2 " <<  DeltaReco1HLT2 << endl;
	     cout << "deltaR Reco2-HLT1 " <<  DeltaReco2HLT1 << endl;
	     cout << "deltaR Reco2-HLT2 " <<  DeltaReco2HLT2 << endl;
	   }       
	   
	   vector<int> IdxMatching;
	   if( DeltaReco1HLT1 <= DeltaReco2HLT1 )
	     {
	       IdxMatching.push_back( 0 );
	       IdxMatching.push_back( 1 );
	     }else
	     {
	       IdxMatching.push_back( 1 );
	       IdxMatching.push_back( 0 );
	     }
	   
	   //////////// event variables	 

	   CreateAndFillUserTH1D("metHLT", 500, 0, 5000, met );
	   CreateAndFillUserTH1D("metReco", 500, 0, 5000, metreco );
	   CreateAndFillUserTH2D("met", 500, 0, 5000, 500, 0, 5000, metreco, met );
	   
	   CreateAndFillUserTH1D("metSigHLT", 100, 0, 2, metSig );
	   CreateAndFillUserTH1D("metSigReco", 100, 0, 2, metrecoSig );
	   CreateAndFillUserTH2D("metSig", 100, 0, 2, 100, 0, 2, metrecoSig, metSig );

	   CreateAndFillUserTH1D("NVertexHLT", 50, 0, 50, nVtx );
	   CreateAndFillUserTH1D("NVertexReco", 50, 0, 50, nVtxreco );
	   CreateAndFillUserTH2D("NVertex", 50, 0, 50, 50, 0, 50, nVtxreco, nVtx );
	   
	   CreateAndFillUserTH1D("AK4_deltaETAjjHLT", 100, 0, 3, deltaETAjjAK4 );
	   CreateAndFillUserTH1D("AK4_deltaETAjjReco", 100, 0, 3, deltaETAjjAK4reco );
	   CreateAndFillUserTH2D("AK4_deltaETAjj", 100, 0, 3, 100, 0, 3, deltaETAjjAK4reco, deltaETAjjAK4 );
	   
	   CreateAndFillUserTH1D("deltaETAjjHLT", 100, 0, 3, deltaETAjj );
	   CreateAndFillUserTH1D("deltaETAjjReco", 100, 0, 3, deltaETAjjreco );
	   CreateAndFillUserTH2D("deltaETAjj", 100, 0, 3, 100, 0, 3, deltaETAjjreco, deltaETAjj );
	   
	   double dR_min;
	   // HLT-Reco Matching : AK4
	   for(int ii=0; ii<2; ii++){
	     dR_min= AK4recojets[IdxMatching.at(ii)].DeltaR(AK4jets[ii]);
	     if(verbose) cout << "dR_min " << dR_min << endl;
	     CreateAndFillUserTH1D("AK4_deltaR_minimum", 500, 0, 5, dR_min);
	     
	     if(dR_min< DeltaR_ ){   
	       CreateAndFillUserTH1D("AK4_deltaR_minimum_AfterMatch", 500, 0, 5, dR_min);
	       if(verbose){
		 std::cout<<"AK4"<<std::endl;
		 std::cout<<"Variabili del jets # "<< ii <<std::endl;
		 std::cout<<"Reco"<<std::endl;
		 std::cout<<"Pt= "<< AK4recojets[IdxMatching.at(ii)].Pt()<<" Eta= "<<AK4recojets[IdxMatching.at(ii)].Eta()<< " Phi= "<<AK4recojets[IdxMatching.at(ii)].Phi()  <<std::endl;
		 std::cout<<"HLT" <<std::endl;
		 std::cout<<"Pt= "<< AK4jets[ii].Pt()<<" Eta= "<<AK4jets[ii].Eta()<< " Phi= "<<AK4jets[ii].Phi()  <<std::endl;
	       }	 
	       
	       CreateAndFillUserTH1D("chargedHadEnFracReco", 100, 0, 1, chargedHadEnFracReco[IdxMatching.at(ii)] );
	       CreateAndFillUserTH1D("neutrHadEnFracReco", 100, 0, 1, neutrHadEnFracReco[IdxMatching.at(ii)] );
	       CreateAndFillUserTH1D("chargedElectromFracReco", 100, 0, 1, chargedElectromFracReco[IdxMatching.at(ii)] );
	       CreateAndFillUserTH1D("neutrElectromFracReco", 100, 0, 1, neutrElectromFracReco[IdxMatching.at(ii)] );
	       CreateAndFillUserTH1D("muEnFracReco", 100, 0, 1, muEnFracReco[IdxMatching.at(ii)] );
	       CreateAndFillUserTH1D("chargedMultReco", 100, 0, 100, chargedMultReco[IdxMatching.at(ii)] );
	       CreateAndFillUserTH1D("neutrMultReco", 100, 0, 100, neutrMultReco[IdxMatching.at(ii)] );
	       CreateAndFillUserTH1D("photonMultReco", 100, 0, 100, photonMultReco[IdxMatching.at(ii)] );
	       CreateAndFillUserTH1D("jetCSVAK4Reco", 100, 0, 2, jetCSVAK4Reco[IdxMatching.at(ii)] );
	     
	       CreateAndFillUserTH1D("chargedHadEnFracHLT", 100, 0, 1, chargedHadEnFrac[ii] );
	       CreateAndFillUserTH1D("neutrHadEnFracHLT", 100, 0, 1, neutrHadEnFrac[ii] );
	       CreateAndFillUserTH1D("chargedElectromFracHLT", 100, 0, 1, chargedElectromFrac[ii] );
	       CreateAndFillUserTH1D("neutrElectromFracHLT", 100, 0, 1, neutrElectromFrac[ii] );
	       CreateAndFillUserTH1D("muEnFracHLT", 100, 0, 1, muEnFrac[ii] );
	       CreateAndFillUserTH1D("chargedMultHLT", 100, 0, 100, chargedMult[ii] );
	       CreateAndFillUserTH1D("neutrMultHLT", 100, 0, 100, neutrMult[ii] );
	       CreateAndFillUserTH1D("photonMultHLT", 100, 0, 100, photonMult[ii] );
	       CreateAndFillUserTH1D("jetCSVAK4HLT", 100, 0, 2, jetCSVAK4[ii] );
	     
	       CreateAndFillUserTH2D("chargedHadEnFrac", 100, 0, 1, 100, 0, 1, chargedHadEnFracReco[IdxMatching.at(ii)], chargedHadEnFrac[ii] );
	       CreateAndFillUserTH2D("neutrHadEnFrac", 100, 0, 1, 100, 0, 1, neutrHadEnFracReco[IdxMatching.at(ii)], neutrHadEnFrac[ii] );
	       CreateAndFillUserTH2D("chargedElectromFrac", 100, 0, 1, 100, 0, 1, chargedElectromFracReco[IdxMatching.at(ii)], chargedElectromFrac[ii] );
	       CreateAndFillUserTH2D("neutrElectromFrac", 100, 0, 1, 100, 0, 1, neutrElectromFracReco[IdxMatching.at(ii)], neutrElectromFrac[ii] );
	       CreateAndFillUserTH2D("muEnFrac", 100, 0, 1, 100, 0, 1, muEnFracReco[IdxMatching.at(ii)], muEnFrac[ii] );
	       CreateAndFillUserTH2D("chargedMult", 100, 0, 100, 100, 0, 100, chargedMultReco[IdxMatching.at(ii)], chargedMult[ii] );
	       CreateAndFillUserTH2D("neutrMult", 100, 0, 100, 100, 0, 100, neutrMultReco[IdxMatching.at(ii)], neutrMult[ii] );
	       CreateAndFillUserTH2D("photonMult", 100, 0, 100, 100, 0, 100, photonMultReco[IdxMatching.at(ii)], photonMult[ii] );
	       CreateAndFillUserTH2D("jetCSVAK4", 100, 0, 2, 100, 0, 2, jetCSVAK4Reco[IdxMatching.at(ii)], jetCSVAK4[ii] );

	       CreateAndFillUserTH2D("AK4_pt_Jets", 300, 0, 3000, 300, 0, 3000, AK4recojets[IdxMatching.at(ii)].Pt(), AK4jets[ii].Pt() );

	       //pT(Reco) - pT(HLT) / pT(Reco)	   
	       double pTResolution = ( AK4recojets[IdxMatching.at(ii)].Pt() - AK4jets[ii].Pt() )/ AK4recojets[IdxMatching.at(ii)].Pt(); 
	       CreateAndFillUserTH1D("AK4_pTResolution", 300, -2, 2, pTResolution);

	       ////// TProfile: pTResolution vs nVtx

	       AK4_pTResolution_vs_Nvtx->Fill(nVtxreco, pTResolution);   	     
	       AK4_pT_vs_NvtxReco->Fill(nVtxreco, AK4recojets[IdxMatching.at(ii)].Pt() );   	     
	       AK4_pT_vs_NvtxHLT ->Fill(nVtx, AK4jets[ii].Pt() );   	     

	       //////////////////// analysis in bin di pT e Eta	     
	       EtaBinning mEtaBinning;
	       PtBinning mPtBinning;
	     
	       int etaBin = mEtaBinning.getBin( fabs(AK4recojets[IdxMatching.at(ii)].Eta() ) );
	       int ptBin = mPtBinning.getPtBin( AK4recojets[IdxMatching.at(ii)].Pt() );
	     
	       std::pair<float, float> ptBins = mPtBinning.getBinValue(ptBin);	     
	     
	       if(verbose){
		 std::pair<float, float> etaBins = mEtaBinning.getBinValue(etaBin);
		 cout<<"etaReco  "<< fabs(AK4recojets[IdxMatching.at(ii)].Eta() )  << "    pT Reco   "<<  AK4recojets[IdxMatching.at(ii)].Pt() <<endl;
		 cout<<"etaBin  "<< etaBin << "    pTBin   "<<ptBin<<endl;
		 cout<<"etaBin.first  "<< etaBins.first << "    etaBin.second   "<<etaBins.second<<endl;
		 cout<<"ptBin.first  "<< ptBins.first << "    ptBin.second   "<<ptBins.second<<endl;
	       }
	     
	       std::string etaName = mEtaBinning.getBinName(etaBin);
	     
	       std::string HistoName = TString::Format("AK4_pTResolution_%s_pT_%i_%i", etaName.c_str(), (int) ptBins.first, (int) ptBins.second ).Data();
	       std::string HistoName2 = TString::Format("AK4_pt_Jets_%s_pT_%i_%i", etaName.c_str(), (int) ptBins.first, (int) ptBins.second ).Data();

	       if(verbose) std::cout << HistoName.c_str()<< std::endl;

	       CreateAndFillUserTH1D(HistoName.c_str(), 300, -2, 2, pTResolution);
	       CreateAndFillUserTH2D(HistoName2.c_str(), 300, 0, 3000, 300, 0, 3000, AK4recojets[IdxMatching.at(ii)].Pt(), AK4jets[ii].Pt() );

	        
	     }//if matching
	   }// for on number of jets
	 }// if number of jets != 0
       
	 ////// End AK4 analysis
	 ////// Start WIDEJETS
	 vector<TLorentzVector> widejetsReco;
	 TLorentzVector wj1Reco, wj2Reco, wdijetReco; 
       
	 vector<TLorentzVector> widejets;
	 TLorentzVector wj1, wj2, wdijet; 
       
	 //Reco
	 wj1Reco.SetPtEtaPhiM(pTWJ_recoj1, etaWJ_recoj1, phiWJ_recoj1, massWJ_recoj1);
	 wj2Reco.SetPtEtaPhiM(pTWJ_recoj2, etaWJ_recoj2, phiWJ_recoj2, massWJ_recoj2);
       
	 if( wj1Reco.Pt()>0 && wj2Reco.Pt()>0 )
	   {
	     // Put widejets in the container
	     widejetsReco.push_back( wj1Reco );
	     widejetsReco.push_back( wj2Reco );
	   }
       
	 //HLT
	 wj1.SetPtEtaPhiM(pTWJ_j1, etaWJ_j1, phiWJ_j1, massWJ_j1);
	 wj2.SetPtEtaPhiM(pTWJ_j2, etaWJ_j2, phiWJ_j2, massWJ_j2);
       
	 if( wj1.Pt()>0 && wj2.Pt()>0 )
	   {
	     // Put widejets in the container
	     widejets.push_back( wj1 );
	     widejets.push_back( wj2 );
	   }
       
	 // Compare with WIDEJETS
	 if( widejetsReco.size() !=0  && widejets.size() !=0 ) {
	 
	   //////////// matching	 
	   double DeltaReco1HLT1 = widejetsReco[0].DeltaR(widejets[0]);
	   double DeltaReco1HLT2 = widejetsReco[0].DeltaR(widejets[1]);
	   double DeltaReco2HLT1 = widejetsReco[1].DeltaR(widejets[0]);
	   double DeltaReco2HLT2 = widejetsReco[1].DeltaR(widejets[1]);
	 
	   if(verbose){
	     //	   cout << endl;
	     cout << "widejets " << endl;
	     cout << "deltaR Reco1-HLT1 " <<  DeltaReco1HLT1 << endl;
	     cout << "deltaR Reco1-HLT2 " <<  DeltaReco1HLT2 << endl;
	     cout << "deltaR Reco2-HLT1 " <<  DeltaReco2HLT1 << endl;
	     cout << "deltaR Reco2-HLT2 " <<  DeltaReco2HLT2 << endl;
	   }       
	 
	   vector<TLorentzVector> widejetsReco_SortedMatching;
	   if( DeltaReco1HLT1 <= DeltaReco2HLT1 )
	     {
	       widejetsReco_SortedMatching.push_back( widejetsReco[0] );
	       widejetsReco_SortedMatching.push_back( widejetsReco[1] );	 
	     }else
	     {
	       widejetsReco_SortedMatching.push_back( widejetsReco[1] );
	       widejetsReco_SortedMatching.push_back( widejetsReco[0] );
	     }
	 
	   //       std::cout<<"DeltaR_ "<< DeltaR_ <<std::endl; // test - works
	 
	   double dR_min;
	   //il matching e' molto buono
	   for(int ii=0; ii<2; ii++){
	     dR_min= widejetsReco_SortedMatching[ii].DeltaR(widejets[ii]);
	     if(verbose) cout << "dR_min " << dR_min << endl;
	     CreateAndFillUserTH1D("deltaR_minimum", 500, 0, 5, dR_min);
	   
	     if(dR_min< DeltaR_ ){
	       CreateAndFillUserTH1D("deltaR_minimum_AfterMatch", 500, 0, 5, dR_min);
	       if(verbose){
		 std::cout<<"WIDEJETS" <<std::endl;
		 std::cout<<"Variabili del jets # "<< ii <<std::endl;
		 std::cout<<"Reco"<<std::endl;
		 std::cout<<"Pt= "<< widejetsReco_SortedMatching[ii].Pt()<<" Eta= "<<widejetsReco_SortedMatching[ii].Eta()<< " Phi= "<<widejetsReco_SortedMatching[ii].Phi()  <<std::endl;
		 std::cout<<"HLT" <<std::endl;
		 std::cout<<"Pt= "<< widejets[ii].Pt()<<" Eta= "<<widejets[ii].Eta()<< " Phi= "<<widejets[ii].Phi()  <<std::endl;
	       }
	      
	       // no bins, all jets matched	     
	       CreateAndFillUserTH1D("pt_WideJetHLT", 300, 0, 3000, widejets[ii].Pt()  );
	       CreateAndFillUserTH1D("eta_WideJetHLT", 100, -2.5, 2.5, widejets[ii].Eta()  );
	       CreateAndFillUserTH1D("phi_WideJetHLT", 100, -3.5, 3.5, widejets[ii].Phi()  );
	       CreateAndFillUserTH1D("MJJ_WideJetHLT", 300, 0, 3000, mjj );

	       CreateAndFillUserTH1D("pt_WideJetReco", 300, 0, 3000, widejetsReco_SortedMatching[ii].Pt()  );
	       CreateAndFillUserTH1D("eta_WideJetReco", 100, -2.5, 2.5, widejetsReco_SortedMatching[ii].Eta() );
	       CreateAndFillUserTH1D("phi_WideJetReco", 100, -3.5, 3.5, widejetsReco_SortedMatching[ii].Phi()  );
	       CreateAndFillUserTH1D("MJJ_WideJetReco", 300, 0, 3000, mjjreco );
	     
	       CreateAndFillUserTH2D("pt_WideJet", 300, 0, 3000, 300, 0, 3000, widejetsReco_SortedMatching[ii].Pt(), widejets[ii].Pt() );
	       CreateAndFillUserTH2D("eta_WideJet", 100, -2.5, 2.5, 100, -2.5, 2.5, widejetsReco_SortedMatching[ii].Eta(), widejets[ii].Eta() );
	       CreateAndFillUserTH2D("phi_WideJet", 100, -3.5, 3.5, 100, -3.5, 3.5, widejetsReco_SortedMatching[ii].Phi(), widejets[ii].Phi() );
	       CreateAndFillUserTH2D("MJJ_WideJet", 300, 0, 3000, 300, 0, 3000, mjjreco, mjj );
	       CreateAndFillUserTH2D("FabsPhi_WideJet", 50, 0, 3.5, 50, 0, 3.5, fabs(widejetsReco_SortedMatching[ii].Phi()), fabs(widejets[ii].Phi()) );
	     
	     
	       //pT(Reco) - pT(HLT) / pT(Reco)	   
	       double pTResolution = ( widejetsReco_SortedMatching[ii].Pt() - widejets[ii].Pt() )/ widejetsReco_SortedMatching[ii].Pt(); 
	       CreateAndFillUserTH1D("pTResolution", 300, -2, 2, pTResolution);
	  
	       ////// TProfile: pTResolution vs nVtx

	       pTResolution_vs_Nvtx->Fill(nVtxreco, pTResolution);   	     
	       pT_vs_NvtxReco->Fill(nVtxreco, widejetsReco_SortedMatching[ii].Pt() );   	     
	       pT_vs_NvtxHLT ->Fill(nVtx, widejets[ii].Pt() );   	     

	       //////////// divide barrel and endcap

	       if(fabs(widejetsReco_SortedMatching[ii].Eta())<= 1.5 ){
		 CreateAndFillUserTH1D("pTResolution_Barrel", 300, -2, 2, pTResolution);
		 pTResolution_vs_Nvtx_Barrel->Fill(nVtxreco, pTResolution);   	     
		 pT_vs_Nvtx_BarrelReco->Fill(nVtxreco, widejetsReco_SortedMatching[ii].Pt() );   	     
		 pT_vs_Nvtx_BarrelHLT ->Fill(nVtx, widejets[ii].Pt() );   	     
	       }else if( fabs(widejetsReco_SortedMatching[ii].Eta())> 1.5 && fabs(widejetsReco_SortedMatching[ii].Eta())<= 2.5  ){
		 CreateAndFillUserTH1D("pTResolution_Endcap", 300, -2, 2, pTResolution);
		 pTResolution_vs_Nvtx_Endcap->Fill(nVtxreco, pTResolution);   	     
		 pT_vs_Nvtx_EndcapReco->Fill(nVtxreco, widejetsReco_SortedMatching[ii].Pt() );   	     
		 pT_vs_Nvtx_EndcapHLT ->Fill(nVtx, widejets[ii].Pt() );   	     
	       }

	       //////////////////// analysis in bin di pT e Eta	     
	       EtaBinning mEtaBinning;
	       PtBinning mPtBinning;
	     
	       int etaBin = mEtaBinning.getBin( fabs(widejetsReco_SortedMatching[ii].Eta() ) );
	       int ptBin = mPtBinning.getPtBin( widejetsReco_SortedMatching[ii].Pt() );
	     
	       std::pair<float, float> ptBins = mPtBinning.getBinValue(ptBin);	     
	     
	       if(verbose){
		 std::pair<float, float> etaBins = mEtaBinning.getBinValue(etaBin);
		 cout<<"etaReco  "<< fabs(widejetsReco_SortedMatching[ii].Eta() )  << "    pT Reco   "<<  widejetsReco_SortedMatching[ii].Pt() <<endl;
		 cout<<"etaBin  "<< etaBin << "    pTBin   "<<ptBin<<endl;
		 cout<<"etaBin.first  "<< etaBins.first << "    etaBin.second   "<<etaBins.second<<endl;
		 cout<<"ptBin.first  "<< ptBins.first << "    ptBin.second   "<<ptBins.second<<endl;
	       }
	     
	       std::string etaName = mEtaBinning.getBinName(etaBin);
	     
	       std::string HistoName = TString::Format("pTResolution_%s_pT_%i_%i", etaName.c_str(), (int) ptBins.first, (int) ptBins.second ).Data();
	       std::string HistoName2 = TString::Format("pt_WideJet_%s_pT_%i_%i", etaName.c_str(), (int) ptBins.first, (int) ptBins.second ).Data();
	     
	       if(verbose) std::cout << HistoName.c_str()<< std::endl;
	     
	       CreateAndFillUserTH1D(HistoName.c_str(), 300, -2, 2, pTResolution);
	       CreateAndFillUserTH2D(HistoName2.c_str(), 300, 0, 3000, 300, 0, 3000, widejetsReco_SortedMatching[ii].Pt(), widejets[ii].Pt() );
	     	     
	     }// if matching
	   }// for on jets number
	 }// numbers of jets != 0
       }// mjj > 500
     } // pass Json
     /////////////////////////////////////////////////////////////////
     
     //== Fill Variables ==
     fillVariableWithValue("PassJSON", PassJSON);     
     fillVariableWithValue("MJJ_reco", mjjreco);     
     //     fillVariableWithValue("dR_min", deltaR[0]);          
    
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
   pTResolution_vs_Nvtx -> Write();
   pT_vs_NvtxReco          -> Write();
   pT_vs_NvtxHLT           -> Write();

   pTResolution_vs_Nvtx_Barrel -> Write();
   pT_vs_Nvtx_BarrelReco          -> Write();
   pT_vs_Nvtx_BarrelHLT           -> Write();

   pTResolution_vs_Nvtx_Endcap -> Write();
   pT_vs_Nvtx_EndcapReco          -> Write();
   pT_vs_Nvtx_EndcapHLT           -> Write();

   AK4_pTResolution_vs_Nvtx -> Write();
   AK4_pT_vs_NvtxReco          -> Write();
   AK4_pT_vs_NvtxHLT            -> Write();

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
