#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>

bool verbose = false;
int Counter = 0;   // debugging -- can be erased    
int CounterPtMore300 = 0;   // debugging -- can be erased    
int CounterPtLess300 = 0;   // debugging -- can be erased    
double DeltaR_;

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile)
{
  std::cout << "analysisClass::analysisClass(): begins " << std::endl;

  //  std::string jetAlgo = getPreCutString1("jetAlgo");
  //  double rParam = getPreCutValue1("DeltaR");

  DeltaR_ = getPreCutValue1("DeltaR");

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
     std::cout<<"Widejet 1"<< std::endl;
     std::cout<<"Pt: "<<pTWJ_recoj1<<" eta: "<<etaWJ_recoj1<<" phi: "<<phiWJ_recoj1<< std::endl;
     std::cout<<"Widejet 2"<< std::endl;
     std::cout<<"Pt: "<<pTWJ_recoj2<<" eta: "<<etaWJ_recoj2<<" phi: "<<phiWJ_recoj2<< std::endl;
     std::cout<<"HLT"<< std::endl;
     std::cout<<"Widejet 1"<< std::endl;
     std::cout<<"Pt: "<<pTWJ_j1<<" eta: "<<etaWJ_j1<<" phi: "<<phiWJ_j1<< std::endl;
     std::cout<<"Widejet 2"<< std::endl;
     std::cout<<"Pt: "<<pTWJ_j2<<" eta: "<<etaWJ_j2<<" phi: "<<phiWJ_j2<< std::endl;
     }


     if(!PassJSON){

       vector<TLorentzVector> widejetsReco;
       TLorentzVector wj1Reco, wj2Reco, wdijetReco; 
       
       vector<TLorentzVector> widejets;
       TLorentzVector wj1, wj2, wdijet; 
       
       //Reco
       wj1Reco.SetPtEtaPhiM(pTWJ_recoj1, etaWJ_recoj1, phiWJ_recoj1, massWJ_recoj1);
       wj2Reco.SetPtEtaPhiM(pTWJ_recoj2, etaWJ_recoj2, phiWJ_recoj2, massWJ_recoj2);
       
       //     double MJJWideReco = 0; 
       //     double DeltaEtaJJWideReco = 0;
       //     double DeltaPhiJJWideReco = 0;
       if( wj1Reco.Pt()>0 && wj2Reco.Pt()>0 )
       {
	 // Create dijet system
	 //	 wdijetReco = wj1Reco + wj2Reco;
	 //	 MJJWideReco = wdijetReco.M();
	 //	 DeltaEtaJJWideReco = fabs(wj1Reco.Eta()-wj2Reco.Eta());
	 //	 DeltaPhiJJWideReco = fabs(wj1Reco.DeltaPhi(wj2Reco));
	 
	 // Put widejets in the container
	 widejetsReco.push_back( wj1Reco );
	 widejetsReco.push_back( wj2Reco );
       }

       // confronto mjj ricalcolata e mjj nell'ntupla 
       //     std::cout<<"MJJWideReco: "<<MJJWideReco<<" mjjreco: "<<mjjreco<< std::endl; 
       // sono uguali
       
       //HLT
       wj1.SetPtEtaPhiM(pTWJ_j1, etaWJ_j1, phiWJ_j1, massWJ_j1);
       wj2.SetPtEtaPhiM(pTWJ_j2, etaWJ_j2, phiWJ_j2, massWJ_j2);
       
       //     double MJJWide = 0; 
       //     double DeltaEtaJJWide = 0;
       //     double DeltaPhiJJWide = 0;
       if( wj1.Pt()>0 && wj2.Pt()>0 )
	 {
	   // Create dijet system
	   //	 wdijet = wj1 + wj2;
	   //	 MJJWide = wdijet.M();
	   //	 DeltaEtaJJWide = fabs(wj1.Eta()-wj2.Eta());
	   //	 DeltaPhiJJWide = fabs(wj1.DeltaPhi(wj2));
	   
	   // Put widejets in the container
	   widejets.push_back( wj1 );
	   widejets.push_back( wj2 );
	 }
       
       // confronto mjj ricalcolata e mjj nell'ntupla 
       //     std::cout<<"MJJWide: "<<MJJWide<<" mjj: "<<mjj<< std::endl;
       // sono uguali
       
       // Compare with WIDEJETS
       if( widejetsReco.size() !=0  && widejets.size() !=0 ) { // forse non dovrebbe essere piu necessario -- variabili riempite nell'analizer se vector size >1 || >2
	 //       double R_min = 1000;
	 
	 double DeltaReco1HLT1 = widejetsReco[0].DeltaR(widejets[0]);
	 double DeltaReco1HLT2 = widejetsReco[0].DeltaR(widejets[1]);
	 double DeltaReco2HLT1 = widejetsReco[1].DeltaR(widejets[0]);
	 double DeltaReco2HLT2 = widejetsReco[1].DeltaR(widejets[1]);
       
	 if(verbose){
	   cout << endl;
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
	     if(verbose){
	       std::cout<<"Variabili del jets # "<< ii <<std::endl;
	       std::cout<<"Reco"<<std::endl;
	       std::cout<<"Pt= "<< widejetsReco_SortedMatching[ii].Pt()<<" Eta= "<<widejetsReco_SortedMatching[ii].Eta()<< " Phi= "<<widejetsReco_SortedMatching[ii].Phi()  <<std::endl;
	       std::cout<<"HLT" <<std::endl;
	       std::cout<<"Pt= "<< widejets[ii].Pt()<<" Eta= "<<widejets[ii].Eta()<< " Phi= "<<widejets[ii].Phi()  <<std::endl;
	     }
	     
	     CreateAndFillUserTH2D("ptWJ_Reco_vs_HLT", 300, 0, 3000, 300, 0, 3000, widejetsReco_SortedMatching[ii].Pt(), widejets[ii].Pt() );
	     CreateAndFillUserTH2D("EtaWJ_Reco_vs_HLT", 100, -2.5, 2.5, 100, -2.5, 2.5, widejetsReco_SortedMatching[ii].Eta(), widejets[ii].Eta() );
	     CreateAndFillUserTH2D("PhiWJ_Reco_vs_HLT", 100, -3.5, 3.5, 100, -3.5, 3.5, widejetsReco_SortedMatching[ii].Phi(), widejets[ii].Phi() );
	     CreateAndFillUserTH2D("FabsPhiWJ_Reco_vs_HLT", 50, 0, 3.5, 50, 0, 3.5, fabs(widejetsReco_SortedMatching[ii].Phi()), fabs(widejets[ii].Phi()) );
	     CreateAndFillUserTH2D("MJJWide_Reco_vs_HLT", 300, 0, 3000, 300, -0, 3000, mjjreco, mjj ); // cambia variabile
	     // CreateAndFillUserTH2D("Nvtx_Reco_vs_HLT", 51, 0, 50, 51, 0, 50, nvtxreco, nvtx ); //cambia variabile
	     //pT(Reco) - pT(HLT) / pT(Reco)	   
	     double pTResolution = ( widejetsReco_SortedMatching[ii].Pt() - widejets[ii].Pt() )/ widejetsReco_SortedMatching[ii].Pt(); 
	     CreateAndFillUserTH1D("pTResolution", 300, -2, 2, pTResolution);
	     
	     if(verbose){
	       if( widejetsReco_SortedMatching[ii].Pt() > 300){
		 CounterPtMore300++;
		 std::cout<<"numero di eventi con pT > 300 = "<<CounterPtMore300<<std::endl;
	       }else{
		 CounterPtLess300++;
		 std::cout<<"numero di eventi con pT < 300 = "<<CounterPtLess300<<std::endl;
	       }	   
	     }
	     
	     if( widejetsReco_SortedMatching[ii].Pt() > 300 && fabs(widejetsReco_SortedMatching[ii].Eta()) < 1.5 ){
	       CreateAndFillUserTH2D("ptWJ_Reco_vs_HLT_Cut1", 300, 0, 3000, 300, 0, 3000, widejetsReco_SortedMatching[ii].Pt(), widejets[ii].Pt() );
	       CreateAndFillUserTH2D("EtaWJ_Reco_vs_HLT_Cut1", 100, -2.5, 2.5, 100, -2.5, 2.5, widejetsReco_SortedMatching[ii].Eta(), widejets[ii].Eta() );
	       CreateAndFillUserTH2D("PhiWJ_Reco_vs_HLT_Cut1", 100, -3.5, 3.5, 100, -3.5, 3.5, widejetsReco_SortedMatching[ii].Phi(), widejets[ii].Phi() );
	       CreateAndFillUserTH2D("FabsPhiWJ_Reco_vs_HLT_Cut1", 50, 0, 3.5, 50, 0, 3.5, fabs(widejetsReco_SortedMatching[ii].Phi()), fabs(widejets[ii].Phi()) );
	       CreateAndFillUserTH2D("MJJWide_Reco_vs_HLT_Cut1", 300, 0, 3000, 300, -0, 3000, mjjreco, mjj );
	       //   CreateAndFillUserTH2D("Nvtx_Reco_vs_HLT_Cut1", 51, 0, 50, 51, 0, 50, nvtxreco, nvtx );
	       CreateAndFillUserTH1D("pTResolution_Cut1", 300, -2, 2, pTResolution);
	     }
	     
	     if( widejetsReco_SortedMatching[ii].Pt() > 300 && fabs(widejetsReco_SortedMatching[ii].Eta()) >= 1.5 && fabs(widejetsReco_SortedMatching[ii].Eta()) < 2.5 ){
	       CreateAndFillUserTH2D("ptWJ_Reco_vs_HLT_Cut2", 300, 0, 3000, 300, 0, 3000, widejetsReco_SortedMatching[ii].Pt(), widejets[ii].Pt() );
	       CreateAndFillUserTH2D("EtaWJ_Reco_vs_HLT_Cut2", 100, -2.5, 2.5, 100, -2.5, 2.5, widejetsReco_SortedMatching[ii].Eta(), widejets[ii].Eta() );
	       CreateAndFillUserTH2D("PhiWJ_Reco_vs_HLT_Cut2", 100, -3.5, 3.5, 100, -3.5, 3.5, widejetsReco_SortedMatching[ii].Phi(), widejets[ii].Phi() );
	       CreateAndFillUserTH2D("FabsPhiWJ_Reco_vs_HLT_Cut2", 50, 0, 3.5, 50, 0, 3.5, fabs(widejetsReco_SortedMatching[ii].Phi()), fabs(widejets[ii].Phi()) );
	       CreateAndFillUserTH2D("MJJWide_Reco_vs_HLT_Cut2", 300, 0, 3000, 300, -0, 3000, mjjreco, mjj );
	       //	   CreateAndFillUserTH2D("Nvtx_Reco_vs_HLT_Cut2", 51, 0, 50, 51, 0, 50, nvtxreco, nvtx );
	   CreateAndFillUserTH1D("pTResolution_Cut2", 300, -2, 2, pTResolution);
	     }
	     
	   }
	   /*else{//dR cut
	     std::cout<<"dR minimo troppo grande = "<<dR_min<<std::endl;
	     Counter++;
	     std::cout<<"numero di eventi = "<<Counter<<std::endl;
	     }*/ // debugging -- ok

	 }// for on jets number
       }// numbers of jets != 0
     }
     /////////////////////////////////////////////////////////////////
     
     //== Fill Variables ==
     fillVariableWithValue("PassJSON", PassJSON);     
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
