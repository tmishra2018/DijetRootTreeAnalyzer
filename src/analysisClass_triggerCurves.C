#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include "TGraphAsymmErrors.h"


analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile)
{
  std::cout << "analysisClass::analysisClass(): begins " << std::endl;

  std::string jetAlgo = getPreCutString1("jetAlgo");
  double rParam = getPreCutValue1("DeltaR");

  if( jetAlgo == "AntiKt" )
    fjJetDefinition = JetDefPtr( new fastjet::JetDefinition(fastjet::antikt_algorithm, rParam) );
  else if( jetAlgo == "Kt" )
    fjJetDefinition = JetDefPtr( new fastjet::JetDefinition(fastjet::kt_algorithm, rParam) );
  else 
    fjJetDefinition = JetDefPtr( new fastjet::JetDefinition(fastjet::cambridge_algorithm, rParam) );

  // For JECs
  if( int(getPreCutValue1("useJECs"))==1 )
  {
    std::string L1Path = "/afs/cern.ch/user/f/ferencek/public/JEC_txt_files/PHYS14_25_V2/PHYS14_25_V2_L1FastJet_AK4PFchs.txt";
    std::string L2Path = "/afs/cern.ch/user/f/ferencek/public/JEC_txt_files/PHYS14_25_V2/PHYS14_25_V2_L2Relative_AK4PFchs.txt"; 
    std::string L3Path = "/afs/cern.ch/user/f/ferencek/public/JEC_txt_files/PHYS14_25_V2/PHYS14_25_V2_L3Absolute_AK4PFchs.txt";
    
    L1Par = new JetCorrectorParameters(L1Path);
    L2Par = new JetCorrectorParameters(L2Path);
    L3Par = new JetCorrectorParameters(L3Path);

    std::vector<JetCorrectorParameters> vPar;
    vPar.push_back(*L1Par);
    vPar.push_back(*L2Par);
    vPar.push_back(*L3Par);

    JetCorrector = new FactorizedJetCorrector(vPar);
  }

///////////
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
   //giulia -- method to calculate uncertainties on trigger eff ?
   int MODE = 1;
   double scale(1.),prescale(1.);
    
   if (fChain == 0) return;
   
   //////////book histos here

   // TH1F *h_nJetFinal = new TH1F ("h_nJetFinal","",10,0,10);
   // h_nJetFinal->Sumw2();      
   // TH1F *h_nVtx = new TH1F ("h_nVtx","",30,0,30);
   // h_nVtx->Sumw2(); 
   // TH1F *h_trueVtx = new TH1F ("h_trueVtx","",40,0,40);
   // h_trueVtx->Sumw2();  
   // TH1F *h_pT1stJet = new TH1F ("h_pT1stJet","",100,0,3000);
   // h_pT1stJet->Sumw2();
   // TH1F *h_pT2ndJet = new TH1F ("h_pT2ndJet","",100,0,3000);
   // h_pT2ndJet->Sumw2();
   // TH1F *h_eta1stJet = new TH1F ("h_eta1stJet","",5,-2.5,2.5);
   // h_eta1stJet->Sumw2();
   // TH1F *h_eta2ndJet = new TH1F ("h_eta2ndJet","",5,-2.5,2.5);
   // h_eta2ndJet->Sumw2();
   // TH1F *h_DijetMass = new TH1F ("h_DijetMass","",600,0,6000);
   // h_DijetMass->Sumw2();
   // TH1F *h_DeltaETAjj = new TH1F ("h_DeltaETAjj","",120,0,3.);
   // h_DeltaETAjj->Sumw2();
   
   // variable binning for mjj trigger efficiency plots
   const int nMassBins = 103;

   double massBoundaries[nMassBins+1] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325,
     354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687,
     1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509,
     4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 
     10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000};

   //for TGraphAsymErr plots of efficiency
   double vx[nMassBins][3],vy[nMassBins][3],vexl[nMassBins][3],vexh[nMassBins][3],veyl[nMassBins][3],veyh[nMassBins][3];



   char* HLTname[50] = {"PFHT350","PFHT900","PFHT650MJJ900","PFHT900_AND_PFHT650MJJ900"};
   TH1F* h_mjj_HLTpass[4];
   char name_histoHLT[50];
   for (int i=0; i<4; i++){  
     sprintf(name_histoHLT,"h_mjj_HLTpass_%s",HLTname[i]);
     h_mjj_HLTpass[i]= new TH1F(name_histoHLT,"",103,massBoundaries);
   }
  
   /////////initialize variables

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<50;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     // if (Cut(ientry) < 0) continue;

     ////////////////////// User's code starts here ///////////////////////

     ///Stuff to be done for every event
 
     size_t no_jets_ak4=jetPtAK4->size();
     vector<TLorentzVector> widejets;
     TLorentzVector wj1, wj2, wdijet;

     resetCuts();

     std::vector<double> jecFactors;
     // new JECs could change the jet pT ordering. the vector below
     // holds sorted jet indices after the new JECs had been applied
     std::vector<unsigned> sortedJetIdx;
     if( int(getPreCutValue1("useJECs"))==1 )
     {
       // sort jets by increasing pT
       std::multimap<double, unsigned> sortedJets;
       for(size_t j=0; j<no_jets_ak4; ++j)
       {
	 JetCorrector->setJetEta(jetEtaAK4->at(j));
	 JetCorrector->setJetPt(jetPtAK4->at(j)/jetJecAK4->at(j));
	 JetCorrector->setJetA(jetAreaAK4->at(j));
	 JetCorrector->setRho(rho);

	 double correction = JetCorrector->getCorrection();

	 jecFactors.push_back(correction);
	 sortedJets.insert(std::make_pair((jetPtAK4->at(j)/jetJecAK4->at(j))*correction, j));
       }
       // get jet indices in decreasing pT order
       for(std::multimap<double, unsigned>::const_reverse_iterator it = sortedJets.rbegin(); it != sortedJets.rend(); ++it)
	 sortedJetIdx.push_back(it->second);
     }
     else
     {
       for(size_t j=0; j<no_jets_ak4; ++j)
       {
	 jecFactors.push_back(jetJecAK4->at(j));
	 sortedJetIdx.push_back(j);
       }
     }

     // if(no_jets_ak4>=2){
     //  if(!(fabs(jetEtaAK4->at(0)) < getPreCutValue1("jetFidRegion") && idTAK4->at(0) == getPreCutValue1("tightJetID"))){
     //    std::cout << " JET 0 FAIL " << jetEtaAK4->at(0) << " JET 0  ID " << idTAK4->at(0) << std::endl;
     //  }
     //  if(!(fabs(jetEtaAK4->at(1)) < getPreCutValue1("jetFidRegion") && idTAK4->at(1) == getPreCutValue1("tightJetID"))){
     //    std::cout << " JET 1 FAIL " << jetEtaAK4->at(1) << " JET 1  ID " << idTAK4->at(1) << std::endl;
     //  }  
     // }

     if( int(getPreCutValue1("useFastJet"))==1 )
     {
       // vector of ak4 jets used for wide jet clustering
       std::vector<fastjet::PseudoJet> fjInputs;

       for(size_t j=0; j<no_jets_ak4; ++j)
       {
	 if( !(jetEtaAK4->at(j) < getPreCutValue1("jetFidRegion")
	       && idTAK4->at(j) == getPreCutValue1("tightJetID")) ) continue;

	 if( j==0 && !((jecFactors[j]/jetJecAK4->at(sortedJetIdx[j]))*jetPtAK4->at(sortedJetIdx[j]) > getPreCutValue1("pt0Cut")) ) continue;
	 else if( j==1 && !((jecFactors[j]/jetJecAK4->at(sortedJetIdx[j]))*jetPtAK4->at(sortedJetIdx[j]) > getPreCutValue1("pt1Cut")) ) continue;
	 else if( !((jecFactors[j]/jetJecAK4->at(sortedJetIdx[j]))*jetPtAK4->at(sortedJetIdx[j]) > getPreCutValue1("ptCut")) ) continue;

	 TLorentzVector tempJet;
	 tempJet.SetPtEtaPhiM((jecFactors[j]/jetJecAK4->at(sortedJetIdx[j]))*jetPtAK4->at(sortedJetIdx[j]),jetEtaAK4->at(j),jetPhiAK4->at(j),jetMassAK4->at(j));

	 fjInputs.push_back(fastjet::PseudoJet(tempJet.Px(),tempJet.Py(),tempJet.Pz(),tempJet.E()));
       }

       fjClusterSeq = ClusterSequencePtr( new fastjet::ClusterSequence( fjInputs, *fjJetDefinition ) );

       std::vector<fastjet::PseudoJet> inclusiveWideJets = fastjet::sorted_by_pt( fjClusterSeq->inclusive_jets(0.) );

       if( inclusiveWideJets.size()>1 )
       {
	 wj1.SetPxPyPzE(inclusiveWideJets.at(0).px(), inclusiveWideJets.at(0).py(), inclusiveWideJets.at(0).pz(), inclusiveWideJets.at(0).e());
	 wj2.SetPxPyPzE(inclusiveWideJets.at(1).px(), inclusiveWideJets.at(1).py(), inclusiveWideJets.at(1).pz(), inclusiveWideJets.at(1).e());
       }
     }
     else
     {
       TLorentzVector wj1_tmp, wj2_tmp;
       double wideJetDeltaR_ = getPreCutValue1("DeltaR");

       if(no_jets_ak4>=2)
	 {
	   if(fabs(jetEtaAK4->at(0)) < getPreCutValue1("jetFidRegion") 
	      && (jecFactors[0]/jetJecAK4->at(sortedJetIdx[0]))*jetPtAK4->at(sortedJetIdx[0]) > getPreCutValue1("pt0Cut"))
	     {
	       if(fabs(jetEtaAK4->at(1)) < getPreCutValue1("jetFidRegion") 
		  && (jecFactors[1]/jetJecAK4->at(sortedJetIdx[1]))*jetPtAK4->at(sortedJetIdx[1]) > getPreCutValue1("pt1Cut"))
		 {
		   TLorentzVector jet1, jet2;
		   jet1.SetPtEtaPhiM((jecFactors[0]/jetJecAK4->at(sortedJetIdx[0]))*jetPtAK4->at(sortedJetIdx[0]),jetEtaAK4->at(0),jetPhiAK4->at(0),jetMassAK4->at(0));
		   jet2.SetPtEtaPhiM((jecFactors[1]/jetJecAK4->at(sortedJetIdx[1]))*jetPtAK4->at(sortedJetIdx[1]),jetEtaAK4->at(1),jetPhiAK4->at(1),jetMassAK4->at(1));
		   
		   for(Long64_t ijet=0; ijet<no_jets_ak4; ijet++)
		     { //jet loop for ak4
		       TLorentzVector currentJet;
		       
		       if(fabs(jetEtaAK4->at(ijet)) < getPreCutValue1("jetFidRegion") 
			  && idTAK4->at(ijet) == getPreCutValue1("tightJetID") 
			  && (jecFactors[ijet]/jetJecAK4->at(sortedJetIdx[ijet]))*jetPtAK4->at(sortedJetIdx[ijet]) > getPreCutValue1("ptCut"))
			 {
			   TLorentzVector currentJet;
			   currentJet.SetPtEtaPhiM((jecFactors[ijet]/jetJecAK4->at(sortedJetIdx[ijet]))*jetPtAK4->at(sortedJetIdx[ijet]),jetEtaAK4->at(ijet),jetPhiAK4->at(ijet),jetMassAK4->at(ijet));   
			   
			   double DeltaR1 = currentJet.DeltaR(jet1);
			   double DeltaR2 = currentJet.DeltaR(jet2);
			   
			   if(DeltaR1 < DeltaR2 && DeltaR1 < wideJetDeltaR_)
			     {
			       wj1_tmp += currentJet;
			     }
			   else if(DeltaR2 < wideJetDeltaR_)
			     {
			       wj2_tmp += currentJet;
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
	 }
       else
	 {
	   wj1 = wj2_tmp;
	   wj2 = wj1_tmp;
	 }
     }

     double MJJWide = 0; 
     double DeltaEtaJJWide = 0;
     double DeltaPhiJJWide = 0;
     if( wj1.Pt()>0 && wj2.Pt()>0 )
     {
       // Create dijet system
       wdijet = wj1 + wj2;
       MJJWide = wdijet.M();
       DeltaEtaJJWide = fabs(wj1.Eta()-wj2.Eta());
       DeltaPhiJJWide = fabs(wj1.DeltaPhi(wj2));

       // Put widejets in the container
       widejets.push_back( wj1 );
       widejets.push_back( wj2 );
     }


     //== Fill Variables ==

     fillVariableWithValue("run",runNo);     
     fillVariableWithValue("event",evtNo);     
     fillVariableWithValue("lumi",lumi);     
     fillVariableWithValue("nVtx",nvtx);     
     fillVariableWithValue("nJet",widejets.size());
     fillVariableWithValue("metSig",metSig);

     // Trigger
     int NtriggerBits = triggerResult->size();
     if( NtriggerBits > 0)
       fillVariableWithValue("passHLT",triggerResult->at(0));// HLT_PFHT900_v*    

     if( no_jets_ak4 >=1 ){
       fillVariableWithValue( "IdTight_j1",idTAK4->at(0));
       fillVariableWithValue( "pTAK4_j1", jetPtAK4->at(0) );
       fillVariableWithValue( "etaAK4_j1", jetEtaAK4->at(0));
       fillVariableWithValue( "phiAK4_j1", jetPhiAK4->at(0));
     }
     if( no_jets_ak4 >=2 ){
       fillVariableWithValue( "IdTight_j2",idTAK4->at(1));
       fillVariableWithValue( "pTAK4_j2", jetPtAK4->at(1) );
       fillVariableWithValue( "etaAK4_j2", jetEtaAK4->at(1));
       fillVariableWithValue( "phiAK4_j2", jetPhiAK4->at(1));
       fillVariableWithValue( "Dijet_MassAK4", mjjAK4) ; 
       fillVariableWithValue( "CosThetaStarAK4", TMath::TanH( (jetEtaAK4->at(0)-jetEtaAK4->at(1))/2 )); 
     }

     if( widejets.size() >= 1 ){
         fillVariableWithValue( "pTWJ_j1", widejets[0].Pt() );
         fillVariableWithValue( "etaWJ_j1", widejets[0].Eta());

	 //no cuts on these variables, just to store in output
         fillVariableWithValue( "phiWJ_j1", widejets[0].Phi());
         fillVariableWithValue( "neutrHadEnFrac_j1", jetNhfAK4->at(0));
         fillVariableWithValue( "chargedHadEnFrac_j1", jetChfAK4->at(0));
         fillVariableWithValue( "photonEnFrac_j1", jetPhfAK4->at(0));
         fillVariableWithValue( "eleEnFract_j1", jetElfAK4->at(0));
         fillVariableWithValue( "muEnFract_j1", jetMufAK4->at(0));
         fillVariableWithValue( "chargedMult_j1", chMultAK4->at(0));
         fillVariableWithValue( "neutrMult_j1", neMultAK4->at(0));
         fillVariableWithValue( "photonMult_j1", phoMultAK4->at(0));
       }

     if( widejets.size() >= 2 ){
         fillVariableWithValue( "pTWJ_j2", widejets[1].Pt() );
         fillVariableWithValue( "etaWJ_j2", widejets[1].Eta());
	 fillVariableWithValue( "deltaETAjj", DeltaEtaJJWide ) ;
         fillVariableWithValue( "mjj", MJJWide ) ;
         fillVariableWithValue( "CosThetaStarWJ", TMath::TanH( (widejets[0].Eta()-widejets[1].Eta())/2 )); 

	 //no cuts on these variables, just to store in output
         fillVariableWithValue( "phiWJ_j2", widejets[1].Phi());	
         fillVariableWithValue( "neutrHadEnFrac_j2", jetNhfAK4->at(1));
         fillVariableWithValue( "chargedHadEnFrac_j2", jetChfAK4->at(1));
         fillVariableWithValue( "photonEnFrac_j2", jetPhfAK4->at(1));
         fillVariableWithValue( "eleEnFract_j2", jetElfAK4->at(1));
         fillVariableWithValue( "muEnFract_j2", jetMufAK4->at(1));
         fillVariableWithValue( "chargedMult_j2", chMultAK4->at(1));
         fillVariableWithValue( "neutrMult_j2", neMultAK4->at(1));
         fillVariableWithValue( "photonMult_j2", phoMultAK4->at(1));
	 fillVariableWithValue( "deltaPHIjj", DeltaPhiJJWide ) ;

	 //fillVariableWithValue( "Dijet_MassAK8", mjjAK8 ) ;  
	 //fillVariableWithValue( "Dijet_MassC", mjjCA8 ) ;
	 // if(wdijet.M()<1){
	 //    std::cout << " INV MASS IS " << wdijet.M() << std::endl;
	 //    std::cout << " Delta Eta IS " << DeltaEtaJJWide << " n is  " << widejets.size() << std::endl;
	 //    std::cout << " INV MASS FROM NTUPLE AK8 " << mjjAK8 << std::endl;
	 //    //std::cout << " INV MASS FROM NTUPLE CA8 " << mjjCA8 << std::endl;
       }

     //find intime BX
     int idx_InTimeBX=-1;
     for(size_t j=0; j<PileupOriginBX->size(); ++j)
       {
	 //cout << PileupOriginBX->at(j) << endl;	 
	 if(PileupOriginBX->at(j)==0)
	   {
	     idx_InTimeBX = j;
	     //cout << "idx_InTimeBX: " << idx_InTimeBX << endl; 
	   }
       }    

     //no cuts on these variables, just to store in output
     if(idx_InTimeBX > -1 )
       fillVariableWithValue("trueVtx",PileupInteractions->at(idx_InTimeBX));
     else
       fillVariableWithValue("trueVtx",999);
     fillVariableWithValue("MET",met);
     double METoverHTAK4=double(met/htAK4);
     fillVariableWithValue("METoverHTAK4",METoverHTAK4);
     fillVariableWithValue("HTAK4",htAK4);
     fillVariableWithValue("ptHat",ptHat);

     // Evaluate cuts (but do not apply them)
     evaluateCuts();
     
     // optional call to fill a skim with the full content of the input roottuple
     //if( passedCut("nJetFinal") ) fillSkimTree();
     
     if( passedCut("nVtx") 
	 && passedCut("IdTight_j1")
	 && passedCut("IdTight_j2")
	 && passedCut("nJet")
	 && passedCut("pTAK4_j1")
	 && passedCut("etaAK4_j1")
	 && passedCut("pTAK4_j2")
	 && passedCut("etaAK4_j2")
	 && passedCut("mjj") 
	 && passedCut("deltaETAjj") ){

       if(triggerResult->at(3)) h_mjj_HLTpass[0] -> Fill(MJJWide); //PFHT350
       if(triggerResult->at(0)) h_mjj_HLTpass[1] -> Fill(MJJWide); //PFHT900
       if(triggerResult->at(5)) h_mjj_HLTpass[2] -> Fill(MJJWide); //PFHT650MJJ900
       if(triggerResult->at(0) && triggerResult->at(5)) h_mjj_HLTpass[3] -> Fill(MJJWide); //PFHT900 && PFHT650MJJ900

       //std::cout << "triggerResult->at(3) = " << triggerResult->at(3) << "  triggerResult->at(0) = " << triggerResult->at(0) << "  triggerResult->at(5) = " << triggerResult->at(5) << std::endl;

     } 
     
     
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
	 //     CreateAndFillUserTH1D("h_htak4_mjjgt4000", 1000, 0, 10000, getVariableValue("HTAK4"));
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
   for (int i=0; i<4; i++){
     h_mjj_HLTpass[i]->Write();
   }
//   for(int ii=0; ii<3; ii++){
//     for(int i=0; i<nMassBins; i++)
//       {
//         double N1 = h_mjj_HLTpass[0]->GetBinContent(i+1);
//         double N2 = h_mjj_HLTpass[ii]->GetBinContent(i+1);
//	 double p  = 0;
//         double eU = 0;
//         double eL = 0;
//         double n,w;
//         if (N1 > 0)
//           { 
//             p = N2/N1;
//             n = N1+N2;
//             w = N2/n;
//             if (MODE==1) // Wilson for binomial
//               {
//        	 scale = 1.0; // makes sense only for the unprescaled trigger 
//        	 double d = sqrt(p*(1-p)/N1+0.25/(N1*N1));
//        	 eU = (p+0.5/N1+d)/(1+1/N1)-p;
//        	 eL = p-(p+0.5/N1-d)/(1+1/N1);
//               } 
//             else // Wilson for Poisson ratio
//               {
//        	 double d  = sqrt(w*(1-w)/n+0.25/(n*n));
//        	 double UB = (w+0.5/n+d)/(1+1/n);
//        	 double LB = (w+0.5/n-d)/(1+1/n);
//        	 eU = UB/(1-UB)-p;
//        	 eL = p-LB/(1-LB);     
//               }   
//           }
//         cout<<N1<<" "<<N2<<" "<<p<<" "<<eL<<" "<<eU<<endl;
//         vx[i][ii]   = h_mjj_HLTpass[0]->GetBinCenter(i+1);
//         vy[i][ii]   = p*scale;
//         vexl[i][ii] = h_mjj_HLTpass[0]->GetBinWidth(i+1)/2;
//         vexh[i][ii] = h_mjj_HLTpass[0]->GetBinWidth(i+1)/2;      
//         veyl[i][ii] = eL*scale;
//         veyh[i][ii] = eU*scale;     
//       }
//   }
//
//   TGraphAsymmErrors* g_eff_PFHT900 = new TGraphAsymmErrors(nMassBins,vx[0],vy[0],vexl[0],vexh[0],veyl[0],veyh[0]);  
//   TGraphAsymmErrors* g_eff_PFHT650MJJ900 = new TGraphAsymmErrors(nMassBins,vx[1],vy[1],vexl[1],vexh[1],veyl[1],veyh[1]);  
//   TGraphAsymmErrors* g_eff_PFHT900_AND_MJJ900 = new TGraphAsymmErrors(nMassBins,vx[2],vy[2],vexl[2],vexh[2],veyl[2],veyh[2]);  
//   
//   g_eff_PFHT900->SetName("g_eff_PFHT900");
//   g_eff_PFHT650MJJ900->SetName("g_eff_PFHT650MJJ900");
//   g_eff_PFHT900_AND_MJJ900->SetName("g_eff_PFHT900_AND_MJJ900");
//   g_eff_PFHT900->Write();
//   g_eff_PFHT650MJJ900->Write();
//   g_eff_PFHT900_AND_MJJ900->Write();
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
