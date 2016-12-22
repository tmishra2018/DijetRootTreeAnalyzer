#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TParameter.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>

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
    std::cout << "Reapplying JECs on the fly" << std::endl;

    std::string L1Path = "data/Spring16_25nsV8BCD_MC/Spring16_25nsV8BCD_MC_L1FastJet_AK4PFchs.txt";
    std::string L2Path = "data/Spring16_25nsV8BCD_MC/Spring16_25nsV8BCD_MC_L2Relative_AK4PFchs.txt";
    std::string L3Path = "data/Spring16_25nsV8BCD_MC/Spring16_25nsV8BCD_MC_L3Absolute_AK4PFchs.txt";
    std::string L1DATAPath = "data/Spring16_25nsV8BCD_DATA/Spring16_25nsV8BCD_DATA_L1FastJet_AK4PFchs.txt";
    std::string L2DATAPath = "data/Spring16_25nsV8BCD_DATA/Spring16_25nsV8BCD_DATA_L2Relative_AK4PFchs.txt"; 
    std::string L3DATAPath = "data/Spring16_25nsV8BCD_DATA/Spring16_25nsV8BCD_DATA_L3Absolute_AK4PFchs.txt";
    std::string L2L3ResidualPath = "data/Spring16_25nsV8BCD_DATA/Spring16_25nsV8BCD_DATA_L2Residual_AK4PFchs.txt"; //don't forget to put l2l3res for the analysis
    
    
    std::string L1RCcorrDATAPath = "data/Spring16_25nsV8BCD_DATA/Spring16_25nsV8BCD_DATA_L1RC_AK4PFchs.txt";
    std::string L1RCcorrMCPath = "data/Spring16_25nsV8BCD_MC/Spring16_25nsV8BCD_MC_L1RC_AK4PFchs.txt";

    //uncertainty
    unc = new JetCorrectionUncertainty("data/Spring16_25nsV8BCD_DATA/Spring16_25nsV8BCD_DATA_Uncertainty_AK4PFchs.txt");
        
    L1Par = new JetCorrectorParameters(L1Path);
    L2Par = new JetCorrectorParameters(L2Path);
    L3Par = new JetCorrectorParameters(L3Path);
    L1DATAPar = new JetCorrectorParameters(L1DATAPath);
    L2DATAPar = new JetCorrectorParameters(L2DATAPath);
    L3DATAPar = new JetCorrectorParameters(L3DATAPath);
    L2L3Residual = new JetCorrectorParameters(L2L3ResidualPath);
    
    L1JetParForTypeI =new JetCorrectorParameters(L1RCcorrDATAPath);
    L1JetParForTypeIMC =new JetCorrectorParameters(L1RCcorrMCPath);


    std::vector<JetCorrectorParameters> vPar;
    std::vector<JetCorrectorParameters> vPar_data;
    vPar.push_back(*L1Par);
    vPar.push_back(*L2Par);
    vPar.push_back(*L3Par);
    
     std::vector<JetCorrectorParameters> vParTypeI;
    vParTypeI.push_back(*L1JetParForTypeI);
    std::vector<JetCorrectorParameters> vParTypeIL123;
    vParTypeIL123.push_back(*L1DATAPar);
    vParTypeIL123.push_back(*L2DATAPar);
    vParTypeIL123.push_back(*L3DATAPar);
    vParTypeIL123.push_back(*L2L3Residual);
    
    std::vector<JetCorrectorParameters> vParTypeIMC;
    vParTypeIMC.push_back(*L1JetParForTypeIMC);
    std::vector<JetCorrectorParameters> vParTypeIL123MC;
    vParTypeIL123MC.push_back(*L1Par);
    vParTypeIL123MC.push_back(*L2Par);
    vParTypeIL123MC.push_back(*L3Par);
    
   
    //residuals are applied only to data
    vPar_data.push_back(*L1DATAPar);
    vPar_data.push_back(*L2DATAPar);
    vPar_data.push_back(*L3DATAPar);
    vPar_data.push_back(*L2L3Residual);

    JetCorrector = new FactorizedJetCorrector(vPar); assert(JetCorrector);
    JetCorrector_data = new FactorizedJetCorrector(vPar_data); assert(JetCorrector_data);
    
    JetCorrectortypI      = new FactorizedJetCorrector(vParTypeI); assert(JetCorrectortypI);
    JetCorrectortypIL123  = new FactorizedJetCorrector(vParTypeIL123); assert(JetCorrectortypIL123);
    JetCorrectortypIMC      = new FactorizedJetCorrector(vParTypeIMC); assert(JetCorrectortypIMC);
    JetCorrectortypIL123MC  = new FactorizedJetCorrector(vParTypeIL123MC); assert(JetCorrectortypIL123MC);
    
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
   
   
    Double_t genweight = 0.;
    Float_t trueInteractionall =0;
    
    TTree *PUvariable = new TTree("puvariable","puvariable");
    PUvariable->Branch("Generatorweight", &genweight, "genweight/D");
    PUvariable->Branch("TrueInteractionall", &trueInteractionall, "trueInteractionall/F");    
     
     
    TH2F *DeltaPhiAlpha = new TH2F("DeltaPhi_vs_alpha", "DeltaPhi_vs_alpha", 50, 0, 1. ,50, 0, 4  ) ;
      
    TH1F *SumWeight  = new TH1F("h_sumW", "h_sumW", 1, -0.5, 5.5) ;    
    SumWeight->Sumw2();
    TParameter<float> *totalluminosityP = new  TParameter<float>("totallumi", 0.);
    float storelumi = -15;
    float totalluminosity = 0.; 
   int nselectedevent = 0;
   /////////initialize variables

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;      
   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
  Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<2000;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     
     nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     // if (Cut(ientry) < 0) continue;

     ////////////////////// User's code starts here ///////////////////////

     ///Stuff to be done for every event
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
       bool isData = 0;
     if(idx_InTimeBX > -1 ) isData = 0;
     else isData = 1;
     
     
     if(!isData)
     {
         SumWeight->Fill(0.,weight);
         genweight = weight;
         trueInteractionall = npu->at(1)/*PileupInteractions->at(idx_InTimeBX)*/;
         PUvariable->Fill();
         
     }else
     {
          SumWeight->Fill(0.,1);
          weight = 1 ;  
          trueInteractionall = 999;   
          PUvariable->Fill();
     }
    
    
     
     
     
     if(nPhotonTight==1 && nPhoton == 1  && !HaspixelSeed->at(0) && PhotonLoosePt->at(0)>40.  && nMuonsLoose == -999) 
{
     size_t no_jets_ak4=jetPtAK4->size();
    
     
     TLorentzVector gamma1      ;
     TLorentzVector gamma1smear ;
       gamma1.SetPtEtaPhiE(PhotonLoosePt->at(0),PhotonLooseEta->at(0),PhotonLoosePhi->at(0),PhotonLooseEnergy->at(0));
       gamma1smear.SetPtEtaPhiE(PhotonsmearPt->at(0),PhotonsmearEta->at(0),PhotonsmearPhi->at(0),PhotonsmearEnergy->at(0)); 
       
       
       
     size_t neletron =  electronPt->size();
     bool keepEvent = true;
    for (size_t j = 0; j < neletron; j++) {
      double deltaR = std::hypot((electronEta->at(j)-gamma1.Eta()),(electronPhi->at(j)-gamma1.Phi()));
      if (deltaR < 0.13) {
        keepEvent = false;
        break;
      }
    }
    
     if(!keepEvent) continue;
     vector<TLorentzVector> AK4jets;
     TLorentzVector ak4j1, ak4j2;
     
     TLorentzVector rawMet;
     rawMet.SetPtEtaPhiE(metPt,metEta,metPhi,metEnergy);
     TLorentzVector MetTypeI;      

     resetCuts();

     //find intime BX
    

     std::vector<double> jecFactors;
     std::vector<double> jecUncertainty;
     // new JECs could change the jet pT ordering. the vector below
     // holds sorted jet indices after the new JECs had been applied
     std::vector<unsigned> sortedJetIdx;
     
     
     int nfakejet = 0;
     int fakejetIdx = -999;
     double fakecorrection = 0. ;

     if( int(getPreCutValue1("useJECs"))==1 )
       {
	 // sort jets by increasing pT

	 std::multimap<double, unsigned> sortedJets;
	 for(size_t j=0; j<no_jets_ak4; ++j)
	   {
	      
	     JetCorrector->setJetEta(jetEtaAK4->at(j));
	     JetCorrector->setJetPt(jetPtAK4->at(j)/jetJecAK4->at(j)); //pTraw
	     JetCorrector->setJetA(jetAreaAK4->at(j));
	     JetCorrector->setRho(rho);

  	     JetCorrector_data->setJetEta(jetEtaAK4->at(j));
	     JetCorrector_data->setJetPt(jetPtAK4->at(j)/jetJecAK4->at(j)); //pTraw
	     JetCorrector_data->setJetA(jetAreaAK4->at(j));
	     JetCorrector_data->setRho(rho);


  	     //nominal value of JECs
	     double correction;//, old_correction, nominal_correction;
	     
	     if (isData == 1) correction = JetCorrector_data->getCorrection();
	     else correction = JetCorrector->getCorrection();
	     
	     //JEC uncertainties
	     unc->setJetEta(jetEtaAK4->at(j));
	     unc->setJetPt(jetPtAK4->at(j)/jetJecAK4->at(j)*correction);
	     double uncertainty = unc->getUncertainty(true);
	     jecUncertainty.push_back(uncertainty); 

	     
	     //use "shifted" JECs for study of systematic uncertainties 
	     if( int(getPreCutValue1("shiftJECs"))==1 ){
	       
	       
	       //shift of the corresponding unc
	       correction = correction + getPreCutValue2("shiftJECs")*uncertainty*correction;
	       
	       
	   }

	 jecFactors.push_back(correction);
	 sortedJets.insert(std::make_pair((jetPtAK4->at(j)/jetJecAK4->at(j))*correction, j));
	 if(std::hypot((jetEtaAK4->at(j)-gamma1.Eta()),(jetPhiAK4->at(j)-gamma1.Phi()))<0.4)
	      { 
	         nfakejet++ ;
	         fakejetIdx = j;
	         fakecorrection = correction;
	        
	      }

       }
       int size_all = sortedJets.size();
       if(nfakejet !=0)
       {
          sortedJets.erase((jetPtAK4->at(fakejetIdx)/jetJecAK4->at(fakejetIdx))*fakecorrection);

        }
        int size_true = sortedJets.size();
        
       if(size_all - size_true > 1) std::cout<<"Warning : more than two jet removed because of photon misidentification"<<std::endl;
     // get jet indices in decreasing pT order
     for(std::multimap<double, unsigned>::const_reverse_iterator it = sortedJets.rbegin(); it != sortedJets.rend(); ++it)
	 sortedJetIdx.push_back(it->second);

     }
     else if( int(getPreCutValue1("noJECs"))==1  )
       {
	 // sort jets by increasing pT
	 std::multimap<double, unsigned> sortedJets;
	 for(size_t j=0; j<no_jets_ak4; ++j) //same ordering of original root trees
	   {
	     jecUncertainty.push_back(0.); 
	     jecFactors.push_back(1.);
	     sortedJets.insert(std::make_pair((jetPtAK4->at(j)/jetJecAK4->at(j)), j)); //raw
	   }       
	 // get jet indices in decreasing pT order
	 for(std::multimap<double, unsigned>::const_reverse_iterator it = sortedJets.rbegin(); it != sortedJets.rend(); ++it)
	   sortedJetIdx.push_back(it->second);
       }
     else
       {
	 for(size_t j=0; j<no_jets_ak4; ++j) //same ordering of original root trees
	   {
	     jecFactors.push_back(jetJecAK4->at(j));
	     jecUncertainty.push_back(0.); 
	     sortedJetIdx.push_back(j);
	   }
       }


     //#############################################################
     //########## NOTE: from now on sortedJetIdx[ijet] should be used
     //#############################################################

    
     //count ak4 jets passing pt threshold and id criteria
     int Nak4 = 0;
     double HTak4 = 0;
     for(size_t ijet=0; ijet<no_jets_ak4-nfakejet; ++ijet)
       {	 
	  
	 if(fabs(jetEtaAK4->at(sortedJetIdx[ijet])) < getPreCutValue1("jetFidRegion")
	    && idLAK4->at(sortedJetIdx[ijet]) == getPreCutValue1("tightJetID")
	    && (jecFactors[sortedJetIdx[ijet]]/jetJecAK4->at(sortedJetIdx[ijet]))*jetPtAK4->at(sortedJetIdx[ijet]) > getPreCutValue1("ptCut"))
	   {
	     Nak4 += 1;
	     HTak4 += (jecFactors[sortedJetIdx[ijet]]/jetJecAK4->at(sortedJetIdx[ijet]))*jetPtAK4->at(sortedJetIdx[ijet]);
	   }
       }




     //AK4 jets
    if(no_jets_ak4-nfakejet>=1) 
      { 

     
	 if(fabs(jetEtaAK4->at(sortedJetIdx[0])) < getPreCutValue1("jetFidRegion") 
	    && (jecFactors[sortedJetIdx[0]]/jetJecAK4->at(sortedJetIdx[0]))*jetPtAK4->at(sortedJetIdx[0]) > getPreCutValue1("pt0Cut") && idLAK4->at(sortedJetIdx[0]) == getPreCutValue1("tightJetID"))
	   {
	     
		 //cout << "filling ak4j1 and ak4j2" << endl;
		 //cout << "pt ak4 j1 = " << (jecFactors[sortedJetIdx[0]]/jetJecAK4->at(sortedJetIdx[0])) *jetPtAK4->at(sortedJetIdx[0]) << endl;
		 ak4j1.SetPtEtaPhiM( (jecFactors[sortedJetIdx[0]]/jetJecAK4->at(sortedJetIdx[0])) *jetPtAK4->at(sortedJetIdx[0])
				     ,jetEtaAK4->at(sortedJetIdx[0])
				     ,jetPhiAK4->at(sortedJetIdx[0])
				     , (jecFactors[sortedJetIdx[0]]/jetJecAK4->at(sortedJetIdx[0])) *jetMassAK4->at(sortedJetIdx[0]));
		if(no_jets_ak4-nfakejet >=2 && fabs(jetEtaAK4->at(sortedJetIdx[1])) < getPreCutValue1("jetFidRegion") 
		&& (jecFactors[sortedJetIdx[1]]/jetJecAK4->at(sortedJetIdx[1]))*jetPtAK4->at(sortedJetIdx[1]) > getPreCutValue1("pt1Cut") && idLAK4->at(sortedJetIdx[1]) == getPreCutValue1("tightJetID"))
	       {
		      ak4j2.SetPtEtaPhiM( (jecFactors[sortedJetIdx[1]]/jetJecAK4->at(sortedJetIdx[1])) *jetPtAK4->at(sortedJetIdx[1])
				     ,jetEtaAK4->at(sortedJetIdx[1])
				     ,jetPhiAK4->at(sortedJetIdx[1])
				     , (jecFactors[sortedJetIdx[1]]/jetJecAK4->at(sortedJetIdx[1])) *jetMassAK4->at(sortedJetIdx[1]));
	       }
	   }
   
     }


 //----------------------TYPE I MET computation------------------
 // std::cout<<"rawMet: "<< rawMet.Pt() << " "<< rawMet.Eta() << " " << rawMet.Phi() << " " << rawMet.Et() <<std::endl; 
  
               

  
  double deltaPx = 0., deltaPy = 0.;
  TLorentzVector  jetRC,corrJet; 
  
  
  for (Long64_t it=0; it< no_jets_ak4-nfakejet ; it++) {
    if(fakejetIdx==it)continue;
    
    
    double corrs = 1.;
    double corrsForTypeI = 1.;
    
    
      if(isData)
      {
      JetCorrectortypI->setJetEta(jetEtaAK4RC->at(it));
      JetCorrectortypI->setJetPt(jetPtAK4RC->at(it));
      JetCorrectortypI->setJetA(jetAreaAK4->at(it));
      JetCorrectortypI->setRho(rho);
      corrsForTypeI = JetCorrectortypI->getCorrection(); //only RC

      JetCorrectortypIL123 ->setJetEta(jetEtaAK4RC->at(it));
      JetCorrectortypIL123 ->setJetPt(jetPtAK4RC->at(it));
      JetCorrectortypIL123 ->setJetA(jetAreaAK4->at(it));
      JetCorrectortypIL123 ->setRho(rho);
      corrs = JetCorrectortypIL123->getCorrection(); // L1L2L3
      
      }else{

      JetCorrectortypIMC->setJetEta(jetEtaAK4RC->at(it));
      JetCorrectortypIMC->setJetPt(jetPtAK4RC->at(it));
      JetCorrectortypIMC->setJetA(jetAreaAK4->at(it));
      JetCorrectortypIMC->setRho(rho);
      corrsForTypeI = JetCorrectortypIMC->getCorrection(); //only RC
      JetCorrectortypIL123MC ->setJetEta(jetEtaAK4RC->at(it));
      JetCorrectortypIL123MC ->setJetPt(jetPtAK4RC->at(it));
      JetCorrectortypIL123MC ->setJetA(jetAreaAK4->at(it));
      JetCorrectortypIL123MC ->setRho(rho);
      corrs = JetCorrectortypIL123MC->getCorrection(); // L1L2L3

      }

    jetRC.SetPtEtaPhiE((jetPtAK4RC->at(it))*corrsForTypeI,jetEtaAK4RC->at(it),jetPhiAK4RC->at(it),(jetEnergyAK4RC->at(it))*corrsForTypeI);
     // only RC
      
    corrJet.SetPtEtaPhiE((jetPtAK4RC->at(it)*corrs),jetEtaAK4RC->at(it),jetPhiAK4RC->at(it),(jetEnergyAK4RC->at(it))*corrs);;
      
    double dR = gamma1.DeltaR(corrJet);
   //   std::cout<<"pt of ak4jets "<<corrJet.Pt()<<" dr "<<dR<<std::endl;
    if(corrJet.Pt() > 15 && dR > 0.25) {
	
      double emEnergyFraction = jetNemfAK4->at(it) + jetCemfAK4->at(it);
      if (emEnergyFraction > 0.90)
	continue;
	
      deltaPx += (corrJet.Px() - jetRC.Px());
      deltaPy += (corrJet.Py() - jetRC.Py());
    } // jet.pt() && dR
  }//loop over jets
    
  double correctedMetPx = rawMet.Px() - deltaPx;
  double correctedMetPy = rawMet.Py() - deltaPy;
  double correctedMetPt = sqrt(correctedMetPx * correctedMetPx + correctedMetPy * correctedMetPy);
  
  MetTypeI.SetPxPyPzE(correctedMetPx,correctedMetPy, rawMet.Pz(), std::hypot(correctedMetPx,correctedMetPy));  

  // std::coulst<<"MET type I: "<< MetTypeI.Pt() << " "<< MetTypeI.Eta() << " " << MetTypeI.Phi() << " " << MetTypeI.Et() <<std::endl; 

//------------------------END met calculation------------------


  double Rbal = -999. ;
  double Rmpf = -999. ;
  double Rmpfraw = -999. ;

  Rbal = (ak4j1.Pt()/gamma1.Pt());
  double deltPHIgj = gamma1.DeltaPhi(ak4j1);
  double deltaphiPhomet = (gamma1.Phi()- metPhi);

  Rmpf = 1. + (MetTypeI.Px()*gamma1.Px()+MetTypeI.Py()*gamma1.Py())/std::pow(gamma1.Pt(),2);
  Rmpfraw = 1. + metPt*gamma1.Pt()*cos(deltaphiPhomet)/std::pow(gamma1.Pt(),2);
  double alpha = -999. ;
  
  alpha = (ak4j2.Pt()/gamma1.Pt());

 DeltaPhiAlpha->Fill(alpha,deltPHIgj);

     if(deltPHIgj>=2.8 && alpha<0.3)
     {
     if( ak4j1.Pt()>0   )
     {
       
       AK4jets.push_back( ak4j1 );
       if(ak4j2.Pt()>0)
       {
          AK4jets.push_back( ak4j2 );
       }
     }
     
     
      //== Fill Variables ==
     if(ak4j1.Pt()>0 && idLAK4->at(sortedJetIdx[0] == 1)) 
    {
     nselectedevent++;
     
     if(lumi != storelumi)
      {
        totalluminosity += lumi ;
      }
      storelumi = lumi;

     
     
     
     fillVariableWithValue("N_photon", nPhoton);
     fillVariableWithValue("run",runNo);     
     fillVariableWithValue("event",evtNo);     
     fillVariableWithValue("lumi",lumi);     
     fillVariableWithValue("nVtx",nvtx);     
     fillVariableWithValue("nJet",no_jets_ak4-nfakejet);
     fillVariableWithValue("metSig",metSig);
     fillVariableWithValue("Nak4",Nak4);
     fillVariableWithValue( "PassJSON", passJSON (runNo, lumi, isData));
     fillVariableWithValue("rho",rho);
     
     fillVariableWithValue("Pt_photon", gamma1.Pt());
     fillVariableWithValue("Eta_photon", gamma1.Eta());
     fillVariableWithValue("Phi_photon", gamma1.Phi());
     fillVariableWithValue("Energy_photon", gamma1.Energy());
     fillVariableWithValue("Rbalancing", Rbal);
     fillVariableWithValue("RMPF", Rmpf);
     fillVariableWithValue("RMPFRAW", Rmpfraw);
     fillVariableWithValue("alpha", alpha);
     fillVariableWithValue("deltaPHIgj", deltPHIgj);
     fillVariableWithValue("sigmaietaieta_photon", Photonfull5x5SigmaIEtaIEtaMapToken->at(0));
     fillVariableWithValue("CHiso_photon", PhotonphoChargedIsolationToken->at(0));
     fillVariableWithValue("NHiso_photon", PhotonphoNeutralHadronIsolationToken->at(0));
     fillVariableWithValue("Photoniso_photon", PhotonphoPhotonIsolationToken->at(0));
     fillVariableWithValue("Pt_photonSC",PhotonSCPt->at(0));
     fillVariableWithValue("Eta_photonSC",PhotonSCEta->at(0));
     fillVariableWithValue("Phi_photonSC",PhotonSCPhi->at(0));
     fillVariableWithValue("Energy_photonSC",PhotonSCEnergy->at(0));
     fillVariableWithValue("hadTowOverEm",hadTowOverEm->at(0));
     if(isData)
       {
          fillVariableWithValue("weight", 1);
          
         //  SumWeight->Fill(0.,weight);
        }
       if(!isData)
       {
         fillVariableWithValue("Pt_photonGEN",photonPtGen->at(0));
	 fillVariableWithValue("Eta_photonGEN",photonEtaGen->at(0));
	 fillVariableWithValue("Phi_photonGEN",photonPhiGen->at(0));
	 fillVariableWithValue("Energy_photonGEN",photonEnergyGen->at(0));
	 
	 fillVariableWithValue("weight",weight);
	// SumWeight->Fill(0.,weight);
	 
       }//MC
     

     if( AK4jets.size() >=1 )
     {
       //cout << "AK4jets.size() " <<  AK4jets.size() << endl;
       //cout << "IdTight_j1 : " << idTAK4->at(sortedJetIdx[0]) << endl;
       fillVariableWithValue( "IdTight_j1",idLAK4->at(sortedJetIdx[0]));
       fillVariableWithValue( "pTAK4_j1", AK4jets[0].Pt());
       fillVariableWithValue( "etaAK4_j1", AK4jets[0].Eta());
       fillVariableWithValue( "phiAK4_j1", AK4jets[0].Phi());
       
       //fillVariableWithValue( "jetPtAK4matchCaloJet_j1", jetPtAK4matchCaloJet->at(sortedJetIdx[0]));
       
       fillVariableWithValue( "jetJecAK4_j1", jecFactors[sortedJetIdx[0]] );
       fillVariableWithValue( "jetJecUncAK4_j1", jecUncertainty[sortedJetIdx[0]] );
       fillVariableWithValue( "jetCSVAK4_j1", jetCSVAK4->at(sortedJetIdx[0]) );
       //jetID
       fillVariableWithValue( "neutrHadEnFrac_j1", jetNhfAK4->at(sortedJetIdx[0]));
       fillVariableWithValue( "chargedHadEnFrac_j1", jetChfAK4->at(sortedJetIdx[0]));
       fillVariableWithValue( "photonEnFrac_j1", jetPhfAK4->at(sortedJetIdx[0]));
       fillVariableWithValue( "eleEnFract_j1", jetElfAK4->at(sortedJetIdx[0]));
       fillVariableWithValue( "muEnFract_j1", jetMufAK4->at(sortedJetIdx[0]));
       fillVariableWithValue( "neutrElectromFrac_j1", jetNemfAK4->at(sortedJetIdx[0]));
       fillVariableWithValue( "chargedElectromFrac_j1", jetCemfAK4->at(sortedJetIdx[0]));
       fillVariableWithValue( "chargedMult_j1", chMultAK4->at(sortedJetIdx[0]));
       fillVariableWithValue( "neutrMult_j1", neMultAK4->at(sortedJetIdx[0]));
       fillVariableWithValue( "photonMult_j1", phoMultAK4->at(sortedJetIdx[0]));
       if(!isData )
       {
         fillVariableWithValue( "pTAK4_j1GEN", jetPtGenAK4->at(sortedJetIdx[0]));
	 fillVariableWithValue( "etaAK4_j1GEN", jetEtaGenAK4->at(sortedJetIdx[0]));	       
	 fillVariableWithValue("PDGIDAK4_j1",jetpdgIDGenAK4->at(sortedJetIdx[0]));
       }
     }
     
     if( AK4jets.size() >=2){

       fillVariableWithValue( "IdTight_j2",idLAK4->at(sortedJetIdx[1]));
       fillVariableWithValue( "pTAK4_j2", AK4jets[1].Pt() );
       fillVariableWithValue( "etaAK4_j2", AK4jets[1].Eta());
       fillVariableWithValue( "phiAK4_j2", AK4jets[1].Phi());
       //fillVariableWithValue( "jetPtAK4matchCaloJet_j2", jetPtAK4matchCaloJet->at(sortedJetIdx[1]));
       fillVariableWithValue( "jetJecAK4_j2", jecFactors[sortedJetIdx[1]]); 
       fillVariableWithValue( "jetJecUncAK4_j2", jecUncertainty[sortedJetIdx[1]] );
       fillVariableWithValue( "jetCSVAK4_j2", jetCSVAK4->at(sortedJetIdx[1]) );
       //jetID
       fillVariableWithValue( "neutrHadEnFrac_j2", jetNhfAK4->at(sortedJetIdx[1]));
       fillVariableWithValue( "chargedHadEnFrac_j2", jetChfAK4->at(sortedJetIdx[1]));
       fillVariableWithValue( "photonEnFrac_j2", jetPhfAK4->at(sortedJetIdx[1]));
       fillVariableWithValue( "eleEnFract_j2", jetElfAK4->at(sortedJetIdx[1]));
       fillVariableWithValue( "muEnFract_j2", jetMufAK4->at(sortedJetIdx[1]));
       fillVariableWithValue( "neutrElectromFrac_j2", jetNemfAK4->at(sortedJetIdx[1]));
       fillVariableWithValue( "chargedElectromFrac_j2", jetCemfAK4->at(sortedJetIdx[1]));
       fillVariableWithValue( "chargedMult_j2", chMultAK4->at(sortedJetIdx[1]));
       fillVariableWithValue( "neutrMult_j2", neMultAK4->at(sortedJetIdx[1]));
       fillVariableWithValue( "photonMult_j2", phoMultAK4->at(sortedJetIdx[1]));
       //dijet
       fillVariableWithValue( "CosThetaStarAK4", TMath::TanH( (AK4jets[0].Eta()-AK4jets[1].Eta())/2 ));
       if(!isData )
       {
          

         fillVariableWithValue( "pTAK4_j2GEN", jetPtGenAK4->at(sortedJetIdx[1]));
	 fillVariableWithValue( "etaAK4_j2GEN", jetEtaGenAK4->at(sortedJetIdx[1]));
 	 fillVariableWithValue("PDGIDAK4_j2",jetpdgIDGenAK4->at(sortedJetIdx[1]));
	 

       } 
     }
     //no cuts on these variables, just to store in output
     if(!isData)
       fillVariableWithValue("trueInteraction",npu->at(1)/*PileupInteractions->at(idx_InTimeBX)*/);
     else if(isData)
       fillVariableWithValue("trueInteraction",999);     

     fillVariableWithValue("MET",MetTypeI.Et());
     fillVariableWithValue("MET_Pt",MetTypeI.Pt());
     fillVariableWithValue("MET_Eta",MetTypeI.Eta());
     fillVariableWithValue("MET_Phi",MetTypeI.Phi());
     fillVariableWithValue("METRAW",metEnergy);
     if(!isData && metEnergyGen)
     {
       fillVariableWithValue("METGEN",metEnergyGen);
       fillVariableWithValue("METGEN_Pt",metPtGen);
       fillVariableWithValue("METGEN_Eta",metEtaGen);
       fillVariableWithValue("METGEN_Phi",metPhiGen);
     }
     
     //double METoverHTAK4=double(met/htAK4);
     double METoverHTAK4=double(MetTypeI.Et()/HTak4);
     fillVariableWithValue("METoverHTAK4",METoverHTAK4);
     //fillVariableWithValue("HTAK4",htAK4);
     fillVariableWithValue("HTAK4",HTak4);
     fillVariableWithValue("ptHat",ptHat);

     // Trigger
     int NtriggerBits = triggerResult->size();
     if( NtriggerBits > 0 && isData)
     {
       fillVariableWithValue("passHLT_Photon30",triggerResult->at(0));// 
     //  fillVariableWithValue("prescaletrigger_Photon30",triggerPrescale->at(0));// 

     }
    if( NtriggerBits > 1 && isData)
    {
       fillVariableWithValue("passHLT_Photon50",triggerResult->at(1));//

     }
     if( NtriggerBits > 2 && isData)
     {
       fillVariableWithValue("passHLT_Photon75",triggerResult->at(2));// 

     }
     if( NtriggerBits > 3 && isData)
     {
       fillVariableWithValue("passHLT_Photon90",triggerResult->at(3));//

     }
     if( NtriggerBits > 4 && isData)
     {
       fillVariableWithValue("passHLT_Photon120",triggerResult->at(4));//

     }
     if( NtriggerBits > 5 && isData)
     {
       fillVariableWithValue("passHLT_Photon165",triggerResult->at(5));//

     }
     
     

     // Evaluate cuts (but do not apply them)
     evaluateCuts();
     
     

     // optional call to fill a skim with a subset of the variables defined in the cutFile (use flag SAVE)
     if( passedCut("Pt_photon") ) 
       {
	 fillReducedSkimTree();

       }

}//end dphi cut
}//end of pt j1 cut
     ////////////////////// User's code ends here ///////////////////////

}//end photon condition
 
   } // End loop over events
   totalluminosityP->SetVal(totalluminosity);
  
   SumWeight->Write();
   totalluminosityP->Write();
   PUvariable->Write();
   DeltaPhiAlpha->Write();
   std::cout<<" nb of selected event " << nselectedevent<<std::endl;
   std::cout << "analysisClass::Loop() ends" <<std::endl; 
   
    delete  PUvariable;

   }
