#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TParameter.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <vector>
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "CondFormats/DataRecord/interface/JetResolutionRcd.h"
#include "CondFormats/DataRecord/interface/JetResolutionScaleFactorRcd.h"
#define RESET_COLOR "\033[m"
#define MAKE_RED "\033[31m"
#define MAKE_BLUE "\033[34m"

/*

bool isNewonOldValidJetTight(const float& Eta_ak4, const float& nhf, const float& neMult, const float& chMult, const float& nemf, const float& muF, const float& chf, const float& cemF, const bool& isoldvalid, const bool& isDATA, const bool& hasgenjet, const bool& isFirst)
{
   
    int idL = -999 ;
    if(isoldvalid) {
    return true ;  
    }else{return false;}   
  //  if(!isDATA && !hasgenjet) return false;
    
    if(fabs(Eta_ak4) > 2.7  && fabs(Eta_ak4) <= 3.0)
    {
       idL = ( nemf>0.01 && nhf<0.98 && neMult>2)  ;
       
       
    }else{
       if(isoldvalid){
       idL =1;
       }else{idL =0;}
       
    } 
    
 return idL;

    
    
}*/

/*bool IsTight_photon(const float& hadronicOverEm, const float& full5x5SigmaIEtaIEta, const float& Photon_pt, ){

    bool isValid = true;
   // #1: H/E
    isValid &= hadronicOverEm < 0.0269;
    if (! isValid)
    return false;
    //#2: sigma ietaieta
    isValid &= full5x5SigmaIEtaIEta < 0.00994; //Official    
    if (! isValid)
    return false;

    
    isValid &= getCorrectedPFIsolation((*phoChargedIsolationMap)[photonRef], rho, photonRef->eta(), IsolationType::CHARGED_HADRONS) < 0.202;
    isValid &= getCorrectedPFIsolation((*phoNeutralHadronIsolationMap)[photonRef], rho, photonRef->eta(), IsolationType::NEUTRAL_HADRONS) < (0.264 + 0.0148*Photon_pt+0.000017*(Photon_pt*Photon_pt ) );
    isValid &= getCorrectedPFIsolation((*phoPhotonIsolationMap)[photonRef], rho, photonRef->eta(), IsolationType::PHOTONS) < (2.362+0.0047*Photon_pt);
    isValid &= photonRef->passElectronVeto();
    if (! isValid)
    return false;
    isValid &= photonRef->r9() >0.90;
    if (! isValid)
    return false;
return isValid;


}*/

// Jet selection Cut ID from this twiki : https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
bool isNewonOldValidJetTight(const float& Eta_ak4, const float& nhf, const float& neMult, const float& chMult, const float& nemf, const float& muF, const float& chf, const float& cemF, const bool& isoldvalid, const bool& isDATA, const bool& hasgenjet, const bool& isFirst)
{
   
    int idL = -999 ;
    
  
  if(fabs(Eta_ak4) > 3.0){
    idL = ( neMult > 10 && nhf > 0.02 && nemf < 0.9 );
  
  }
    
    if(fabs(Eta_ak4) > 2.7  && fabs(Eta_ak4) <= 3.0)
    {
       idL = ( nemf < 0.99 && neMult>2)  ;
      
      
      }
      
      if(fabs(Eta_ak4) >= 2.4  && fabs(Eta_ak4) <= 2.7){
        
        idL = (nhf < 0.9 && nemf < 0.9  && neMult >1 && muF < 0.8);
        
        
        
      }
      
       if( fabs(Eta_ak4) < 2.4 ){
  
  
          idL = (nhf < 0.9 && nemf < 0.9  && neMult >1 && muF < 0.8 && chMult > 0 && chf > 0 && cemF < 0.8) ;
      
      }
     
    
 return idL;

    
    
}




analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile)
{
  std::cout << "analysisClass::analysisClass(): begins " << std::endl;
    
  // For JECs
  if( int(getPreCutValue1("useJECs"))==1 )
  {
    std::cout << "Reapplying JECs on the fly" << std::endl;
    // getting the JEC textfile from config file
    std::string L1Path =  getPreCutString1("L1MC");
    std::string L2Path = getPreCutString1("L2MC");
    std::string L3Path = getPreCutString1("L3MC");
    std::string L1DATAPath = getPreCutString1("L1Dat");
    std::string L2DATAPath = getPreCutString1("L2Dat");
    std::string L3DATAPath = getPreCutString1("L3Dat");
    std::string L2L3ResidualPath = getPreCutString1("L2L3Dat");

    
    std::string L1RCcorrDATAPath = getPreCutString1("L1RCDat");
    std::string L1RCcorrMCPath = getPreCutString1("L1RCMC");
    
    std::string JER_resopath = getPreCutString1("JERreso");
    std::string JER_SFpath = getPreCutString1("JERSF");
    //uncertainty
    unc = new JetCorrectionUncertainty(getPreCutString1("JECuncDat"));
    unc_TI = new JetCorrectionUncertainty(getPreCutString1("JECuncDat"));
    unc_RC = new JetCorrectionUncertainty(getPreCutString1("JECuncDat"));
    unc_MC = new JetCorrectionUncertainty(getPreCutString1("JECuncMC"));
    unc_MC_RC = new JetCorrectionUncertainty(getPreCutString1("JECuncMC"));
    unc_MC_TI = new JetCorrectionUncertainty(getPreCutString1("JECuncMC"));    
    L1Par = new JetCorrectorParameters(L1Path);

    L2Par = new JetCorrectorParameters(L2Path);

    L3Par = new JetCorrectorParameters(L3Path);

    L1DATAPar = new JetCorrectorParameters(L1DATAPath);

    L2DATAPar = new JetCorrectorParameters(L2DATAPath);

    L3DATAPar = new JetCorrectorParameters(L3DATAPath);

    L2L3Residual = new JetCorrectorParameters(L2L3ResidualPath);

    
    L1JetParForTypeI =new JetCorrectorParameters(L1RCcorrDATAPath);

    L1JetParForTypeIMC =new JetCorrectorParameters(L1RCcorrMCPath);
    
    resolution = new JME::JetResolution(JER_resopath);
    resolution_sf = new JME::JetResolutionScaleFactor(JER_SFpath);
    
    resolution_TI = new JME::JetResolution(JER_resopath);
    resolution_sf_TI = new JME::JetResolutionScaleFactor(JER_SFpath);

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
    if(int(getPreCutValue1("useResidual"))==1){ 
    std::cout<<" apply residual "<<std::endl;
    vParTypeIL123.push_back(*L2L3Residual);}
    
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
    if(int(getPreCutValue1("useResidual"))==1){ 
    std::cout<<" apply residual "<<std::endl;
    vPar_data.push_back(*L2L3Residual);}

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
    
   
    
    
    
    
     
    TH2F *DeltaPhiAlpha = new TH2F("DeltaPhi_vs_alpha", "DeltaPhi_vs_alpha", 50, 0, 1. ,50, 0, 5  ) ;
    TH2F *EtaPhiJet_beforehot = new TH2F("EtaPhioccupency_before", "EtaPhioccupency_before", 100, -5.2, 5.2 ,100, -3.5, 3.5  ) ;
    TH2F *EtaPhiJet_afterhot = new TH2F("EtaPhioccupency_after", "EtaPhioccupency_after", 100, -5.2, 5.2 ,100, -3.5, 3.5  ) ;
      
    TH1F *SumWeight  = new TH1F("h_sumW", "h_sumW", 1, -0.5, 5.5) ;    
    SumWeight->Sumw2();
    TParameter<float> *totalluminosityP = new  TParameter<float>("totallumi", 0.);
    float storelumi = -15;
    float totalluminosity = 0.; 
    int nselectedevent = 0;
    int ncut_photon = 0;
    int ncut_photonpt = 0;
    int ncut_nophoton = 0;
    int ncut_deltaphi = 0;
    int ncut_alpha = 0;
    int ncut_pixelseed = 0;
    int ncut_muons = 0;
    int ncut_electron = 0;
    int ncut_jet = 0 ;
    int ncut_ptjet = 0;
    int Vtxcut = 0 ;
    int testcount = 0 ;
    int Is_PU = 0 ;
    int is_hot_area = 0 ;
   /////////initialize variables

   Long64_t nentries = fChain->GetEntriesFast();//10000;//fChain->GetEntriesFast();//10000;//10000;//1000000; //
   
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;      
   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
  Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     if(jentry > 0){            
             pt_jets_    ->clear();
            phi_jets_   ->clear();
            eta_jets_   ->clear();
            mass_jets_  ->clear();
            emF_jets_   ->clear();
            IsID_jets_  ->clear();
            }
   //for (Long64_t jentry=0; jentry<2000;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     
     nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry < 10 || jentry%100000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
    
   
     ////////////////////// User's code starts here ///////////////////////

     ///Stuff to be done for every event
     int idx_InTimeBX=-1;
     for(size_t j=0; j<PileupOriginBX->size(); ++j)
       {

	 if(PileupOriginBX->at(j)==0)
	   {
	     idx_InTimeBX = j;

	   }
       }
       bool isData = 0;
     if(idx_InTimeBX > -1 ) isData = 0;
     else isData = 1;
     
     //start selection workflow 
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
          trueInteractionall = -999;   
          PUvariable->Fill();
     }

    if(!goodPVtx){
       
       Vtxcut++;
       continue;
    } 
     if(nPhoton > 0){
     
   

          if( nPhotonTight == 1){

          testcount ++; 
             size_t nb_photon =  PhotonLoosePt->size();
             int indexgoodpho = 0  ; 
             for(size_t pho = 0 ; pho <  nb_photon ; pho ++)
          {
              if(isPhotonTight->at(pho)){
               indexgoodpho = pho ;
               break ; 
              }
          }
          
               if(!HaspixelSeed->at(indexgoodpho)){
                     if(PhotonLoosePt->at(indexgoodpho)>= 40.){
                        bool keepmuon = true ;
                        if(nMuonsLoose != 0){ 
                        size_t nb_muons = muonPt->size();
                       
                        for(size_t mu = 0 ; mu < nb_muons ; mu++)
                        {
                          if(muonPt->at(mu) > 5.){
                           keepmuon = false ;
                           break;
                          }
                        }
                        }
                        
                          if(nMuonsLoose == 0  || keepmuon){


     size_t no_jets_ak4=jetPtAK4->size();


     TLorentzVector gamma1     ;
     TLorentzVector gammaloose ;
     gammaloose.SetPtEtaPhiE(PhotonLoosePt->at(indexgoodpho),PhotonLooseEta->at(indexgoodpho),PhotonLoosePhi->at(indexgoodpho),PhotonLooseEnergy->at(indexgoodpho));
     gamma1.SetPtEtaPhiE(PhotonLoosePt->at(indexgoodpho),PhotonLooseEta->at(indexgoodpho),PhotonLoosePhi->at(indexgoodpho),PhotonLooseEnergy->at(indexgoodpho));
       
     size_t neletron =  electronPt->size();
     bool keepEvent = true;
    for (size_t j = 0; j < neletron; j++) {
      double deltaR = std::hypot((electronEta->at(j)-gamma1.Eta()),(electronPhi->at(j)-gamma1.Phi()));
      if (deltaR < 0.13) {
        keepEvent = false;
        break;
      }
    }
    
     if(!keepEvent){
      ncut_electron++;
      continue;
      
      }
     vector<TLorentzVector> AK4jets;
     TLorentzVector ak4j1, ak4j2;
     
     TLorentzVector rawMet;
     rawMet.SetPtEtaPhiE(metPt,metEta,metPhi,metEnergy);
     TLorentzVector MetTypeI;      

     resetCuts();
     std::vector<double> jecFactors;
     std::vector<double> jecUncertainty;
     std::vector<double> jerFactors;
     std::vector<double> CorrFactors;
     // new JECs could change the jet pT ordering. the vector below
     // holds sorted jet indices after the new JECs had been applied
     std::vector<unsigned> sortedJetIdx;
     
     
     int nfakejet = 0;
     int fakejetIdx = -999;
     double fakecorrection = 0. ;
     double fakejer = 0;
     bool isthirdjet = false ;
     

     TLorentzVector Tmpjet;
     
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
	     double correction = 1.;//, old_correction, nominal_correction;
	     double correction_SF = 1.;
	     if (isData == 1) correction = JetCorrector_data->getCorrection();
	     else correction = JetCorrector->getCorrection();
	     
	     
	     //JEC uncertainties
	     double uncertainty = 0 ;
	     if(isData){
	     unc->setJetEta(jetEtaAK4->at(j));
	     unc->setJetPt(jetPtAK4->at(j)/jetJecAK4->at(j)*correction);
	     uncertainty = unc->getUncertainty(true);
	     }else{
	     unc_MC->setJetEta(jetEtaAK4->at(j));
	     unc_MC->setJetPt(jetPtAK4->at(j)/jetJecAK4->at(j)*correction);
	     uncertainty = unc_MC->getUncertainty(true);
	     
	     }
	     jecUncertainty.push_back(uncertainty); 

	     
	     //use "shifted" JECs for study of systematic uncertainties 
	     if( int(getPreCutValue1("shiftJECs"))==1 ){
	       
	       
	       //shift of the corresponding unc
	       correction = correction + getPreCutValue2("shiftJECs")*uncertainty;
	       
	       
	   }
	   
	   //Applying JER onto MC
	   if(!isData && int(getPreCutValue1("useJERs"))==1){
	       
	       JME::JetParameters parameters_1 = {{JME::Binning::JetPt, correction*jetPtAK4->at(j)/jetJecAK4->at(j)},{JME::Binning::JetEta, jetEtaAK4->at(j)}, {JME::Binning::Rho, rho}};
	      

	       float res = resolution->getResolution(parameters_1);
	       
	       
	       JME::JetParameters parameters = {{JME::Binning::JetEta, jetEtaAK4->at(j)}, {JME::Binning::Rho, rho}};
	       
	       
               float scalefactor = resolution_sf->getScaleFactor(parameters);;
	   
	       if (jetPtGenAK4->at(j)>=0.)
        {
            // Smearing is done as here [1]
            //[1] https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_8/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#L236-L237
            double const jerFactor = 1. + (scalefactor - 1.) * ((correction*jetPtAK4->at(j)/jetJecAK4->at(j)) - jetPtGenAK4->at(j)) / (correction*jetPtAK4->at(j)/jetJecAK4->at(j));
            
            // Different definition for debugging with PEC tuples 3.1.0
            // double const jecCorrE = jet.RawP4().E() * jecCorrPt / jet.RawP4().Pt();
            // double const jerFactor = 1. + (jerSF - 1.) * (jecCorrE - genJet->E()) / jecCorrE;
            
             correction_SF = jerFactor;
        }else 
        {
            // Follow the same approach as here [1]
            //[1] https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_8/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#L244_L250
            double const ptResolution = res;
            TF1* rGen = new TF1("Gaus_jer","gaus",-2.,2.);  
            rGen -> SetParameters(1,0.,ptResolution);
            double randomnumber = rGen->GetRandom();
            double const jerFactor = 1. + randomnumber * std::sqrt(std::max(std::pow(scalefactor, 2) - 1., 0.));
           // std::cout<<" random number "<<randomnumber<<" jer "<< ptResolution<<" sigma gaussian "<< rGen->GetParameter(2) <<std::endl;
            correction_SF = jerFactor;
	   
	   
	   delete rGen;
	   
	   
	   }
	   
	 }
         
	         jecFactors.push_back(correction);
	         jerFactors.push_back(correction_SF);
	         if(!isData && int(getPreCutValue1("useJERs"))==1) {
	         CorrFactors.push_back(correction*correction_SF);
	         }else{
	         CorrFactors.push_back(correction);
	         }
	  
	  Tmpjet.SetPtEtaPhiM(jetPtAK4->at(j)/jetJecAK4->at(j)*correction, jetEtaAK4->at(j), jetPhiAK4->at(j), jetMassAK4->at(j)/jetJecAK4->at(j)*correction);
	  
	 sortedJets.insert(std::make_pair((jetPtAK4->at(j)/jetJecAK4->at(j))*correction*correction_SF, j));
	 if( gamma1.DeltaR(Tmpjet) < 0.4)
	      { 
	         nfakejet++ ;
	         fakejetIdx = j;
	         fakecorrection = correction;
	         fakejer        = correction_SF;
	        
	      }
	      

       }
       int size_all = sortedJets.size();
       // remove out tight  photon from jet collection
       if(nfakejet !=0 )
       {
          
          sortedJets.erase((jetPtAK4->at(fakejetIdx)/jetJecAK4->at(fakejetIdx))*fakecorrection*fakejer);


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
	    && (CorrFactors[sortedJetIdx[ijet]]/jetJecAK4->at(sortedJetIdx[ijet]))*jetPtAK4->at(sortedJetIdx[ijet]) >= getPreCutValue1("ptCut"))
	   {
	     Nak4 += 1;
	     HTak4 += (CorrFactors[sortedJetIdx[ijet]]/jetJecAK4->at(sortedJetIdx[ijet]))*jetPtAK4->at(sortedJetIdx[ijet]);
	   }
       }



    int sndjetidx = -10 ;
    bool hasgen = false;
    bool hasgen2 = false;
     //AK4 jets , Select the two first jet for Photon + Jets analysis
    if(no_jets_ak4-nfakejet>=1  ) 
      { 
     
      
      if(!isData && jetPtGenAK4->at(sortedJetIdx[0])>=0.){ hasgen = 1;}
     
         
         
	 if((CorrFactors[sortedJetIdx[0]]/jetJecAK4->at(sortedJetIdx[0]))*jetPtAK4->at(sortedJetIdx[0]) >=  15.  && isNewonOldValidJetTight(jetEtaAK4->at(sortedJetIdx[0]), jetNhfAK4->at(sortedJetIdx[0]), neMultAK4->at(sortedJetIdx[0]), chMultAK4->at(sortedJetIdx[0]), jetNemfAK4->at(sortedJetIdx[0]), jetMufAK4->at(sortedJetIdx[0]), jetChfAK4->at(sortedJetIdx[0]), jetCemfAK4->at(sortedJetIdx[0]),idLAK4->at(sortedJetIdx[0]), isData, hasgen,true))
	   {
		 ak4j1.SetPtEtaPhiM( (CorrFactors[sortedJetIdx[0]]/*1.*//jetJecAK4->at(sortedJetIdx[0])) *jetPtAK4->at(sortedJetIdx[0]) ,jetEtaAK4->at(sortedJetIdx[0])
				     ,jetPhiAK4->at(sortedJetIdx[0])
				     , (CorrFactors[sortedJetIdx[0]]/*1.*//jetJecAK4->at(sortedJetIdx[0])) *jetMassAK4->at(sortedJetIdx[0]));
		
		
		for(size_t secjet = 1 ; secjet < no_jets_ak4-nfakejet ; secjet++ ){	
		
                if(!isData && jetPtGenAK4->at(sortedJetIdx[secjet] > 0.)){ hasgen2 = 1;}	     
		if(no_jets_ak4-nfakejet >= secjet + 1 && (CorrFactors[sortedJetIdx[secjet]]/jetJecAK4->at(sortedJetIdx[secjet]))*jetPtAK4->at(sortedJetIdx[secjet]) >= 10 && isNewonOldValidJetTight(jetEtaAK4->at(sortedJetIdx[secjet]), jetNhfAK4->at(sortedJetIdx[secjet]), neMultAK4->at(sortedJetIdx[secjet]), chMultAK4->at(sortedJetIdx[secjet]), jetNemfAK4->at(sortedJetIdx[secjet]), jetMufAK4->at(sortedJetIdx[secjet]), jetChfAK4->at(sortedJetIdx[secjet]), jetCemfAK4->at(sortedJetIdx[secjet]),idLAK4->at(sortedJetIdx[0]), isData,hasgen2,false))
	       {
		      ak4j2.SetPtEtaPhiM( (CorrFactors[sortedJetIdx[secjet]]/*1.*//jetJecAK4->at(sortedJetIdx[secjet])) *jetPtAK4->at(sortedJetIdx[secjet])
				     ,jetEtaAK4->at(sortedJetIdx[secjet])
				     ,jetPhiAK4->at(sortedJetIdx[secjet])
				     , (CorrFactors[sortedJetIdx[secjet]]/*1.*//jetJecAK4->at(sortedJetIdx[secjet])) *jetMassAK4->at(sortedJetIdx[secjet]));
				     sndjetidx = secjet;
				     break;
				     
	       }ak4j2.SetPtEtaPhiM(0.,0.,0.,0.);}
	   }else{
	   ncut_ptjet++;
	   continue;
	   }
   
     }else{ncut_jet++;
     continue;}
     bool keep_event = true ; 
     //2017 hot zone cleaning 
     if(isData){
         EtaPhiJet_beforehot->Fill(ak4j1.Eta(),ak4j1.Phi());
         EtaPhiJet_beforehot->Fill(ak4j2.Eta(),ak4j2.Phi());
     }
     //cleaning needed for 2017, have to be removed for 2016 
     // if (isData && (h_hotjets.GetBinContent(h_hotjets.FindBin(ak4j1.Eta(), ak4j1.Phi())) >= 10 || h_hotjets.GetBinContent(h_hotjets.FindBin(ak4j1.Eta(), ak4j1.Phi())) >= 10 )) keep_event=false; //
     if(!keep_event){      
        is_hot_area++;
        continue;
     }
     if (isData && (ak4j1.Eta() >= 2.853 &&  ak4j1.Eta() <= 2.964 && ak4j1.Phi()>= 0.6 && ak4j1.Phi() <= 1.)){
        is_hot_area++;
        continue;
     }
     if (isData && (ak4j1.Eta() >= 2.853 &&  ak4j1.Eta() <= 2.964 && ak4j1.Phi()>= 2.2 && ak4j1.Phi() <= 2.6)){
        is_hot_area++;
        continue;
     }
     if (isData && (ak4j1.Eta() >= 2.853 &&  ak4j1.Eta() <= 2.964 && ak4j1.Phi()>= -2.6 && ak4j1.Phi() <= -2.2)){
        is_hot_area++;
        continue;
     }
     if (isData && (ak4j1.Eta() >= 2.853 &&  ak4j1.Eta() <= 2.964 && ak4j1.Phi()>= 2.9 && ak4j1.Phi() <= 3.1)){
        is_hot_area++;
        continue;
     }
     
     if (isData && (ak4j1.Eta() >= -2.964 &&  ak4j1.Eta() <= -2.853 && ak4j1.Phi()>= -2.6 && ak4j1.Phi() <= -2.2)){
        is_hot_area++;
        continue;
     }
     if (isData && (ak4j1.Eta() >= -2.964 &&  ak4j1.Eta() <= -2.853 && ak4j1.Phi()>= 0.6 && ak4j1.Phi() <= 1.)){
        is_hot_area++;
        continue;
     }
     if (isData && (ak4j1.Eta() >= -2.964 &&  ak4j1.Eta() <= -2.853 && ak4j1.Phi()>= 2.2 && ak4j1.Phi() <= 2.6)){
        is_hot_area++;
        continue;
     }
     if (isData && (ak4j1.Eta() >= -2.964 &&  ak4j1.Eta() <= -2.853 && ak4j1.Phi()>= 2.9 && ak4j1.Phi() <= 3.1)){
        is_hot_area++;
        continue;
     }
     
     // jet2
     
     if (isData && (ak4j2.Eta() >= 2.853 &&  ak4j2.Eta() <= 2.964 && ak4j2.Phi()>= 0.6 && ak4j2.Phi() <= 1.)){
        is_hot_area++;
        continue;
     }
     if (isData && (ak4j2.Eta() >= 2.853 &&  ak4j2.Eta() <= 2.964 && ak4j2.Phi()>= 2.2 && ak4j2.Phi() <= 2.6)){
        is_hot_area++;
        continue;
     }
     if (isData && (ak4j2.Eta() >= 2.853 &&  ak4j2.Eta() <= 2.964 && ak4j2.Phi()>= -2.6 && ak4j2.Phi() <= -2.2)){
        is_hot_area++;
        continue;
     }
     if (isData && (ak4j2.Eta() >= 2.853 &&  ak4j2.Eta() <= 2.964 && ak4j2.Phi()>= 2.9 && ak4j2.Phi() <= 3.1)){
        is_hot_area++;
        continue;
     }
     
     if (isData && (ak4j2.Eta() >= -2.964 &&  ak4j2.Eta() <= -2.853 && ak4j2.Phi()>= -2.6 && ak4j2.Phi() <= -2.2)){
        is_hot_area++;
        continue;
     }
     if (isData && (ak4j2.Eta() >= -2.964 &&  ak4j2.Eta() <= -2.853 && ak4j2.Phi()>= 0.6 && ak4j2.Phi() <= 1.)){
        is_hot_area++;
        continue;
     }
     if (isData && (ak4j2.Eta() >= -2.964 &&  ak4j2.Eta() <= -2.853 && ak4j2.Phi()>= 2.2 && ak4j2.Phi() <= 2.6)){
        is_hot_area++;
        continue;
     }
     if (isData && (ak4j2.Eta() >= -2.964 &&  ak4j2.Eta() <= -2.853 && ak4j2.Phi()>= 2.9 && ak4j2.Phi() <= 3.1)){
        is_hot_area++;
        continue;
     }
     
     if(isData){
         EtaPhiJet_afterhot->Fill(ak4j1.Eta(),ak4j1.Phi());
         EtaPhiJet_afterhot->Fill(ak4j2.Eta(),ak4j2.Phi());
     }
     //store all corrected jets of the events 
     for(int i = 0 ; i < no_jets_ak4-nfakejet; i++){
     
        if(!isData && jetPtGenAK4->at(sortedJetIdx[0]) > 0.){ hasgen = 1;}
        if((CorrFactors[sortedJetIdx[i]]/jetJecAK4->at(sortedJetIdx[i]))*jetPtAK4->at(sortedJetIdx[i]) >=  2. )
	   {
		 pt_jets_   ->push_back( (CorrFactors[sortedJetIdx[i]]/jetJecAK4->at(sortedJetIdx[i])) *jetPtAK4->at(sortedJetIdx[i])) ;
		 eta_jets_  ->push_back(jetEtaAK4->at(sortedJetIdx[i]));
		 phi_jets_  ->push_back(jetPhiAK4->at(sortedJetIdx[i]));
		 mass_jets_ ->push_back((CorrFactors[sortedJetIdx[i]]/jetJecAK4->at(sortedJetIdx[i])) *jetMassAK4->at(sortedJetIdx[i]));
                 emF_jets_  ->push_back(jetNemfAK4->at(sortedJetIdx[i]) + jetCemfAK4->at(sortedJetIdx[i]));
                 bool isID = false;
                if(isNewonOldValidJetTight(jetEtaAK4->at(sortedJetIdx[i]), jetNhfAK4->at(sortedJetIdx[i]), neMultAK4->at(sortedJetIdx[i]), chMultAK4->at(sortedJetIdx[i]), jetNemfAK4->at(sortedJetIdx[i]), jetMufAK4->at(sortedJetIdx[i]), jetChfAK4->at(sortedJetIdx[i]), jetCemfAK4->at(sortedJetIdx[i]),idLAK4->at(sortedJetIdx[i]), isData, hasgen,false)) isID = true; 
                 IsID_jets_ ->push_back(isID);
           }
      }
     
     




 //----------------------TYPE I MET computation------------------
 
  double deltaPx = 0., deltaPy = 0.;
  TLorentzVector  jetRC,corrJet; 
  
  if(nPhotonTight == 1){
  for (Long64_t it=0; it< no_jets_ak4-nfakejet ; it++) {
    if(fakejetIdx==it)continue;
    
    
    double corrs = 1.;
    double corrsForTypeI = 1.;
    double correction_SF = 1.;
    
    
   
      if(isData)
      {
      JetCorrectortypI->setJetEta(jetEtaAK4->at(it));
      JetCorrectortypI->setJetPt(jetPtAK4->at(it)/jetJecAK4->at(it));
      JetCorrectortypI->setJetA(jetAreaAK4->at(it));
      JetCorrectortypI->setRho(rho);
      corrsForTypeI = JetCorrectortypI->getCorrection(); //only RC
      
      
      //JEC uncertainties
      unc_RC->setJetEta(jetEtaAK4->at(it));
      unc_RC->setJetPt(jetPtAK4->at(it)/jetJecAK4->at(it)*corrsForTypeI);
      double uncertainty = unc_RC->getUncertainty(true);

	     
      //use "shifted" JECs for study of systematic uncertainties 
      if( int(getPreCutValue1("shiftJECs"))==1 ){
	       
	       
        //shift of the corresponding unc
        corrsForTypeI = corrsForTypeI + getPreCutValue2("shiftJECs")*uncertainty;
	       
	       
     }

      JetCorrectortypIL123 ->setJetEta(jetEtaAK4->at(it));
      JetCorrectortypIL123 ->setJetPt(jetPtAK4->at(it)/jetJecAK4->at(it));
      JetCorrectortypIL123 ->setJetA(jetAreaAK4->at(it));
      JetCorrectortypIL123 ->setRho(rho);
      corrs = JetCorrectortypIL123->getCorrection(); // L1L2L3
      
      //JEC uncertainties
      unc_TI->setJetEta(jetEtaAK4->at(it));
      unc_TI->setJetPt(jetPtAK4->at(it)/jetJecAK4->at(it)*corrs);
      double uncertainty_TI = unc_TI->getUncertainty(true);

	     
      //use "shifted" JECs for study of systematic uncertainties 
      if( int(getPreCutValue1("shiftJECs"))==1 ){
	       
	       
        //shift of the corresponding unc
        corrs = corrs + getPreCutValue2("shiftJECs")*uncertainty_TI;
	       
	       
     }
      
      
      }else{

      JetCorrectortypIMC->setJetEta(jetEtaAK4->at(it));
      JetCorrectortypIMC->setJetPt(jetPtAK4->at(it)/jetJecAK4->at(it));
      JetCorrectortypIMC->setJetA(jetAreaAK4->at(it));
      JetCorrectortypIMC->setRho(rho);
      corrsForTypeI = JetCorrectortypIMC->getCorrection(); //only RC
      
      //JEC uncertainties
      unc_MC_RC->setJetEta(jetEtaAK4->at(it));
      unc_MC_RC->setJetPt(jetPtAK4->at(it)/jetJecAK4->at(it)*corrsForTypeI);
      double uncertainty = unc_MC_RC->getUncertainty(true);

	     
      //use "shifted" JECs for study of systematic uncertainties 
      if( int(getPreCutValue1("shiftJECs"))==1 ){
	       
	       
        //shift of the corresponding unc
        corrsForTypeI = corrsForTypeI + getPreCutValue2("shiftJECs")*uncertainty;
	       
	       
     }
      
      
      JetCorrectortypIL123MC ->setJetEta(jetEtaAK4->at(it));
      JetCorrectortypIL123MC ->setJetPt(jetPtAK4->at(it)/jetJecAK4->at(it));
      JetCorrectortypIL123MC ->setJetA(jetAreaAK4->at(it));
      JetCorrectortypIL123MC ->setRho(rho);
      corrs = JetCorrectortypIL123MC->getCorrection(); // L1L2L3
      //JEC uncertainties
      unc_MC_TI->setJetEta(jetEtaAK4->at(it));
      unc_MC_TI->setJetPt(jetPtAK4->at(it)/jetJecAK4->at(it)*corrs);
      double uncertainty_TI = unc_MC_TI->getUncertainty(true);

	     
      //use "shifted" JECs for study of systematic uncertainties 
      if( int(getPreCutValue1("shiftJECs"))==1 ){
	       
	       
        //shift of the corresponding unc
        corrs = corrs + getPreCutValue2("shiftJECs")*uncertainty_TI;
	       
	       
     }
      
      //Applying JER onto MC
	  
             if(int(getPreCutValue1("useJERs"))==1 ){
	       JME::JetParameters parameters_1 = {{JME::Binning::JetPt, corrs*jetPtAK4->at(it)/jetJecAK4->at(it)},{JME::Binning::JetEta, jetEtaAK4->at(it)}, {JME::Binning::Rho, rho}};

	       float res = resolution_TI->getResolution(parameters_1);
	       
	       
	       JME::JetParameters parameters = {{JME::Binning::JetEta, jetEtaAK4->at(it)}, {JME::Binning::Rho, rho}};
	       
	       
               float scalefactor = resolution_sf_TI->getScaleFactor(parameters);;
	   
	       if (jetPtGenAK4->at(it)>=0.)
        {
            // Smearing is done as here [1]
            //[1] https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_8/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#L236-L237
            double const jerFactor = 1. + (scalefactor - 1.) * ((corrs*jetPtAK4->at(it)/jetJecAK4->at(it)) - jetPtGenAK4->at(it)) / (corrs*jetPtAK4->at(it)/jetJecAK4->at(it));
            
            
             correction_SF = jerFactor;
        }else 
        {
            // Follow the same approach as here [1]
            //[1] https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_8/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#L244_L250
            double const ptResolution = res;
            TF1* rGen = new TF1("Gaus_jer","gaus",-2.,2.);  
            rGen -> SetParameters(1,0.,ptResolution);
            double randomnumber = rGen->GetRandom();
            double const jerFactor = 1. + randomnumber * std::sqrt(std::max(std::pow(scalefactor, 2) - 1., 0.));
           // std::cout<<" random number "<<randomnumber<<" jer "<< ptResolution<<" sigma gaussian "<< rGen->GetParameter(2) <<std::endl;
            correction_SF = jerFactor;
	   
	   
	   delete rGen;
	   
	   
	   }}
	   
	 

      }

    jetRC.SetPtEtaPhiE((jetPtAK4->at(it)/jetJecAK4->at(it))*corrsForTypeI,jetEtaAK4->at(it),jetPhiAK4->at(it),(jetEnergyAK4->at(it)/jetJecAK4->at(it))*corrsForTypeI);
     // only RC
      
    corrJet.SetPtEtaPhiE((jetPtAK4->at(it)/jetJecAK4->at(it))*corrs*correction_SF,jetEtaAK4->at(it),jetPhiAK4->at(it),(jetEnergyAK4->at(it)/jetJecAK4->at(it))*corrs*correction_SF);
     
     //Type I modified MET for 2017
     if(corrJet.Pt() < 75. && corrJet.Eta() > 2.605 && corrJet.Eta() < 3.139) continue ; 
     
    double dR = gamma1.DeltaR(corrJet);
   
    
    if(corrJet.Pt() > 15. && dR > 0.25) {
	
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
  }else{
  MetTypeI.SetPxPyPzE(0,0,0,0);
  
  } 
 
 

//------------------------END met calculation------------------


  double Rbal = -999. ;
  double Rmpf = -999. ;
  double Rmpfraw = -999. ;

  Rbal = (ak4j1.Pt()/gamma1.Pt());
  double deltPHIgj = fabs(gamma1.DeltaPhi(ak4j1));
  double deltaphiPhomet = (gamma1.Phi()- metPhi);

  Rmpf = 1. + (MetTypeI.Px()*gamma1.Px()+MetTypeI.Py()*gamma1.Py())/std::pow(gamma1.Pt(),2);
  Rmpfraw = 1. + metPt*gamma1.Pt()*cos(deltaphiPhomet)/std::pow(gamma1.Pt(),2);
  double alpha = -999. ;
  
  alpha = (ak4j2.Pt()/gamma1.Pt());

 DeltaPhiAlpha->Fill(alpha,deltPHIgj);

 
 
     if(deltPHIgj>=2.8 )
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
      bool hasgen = false;
      if(!isData && jetPtGenAK4->at(sortedJetIdx[0])){ hasgen = 1;}
     if(ak4j1.Pt()>0 && isNewonOldValidJetTight(jetEtaAK4->at(sortedJetIdx[0]), jetNhfAK4->at(sortedJetIdx[0]), neMultAK4->at(sortedJetIdx[0]), chMultAK4->at(sortedJetIdx[0]), jetNemfAK4->at(sortedJetIdx[0]), jetMufAK4->at(sortedJetIdx[0]), jetChfAK4->at(sortedJetIdx[0]), jetCemfAK4->at(sortedJetIdx[0]),idLAK4->at(sortedJetIdx[0]), isData, hasgen,false)) 
    {
     nselectedevent++;
     
     if(lumi != storelumi)
      {
        totalluminosity += lumi ;
      }
      storelumi = lumi;

     
 //store everything variable needed for the final computation    
     
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
     fillVariableWithValue("sigmaietaieta_photon", Photonfull5x5SigmaIEtaIEtaMapToken->at(indexgoodpho));
     fillVariableWithValue("CHiso_photon", PhotonphoChargedIsolationToken->at(indexgoodpho));
     fillVariableWithValue("NHiso_photon", PhotonphoNeutralHadronIsolationToken->at(indexgoodpho));
     fillVariableWithValue("Photoniso_photon", PhotonphoPhotonIsolationToken->at(indexgoodpho));
     fillVariableWithValue("Pt_photonSC",PhotonSCPt->at(indexgoodpho));
     fillVariableWithValue("Eta_photonSC",PhotonSCEta->at(indexgoodpho));
     fillVariableWithValue("Phi_photonSC",PhotonSCPhi->at(indexgoodpho));
     fillVariableWithValue("Energy_photonSC",PhotonSCEnergy->at(indexgoodpho));
     fillVariableWithValue("hadTowOverEm",hadTowOverEm->at(indexgoodpho));

     if(isData)
       {
          fillVariableWithValue("weight", 1);
          
        }
       if(!isData)
       {
         fillVariableWithValue("Pt_photonGEN",photonPtGen->at(indexgoodpho));
	 fillVariableWithValue("Eta_photonGEN",photonEtaGen->at(indexgoodpho));
	 fillVariableWithValue("Phi_photonGEN",photonPhiGen->at(indexgoodpho));
	 fillVariableWithValue("Energy_photonGEN",photonEnergyGen->at(indexgoodpho));
	 fillVariableWithValue("PassGenmatching",isGenMatch->at(indexgoodpho));

	 fillVariableWithValue("weight",weight);
	 
       }//MC
     

     if( AK4jets.size() >=1 )
     {


       fillVariableWithValue( "IdTight_j1",idLAK4->at(sortedJetIdx[0]));
       fillVariableWithValue( "pTAK4_j1", AK4jets[0].Pt());
       fillVariableWithValue( "etaAK4_j1", AK4jets[0].Eta());
       fillVariableWithValue( "phiAK4_j1", AK4jets[0].Phi());
       

       
       fillVariableWithValue( "jetJecAK4_j1", jecFactors[sortedJetIdx[0]] );
       fillVariableWithValue( "jetJecUncAK4_j1", jecUncertainty[sortedJetIdx[0]] );
       fillVariableWithValue( "jetCSVAK4_j1", jetCSVAK4->at(sortedJetIdx[0]) );
       fillVariableWithValue( "jetQGDAK4_j1", jetQGDAK4->at(sortedJetIdx[0]) );
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
       if(!isData){
       
       fillVariableWithValue( "jetJerAK4_j1", jerFactors[sortedJetIdx[0]] );
       
       }
       if(!isData && hasgen )
       {
         fillVariableWithValue( "pTAK4_j1GEN", jetPtGenAK4->at(sortedJetIdx[0]));
	 fillVariableWithValue( "etaAK4_j1GEN", jetEtaGenAK4->at(sortedJetIdx[0]));	       
	 fillVariableWithValue("PDGIDAK4_j1",jetpdgIDGenAK4->at(sortedJetIdx[0]));

       }
     }
     
     if( AK4jets.size() >=2){

       fillVariableWithValue( "IdTight_j2",idLAK4->at(sortedJetIdx[sndjetidx]));
       fillVariableWithValue( "pTAK4_j2", AK4jets[1].Pt() );
       fillVariableWithValue( "etaAK4_j2", AK4jets[1].Eta());
       fillVariableWithValue( "phiAK4_j2", AK4jets[1].Phi());

       fillVariableWithValue( "jetJecAK4_j2", jecFactors[sortedJetIdx[sndjetidx]]); 
       fillVariableWithValue( "jetJecUncAK4_j2", jecUncertainty[sortedJetIdx[sndjetidx]] );
       fillVariableWithValue( "jetCSVAK4_j2", jetCSVAK4->at(sortedJetIdx[sndjetidx]) );
       fillVariableWithValue( "jetQGDAK4_j2", jetQGDAK4->at(sortedJetIdx[sndjetidx]) );
       //jetID
       fillVariableWithValue( "neutrHadEnFrac_j2", jetNhfAK4->at(sortedJetIdx[sndjetidx]));
       fillVariableWithValue( "chargedHadEnFrac_j2", jetChfAK4->at(sortedJetIdx[sndjetidx]));
       fillVariableWithValue( "photonEnFrac_j2", jetPhfAK4->at(sortedJetIdx[sndjetidx]));
       fillVariableWithValue( "eleEnFract_j2", jetElfAK4->at(sortedJetIdx[sndjetidx]));
       fillVariableWithValue( "muEnFract_j2", jetMufAK4->at(sortedJetIdx[sndjetidx]));
       fillVariableWithValue( "neutrElectromFrac_j2", jetNemfAK4->at(sortedJetIdx[sndjetidx]));
       fillVariableWithValue( "chargedElectromFrac_j2", jetCemfAK4->at(sortedJetIdx[sndjetidx]));
       fillVariableWithValue( "chargedMult_j2", chMultAK4->at(sortedJetIdx[sndjetidx]));
       fillVariableWithValue( "neutrMult_j2", neMultAK4->at(sortedJetIdx[sndjetidx]));
       fillVariableWithValue( "photonMult_j2", phoMultAK4->at(sortedJetIdx[sndjetidx]));
       //dijet
       fillVariableWithValue( "CosThetaStarAK4", TMath::TanH( (AK4jets[0].Eta()-AK4jets[1].Eta())/2 ));
       if(!isData)fillVariableWithValue( "jetJerAK4_j2", jerFactors[sortedJetIdx[sndjetidx]] );
       if(!isData && hasgen2)
       {
          

         fillVariableWithValue( "pTAK4_j2GEN", jetPtGenAK4->at(sortedJetIdx[sndjetidx]));
	 fillVariableWithValue( "etaAK4_j2GEN", jetEtaGenAK4->at(sortedJetIdx[sndjetidx]));
 	 fillVariableWithValue("PDGIDAK4_j2",jetpdgIDGenAK4->at(sortedJetIdx[sndjetidx]));
	 

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
       fillVariableWithValue("METGEN",metEnergyPUPPIGen);
       fillVariableWithValue("METGEN_Pt",metPtPUPPIGen);
       fillVariableWithValue("METGEN_Eta",metEtaPUPPIGen);
       fillVariableWithValue("METGEN_Phi",metPhiPUPPIGen);
     }
     

     double METoverHTAK4=double(MetTypeI.Et()/HTak4);
     fillVariableWithValue("METoverHTAK4",METoverHTAK4);

     fillVariableWithValue("HTAK4",HTak4);
     fillVariableWithValue("ptHat",ptHat);

     // Trigger

     if( isMatch30->size() > 0 && indexgoodpho < isMatch30->size() )//&& isData)
     {
       fillVariableWithValue("phomatchHLT_Photon30",isMatch30->at(indexgoodpho));// 
       
     }
    if( isMatch50->size() > 0 && indexgoodpho < isMatch50->size() )//&& isData)
    {
       fillVariableWithValue("phomatchHLT_Photon50",isMatch50->at(indexgoodpho));//

     }
     if( isMatch75->size() > 0 && indexgoodpho < isMatch75->size() )//&& isData)
     { 
       fillVariableWithValue("phomatchHLT_Photon75",isMatch75->at(indexgoodpho));// 

     }
     if( isMatch90->size() > 0 && indexgoodpho < (isMatch90->size() )  )// && isData)
     {
       fillVariableWithValue("phomatchHLT_Photon90",isMatch90->at(indexgoodpho));//

     }
     if( isMatch120->size() > 0 && indexgoodpho < isMatch120->size() )//&& isData)
     {
       fillVariableWithValue("phomatchHLT_Photon120",isMatch120->at(indexgoodpho));//

     }
     if( isMatch165->size() > 0 && indexgoodpho < isMatch165->size() )// && isData)
     {
       fillVariableWithValue("phomatchHLT_Photon165",isMatch165->at(indexgoodpho));//
     }

     
          int NtriggerBits = triggerResult->size();
     if( NtriggerBits > 0 )//&& isData)
     {
       fillVariableWithValue("passHLT_Photon30",triggerResult->at(0));// 
     }
    if( NtriggerBits > 1 )//&& isData)
    {
       fillVariableWithValue("passHLT_Photon50",triggerResult->at(1));//

     }
     if( NtriggerBits > 2 )//&& isData)
     {
       fillVariableWithValue("passHLT_Photon75",triggerResult->at(2));// 

     }
     if( NtriggerBits > 3)// && isData)
     {
       fillVariableWithValue("passHLT_Photon90",triggerResult->at(3));//

     }
     if( NtriggerBits > 4 )//&& isData)
     {
       fillVariableWithValue("passHLT_Photon120",triggerResult->at(4));//

     }
     if( NtriggerBits > 5)// && isData)
     {
       fillVariableWithValue("passHLT_Photon165",triggerResult->at(5));//

     }
     if( NtriggerBits > 6)// && isData)
     {
       fillVariableWithValue("passHLT_Photon200",triggerResult->at(6));//

     }


     // Evaluate cuts (but do not apply them)
     evaluateCuts();
     
     

     // optional call to fill a skim with a subset of the variables defined in the cutFile (use flag SAVE)
     if( passedCut("Pt_photon") ) 
       {
	 fillReducedSkimTree();

       }
}
//}else{ncut_alpha++;
//continue;}
}else{ncut_deltaphi++;
continue;}//end dphi cut
//end of pt j1 cut
     ////////////////////// User's code ends here ///////////////////////
     
     
}else{ncut_muons++;
continue;}
}else{ncut_photonpt++;
continue;}
}else{ncut_pixelseed++;
continue;}
}else{ncut_photon++;
 continue;}
//end photon condition
}else{
   ncut_nophoton++;
   continue;
}//no photon 
   } // End loop over events
   
  // reduced_skim_file_->cd();
   totalluminosityP->SetVal(totalluminosity);
  
   SumWeight->Write();
   totalluminosityP->Write();
   PUvariable->Write();
   DeltaPhiAlpha->Write();
   EtaPhiJet_beforehot->Write();
   EtaPhiJet_afterhot->Write();
   
   
    std::cout << std::endl;
  std::cout << "Absolute efficiency : related to initial number of event =  " << testcount<< std::endl;
  std::cout << "Nevent after photon presence cut: " << MAKE_RED << testcount - (double) Vtxcut -  (double)  ncut_nophoton    << RESET_COLOR << std::endl;
  std::cout << "Nevent after photon cut: " << MAKE_RED <<  testcount - (double) Vtxcut - (double)  ncut_nophoton - (double) ncut_photon  << RESET_COLOR << std::endl;
  std::cout << "Nevent after pixel seed veto cut: " << MAKE_RED <<  testcount /*- (double) Vtxcut - (double)  ncut_nophoton - (double) ncut_photon */-(double) ncut_pixelseed  << RESET_COLOR << std::endl;    
  std::cout << "Nevent after  pt photon cut: " << MAKE_RED <<   testcount /*- (double) Vtxcut - (double)  ncut_nophoton - (double) ncut_photon */-(double) ncut_pixelseed  - (double) ncut_photonpt  << RESET_COLOR << std::endl;
  std::cout << "Nevent after muons cut: " << MAKE_RED << testcount /*- (double) Vtxcut - (double)  ncut_nophoton - (double) ncut_photon */-(double) ncut_pixelseed  - (double) ncut_photonpt - (double) ncut_muons  << RESET_COLOR << std::endl; 
  std::cout << "Nevent after electrons cut: " << MAKE_RED << testcount  /*- (double) Vtxcut - (double)  ncut_nophoton - (double) ncut_photon */-(double) ncut_pixelseed  - (double) ncut_photonpt - (double) ncut_muons - (double) ncut_electron << RESET_COLOR << std::endl; 
  std::cout << "Nevent after n jet cut: " << MAKE_RED <<testcount  /*- (double) Vtxcut - (double)  ncut_nophoton - (double) ncut_photon */-(double) ncut_pixelseed  - (double) ncut_photonpt - (double) ncut_muons - (double) ncut_electron - (double) ncut_jet << RESET_COLOR << std::endl;  
  std::cout << "Nevent after Pt(j1) cut: " << MAKE_RED << testcount /*- (double) Vtxcut - (double)  ncut_nophoton - (double) ncut_photon */ -(double) ncut_pixelseed  - (double) ncut_photonpt - (double) ncut_muons - (double) ncut_electron - (double) ncut_jet - (double)  ncut_ptjet<< RESET_COLOR << std::endl;   
  std::cout << "Nevent after Δφ cut: " << MAKE_RED << testcount /*- (double) Vtxcut - (double)  ncut_nophoton - (double) ncut_photon */-(double) ncut_pixelseed  - (double) ncut_photonpt - (double) ncut_muons - (double) ncut_electron - (double) ncut_jet - (double)  ncut_ptjet - (double) ncut_deltaphi  << RESET_COLOR << std::endl;
  std::cout << "Nevent after α cut: " << MAKE_RED << testcount /*- (double) Vtxcut - (double)  ncut_nophoton - (double) ncut_photon */ -(double) ncut_pixelseed  - (double) ncut_photonpt - (double) ncut_muons - (double) ncut_electron - (double) ncut_jet - (double)  ncut_ptjet - (double) ncut_deltaphi - (double) ncut_alpha  << RESET_COLOR << std::endl;  
  
  std::cout << "Event rejected because of PU: " << MAKE_RED << Is_PU << RESET_COLOR << std::endl;
  std::cout << "Event rejected because of Hot tower: " << MAKE_RED << is_hot_area << RESET_COLOR << std::endl;
  std::cout << "Selection efficiency: " << MAKE_RED << (double) nselectedevent  << RESET_COLOR << std::endl;
  std::cout << std::endl;
   
   std::cout<<" nb of selected event " << nselectedevent<<std::endl;
   
   
   
   
   
   
   
   
   
   std::cout << "analysisClass::Loop() ends" <<std::endl; 
   
    delete  PUvariable;
    delete pt_jets_;
    delete phi_jets_;
    delete eta_jets_;
    delete mass_jets_;
    delete emF_jets_;
    delete IsID_jets_;
    delete DeltaPhiAlpha;
    delete SumWeight;
    delete EtaPhiJet_beforehot;
    delete EtaPhiJet_afterhot;
    //delete h_hotjets;
     
    
     

   }
