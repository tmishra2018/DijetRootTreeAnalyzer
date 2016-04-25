#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TF1.h>

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile)
{
  std::cout << "analysisClass::analysisClass(): begins " << std::endl;
  
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

   // TH1F* h_mjj_HLTpass[8];
   // char name_histoHLT[50];
   // for (int i=0; i<8; i++){  
   //   sprintf(name_histoHLT,"h_mjj_HLTpass_%s",HLTname[i]);
   //   h_mjj_HLTpass[i]= new TH1F(name_histoHLT,"",103,massBoundaries);
   // }

   //TH1F* h_mjj_fullSel_varBin = new TH1F("h_mjj_fullSel_varBin","",103,massBoundaries);
   //TH1F* h_mjj_fullSel_fixBin = new TH1F("h_mjj_fullSel_fixBin","",getHistoNBins("mjj"),getHistoMin("mjj"),getHistoMax("mjj"));

   //For trigger efficiency measurements

   //var bin
   TH1F* h_mjj_NoTrigger = new TH1F("h_mjj_NoTrigger","",103,massBoundaries);
   //L1
   TH1F* h_mjj_HLTpass_ZeroBias = new TH1F("h_mjj_HLTpass_ZeroBias","",103,massBoundaries);
   TH1F* h_mjj_HLTpass_ZeroBias_L1HTT150 = new TH1F("h_mjj_HLTpass_ZeroBias_L1HTT150","",103,massBoundaries);  
   TH1F* h_mjj_HLTpass_CaloJet40 = new TH1F("h_mjj_HLTpass_CaloJet40","",14000,0,14000);
   TH1F* h_mjj_HLTpass_CaloJet40_L1HTT150 = new TH1F("h_mjj_HLTpass_CaloJet40_L1HTT150","",14000,0,14000);
   //HLT
   TH1F* h_mjj_HLTpass_L1HTT150 = new TH1F("h_mjj_HLTpass_L1HTT150","",103,massBoundaries);
   TH1F* h_mjj_HLTpass_L1HTT150_HT450 = new TH1F("h_mjj_HLTpass_L1HTT150_HT450","",103,massBoundaries);
   TH1F* h_mjj_HLTpass_L1HTT150_p1 = new TH1F("h_mjj_HLTpass_L1HTT150_p1","",103,massBoundaries);
   TH1F* h_mjj_HLTpass_L1HTT150_HT450_p1 = new TH1F("h_mjj_HLTpass_L1HTT150_HT450_p1","",103,massBoundaries);
   TH1F* h_mjj_HLTpass_L1HTT150_p2 = new TH1F("h_mjj_HLTpass_L1HTT150_p2","",103,massBoundaries);
   TH1F* h_mjj_HLTpass_L1HTT150_HT450_p2 = new TH1F("h_mjj_HLTpass_L1HTT150_HT450_p2","",103,massBoundaries);
   TH1F* h_mjj_HLTpass_L1HTT150_p3 = new TH1F("h_mjj_HLTpass_L1HTT150_p3","",103,massBoundaries);
   TH1F* h_mjj_HLTpass_L1HTT150_HT450_p3 = new TH1F("h_mjj_HLTpass_L1HTT150_HT450_p3","",103,massBoundaries);
   TH1F* h_mjj_HLTpass_L1HTT150_p4 = new TH1F("h_mjj_HLTpass_L1HTT150_p4","",103,massBoundaries);
   TH1F* h_mjj_HLTpass_L1HTT150_HT450_p4 = new TH1F("h_mjj_HLTpass_L1HTT150_HT450_p4","",103,massBoundaries);

   //1 GeV bin
   TH1F* h_mjj_NoTrigger_1GeVbin = new TH1F("h_mjj_NoTrigger_1GeVbin","",14000,0,14000);
   //L1
   TH1F* h_mjj_HLTpass_ZeroBias_1GeVbin = new TH1F("h_mjj_HLTpass_ZeroBias_1GeVbin","",14000,0,14000);
   TH1F* h_mjj_HLTpass_ZeroBias_L1HTT150_1GeVbin = new TH1F("h_mjj_HLTpass_ZeroBias_L1HTT150_1GeVbin","",14000,0,14000); 
   TH1F* h_mjj_HLTpass_CaloJet40_1GeVbin = new TH1F("h_mjj_HLTpass_CaloJet40_1GeVbin","",14000,0,14000);
   TH1F* h_mjj_HLTpass_CaloJet40_L1HTT150_1GeVbin = new TH1F("h_mjj_HLTpass_CaloJet40_L1HTT150_1GeVbin","",14000,0,14000);

   //HLT
   TH1F* h_mjj_HLTpass_L1HTT150_1GeVbin = new TH1F("h_mjj_HLTpass_L1HTT150_1GeVbin","",14000,0,14000);
   TH1F* h_mjj_HLTpass_L1HTT150_HT450_1GeVbin = new TH1F("h_mjj_HLTpass_L1HTT150_HT450_1GeVbin","",14000,0,14000);
   TH1F* h_mjj_HLTpass_L1HTT150_1GeVbin_p1 = new TH1F("h_mjj_HLTpass_L1HTT150_1GeVbin_p1","",14000,0,14000);
   TH1F* h_mjj_HLTpass_L1HTT150_HT450_1GeVbin_p1 = new TH1F("h_mjj_HLTpass_L1HTT150_HT450_1GeVbin_p1","",14000,0,14000);
   TH1F* h_mjj_HLTpass_L1HTT150_1GeVbin_p2 = new TH1F("h_mjj_HLTpass_L1HTT150_1GeVbin_p2","",14000,0,14000);
   TH1F* h_mjj_HLTpass_L1HTT150_HT450_1GeVbin_p2 = new TH1F("h_mjj_HLTpass_L1HTT150_HT450_1GeVbin_p2","",14000,0,14000);
   TH1F* h_mjj_HLTpass_L1HTT150_1GeVbin_p3 = new TH1F("h_mjj_HLTpass_L1HTT150_1GeVbin_p3","",14000,0,14000);
   TH1F* h_mjj_HLTpass_L1HTT150_HT450_1GeVbin_p3 = new TH1F("h_mjj_HLTpass_L1HTT150_HT450_1GeVbin_p3","",14000,0,14000);
   TH1F* h_mjj_HLTpass_L1HTT150_1GeVbin_p4 = new TH1F("h_mjj_HLTpass_L1HTT150_1GeVbin_p4","",14000,0,14000);
   TH1F* h_mjj_HLTpass_L1HTT150_HT450_1GeVbin_p4 = new TH1F("h_mjj_HLTpass_L1HTT150_HT450_1GeVbin_p4","",14000,0,14000);

   //mjj correction factor (https://indico.cern.ch/event/515842/contribution/8/attachments/1251169/1845106/Aprile-01-2016_-_Meeting.pdf slide 6)
   TF1 *f_mjjCorr = new TF1("f_mjjCorr","-0.00119706*pow(log(x),3)+0.00910685*pow(log(x),2)+0.0547697*log(x)-0.402532", 0, 14000);
   //TF1 *f_mjjCorr = new TF1("f_mjjCorr","0.", 0, 14000);
   TH2F* h2_f_mjjCorr_vs_mjj = new TH2F("h2_f_mjjCorr_vs_mjj","",14000,0,14000,1000,0.5,1.5);

   /////////initialize variables

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<200000;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%100000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     // if (Cut(ientry) < 0) continue;

     ////////////////////// User's code starts here ///////////////////////

     ///Stuff to be done for every event

     resetCuts();
    
     //== Fill Variables ==
     fillVariableWithValue("PassJSON",PassJSON);     
     fillVariableWithValue("deltaETAjj",deltaETAjj);     
     fillVariableWithValue("mjjCorr", mjj*(1+f_mjjCorr->Eval(mjj)) );     
     fillVariableWithValue("mjj",mjj);     
     fillVariableWithValue("massCorrection", 1+f_mjjCorr->Eval(mjj) );     

     // Evaluate cuts (but do not apply them)
     evaluateCuts();

     if ( passedCut("all") )
       {

	 //correction
	 h2_f_mjjCorr_vs_mjj->Fill(getVariableValue("mjj"),getVariableValue("massCorrection"));

	 //trigger
	 h_mjj_NoTrigger -> Fill(getVariableValue("mjjCorr")); 
	 h_mjj_NoTrigger_1GeVbin -> Fill(getVariableValue("mjjCorr")); 
	 
	 if(passHLT_ZeroBias_BtagSeq || passHLT_ZeroBias)
	   {
	     h_mjj_HLTpass_ZeroBias -> Fill(getVariableValue("mjjCorr"));  
	     h_mjj_HLTpass_ZeroBias_1GeVbin -> Fill(getVariableValue("mjjCorr"));  
	   }
	 
	 if( (passHLT_ZeroBias_BtagSeq || passHLT_ZeroBias) 
	     && (passHLT_L1HTT150_BtagSeq || passHLT_L1HTT150) )
	   {
	     h_mjj_HLTpass_ZeroBias_L1HTT150 -> Fill(getVariableValue("mjjCorr"));  
	     h_mjj_HLTpass_ZeroBias_L1HTT150_1GeVbin -> Fill(getVariableValue("mjjCorr"));  
	   }

	 if(passHLT_CaloJet40_BtagSeq || passHLT_CaloJet40)
	   {
	     h_mjj_HLTpass_CaloJet40 -> Fill(getVariableValue("mjjCorr"));  
	     h_mjj_HLTpass_CaloJet40_1GeVbin -> Fill(getVariableValue("mjjCorr"));  
	   }

	 if( (passHLT_CaloJet40_BtagSeq || passHLT_CaloJet40) 
	     && (passHLT_L1HTT150_BtagSeq || passHLT_L1HTT150) )
	   {
	     h_mjj_HLTpass_CaloJet40_L1HTT150 -> Fill(getVariableValue("mjjCorr"));  
	     h_mjj_HLTpass_CaloJet40_L1HTT150_1GeVbin -> Fill(getVariableValue("mjjCorr"));  
	   }

	 if( (passHLT_L1HTT150_BtagSeq || passHLT_L1HTT150) )
	   {
	     h_mjj_HLTpass_L1HTT150 -> Fill(getVariableValue("mjjCorr"));  
	     h_mjj_HLTpass_L1HTT150_1GeVbin -> Fill(getVariableValue("mjjCorr"));  

	     if( run >= 257968 && run <=258432)// 415.879 pb-1
	       {
		 h_mjj_HLTpass_L1HTT150_p1 -> Fill(getVariableValue("mjjCorr"));  
		 h_mjj_HLTpass_L1HTT150_1GeVbin_p1 -> Fill(getVariableValue("mjjCorr"));  
	       }

	     if( run >= 258434 && run <=258745)// 460.903 pb-1
	       {
		 h_mjj_HLTpass_L1HTT150_p2 -> Fill(getVariableValue("mjjCorr"));  
		 h_mjj_HLTpass_L1HTT150_1GeVbin_p2 -> Fill(getVariableValue("mjjCorr"));  
	       }

	     if( run >= 258749 && run <=260425)// 457.316 pb-1
	       {
		 h_mjj_HLTpass_L1HTT150_p3 -> Fill(getVariableValue("mjjCorr"));  
		 h_mjj_HLTpass_L1HTT150_1GeVbin_p3 -> Fill(getVariableValue("mjjCorr"));  
	       }

	     if( run >= 260426 && run <=260627)// 470.504 pb-1
	       {
		 h_mjj_HLTpass_L1HTT150_p4 -> Fill(getVariableValue("mjjCorr"));  
		 h_mjj_HLTpass_L1HTT150_1GeVbin_p4 -> Fill(getVariableValue("mjjCorr"));  
	       }
	   }	
 
	 if( (passHLT_L1HTT150_BtagSeq || passHLT_L1HTT150) 
	     && (passHLT_HT450_BtagSeq || passHLT_HT450) )
	   {
	     h_mjj_HLTpass_L1HTT150_HT450 -> Fill(getVariableValue("mjjCorr"));  
	     h_mjj_HLTpass_L1HTT150_HT450_1GeVbin -> Fill(getVariableValue("mjjCorr"));  

	     if( run >= 257968 && run <=258432)// 415.879 pb-1
	       {
		 h_mjj_HLTpass_L1HTT150_HT450_p1 -> Fill(getVariableValue("mjjCorr"));  
		 h_mjj_HLTpass_L1HTT150_HT450_1GeVbin_p1 -> Fill(getVariableValue("mjjCorr"));  
	       }

	     if( run >= 258434 && run <=258745)// 460.903 pb-1
	       {
		 h_mjj_HLTpass_L1HTT150_HT450_p2 -> Fill(getVariableValue("mjjCorr"));  
		 h_mjj_HLTpass_L1HTT150_HT450_1GeVbin_p2 -> Fill(getVariableValue("mjjCorr"));  
	       }

	     if( run >= 258749 && run <=260425)// 457.316 pb-1
	       {
		 h_mjj_HLTpass_L1HTT150_HT450_p3 -> Fill(getVariableValue("mjjCorr"));  
		 h_mjj_HLTpass_L1HTT150_HT450_1GeVbin_p3 -> Fill(getVariableValue("mjjCorr"));  
	       }

	     if( run >= 260426 && run <=260627)// 470.504 pb-1
	       {		 
		 h_mjj_HLTpass_L1HTT150_HT450_p4 -> Fill(getVariableValue("mjjCorr"));  
		 h_mjj_HLTpass_L1HTT150_HT450_1GeVbin_p4 -> Fill(getVariableValue("mjjCorr"));  
	       }

	   }

       }     

     /*
     // optional call to fill a skim with a subset of the variables defined in the cutFile (use flag SAVE)
     if( passedAllPreviousCuts("mjj") && passedCut("mjj") ) 
       {
	 fillReducedSkimTree();
       }
     */

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

   f_mjjCorr -> Write();
   h2_f_mjjCorr_vs_mjj -> Write();

   h_mjj_NoTrigger -> Write();
   h_mjj_HLTpass_ZeroBias -> Write();
   h_mjj_HLTpass_ZeroBias_L1HTT150 -> Write();
   h_mjj_HLTpass_L1HTT150 -> Write();
   h_mjj_HLTpass_L1HTT150_HT450 -> Write();
   h_mjj_HLTpass_L1HTT150_p1 -> Write();
   h_mjj_HLTpass_L1HTT150_HT450_p1 -> Write();
   h_mjj_HLTpass_L1HTT150_p2 -> Write();
   h_mjj_HLTpass_L1HTT150_HT450_p2 -> Write();
   h_mjj_HLTpass_L1HTT150_p3 -> Write();
   h_mjj_HLTpass_L1HTT150_HT450_p3 -> Write();
   h_mjj_HLTpass_L1HTT150_p4 -> Write();
   h_mjj_HLTpass_L1HTT150_HT450_p4 -> Write();

   h_mjj_NoTrigger_1GeVbin -> Write();
   h_mjj_HLTpass_ZeroBias_1GeVbin -> Write();
   h_mjj_HLTpass_ZeroBias_L1HTT150_1GeVbin -> Write();
   h_mjj_HLTpass_L1HTT150_1GeVbin -> Write();
   h_mjj_HLTpass_L1HTT150_HT450_1GeVbin -> Write();
   h_mjj_HLTpass_L1HTT150_1GeVbin_p1 -> Write();
   h_mjj_HLTpass_L1HTT150_HT450_1GeVbin_p1 -> Write();
   h_mjj_HLTpass_L1HTT150_1GeVbin_p2 -> Write();
   h_mjj_HLTpass_L1HTT150_HT450_1GeVbin_p2 -> Write();
   h_mjj_HLTpass_L1HTT150_1GeVbin_p3 -> Write();
   h_mjj_HLTpass_L1HTT150_HT450_1GeVbin_p3 -> Write();
   h_mjj_HLTpass_L1HTT150_1GeVbin_p4 -> Write();
   h_mjj_HLTpass_L1HTT150_HT450_1GeVbin_p4 -> Write();

   h_mjj_HLTpass_CaloJet40 -> Write();
   h_mjj_HLTpass_CaloJet40_1GeVbin -> Write();
   h_mjj_HLTpass_CaloJet40_L1HTT150 -> Write();
   h_mjj_HLTpass_CaloJet40_L1HTT150_1GeVbin -> Write();

   //h_mjj_fullSel_varBin->Write();
   //h_mjj_fullSel_fixBin->Write();

   // //pT of both jets, to be built using the histograms produced automatically by baseClass
   // TH1F * h_pTJets = new TH1F ("h_pTJets","", getHistoNBins("pT1stJet"), getHistoMin("pT1stJet"), getHistoMax("pT1stJet"));
   // h_pTJets->Add( & getHisto_noCuts_or_skim("pT1stJet") ); // all histos can be retrieved, see other getHisto_xxxx methods in baseClass.h
   // h_pTJets->Add( & getHisto_noCuts_or_skim("pT2ndJet") );
   // //one could also do:  *h_pTJets = getHisto_noCuts_or_skim("pT1stJet") + getHisto_noCuts_or_skim("pT2ndJet");
   // h_pTJets->Write();
   // //one could also do:   const TH1F& h = getHisto_noCuts_or_skim// and use h

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
   }
