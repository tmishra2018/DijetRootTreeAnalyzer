//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jan 30 19:27:01 2016 by ROOT version 6.02/05
// from TChain rootTupleTree/tree/
//////////////////////////////////////////////////////////

#ifndef rootNtupleClass_h
#define rootNtupleClass_h

//// Lines added by make_rootNtupleClass.sh - BEGIN 
#include <vector> 
using namespace std; 
//// Lines added by make_rootNtupleClass.sh - END 

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class rootNtupleClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        CosThetaStarAK4;
   Double_t        CosThetaStarAK4reco;
   Double_t        CosThetaStarWJ;
   Double_t        CosThetaStarWJreco;
   Double_t        Dijet_MassAK4;
   Double_t        Dijet_MassAK4reco;
   Double_t        IdTight_j1;
   Double_t        IdTight_j2;
   Double_t        IdTight_recoj1;
   Double_t        IdTight_recoj2;
   Double_t        Nak4;
   Double_t        Nak4reco;
   Double_t        PassJSON;
   Double_t        chargedElectromFrac_j1;
   Double_t        chargedElectromFrac_j2;
   Double_t        chargedElectromFrac_recoj1;
   Double_t        chargedElectromFrac_recoj2;
   Double_t        chargedHadEnFrac_j1;
   Double_t        chargedHadEnFrac_j2;
   Double_t        chargedHadEnFrac_recoj1;
   Double_t        chargedHadEnFrac_recoj2;
   Double_t        chargedMult_j1;
   Double_t        chargedMult_j2;
   Double_t        chargedMult_recoj1;
   Double_t        chargedMult_recoj2;
   Double_t        deltaETAjj;
   Double_t        deltaETAjjAK4;
   Double_t        deltaETAjjAK4reco;
   Double_t        deltaETAjjreco;
   Double_t        deltaPHIjj;
   Double_t        deltaPHIjjAK4;
   Double_t        deltaPHIjjAK4reco;
   Double_t        deltaPHIjjreco;
   Double_t        eleEnFract_j1;
   Double_t        eleEnFract_j2;
   Double_t        eleEnFract_recoj1;
   Double_t        eleEnFract_recoj2;
   Double_t        etaAK4_j1;
   Double_t        etaAK4_j2;
   Double_t        etaAK4_recoj1;
   Double_t        etaAK4_recoj2;
   Double_t        etaWJ_j1;
   Double_t        etaWJ_j2;
   Double_t        etaWJ_recoj1;
   Double_t        etaWJ_recoj2;
   Double_t        event;
   Double_t        htAK4;
   Double_t        htAK4reco;
   Double_t        isData;
   Double_t        jetCSVAK4_j1;
   Double_t        jetCSVAK4_j2;
   Double_t        jetCSVAK4_recoj1;
   Double_t        jetCSVAK4_recoj2;
   Double_t        jetJecAK4_j1;
   Double_t        jetJecAK4_j2;
   Double_t        jetJecAK4_recoj1;
   Double_t        jetJecAK4_recoj2;
   Double_t        jetJecUncAK4_j1;
   Double_t        jetJecUncAK4_j2;
   Double_t        jetJecUncAK4_recoj1;
   Double_t        jetJecUncAK4_recoj2;
   Double_t        lumi;
   Double_t        massWJ_j1;
   Double_t        massWJ_j2;
   Double_t        massWJ_recoj1;
   Double_t        massWJ_recoj2;
   Double_t        met;
   Double_t        metSig;
   Double_t        metreco;
   Double_t        metrecoSig;
   Double_t        mhtAK4;
   Double_t        mhtAK4Sig;
   Double_t        mhtAK4reco;
   Double_t        mhtAK4recoSig;
   Double_t        mjj;
   Double_t        mjj_shiftJEC;
   Double_t        mjjreco;
   Double_t        mjjreco_shiftJEC;
   Double_t        muEnFract_j1;
   Double_t        muEnFract_j2;
   Double_t        muEnFract_recoj1;
   Double_t        muEnFract_recoj2;
   Double_t        nJet;
   Double_t        nJetreco;
   Double_t        nVtx;
   Double_t        nVtxreco;
   Double_t        neutrElectromFrac_j1;
   Double_t        neutrElectromFrac_j2;
   Double_t        neutrElectromFrac_recoj1;
   Double_t        neutrElectromFrac_recoj2;
   Double_t        neutrHadEnFrac_j1;
   Double_t        neutrHadEnFrac_j2;
   Double_t        neutrHadEnFrac_recoj1;
   Double_t        neutrHadEnFrac_recoj2;
   Double_t        neutrMult_j1;
   Double_t        neutrMult_j2;
   Double_t        neutrMult_recoj1;
   Double_t        neutrMult_recoj2;
   Double_t        offMet;
   Double_t        offMetSig;
   Double_t        pTAK4_j1;
   Double_t        pTAK4_j2;
   Double_t        pTAK4_recoj1;
   Double_t        pTAK4_recoj2;
   Double_t        pTWJ_j1;
   Double_t        pTWJ_j2;
   Double_t        pTWJ_recoj1;
   Double_t        pTWJ_recoj2;
   Double_t        passHLT_CaloJet40;
   Double_t        passHLT_CaloJet40_BtagSeq;
   Double_t        passHLT_HT450;
   Double_t        passHLT_HT450_BtagSeq;
   Double_t        passHLT_L1DoubleMu;
   Double_t        passHLT_L1DoubleMu_BtagSeq;
   Double_t        passHLT_L1HTT150;
   Double_t        passHLT_L1HTT150_BtagSeq;
   Double_t        passHLT_PFHT650MJJ900;
   Double_t        passHLT_PFHT650MJJ950;
   Double_t        passHLT_PFHT800;
   Double_t        passHLT_ZeroBias;
   Double_t        passHLT_ZeroBias_BtagSeq;
   Double_t        phiAK4_j1;
   Double_t        phiAK4_j2;
   Double_t        phiAK4_recoj1;
   Double_t        phiAK4_recoj2;
   Double_t        phiWJ_j1;
   Double_t        phiWJ_j2;
   Double_t        phiWJ_recoj1;
   Double_t        phiWJ_recoj2;
   Double_t        photonEnFrac_j1;
   Double_t        photonEnFrac_j2;
   Double_t        photonEnFrac_recoj1;
   Double_t        photonEnFrac_recoj2;
   Double_t        photonMult_j1;
   Double_t        photonMult_j2;
   Double_t        photonMult_recoj1;
   Double_t        photonMult_recoj2;
   Double_t        run;

   // List of branches
   TBranch        *b_CosThetaStarAK4;   //!
   TBranch        *b_CosThetaStarAK4reco;   //!
   TBranch        *b_CosThetaStarWJ;   //!
   TBranch        *b_CosThetaStarWJreco;   //!
   TBranch        *b_Dijet_MassAK4;   //!
   TBranch        *b_Dijet_MassAK4reco;   //!
   TBranch        *b_IdTight_j1;   //!
   TBranch        *b_IdTight_j2;   //!
   TBranch        *b_IdTight_recoj1;   //!
   TBranch        *b_IdTight_recoj2;   //!
   TBranch        *b_Nak4;   //!
   TBranch        *b_Nak4reco;   //!
   TBranch        *b_PassJSON;   //!
   TBranch        *b_chargedElectromFrac_j1;   //!
   TBranch        *b_chargedElectromFrac_j2;   //!
   TBranch        *b_chargedElectromFrac_recoj1;   //!
   TBranch        *b_chargedElectromFrac_recoj2;   //!
   TBranch        *b_chargedHadEnFrac_j1;   //!
   TBranch        *b_chargedHadEnFrac_j2;   //!
   TBranch        *b_chargedHadEnFrac_recoj1;   //!
   TBranch        *b_chargedHadEnFrac_recoj2;   //!
   TBranch        *b_chargedMult_j1;   //!
   TBranch        *b_chargedMult_j2;   //!
   TBranch        *b_chargedMult_recoj1;   //!
   TBranch        *b_chargedMult_recoj2;   //!
   TBranch        *b_deltaETAjj;   //!
   TBranch        *b_deltaETAjjAK4;   //!
   TBranch        *b_deltaETAjjAK4reco;   //!
   TBranch        *b_deltaETAjjreco;   //!
   TBranch        *b_deltaPHIjj;   //!
   TBranch        *b_deltaPHIjjAK4;   //!
   TBranch        *b_deltaPHIjjAK4reco;   //!
   TBranch        *b_deltaPHIjjreco;   //!
   TBranch        *b_eleEnFract_j1;   //!
   TBranch        *b_eleEnFract_j2;   //!
   TBranch        *b_eleEnFract_recoj1;   //!
   TBranch        *b_eleEnFract_recoj2;   //!
   TBranch        *b_etaAK4_j1;   //!
   TBranch        *b_etaAK4_j2;   //!
   TBranch        *b_etaAK4_recoj1;   //!
   TBranch        *b_etaAK4_recoj2;   //!
   TBranch        *b_etaWJ_j1;   //!
   TBranch        *b_etaWJ_j2;   //!
   TBranch        *b_etaWJ_recoj1;   //!
   TBranch        *b_etaWJ_recoj2;   //!
   TBranch        *b_event;   //!
   TBranch        *b_htAK4;   //!
   TBranch        *b_htAK4reco;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_jetCSVAK4_j1;   //!
   TBranch        *b_jetCSVAK4_j2;   //!
   TBranch        *b_jetCSVAK4_recoj1;   //!
   TBranch        *b_jetCSVAK4_recoj2;   //!
   TBranch        *b_jetJecAK4_j1;   //!
   TBranch        *b_jetJecAK4_j2;   //!
   TBranch        *b_jetJecAK4_recoj1;   //!
   TBranch        *b_jetJecAK4_recoj2;   //!
   TBranch        *b_jetJecUncAK4_j1;   //!
   TBranch        *b_jetJecUncAK4_j2;   //!
   TBranch        *b_jetJecUncAK4_recoj1;   //!
   TBranch        *b_jetJecUncAK4_recoj2;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_massWJ_j1;   //!
   TBranch        *b_massWJ_j2;   //!
   TBranch        *b_massWJ_recoj1;   //!
   TBranch        *b_massWJ_recoj2;   //!
   TBranch        *b_met;   //!
   TBranch        *b_metSig;   //!
   TBranch        *b_metreco;   //!
   TBranch        *b_metrecoSig;   //!
   TBranch        *b_mhtAK4;   //!
   TBranch        *b_mhtAK4Sig;   //!
   TBranch        *b_mhtAK4reco;   //!
   TBranch        *b_mhtAK4recoSig;   //!
   TBranch        *b_mjj;   //!
   TBranch        *b_mjj_shiftJEC;   //!
   TBranch        *b_mjjreco;   //!
   TBranch        *b_mjjreco_shiftJEC;   //!
   TBranch        *b_muEnFract_j1;   //!
   TBranch        *b_muEnFract_j2;   //!
   TBranch        *b_muEnFract_recoj1;   //!
   TBranch        *b_muEnFract_recoj2;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_nJetreco;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_nVtxreco;   //!
   TBranch        *b_neutrElectromFrac_j1;   //!
   TBranch        *b_neutrElectromFrac_j2;   //!
   TBranch        *b_neutrElectromFrac_recoj1;   //!
   TBranch        *b_neutrElectromFrac_recoj2;   //!
   TBranch        *b_neutrHadEnFrac_j1;   //!
   TBranch        *b_neutrHadEnFrac_j2;   //!
   TBranch        *b_neutrHadEnFrac_recoj1;   //!
   TBranch        *b_neutrHadEnFrac_recoj2;   //!
   TBranch        *b_neutrMult_j1;   //!
   TBranch        *b_neutrMult_j2;   //!
   TBranch        *b_neutrMult_recoj1;   //!
   TBranch        *b_neutrMult_recoj2;   //!
   TBranch        *b_offMet;   //!
   TBranch        *b_offMetSig;   //!
   TBranch        *b_pTAK4_j1;   //!
   TBranch        *b_pTAK4_j2;   //!
   TBranch        *b_pTAK4_recoj1;   //!
   TBranch        *b_pTAK4_recoj2;   //!
   TBranch        *b_pTWJ_j1;   //!
   TBranch        *b_pTWJ_j2;   //!
   TBranch        *b_pTWJ_recoj1;   //!
   TBranch        *b_pTWJ_recoj2;   //!
   TBranch        *b_passHLT_CaloJet40;   //!
   TBranch        *b_passHLT_CaloJet40_BtagSeq;   //!
   TBranch        *b_passHLT_HT450;   //!
   TBranch        *b_passHLT_HT450_BtagSeq;   //!
   TBranch        *b_passHLT_L1DoubleMu;   //!
   TBranch        *b_passHLT_L1DoubleMu_BtagSeq;   //!
   TBranch        *b_passHLT_L1HTT150;   //!
   TBranch        *b_passHLT_L1HTT150_BtagSeq;   //!
   TBranch        *b_passHLT_PFHT650MJJ900;   //!
   TBranch        *b_passHLT_PFHT650MJJ950;   //!
   TBranch        *b_passHLT_PFHT800;   //!
   TBranch        *b_passHLT_ZeroBias;   //!
   TBranch        *b_passHLT_ZeroBias_BtagSeq;   //!
   TBranch        *b_phiAK4_j1;   //!
   TBranch        *b_phiAK4_j2;   //!
   TBranch        *b_phiAK4_recoj1;   //!
   TBranch        *b_phiAK4_recoj2;   //!
   TBranch        *b_phiWJ_j1;   //!
   TBranch        *b_phiWJ_j2;   //!
   TBranch        *b_phiWJ_recoj1;   //!
   TBranch        *b_phiWJ_recoj2;   //!
   TBranch        *b_photonEnFrac_j1;   //!
   TBranch        *b_photonEnFrac_j2;   //!
   TBranch        *b_photonEnFrac_recoj1;   //!
   TBranch        *b_photonEnFrac_recoj2;   //!
   TBranch        *b_photonMult_j1;   //!
   TBranch        *b_photonMult_j2;   //!
   TBranch        *b_photonMult_recoj1;   //!
   TBranch        *b_photonMult_recoj2;   //!
   TBranch        *b_run;   //!

   rootNtupleClass(TTree *tree=0);
   virtual ~rootNtupleClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef rootNtupleClass_cxx
rootNtupleClass::rootNtupleClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("rootTupleTree/tree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("rootTupleTree/tree","");
      chain->Add("output/rootTest_reduced_skim.root/rootTupleTree/tree");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

rootNtupleClass::~rootNtupleClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t rootNtupleClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t rootNtupleClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void rootNtupleClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("CosThetaStarAK4", &CosThetaStarAK4, &b_CosThetaStarAK4);
   fChain->SetBranchAddress("CosThetaStarAK4reco", &CosThetaStarAK4reco, &b_CosThetaStarAK4reco);
   fChain->SetBranchAddress("CosThetaStarWJ", &CosThetaStarWJ, &b_CosThetaStarWJ);
   fChain->SetBranchAddress("CosThetaStarWJreco", &CosThetaStarWJreco, &b_CosThetaStarWJreco);
   fChain->SetBranchAddress("Dijet_MassAK4", &Dijet_MassAK4, &b_Dijet_MassAK4);
   fChain->SetBranchAddress("Dijet_MassAK4reco", &Dijet_MassAK4reco, &b_Dijet_MassAK4reco);
   fChain->SetBranchAddress("IdTight_j1", &IdTight_j1, &b_IdTight_j1);
   fChain->SetBranchAddress("IdTight_j2", &IdTight_j2, &b_IdTight_j2);
   fChain->SetBranchAddress("IdTight_recoj1", &IdTight_recoj1, &b_IdTight_recoj1);
   fChain->SetBranchAddress("IdTight_recoj2", &IdTight_recoj2, &b_IdTight_recoj2);
   fChain->SetBranchAddress("Nak4", &Nak4, &b_Nak4);
   fChain->SetBranchAddress("Nak4reco", &Nak4reco, &b_Nak4reco);
   fChain->SetBranchAddress("PassJSON", &PassJSON, &b_PassJSON);
   fChain->SetBranchAddress("chargedElectromFrac_j1", &chargedElectromFrac_j1, &b_chargedElectromFrac_j1);
   fChain->SetBranchAddress("chargedElectromFrac_j2", &chargedElectromFrac_j2, &b_chargedElectromFrac_j2);
   fChain->SetBranchAddress("chargedElectromFrac_recoj1", &chargedElectromFrac_recoj1, &b_chargedElectromFrac_recoj1);
   fChain->SetBranchAddress("chargedElectromFrac_recoj2", &chargedElectromFrac_recoj2, &b_chargedElectromFrac_recoj2);
   fChain->SetBranchAddress("chargedHadEnFrac_j1", &chargedHadEnFrac_j1, &b_chargedHadEnFrac_j1);
   fChain->SetBranchAddress("chargedHadEnFrac_j2", &chargedHadEnFrac_j2, &b_chargedHadEnFrac_j2);
   fChain->SetBranchAddress("chargedHadEnFrac_recoj1", &chargedHadEnFrac_recoj1, &b_chargedHadEnFrac_recoj1);
   fChain->SetBranchAddress("chargedHadEnFrac_recoj2", &chargedHadEnFrac_recoj2, &b_chargedHadEnFrac_recoj2);
   fChain->SetBranchAddress("chargedMult_j1", &chargedMult_j1, &b_chargedMult_j1);
   fChain->SetBranchAddress("chargedMult_j2", &chargedMult_j2, &b_chargedMult_j2);
   fChain->SetBranchAddress("chargedMult_recoj1", &chargedMult_recoj1, &b_chargedMult_recoj1);
   fChain->SetBranchAddress("chargedMult_recoj2", &chargedMult_recoj2, &b_chargedMult_recoj2);
   fChain->SetBranchAddress("deltaETAjj", &deltaETAjj, &b_deltaETAjj);
   fChain->SetBranchAddress("deltaETAjjAK4", &deltaETAjjAK4, &b_deltaETAjjAK4);
   fChain->SetBranchAddress("deltaETAjjAK4reco", &deltaETAjjAK4reco, &b_deltaETAjjAK4reco);
   fChain->SetBranchAddress("deltaETAjjreco", &deltaETAjjreco, &b_deltaETAjjreco);
   fChain->SetBranchAddress("deltaPHIjj", &deltaPHIjj, &b_deltaPHIjj);
   fChain->SetBranchAddress("deltaPHIjjAK4", &deltaPHIjjAK4, &b_deltaPHIjjAK4);
   fChain->SetBranchAddress("deltaPHIjjAK4reco", &deltaPHIjjAK4reco, &b_deltaPHIjjAK4reco);
   fChain->SetBranchAddress("deltaPHIjjreco", &deltaPHIjjreco, &b_deltaPHIjjreco);
   fChain->SetBranchAddress("eleEnFract_j1", &eleEnFract_j1, &b_eleEnFract_j1);
   fChain->SetBranchAddress("eleEnFract_j2", &eleEnFract_j2, &b_eleEnFract_j2);
   fChain->SetBranchAddress("eleEnFract_recoj1", &eleEnFract_recoj1, &b_eleEnFract_recoj1);
   fChain->SetBranchAddress("eleEnFract_recoj2", &eleEnFract_recoj2, &b_eleEnFract_recoj2);
   fChain->SetBranchAddress("etaAK4_j1", &etaAK4_j1, &b_etaAK4_j1);
   fChain->SetBranchAddress("etaAK4_j2", &etaAK4_j2, &b_etaAK4_j2);
   fChain->SetBranchAddress("etaAK4_recoj1", &etaAK4_recoj1, &b_etaAK4_recoj1);
   fChain->SetBranchAddress("etaAK4_recoj2", &etaAK4_recoj2, &b_etaAK4_recoj2);
   fChain->SetBranchAddress("etaWJ_j1", &etaWJ_j1, &b_etaWJ_j1);
   fChain->SetBranchAddress("etaWJ_j2", &etaWJ_j2, &b_etaWJ_j2);
   fChain->SetBranchAddress("etaWJ_recoj1", &etaWJ_recoj1, &b_etaWJ_recoj1);
   fChain->SetBranchAddress("etaWJ_recoj2", &etaWJ_recoj2, &b_etaWJ_recoj2);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("htAK4", &htAK4, &b_htAK4);
   fChain->SetBranchAddress("htAK4reco", &htAK4reco, &b_htAK4reco);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("jetCSVAK4_j1", &jetCSVAK4_j1, &b_jetCSVAK4_j1);
   fChain->SetBranchAddress("jetCSVAK4_j2", &jetCSVAK4_j2, &b_jetCSVAK4_j2);
   fChain->SetBranchAddress("jetCSVAK4_recoj1", &jetCSVAK4_recoj1, &b_jetCSVAK4_recoj1);
   fChain->SetBranchAddress("jetCSVAK4_recoj2", &jetCSVAK4_recoj2, &b_jetCSVAK4_recoj2);
   fChain->SetBranchAddress("jetJecAK4_j1", &jetJecAK4_j1, &b_jetJecAK4_j1);
   fChain->SetBranchAddress("jetJecAK4_j2", &jetJecAK4_j2, &b_jetJecAK4_j2);
   fChain->SetBranchAddress("jetJecAK4_recoj1", &jetJecAK4_recoj1, &b_jetJecAK4_recoj1);
   fChain->SetBranchAddress("jetJecAK4_recoj2", &jetJecAK4_recoj2, &b_jetJecAK4_recoj2);
   fChain->SetBranchAddress("jetJecUncAK4_j1", &jetJecUncAK4_j1, &b_jetJecUncAK4_j1);
   fChain->SetBranchAddress("jetJecUncAK4_j2", &jetJecUncAK4_j2, &b_jetJecUncAK4_j2);
   fChain->SetBranchAddress("jetJecUncAK4_recoj1", &jetJecUncAK4_recoj1, &b_jetJecUncAK4_recoj1);
   fChain->SetBranchAddress("jetJecUncAK4_recoj2", &jetJecUncAK4_recoj2, &b_jetJecUncAK4_recoj2);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("massWJ_j1", &massWJ_j1, &b_massWJ_j1);
   fChain->SetBranchAddress("massWJ_j2", &massWJ_j2, &b_massWJ_j2);
   fChain->SetBranchAddress("massWJ_recoj1", &massWJ_recoj1, &b_massWJ_recoj1);
   fChain->SetBranchAddress("massWJ_recoj2", &massWJ_recoj2, &b_massWJ_recoj2);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("metSig", &metSig, &b_metSig);
   fChain->SetBranchAddress("metreco", &metreco, &b_metreco);
   fChain->SetBranchAddress("metrecoSig", &metrecoSig, &b_metrecoSig);
   fChain->SetBranchAddress("mhtAK4", &mhtAK4, &b_mhtAK4);
   fChain->SetBranchAddress("mhtAK4Sig", &mhtAK4Sig, &b_mhtAK4Sig);
   fChain->SetBranchAddress("mhtAK4reco", &mhtAK4reco, &b_mhtAK4reco);
   fChain->SetBranchAddress("mhtAK4recoSig", &mhtAK4recoSig, &b_mhtAK4recoSig);
   fChain->SetBranchAddress("mjj", &mjj, &b_mjj);
   fChain->SetBranchAddress("mjj_shiftJEC", &mjj_shiftJEC, &b_mjj_shiftJEC);
   fChain->SetBranchAddress("mjjreco", &mjjreco, &b_mjjreco);
   fChain->SetBranchAddress("mjjreco_shiftJEC", &mjjreco_shiftJEC, &b_mjjreco_shiftJEC);
   fChain->SetBranchAddress("muEnFract_j1", &muEnFract_j1, &b_muEnFract_j1);
   fChain->SetBranchAddress("muEnFract_j2", &muEnFract_j2, &b_muEnFract_j2);
   fChain->SetBranchAddress("muEnFract_recoj1", &muEnFract_recoj1, &b_muEnFract_recoj1);
   fChain->SetBranchAddress("muEnFract_recoj2", &muEnFract_recoj2, &b_muEnFract_recoj2);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("nJetreco", &nJetreco, &b_nJetreco);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("nVtxreco", &nVtxreco, &b_nVtxreco);
   fChain->SetBranchAddress("neutrElectromFrac_j1", &neutrElectromFrac_j1, &b_neutrElectromFrac_j1);
   fChain->SetBranchAddress("neutrElectromFrac_j2", &neutrElectromFrac_j2, &b_neutrElectromFrac_j2);
   fChain->SetBranchAddress("neutrElectromFrac_recoj1", &neutrElectromFrac_recoj1, &b_neutrElectromFrac_recoj1);
   fChain->SetBranchAddress("neutrElectromFrac_recoj2", &neutrElectromFrac_recoj2, &b_neutrElectromFrac_recoj2);
   fChain->SetBranchAddress("neutrHadEnFrac_j1", &neutrHadEnFrac_j1, &b_neutrHadEnFrac_j1);
   fChain->SetBranchAddress("neutrHadEnFrac_j2", &neutrHadEnFrac_j2, &b_neutrHadEnFrac_j2);
   fChain->SetBranchAddress("neutrHadEnFrac_recoj1", &neutrHadEnFrac_recoj1, &b_neutrHadEnFrac_recoj1);
   fChain->SetBranchAddress("neutrHadEnFrac_recoj2", &neutrHadEnFrac_recoj2, &b_neutrHadEnFrac_recoj2);
   fChain->SetBranchAddress("neutrMult_j1", &neutrMult_j1, &b_neutrMult_j1);
   fChain->SetBranchAddress("neutrMult_j2", &neutrMult_j2, &b_neutrMult_j2);
   fChain->SetBranchAddress("neutrMult_recoj1", &neutrMult_recoj1, &b_neutrMult_recoj1);
   fChain->SetBranchAddress("neutrMult_recoj2", &neutrMult_recoj2, &b_neutrMult_recoj2);
   fChain->SetBranchAddress("offMet", &offMet, &b_offMet);
   fChain->SetBranchAddress("offMetSig", &offMetSig, &b_offMetSig);
   fChain->SetBranchAddress("pTAK4_j1", &pTAK4_j1, &b_pTAK4_j1);
   fChain->SetBranchAddress("pTAK4_j2", &pTAK4_j2, &b_pTAK4_j2);
   fChain->SetBranchAddress("pTAK4_recoj1", &pTAK4_recoj1, &b_pTAK4_recoj1);
   fChain->SetBranchAddress("pTAK4_recoj2", &pTAK4_recoj2, &b_pTAK4_recoj2);
   fChain->SetBranchAddress("pTWJ_j1", &pTWJ_j1, &b_pTWJ_j1);
   fChain->SetBranchAddress("pTWJ_j2", &pTWJ_j2, &b_pTWJ_j2);
   fChain->SetBranchAddress("pTWJ_recoj1", &pTWJ_recoj1, &b_pTWJ_recoj1);
   fChain->SetBranchAddress("pTWJ_recoj2", &pTWJ_recoj2, &b_pTWJ_recoj2);
   fChain->SetBranchAddress("passHLT_CaloJet40", &passHLT_CaloJet40, &b_passHLT_CaloJet40);
   fChain->SetBranchAddress("passHLT_CaloJet40_BtagSeq", &passHLT_CaloJet40_BtagSeq, &b_passHLT_CaloJet40_BtagSeq);
   fChain->SetBranchAddress("passHLT_HT450", &passHLT_HT450, &b_passHLT_HT450);
   fChain->SetBranchAddress("passHLT_HT450_BtagSeq", &passHLT_HT450_BtagSeq, &b_passHLT_HT450_BtagSeq);
   fChain->SetBranchAddress("passHLT_L1DoubleMu", &passHLT_L1DoubleMu, &b_passHLT_L1DoubleMu);
   fChain->SetBranchAddress("passHLT_L1DoubleMu_BtagSeq", &passHLT_L1DoubleMu_BtagSeq, &b_passHLT_L1DoubleMu_BtagSeq);
   fChain->SetBranchAddress("passHLT_L1HTT150", &passHLT_L1HTT150, &b_passHLT_L1HTT150);
   fChain->SetBranchAddress("passHLT_L1HTT150_BtagSeq", &passHLT_L1HTT150_BtagSeq, &b_passHLT_L1HTT150_BtagSeq);
   fChain->SetBranchAddress("passHLT_PFHT650MJJ900", &passHLT_PFHT650MJJ900, &b_passHLT_PFHT650MJJ900);
   fChain->SetBranchAddress("passHLT_PFHT650MJJ950", &passHLT_PFHT650MJJ950, &b_passHLT_PFHT650MJJ950);
   fChain->SetBranchAddress("passHLT_PFHT800", &passHLT_PFHT800, &b_passHLT_PFHT800);
   fChain->SetBranchAddress("passHLT_ZeroBias", &passHLT_ZeroBias, &b_passHLT_ZeroBias);
   fChain->SetBranchAddress("passHLT_ZeroBias_BtagSeq", &passHLT_ZeroBias_BtagSeq, &b_passHLT_ZeroBias_BtagSeq);
   fChain->SetBranchAddress("phiAK4_j1", &phiAK4_j1, &b_phiAK4_j1);
   fChain->SetBranchAddress("phiAK4_j2", &phiAK4_j2, &b_phiAK4_j2);
   fChain->SetBranchAddress("phiAK4_recoj1", &phiAK4_recoj1, &b_phiAK4_recoj1);
   fChain->SetBranchAddress("phiAK4_recoj2", &phiAK4_recoj2, &b_phiAK4_recoj2);
   fChain->SetBranchAddress("phiWJ_j1", &phiWJ_j1, &b_phiWJ_j1);
   fChain->SetBranchAddress("phiWJ_j2", &phiWJ_j2, &b_phiWJ_j2);
   fChain->SetBranchAddress("phiWJ_recoj1", &phiWJ_recoj1, &b_phiWJ_recoj1);
   fChain->SetBranchAddress("phiWJ_recoj2", &phiWJ_recoj2, &b_phiWJ_recoj2);
   fChain->SetBranchAddress("photonEnFrac_j1", &photonEnFrac_j1, &b_photonEnFrac_j1);
   fChain->SetBranchAddress("photonEnFrac_j2", &photonEnFrac_j2, &b_photonEnFrac_j2);
   fChain->SetBranchAddress("photonEnFrac_recoj1", &photonEnFrac_recoj1, &b_photonEnFrac_recoj1);
   fChain->SetBranchAddress("photonEnFrac_recoj2", &photonEnFrac_recoj2, &b_photonEnFrac_recoj2);
   fChain->SetBranchAddress("photonMult_j1", &photonMult_j1, &b_photonMult_j1);
   fChain->SetBranchAddress("photonMult_j2", &photonMult_j2, &b_photonMult_j2);
   fChain->SetBranchAddress("photonMult_recoj1", &photonMult_recoj1, &b_photonMult_recoj1);
   fChain->SetBranchAddress("photonMult_recoj2", &photonMult_recoj2, &b_photonMult_recoj2);
   fChain->SetBranchAddress("run", &run, &b_run);
   Notify();
}

Bool_t rootNtupleClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void rootNtupleClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t rootNtupleClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef rootNtupleClass_cxx
