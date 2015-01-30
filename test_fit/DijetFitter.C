#include "string"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TChain.h"
#include "TTree.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "RooCurve.h" 
#include "TLatex.h"


#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooBinning.h"
#include "RooWorkspace.h"
#include "RooHistFunc.h"
#include "RooPolynomial.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooLandau.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "RooProdPdf.h"
#include "RooExtendPdf.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooRealConstant.h"
#include "RooCategory.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooProduct.h"
#include "RooMinimizer.h"
#include "RooGenericPdf.h"
#include "RooStats/RooStatsUtils.h"
#include "RooStats/HLFactory.h"
#include "RooRandom.h"
#include "RooHist.h"
#include "RooHistPdf.h"

#include "DijetFitter.h"
#include "GetMinMAxFit.h"
//#include "Graphics.h"
#include "WriteWorkspace.h"
#include "BkgModelFunctions.h"

//#include "CMS_lumi.C"
//#include "SigModelFunctions.h"
//#include "WriteDatacard.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

//void AddSigData(RooWorkspace*, Float_t);
void AddBkgData(RooWorkspace*, Float_t);

//void SetConstantParams(const RooArgSet* params);
//void SigModelResponseFcnFit(RooWorkspace*);
//void SigModelFitConvBW(RooWorkspace*, Float_t, Double_t);

//RooFitResult* BkgModelFitExpPARFunc(RooWorkspace*, Bool_t, Float_t, Int_t, bool); 
RooFitResult*  BkgModelFitDiJetFunc(RooWorkspace*, Bool_t, Float_t, Int_t, bool);
//RooFitResult* BkgModelFitDiJetEXPFunc(RooWorkspace*, Bool_t, Float_t, Int_t, bool);
//RooFitResult* BkgModelFitDiJetEXPOLFunc(RooWorkspace*, Bool_t, Float_t, Int_t, bool);
//RooFitResult* BkgModelFitExpolFunc(RooWorkspace*, Bool_t, Float_t, bool);

RooArgSet* defineVariables(); 

// RooHistFunc* getRooHistFunc(int, RooRealVar* );
// void makeRooHistPlot(RooHistFunc*,RooHistFunc* ,TF1*, TF1*, RooRealVar*);
// RooHistFunc* getRooHistFuncFitMAX(int, RooRealVar*);
// RooHistFunc* getRooHistFuncFitMIN(int, RooRealVar*);
// TF1* getFuncFitMAX(int, RooRealVar* );
// TF1* getFuncFitMIN(int, RooRealVar* );

RooDataHist* binnedDataset(RooWorkspace*);

/////////////////////////////////////////////////////////////
//
//                     / _(_) |      
//  _ __ _   _ _ __   | |_ _| |_ ___ 
// | '__| | | | '_ \  |  _| | __/ __|
// | |  | |_| | | | | | | | | |_\__\
// |_|   \__,_|_| |_| |_| |_|\__|___/
//
///////////////////////////////////////////////////////////                 

void runfits(const Float_t mass=5000., Bool_t dobands = false,  std::string model = "") {
  
  //******************************************************************//
  //  Running mode  corresponds to the following cases
  //         - full run set:
  //         - create signal and background data sets 
  //         - make and fit signal and background  models 
  //         - write signal and background workspaces in root files
  //         - write data card
  //*******************************************************************//

  TString fileBaseName("DijetMass");    
  TString fileBkgName("DijetMass.inputbkg");
  
  //RooWorkspace* w = new RooWorkspace("w", "workspace");

  //initialization ??
  TString card_name("DijetMass_Bkg_13TeV_initialization.rs");
  HLFactory hlf("HLFactory", card_name, false);
  RooWorkspace* w = hlf.GetWs();

  RooFitResult* fitresults;

  // Luminosity:
  Float_t Lum = 1000.0;  
  RooRealVar lumi("lumi","lumi",Lum);
  w->import(lumi); 
  w->Print();
 

  //what?
  //compute roohistfunc with min and max for the fit range 
  //RooRealVar* var = new RooRealVar("var","var", 1100, 6000 );
  
  
  // RooHistFunc* rooFitMin = getRooHistFuncFitMIN(0, var);
  // RooHistFunc* rooFitMax = getRooHistFuncFitMAX(0, var);
  // TF1* fFitMin = getFuncFitMIN(0, var);
  // TF1* fFitMax = getFuncFitMAX(0, var);

  // std::cout<<fFitMin->Eval(mass)<<"    "<<fFitMax->Eval(mass)<<std::endl;
  // makeRooHistPlot(rooFitMin,rooFitMax,fFitMin,fFitMax,var); 
  
  // var->setVal(mass);
  // double newmin = fFitMin->Eval(mass);
  // double newmax = fFitMax->Eval(mass);

  // if(mass==5000){
  //   newmin=1100;
  //   newmax=6000;
  // }

  // MINmass=newmin;
  // MAXmass=newmax;
  // std::cout<<"  MIN MASS: "<<MINmass<<"   MAX MASS: "<<MAXmass<<std::endl;
  
  //giulia -  not for now
  //cout << endl; cout << "Now AddSigData" << endl;
  //AddSigData(w, mass, model);   
  //cout << endl; cout << "Now SigModelFit" << endl;  
  //SigModelResponseFcnFit(w, mass, model); 
  /* w->var("PhotonsMass")->setMin(MINmass);
     w->var("PhotonsMass")->setMax(MAXmass);*/

  cout << endl; cout << "Now AddBkgData" << endl;
  AddBkgData(w, mass);         
  RooDataHist* dataBinned = binnedDataset(w);
  w->import(*dataBinned);


  // //fit sig
  // cout << endl; cout << "Now SigModelFit" << endl;    
  // // SigModelFitConvBW(w, mass, width, model);      
  // //SigModelFitConvCPS(w, mass, width, model);      
  // //SigModelFitConvRelBW(w, mass, width, model); 
  

  cout << endl; cout << "Now BkgModelFit" << endl;    
  bool blind = false; 
  bool dobandsHere= false;

  //do the fit and write result on a file
  TFile* file_result;
  for(int c=0; c<NCAT;c++){
    fitresults = BkgModelFitDiJetFunc(w, dobandsHere, mass, c, blind);
    file_result = new TFile(TString::Format("roofit_result_cat%d.root",c),"recreate");
    file_result->cd();
    fitresults->Write();
    file_result->Close();
  }

  //write workspace
  cout << "now write workspace" << endl; 
  MakeBkgWS(w,fileBaseName,mass); 



  return;
}

/////////////////////////////////////////////////////////

//                | |                     / _(_) |      
//   ___ _ __   __| |  _ __ _   _ _ __   | |_ _| |_ ___ 
//  / _ \ '_ \ / _` | | '__| | | | '_ \  |  _| | __/ __|
// |  __/ | | | (_| | | |  | |_| | | | | | | | | |_\__		
//  \___|_| |_|\__,_| |_|   \__,_|_| |_| |_| |_|\__|___/
                                                     
////////////////////////////////////////////////////////////    

RooArgSet* defineVariables() {
  
  // define variables of the input ntuple 
  
  RooRealVar* mjj  = new RooRealVar("mjj", "M(WjWj)",1, 14000,"GeV");
  RooRealVar* pt_j1 = new RooRealVar("pt_j1", "pt(Wj1)",500,5000,"GeV");
  RooRealVar* pt_j2 = new RooRealVar("pt_j2", "pt(Wj2)",500,5000,"GeV");
  RooRealVar* eta_j1 = new RooRealVar("eta_j1", "eta(Wj1)",-3.,3.,"");
  RooRealVar* eta_j2 = new RooRealVar("eta_j2", "eta(Wj2)",-3.,3.,"");
  RooRealVar* deltaETAjj = new RooRealVar("deltaETAjj", "dEta(WjWj)",-3.,3.,"");
  
  //RooRealVar* evweight = new RooRealVar("evweight","weightings",0,1000,"");
  //RooRealVar* nvtx = new RooRealVar("nvtx", "nvtx", 0, 50, "");
  RooArgSet* ntplVars = new RooArgSet(*mjj, *pt_j1, *pt_j2, *eta_j1, *eta_j2, *deltaETAjj );
  
  return ntplVars;
}


//giulia - da sistemare, per ora non serve
// // Signal Data Set
// void AddSigData(RooWorkspace* w, Float_t mass, std::string model) {

//   Int_t ncat = NCAT;
//   TString inDir = "";

//   Float_t MASS(mass);


//   // Variables
//   RooArgSet* ntplVars = defineVariablesM150();
//   ntplVars->add(*w->var("PhotonsMass"));
//   (*w->var("PhotonsMass")).setRange(130,230);
//   ntplVars->add(*w->var("PhotonsMassTrue"));
//   (*w->var("PhotonsMassTrue")).setRange(130,230);
  
//   int iMass = abs(mass);      

  
//   //chain summing up all production modes
//   TChain* sigTree1  = new TChain();
//   /*  sigTree1->Add(TString::Format("histograms_CMS-HGG_19032013.root/ggh_m250_8TeV", iMass));
//   sigTree1->Add(TString::Format("histograms_CMS-HGG_19032013.root/vbf_m250_8TeV", iMass));
//   sigTree1->Add(TString::Format("histograms_CMS-HGG_19032013.root/wzh_m250_8TeV", iMass));
//   sigTree1->Add(TString::Format("histograms_CMS-HGG_19032013.root/tth_m250_8TeV", iMass));*/


//   //histograms_CMS-HGG_17042014_MC.root

//   sigTree1->Add("histograms_CMS-HGG_08052014_MC.root/ggh_m150_8TeV");
// /*  if(model!="GGH"){
//   sigTree1->Add(TString::Format("histograms_CMS-HGG_17042014_MC.root/vbf_m150_8TeV", iMass));
//   if(model!="VBF"){		 				   
//   sigTree1->Add(TString::Format("histograms_CMS-HGG_17042014_MC.root/wzh_m150_8TeV", iMass));
//   sigTree1->Add(TString::Format("histograms_CMS-HGG_17042014_MC.root/tth_m150_8TeV", iMass));
//   }
//   }*/

//   sigTree1->SetTitle("sigTree1");
//   sigTree1->SetName("sigTree1");


//   // common preselection cut
//   TString mainCut = "PhotonsMass>=130 && PhotonsMass<=1000";   // livia
  
  
//   // Create signal dataset composed with different productions, the weight is already applied in our ntuples
//   RooDataSet sigWeighted("sigWeighted","dataset",sigTree1,*ntplVars,mainCut,"evweight");
//   cout << endl;
//   cout << "sigWeighted" << endl;
//   sigWeighted.Print("v");
//   cout << "---- nX:  " << sigWeighted.sumEntries() << endl; 
  
//   // apply a common preselection cut; split in categories
//   cout << endl;
//   RooDataSet* signal[NCAT];
//   for (int c=0; c<ncat; ++c) {

  
   
//     // 1)  prime 4 cat livia
//     if (c==0) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("PhotonsMass"),mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
//     if (c==1) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("PhotonsMass"),mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
//     if (c==2) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("PhotonsMass"),mainCut+TString::Format("&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
//     if (c==3) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("PhotonsMass"),mainCut+TString::Format("&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));

//     w->import(*signal[c],Rename(TString::Format("SigWeight_cat%d",c)));
    
//     cout << "cat " << c << ", signal[c]: " << endl;
//     signal[c]->Print("v");
//     cout << "---- for category " << c << ", nX for signal[c]:  " << signal[c]->sumEntries() << endl; 
//     cout << endl;
//   }

//   // Create full weighted signal data set without categorization
//   RooDataSet* signalAll = (RooDataSet*) sigWeighted.reduce(*w->var("PhotonsMass"),mainCut);
//   w->import(*signalAll, Rename("SigWeight"));
//   cout << "now signalAll" << endl;
//   signalAll->Print("v");
//   cout << "---- nX for signalAll:  " << signalAll->sumEntries() << endl; 
//   cout << endl;
// }

RooDataHist* binnedDataset(RooWorkspace* w){
  // Create binning object with range (-10,10)
  const int nMassBins = 103;

  double massBoundaries[nMassBins+1] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000};
  RooBinning abins(nMassBins,massBoundaries,"mjj_binning") ;
  RooDataSet* data = (RooDataSet*) w->data(TString::Format("Data_cat0"));
  RooRealVar* mjj = (RooRealVar*) data->get()->find("mjj");

  mjj->setBinning(abins) ;
  RooDataHist* dataBinned = new RooDataHist(TString::Format("data_obs_cat0"),"",RooArgSet(*mjj),*data) ;
  return dataBinned;
}

// Data dataset
void AddBkgData(RooWorkspace* w, Float_t mass) {

  // initializations
  Int_t ncat = NCAT;

  Float_t minMassFit, maxMassFit;
  minMassFit = MINmass;
  maxMassFit = MAXmass;

  // retrieve the data tree; no common preselection cut applied yet; 
  string inDir = "/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/output/";
  TChain* dataChain[9];
  for (int i=0; i<9; ++i) {
    dataChain[i]  = new TChain("rootTupleTree/tree");
  }
  dataChain[0]->Add((inDir+"rootfile_QCD_Pt-300to470_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_trkalmb_castor_PHYS14_25_V1-v2__MINIAODSIM_reduced_skim.root").c_str());
  dataChain[1]->Add((inDir+"rootfile_QCD_Pt-470to600_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_trkalmb_castor_PHYS14_25_V1-v2__MINIAODSIM_reduced_skim.root").c_str());
  dataChain[2]->Add((inDir+"rootfile_QCD_Pt-600to800_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_trkalmb_castor_PHYS14_25_V1-v1__MINIAODSIM_reduced_skim.root").c_str());
  dataChain[3]->Add((inDir+"rootfile_QCD_Pt-800to1000_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_trkalmb_castor_PHYS14_25_V1-v2__MINIAODSIM_reduced_skim.root").c_str());
  dataChain[4]->Add((inDir+"rootfile_QCD_Pt-1000to1400_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_trkalmb_castor_PHYS14_25_V1-v1__MINIAODSIM_reduced_skim.root").c_str());
  dataChain[5]->Add((inDir+"rootfile_QCD_Pt-1400to1800_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_trkalmb_castor_PHYS14_25_V1-v1__MINIAODSIM_reduced_skim.root").c_str());
  dataChain[6]->Add((inDir+"rootfile_QCD_Pt-1800to2400_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_trkalmb_PHYS14_25_V1-v2__MINIAODSIM_reduced_skim.root").c_str());
  dataChain[7]->Add((inDir+"rootfile_QCD_Pt-2400to3200_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_trkalmb_PHYS14_25_V1-v1__MINIAODSIM_reduced_skim.root").c_str());
  dataChain[8]->Add((inDir+"rootfile_QCD_Pt-3200_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_trkalmb_PHYS14_25_V1-v1__MINIAODSIM_reduced_skim.root").c_str());

  // Variables
  RooArgSet* ntplVars = defineVariables();
  
  //giulia - non serve aggiungerla, gia` c'e`
  //ntplVars->add(*w->var("Dijet_MassW"));
  //(*w->var("Dijet_MassW")).setRange(MINmass, MAXmass);
  
  // common preselection cut
  TString mainCut = TString::Format(" mjj > 1100 && TMath::Abs(eta_j1)<2.5 && TMath::Abs(eta_j2)<2.5 && pt_j1>30 && pt_j2>30 ");  

  // Create dataset
  RooDataSet* Data[9];
  RooFormulaVar* wFunc[9]; 
  RooRealVar* lumi = (RooRealVar*)w->var("lumi");
  RooRealVar* xsec[9];
  xsec[0] = new RooRealVar(TString::Format("xsec"),"",7475.); 
  xsec[1] = new RooRealVar(TString::Format("xsec"),"",587.); 
  xsec[2] = new RooRealVar(TString::Format("xsec"),"",167.); 
  xsec[3] = new RooRealVar(TString::Format("xsec"),"",28.25); 
  xsec[4] = new RooRealVar(TString::Format("xsec"),"",8.195); 
  xsec[5] = new RooRealVar(TString::Format("xsec"),"",0.7346); 
  xsec[6] = new RooRealVar(TString::Format("xsec"),"",0.102); 
  xsec[7] = new RooRealVar(TString::Format("xsec"),"",0.00644); 
  xsec[8] = new RooRealVar(TString::Format("xsec"),"",0.000163); 
    
  RooRealVar* Ngen[9];
  Ngen[0] = new RooRealVar(TString::Format("Ngen"),"",1986177);
  Ngen[1] = new RooRealVar(TString::Format("Ngen"),"",2001071);
  Ngen[2] = new RooRealVar(TString::Format("Ngen"),"",1997744);
  Ngen[3] = new RooRealVar(TString::Format("Ngen"),"",1000065);
  Ngen[4] = new RooRealVar(TString::Format("Ngen"),"",500500);
  Ngen[5] = new RooRealVar(TString::Format("Ngen"),"",199627);
  Ngen[6] = new RooRealVar(TString::Format("Ngen"),"",200079);
  Ngen[7] = new RooRealVar(TString::Format("Ngen"),"",200453);
  Ngen[8] = new RooRealVar(TString::Format("Ngen"),"",200200);

  RooRealVar* eff[9];
  eff[0] = new RooRealVar(TString::Format("eff"),"",0.529013);
  eff[1] = new RooRealVar(TString::Format("eff"),"",0.5813171);
  eff[2] = new RooRealVar(TString::Format("eff"),"",0.612691);
  eff[3] = new RooRealVar(TString::Format("eff"),"",0.650460);
  eff[4] = new RooRealVar(TString::Format("eff"),"",0.693160);
  eff[5] = new RooRealVar(TString::Format("eff"),"",0.760388);
  eff[6] = new RooRealVar(TString::Format("eff"),"",0.827733);
  eff[7] = new RooRealVar(TString::Format("eff"),"",0.905659);
  eff[8] = new RooRealVar(TString::Format("eff"),"",0.96965);
  
  RooRealVar* weight[9];
  for (int i=0; i<9; ++i){
    Data[i] = new RooDataSet(TString::Format("Data"), TString::Format("Data_%d",i), *ntplVars, /*WeightVar("evweight"),*/ RooFit::Import(*dataChain[i]));
    //RooDataSet Data("Data","dataset",dataChain,*ntplVars,mainCut,"evweight");
    cout << endl;
    cout << "Data, everything: " << endl;
    Data[i]->Print("v");
    cout << "---- nX:  " << Data[i]->sumEntries() << endl;
    cout << endl;
  
    // Construct formula to calculate weight for events
    wFunc[i] = new RooFormulaVar(TString::Format("w"),"luminosity weight","@0*@1*@2*1.3/@3",RooArgSet(*lumi,*xsec[i],*eff[i],*Ngen[i])) ;

    //
    // Add column with variable w to previously generated dataset
    weight[i] = (RooRealVar*) Data[i]->addColumn(*wFunc[i]) ;
    weight[i]->SetName("lumi_weight");
  }
  cout << "debug: added weight to datasets" << endl;
  RooDataSet* Data_sum = (RooDataSet*) Data[0]->Clone("Data_sum");
  cout << "now merging the datasets" << endl;
  for (int i=1; i<9; i++){
    Data_sum ->append(*Data[i]); 
  }
  cout << "done with merging datasets" << endl;
  cout << "create weighted dataset" << endl;
  RooDataSet Data_weighted("Data_weighted","",Data_sum,*Data_sum->get(),0,"lumi_weight") ;
  cout << "is weighted ? " <<  Data_weighted.isWeighted() << endl;
  cout << "is non Poisson weighted? " << Data_weighted.isNonPoissonWeighted() << endl;
  cout << "done with weighted " << endl;

  // split into NCAT categories
  RooDataSet* dataToFit[NCAT];  
  for (int c=0; c<ncat; ++c) {
   
    cout << "Reducing dataset for cat " << c << endl;  

    // 1) just 1 category - giulia
    //if (c==0) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("Dijet_MassW"),mainCut+TString::Format("&& TMath::Abs(deltaETAjj)<1.3"));
    if (c==0) dataToFit[c] = (RooDataSet*) Data_weighted.reduce(mainCut+TString::Format("&& TMath::Abs(deltaETAjj)<1.3"));
  
    cout << endl; cout << "for category = " << c << endl;
    dataToFit[c]->Print("v");
    cout << "---- nX:  " << dataToFit[c]->sumEntries() << endl;
  
    cout << "Importing in the RooDataSet";
    w->import(*dataToFit[c],Rename(TString::Format("Data_cat%d",c)));
  }
  
  // cout << "data, no split" << endl;
  // // Create full data set without categorization
  // RooDataSet* data = (RooDataSet*) Data.reduce(*w->var("PhotonsMass"),mainCut);
  // w->import(*data, Rename("Data"));
  // cout << endl;
 
  // data->Print("v");
  // cout << "---- nX:  " << data->sumEntries() << endl; 
}



void makeRooHistPlot(RooHistFunc* rooFitMin,RooHistFunc* rooFitMax ,TF1* fmin, TF1* fmax, RooRealVar* var){

  TCanvas* cfit = new TCanvas("cfit", "cfit", 1);
  cfit->cd();
  
  RooPlot* p = var->frame();
  rooFitMin->plotOn(p, LineColor(kAzure+9));
  rooFitMax->plotOn(p, LineColor(kViolet+9));
  // p->Draw();
  fmax->SetLineColor(kViolet+9);
  fmax->SetLineWidth(2);
  fmin->SetLineWidth(2);
  //  fmax->SetLineStyle(kDashed);
  fmin->SetLineColor(kAzure+9);
  //fmin->SetLineStyle(kDashed);
  fmax->GetYaxis()->SetRangeUser(0, 1100);
  fmax->Draw("L");
  fmin->Draw("Lsame");
  fmax->GetYaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
  fmax->GetXaxis()->SetTitle("m_{X} [GeV]");

  TLegend* leg = new TLegend(0.52, 0.16, 0.85, 0.35, "", "brNDC");
  //  leg->SetTextSize(0.0206044);  
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(fmax, "Fit Range Maximum", "L");
  leg->AddEntry(fmin, "Fit Range Minimum", "L");
  leg->Draw("same");
  

  //TPaveText* label_cms = get_labelCMS(0, "2014", false);
  //TPaveText* label_sqrt = get_labelSqrt(0);
  //  label_cms->Draw("same");
  //label_sqrt->Draw("same");

  int iPos=11 ;
  //CMS_lumi( cfit,true,iPos );
  cfit->SaveAs("plots/bkg_fitRange.png");
  cfit->SaveAs("plots/bkg_fitRange.pdf");
  cfit->SaveAs("plots/BkgFit/bkg_fitRange.pdf");
  cfit->SaveAs("plots/BkgFit/bkg_fitRange.png");

  // cfit->SaveAs("~/www/plotsPAS/bkg_fitRange.pdf");
  // cfit->SaveAs("~/www/plotsPAS/bkg_fitRange.png");
  // cfit->SaveAs("~/www/plotsPAS/bkg_fitRange.C");
  // cfit->SaveAs("~/www/plotsPAS/bkg_fitRange.root");


}

/////////////////////////////////////////////////
// not useful for now
/////////////////////////////////////////////////


// void SetConstantParams(const RooArgSet* params) {

//   cout << endl; cout << "Entering SetConstantParams" << endl;
//   TIterator* iter(params->createIterator());
//   for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
//     RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
//     if (rrv) { rrv->setConstant(true); std::cout << " " << rrv->GetName(); }
//   }  
// }



// Double_t effSigma(TH1 *hist) {
  
//   TAxis *xaxis = hist->GetXaxis();
//   Int_t nb = xaxis->GetNbins();
//   if(nb < 10) {
//     std::cout << "effsigma: Not a valid histo. nbins = " << nb << std::endl;
//     return 0.;
//   }
  
//   Double_t bwid = xaxis->GetBinWidth(1);
//   if(bwid == 0) {
//     std::cout << "effsigma: Not a valid histo. bwid = " << bwid << std::endl;
//     return 0.;
//   }
//   Double_t xmax = xaxis->GetXmax();
//   Double_t xmin = xaxis->GetXmin();
//   Double_t ave = hist->GetMean();
//   Double_t rms = hist->GetRMS();

//   Double_t total=0.;
//   for(Int_t i=0; i<nb+2; i++) {
//     total+=hist->GetBinContent(i);
//   }
//   if(total < 100.) {
//     std::cout << "effsigma: Too few entries " << total << std::endl;
//     return 0.;
//   }
//   Int_t ierr=0;
//   Int_t ismin=999;

//   Double_t rlim=0.683*total;
//   Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
//   if(nrms > nb/10) nrms=nb/10; // Could be tuned...

//   Double_t widmin=9999999.;
//   for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
//     Int_t ibm=(ave-xmin)/bwid+1+iscan;
//     Double_t x=(ibm-0.5)*bwid+xmin;
//     Double_t xj=x;
//     Double_t xk=x;
//     Int_t jbm=ibm;
//     Int_t kbm=ibm;
//     Double_t bin=hist->GetBinContent(ibm);
//     total=bin;
//     for(Int_t j=1;j<nb;j++){
//       if(jbm < nb) {
//         jbm++;
//         xj+=bwid;
//         bin=hist->GetBinContent(jbm);
//         total+=bin;
//         if(total > rlim) break;
//       }
//       else ierr=1;
//       if(kbm > 0) {
//         kbm--;
//         xk-=bwid;
//         bin=hist->GetBinContent(kbm);
//         total+=bin;
//         if(total > rlim) break;
//       }
//       else ierr=1;
//     }
//     Double_t dxf=(total-rlim)*bwid/bin;
//     Double_t wid=(xj-xk+bwid-dxf)*0.5;
//     if(wid < widmin) {
//       widmin=wid;
//       ismin=iscan;
//     }
//   }
//   if(ismin == nrms || ismin == -nrms) ierr=3;
//   if(ierr != 0) std::cout << "effsigma: Error of type " << ierr << std::endl;

//   return widmin;
// }

// void R2JJFitter(double mass, std::string postfix="")
// {
//     filePOSTfix=postfix;
//     if(postfix!="")
//     {
//       // for optimization studies
//       MMIN=1000;
//       if(mass==1000)
//          signalScaler=0.034246;
//       if((mass>1000)&&(mass<2000))
//          signalScaler=0.02469;
//       if(mass==2000)
//          signalScaler=2.0;
//     };
//     runfits(mass, 1);
//     if(postfix!="")
//     {
//       // for optimization studies
//       MMIN=1000;
//       if(mass==1000)
//          signalScaler=0.033500;
//       if((mass>1000)&&(mass<2000))
//          signalScaler=0.02016;
//       if(mass==2000)
//          signalScaler=2.22222;
//     };
//     runfits(mass, 0);
//     runfits(mass, 2);
// }



// RooHistFunc* getNorm2D(int cat,RooRealVar* var,double mass, double width, std::string model ){
  
//   double norm;
//   TFile* f;
//   if(model=="GRAVITON") f= new TFile("effXacc2D_GRAVITON.root", "READ");
//   else f= new TFile("effXacc2D.root", "READ");
//   f->Print();
//   TH2D* h2 = (TH2D*) f->Get(TString::Format("h2_cat%d",cat));
//   int bin = h2->GetYaxis()->FindBin(width);
//   std::cout<<bin<<std::endl;
 
//   TH1D* h1 = (TH1D*)h2->ProjectionX("py", bin, bin);

//   double eff;
//   eff = h1->GetBinContent(h1->FindBin(mass));
//   norm = 1.87004999e-02*eff/100.;
//   for(int i = 0; i<h1->GetNbinsX()+1; i++){
//     std::cout<< i<<" M: "<<h1->GetBinCenter(i)<<"   "<<h1->GetBinContent(i)<<std::endl;
//     h1->SetBinContent(i,1.87004999e-02*h1->GetBinContent(i)/100);
//     //std::cout<<h1->GetBinCenter(i)<<"   "<<h1->GetBinContent(i)<<std::endl;
//   }
//   std::cout<<norm<<std::endl;
//   RooDataHist* rooData_all = new RooDataHist("rooData_all", "roData_all",*var,h1);
//   RooHistFunc* rooFunc_all = new RooHistFunc(TString::Format("PhotonsMassSig_cat%d_norm",cat),TString::Format("PhotonsMassSig_cat%d_norm",cat), *var,*rooData_all, 3);
//   // std::cout<<rooFunc_all->getVal(mass)<<std::endl;
//   TCanvas* c = new TCanvas("c", "c", 1);
//   c->cd();
//   //  TPaveText* label_cms = get_labelCMS(0, false);
//   //TPaveText* label_sqrt = get_labelSqrt(0);  
//   // h1->GetYaxis()->SetRangeUser(0.0001, 0.1);
//   h1->GetXaxis()->SetTitle("m_{X} [GeV]");
//   h1->GetYaxis()->SetTitle("Signal Yield/1pb");
//   h1->Draw();
//   RooPlot* plot = var->frame();
//   rooData_all->plotOn(plot);
//   rooFunc_all->plotOn(plot);
//   plot->Draw("same");
//   c->SetLogy();

//   c->SaveAs(TString::Format("~/www/signalYield_cat%d.png",cat));
//   c->SaveAs(TString::Format("~/www/signalYield_cat%d.pdf",cat));

//   return rooFunc_all;

// }


//////////////////////////

