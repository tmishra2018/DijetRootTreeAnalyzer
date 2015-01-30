// Signal Data Set
void MakeRooKeysPDFMCBkg(RooWorkspace* w, Float_t mass, Bool_t isMirror) {


  TString wsDir = "BiasStudy/workspaces/"+filePOSTfix;
  
  RooWorkspace *wBias = new RooWorkspace("w_bias","w_bias");  


  Int_t ncat = NCAT;
  TString inDir = "";

  Float_t MASS(mass);

  // Luminosity:
  Float_t Lum = 19500.0;  
  RooRealVar lumi("lumi","lumi",Lum);
  wBias->import(lumi); 
  
  // Variables
  RooArgSet* ntplVars = defineVariables();
  ntplVars->Print("V");
  RooArgSet* ntplVars_newweight = defineVariables_newWeight();
  int iMass = abs(mass);  
  RooRealVar* PhotonsMass = w->var("PhotonsMass");  
  PhotonsMass->setRange(MINmass, MAXmass);  
  TFile sigFile1("histograms_CMS-HGG_24072013.root");   //ggh prod mode tree livia
  
  // common preselection cut
  //TString mainCut = "PhotonsMass>=100 && PhotonsMass<=1000";   // livia
  TString mainCut = "mjjAK8>=100 && mjjAK8<=1000";   //giulia
  
  //get sumEntries of QCD 
  TChain* qcdTree=new TChain();
  //qcdTree->Add("histograms_CMS-HGG_24072013.root/qcd_30_8TeV_pf");
  //qcdTree->Add("histograms_CMS-HGG_24072013.root/qcd_40_8TeV_pf");


  //giulia : VA AGGIUNTO IL PESO!

  qcdTree->Add("/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_1_0_pre9_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/data/output/rootFile_QCD_Pt-300to470_reduced_skim_withWeight.root/tree");
  qcdTree->Add("/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_1_0_pre9_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/data/output/rootFile_QCD_Pt-470to600_reduced_skim_withWeight.root/tree");
  qcdTree->Add("/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_1_0_pre9_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/data/output/rootFile_QCD_Pt-600to800_reduced_skim_withWeight.root/tree");
  qcdTree->Add("/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_1_0_pre9_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/data/output/rootFile_QCD_Pt-800to1000_reduced_skim_withWeight.root/tree");
  qcdTree->Add("/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_1_0_pre9_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/data/output/rootFile_QCD_Pt-1000to1400_reduced_skim_withWeight.root/tree");
  qcdTree->Add("/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_1_0_pre9_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/data/output/rootFile_QCD_Pt-1400to1800_reduced_skim_withWeight.root/tree");
  qcdTree->Add("/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_1_0_pre9_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/data/output/rootFile_QCD_Pt-1800to2400_reduced_skim_withWeight.root/tree");
  qcdTree->Add("/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_1_0_pre9_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/data/output/rootFile_QCD_Pt-2400to3200_reduced_skim_withWeight.root/tree");
  qcdTree->Add("/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_1_0_pre9_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/data/output/rootFile_QCD_Pt-3300_reduced_skim_withWeight.root/tree");

 
  RooDataSet qcdMCWeighted("qcdMCWeighted","MC qcd weighted",qcdTree,*ntplVars,mainCut,"evweight");
  cout << endl;
  cout << "qcdMCWeighted" << endl;
  qcdMCWeighted.Print("v");
  Double_t qcdInt_ = qcdMCWeighted.sumEntries();
  cout << "---- nX: qcd Int " << qcdInt_ << endl; 

  //giulia gjets non serve
  //get sumEntries of QCD 
  // TChain* gjTree=new TChain();
  // gjTree->Add("histograms_CMS-HGG_24072013.root/gjet_20_8TeV_pf");
  // gjTree->Add("histograms_CMS-HGG_24072013.root/gjet_40_8TeV_pf");

  //  RooDataSet gjMCWeighted("gjMCWeighted","MC gj weighted",gjTree,*ntplVars,mainCut,"evweight");
  // cout << endl;
  // cout << "gjMCWeighted" << endl;
  // gjMCWeighted.Print("v");
  // Double_t gjInt_ = gjMCWeighted.sumEntries();
  // cout << "---- nX: gj Int " << gjInt_ << endl; 

 
  //chain summing up all production modes
  TChain* sigTree1  = new TChain();
  // sigTree1->Add("histograms_CMS-HGG_24072013.root/gjet_20_8TeV_pf");
  // sigTree1->Add("histograms_CMS-HGG_24072013.root/gjet_40_8TeV_pf");
  // sigTree1->Add("histograms_CMS-HGG_24072013.root/diphojet_8TeV");
  // sigTree1->Add("histograms_CMS-HGG_24072013.root/dipho_Box_25_8TeV");
  // sigTree1->Add("histograms_CMS-HGG_24072013.root/dipho_Box_250_8TeV");
  

  sigTree1->Add("/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_1_0_pre9_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/data/output/rootFile_RSGravToJJ_kMpl01_M-5000_reduced_skim_withWeight.root/tree");

  sigTree1->SetTitle("sigTree1");
  sigTree1->SetName("sigTree1");
  
  
  // Create signal dataset composed with different productions, the weight is already applied in our ntuples
  RooDataSet BkgMCWeighted("BkgMCWeighted","MC BKG weighted",sigTree1,*ntplVars,mainCut,"evweight");
  cout << endl;
  cout << "BkgMCWeighted" << endl;
  BkgMCWeighted.Print("v");
  cout << "---- nX:  " << BkgMCWeighted.sumEntries() << endl; 
  
  // apply a common preselection cut; split in categories
  cout << endl;
 
  RooDataSet* data[NCAT];
  RooDataSet* qcdMC[NCAT];
  //RooDataSet* gjMC[NCAT];
  TTree* gjAndQcdTree[NCAT];
  RooDataSet* QcdToGjMC[NCAT];

  RooDataSet* BkgMC[NCAT];
  TTree* BkgMCcopyTree[NCAT];
  RooDataSet* BkgMCcopy[NCAT];
  RooKeysPdf* BkgMCKeyPdf[NCAT];
  RooKeysPdf* BkgMCKeyPdf_bw2[NCAT];
  RooKeysPdf* BkgMCKeyPdf_bw3[NCAT];
  RooKeysPdf* BkgMCKeyPdf_bw4[NCAT];
  RooKeysPdf* BkgMCKeyPdf_bw2_noMirr[NCAT];


  RooDerivative* BkgMCKeyPdf_D1[NCAT];
  RooDerivative* BkgMCKeyPdf_bw2_D1[NCAT];
  RooDerivative* BkgMCKeyPdf_bw3_D1[NCAT];
  RooDerivative* BkgMCKeyPdf_bw4_D1[NCAT];
  RooDerivative* BkgMCKeyPdf_bw2_D1_noMirr[NCAT];

  RooDerivative* BkgMCKeyPdf_D2[NCAT];
  RooDerivative* BkgMCKeyPdf_bw2_D2[NCAT];
  RooDerivative* BkgMCKeyPdf_bw3_D2[NCAT];
  RooDerivative* BkgMCKeyPdf_bw4_D2[NCAT];
  RooDerivative* BkgMCKeyPdf_bw2_D2_noMirr[NCAT];

  RooPlot* plotPhotonsMassBkgMC[NCAT];
  RooPlot* plotPhotonsMassBkgMC_D1[NCAT];
  RooPlot* plotPhotonsMassBkgMC_D2[NCAT];
  TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
  Int_t nBinsMass(200);
  Double_t  minMassFit = MINmass;
  Double_t  maxMassFit = MAXmass;

  //  RooArgSet* argset_ = new RooArgSet(*w->var("PhotonsMass"), *w->var("evweight"));

  for (int c=1; c<2; ++c) {

    // 0) chiara: 1cat only
    // signal[c] =  (RooDataSet*) sigWeighted.reduce(*w->var("massggnewvtx"),mainCut);   

    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));

    // 1)  prime 4 cat livia
  
   
    //reduce QCD dataset
    if (c==0) qcdMC[c] = (RooDataSet*) qcdMCWeighted.reduce(*ntplVars, mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==1) qcdMC[c] = (RooDataSet*) qcdMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
    if (c==2) qcdMC[c] = (RooDataSet*) qcdMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==3) qcdMC[c] = (RooDataSet*) qcdMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9>0.94 || ph2_r9>0.94 ) "));

 
    //reduce gj dataset
    if (c==0) gjMC[c] = (RooDataSet*) gjMCWeighted.reduce(*ntplVars, mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==1) gjMC[c] = (RooDataSet*) gjMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
    if (c==2) gjMC[c] = (RooDataSet*) gjMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==3) gjMC[c] = (RooDataSet*) gjMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9>0.94 || ph2_r9>0.94 ) "));
 
    gjMC[c]->Print();



    Float_t qcdInt = qcdMC[c]->sumEntries();
    Float_t gjInt = gjMC[c]->sumEntries();
    Float_t gjEntries = gjMC[c]->numEntries();
    std::cout<<"qcd: "<<qcdInt<<" gj: "<<gjInt<<" gjEntries: "<<gjEntries<<std::endl;
    Float_t qcdWeight = qcdInt/gjInt;
 
    
    gjAndQcdTree[c] = (TTree*)dataset2tree(gjMC[c], ntplVars, qcdWeight);
    
   
    QcdToGjMC[c] = new RooDataSet(gjMC[c]->GetName(),gjMC[c]->GetTitle(),gjAndQcdTree[c],*ntplVars_newweight, mainCut, "newweight"  );  
    QcdToGjMC[c]->Print("");
  
    


    if (c==0) BkgMC[c] = (RooDataSet*) BkgMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==1) BkgMC[c] = (RooDataSet*) BkgMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
    if (c==2) BkgMC[c] = (RooDataSet*) BkgMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==3) BkgMC[c] = (RooDataSet*) BkgMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9>0.94 || ph2_r9>0.94 ) "));


    BkgMC[c]->Print();   
    BkgMCcopyTree[c] = (TTree*) dataset2tree(BkgMC[c], ntplVars, 1.);

    BkgMCcopy[c] =  new RooDataSet(BkgMC[c]->GetName(),BkgMC[c]->GetTitle(),BkgMCcopyTree[c] ,*ntplVars_newweight, mainCut, "newweight");     
    BkgMCcopy[c]->append(*QcdToGjMC[c]);
    BkgMCcopy[c]->Print();


   

    wBias->import(*BkgMCcopy[c],Rename(TString::Format("BkgMCWeight_cat%d",c)));
    
    cout << "cat " << c << ", BkgMC[c]: " << endl;
    BkgMCcopy[c]->Print("v");
    cout << "---- for category " << c << ", nX for [c]:  " << BkgMCcopy[c]->sumEntries() << endl; 
    cout << endl;


    isMirror = true;
    if(isMirror){
    
    //create the rookeyspdf
      //   BkgMCKeyPdf[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_cat%d",c),TString::Format("BkgMCKeyPdf_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::MirrorBoth );
      BkgMCKeyPdf_bw2[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw2_cat%d",c),TString::Format("BkgMCKeyPdf_bw2_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::MirrorBoth,2 );
      BkgMCKeyPdf_bw3[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw3_cat%d",c),TString::Format("BkgMCKeyPdf_bw3_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::MirrorBoth,3 );
      BkgMCKeyPdf_bw4[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw4_cat%d",c),TString::Format("BkgMCKeyPdf_bw4_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::MirrorBoth,4 );
      
    }else{

      //create the rookeyspdf
      BkgMCKeyPdf[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_cat%d",c),TString::Format("BkgMCKeyPdf_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::NoMirror );
      BkgMCKeyPdf_bw2[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw2_cat%d",c),TString::Format("BkgMCKeyPdf_bw2_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::NoMirror,2 );
      BkgMCKeyPdf_bw3[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw3_cat%d",c),TString::Format("BkgMCKeyPdf_bw3_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::NoMirror,3 );
      BkgMCKeyPdf_bw4[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw4_cat%d",c),TString::Format("BkgMCKeyPdf_bw4_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::NoMirror,4 );
    

    }


    wBias->import(*BkgMCKeyPdf_bw2[c]);
    //   wBias->import(*BkgMCKeyPdf_bw[c]);
    wBias->import(*BkgMCKeyPdf_bw3[c]);
    wBias->import(*BkgMCKeyPdf_bw4[c]);
    
    plotPhotonsMassBkgMC[c] = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
    data[c]->plotOn(plotPhotonsMassBkgMC[c]);
    //BkgMC[c]->plotOn(plotPhotonsMassBkgMC[c],"PE", MarkerColor(kRed), LineColor(kRed), MarkerSize(1.));
    //BkgMCcopy[c]->plotOn(plotPhotonsMassBkgMC[c],"PE", MarkerColor(kGreen), LineColor(kGreen), MarkerSize(0.5));
    //   BkgMCKeyPdf[c]->plotOn(plotPhotonsMassBkgMC[c],"L", LineColor(kBlue), LineWidth(2));
    BkgMCKeyPdf_bw2[c]->plotOn(plotPhotonsMassBkgMC[c],"L", LineColor(kBlack), LineWidth(2));
    BkgMCKeyPdf_bw3[c]->plotOn(plotPhotonsMassBkgMC[c],"L", LineColor(8), LineWidth(2));
    BkgMCKeyPdf_bw4[c]->plotOn(plotPhotonsMassBkgMC[c],"L", LineColor(kOrange), LineWidth(2));

    plotPhotonsMassBkgMC[c]->SetAxisRange(0.001,plotPhotonsMassBkgMC[c]->GetMaximum()*30.,"Y");
    
    /*
    //create the rookeyspdf_D1

    BkgMCKeyPdf_D1[c] =   (RooDerivative*) BkgMCKeyPdf[c]->derivative(*PhotonsMass, 1, 0.01);
    BkgMCKeyPdf_bw2_D1[c] =  (RooDerivative*)   BkgMCKeyPdf_bw2[c]->derivative(*PhotonsMass, 1, 0.01);
    BkgMCKeyPdf_bw3_D1[c] =   (RooDerivative*) BkgMCKeyPdf_bw3[c]->derivative(*PhotonsMass, 1, 0.01);
    BkgMCKeyPdf_bw4_D1[c] =   (RooDerivative*) BkgMCKeyPdf_bw4[c]->derivative(*PhotonsMass, 1, 0.01);


    //create the rookeyspdf_D2
    BkgMCKeyPdf_D2[c] =   (RooDerivative*) BkgMCKeyPdf[c]->derivative(*PhotonsMass, 2, 0.01);
    BkgMCKeyPdf_bw2_D2[c] = (RooDerivative*)  BkgMCKeyPdf_bw2[c]->derivative(*PhotonsMass, 2,0.01);
    BkgMCKeyPdf_bw3_D2[c] = (RooDerivative*)   BkgMCKeyPdf_bw3[c]->derivative(*PhotonsMass, 2,0.01);
    BkgMCKeyPdf_bw4_D2[c] =  (RooDerivative*)  BkgMCKeyPdf_bw4[c]->derivative(*PhotonsMass, 2,0.01);

    
    plotPhotonsMassBkgMC_D1[c] = PhotonsMass->frame(minMassFit, maxMassFit);
    BkgMCKeyPdf_D1[c]->plotOn(plotPhotonsMassBkgMC_D1[c],"L", LineColor(kBlue), LineWidth(2));
    BkgMCKeyPdf_bw2_D1[c]->plotOn(plotPhotonsMassBkgMC_D1[c],"L", LineColor(kBlack), LineWidth(2));
    BkgMCKeyPdf_bw3_D1[c]->plotOn(plotPhotonsMassBkgMC_D1[c],"L", LineColor(8), LineWidth(2));
    BkgMCKeyPdf_bw4_D1[c]->plotOn(plotPhotonsMassBkgMC_D1[c],"L", LineColor(kOrange), LineWidth(2));

    plotPhotonsMassBkgMC_D1[c]->SetAxisRange(plotPhotonsMassBkgMC_D1[c]->GetMinimum()*1.3,0.15,"Y");   
    
   
    plotPhotonsMassBkgMC_D2[c] = PhotonsMass->frame(minMassFit, maxMassFit);
    BkgMCKeyPdf_D2[c]->plotOn(plotPhotonsMassBkgMC_D2[c],"", LineColor(kBlue), LineWidth(2));
    BkgMCKeyPdf_bw2_D2[c]->plotOn(plotPhotonsMassBkgMC_D2[c],"", LineColor(kBlack), LineWidth(2));
    BkgMCKeyPdf_bw3_D2[c]->plotOn(plotPhotonsMassBkgMC_D2[c],"", LineColor(8), LineWidth(2));
    BkgMCKeyPdf_bw4_D2[c]->plotOn(plotPhotonsMassBkgMC_D2[c],"", LineColor(kOrange), LineWidth(2));

    plotPhotonsMassBkgMC_D2[c]->SetAxisRange(plotPhotonsMassBkgMC_D2[c]->GetMinimum()*1.3,0.04,"Y");
    
    */
 
    TLegend *leg1 = new TLegend(0.2175,0.8236441,0.6575,0.9240678, TString::Format("Category %d",c), "brNDC");
    leg1->AddEntry(plotPhotonsMassBkgMC[c]->getObject(0),"Bkg MC","LPE");

    TLegend *leg2 = new TLegend(0.4375,0.7236441,0.85,0.9240678, TString::Format("RooKeysPdf",c), "brNDC");

    TLegend *leg3 = new TLegend(0.2175,0.8236441,0.6575,0.9240678, TString::Format("Category %d",c), "brNDC");
    if(isMirror){
      leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(1),"Default Bw","L");
      leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(2),"Bw x 2","L");
      leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(3),"Bw x 3","L");
      /* leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(4),"Bw x 4","L");
      leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(4),"Bw x 4","L");
      leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(4),"Bw x 4","L");*/
    }else{
      leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(1),"Default bw NoMirr","L");
      leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(2),"Bw x 2 NoMirr","L");
      leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(3),"Bw x 3 NoMirr","L");
      leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(4),"Bw x 4 NoMirr","L");
    }
    
    leg1->SetTextSize(0.035);
    leg1->SetTextFont(42);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
   
    leg2->SetTextSize(0.035);
    leg2->SetTextFont(42);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);

    leg3->SetTextSize(0.035);
    leg3->SetTextFont(42);
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
   
   
    TPaveText* label_cms = get_labelCMS(0, "2012", true);
    TPaveText* label_sqrt = get_labelSqrt(0);

  
    ctmp->cd();
    plotPhotonsMassBkgMC[c]->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    plotPhotonsMassBkgMC[c]->Draw();  
    leg1->Draw("same");
    leg2->Draw("same");
    
    label_cms->Draw("same");
    label_sqrt->Draw("same");

   
    ctmp->SetLogy(1);
    if(isMirror){
    ctmp->SaveAs("plots/BKG_MC_rookeyspdf_"+TString::Format("cat%d_LOG.png", c));
    ctmp->SaveAs("plots/BKG_MC_rookeyspdf_"+TString::Format("cat%d_LOG.pdf", c));
    ctmp->SaveAs("plots/BKG_MC_rookeyspdf_"+TString::Format("cat%d_LOG.root", c));
    }else{
      ctmp->SaveAs("plots/BKG_MC_rookeyspdf_NoMirr_"+TString::Format("cat%d_LOG.png", c));
      ctmp->SaveAs("plots/BKG_MC_rookeyspdf_NoMirr_"+TString::Format("cat%d_LOG.pdf", c));
      ctmp->SaveAs("plots/BKG_MC_rookeyspdf_NoMirr_"+TString::Format("cat%d_LOG.root", c));
  
    }
       ctmp->SetLogy(0);
       /*  //plot D1
       plotPhotonsMassBkgMC_D1[c]->GetYaxis()->SetTitle("First Derivative ");
       plotPhotonsMassBkgMC_D1[c]->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
       plotPhotonsMassBkgMC_D1[c]->Draw();  
       leg2->SetHeader("RooKeysPDF - 1st Derivative");
       leg2->Draw("same");
       leg3->Draw("same");
       
       label_cms->Draw("same");
       label_sqrt->Draw("same");
       if(isMirror){
	 ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D1_"+TString::Format("cat%d.png", c));
	 ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D1_"+TString::Format("cat%d.pdf", c));
	 ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D1_"+TString::Format("cat%d.root", c));
       }else{
	 ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D1_NoMirr_"+TString::Format("cat%d_LOG.png", c));
	 ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D1_NoMirr_"+TString::Format("cat%d_LOG.pdf", c));
	 ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D1_NoMirr_"+TString::Format("cat%d_LOG.root", c));
	 
       }  
       ctmp->SetLogy(0);
       ctmp->Clear();

    //plot D2
    plotPhotonsMassBkgMC_D2[c]->GetYaxis()->SetTitle("Second Derivative");
    plotPhotonsMassBkgMC_D2[c]->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    plotPhotonsMassBkgMC_D2[c]->Draw();  
    leg2->SetHeader("RooKeysPDF - 2nd Derivative");
    leg2->Draw("same");
    leg3->Draw("same");
	
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    if(isMirror){
    ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D2_"+TString::Format("cat%d.png", c));
    ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D2_"+TString::Format("cat%d.pdf", c));
    ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D2_"+TString::Format("cat%d.root", c));
    }else{
      ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D2_NoMirr_"+TString::Format("cat%d_LOG.png", c));
      ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D2_NoMirr_"+TString::Format("cat%d_LOG.pdf", c));
      ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D2_NoMirr_"+TString::Format("cat%d_LOG.root", c));
      
      } */
  }

  
    
  std::cout << "done with importing MC background datasets and RooKeysPdfs" << std::endl;
  

  TString filename(wsDir+"HighMass-hgg.RooKeysPdfMCBkg_8TeV_cat123.root");
  wBias->writeToFile(filename);
  cout << "Write background RooKeys workspace in: " << filename << " file" << endl;

  return;

}

