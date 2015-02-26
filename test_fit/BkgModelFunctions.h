RooFitResult* BkgModelFitDiJetFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, Int_t c,  bool blind) {
//void BkgModelFitDiJetFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, Int_t c,  bool blind) {

  //Int_t ncat = NCAT;
  std::cout<<"isBlind: "<<blind<<std::endl;
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data;
  RooDataHist* dataBinned;
  RooFitResult* fitresult;
  
  RooPlot* plotDijetMassBkg;

  // dobands and dosignal
  RooDataSet* signal;
  RooAbsPdf* DijetMassSig;
  
  Float_t minMassFit, maxMassFit;
  minMassFit = MINmass;
  maxMassFit = MAXmass;
   
  // Fit data with background pdf for data limit
  RooRealVar* mjj = w->var("mjj");  
  mjj->setUnit("GeV");
 
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  

  data = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
  dataBinned = (RooDataHist*) w->data(TString::Format("data_obs_cat%d",c));

 
  //giulia test: binning 1 GeV
  //RooBinning oneGeVbins(13999,1.,14000.) ;
  //mjj->setBinning(oneGeVbins) ;
  //RooDataHist* dataBinned_1GeV = new RooDataHist(TString::Format("data_1GeV_cat0"),"",RooArgSet(*mjj),*data) ; 
  
  
  cout << "start debug" << endl;
  //Debug
  w->Print();
  RooRealVar* debug1 = (RooRealVar*)(w->var(TString::Format("DijetsMass_bkg_13TeV_norm_cat%d",c))); 
  RooRealVar* debug2 =  (RooRealVar*)(w->var(TString::Format("DijetsMass_bkg_13TeV_slope1_cat%d",c))); 
  RooRealVar* debug3 =  (RooRealVar*)(w->var(TString::Format("DijetsMass_bkg_13TeV_slope2_cat%d",c))); 
  RooRealVar* debug4 =  (RooRealVar*)(w->var(TString::Format("DijetsMass_bkg_13TeV_slope3_cat%d",c))); 
  debug1->Print();
  debug2->Print();
  debug3->Print();
  debug4->Print();
  cout << "end debug" << endl;

  // fit a la dijets
  RooFormulaVar *p0mod = new RooFormulaVar(TString::Format("par0DiJet_cat%d",c),"","@0",*w->var(TString::Format("DijetsMass_bkg_13TeV_norm_cat%d",c)));
  RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("par1DiJet_cat%d",c),"","@0",*w->var(TString::Format("DijetsMass_bkg_13TeV_slope1_cat%d",c)));
  RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("par2DiJet_cat%d",c),"","@0",*w->var(TString::Format("DijetsMass_bkg_13TeV_slope2_cat%d",c)));
  RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("par3DiJet_cat%d",c),"","@0",*w->var(TString::Format("DijetsMass_bkg_13TeV_slope3_cat%d",c)));
  mjj->setRange("bkg range", MINmass,MAXmass);
  RooFormulaVar *x     = new RooFormulaVar(TString::Format("xDiJet_cat%d",c),"","@0/13000.",*w->var("mjj"));

  
  RooAbsPdf* DijetMassBkgTmp0 = new RooGenericPdf(TString::Format("DijetMassBkg_cat%d",c), "pow(1-@0, @1)/pow(@0, @2+@3*log(@0))", RooArgList(*x, *p1mod, *p2mod,*p3mod));
  
  //fitresult = DijetMassBkgTmp0->fitTo(*data, RooFit::Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), RooFit::SumW2Error(kTRUE), RooFit::Save(kTRUE));
  //fitresult = DijetMassBkgTmp0->fitTo(*dataBinned, RooFit::Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), RooFit::SumW2Error(kTRUE), RooFit::Save(kTRUE));
  //fitresult = DijetMassBkgTmp0->fitTo(*dataBinned_1GeV, RooFit::Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), RooFit::SumW2Error(kTRUE), RooFit::Save(kTRUE));


  std::cout<<TString::Format("******************************** Background Fit results mass %f cat %d***********************************", mass, c)<<std::endl;
  //fitresult->Print("V");
   
  //Generate Pseudodataset and fit again with the dijet function
  //WARNING: at the moment number of events for 1 fb-1 by hand!!!! 
  //TO BE CHANGED	 

  //
  TH1F* h_data_binned_mjj = (TH1F*)dataBinned->createHistogram("mjj");
  int nbins = h_data_binned_mjj->GetNbinsX();
  double Xmin = h_data_binned_mjj->GetBinLowEdge(1);
  double Xmax = h_data_binned_mjj->GetBinLowEdge(nbins+1);
  cout << "#######################" << endl;
  cout << "h_data_binned_mjj : nbins="<<nbins << " Xmin="<< Xmin<<"  Xmax="<<Xmax << endl;
  TH1F* h_data_binned_w_mjj = new TH1F("h_data_binned_w_mjj","", nbins, Xmin, Xmax); 
  for(int i=1; i<=h_data_binned_mjj->GetNbinsX(); i++){
    h_data_binned_w_mjj->SetBinContent(i,h_data_binned_mjj->GetBinContent(i)/h_data_binned_w_mjj->GetBinWidth(i));
  }

  TCanvas c1;
  c1.cd();
  c1.SetLogy(1);
  h_data_binned_mjj->Draw()  ;
  c1.SaveAs("databinned_debug.png");



  TH1F* h_pseudodata_binned_mjj = new TH1F("h_pseudodata_binned_mjj","", nbins, Xmin, Xmax);
  for(int i=1; i<=h_data_binned_mjj->GetNbinsX(); i++){
    h_pseudodata_binned_mjj->SetBinContent(i,RooRandom::randomGenerator()->Poisson(h_data_binned_mjj->GetBinContent(i)));
    h_pseudodata_binned_mjj->SetBinError(i,TMath::Sqrt(h_data_binned_mjj->GetBinContent(i)));
    cout << "bin" << i << "content : " << h_pseudodata_binned_mjj->GetBinContent(i) <<  "  error : " << h_pseudodata_binned_mjj->GetBinError(i) << endl; 
  }
  TCanvas c2;
  c2.SetLogy(1);
  c2.cd();
  h_pseudodata_binned_mjj->Draw();
  c2.SaveAs("pseudodatabinned_debug.png");


  TH1F* h_pseudodata_binned_w_mjj = new TH1F("h_pseudodata_binned_w_mjj","", nbins, Xmin, Xmax);
  for(int i=1; i<=h_data_binned_mjj->GetNbinsX(); i++){
    h_pseudodata_binned_w_mjj->SetBinContent(i,h_pseudodata_binned_mjj->GetBinContent(i)/h_pseudodata_binned_w_mjj->GetBinWidth(i));
    h_pseudodata_binned_w_mjj->SetBinError(i,TMath::Sqrt(h_pseudodata_binned_mjj->GetBinContent(i)/h_pseudodata_binned_w_mjj->GetBinWidth(i)));
    cout << "random bin " << i << " content : " << h_pseudodata_binned_w_mjj->GetBinContent(i) <<"  error : " << h_pseudodata_binned_w_mjj->GetBinError(i) << endl;
  }
  TCanvas c3;
  c3.SetLogy(1);
  c3.cd();
  h_pseudodata_binned_w_mjj->Draw(); 
  c3.SaveAs("pseudodatabinned_w_debug.png");

  RooDataHist* pseudoDataBinned = new RooDataHist(TString::Format("pseudodata_cat0"),"",RooArgSet(*mjj),h_pseudodata_binned_mjj) ;
  RooDataHist* pseudoDataBinned_w = new RooDataHist(TString::Format("pseudodata_w_cat0"),"",RooArgSet(*mjj),h_pseudodata_binned_w_mjj) ;
  fitresult = DijetMassBkgTmp0->fitTo(*pseudoDataBinned, RooFit::Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), RooFit::SumW2Error(kTRUE), RooFit::Save(kTRUE));
  fitresult->Print("V");
  w->import(*DijetMassBkgTmp0);

  //*****************plot histogram rebinned **************************
  //debug
  //  const int nMassBins = 103;
//  double massBoundaries[nMassBins+1] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000};
//  TH1F* h_data_reBinned_mjj = (TH1F*)h_pseudodata_binned_mjj->Clone("h_data_reBinned_mjj");
//  h_data_reBinned_mjj->Rebin(nMassBins,"h_data_reBinned_mjj",massBoundaries);
//  int nbins_rebin = h_data_reBinned_mjj->GetNbinsX();
//  double Xmin_rebin = h_data_reBinned_mjj->GetBinLowEdge(1);
//  double Xmax_rebin = h_data_reBinned_mjj->GetBinLowEdge(nbins_rebin+1);
//
//  TH1F* h_data_reBinned_w_mjj = new TH1F("h_data_reBinned_w_mjj","",nbins_rebin,Xmin_rebin,Xmax_rebin);
//  for(int i=1; i<=h_data_reBinned_mjj->GetNbinsX(); i++){
//    h_data_reBinned_w_mjj->SetBinContent(i,h_data_reBinned_mjj->GetBinContent(i)/h_data_reBinned_w_mjj->GetBinWidth(i));
//    h_data_reBinned_w_mjj->SetBinError(i,TMath::Sqrt(h_data_reBinned_mjj->GetBinContent(i)/h_data_reBinned_w_mjj->GetBinWidth(i)));
//    cout << "new bin " << i << " content : " << h_pseudodata_binned_w_mjj->GetBinContent(i) <<endl;
//  }  
//  
  //RooBinning abins(nMassBins,massBoundaries,"mjj_binning") ;
  //RooRealVar* mjj_binned = new RooRealVar("mjj_binned","",1,14000);
  //mjj_binned->setBinning(abins) ;
  //RooDataHist* pseudodata_reBinned_w_mjj = new RooDataHist(TString::Format("pseudodata_w_cat0"),"",RooArgSet(*mjj_binned),h_data_reBinned_w_mjj) ;
  //RooHist* roo_h_pseudodata_reBinned_w_mjj = new RooHist(h_data_reBinned_w_mjj);


  //RooDataSet* pseudoDataset = DijetMassBkgTmp0->generate(*mjj,RooRandom::randomGenerator()->Poisson(275814));
//
  //************************************************
  // Plot DijetMass background fit results per categories 
  TCanvas* ctmp = new TCanvas("ctmp","DijetMass Background Categories",0,0,500,500);
  //Int_t nBinsMass(103);
  //RooRealVar* mjj_binned = (RooRealVar*)(dataBinned->get()->find("mjj")); 
  plotDijetMassBkg = mjj->frame(RooFit::Range(minMassFit,maxMassFit));//(minMassFit, maxMassFit,nBinsMass);
  
  //data->plotOn(plotDijetMassBkg,RooFit::MarkerColor(kPink), RooFit::Binning(445,1100.,10000.)/*,RooFit::Invisible()*/);    
  //dataBinned_1GeV->plotOn(plotDijetMassBkg);      
  //dataBinned->plotOn(plotDijetMassBkg,RooFit::MarkerColor(kPink));     
  //pseudoDataBinned->plotOn(plotDijetMassBkg);
  pseudoDataBinned_w->plotOn(plotDijetMassBkg, RooFit::Range(minMassFit,maxMassFit));
  //pseudodata_reBinned_w_mjj->plotOn(plotDijetMassBkg,RooFit::Range(minMassFit,maxMassFit));
  //plotDijetMassBkg->addPlotable(roo_h_pseudodata_reBinned_w_mjj,"P");

//  //giulia : example to adapt
//  //generate pseudodataset with correct number of events to preserve the shape
//  const int nMassBins = 103;
//
//  double massBoundaries[nMassBins+1] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000};
//
// 
//  ///////////////////
//  RooBinning abins(nMassBins,massBoundaries,"mjj_binning") ;
//
//  mjj->setBinning(abins) ;
//  RooDataHist* pseudoDataBinned = new RooDataHist(TString::Format("pseudodata_cat0"),"",RooArgSet(*mjj),*pseudoDataset) ;
//  //pseudoDataBinned->plotOn(plotDijetMassBkg);
//
//
//  ///////////
//  RooHistPdf* DijetMassBkgTmp0_hist = new RooHistPdf(TString::Format("DijetMassBkg_hist_cat%d",c),"",RooArgSet(*mjj),*dataBinned_1GeV); 
// 
//  mjj->setBinning(oneGeVbins) ; 
//  RooDataSet* pseudoData_2 = DijetMassBkgTmp0_hist->generate(*mjj,RooRandom::randomGenerator()->Poisson(275814));
//  //mjj->setBinning(abins) ;
//  //RooDataHist* pseudoDataBinned_2 = new RooDataHist(TString::Format("pseudodata_cat0"),"",RooArgSet(*mjj),*pseudoData_2) ;
//  RooFitResult*  fitresult_2 = DijetMassBkgTmp0->fitTo(*pseudoData_2, RooFit::Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), RooFit::SumW2Error(kTRUE), RooFit::Save(kTRUE));
//
//  fitresult_2->Print();

//  TH1F* h_pseudodata_binned_2 = (TH1F*)pseudoDataBinned_2->createHistogram("mjj");
//  int nbins = h_pseudodata_binned_2->GetNbinsX();   
//  double Xmin = h_pseudodata_binned_2->GetBinLowEdge(1);
//  double Xmax = h_pseudodata_binned_2->GetBinLowEdge(nbins+1);
//  TH1F* h_pseudodata_binned_w_2 = new TH1F("h_pseudodata_binned_w_2","", nbins, Xmin, Xmax);
//  for(int i=1; i<=h_pseudodata_binned_2->GetNbinsX(); i++){
//    h_pseudodata_binned_w_2->SetBinContent(i,h_pseudodata_binned_2->GetBinContent(i)/h_pseudodata_binned_2->GetBinWidth(i));
//  }
//
//  RooDataHist* pseudoDataBinned_w_2 = new RooDataHist(TString::Format("pseudodata_w_2_cat0"),"",RooArgSet(*mjj),h_pseudodata_binned_w_2) ;
//pseudoDataBinned_w_2->plotOn(plotDijetMassBkg);
 // pseudoData_2->plotOn(plotDijetMassBkg);
  DijetMassBkgTmp0->plotOn(plotDijetMassBkg,RooFit::LineColor(kBlue),RooFit::Range(minMassFit,maxMassFit)/*,RooFit::NormRange("bkg range")*/); 
  double chi2 = plotDijetMassBkg->chiSquare(3);
  Int_t ndof = (Int_t)plotDijetMassBkg->GetNbinsX() -3;
  std::cout<<"------> "<< ndof<<std::endl;
  double prob = TMath::Prob(chi2*ndof, ndof);
  std::cout<<prob<<std::endl;
 
  //blind = true;
  /* 
     if( blind ) {

     RooDataSet* data_down = (RooDataSet*) data.reduce(*w->var("Dijet_MassW"),"Dijet_MassW < 178.");
      RooDataSet* data_up = (RooDataSet*) data.reduce(*w->var("Dijet_MassW"),"Dijet_MassW >402");

      data_up->plotOn(plotDijetMassBkg);    
      data_down->plotOn(plotDijetMassBkg); 


   
    } else {
      data->plotOn(plotDijetMassBkg);    
      } 
       
    */

  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist* hpull = plotDijetMassBkg->pullHist() ;

  // Create a new frame to draw the pull distribution and add the distribution to the frame
  RooPlot* pull_frame = mjj->frame(RooFit::Title("Pull Distribution"),RooFit::Range(minMassFit,maxMassFit)) ;
  pull_frame->addPlotable(hpull,"P");

  ctmp->Divide(1,2);
  TPad *p1 = (TPad*)ctmp->cd(1) ; 
  gPad->SetLeftMargin(0.15) ; 
  
  plotDijetMassBkg->GetYaxis()->SetTitleOffset(1.6) ; 
  plotDijetMassBkg->GetXaxis()->SetTitle("m(jj) [GeV]");
  plotDijetMassBkg->SetAxisRange(0.001,plotDijetMassBkg->GetMaximum()*1.5,"Y");
  plotDijetMassBkg->Draw();  

  TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Category %d",c), "brNDC");
  //legdata->AddEntry(plotDijetMassBkg->getObject(3),"Data","LPE");
  legdata->AddEntry(plotDijetMassBkg->getObject(2),"PseudoData","LPE");
  legdata->AddEntry(plotDijetMassBkg->getObject(1),"Parametric Model: DiJet","L");
  legdata->SetTextSize(0.035);
  legdata->SetTextFont(42);
  // legdata->SetTextAlign(31);
  legdata->SetBorderSize(0);
  legdata->SetFillStyle(0);
  legdata->Draw("same");

  //TPaveText* label_cms = get_labelCMS(0, "2014", false);
  //TPaveText* label_sqrt = get_labelSqrt(0);
  //label_cms->Draw("same");
  //label_sqrt->Draw("same");

  //write down the chi2 of the fit on the
 
  TPaveText* label_chi2 = new TPaveText(0.5524194,0.6419492,0.796371,0.7690678, "brNDC");
  label_chi2->SetFillColor(kWhite);
  label_chi2->SetTextSize(0.035);
  label_chi2->SetTextFont(42);
  label_chi2->SetTextAlign(31); // align right
  label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", chi2));
  label_chi2->AddText(TString::Format("Chi square Prob = %.3f", prob));
  label_chi2->Draw("same");

  ctmp->cd(2) ; 
  gPad->SetLeftMargin(0.15) ; 
  pull_frame->GetYaxis()->SetTitleOffset(1.6) ; 
  pull_frame->Draw() ;


  //********************************************************************************// 
  //to be implemented
  //if (dobands) {
  //   
  //   RooAbsPdf *cpdf; cpdf = DijetMassBkgTmp0;
  //   TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
  //   TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
  //   
  //   RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
  //   nlim->removeRange();
  //   
  //   RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotDijetMassBkg->getObject(1));
  //   
  //   for (int i=1; i<(plotDijetMassBkg->GetXaxis()->GetNbins()+1); ++i) {
  //     double lowedge = plotDijetMassBkg->GetXaxis()->GetBinLowEdge(i);
  //     double upedge  = plotDijetMassBkg->GetXaxis()->GetBinUpEdge(i);
  //     double center  = plotDijetMassBkg->GetXaxis()->GetBinCenter(i);
  //     
  //     double nombkg = nomcurve->interpolate(center);
  //     nlim->setVal(nombkg);
  //     DijetMass->setRange("errRange",lowedge,upedge);
  //     RooAbsPdf *epdf = 0;
  //     epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
  //     
  //     RooAbsReal *nll = epdf->createNLL(*(data),RooFit::Extended());
  //     RooMinimizer minim(*nll);
  //     minim.setStrategy(0);
  //     double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
  //     double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
  //     
  //     minim.migrad();
  //     minim.minos(*nlim);
  //     // printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
  //     
  //     onesigma->SetPoint(i-1,center,nombkg);
  //     onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
  //     
  //     minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
  //     // eventually if cl = 0.95 this is the usual 1.92!      
  //     
  //     minim.migrad();
  //     minim.minos(*nlim);
  //     
  //     twosigma->SetPoint(i-1,center,nombkg);
  //     twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
  //     
  //     delete nll;
  //     delete epdf;
  //   }
  //   
  //   DijetMass->setRange("errRange",minMassFit,maxMassFit);
  //   
  //   twosigma->SetLineColor(kGreen);
  //   twosigma->SetFillColor(kGreen);
  //   twosigma->SetMarkerColor(kGreen);
  //   twosigma->Draw("L3 SAME");
  //   
  //   onesigma->SetLineColor(kYellow);
  //   onesigma->SetFillColor(kYellow);
  //   onesigma->SetMarkerColor(kYellow);
  //   onesigma->Draw("L3 SAME");
  //   
  //   plotDijetMassBkg->Draw("SAME"); 
  // }

  int massI(mass);
  ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_M%d.png",c,massI));
  ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_M%d.pdf",c,massI));
  ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_M%d.root",c,massI));

  
  p1->SetLogy();
  plotDijetMassBkg->SetAxisRange(0.001,plotDijetMassBkg->GetMaximum()*1.5,"Y");
  ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_LOG_M%d.png",c,massI));
  ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_LOG_M%d.pdf",c,massI));
  ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_LOG_M%d.root",c,massI));

  cout << "*************************************" << endl;
  cout << "return result:" << endl; 
  fitresult->Print("V");

return fitresult;
}
