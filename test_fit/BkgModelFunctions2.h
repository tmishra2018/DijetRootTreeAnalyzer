RooFitResult* BkgModelFitExpPARFunc(RooWorkspace* w, Bool_t dobands, Float_t mass,Int_t c,  bool blind) {

  Int_t ncat = NCAT;

  RooDataSet* data;
 
  RooFitResult* fitresult;


  RooPlot* plotPhotonsMassBkg;

  Float_t minMassFit, maxMassFit;
 
    minMassFit = MINmass;
    maxMassFit = MAXmass;
  
  // Fit data with background pdf for data limit
  RooRealVar* PhotonsMass = w->var("PhotonsMass");  
  PhotonsMass->setUnit("GeV");
 
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  
  
    data = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    PhotonsMass->setRange("bkg range", MINmass, MAXmass);
    
    // fit con expol 
    RooFormulaVar *p1mod= new RooFormulaVar(TString::Format("par1ExpPAR_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_ExpPAR1_cat%d",c)));
    RooFormulaVar *p2mod= new RooFormulaVar(TString::Format("par2ExpPAR_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_ExpPAR2_cat%d",c)));
    
    RooAbsPdf* PhotonsMassBkg = new RooGenericPdf(TString::Format("PhotonsMassBkg_cat%d",c), "exp(-@1*@0)*pow(@0, @2)", RooArgList(*PhotonsMass, *p1mod, *p2mod));
    
    fitresult = PhotonsMassBkg->fitTo(*data, Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), SumW2Error(kTRUE), Save(kTRUE));   
   
    std::cout<<TString::Format("******************************** Background Fit results mass %f cat %d ***********************************", mass, c)<<std::endl;
    fitresult->Print("V");
    w->import(*PhotonsMassBkg);
 

    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background ",0,0,800,800);
    Int_t nBinsMass(60);
    plotPhotonsMassBkg = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
  
    data->plotOn(plotPhotonsMassBkg,RooFit::Invisible());    
   
    PhotonsMassBkg->plotOn(plotPhotonsMassBkg,LineColor(kBlue),Range(minMassFit,maxMassFit),NormRange("bkg range")); 
   
    double chi2 = plotPhotonsMassBkg->chiSquare(2);
    Int_t ndof = nBinsMass-2;
    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;


    blind=false;
    if( blind ) {
     
      RooDataSet* data_down = (RooDataSet*) data.reduce(*w->var("PhotonsMass"),"PhotonsMass < 327.666 ");
      RooDataSet* data_up = (RooDataSet*) data.reduce(*w->var("PhotonsMass"),"PhotonsMass >650");

      data_up->plotOn(plotPhotonsMassBkg);    
      data_down->plotOn(plotPhotonsMassBkg); 


   
    } else {
      data->plotOn(plotPhotonsMassBkg);    
      } 
      TH1F* h_data = new TH1F("h_data","h_data", 60,minMassFit, maxMassFit);
    TH1F* h_pdf = new TH1F("h_pdf","h_pdf", 60,minMassFit, maxMassFit);
    h_data = (TH1F*) data->createHistogram("PhotonsMass", 60, 0, 0);
    h_pdf =  (TH1F*) PhotonsMassBkg->createHistogram("PhotonsMass",60);
    h_pdf ->Scale(h_data->Integral()/h_pdf->Integral());

    //-------pad 1-------//
    //TPad * pad1 = new TPad("pad1", "pad1",0.01256281,0.1304945,0.5741206,1);  
    TPad * pad1 = new TPad("pad1", "pad1",0., 0., 1., 1.);  
      
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();
    pad1->Range(154.1111,-5.730293,650.7778,4.501062);
    pad1->SetLeftMargin(0.1789709);
    pad1->SetRightMargin(0.01565995);
    pad1->SetTopMargin(0.04897314);
    pad1->SetBottomMargin(0.1691167);
   
   
    plotPhotonsMassBkg->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    plotPhotonsMassBkg->GetYaxis()->SetTitle("Events/(6.7)");
    plotPhotonsMassBkg->SetAxisRange(0.001,plotPhotonsMassBkg->GetMaximum()*1.5,"Y");
    plotPhotonsMassBkg->Draw();  

    TLegend *legdata = new TLegend(0.5334677,0.680339,0.8245968,0.8958475, TString::Format("Category %d",c), "brNDC");
    if(mass<450.)legdata = new TLegend(0.2334677,0.300339,0.5645968,0.4958475, TString::Format("Category %d",c), "brNDC");
    legdata->AddEntry(plotPhotonsMassBkg->getObject(2),"Data","LPE");
    legdata->AddEntry(plotPhotonsMassBkg->getObject(1),"Fit Model","L");
  
  

    dobands=false;
    //********************************************************************************/
    /*    if (dobands) {

      RooAbsPdf *cpdf; cpdf = PhotonsMassBkg;
      TGraphAsymmErrors onesigma;
      TGraphAsymmErrors twosigma;
      
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotPhotonsMassBkg->getObject(1));
      
      double el1;
      double eh1;
      double el2;
      double eh2;
  
      int j = 0;
      for (int i=1; i<(plotPhotonsMassBkg->GetXaxis()->GetNbins()+1); ++i) {
	double lowedge = plotPhotonsMassBkg->GetXaxis()->GetBinLowEdge(i);
	double upedge  = plotPhotonsMassBkg->GetXaxis()->GetBinUpEdge(i);
	double center  = plotPhotonsMassBkg->GetXaxis()->GetBinCenter(i);
	
	double nombkg = nomcurve->interpolate(center);


	nlim->setVal(nombkg);
	PhotonsMass->setRange("errRange",lowedge,upedge);
	RooAbsPdf *epdf = 0;
	epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
	
	RooAbsReal *nll = epdf->createNLL(*(data),Extended());
	RooMinimizer minim(*nll);
	minim.setStrategy(0);
	double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
	double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	
	minim.migrad();
	minim.minos(*nlim);
	printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());


	el1 = nlim->getErrorLo();
	eh1= nlim->getErrorHi();
	//	std::cout<<"-----------------------------------------------------------------> "<< minim.minos(*nlim)<<std::endl;
	minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
	minim.migrad();
	minim.minos(*nlim);
	el2 = nlim->getErrorLo();
	eh2= nlim->getErrorHi();


	delete nll;
	delete epdf;
	if( el1 != 0. && eh1 != 0. && el2 != 0. && eh2 != 0. &&  el1 != 1. && eh1 != 1.  && el2 != 1.  && eh2 != 1. ) {
	onesigma.SetPoint(j,center,nombkg);
	twosigma.SetPoint(j,center,nombkg);
	onesigma.SetPointError(j,0.,0.,-el1,eh1);
	twosigma.SetPointError(j,0.,0.,-el2,eh2);
	j++;
	}

	}

      PhotonsMass->setRange("errRange",minMassFit,maxMassFit);
      
      twosigma.SetLineColor(kGreen);
      twosigma.SetFillColor(kGreen);
      twosigma.SetMarkerColor(kGreen);
      twosigma.Draw("C3 SAME");
  
      onesigma.SetLineColor(kYellow);
      onesigma.SetFillColor(kYellow);
      onesigma.SetMarkerColor(kYellow);
      onesigma.Draw("C3 SAME");
   
      legdata->AddEntry(&onesigma, "#pm 1 #sigma", "F" );
      legdata->AddEntry(&twosigma, "#pm 2 #sigma","F" );
      plotPhotonsMassBkg->Draw("SAME"); 
     
      }
    */
    legdata->SetTextSize(0.035);
    legdata->SetTextFont(42);
    legdata->SetBorderSize(0);
    legdata->SetFillStyle(0);
    legdata->Draw("same");


    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    //    label_cms->Draw("same");
    //label_sqrt->Draw("same");
    int iPos=11 ;
    CMS_lumi( pad1,false,iPos );

    //write down the chi2 of the fit on the
      
    TPaveText* label_chi2= new TPaveText(0.5744355,0.6050847,0.80583871,0.65822034,"brNDC");
    if(mass<450. || c >1)label_chi2 = new TPaveText(0.4479505,0.1861729,0.6095114,0.2933837,"brNDC");
    label_chi2->SetFillColor(kWhite);
    label_chi2->SetTextSize(0.035);
    label_chi2->SetTextFont(42);
    label_chi2->SetTextAlign(31); // align right
    label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", chi2));
    label_chi2->AddText(TString::Format("Chi square Prob = %.3f", prob));
    label_chi2->Draw("same");
    

    ctmp->cd();
    /*  //-------pad 2------//
    TPad * pad2 = new TPad("pad2", "pad2",0.01507538,0.01510989,0.5854271,0.260989);
    pad2->SetGrid();
    
    pad2->Draw();
    pad2->Range(157.9171,-7.736842,659.5746,3.568421);
    pad2->SetLeftMargin(0.1696035);
    pad2->SetRightMargin(0.03303965);
    pad2->SetTopMargin(0.05027933);
    pad2->SetBottomMargin(0.4189944);
   
    

   
    
    pad2->cd();
    h_data->Sumw2();
    h_pdf->Sumw2();
    TH1F* pull = new TH1F("pull", "pull", 20., -3., 3.);
    TH1F* h_data1 = h_data->Clone();
    h_data->Add(h_pdf, -1);
    double ndata;
    double nfit;
    double errdata;
    double errfit;

    for(int i = 1; i< h_data->GetNbinsX()+1;i++){
      ndata=h_data1->GetBinContent(i);
      errdata=h_data1->GetBinError(i);
      nfit=h_pdf->GetBinContent(i);
      errfit=h_pdf->GetBinError(i);
      h_data->SetBinContent(i,h_data->GetBinContent(i)/errfit);
      pull->Fill(h_data->GetBinContent(i));
    
    }

    Double_t mean = pull->GetMean();
    Double_t sigma = pull->GetRMS();
    Double_t meanErr =pull->GetMeanError();
    Double_t sigmaErr = pull->GetRMSError();

    std::cout<<"MEAN: "<<pull->GetMean()<<"  RMS: "<<pull->GetRMS()<<std::endl;
    //    for(int i = 0; i< h_data->GetNbinsX();i++) h_data->SetBinError(i,h_data1->GetBinError(i));
    h_data->GetYaxis()->SetRangeUser(-3., 3.);
    h_data->GetYaxis()->SetNdivisions(505);
    h_data->SetMarkerSize(0.4);
    h_data->GetXaxis()->SetTitle("m_{#gamma #gamma}");
    h_data->GetXaxis()->SetTitleSize(0.15);
    h_data->GetXaxis()->SetLabelSize(0.1);
    h_data->GetYaxis()->SetLabelSize(0.1);
    h_data->GetYaxis()->SetTitleSize(0.13);
    h_data->GetYaxis()->SetTitle("#frac{N_{obs}-N_{exp}}{#sqrt{N_{exp}}}");
    h_data->GetYaxis()->SetTitleOffset(0.45);
    h_data->GetXaxis()->SetTitleOffset(0.8);
    for(int i=0; i<h_data->GetNbinsX()+1;i++) h_data->SetBinError(i, 0.);
    h_data->Draw("P");
    TH1F* h_dataCopy = h_data->Clone();
    for (int i = 0; i< h_dataCopy->GetNbinsX()+1;i++) if (h_dataCopy->GetBinContent(i)==0) h_dataCopy->SetBinContent(i, 1);
   
    h_dataCopy->Add(h_dataCopy, -1);
    for (int i = 0; i< h_dataCopy->GetNbinsX()+1;i++) h_dataCopy->SetBinError(i, 1);
    h_dataCopy->SetLineColor(kAzure-2);
    h_dataCopy->SetFillColor(kAzure-2);
    h_dataCopy->SetFillStyle(3002);
    h_dataCopy->SetMarkerSize(0.);
    h_dataCopy->Draw("HISTE23same");
    h_data->Draw("PSAME");


  //-------pad 3------//
    
    ctmp->cd();
    TPad * pad3 = new TPad("pad3", "pad3",0.5703518,0.01785714,0.7738693,0.2554945);
    pad3->SetGrid();
    
    pad3->Draw();
    pad3->cd();
    pad3->Range(-1.052224,-7.652275,12.89804,3.257984);
      pad3->SetLeftMargin(0.07542683);
   pad3->SetRightMargin(0.1718992);
   pad3->SetTopMargin(0.02364605);
   pad3->SetBottomMargin(0.4264129);
   
    pad3->SetFrameFillStyle(0);
    pad3->SetFrameBorderMode(0);
    pad3->SetFrameFillStyle(0);
    pad3->SetFrameBorderMode(0);

    pull->SetFillColor(kBlack);
    pull->Draw("hbar");

   
    ctmp->cd();
    TPaveText* label2 = new TPaveText(0.5745968,0.1779661,0.7681452,0.3919492, "brNDC" );
    label2->SetFillColor(kWhite);
    label2->SetBorderSize(0.);
    label2->SetTextSize(0.0213922);
    label2->SetTextAlign(11);
    label2->SetTextFont(42);
    label2->AddText(TString::Format("#splitline{Mean: %.3f #pm %.3f }{RMS: %.3f #pm %.3f}", mean,meanErr,sigma, sigmaErr));


    label2->Draw("same");


    pad3->Draw();

    */

    int massI(mass);
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPAR_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPAR_M%d.pdf",c,massI));

    ctmp->SetLogy();
    if(mass<650)plotPhotonsMassBkg->SetAxisRange(0.5,plotPhotonsMassBkg->GetMaximum()*10,"Y");
    else plotPhotonsMassBkg->SetAxisRange(0.5,100,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.pdf",c,massI));
   
    // if(mass==350){
    //   ctmp->SaveAs("~/www/plotsNota/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.png",c,massI));
    //   ctmp->SaveAs("~/www/plotsNota/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.pdf",c,massI));
    //   ctmp->SaveAs("~/www/plotsPAS/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.png",c,massI));
    //   ctmp->SaveAs("~/www/plotsPAS/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.pdf",c,massI));
    //   ctmp->SaveAs("~/www/plotsPAS/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.C",c,massI));
    //   ctmp->SaveAs("~/www/plotsPAS/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.root",c,massI));
    // }
    

  RooFitResult* r;
    
  return r;
}




RooFitResult* BkgModelFitDiJetFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, bool blind) {

  Int_t ncat = NCAT;
  std::cout<<"isBlind: "<<blind<<std::endl;
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data[NCAT];
 
  RooFitResult* fitresult[NCAT];;
  RooPlot* plotPhotonsMassBkg[NCAT];

  // dobands and dosignal
  RooDataSet* signal[NCAT];

  RooAbsPdf* PhotonsMassSig[NCAT];
  
  Float_t minMassFit, maxMassFit;
  minMassFit = MINmass;
  maxMassFit = MAXmass;
   
  // Fit data with background pdf for data limit
  RooRealVar* PhotonsMass = w->var("PhotonsMass");  
  PhotonsMass->setUnit("GeV");
 
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  
  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    
  
    // fit a la dijets
    RooFormulaVar *p0mod = new RooFormulaVar(TString::Format("par0DiJet_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_norm_cat%d",c)));
    RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("par1DiJet_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope1_cat%d",c)));
    RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("par2DiJet_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope2_cat%d",c)));
    RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("par3DiJet_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope3_cat%d",c)));
    PhotonsMass->setRange("bkg range", MINmass,MAXmass);
    RooFormulaVar *x     = new RooFormulaVar(TString::Format("xDiJet_cat%d",c),"","@0/8000.",*w->var("PhotonsMass"));

   
    RooAbsPdf* PhotonsMassBkgTmp0 = new RooGenericPdf(TString::Format("PhotonsMassBkg_cat%d",c), "pow(1-@0, @2)/pow(@0, @1+@3*log(@0))", RooArgList(*x, *p1mod, *p2mod,*p3mod));
   

    fitresult[c] = PhotonsMassBkgTmp0->fitTo(*data[c], Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), SumW2Error(kTRUE), Save(kTRUE));
    w->import(*PhotonsMassBkgTmp0);
   
    std::cout<<TString::Format("******************************** Background Fit results mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresult[c]->Print("V");
   

    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
    Int_t nBinsMass(60);
    plotPhotonsMassBkg[c] = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
  
    data[c]->plotOn(plotPhotonsMassBkg[c],RooFit::Invisible());    
   
    PhotonsMassBkgTmp0->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue),Range(minMassFit,maxMassFit),NormRange("bkg range")); 
    double chi2 = plotPhotonsMassBkg[c]->chiSquare(3);
    Int_t ndof = nBinsMass-3;
    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;
    blind = true;
    if( blind ) {
     
      RooDataSet* data_down = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass < 178.");
      RooDataSet* data_up = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass >402");

      data_up->plotOn(plotPhotonsMassBkg[c]);    
      data_down->plotOn(plotPhotonsMassBkg[c]); 


   
    } else {
      data[c]->plotOn(plotPhotonsMassBkg[c]);    
      } 
       
   
    plotPhotonsMassBkg[c]->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");
    plotPhotonsMassBkg[c]->SetAxisRange(0.001,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    plotPhotonsMassBkg[c]->Draw();  

    TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Category %d",c), "brNDC");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(2),"Data","LPE");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(1),"Parametric Model: DiJet","L");
    legdata->SetTextSize(0.035);
    legdata->SetTextFont(42);
    // legdata->SetTextAlign(31);
    legdata->SetBorderSize(0);
    legdata->SetFillStyle(0);
    legdata->Draw("same");

    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");

    //write down the chi2 of the fit on the
 
    TPaveText* label_chi2 = new TPaveText(0.5524194,0.6419492,0.796371,0.7690678, "brNDC");
    label_chi2->SetFillColor(kWhite);
    label_chi2->SetTextSize(0.035);
    label_chi2->SetTextFont(42);
    label_chi2->SetTextAlign(31); // align right
    label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", chi2));
    label_chi2->AddText(TString::Format("Chi square Prob = %.3f", prob));
    label_chi2->Draw("same");

    

    //********************************************************************************//
    if (dobands) {

      RooAbsPdf *cpdf; cpdf = PhotonsMassBkgTmp0;
      TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
      TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
      
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotPhotonsMassBkg[c]->getObject(1));
      
      for (int i=1; i<(plotPhotonsMassBkg[c]->GetXaxis()->GetNbins()+1); ++i) {
	double lowedge = plotPhotonsMassBkg[c]->GetXaxis()->GetBinLowEdge(i);
	double upedge  = plotPhotonsMassBkg[c]->GetXaxis()->GetBinUpEdge(i);
	double center  = plotPhotonsMassBkg[c]->GetXaxis()->GetBinCenter(i);
	
	double nombkg = nomcurve->interpolate(center);
	nlim->setVal(nombkg);
	PhotonsMass->setRange("errRange",lowedge,upedge);
	RooAbsPdf *epdf = 0;
	epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
	
	RooAbsReal *nll = epdf->createNLL(*(data[c]),Extended());
	RooMinimizer minim(*nll);
	minim.setStrategy(0);
	double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
	double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	
	minim.migrad();
	minim.minos(*nlim);
	// printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
	
	onesigma->SetPoint(i-1,center,nombkg);
	onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
	// eventually if cl = 0.95 this is the usual 1.92!      
	
	minim.migrad();
	minim.minos(*nlim);
	
	twosigma->SetPoint(i-1,center,nombkg);
	twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	delete nll;
	delete epdf;
      }

      PhotonsMass->setRange("errRange",minMassFit,maxMassFit);
      
      twosigma->SetLineColor(kGreen);
      twosigma->SetFillColor(kGreen);
      twosigma->SetMarkerColor(kGreen);
      twosigma->Draw("L3 SAME");
      
      onesigma->SetLineColor(kYellow);
      onesigma->SetFillColor(kYellow);
      onesigma->SetMarkerColor(kYellow);
      onesigma->Draw("L3 SAME");
      
      plotPhotonsMassBkg[c]->Draw("SAME"); 
    }

    int massI(mass);
    ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_M%d.png",c,massI));
    ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_M%d.pdf",c,massI));
    ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_M%d.root",c,massI));

    ctmp->SetLogy();
    plotPhotonsMassBkg[c]->SetAxisRange(1.3,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_LOG_M%d.png",c,massI));
    ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_LOG_M%d.pdf",c,massI));
    ctmp->SaveAs("plots/prelimBkg"+TString::Format("_cat%d_DIJET_LOG_M%d.root",c,massI));

  }



  return fitresult;
}


RooFitResult* BkgModelFitDiJetEXPFunc(RooWorkspace* w, Bool_t dobands, Float_t mass,Int_t c, bool blind) {

  Int_t ncat = NCAT;
  std::cout<<"isBlind: "<<blind<<std::endl;
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data;
 
  RooFitResult* fitresult;;
  
  RooPlot* plotPhotonsMassBkg;

  // dobands and dosignal
  RooDataSet* signal;

  RooAbsPdf* PhotonsMassSig;
  
  Float_t minMassFit, maxMassFit;

    minMassFit = MINmass;
    maxMassFit = MAXmass;

  
  // Fit data with background pdf for data limit
  RooRealVar* PhotonsMass = w->var("PhotonsMass");  
  PhotonsMass->setUnit("GeV");
 
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  
 
    data = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    
  
    // fit a la dijets
    
    RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("par1DiJetEXP_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope1_3_cat%d",c)));
    RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("par2DiJetEXP_cat%d",cv),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope2_3_cat%d",c)));
    RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("par3DiJetEXP_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope3_3_cat%d",c)));
    RooFormulaVar *exp1 = new RooFormulaVar(TString::Format("expDiJetEXP_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_exp1DiJetEXP_cat%d",c)));
    RooFormulaVar *pFrac1 = new RooFormulaVar(TString::Format("fracDiJetEXP_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_fracDiJetEXP_cat%d",c)));
   
   
    RooGenericPdf* PhotonsMassBkgTmp0DiJet = new RooGenericPdf(TString::Format("PhotonsMassBkg_DIJETE_truth_cat%d",c), "pow(1-@0/8000., @2)/pow(@0/8000., @1+@3*log(@0/8000.))", RooArgList(*w->var("PhotonsMass"), *p1mod, *p2mod,*p3mod));
  
    RooExponential* PhotonsMassBkgTmp0Exp = new RooExponential(TString::Format("PhotonsMassBkg_EXP_truth_cat%d",c),"", *w->var("PhotonsMass"),  *exp1);
    
   
    RooAddPdf* PhotonsMassBkgTmpAdd = new RooAddPdf(TString::Format("PhotonsMassBkg_cat%d",c),TString::Format("PhotonsMassBkg_cat%d",c) , RooArgList(*PhotonsMassBkgTmp0DiJet, *PhotonsMassBkgTmp0Exp), RooArgList(*pFrac1));
    
    fitresult = PhotonsMassBkgTmpAdd->fitTo(*data,RooFit::FitOptions("MHTR"), Save(kTRUE));//RooFit::FitOptions("MHTER"), Range(minMassFit,maxMassFit),    
    w->import(*PhotonsMassBkgTmpAdd);
  
    std::cout<<TString::Format("******************************** Background DiJetEXP Fit results mass %f cat %d ***********************************", mass, c)<<std::endl;
    fitresult->Print("V");
   

    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
    Int_t nBinsMass(60);
    plotPhotonsMassBkg = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
  
    data->plotOn(plotPhotonsMassBkg,RooFit::Invisible());    
   
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg,LineColor(kBlue)); 
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg,Components(TString::Format("PhotonsMassBkg_EXP_truth_cat%d",c)),LineColor(kViolet),LineStyle(kDashed)); 
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg,Components(TString::Format("PhotonsMassBkg_DIJETE_truth_cat%d",c)),LineColor(kOrange),LineStyle(kDashed));   

    double chi2 = plotPhotonsMassBkg->chiSquare(3);
    Int_t ndof = nBinsMass-5;
    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;
   
    if( blind ) {
 
      RooDataSet* data_down = (RooDataSet*) data.reduce(*w->var("PhotonsMass"),"PhotonsMass < 173.5");
      RooDataSet* data_up = (RooDataSet*) data.reduce(*w->var("PhotonsMass")," PhotonsMass >850");
      TH1F* h_up= new TH1F("h_up", "h_up",nBinsMass, 130, 1000);
      h_up->Sumw2();
      data_up->fillHistogram(h_up, RooArgList(*PhotonsMass));
      TH1F* h_down= new TH1F("h_down", "h_down",nBinsMass, 130, 1000);
      h_down->Sumw2();
      data_down->fillHistogram(h_down, RooArgList(*PhotonsMass));
   
    } else {
      data->plotOn(plotPhotonsMassBkg);    
      } 
       
   
    plotPhotonsMassBkg->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    plotPhotonsMassBkg->SetAxisRange(0.1,10000,"Y");
    plotPhotonsMassBkg->Draw();  
   if( blind ) {
       h_up->Draw("sameP");
       h_down->Draw("sameP");
     }
  
    TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Category %d",c), "brNDC");
    legdata->AddEntry(plotPhotonsMassBkg->getObject(4),"Data","LPE");
    legdata->AddEntry(plotPhotonsMassBkg->getObject(1),"Parametric Model: DiJetEXP","L");  
    legdata->AddEntry(plotPhotonsMassBkg->getObject(2),"Parametric Model: EXP","L");
    legdata->AddEntry(plotPhotonsMassBkg->getObject(3),"Parametric Model: DiJet","L");
  
    legdata->SetTextSize(0.035);
    legdata->SetTextFont(42);
    // legdata->SetTextAlign(31);
    legdata->SetBorderSize(0);
    legdata->SetFillStyle(0);
    legdata->Draw("same");

    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");

    //write down the chi2 of the fit on the
 
    TPaveText* label_chi2 = new TPaveText(0.5524194,0.6419492,0.796371,0.7690678, "brNDC");
    label_chi2->SetFillColor(kWhite);
    label_chi2->SetTextSize(0.035);
    label_chi2->SetTextFont(42);
    label_chi2->SetTextAlign(31); // align right
    label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", chi2));
    label_chi2->AddText(TString::Format("Chi square Prob = %.3f", prob));
    label_chi2->Draw("same");

    

    //********************************************************************************//
    if (dobands) {

      RooAbsPdf *cpdf; cpdf = PhotonsMassBkgTmp0;
      TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
      TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
      
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotPhotonsMassBkg[c]->getObject(1));
      
      for (int i=1; i<(plotPhotonsMassBkg->GetXaxis()->GetNbins()+1); ++i) {
	double lowedge = plotPhotonsMassBkg->GetXaxis()->GetBinLowEdge(i);
	double upedge  = plotPhotonsMassBkg->GetXaxis()->GetBinUpEdge(i);
	double center  = plotPhotonsMassBkg->GetXaxis()->GetBinCenter(i);
	
	double nombkg = nomcurve->interpolate(center);
	nlim->setVal(nombkg);
	PhotonsMass->setRange("errRange",lowedge,upedge);
	RooAbsPdf *epdf = 0;
	epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
	
	RooAbsReal *nll = epdf->createNLL(*(data),Extended());
	RooMinimizer minim(*nll);
	minim.setStrategy(0);
	double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
	double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	
	minim.migrad();
	minim.minos(*nlim);
	// printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
	
	onesigma->SetPoint(i-1,center,nombkg);
	onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
	// eventually if cl = 0.95 this is the usual 1.92!      
	
	minim.migrad();
	minim.minos(*nlim);
	
	twosigma->SetPoint(i-1,center,nombkg);
	twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	delete nll;
	delete epdf;
      }

      PhotonsMass->setRange("errRange",minMassFit,maxMassFit);
      
      twosigma->SetLineColor(kGreen);
      twosigma->SetFillColor(kGreen);
      twosigma->SetMarkerColor(kGreen);
      twosigma->Draw("L3 SAME");
      
      onesigma->SetLineColor(kYellow);
      onesigma->SetFillColor(kYellow);
      onesigma->SetMarkerColor(kYellow);
      onesigma->Draw("L3 SAME");
      
      plotPhotonsMassBkg->Draw("SAME"); 
    }

    int massI(mass);
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_M%d.root",c,massI));

    ctmp->SetLogy();
    //  plotPhotonsMassBkg->SetAxisRange(1.3,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_LOG_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_LOG_M%d.root",c,massI));

 

  return fitresult;
}






RooFitResult* BkgModelFitDiJetEXPOLFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, Int_t c, bool blind) {

  Int_t ncat = NCAT;
 
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data;
 
  RooFitResult* fitresult;;

  RooPlot* plotPhotonsMassBkg;

  // dobands and dosignal
  RooDataSet* signal;

  RooAbsPdf* PhotonsMassSig;
  
  Float_t minMassFit, maxMassFit;

    minMassFit = MINmass;
    maxMassFit = MAXmass;

  
  // Fit data with background pdf for data limit
  RooRealVar* PhotonsMass = w->var("PhotonsMass");  
  PhotonsMass->setUnit("GeV");
 
v  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  
  
    data = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    
  
    // fit a la dijets
    
    RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("par1DiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope1_4_cat%d",c)));
    RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("par2DiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope2_4_cat%d",c)));
    RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("par3DiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope3_4_cat%d",c)));
    RooFormulaVar *expol1 = new RooFormulaVar(TString::Format("expol1DiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol1_4_cat%d",c)));
    RooFormulaVar *expol2 = new RooFormulaVar(TString::Format("expol2DiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol2_4_cat%d",c)));
    RooFormulaVar *pFrac1 = new RooFormulaVar(TString::Format("fracDiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_fracDiJetEXPOL_cat%d",c)));
   
   
    RooGenericPdf* PhotonsMassBkgTmp0DiJet = new RooGenericPdf(TString::Format("PhotonsMassBkg_DiJetEx_cat%d",c), "pow(1-@0/8000., @2)/pow(@0/8000., @1+@3*log(@0/8000.))", RooArgList(*w->var("PhotonsMass"), *p1mod, *p2mod,*p3mod));
    RooGenericPdf* PhotonsMassBkgTmp0Expol = new RooGenericPdf(TString::Format("PhotonsMassBkg_ExpolDiJ_cat%d",c), "exp(-@0/(@1+@2*@0))", RooArgList(*w->var("PhotonsMass"), *expol1, *expol2));
    
   
    RooAddPdf* PhotonsMassBkgTmpAdd = new RooAddPdf(TString::Format("PhotonsMassBkg_cat%d",c),TString::Format("PhotonsMassBkg_cat%d",c) , RooArgList(*PhotonsMassBkgTmp0DiJet, *PhotonsMassBkgTmp0Expol), RooArgList(*pFrac1));
    
    fitresult = PhotonsMassBkgTmpAdd->fitTo(*data,RooFit::FitOptions("MHTR"), Save(kTRUE));//RooFit::FitOptions("MHTER"), Range(minMassFit,maxMassFit),    
    w->import(*PhotonsMassBkgTmpAdd);
  
    std::cout<<TString::Format("******************************** Background DiJetEXPOL Fit results mass %f cat %d ***********************************", mass, c)<<std::endl;
    fitresult->Print("V");
   

    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
    Int_t nBinsMass(60);
    plotPhotonsMassBkg = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
  
    data->plotOn(plotPhotonsMassBkg,RooFit::Invisible());    
   
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg,LineColor(kBlue)); 
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg,Components(TString::Format("PhotonsMassBkg_ExpolDiJ_cat%d",c)),LineColor(kViolet),LineStyle(kDashed)); 
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg,Components(TString::Format("PhotonsMassBkg_DiJetEx_cat%d",c)),LineColor(kOrange),LineStyle(kDashed));   

    double chi2 = plotPhotonsMassBkg->chiSquare(3);
    Int_t ndof = nBinsMass-6;
    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;
    blind = false;
    if( blind ) {
     
      RooDataSet* data_down = (RooDataSet*) data.reduce(*w->var("PhotonsMass"),"PhotonsMass < 178.");
      RooDataSet* data_up = (RooDataSet*) data.reduce(*w->var("PhotonsMass"),"PhotonsMass >402");

      data_up->plotOn(plotPhotonsMassBkg);    
      data_down->plotOn(plotPhotonsMassBkg); 


   
    } else {
      data->plotOn(plotPhotonsMassBkg);    
      } 
       
   
    plotPhotonsMassBkg->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    plotPhotonsMassBkg->SetAxisRange(0.001,plotPhotonsMassBkg->GetMaximum()*1.5,"Y");
    plotPhotonsMassBkg->Draw();  

    TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Category %d",c), "brNDC");
    legdata->AddEntry(plotPhotonsMassBkg->getObject(4),"Data","LPE");
    legdata->AddEntry(plotPhotonsMassBkg->getObject(1),"Parametric Model: DiJetEXPOL","L");  
    legdata->AddEntry(plotPhotonsMassBkg->getObject(2),"Parametric Model: Expol","L");
    legdata->AddEntry(plotPhotonsMassBkg->getObject(3),"Parametric Model: DiJet","L");
  
    legdata->SetTextSize(0.035);
    legdata->SetTextFont(42);
v    // legdata->SetTextAlign(31);
    legdata->SetBorderSize(0);
    legdata->SetFillStyle(0);
    legdata->Draw("same");

    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");

    //write down the chi2 of the fit on the
 
    TPaveText* label_chi2 = new TPaveText(0.5524194,0.6419492,0.796371,0.7690678, "brNDC");
    label_chi2->SetFillColor(kWhite);
    label_chi2->SetTextSize(0.035);
    label_chi2->SetTextFont(42);
    label_chi2->SetTextAlign(31); // align right
    label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", chi2));
    label_chi2->AddText(TString::Format("Chi square Prob = %.3f", prob));
    label_chi2->Draw("same");

    

    //********************************************************************************//
    if (dobands) {

      RooAbsPdf *cpdf; cpdf = PhotonsMassBkgTmp0;
      TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
      TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
      
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotPhotonsMassBkg->getObject(1));
      
      for (int i=1; i<(plotPhotonsMassBkg->GetXaxis()->GetNbins()+1); ++i) {
	double lowedge = plotPhotonsMassBkg->GetXaxis()->GetBinLowEdge(i);
	double upedge  = plotPhotonsMassBkg->GetXaxis()->GetBinUpEdge(i);
	double center  = plotPhotonsMassBkg->GetXaxis()->GetBinCenter(i);
	
	double nombkg = nomcurve->interpolate(center);
	nlim->setVal(nombkg);
	PhotonsMass->setRange("errRange",lowedge,upedge);
	RooAbsPdf *epdf = 0;
	epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
	
	RooAbsReal *nll = epdf->createNLL(*(data),Extended());
	RooMinimizer minim(*nll);
	minim.setStrategy(0);
	double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
	double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	
	minim.migrad();
	minim.minos(*nlim);
	// printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
	
	onesigma->SetPoint(i-1,center,nombkg);
	onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
	// eventually if cl = 0.95 this is the usual 1.92!      
	
	minim.migrad();
	minim.minos(*nlim);
	
	twosigma->SetPoint(i-1,center,nombkg);
	twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	delete nll;
	delete epdf;
      }

      PhotonsMass->setRange("errRange",minMassFit,maxMassFit);
      
      twosigma->SetLineColor(kGreen);
      twosigma->SetFillColor(kGreen);
      twosigma->SetMarkerColor(kGreen);
      twosigma->Draw("L3 SAME");
      
      onesigma->SetLineColor(kYellow);
      onesigma->SetFillColor(kYellow);
      onesigma->SetMarkerColor(kYellow);
      onesigma->Draw("L3 SAME");
      
      plotPhotonsMassBkg->Draw("SAME"); 
    }

    int massI(mass);
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_M%d.root",c,massI));

    ctmp->SetLogy();
    plotPhotonsMassBkg->SetAxisRange(1.3,plotPhotonsMassBkg->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_LOG_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_LOG_M%d.root",c,massI));

  



  return fitresult;
}


RooFitResult* BkgModelFitExpolFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, bool blind) {

  Int_t ncat = NCAT;
  std::cout<<"isBlind: "<<blind<<std::endl;
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data[NCAT];
 
  RooFitResult* fitresult[NCAT];;
  RooPlot* plotPhotonsMassBkg[NCAT];

  // dobands and dosignal
  RooDataSet* signal[NCAT];

  RooAbsPdf* PhotonsMassSig[NCAT];
  
  Float_t minMassFit, maxMassFit;
 
    minMassFit = MINmass;
    maxMassFit = MAXmass;
  
  // Fit data with background pdf for data limit
  RooRealVar* PhotonsMass = w->var("PhotonsMass");  
  PhotonsMass->setUnit("GeV");
 
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  
  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    
  
    // fit con expo pol
   
    RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("par1Expol_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol1_cat%d",c)));
    RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("par2Expol_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol2_cat%d",c)));
    RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("par3Expol_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol3_cat%d",c)));
 
    PhotonsMass->setRange("bkg range", MINmass,MAXmass);
    RooFormulaVar *x     = new RooFormulaVar(TString::Format("xExpol_cat%d",c),"","@0",*w->var("PhotonsMass"));

   
    RooAbsPdf* PhotonsMassBkgTmp0 = new RooGenericPdf(TString::Format("PhotonsMassBkg_cat%d",c), "exp(-@0*@0/(@1+@2*@0+@3*@0*@0))", RooArgList(*x, *p1mod, *p2mod,*p3mod));
   

    fitresult[c] = PhotonsMassBkgTmp0->fitTo(*data[c], Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), SumW2Error(kTRUE), Save(kTRUE));
    w->import(*PhotonsMassBkgTmp0);
   
    std::cout<<TString::Format("******************************** Background Fit results mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresult[c]->Print("V");
   

    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
    Int_t nBinsMass(60);
    plotPhotonsMassBkg[c] = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
  
    data[c]->plotOn(plotPhotonsMassBkg[c],RooFit::Invisible());    
   
    PhotonsMassBkgTmp0->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue),Range(minMassFit,maxMassFit),NormRange("bkg range")); 
    double chi2 = plotPhotonsMassBkg[c]->chiSquare(3);
    Int_t ndof = nBinsMass-2;
    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;
    blind =false;
    if( blind ) {
     
      RooDataSet* data_down = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass < 178.");
      RooDataSet* data_up = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass >402");

      data_up->plotOn(plotPhotonsMassBkg[c]);    
      data_down->plotOn(plotPhotonsMassBkg[c]); 


   
    } else {
      data[c]->plotOn(plotPhotonsMassBkg[c]);    
      } 
       
  
    plotPhotonsMassBkg[c]->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");
    plotPhotonsMassBkg[c]->SetAxisRange(0.001,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    plotPhotonsMassBkg[c]->Draw();  

    TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Category %d",c), "brNDC");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(2),"Data","LPE");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(1),"Parametric Model: Expol","L");
    legdata->SetTextSize(0.035);
    legdata->SetTextFont(42);
    // legdata->SetTextAlign(31);
    legdata->SetBorderSize(0);
    legdata->SetFillStyle(0);
    legdata->Draw("same");

    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");

    //write down the chi2 of the fit on the
      
    TPaveText* label_chi2 = new TPaveText(0.5524194,0.6419492,0.796371,0.7690678, "brNDC");
    label_chi2->SetFillColor(kWhite);
    label_chi2->SetTextSize(0.035);
    label_chi2->SetTextFont(42);
    label_chi2->SetTextAlign(31); // align right
    label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", chi2));
    label_chi2->AddText(TString::Format("Chi square Prob = %.3f", prob));
    label_chi2->Draw("same");

    
    dobands = false;
    //********************************************************************************//
    if (dobands) {

      RooAbsPdf *cpdf; cpdf = PhotonsMassBkgTmp0;
      TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
      TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
      
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotPhotonsMassBkg[c]->getObject(1));
      
      for (int i=1; i<(plotPhotonsMassBkg[c]->GetXaxis()->GetNbins()+1); ++i) {
	double lowedge = plotPhotonsMassBkg[c]->GetXaxis()->GetBinLowEdge(i);
	double upedge  = plotPhotonsMassBkg[c]->GetXaxis()->GetBinUpEdge(i);
	double center  = plotPhotonsMassBkg[c]->GetXaxis()->GetBinCenter(i);
	
	double nombkg = nomcurve->interpolate(center);
	nlim->setVal(nombkg);
	PhotonsMass->setRange("errRange",lowedge,upedge);
	RooAbsPdf *epdf = 0;
	epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
	
	RooAbsReal *nll = epdf->createNLL(*(data[c]),Extended());
	RooMinimizer minim(*nll);
	minim.setStrategy(0);
	double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
	double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	
	minim.migrad();
	minim.minos(*nlim);
	// printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
	
	onesigma->SetPoint(i-1,center,nombkg);
	onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
	// eventually if cl = 0.95 this is the usual 1.92!      
	
	minim.migrad();
	minim.minos(*nlim);
	
	twosigma->SetPoint(i-1,center,nombkg);
	twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	delete nll;
	delete epdf;
      }

      PhotonsMass->setRange("errRange",minMassFit,maxMassFit);
      
      twosigma->SetLineColor(kGreen);
      twosigma->SetFillColor(kGreen);
      twosigma->SetMarkerColor(kGreen);
      twosigma->Draw("L3 SAME");
      
      onesigma->SetLineColor(kYellow);
      onesigma->SetFillColor(kYellow);
      onesigma->SetMarkerColor(kYellow);
      onesigma->Draw("L3 SAME");
      
      plotPhotonsMassBkg[c]->Draw("SAME"); 
    }
    int massI(mass);
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_M%d.root",c,massI));

    ctmp->SetLogy();
    plotPhotonsMassBkg[c]->SetAxisRange(1.3,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_LOG_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_LOG_M%d.root",c,massI));

  }



  return fitresult;
}

