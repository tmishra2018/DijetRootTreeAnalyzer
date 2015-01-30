void SigModelFitCBC(RooWorkspace* w, Float_t mass) {

  Int_t ncat = NCAT;
  RooDataSet* sigToFit[NCAT];
  
  
  Float_t MASS(mass);  
  Float_t minMassFit(mass*0.8);
  Float_t maxMassFit(mass*1.2); 
  
  if(mass==150.) minMassFit = MINmass;
  std::cout<<"----------------------------------------------------------------------------------------"<<std::endl;
  //PhotonsMass->setRange("sigrange",minMassFit-20,maxMassFit+20); 
  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
 
  // Fit Signal 
  for (int c=0; c<ncat; ++c) {
    cout << "---------- Category = " << c << endl;
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));
   
    RooRealVar* PhotonsMass = w->var("PhotonsMass"); 
    //cb
   
    RooFormulaVar CBpos_mean_draw(TString::Format("CBpos_mean_draw_cat%d",c),"","@0+250",*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)));
    RooFormulaVar CBpos_sigma(TString::Format("CBpos_sigma_cat%d",c),"","@0*250",*w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c) ));
    RooFormulaVar CBpos_alphaCB(TString::Format("CBpos_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)) );
    RooFormulaVar CBneg_alphaCB(TString::Format("CBneg_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)) );
    RooFormulaVar CBpos_n(TString::Format("CBpos_n_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_Npos_cat%d",c)) );
    RooFormulaVar CBneg_n(TString::Format("CBneg_n_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c)) );
    RooFormulaVar CBpos_frac(TString::Format("CBpos_frac_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)) );
    
    RooCBShape ResCBpos_draw(TString::Format("ResCBpos_draw_cat%d",c),TString::Format("ResCBpos_draw_cat%d",c) , *PhotonsMass, CBpos_mean_draw, CBpos_sigma, CBpos_alphaCB, CBpos_n) ;
    RooCBShape ResCBneg_draw(TString::Format("ResCBneg_draw_cat%d",c),TString::Format("ResCBneg_draw_cat%d",c) , *PhotonsMass, CBpos_mean_draw, CBpos_sigma, CBneg_alphaCB, CBneg_n) ;
   
    //PhotonsMass->setBins(40000, "cache");  
   

    //add CB neg + Gauss
    RooFormulaVar Gauss_frac(TString::Format("Gauss_frac_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_fracGauss_cat%d",c)));    
    RooFormulaVar Gauss_sigma(TString::Format("Gauss_sigma_cat%d",c),"","@0*250",*w->var(TString::Format("ReducedMass_sig_sigmaGauss_cat%d",c)));
    RooGaussian ResGauss_draw(TString::Format("ResGauss_draw_cat%d",c),TString::Format("ResGauss_draw_cat%d",c),*PhotonsMass, CBpos_mean_draw, Gauss_sigma );
    RooAddPdf ResAddGaussPdf_draw(TString::Format("PhotonsMassSig_cat%d",c),TString::Format("PhotonsMassSig_cat%d",c), RooArgList(ResCBneg_draw, ResGauss_draw), Gauss_frac);
    // w->import(ResAddGaussPdf_draw);
    
    //CBC
    RooFormulaVar CBC_mean(TString::Format("CBC_mean_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_sig_mean_cat%d",c)) );
    RooFormulaVar CBC_sigma(TString::Format("CBC_sigma_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_sig_sigma_cat%d",c)) );
    RooFormulaVar CBC_alphaC(TString::Format("CBC_alphaC_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_sig_alphaC_cat%d",c)) );
    RooFormulaVar CBC_alphaCB(TString::Format("CBC_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_sig_alphaCB_cat%d",c)) );
    RooFormulaVar CBC_n(TString::Format("CBC_n_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_sig_n_cat%d",c)) );

    
    RooCBCrujffPdf ResCBCPdf_draw(TString::Format("PhotonsMassSig_cat%d",c),TString::Format("PhotonsMassSig_cat%d",c) , *PhotonsMass, CBC_mean, CBC_sigma, CBC_alphaC, CBC_alphaCB, CBC_n) ; 
    w->import(ResCBCPdf_draw);


    double width=0.1;
    if(width < 2.){ //if i want to plot the fit
      //RooFitResult* fitresults_Gauss = (RooFitResult* ) ResAddGaussPdf_draw.fitTo(*sigToFit[c], RooFit::Save(kTRUE));
      RooFitResult* fitresults_Gauss = (RooFitResult* ) ResCBCPdf_draw.fitTo(*sigToFit[c],Range(minMassFit,maxMassFit),SumW2Error(kTRUE), RooFit::Save(kTRUE));
      std::cout<<TString::Format("******************************** Signal Fit results Gauss mass %f cat %d***********************************", mass, c)<<std::endl;
      fitresults_Gauss->Print("V");
      
      RooPlot* plotOnlyResPdf = w->var("PhotonsMass")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
      sigToFit[c]->plotOn(plotOnlyResPdf, LineColor(kRed), LineStyle(kDashed));
      double max = plotOnlyResPdf->GetMaximum();
      
      // Plot to verify everything is ok
      RooPlot* plotPhotonsMassAll = w->var("PhotonsMass")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
      sigToFit[c]->plotOn(plotPhotonsMassAll);

      // ResAddGaussPdf_draw.plotOn(plotPhotonsMassAll, LineColor(kRed), LineStyle(kDashed), NormRange("sigrange"));
      //ResAddGaussPdf_draw.plotOn(plotPhotonsMassAll,Components(TString::Format("ResGauss_draw_cat%d",c)), LineColor(kOrange), LineStyle(kDashed));
      //ResAddGaussPdf_draw.plotOn(plotPhotonsMassAll,Components(TString::Format("ResCBneg_draw_cat%d",c)), LineColor(kViolet), LineStyle(kDashed));
      ResCBCPdf_draw.plotOn(plotPhotonsMassAll, LineColor(kRed), LineStyle(kDashed));

      
      TCanvas* c1 = new TCanvas("c1","PhotonsMass",0,0,800,800);
      c1->cd(1);
      
      plotPhotonsMassAll->Draw();  
      
      TLegend *legmc = new TLegend(0.5491457,0.75,0.801457,0.9340659, TString::Format("Category %d",c), "brNDC");
      legmc->AddEntry(plotPhotonsMassAll->getObject(0),"Simulation","LPE");
      //legmc->AddEntry(plotPhotonsMassAll->getObject(2),"CB + Gauss","L");
      //legmc->AddEntry(plotPhotonsMassAll->getObject(3),"Gauss","L");
      //legmc->AddEntry(plotPhotonsMassAll->getObject(4),"CB","L");
      
      legmc->SetTextSize(0.0206044);
      legmc->SetTextFont(42);
      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      legmc->Draw();
      
      label_cms->Draw("same");
      label_sqrt->Draw("same");
      

      
      int massI(mass);
      c1->SetLogy();
      plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.0000000001,max*10. );
      c1->SaveAs("plots/prelimSignalCBGaussCconvBW"+TString::Format("_M%d_cat%d_LOG.png",massI,c));
      c1->SaveAs("plots/prelimSignalCBGaussCconvBW"+TString::Format("_M%d_cat%d_LOG.root",massI,c));
      
      /*  plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.0001, max*1.1 );
      c1->SetLogy(0);
      c1->SaveAs("plots/prelimSignalCBGaussCconvBW"+TString::Format("_M%d_cat%d.png",massI, c));
      c1->SaveAs("plots/prelimSignalCBGaussCconvBW"+TString::Format("_M%d_cat%d.root",massI, c));
      
        c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format("_M%d_cat%d.png",massI, c));
	   c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format("_M%d_cat%d.root",massI, c));
	   
	   c1->SetLogy();
	   c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format("_M%d_cat%d_LOG.png",massI,c));
	   c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format("_M%d_cat%d_LOG.root",massI,c));
      */
      
    }
    
    // IMPORTANT: fix all pdf parameters to constant
    
    /*  w->defineSet(TString::Format("ConvolutedPdfGaussParam_cat%d",c),RooArgSet( *w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)), 
									       *w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)),
									       *w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c)),	   
									       *w->var(TString::Format("ReducedMass_sig_fracGauss_cat%d",c)),  
									       *w->var(TString::Format("ReducedMass_sig_sigmaGauss_cat%d",c)),  
									       *w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),
									       *w->var(TString::Format("sigmaBW_var_cat%d",c))));
    
    SetConstantParams(w->set(TString::Format("ConvolutedPdfGaussParam_cat%d",c)));
    */

    w->defineSet(TString::Format("CBCParam_cat%d",c),RooArgSet(  *w->var(TString::Format("PhotonsMass_sig_alphaCB_cat%d",c)),
										*w->var(TString::Format("PhotonsMass_sig_n_cat%d",c)),	   
										*w->var(TString::Format("PhotonsMass_sig_alphaC_cat%d",c)),  
										*w->var(TString::Format("PhotonsMass_sig_mean_cat%d",c)),  
										*w->var(TString::Format("PhotonsMass_sig_sigma_cat%d",c))));  
										
    
    SetConstantParams(w->set(TString::Format("CBCParam_cat%d",c)));
    //w->Print("V");
    
  }

}



void SigModelResponseFcnFit(RooWorkspace* w, Float_t mass, std::string model) {


  
  //chain summing up all production modes
  TChain* sigTree1  = new TChain();
  //histograms_CMS-HGG_Moriond_M150.root
  
  sigTree1->Add("histograms_CMS-HGG_08052014_MC.root/ggh_m150_8TeV");
  /*  if(model!="GGH"){
  sigTree1->Add("histograms_CMS-HGG_Moriond_M150.root/vbf_m150_8TeV");
  if(model!="VBF"){
  sigTree1->Add("histograms_CMS-HGG_Moriond_M150.root/wzh_m150_8TeV");
  sigTree1->Add("histograms_CMS-HGG_Moriond_M150.root/tth_m150_8TeV");
  }
  }*/
  sigTree1->SetTitle("sigTree1");
  sigTree1->SetName("sigTree1");

  // Variables
  RooArgSet* ntplVars = defineVariablesM150();
  ntplVars->add(*w->var("PhotonsMassTrue"));
  (*w->var("PhotonsMassTrue")).setRange(130,230);
  ntplVars->add(*w->var("PhotonsMass"));
  (*w->var("PhotonsMass")).setRange(130,230);
  

  //TString mainCut1 = TString::Format("PhotonsMass > 130 && PhotonsMass<1000");   // livia
  TString mainCut1 = TString::Format("1");   // giulia
  RooDataSet sigWeighted("sigWeighted","dataset",sigTree1,*ntplVars,mainCut1,"evweight");

  RooRealVar* PhotonsMass = w->var("PhotonsMass");
   
  RooRealVar *mH = new RooRealVar("MH", "MH", MINmass, MAXmass);
  mH->setVal(mass);
  mH->setConstant();
  w->import(*mH);

   RooFormulaVar *massReduced_formula  =  new RooFormulaVar("massReduced_formula","","@0/@1 -1",RooArgList(*w->var("PhotonsMass"),*w->var("PhotonsMassTrue")));
  //RooFormulaVar *massReduced_formula  =  new RooFormulaVar("massReduced_formula","","@0/150 -1",RooArgList(*w->var("PhotonsMass")));
  // = new RooFormulaVar("massReduced_formula","","@0/250 -1",*w->var("PhotonsMass"));
  RooRealVar* massReduced = (RooRealVar*) sigWeighted.addColumn(*massReduced_formula);
  massReduced->SetName("massReduced");
  massReduced->SetTitle("massReduced");
  w->import(*massReduced);  
  massReduced->setRange(-0.5, 0.5);

 // common preselection cut
  TString mainCut = TString::Format("massReduced>-0.5 && massReduced <0.5");   // livia
  
  
  RooDataSet* signal[NCAT];
  RooCBShape* ResponseCBpos[NCAT];
  RooCBShape* ResponseCBneg[NCAT];
  RooGaussian* ResponseGauss[NCAT];
  RooAddPdf* ResponseAddGauss[NCAT];
  RooAddPdf* ResponseAdd[NCAT];
  
  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);

 
  Double_t scaleSyst;
  Double_t smearSyst;
  
  for(int c = 0; c<NCAT; c++){
    
    TLatex *lat  = new TLatex(0.6,0.9,TString::Format("Cat: %d", c));  
    lat->SetTextSize(0.038);
    lat->SetTextAlign(11);
    lat->SetTextFont(42); 
    lat->SetNDC();

    RooRealVar *mShift = new RooRealVar(TString::Format("mShift_cat%d", c), TString::Format("mShift_cat%d",c), -10., 10.);
    mShift->setVal(0.);
    //  mShift->setConstant();
    w->import(*mShift);

    RooRealVar *mSmear = new RooRealVar(TString::Format("mSmear_cat%d", c), TString::Format("mSmear_cat%d",c), -10., 10.);
    mSmear->setVal(0.);
    //  mSmear->setConstant();
    w->import(*mSmear);

    w->Print();

   // 1)  prime 4 cat livia
   if (c==0) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)<1.5 && abs(ph2_eta)<1.5) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
   if (c==1) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)<1.5 && abs(ph2_eta)<1.5) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
   if (c==2) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)>1.5 || abs(ph2_eta)>1.5) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
   if (c==3) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)>1.5 || abs(ph2_eta)>1.5) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));

   // w->import(*signal[c],Rename(TString::Format("SigWeight_cat%d",c)));
   //add cb neg +pos

   if(c==0 || c==1)scaleSyst = 0.005;
   if(c==2 || c==3)scaleSyst = 0.007;
   if(c==0)smearSyst = 0.005;
   if(c==1)smearSyst = 0.0058;
   if(c==2 || c==3)smearSyst = 0.01;

   //cb pos                                                                                                                     
   //  RooFormulaVar CBpos_mean(TString::Format("ReducedMass_CBpos_sig_mean_cat%d",c),"", TString::Format("@0+%f*@1", scaleSyst), RooArgList(*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)), *w->var(TString::Format("mShift_cat%d",c))));
   RooFormulaVar CBpos_mean(TString::Format("ReducedMass_CBpos_sig_mean_cat%d",c),"", "@0", *w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)));
   RooFormulaVar CBpos_sigma(TString::Format("ReducedMass_CBpos_sig_sigma_cat%d",c), "", "sqrt(@0*@0)", *w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)));
   RooFormulaVar CBpos_alphaCB(TString::Format("ReducedMass_CBpos_sig_alphaCB_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)));
   RooFormulaVar CBpos_n(TString::Format("ReducedMass_CBpos_sig_n_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_Npos_cat%d",c)));
     
   //cb neg
   RooFormulaVar CBneg_n(TString::Format("ReducedMass_CBneg_sig_n_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_Nneg_cat%d",c)));
   RooFormulaVar CBneg_alphaCB(TString::Format("ReducedMass_CBneg_sig_alphaCB_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)));
   
   ResponseCBpos[c] =  new RooCBShape(TString::Format("ResponseCBpos_cat%d",c),TString::Format("ResponseCBpos_cat%d",c) , *massReduced, CBpos_mean, CBpos_sigma,  CBpos_alphaCB, CBpos_n) ;
   

   ResponseCBneg[c] =  new RooCBShape(TString::Format("ResponseCBneg_cat%d",c),TString::Format("ResponseCBneg_cat%d",c) , *massReduced, CBpos_mean, CBpos_sigma,  CBneg_alphaCB, CBneg_n) ;
   

   
   RooFormulaVar CB_frac(TString::Format("ReducedMass_CBpos_sig_frac_cat%d",c),"", "@0", *w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)));
   w->import(CB_frac);  
   ResponseAdd[c]= new RooAddPdf(TString::Format("ResponseAddPdf_cat%d",c),TString::Format("ResponseAddPdf_cat%d",c) , RooArgList(*ResponseCBneg[c], *ResponseCBpos[c]), CB_frac);
   w->import(*ResponseAdd[c]);
   
  
  
  
   RooFitResult* fitresults = (RooFitResult* ) ResponseAdd[c]->fitTo(*signal[c],SumW2Error(kTRUE), Range(-1, 1), RooFit::Save(kTRUE));
   std::cout<<TString::Format("******************************** Signal Fit results CB+CB  mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresults->Print("V");
   

     RooPlot* plotG = massReduced->frame(Range(-0.12, 0.12),Title("Mass Reduced"), Bins(60));
    signal[c]->plotOn(plotG);
   

    ResponseAdd[c]->plotOn(plotG, LineColor(kBlue));
    ResponseAdd[c]->plotOn(plotG,Components(TString::Format("ResponseCBneg_cat%d",c)), LineColor(kGreen), LineStyle(kDashed));
    ResponseAdd[c]->plotOn(plotG,Components(TString::Format("ResponseCBpos_cat%d",c)), LineColor(kRed), LineStyle(kDashed));
  
    plotG->GetYaxis()->SetRangeUser(0.01,plotG->GetMaximum()*10 );
    plotG->GetXaxis()->SetTitle("#frac{#Delta m}{m}");
    plotG->GetXaxis()->SetTitleFont(42);
    plotG->GetXaxis()->SetTitleSize(0.05);
    TLegend* legmc = new TLegend(0.58, 0.52, 0.91, 0.89, "", "brNDC");
    legmc->SetTextSize(0.0286044);  
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    // legmc->SetHeader(("Model: "+model).c_str());
    legmc->AddEntry(plotG->getObject(0),"m_{#gamma#gamma} = 150 GeV","LPE");    
    

    legmc->AddEntry(plotG->getObject(1),"Sum of two CB ","L");
    legmc->AddEntry(plotG->getObject(2),"CB 1","L");   
    legmc->AddEntry(plotG->getObject(3),"CB 2","L");   
    plotG->Draw();
    
    lat->Draw("same");
    legmc->Draw("same");
    int iPos=11 ;
    CMS_lumi( c1,true,iPos );
    
    c1->SetLogy();
    
    c1->SaveAs(TString::Format(("plots/responseFcnFitCBCB_cat%d_LOG_"+model+".png").c_str(),c)); 
    c1->SaveAs(TString::Format(("plots/responseFcnFitCBCB_cat%d_LOG_"+model+".pdf").c_str(),c)); 
    // c1->SaveAs(TString::Format("plots/responseFcnFitCBCB_cat%d_LOG_%s.eps",c,model)); 
    
    
    plotG->GetYaxis()->SetRangeUser(0.0001,plotG->GetMaximum()*0.12 );
    
    c1->SetLogy(0);  
    
    //giulia
    // if(c==0){
    //   c1->SaveAs(TString::Format(("~/www/plotsPAS/responseFcnFitCBCB_cat%d_"+model+".png").c_str(),c,model)); 
    //   c1->SaveAs(TString::Format(("~/www/plotsPAS/responseFcnFitCBCB_cat%d_"+model+".pdf").c_str(),c,model)); 
    //   c1->SaveAs(TString::Format(("~/www/plotsPAS/responseFcnFitCBCB_cat%d_"+model+".C").c_str(),c,model)); 
    //   c1->SaveAs(TString::Format(("~/www/plotsPAS/responseFcnFitCBCB_cat%d_"+model+".root").c_str(),c,model));
    //   c1->SetLogy(1);  
    //   plotG->GetYaxis()->SetRangeUser(0.01,100 );
    //   c1->SaveAs(TString::Format(("~/www/plotsPAS/responseFcnFitCBCB_cat%d_"+model+"_LOG.png").c_str(),c,model)); 
    //   c1->SaveAs(TString::Format(("~/www/plotsPAS/responseFcnFitCBCB_cat%d_"+model+"_LOG.pdf").c_str(),c,model)); 
    //   c1->SaveAs(TString::Format(("~/www/plotsPAS/responseFcnFitCBCB_cat%d_"+model+"_LOG.C").c_str(),c,model)); 
    //   c1->SaveAs(TString::Format(("~/www/plotsPAS/responseFcnFitCBCB_cat%d_"+model+"_LOG.root").c_str(),c,model)); 
    // }
    
    w->defineSet(TString::Format("ResponseAddPdfParam_cat%d",c),RooArgSet(*w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)), 
									  *w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_Npos_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c)),	   
									  *w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)),  
									  *w->var(TString::Format("ReducedMass_sig_mean_cat%d",c))));
    SetConstantParams(w->set(TString::Format("ResponseAddPdfParam_cat%d",c)));
  }
  

 

}





// Fit signal with model gauss pdfs
void SigModelFitConvBW(RooWorkspace* w, Float_t mass, Double_t width, std::string model) {

  Int_t ncat = NCAT;
  RooDataSet* sigToFit[NCAT];
  
  
  Float_t MASS(mass);  
  Float_t minMassFit(mass*0.8);
  Float_t maxMassFit(mass*1.2); 
  RooRealVar* PhotonsMass = w->var("PhotonsMass"); 
  Double_t scaleSyst;
  Double_t smearSyst;

  std::cout<<"----------------------------------------------------------------------------------------"<<std::endl;
  //  PhotonsMass->setRange("sigrange",minMassFit-20,maxMassFit+20); 
  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  TFile* f = new TFile("sigShapeCorrections.root", "READ");
  // Fit Signal 
  for (int c=0; c<ncat; ++c) {
    cout << "---------- Category = " << c << endl;
    //  sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));
    //  w->import(*sigToFit[c]);

    //introduce systs 
    if(c==0 || c==1)scaleSyst = 0.005;
    if(c==2 || c==3)scaleSyst = 0.007;
    if(c==0)smearSyst = 0.005;
    if(c==1)smearSyst = 0.0058;
    if(c==2 || c==3)smearSyst = 0.01;

    //get sigma from TF1:   
    TF1* fcn = (TF1*)f->Get(TString::Format("f%d",c));
    Float_t massF = (Float_t) (*w->var("MH")).getVal();
    Float_t sigmaCorr = fcn->Eval(massF);
    if(massF==150)sigmaCorr=1;
    //std::cout<<"Mass: "<<massF<<" corr: "<<sigmaCorr<<std::endl;
    RooRealVar rooSigmaCorr (TString::Format("rooSigmaCorr_cat%d",c), TString::Format("rooSigmaCorr_cat%d",c), sigmaCorr, "");
    rooSigmaCorr.setConstant();
    w->import(rooSigmaCorr);


    ( *w->var(TString::Format("mShift_cat%d",c))).setConstant();
    ( *w->var(TString::Format("mSmear_cat%d",c))).setConstant();
    //cb
    RooFormulaVar CBpos_mean_draw(TString::Format("CBpos_mean_draw_cat%d",c),"","@0+@1",RooArgList(*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),*w->var("MH")));
    RooFormulaVar CBpos_mean(TString::Format("CBpos_mean_cat%d",c),"",TString::Format("@0+%f*@1", scaleSyst),RooArgList(*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),*w->var(TString::Format("mShift_cat%d",c))));
    //RooFormulaVar CBpos_mean(TString::Format("CBpos_mean_cat%d",c),"",TString::Format("@0"),*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)));   
    RooFormulaVar CBpos_sigma(TString::Format("CBpos_sigma_cat%d",c),"",TString::Format("(sqrt(@0*@0*@3*@3+%f*%f*@2)*@1)",smearSyst,smearSyst),RooArgList(*w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)),*w->var("MH"),*w->var(TString::Format("mSmear_cat%d",c)),*w->var(TString::Format("rooSigmaCorr_cat%d",c)) ) );

    //    std::cout<<"-------------------> SIGMA: "<<CBpos_sigma->getVal()<<"    MASS: "<<(*w->var("MH")).getVal()<<std::endl;


    RooFormulaVar CBpos_alphaCB(TString::Format("CBpos_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)) );
    RooFormulaVar CBneg_alphaCB(TString::Format("CBneg_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)) );
    RooFormulaVar CBpos_n(TString::Format("CBpos_n_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_Npos_cat%d",c)) );
    RooFormulaVar CBneg_n(TString::Format("CBneg_n_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c)) );
    RooFormulaVar CBpos_frac(TString::Format("CBpos_frac_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)) );
    
    RooCBShape ResCBpos_draw(TString::Format("ResCBpos_draw_cat%d",c),TString::Format("ResCBpos_draw_cat%d",c) , *PhotonsMass, CBpos_mean_draw, CBpos_sigma, CBpos_alphaCB, CBpos_n) ;
    RooCBShape ResCBneg_draw(TString::Format("ResCBneg_draw_cat%d",c),TString::Format("ResCBneg_draw_cat%d",c) , *PhotonsMass, CBpos_mean_draw, CBpos_sigma, CBneg_alphaCB, CBneg_n) ;
    RooCBShape ResCBpos(TString::Format("ResCBpos_cat%d",c),TString::Format("ResCBpos_cat%d",c) , *PhotonsMass, CBpos_mean, CBpos_sigma,CBpos_alphaCB, CBpos_n) ;
    RooCBShape ResCBneg(TString::Format("ResCBneg_cat%d",c),TString::Format("ResCBneg_cat%d",c) , *PhotonsMass, CBpos_mean, CBpos_sigma,CBneg_alphaCB, CBneg_n) ;
    PhotonsMass->setBins(40000, "cache");  
    //add CB pos + CB neg
    RooAddPdf ResAddPdf(TString::Format("ResAddPdf_cat%d",c),TString::Format("ResAddPdf_cat%d",c) , RooArgList(ResCBneg, ResCBpos), CBpos_frac);
    RooAddPdf ResAddPdf_draw(TString::Format("ResAddPdf_draw_cat%d",c),TString::Format("ResAddPdf_draw_cat%d",c) , RooArgList(ResCBneg_draw, ResCBpos_draw), CBpos_frac);


    
    //BW
    RooFormulaVar meanBW(TString::Format("meanBW_cat%d",c),"","@0",*w->var("MH"));  
    RooRealVar sigmaBW_var(TString::Format("sigmaBW_var_cat%d",c), TString::Format("sigmaBW_var_cat%d",c), width);
    std::cout<<" width:--------> "<<width<<std::endl;
    sigmaBW_var.setConstant();
    w->import(sigmaBW_var);
    
    
    RooFormulaVar* sigmaBW;
    if(width<1)sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0",*w->var(TString::Format("sigmaBW_var_cat%d",c))); 
    else if(width==2) sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.02",*w->var("MH"));   
    else if(width==5) sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.05",*w->var("MH"));   
    else if(width==7) sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.07",*w->var("MH"));   
    else if(width==10) sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.10",*w->var("MH"));   
    RooBreitWigner SigModelBW(TString::Format("BW_cat%d",c),TString::Format("BW_cat%d",c), *PhotonsMass, meanBW, *sigmaBW);
    /*
    // correction to BW
    TFile* file = new TFile("BW_corrections.root","READ");

    TGraph2D* g0 = (TGraph2D*)file->Get(TString::Format("p0_cat%d",c));   
    double var1 = g0->Interpolate(mass, width);  
    RooRealVar v1(TString::Format("v1_cat%d",c),TString::Format("v1_cat%d",c), var1);
    v1.setConstant();
    std::cout<<"++++++++++++++++++++++++++ "<<var1<<std::endl; 
    w->import(v1);
 
    RooFormulaVar f1(TString::Format("f1_cat%d",c), TString::Format("f1_cat%d",c), "@0",*w->var(TString::Format("v1_cat%d",c)));
   
   
    RooGenericPdf pf(TString::Format("pf_cat%d",c),  "exp(@1*(@0-@2))", RooArgList(*PhotonsMass, f1,*w->var("MH")));
    RooProdPdf* prod;*/
    RooFFTConvPdf*  ConvolutedRes_CB;
   
     
      ConvolutedRes_CB = new RooFFTConvPdf(TString::Format("PhotonsMassSig_cat%d",c),TString::Format("PhotonsMassSig_cat%d",c), *PhotonsMass,SigModelBW, ResAddPdf);
      //prod = new RooProdPdf(TString::Format("prod_cat%d",c),TString::Format("prod_cat%d",c), RooArgList(pf, *ConvolutedRes_CB)); 
      //prod->SetName(TString::Format("PhotonsMassSig_cat%d",c));
      //prod->SetTitle(TString::Format("PhotonsMassSig_cat%d",c));
      //w->import(*prod);
    // std::cout<<TString::Format("******************************** Signal Fit results CB mass %f cat %d***********************************", mass, c)<<std::endl;
    // fitresults_CB->Print("V");
    w->import(*ConvolutedRes_CB);
    // std::cout<<".............> "<<c<<std::endl;
    
    RooHistFunc* rooFunc_norm = getRooHistFunc(c,w->var("MH"), model );
    w->import(*rooFunc_norm);

    std::cout<<"SIG NORM ----->"<<rooFunc_norm->getVal(*w->var("MH"));
    // w->Print("V");

    if(width <2. && mass < 150){ //if i want to plot the fit
      sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));

      RooFitResult* fitresults_CB = (RooFitResult* ) ConvolutedRes_CB.fitTo(*sigToFit[c], RooFit::Save(kTRUE));
      fitresults_CB->Print("V");
      
      RooPlot* plotOnlyResPdf = w->var("PhotonsMass")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
      sigToFit[c]->plotOn(plotOnlyResPdf, LineColor(kRed), LineStyle(kDashed));
      double max = plotOnlyResPdf->GetMaximum();
      
      // Plot to verify everything is ok
      RooPlot* plotPhotonsMassAll = w->var("PhotonsMass")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
      sigToFit[c]->plotOn(plotPhotonsMassAll);
      SigModelBW.plotOn(plotPhotonsMassAll, LineColor(kGreen), LineStyle(kDashed));
      //  ResAddPdf_draw.plotOn(plotPhotonsMassAll, LineColor(kRed), LineStyle(kDashed));
      ConvolutedRes_CB.plotOn(plotPhotonsMassAll, LineColor(kBlue));

      
      TCanvas* c1 = new TCanvas("c1","PhotonsMass",0,0,800,800);
      c1->cd(1);
      
      plotPhotonsMassAll->Draw();  
      plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.01, max*1.2);
      plotPhotonsMassAll->GetXaxis()->SetRangeUser(210, 290);

      TLatex *lat  = new TLatex(0.55,0.9,TString::Format("Cat: %d", c));  
      lat->SetTextSize(0.038);
      lat->SetTextAlign(11);
      lat->SetTextFont(42); 
      lat->SetNDC();

      TLegend *legmc = new TLegend(0.55, 0.6, 0.87, 0.88, ("Model: "+model).c_str(), "brNDC");
      legmc->AddEntry(plotPhotonsMassAll->getObject(0),"Simulation","LPE");
      legmc->AddEntry(plotPhotonsMassAll->getObject(1),"BW","L");
      //legmc->AddEntry(plotPhotonsMassAll->getObject(2)," CB + CB ","L");
      legmc->AddEntry(plotPhotonsMassAll->getObject(2),"BW #otimes Resolution","L");
      legmc->SetTextSize(0.0206044);
      legmc->SetTextFont(42);
      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      legmc->Draw();
      lat->Draw("same");
      label_cms->Draw("same");
      label_sqrt->Draw("same");
      

      
      int massI(mass);
      c1->SetLogy();
     
      plotPhotonsMassAll->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");
     
      c1->SetLogy(0);
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_"+model+".png").c_str(),massI, c));
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_"+model+".pdf").c_str(),massI, c));
      
      c1->SetLogy();
      plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.01,max*10. );
      plotPhotonsMassAll->GetXaxis()->SetRangeUser(210, 290);
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+".png").c_str(),massI,c));
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+".pdf").c_str(),massI,c));
      //c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+".eps").c_str(),massI,c));


    }   
      //plot signal model at different widths
      bool plotW = true;
      if(plotW && c==0){
	RooRealVar var_01("var_w01", "var_w01", 0.1);
	var_01.setConstant();	
	RooFormulaVar sigmaBW_01("w01", "w01","@0", var_01);     
	RooBreitWigner SiBW_01("sigBW_01","sigBW_01" , *PhotonsMass, meanBW, sigmaBW_01);
	RooFFTConvPdf  ConvolutedRes_01("conv01", "conv01", *PhotonsMass,SiBW_01, ResAddPdf);

	RooRealVar var_3("var_w3", "var_w3",3);
	var_3.setConstant();	
	RooFormulaVar sigmaBW_3("w3", "w3","@0",  var_3);     
	RooBreitWigner SiBW_3("sigBW_3","sigBW_3" , *PhotonsMass, meanBW, sigmaBW_3);
	RooFFTConvPdf  ConvolutedRes_3("conv3", "conv3", *PhotonsMass,SiBW_3, ResAddPdf);

	RooRealVar var_6("var_w6", "var_w6", 6);
	var_6.setConstant();	
	RooFormulaVar sigmaBW_6("w6", "w6","@0", var_6);     
	RooBreitWigner SiBW_6("sigBW_6","sigBW_6" , *PhotonsMass, meanBW, sigmaBW_6);
	RooFFTConvPdf  ConvolutedRes_6("conv6", "conv6", *PhotonsMass,SiBW_6, ResAddPdf);

	RooRealVar var_10("var_w10", "var_w10", 10);
	var_10.setConstant();	
	RooFormulaVar sigmaBW_10("w10", "w10","@0", var_10);     
	RooBreitWigner SiBW_10("sigBW_10","sigBW_10" , *PhotonsMass, meanBW, sigmaBW_10);
	RooFFTConvPdf  ConvolutedRes_10("conv10", "conv10", *PhotonsMass,SiBW_10, ResAddPdf);

	RooRealVar var_15("var_w15", "var_w15",15);
	var_15.setConstant();	
	RooFormulaVar sigmaBW_15("w15", "w15","@0",  var_15);     
	RooBreitWigner SiBW_15("sigBW_15","sigBW_15" , *PhotonsMass, meanBW, sigmaBW_15);
	RooFFTConvPdf  ConvolutedRes_15("conv15", "conv15", *PhotonsMass,SiBW_15, ResAddPdf);

	RooPlot* plotWidths = w->var("PhotonsMass")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
	ConvolutedRes_15.plotOn( plotWidths, LineColor(kAzure+3));
	ConvolutedRes_10.plotOn( plotWidths, LineColor(kAzure+2));
	ConvolutedRes_6.plotOn( plotWidths, LineColor(kAzure+1));
	ConvolutedRes_3.plotOn( plotWidths, LineColor(kViolet+1));
	ConvolutedRes_01.plotOn( plotWidths, LineColor(kViolet-9));
	plotWidths->Draw();

	label_cms->Draw("same");
	label_sqrt->Draw("same");
      
	TLegend* leg = new TLegend(0.598851,0.6044755,0.84253,0.928252,"", "brNDC");
  
	leg->SetBorderSize(0.);
	leg->SetFillColor(kWhite);
	leg->SetTextFont(42);
	plotWidths->GetYaxis()->SetRangeUser(0.001, 1.);
	plotWidths->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");
	plotWidths->GetYaxis()->SetTitle(" ");
	leg->AddEntry(plotWidths->getObject(0), "Width = 15 GeV", "L");
	leg->AddEntry(plotWidths->getObject(1), "Width = 10 GeV", "L");
	leg->AddEntry(plotWidths->getObject(2),"Width = 6 GeV", "L");
	leg->AddEntry(plotWidths->getObject(3),"Width = 3 GeV", "L");
	leg->AddEntry(plotWidths->getObject(4), "Width = 0.1 GeV", "L");
	leg->Draw("same");

	c1->SaveAs("plots/SignalModels_differentWidths.png");
	c1->SaveAs("plots/SignalModels_differentWidths.pdf");
      
      


      

    }
    
    // IMPORTANT: fix all pdf parameters to constant
    
    w->defineSet(TString::Format("ConvolutedPdfParam_cat%d",c),RooArgSet( *w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)), 
									  *w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_Npos_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c)),	   
									  *w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)),  
									  *w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),
									  //	  *w->var(TString::Format("v1_cat%d",c)),
									  *w->var(TString::Format("sigmaBW_var_cat%d",c))));
        
    SetConstantParams(w->set(TString::Format("ConvolutedPdfParam_cat%d",c)));
    
    //w->Print("V");
    
  }

}



// Fit signal with model gauss pdfs
void SigModelFitConvRelBW(RooWorkspace* w, Float_t mass, Double_t width, std::string model) {

  Int_t ncat = NCAT;
  RooDataSet* sigToFit[NCAT];
  
  
  Float_t MASS(mass);  
  Float_t minMassFit(mass*0.8);
  Float_t maxMassFit(mass*1.2); 
  RooRealVar* PhotonsMass = w->var("PhotonsMass"); 
  Double_t scaleSyst;
  Double_t smearSyst;

  std::cout<<"----------------------------------------------------------------------------------------"<<std::endl;
  //  PhotonsMass->setRange("sigrange",minMassFit-20,maxMassFit+20); 
  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  TFile* f = new TFile("sigShapeCorrections.root", "READ");
  // Fit Signal 
  for (int c=0; c<ncat; ++c) {
    cout << "---------- Category = " << c << endl;
    //  sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));
    //  w->import(*sigToFit[c]);

    //introduce systs 
    if(c==0 || c==1)scaleSyst = 0.005;
    if(c==2 || c==3)scaleSyst = 0.007;
    if(c==0)smearSyst = 0.005;
    if(c==1)smearSyst = 0.0058;
    if(c==2 || c==3)smearSyst = 0.01;

    //get sigma from TF1:   
    TF1* fcn = (TF1*)f->Get(TString::Format("f%d",c));
    Float_t massF = (Float_t) (*w->var("MH")).getVal();
    Float_t sigmaCorr = fcn->Eval(massF);
    if(massF==150)sigmaCorr=1;
    //std::cout<<"Mass: "<<massF<<" corr: "<<sigmaCorr<<std::endl;
    RooRealVar rooSigmaCorr (TString::Format("rooSigmaCorr_cat%d",c), TString::Format("rooSigmaCorr_cat%d",c), sigmaCorr, "");
    rooSigmaCorr.setConstant();
    w->import(rooSigmaCorr);


    ( *w->var(TString::Format("mShift_cat%d",c))).setConstant();
    ( *w->var(TString::Format("mSmear_cat%d",c))).setConstant();
    //cb
    RooFormulaVar CBpos_mean_draw(TString::Format("CBpos_mean_draw_cat%d",c),"","@0+@1",RooArgList(*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),*w->var("MH")));
    RooFormulaVar CBpos_mean(TString::Format("CBpos_mean_cat%d",c),"",TString::Format("@0+%f*@1", scaleSyst),RooArgList(*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),*w->var(TString::Format("mShift_cat%d",c))));
    //RooFormulaVar CBpos_mean(TString::Format("CBpos_mean_cat%d",c),"",TString::Format("@0"),*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)));   
    RooFormulaVar CBpos_sigma(TString::Format("CBpos_sigma_cat%d",c),"",TString::Format("(sqrt(@0*@0*@3*@3+%f*%f*@2)*@1)",smearSyst,smearSyst),RooArgList(*w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)),*w->var("MH"),*w->var(TString::Format("mSmear_cat%d",c)),*w->var(TString::Format("rooSigmaCorr_cat%d",c)) ) );

    //    std::cout<<"-------------------> SIGMA: "<<CBpos_sigma->getVal()<<"    MASS: "<<(*w->var("MH")).getVal()<<std::endl;


    RooFormulaVar CBpos_alphaCB(TString::Format("CBpos_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)) );
    RooFormulaVar CBneg_alphaCB(TString::Format("CBneg_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)) );
    RooFormulaVar CBpos_n(TString::Format("CBpos_n_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_Npos_cat%d",c)) );
    RooFormulaVar CBneg_n(TString::Format("CBneg_n_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c)) );
    RooFormulaVar CBpos_frac(TString::Format("CBpos_frac_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)) );
    
    RooCBShape ResCBpos_draw(TString::Format("ResCBpos_draw_cat%d",c),TString::Format("ResCBpos_draw_cat%d",c) , *PhotonsMass, CBpos_mean_draw, CBpos_sigma, CBpos_alphaCB, CBpos_n) ;
    RooCBShape ResCBneg_draw(TString::Format("ResCBneg_draw_cat%d",c),TString::Format("ResCBneg_draw_cat%d",c) , *PhotonsMass, CBpos_mean_draw, CBpos_sigma, CBneg_alphaCB, CBneg_n) ;
    RooCBShape ResCBpos(TString::Format("ResCBpos_cat%d",c),TString::Format("ResCBpos_cat%d",c) , *PhotonsMass, CBpos_mean, CBpos_sigma,CBpos_alphaCB, CBpos_n) ;
    RooCBShape ResCBneg(TString::Format("ResCBneg_cat%d",c),TString::Format("ResCBneg_cat%d",c) , *PhotonsMass, CBpos_mean, CBpos_sigma,CBneg_alphaCB, CBneg_n) ;
    PhotonsMass->setBins(40000, "cache");  
    //add CB pos + CB neg
    RooAddPdf ResAddPdf(TString::Format("ResAddPdf_cat%d",c),TString::Format("ResAddPdf_cat%d",c) , RooArgList(ResCBneg, ResCBpos), CBpos_frac);
    RooAddPdf ResAddPdf_draw(TString::Format("ResAddPdf_draw_cat%d",c),TString::Format("ResAddPdf_draw_cat%d",c) , RooArgList(ResCBneg_draw, ResCBpos_draw), CBpos_frac);


    
    //BW
    RooFormulaVar meanBW(TString::Format("meanBW_cat%d",c),"","@0",*w->var("MH"));  
    RooRealVar sigmaBW_var(TString::Format("sigmaBW_var_cat%d",c), TString::Format("sigmaBW_var_cat%d",c), width);
    std::cout<<" width:--------> "<<width<<std::endl;
    sigmaBW_var.setConstant();
    w->import(sigmaBW_var);
    
    
    RooFormulaVar* sigmaBW;
    if(width<1)sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0",*w->var(TString::Format("sigmaBW_var_cat%d",c))); 
    else if(width==2) sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.02",*w->var("MH"));   
    else if(width==5) sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.05",*w->var("MH"));   
    else if(width==7) sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.07",*w->var("MH"));   
    else if(width==10) sigmaBW = new RooFormulaVar(TString::Format("sigmaBW_cat%d",c),"","@0*0.10",*w->var("MH"));   


    RooGenericPdf SigModelBW(TString::Format("BW_cat%d",c),"1 / ( TMath::Power( TMath::Power(@0,2) - TMath::Power(@1,2) , 2 ) + TMath::Power(@1,2)*TMath::Power(@2,2) )", RooArgSet(*PhotonsMass, meanBW, *sigmaBW));
    
    RooFFTConvPdf*  ConvolutedRes_CB;
   
     
      ConvolutedRes_CB = new RooFFTConvPdf(TString::Format("PhotonsMassSig_cat%d",c),TString::Format("PhotonsMassSig_cat%d",c), *PhotonsMass,SigModelBW, ResAddPdf);
    w->import(*ConvolutedRes_CB);
    // std::cout<<".............> "<<c<<std::endl;
    
    //RooHistFunc* rooFunc_norm = getRooHistFunc(c,w->var("MH"), model );
    RooHistFunc* rooFunc_norm = getNorm2D(c,w->var("MH"), w->var("MH")->getVal(), width,model);
    w->import(*rooFunc_norm);
    std::cout<<"SIG NORM ----->"<<rooFunc_norm->getVal(*w->var("MH"))<<std::endl;

    /*   RooHistFunc* rooFunc_norm2D_0 = getNorm2D(c,w->var("MH"),w->var("MH")->getVal(), 0);
    RooHistFunc* rooFunc_norm2D_2 = getNorm2D(c,w->var("MH"),w->var("MH")->getVal(),2);
    RooHistFunc* rooFunc_norm2D_5 = getNorm2D(c,w->var("MH"),w->var("MH")->getVal(),5);
    RooHistFunc* rooFunc_norm2D_10 = getNorm2D(c,w->var("MH"),w->var("MH")->getVal(), 10);

    std::cout<<"SIG NORM ----->"<<rooFunc_norm2D_0->getVal(*w->var("MH"))<<"      "<<rooFunc_norm2D_2->getVal(*w->var("MH"))<<"  "<<rooFunc_norm2D_5->getVal(*w->var("MH"))<<"      "<<rooFunc_norm2D_10->getVal(*w->var("MH"))<<std::endl;;*/
    // w->Print("V");

    if(width <2. && mass <= 150){ //if i want to plot the fit
      sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));

      RooFitResult* fitresults_CB = (RooFitResult* ) ConvolutedRes_CB.fitTo(*sigToFit[c], RooFit::Save(kTRUE));
      fitresults_CB->Print("V");
      
      RooPlot* plotOnlyResPdf = w->var("PhotonsMass")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
      sigToFit[c]->plotOn(plotOnlyResPdf, LineColor(kRed), LineStyle(kDashed));
      double max = plotOnlyResPdf->GetMaximum();
      
      // Plot to verify everything is ok
      RooPlot* plotPhotonsMassAll = w->var("PhotonsMass")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
      sigToFit[c]->plotOn(plotPhotonsMassAll);
      SigModelBW.plotOn(plotPhotonsMassAll, LineColor(kGreen), LineStyle(kDashed));
      //  ResAddPdf_draw.plotOn(plotPhotonsMassAll, LineColor(kRed), LineStyle(kDashed));
      ConvolutedRes_CB.plotOn(plotPhotonsMassAll, LineColor(kBlue));

      
      TCanvas* c1 = new TCanvas("c1","PhotonsMass",0,0,800,800);
      c1->cd(1);
      
      plotPhotonsMassAll->Draw();  
      plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.01, max*1.2);
      plotPhotonsMassAll->GetXaxis()->SetRangeUser(210, 290);

      TLatex *lat  = new TLatex(0.55,0.9,TString::Format("Cat: %d", c));  
      lat->SetTextSize(0.038);
      lat->SetTextAlign(11);
      lat->SetTextFont(42); 
      lat->SetNDC();

      TLegend *legmc = new TLegend(0.55, 0.6, 0.87, 0.88);
      legmc->AddEntry(plotPhotonsMassAll->getObject(0),"Simulation","LPE");
      legmc->AddEntry(plotPhotonsMassAll->getObject(1),"BW","L");
      //legmc->AddEntry(plotPhotonsMassAll->getObject(2)," CB + CB ","L");
      legmc->AddEntry(plotPhotonsMassAll->getObject(2),"BW #otimes Resolution","L");
      legmc->SetTextSize(0.0286044);
      legmc->SetTextFont(42);
      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      legmc->Draw();
      lat->Draw("same");
      //label_cms->Draw("same");
      //label_sqrt->Draw("same");
      int iPos=11 ;
      CMS_lumi( c1,true,iPos );

      
      int massI(mass);
      c1->SetLogy();
     
      plotPhotonsMassAll->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");
     
      c1->SetLogy(0);
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_"+model+".png").c_str(),massI, c));
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_"+model+".pdf").c_str(),massI, c));
      
      // if(c==0){
      // 	c1->SaveAs("~/www/plotsPAS/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_"+model+".png").c_str(),massI, c));
      // 	c1->SaveAs("~/www/plotsPAS/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_"+model+".pdf").c_str(),massI, c));
      // 	c1->SaveAs("~/www/plotsPAS/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_"+model+".C").c_str(),massI, c));
      // 	c1->SaveAs("~/www/plotsPAS/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_"+model+".root").c_str(),massI, c));
      // }

      c1->SetLogy();
      plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.01,max*10. );
      plotPhotonsMassAll->GetXaxis()->SetRangeUser(210, 290);
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+".png").c_str(),massI,c));
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+".pdf").c_str(),massI,c));

      // if(c==0){
      // 	c1->SaveAs("~/www/plotsPAS/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+".png").c_str(),massI,c));
      // 	c1->SaveAs("~/www/plotsPAS/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+".pdf").c_str(),massI,c));
      // 	c1->SaveAs("~/www/plotsPAS/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+".C").c_str(),massI,c));
      // 	c1->SaveAs("~/www/plotsPAS/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+".root").c_str(),massI,c));
      // }

    }
//plot signal model at different widths
      bool plotW = true;
      if(plotW && c==0){
	RooRealVar var_01("var_w01", "var_w01", 0.1);
	var_01.setConstant();	
	RooFormulaVar sigmaBW_01("w01", "w01","@0", var_01); 
	//	RooBreitWigner SiBW_01("sigBW_01","sigBW_01" , *PhotonsMass, meanBW, sigmaBW_01);
	RooGenericPdf SiBW_01("sigBW_01","1 / ( TMath::Power( TMath::Power(@0,2) - TMath::Power(@1,2) , 2 ) + TMath::Power(@1,2)*TMath::Power(@2,2) )", RooArgSet(*PhotonsMass, meanBW, sigmaBW_01));    
	RooFFTConvPdf  ConvolutedRes_01("conv01", "conv01", *PhotonsMass,SiBW_01, ResAddPdf);

	RooRealVar var_3("var_w3", "var_w3",3);
	var_3.setConstant();	
	RooFormulaVar sigmaBW_3("w3", "w3","@0",  var_3);     
	//	RooBreitWigner SiBW_3("sigBW_3","sigBW_3" , *PhotonsMass, meanBW, sigmaBW_3);
	RooGenericPdf SiBW_3("sigBW_3","1 / ( TMath::Power( TMath::Power(@0,2) - TMath::Power(@1,2) , 2 ) + TMath::Power(@1,2)*TMath::Power(@2,2) )", RooArgSet(*PhotonsMass, meanBW, sigmaBW_3));    
	RooFFTConvPdf  ConvolutedRes_3("conv3", "conv3", *PhotonsMass,SiBW_3, ResAddPdf);

	RooRealVar var_6("var_w6", "var_w6", 6);
	var_6.setConstant();	
	RooFormulaVar sigmaBW_6("w6", "w6","@0", var_6);     
	//RooBreitWigner SiBW_6("sigBW_6","sigBW_6" , *PhotonsMass, meanBW, sigmaBW_6);
	RooGenericPdf SiBW_6("sigBW_6","1 / ( TMath::Power( TMath::Power(@0,2) - TMath::Power(@1,2) , 2 ) + TMath::Power(@1,2)*TMath::Power(@2,2) )", RooArgSet(*PhotonsMass, meanBW, sigmaBW_3));    
	RooFFTConvPdf  ConvolutedRes_6("conv6", "conv6", *PhotonsMass,SiBW_6, ResAddPdf);

	RooRealVar var_10("var_w10", "var_w10", 10);
	var_10.setConstant();	
	RooFormulaVar sigmaBW_10("w10", "w10","@0", var_10);     
	//	RooBreitWigner SiBW_10("sigBW_10","sigBW_10" , *PhotonsMass, meanBW, sigmaBW_10);
	RooGenericPdf SiBW_10("sigBW_10","1 / ( TMath::Power( TMath::Power(@0,2) - TMath::Power(@1,2) , 2 ) + TMath::Power(@1,2)*TMath::Power(@2,2) )", RooArgSet(*PhotonsMass, meanBW, sigmaBW_10));    
	RooFFTConvPdf  ConvolutedRes_10("conv10", "conv10", *PhotonsMass,SiBW_10, ResAddPdf);

	RooRealVar var_15("var_w15", "var_w15",15);
	var_15.setConstant();	
	RooFormulaVar sigmaBW_15("w15", "w15","@0",  var_15);     
	//	RooBreitWigner SiBW_15("sigBW_15","sigBW_15" , *PhotonsMass, meanBW, sigmaBW_15);
	RooGenericPdf SiBW_15("sigBW_15","1 / ( TMath::Power( TMath::Power(@0,2) - TMath::Power(@1,2) , 2 ) + TMath::Power(@1,2)*TMath::Power(@2,2) )", RooArgSet(*PhotonsMass, meanBW, sigmaBW_15));    
	RooFFTConvPdf  ConvolutedRes_15("conv15", "conv15", *PhotonsMass,SiBW_15, ResAddPdf);

	RooPlot* plotWidths = w->var("PhotonsMass")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
	ConvolutedRes_15.plotOn( plotWidths, LineColor(kAzure+3));
	ConvolutedRes_10.plotOn( plotWidths, LineColor(kAzure+2));
	ConvolutedRes_6.plotOn( plotWidths, LineColor(kAzure+1));
	ConvolutedRes_3.plotOn( plotWidths, LineColor(kViolet+1));
	ConvolutedRes_01.plotOn( plotWidths, LineColor(kViolet-9));
	plotWidths->Draw();

	label_cms->Draw("same");
	label_sqrt->Draw("same");
      
	TLegend* leg = new TLegend(0.598851,0.6044755,0.84253,0.928252,"", "brNDC");
	std::cout<<meanBW->getVal()<<"   --------"<<std::endl;
	leg->SetBorderSize(0.);
	leg->SetFillColor(kWhite);
	leg->SetTextFont(42);
	plotWidths->GetYaxis()->SetRangeUser(0.001, 1.);
	plotWidths->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");
	plotWidths->GetYaxis()->SetTitle(" ");
	leg->AddEntry(plotWidths->getObject(0), "Width = 15 GeV", "L");
	leg->AddEntry(plotWidths->getObject(1), "Width = 10 GeV", "L");
	leg->AddEntry(plotWidths->getObject(2),"Width = 6 GeV", "L");
	leg->AddEntry(plotWidths->getObject(3),"Width = 3 GeV", "L");
	leg->AddEntry(plotWidths->getObject(4), "Width = 0.1 GeV", "L");
	leg->Draw("same");

	c1->SaveAs("plots/SignalModels_differentWidths.png");
	c1->SaveAs("plots/SignalModels_differentWidths.pdf");

	// c1->SaveAs("~/www/plotsNota/SignalModels_differentWidths.pdf");
	// c1->SaveAs("~/www/plotsNota/SignalModels_differentWidths.png");
      
      }

    
    
    // IMPORTANT: fix all pdf parameters to constant
    
    w->defineSet(TString::Format("ConvolutedPdfParam_cat%d",c),RooArgSet( *w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)), 
									  *w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_Npos_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c)),	   
									  *w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)),  
									  *w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),
									  //	  *w->var(TString::Format("v1_cat%d",c)),
									  *w->var(TString::Format("sigmaBW_var_cat%d",c))));
        
    SetConstantParams(w->set(TString::Format("ConvolutedPdfParam_cat%d",c)));
    
    //w->Print("V");
    
  }

}



// Fit signal with model gauss pdfs
void SigModelFitConvCPS(RooWorkspace* w, Float_t mass, Double_t width, std::string model) {

  Int_t ncat = NCAT;
  RooDataSet* sigToFit[NCAT];
  
  
  Float_t MASS(mass);  
  Float_t minMassFit(mass*0.8);
  Float_t maxMassFit(mass*1.2); 
  RooRealVar* PhotonsMass = w->var("PhotonsMass"); 
  Double_t scaleSyst;
  Double_t smearSyst;

  std::cout<<"----------------------------------------------------------------------------------------"<<std::endl;
  //  PhotonsMass->setRange("sigrange",minMassFit-20,maxMassFit+20); 
  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  TFile* f = new TFile("sigShapeCorrections.root", "READ");
  // Fit Signal 
  for (int c=0; c<ncat; ++c) {
    cout << "---------- Category = " << c << endl;
   
    //introduce systs 
    if(c==0 || c==1)scaleSyst = 0.005;
    if(c==2 || c==3)scaleSyst = 0.007;
    if(c==0)smearSyst = 0.005;
    if(c==1)smearSyst = 0.0058;
    if(c==2 || c==3)smearSyst = 0.01;

    //get sigma from TF1:   
    TF1* fcn = (TF1*)f->Get(TString::Format("f%d",c));
    Float_t massF = (Float_t) (*w->var("MH")).getVal();
    Float_t sigmaCorr = fcn->Eval(massF);
    if(massF==150)sigmaCorr=1;
    //std::cout<<"Mass: "<<massF<<" corr: "<<sigmaCorr<<std::endl;
    RooRealVar rooSigmaCorr (TString::Format("rooSigmaCorr_cat%d",c), TString::Format("rooSigmaCorr_cat%d",c), sigmaCorr, "");
    rooSigmaCorr.setConstant();
    w->import(rooSigmaCorr);


    ( *w->var(TString::Format("mShift_cat%d",c))).setConstant();
    ( *w->var(TString::Format("mSmear_cat%d",c))).setConstant();
    //cb
    RooFormulaVar CBpos_mean_draw(TString::Format("CBpos_mean_draw_cat%d",c),"","@0+@1",RooArgList(*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),*w->var("MH")));
    RooFormulaVar CBpos_mean(TString::Format("CBpos_mean_cat%d",c),"",TString::Format("@0+%f*@1", scaleSyst),RooArgList(*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),*w->var(TString::Format("mShift_cat%d",c))));
    //RooFormulaVar CBpos_mean(TString::Format("CBpos_mean_cat%d",c),"",TString::Format("@0"),*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)));   
    RooFormulaVar CBpos_sigma(TString::Format("CBpos_sigma_cat%d",c),"",TString::Format("(sqrt(@0*@0*@3*@3+%f*%f*@2)*@1)",smearSyst,smearSyst),RooArgList(*w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)),*w->var("MH"),*w->var(TString::Format("mSmear_cat%d",c)),*w->var(TString::Format("rooSigmaCorr_cat%d",c)) ) );

    //    std::cout<<"-------------------> SIGMA: "<<CBpos_sigma->getVal()<<"    MASS: "<<(*w->var("MH")).getVal()<<std::endl;


    RooFormulaVar CBpos_alphaCB(TString::Format("CBpos_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)) );
    RooFormulaVar CBneg_alphaCB(TString::Format("CBneg_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)) );
    RooFormulaVar CBpos_n(TString::Format("CBpos_n_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_Npos_cat%d",c)) );
    RooFormulaVar CBneg_n(TString::Format("CBneg_n_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c)) );
    RooFormulaVar CBpos_frac(TString::Format("CBpos_frac_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)) );
    
    RooCBShape ResCBpos_draw(TString::Format("ResCBpos_draw_cat%d",c),TString::Format("ResCBpos_draw_cat%d",c) , *PhotonsMass, CBpos_mean_draw, CBpos_sigma, CBpos_alphaCB, CBpos_n) ;
    RooCBShape ResCBneg_draw(TString::Format("ResCBneg_draw_cat%d",c),TString::Format("ResCBneg_draw_cat%d",c) , *PhotonsMass, CBpos_mean_draw, CBpos_sigma, CBneg_alphaCB, CBneg_n) ;
    RooCBShape ResCBpos(TString::Format("ResCBpos_cat%d",c),TString::Format("ResCBpos_cat%d",c) , *PhotonsMass, CBpos_mean, CBpos_sigma,CBpos_alphaCB, CBpos_n) ;
    RooCBShape ResCBneg(TString::Format("ResCBneg_cat%d",c),TString::Format("ResCBneg_cat%d",c) , *PhotonsMass, CBpos_mean, CBpos_sigma,CBneg_alphaCB, CBneg_n) ;
    PhotonsMass->setBins(40000, "cache");  
    //add CB pos + CB neg
    RooAddPdf ResAddPdf(TString::Format("ResAddPdf_cat%d",c),TString::Format("ResAddPdf_cat%d",c) , RooArgList(ResCBneg, ResCBpos), CBpos_frac);
    RooAddPdf ResAddPdf_draw(TString::Format("ResAddPdf_draw_cat%d",c),TString::Format("ResAddPdf_draw_cat%d",c) , RooArgList(ResCBneg_draw, ResCBpos_draw), CBpos_frac);


    
    //BW
    RooFormulaVar meanCPS(TString::Format("meanCPS_cat%d",c),"","@0",*w->var("MH"));  
    double cprime;
    if(w->var("MH")->getVal()==400. && width == 10 ) cprime=1.33;
    if(w->var("MH")->getVal()==400. && width == 5 ) cprime=0.66;
    if(w->var("MH")->getVal()==400. && width == 2 ) cprime=0.26;
    if(w->var("MH")->getVal()==800. && width == 10 ) cprime=0.26;
    if(w->var("MH")->getVal()==800. && width == 5 ) cprime=0.13;
    if(w->var("MH")->getVal()==800. && width == 2 ) cprime=0.053;
    if(w->var("MH")->getVal()==600. && width == 10 ) cprime=0.5;
    if(w->var("MH")->getVal()==600. && width == 5 ) cprime=0.25;

    RooRealVar cp_var(TString::Format("cp_var_cat%d",c), TString::Format("cp_var_cat%d",c), cprime);
    std::cout<<" width:--------> "<<width<<"    cprime: "<<cp_var.getVal()<<std::endl;
    cp_var.setConstant();
    w->import(cp_var);    
    RooFormulaVar cp(TString::Format("cp_cat%d",c),"","@0",*w->var(TString::Format("cp_var_cat%d",c)));   
    RooRealVar br_var(TString::Format("br_var_cat%d",c), TString::Format("br_var_cat%d",c),0.);
    br_var.setVal(0.);
    br_var.setConstant();
    w->import(br_var);   
    RooFormulaVar br(TString::Format("br_cat%d",c),"","@0",br_var);
    Bool_t is8 =  true;
   
    RooCPSHighMassGGHNoInterf cpsPdf(TString::Format("CPSpdf_cat%d",c),TString::Format("CPSpdf_cat%d",c),*PhotonsMass, meanCPS, cp,br,is8);
   
    RooFFTConvPdf*  ConvolutedRes_CB;
   
    
    ConvolutedRes_CB = new RooFFTConvPdf(TString::Format("PhotonsMassSig_cat%d",c),TString::Format("PhotonsMassSig_cat%d",c), *PhotonsMass,cpsPdf, ResAddPdf);    
    w->import(*ConvolutedRes_CB);
    
    
    RooHistFunc* rooFunc_norm = getRooHistFunc(c,w->var("MH"), model );
    w->import(*rooFunc_norm);

    std::cout<<"SIG NORM ----->"<<rooFunc_norm->getVal(*w->var("MH"))<<"      "<<getNorm2D(c,w->var("MH")->getVal(),width, model));
    // w->Print("V");

    if(width <2. && mass < 150){ //if i want to plot the fit
      sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));

      RooFitResult* fitresults_CB = (RooFitResult* ) ConvolutedRes_CB.fitTo(*sigToFit[c], RooFit::Save(kTRUE));
      fitresults_CB->Print("V");
      
      RooPlot* plotOnlyResPdf = w->var("PhotonsMass")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
      sigToFit[c]->plotOn(plotOnlyResPdf, LineColor(kRed), LineStyle(kDashed));
      double max = plotOnlyResPdf->GetMaximum();
      
      // Plot to verify everything is ok
      RooPlot* plotPhotonsMassAll = w->var("PhotonsMass")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
      sigToFit[c]->plotOn(plotPhotonsMassAll);
      SigModelBW.plotOn(plotPhotonsMassAll, LineColor(kGreen), LineStyle(kDashed));
      //  ResAddPdf_draw.plotOn(plotPhotonsMassAll, LineColor(kRed), LineStyle(kDashed));
      ConvolutedRes_CB.plotOn(plotPhotonsMassAll, LineColor(kBlue));

      
      TCanvas* c1 = new TCanvas("c1","PhotonsMass",0,0,800,800);
      c1->cd(1);
      
      plotPhotonsMassAll->Draw();  
      plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.01, max*1.2);
      plotPhotonsMassAll->GetXaxis()->SetRangeUser(210, 290);

      TLatex *lat  = new TLatex(0.55,0.9,TString::Format("Cat: %d", c));  
      lat->SetTextSize(0.038);
      lat->SetTextAlign(11);
      lat->SetTextFont(42); 
      lat->SetNDC();

      TLegend *legmc = new TLegend(0.55, 0.6, 0.87, 0.88, ("Model: "+model).c_str(), "brNDC");
      legmc->AddEntry(plotPhotonsMassAll->getObject(0),"Simulation","LPE");
      legmc->AddEntry(plotPhotonsMassAll->getObject(1),"BW","L");
      //legmc->AddEntry(plotPhotonsMassAll->getObject(2)," CB + CB ","L");
      legmc->AddEntry(plotPhotonsMassAll->getObject(2),"BW #otimes Resolution","L");
      legmc->SetTextSize(0.0206044);
      legmc->SetTextFont(42);
      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      legmc->Draw();
      lat->Draw("same");
      label_cms->Draw("same");
      label_sqrt->Draw("same");
      

      
      int massI(mass);
      c1->SetLogy();
     
      plotPhotonsMassAll->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");
     
      c1->SetLogy(0);
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_"+model+".png").c_str(),massI, c));
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_"+model+".pdf").c_str(),massI, c));
      
      c1->SetLogy();
      plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.01,max*10. );
      plotPhotonsMassAll->GetXaxis()->SetRangeUser(210, 290);
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+".png").c_str(),massI,c));
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+".pdf").c_str(),massI,c));
      //c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format(("_M%d_cat%d_LOG_"+model+".eps").c_str(),massI,c));


    }
    
    // IMPORTANT: fix all pdf parameters to constant
    
    w->defineSet(TString::Format("ConvolutedPdfParam_cat%d",c),RooArgSet( *w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)), 
									  *w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_Npos_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c)),	   
									  *w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)),  
									  *w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),
									  *w->var(TString::Format("br_var_cat%d",c)),
									  *w->var(TString::Format("cp_var_cat%d",c))));
        
    SetConstantParams(w->set(TString::Format("ConvolutedPdfParam_cat%d",c)));
    
    //w->Print("V");
    
  }

}


