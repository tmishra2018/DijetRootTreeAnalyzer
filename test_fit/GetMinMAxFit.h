RooHistFunc* getRooHistFuncFitMIN(int cat, RooRealVar* var){

  double mass[1] = {2000.};
  double c0[1] = {2000.};
  double c1[1] = {2000.};
  double c2[1] = {2000.};
  double c3[1] = {2000};

  
  TH1F* h_all = new TH1F("h_all", "h_all", 1, 0, 10000);
  for(int i=0;i<1;i++){
    std::cout<<cat<<std::endl;
    std::cout<<mass[i]<<"  "<<h_all->FindBin(mass[i])<<"   "<<c0[i]<<std::endl;
    if(cat==0) h_all->SetBinContent(h_all->FindBin(mass[i]),c0[i]);
    if(cat==1) h_all->SetBinContent(h_all->FindBin(mass[i]),c1[i]);
    if(cat==2) h_all->SetBinContent(h_all->FindBin(mass[i]),c2[i]);
    if(cat==3) h_all->SetBinContent(h_all->FindBin(mass[i]),c3[i]);
      
  }
  /* TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();
  //TPaveText* label_cms = get_labelCMS(0,"2014", false);
  //TPaveText* label_sqrt = get_labelSqrt(0);  

  h_all->GetXaxis()->SetTitle("m_{X} [GeV]");
  h_all->GetYaxis()->SetTitle("Fit Range: Minimum");
  h_all->Draw("");*/
 
  RooDataHist* rooData_all = new RooDataHist("rooData_all", "roData_all",*var,h_all);
  RooHistFunc* rooFunc_all = new RooHistFunc(TString::Format("DijetMassBkg_cat%d_min",cat),TString::Format("DijetMassBkg_cat%d_min",cat), *var,*rooData_all, 3);
  /* RooPlot* plot = var->frame();
  rooData_all->plotOn(plot);
  rooFunc_all->plotOn(plot);
  plot->Draw("same");
  // c->SetLogy();

  c->SaveAs("plots/.png");
  c->SaveAs("plots/Bkg_fitMinimum.pdf");
*/
  return rooFunc_all;

}



TF1* getFuncFitMIN(int cat, RooRealVar* var){

  double mass[1] = {2000.};
  double c0[1] = {2000.};
  double c1[1] = {2000.};
  double c2[1] = {2000.};
  double c3[1] = {2000};
 
 
  TF1* f = new TF1("f", "[0]+[1]*x", 150, 850);
  f->SetParameter(0,1000 );
  f->SetParameter(1, -0.5);
  f->SetParameter(2, -1000.);
  
  
  TH1F* h_all = new TH1F("h_all", "h_all", 1, 0, 10000);
  for(int i=0;i<1;i++){
    if(cat==0) h_all->SetBinContent(h_all->FindBin(mass[i]),c0[i]);
    if(cat==1) h_all->SetBinContent(h_all->FindBin(mass[i]),c1[i]);
    if(cat==2) h_all->SetBinContent(h_all->FindBin(mass[i]),c2[i]);
    if(cat==3) h_all->SetBinContent(h_all->FindBin(mass[i]),c3[i]);
      
  }
  TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();
  //TPaveText* label_cms = get_labelCMS(0,"2014", false);
  //TPaveText* label_sqrt = get_labelSqrt(0);  

  h_all->GetXaxis()->SetTitle("m_{X} [GeV]");
  h_all->GetYaxis()->SetTitle("Fit Range: Minimum");
  h_all->Draw("");
  h_all->Fit("f");
  f->Draw("Lsame");
  for(int i=0;i<1;i++){
    std::cout<<cat<<std::endl;
    std::cout<<mass[i]<<"  "<<h_all->FindBin(mass[i])<<"   "<<c0[i]<< "    "<<f->Eval(mass[i])<<std::endl;
       
  }
 
  return f;

}






RooHistFunc* getRooHistFuncFitMAX(int cat, RooRealVar* var){


  double mass[1] = {2000.};
  double c0[1] = {2000.};
  double c1[1] = {2000.};
  double c2[1] = {2000.};
  double c3[1] = {2000};


  TH1F* h_all = new TH1F("h_all", "h_all", 1, 0, 10000);
  for(int i=0;i<1;i++){
     std::cout<<cat<<std::endl;
     std::cout<<mass[i]<<"  "<<h_all->FindBin(mass[i])<<"   "<<c0[i]<<std::endl;
     if(cat==0) h_all->SetBinContent(h_all->FindBin(mass[i]),c0[i]);
     if(cat==1) h_all->SetBinContent(h_all->FindBin(mass[i]),c1[i]);
     if(cat==2) h_all->SetBinContent(h_all->FindBin(mass[i]),c2[i]);
     if(cat==3) h_all->SetBinContent(h_all->FindBin(mass[i]),c3[i]);
      
  }
  /* TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();
  //TPaveText* label_cms = get_labelCMS(0,"2014", false);
  //TPaveText* label_sqrt = get_labelSqrt(0);  

  h_all->GetXaxis()->SetTitle("m_{X} [GeV]");
  h_all->GetYaxis()->SetTitle("Fit Range: Maximum");
  h_all->Draw();
  */

  RooDataHist* rooData_all = new RooDataHist("rooData_all", "roData_all",*var,h_all);
  RooHistFunc* rooFunc_all = new RooHistFunc(TString::Format("DijetMassBkg_cat%d_max",cat),TString::Format("DijetMassBkg_cat%d_max",cat), *var,*rooData_all, 7);

  /* RooPlot* plot = var->frame();
  rooData_all->plotOn(plot);
  rooFunc_all->plotOn(plot);
  plot->Draw("same");
  c->SetLogy();

  c->SaveAs("plots/.png");
  c->SaveAs("plots/Bkg_fitMaximum.pdf");*/

  return rooFunc_all;

}



TF1* getFuncFitMAX(int cat, RooRealVar* var){

  double mass[1] = {2000.};
  double c0[1] = {2000.};
  double c1[1] = {2000.};
  double c2[1] = {2000.};
  double c3[1] = {2000};

  TF1* f = new TF1("f", "[0]+[2]*pow(x,[1])", 150, 850);
  f->SetParameter(0,1000 );
  f->SetParameter(1, -0.5);
  f->SetParameter(2, -1000.);
  

  TH1F* h_all = new TH1F("h_all", "h_all", 1, 0, 10000);
  for(int i=0;i<1;i++){
    if(cat==0) h_all->SetBinContent(h_all->FindBin(mass[i]),c0[i]);
    if(cat==1) h_all->SetBinContent(h_all->FindBin(mass[i]),c1[i]);
    if(cat==2) h_all->SetBinContent(h_all->FindBin(mass[i]),c2[i]);
    if(cat==3) h_all->SetBinContent(h_all->FindBin(mass[i]),c3[i]);
      
  }
  TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();
  //TPaveText* label_cms = get_labelCMS(0,"2014", false);
  //TPaveText* label_sqrt = get_labelSqrt(0);  

  h_all->GetXaxis()->SetTitle("m_{X} [GeV]");
  h_all->GetYaxis()->SetTitle("Fit Range: Maximum");
  h_all->Draw();
  h_all->Fit("f");
  f->Draw("Lsame");
  for(int i=0;i<1;i++){
    std::cout<<cat<<std::endl;
    std::cout<<mass[i]<<"  "<<h_all->FindBin(mass[i])<<"   "<<c0[i]<< "    "<<f->Eval(mass[i])<<std::endl;
       
  }
  

 
  return f;

}


