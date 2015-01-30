
RooHistFunc* getRooHistFunc(int cat, RooRealVar* var, std::string model){


  double mass[8] = {2000., 3000., 4000., 5000., 6000., 7000., 8000., 9000.};

  double c0[8];
  double c1[8];
  double c2[8];
  double c3[8];
  double all[8];
  
  TH1F* h_all = new TH1F("h_all", "h_all", 1, 1100., 6000.);
  
  for(int i=0;i<5;i++){
    std::cout<<"i: "<<i<<" mass: "<<mass[i]<<std::endl;
    if(cat==0) h_all->SetBinContent(h_all->FindBin(mass[i]),c0[i]);
    if(cat==1) h_all->SetBinContent(h_all->FindBin(mass[i]),c1[i]);
    if(cat==2) h_all->SetBinContent(h_all->FindBin(mass[i]),c2[i]);
    if(cat==3) h_all->SetBinContent(h_all->FindBin(mass[i]),c3[i]);
    if(cat==4) h_all->SetBinContent(h_all->FindBin(mass[i]),all[i]);
   
  }


  /*  
  for(int i=0;i<16;i++){
    std::cout<<cat<<std::endl;
    // std::cout<<mass[i]<<"  "<<h_all->FindBin(mass[i])<<"   "<<c0[i]<<std::endl;
    std::cout<<i<<"  "<<h_all->GetBinCenter(i)<<"  "<<h_all->GetBinContent(i)<<std::endl;

    }
  */

  TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();
  //  TPaveText* label_cms = get_labelCMS(0, false);
  //TPaveText* label_sqrt = get_labelSqrt(0);  
  h_all->GetYaxis()->SetRangeUser(0.0001, 0.1);
  h_all->GetXaxis()->SetTitle("m_{X} [GeV]");
  h_all->GetYaxis()->SetTitle("Signal Yield");
  h_all->Draw();

  RooDataHist* rooData_all = new RooDataHist("rooData_all", "roData_all",*var,h_all);
  RooHistFunc* rooFunc_all = new RooHistFunc(TString::Format("DijetMassSig_cat%d_norm",cat),TString::Format("DijetMassSig_cat%d_norm",cat), *var,*rooData_all, 3);
  /*  RooPlot* plot = var->frame();
  rooData_all->plotOn(plot);
  rooFunc_all->plotOn(plot);
  plot->Draw("same");
  c->SetLogy();

  c->SaveAs("plots/signalYield.png");
  c->SaveAs("plots/signalYield.pdf");*/
  
  return rooFunc_all;

}

