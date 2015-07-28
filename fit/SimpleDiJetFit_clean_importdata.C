#include <algorithm>

void SimpleDiJetFitV1()
{

  gROOT->Reset();

  // Set TDR Style	
  gROOT->ProcessLine(".L tdrStyle.C");
  gROOT->ProcessLine("setTDRStyle()");
  gROOT->LoadMacro("CMS_lumi.C");

  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_13TeV = "37 pb^{-1}" ;
  lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"

  int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV  
  //   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)
  int iPos = 10;
  //gStyle->SetOptFit(0); 
  //gROOT->ForceStyle();

  //gRandom = new TRandom3(0);
  //gRandom->SetSeed(0);

  double lumi = 37;
  double sigmaQstar = 0.2827E-01;

  //########################
  //##### User Options #####
  //########################
   //#### signal #######
  //TFile* signal_sample_giulia = new TFile("histo_signal_qstar_5000.root","read"); 
  //TH1F* h_w_giulia = signal_sample_giulia->Get("h_w");
  //h_w_giulia->Scale(0.0182*1000/100000);
  //h_w_giulia->SetName("h_w_giulia");
  //TFile* signal_sample = new TFile("/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/fit/Resonance_Shapes_qg_13TeV_newJEC.root", "read");
  TFile* signal_sample = new TFile("/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/test_fit/Resonance_Shapes_qg_PU20_13TeV_newJEC.root", "read");
  TH1F* h_sig = (TH1F*)signal_sample->Get("h_qg_4500");
  const Int_t nbins_sig = h_sig->GetNbinsX();
  double integral = h_sig->Integral();
  h_sig->Scale( sigmaQstar / integral);
  //cout << "n bins x sig : " << nbins_sig << endl;
  double massBinsSig[nbins_sig+1];
//-------------
  for (int i=0; i<nbins_sig+1; i++ ){
    massBinsSig[i]=h_sig->GetXaxis()->GetBinLowEdge(i+1);  
  }

  TH1D* h_w = new TH1D("h_w","", h_sig->GetNbinsX(), massBinsSig);
  for (int i=1; i<nbins_sig+1; i++ ){
    double bincontent = h_sig->GetBinContent(i) / h_sig->GetBinWidth(i);
    //cout << "content bin " << i << " = " << bincontent << endl;
    h_w->SetBinContent(i, bincontent );
  } 
  h_w->Print();
  //cout << "maximum bin h_sig : " << h_sig->GetMaximumBin() << "  maximum bin h_w: " << h_w->GetMaximumBin() << endl; 
  //cout << "maximum bin h_sig center: " << h_sig->GetBinCenter(h_sig->GetMaximumBin()) << "  maximum bin h_w center: " << h_w->GetBinCenter(h_w->GetMaximumBin()) << endl; 
  
  //pseudodata
  //char input_root_file[500] = "/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/test_fit/dijetFitResults_FuncType0_nParFit4_MC_1fb-1_Dinko.root";
  
  // data all  
  //char input_root_file[500] = "/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_run245155/histo_data_mjj_fromTree.root";
  //char input_root_file[500] = "../scripts/histo_data_mjj_fromTree_run246908_247398.root";
  //char input_root_file[500] = "/cmshome/fpreiato/CMSSW_7_4_3/src/CMSDIJET/DijetRootTreeAnalyzer/output/Data2015/UnstableBeam/Plot/histo_data_mjj_fromTree.root";
  //char input_root_file[500] = "../scripts/histo_data_mjj_fromTree_run246908_247398.root";
  //char input_root_file[500] = "/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_Complete0T/histo_data_mjj_fromTree.root";
char input_root_file[500] = "/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_withSF_19_07_15/all/histo_data_mjj_fromTree.root";

  char fileNameSuffix[500] = "data_nosig"; //i.e. run period
  //char fileNameSuffix[500] = "MC_10fb-1"; //i.e. run period

  //char input_1Dhistogram[500] = "hist_allCutsQCD";
  //char input_1Dhistogram[500] = "hist_mass_1GeV";
  char input_1Dhistogram[500] = "h_dat";
  double minX_mass = 1118.;
  //1 fb -1
  //double maxX_mass = 6099.;
  //10 fb-1
  //double maxX_mass = 7000.;
  ////////
  double maxX_mass = 5253.;

  //Fit functions
  // 0: DEFAULT (4 par.) - "( [0]*TMath::Power(1-x/8000,[1]) ) / ( TMath::Power(x/8000,[2]+[3]*log(x/8000)) )" 

  const int FunctionType = 0;

  const int number_of_variableWidth_bins = 103;

  double massBins[number_of_variableWidth_bins+1] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325,
    354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687,
    1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509,
    4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 
    10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000};

  //================================================================================================================

  int nPar=-1; 

  //### input file and 1D histo
  TFile *file0=TFile::Open( input_root_file );
  TH1D* hist_mass_original = (TH1D*) file0->Get( input_1Dhistogram );  

  float bin_width_X =  hist_mass_original->GetBinWidth(1);
  cout << "bin_width_X: " << bin_width_X << endl;

  //variable binning random dataset
  hist_binned= (TH1F*)hist_mass_original->Rebin(number_of_variableWidth_bins,"hist_binned",massBins); 

  //variable binning pseudodataset reweighted by the bin width
  TH1F* hist_mass = (TH1F*)hist_binned->Clone("hist_mass");
  hist_mass->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  for (int i=1; i<= number_of_variableWidth_bins; i++){
    bincontent = hist_binned->GetBinContent(i);
    binwidth = hist_binned->GetBinWidth(i);
    binerror = hist_binned->GetBinError(i);
    hist_mass->SetBinContent(i,bincontent/(binwidth*lumi));   
    hist_mass->SetBinError(i,binerror/(binwidth*lumi));   
    cout << "content bin " << i << " = " <<  hist_mass->GetBinContent(i) << endl;
  }

  //hist_mass->Draw();


  //### fit mass histogram with background function
  TF1 *M1Bkg;

  // 0: DEFAULT (4 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )" 
  if( FunctionType==0 )    
    {
      nPar=4;
      M1Bkg = new TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )",minX_mass,14000);
      // M1Bkg->SetParameter(0,0.05);
      // M1Bkg->SetParameter(1,7.1);
      // M1Bkg->SetParameter(2,5.9);
      // M1Bkg->SetParameter(3,0.2);
      //1 fb-1
      //M1Bkg->SetParameter(0,0.00073);
      //10 fb-1
      //M1Bkg->SetParameter(0,0.0073);
      // pb-1
      M1Bkg->SetParameter(0,0.08);
      M1Bkg->SetParameter(1,28);
      M1Bkg->SetParameter(2,2.8);
      M1Bkg->SetParameter(3,0.2595);

      //1 fb-1
      M1Bkg->SetParLimits(0,0.,10);
      //10 fb-1
      //M1Bkg->SetParLimits(0,0,1.);
      M1Bkg->SetParLimits(1,0.,10000);
      M1Bkg->SetParLimits(2,1.,10000);
      M1Bkg->SetParLimits(3,0.,10000);

      //M1Bkg->FixParameter(3,0.);
    }



  //////////////////////////////////////////////////////////////////////////////////////////////// 
  //  ______ _ _                     _       _     _        _     _             _               
  //  |  ___(_) |                   (_)     | |   | |      | |   (_)           (_)            
  //  | |_   _| |_  __   ____ _ _ __ _  __ _| |__ | | ___  | |__  _ _ __  _ __  _ _ __   __ _ 
  //  |  _| | | __| \ \ / / _` | '__| |/ _` | '_ \| |/ _ \ | '_ \| | '_ \| '_ \| | '_ \ / _` |
  //  | |   | | |_   \ V / (_| | |  | | (_| | |_) | |  __/ | |_) | | | | | | | | | | | | (_| |
  //  \_|   |_|\__|   \_/ \__,_|_|  |_|\__,_|_.__/|_|\___| |_.__/|_|_| |_|_| |_|_|_| |_|\__, |
  //                                                                                     __/ |
  //                                                                                     |___/ 
  // 
  /////////////////////////////////////////////////////////////////////////////////////////////////  

  cout << hist_mass_original->Integral(hist_mass_original->FindBin(3000),hist_mass_original->FindBin(13000)) << endl;
  TFitResultPtr r_bin;
  r_bin->SetName("fit_result_M1Bkg_binned");
  //TFile* fitresult_binned = new TFile("fitresult_binned_1fb-1.root","recreate");
  //fitresult_binned->cd();
  int stopProgram=1;
  for( int loop=0; loop<10; loop++)
    {
      //r = hist_mass->Fit("M1Bkg","LS","",minX_mass,maxX_mass);      
      r_bin = hist_mass->Fit("M1Bkg","ILSR","",minX_mass,maxX_mass);      
      Int_t fitStatus = r_bin;
      if(fitStatus==0)
	{
	  stopProgram=0;
	  r_bin->Print("V"); 
	  //r_bin->Write();
	  break;
	}
    }
  //fitresult_binned->Close();

  if(stopProgram==1)
    {
      cout << "######################" << endl;
      cout << "######################" << endl;
      cout << "ERROR : Fit failed!!!!" << endl;
      cout << "######################" << endl;
      cout << "######################" << endl;
      break;
    }

//  // fit residuals and chi2 from histogram
//  TH1D* hist_fit_residual_vsMass = new TH1D("hist_fit_residual_vsMass","hist_fit_residual_vsMass",number_of_variableWidth_bins,massBins);
//  TH1D* hist_fit_residual = new TH1D("hist_fit_residual","hist_fit_residual",10,-5,5);
//  int NumberOfObservations_VarBin = 0;
//  double chi2_VarBin = 0.;
//
//  //cout << hist_mass_varbin->GetNbinsX() << endl;  
//  for(int bin=1; bin<number_of_variableWidth_bins; bin++)
//    {
//      if( hist_mass->GetXaxis()->GetBinLowEdge(bin)>=minX_mass 
// 	  && hist_mass->GetXaxis()->GetBinUpEdge(bin)<=maxX_mass )
// 	{	  
//	  NumberOfObservations_VarBin++;
//
// 	  //cout << hist_mass_varbin->GetXaxis()->GetBinLowEdge(bin) << endl;
//	  double data = hist_mass->GetBinContent(bin);
//	  double err_data = hist_mass->GetBinError(bin);
//	  if( data == 0 )
//	    {
//	      err_data = 1.8 / (hist_mass->GetBinWidth(bin) * lumi) ; //giulia : perche`??
//	      hist_mass->SetBinError(bin, err_data);
//	    }
//	  //double fit = M1Bkg->Eval(hist_mass->GetBinCenter(bin));
// 	  double fit = M1Bkg->Integral(hist_mass->GetXaxis()->GetBinLowEdge(bin) 
// 	  			       , hist_mass->GetXaxis()->GetBinUpEdge(bin) ); 
//	  fit = fit / ( hist_mass->GetBinWidth(bin) );
//	  
//	  double err_tot = err_data;	  
//	  double fit_residual = (data - fit) / err_tot;
//	  double err_fit_residual = 1;
//	    	  
//	  chi2_VarBin += pow( (data - fit) , 2 ) / pow( err_data , 2 );	 
//
//	  // 	  cout << "data, err_data, fit: " << data << ", " << err_data << ", " << fit << endl;
//	  // 	  cout << "bin, fit residual : " << bin << ", " <<fit_residual << endl;	  
//	  hist_fit_residual_vsMass->SetBinContent(bin,fit_residual);
//	  hist_fit_residual_vsMass->SetBinError(bin,err_fit_residual);
//	  hist_fit_residual->Fill(fit_residual);
// 	}
//    }
//  int ndf_VarBin = NumberOfObservations_VarBin - nPar  ;
//  cout << "============================" << endl;
//  cout << "NumberOfObservations_VarBin: " << NumberOfObservations_VarBin << endl;
//  cout << "ndf_VarBin: " << ndf_VarBin << endl;
//  cout << "chi2_VarBin: " << chi2_VarBin << endl;
//  cout << "============================" << endl;  


//#######################################################
// data in TGraph format (hist binned)
  const double alpha = 1 - 0.6827;
  double vx[number_of_variableWidth_bins];
  double vy[number_of_variableWidth_bins];
  double vexl[number_of_variableWidth_bins];
  double vexh[number_of_variableWidth_bins];
  double veyl[number_of_variableWidth_bins];
  double veyh[number_of_variableWidth_bins];
  
  for(int i=0; i<number_of_variableWidth_bins; ++i)
    {
      double n    = hist_binned->GetBinContent(i+1);
      double dm   = hist_binned->GetBinWidth(i+1);
      double mass = hist_binned->GetBinCenter(i+1);
      double xl   = hist_binned->GetBinLowEdge(i+1);
      double xh   = xl+dm;
      vx[i]   = (xl+xh)/2.;
      vexl[i] = dm/2.;
      vexh[i] = dm/2.;
      vy[i]   = n / (dm*lumi);

      double l = 0.5*TMath::ChisquareQuantile(alpha/2,2*n);
      //double h = (n==0) ? ( 0.5*TMath::ChisquareQuantile(1-alpha,2*(n+1)) ) : ( 0.5*TMath::ChisquareQuantile(1-alpha/2,2*(n+1)) );
      double h = 0.5*TMath::ChisquareQuantile(1-alpha/2,2*(n+1));

      veyl[i] = (n-l)/(lumi*dm);
      veyh[i] = (h-n)/(lumi*dm);
    }
  //output txt file
  char output_txt_file[500];
  sprintf(output_txt_file,"DijetMassSpectrum_CMS_13TeV_36pb-1.txt"); 
  ofstream myfile;
  myfile.open(output_txt_file);
  myfile << "########################################################################################" << endl;
  myfile << "#binLow(GeV)    binHigh(GeV)    dsdm(pb/GeV)    errLow_dsdm(pb/GeV)    errUp_dsdm(pb/GeV)" << endl;
  myfile << "########################################################################################" << endl;

  // data in the graph format
  TGraphAsymmErrors *g = new TGraphAsymmErrors(number_of_variableWidth_bins,vx,vy,vexl,vexh,veyl,veyh);

//#######################################################
  // fit residuals and chi2 from TGraphAsymmErrors
  TH1D* hist_fit_residual_vsMass = new TH1D("hist_fit_residual_vsMass","hist_fit_residual_vsMass",number_of_variableWidth_bins,massBins);
  TH1D* hist_fit_residual = new TH1D("hist_fit_residual","hist_fit_residual",10,-5,5);
  int NumberOfObservations_VarBin = 0;
  double chi2_VarBin = 0.;
  double integral_4p2TeV = M1Bkg->Integral(4200, 13000 );
  double integral_3TeV = M1Bkg->Integral(3000, 13000 );
  double integral_3p5TeV = M1Bkg->Integral(3500, 13000 );
  double integral_4TeV = M1Bkg->Integral(4000, 13000 );
  double integral_3TeV_sig= 0;
  double integral_3p5TeV_sig=0; 
  double integral_4TeV_sig =0; 
  double integral_3TeV_data= 0;
  double integral_3p5TeV_data=0; 
  double integral_4TeV_data =0; 
  

  for(int bin=1; bin<number_of_variableWidth_bins; bin++)
    {
      if (hist_mass->GetXaxis()->GetBinLowEdge(bin)>3000 ){ 
	integral_3TeV_sig+= h_w->GetBinContent(bin)*h_w->GetBinWidth(bin);
	integral_3TeV_data+= hist_mass->GetBinContent(bin)*h_w->GetBinWidth(bin);
      }
  	if (hist_mass->GetXaxis()->GetBinLowEdge(bin)>3500 ) {
	integral_3p5TeV_sig+= h_w->GetBinContent(bin)*h_w->GetBinWidth(bin);
	integral_3p5TeV_data+= hist_mass->GetBinContent(bin)*h_w->GetBinWidth(bin);
	}
      if (hist_mass->GetXaxis()->GetBinLowEdge(bin)>4000 )  {
  	integral_4TeV_sig+=h_w->GetBinContent(bin)*h_w->GetBinWidth(bin);
	integral_4TeV_data+= hist_mass->GetBinContent(bin)*h_w->GetBinWidth(bin);
      }
      if( hist_mass->GetXaxis()->GetBinLowEdge(bin)>=minX_mass 
 	  && hist_mass->GetXaxis()->GetBinUpEdge(bin)<=maxX_mass )
 	{	  
	  NumberOfObservations_VarBin++;
	  cout << "bin content = " << hist_mass->GetBinContent(bin) << "   graph y = " << data[bin-1] << "  error y low = " << g->GetErrorYlow(bin-1) << endl ;
	  double data = hist_mass->GetBinContent(bin);
	  double err_data_low = g->GetErrorYlow(bin-1);
	  double err_data_high= g->GetErrorYhigh(bin-1);

  	  double fit = M1Bkg->Integral(hist_mass->GetXaxis()->GetBinLowEdge(bin) 
				       , hist_mass->GetXaxis()->GetBinUpEdge(bin) ); 
	  fit = fit / ( hist_mass->GetBinWidth(bin) );
	  double err_tot;
	  if(fit > data) err_tot = err_data_high;
	  else err_tot = err_data_low; 

	  double fit_residual = (data - fit) / err_tot;
	  double err_fit_residual = 1;
	    	  
	  chi2_VarBin += pow( (data - fit) , 2 ) / pow( err_tot , 2 );	 

	  // 	  cout << "data, err_tot, fit: " << data << ", " << err_tot << ", " << fit << endl;
	  // 	  cout << "bin, fit residual : " << bin << ", " <<fit_residual << endl;	  
	  hist_fit_residual_vsMass->SetBinContent(bin,fit_residual);
	  hist_fit_residual_vsMass->SetBinError(bin,err_fit_residual);
	  hist_fit_residual->Fill(fit_residual);

	  //write on txt file
	  myfile << hist_mass->GetXaxis()->GetBinLowEdge(bin) << "      " << hist_mass->GetXaxis()->GetBinUpEdge(bin) << "      " << data << "      " << err_data_low << "       " << err_data_high  << endl;
    
	}
  }
  myfile.close();
  
  int ndf_VarBin = NumberOfObservations_VarBin - nPar  ;
  cout << "============================" << endl;
  cout << "NumberOfObservations_VarBin: " << NumberOfObservations_VarBin << endl;
  cout << "ndf_VarBin: " << ndf_VarBin << endl;
  cout << "chi2_VarBin: " << chi2_VarBin << endl;
  cout << "integral fit func > 4.2 TeV " << integral_4p2TeV*lumi << endl;
  cout << "integral fit func > 3.0 TeV " << integral_3TeV*lumi << endl;
  cout << "integral fit func > 3.5 TeV " << integral_3p5TeV*lumi << endl;
  cout << "integral fit func > 4.0 TeV " << integral_4TeV*lumi << endl;
  cout << "integral sig  > 3.0 TeV " << integral_3TeV_sig*lumi << endl;
  cout << "integral sig  > 3.5 TeV " << integral_3p5TeV_sig*lumi << endl;
  cout << "integral sig > 4.0 TeV " << integral_4TeV_sig*lumi << endl;
  cout << "integral data > 3.0 TeV " << integral_3TeV_data*lumi << endl;
  cout << "integral data  > 3.5 TeV " << integral_3p5TeV_data*lumi << endl;
  cout << "integral data > 4.0 TeV " << integral_4TeV_data*lumi << endl;
  cout << "============================" << endl;  


  //### Draw plots
//
//  TCanvas *Canvas0 = new TCanvas("Canvas0","Canvas0",11,51,700,500);
//  Canvas0->cd();
//  Canvas0->SetLogy();
//  hist_mass_original->GetYaxis()->SetTitle("Events");
//  hist_mass_original->Draw();  
//  M1Bkg->SetLineColor(2);
//  M1Bkg->Draw("same");     
//
//  TCanvas *Canvas1 = new TCanvas("Canvas1","Canvas1",11,51,700,500);
//  Canvas1->cd();
//  Canvas1->SetLogy();
//  //hist_mass->GetYaxis()->SetTitle("Events / bin width");
//  hist_mass->GetYaxis()->SetTitle("d#sigma / dm (pb/GeV)");
//  hist_mass->Draw("PE0");  
//  M1Bkg->SetLineColor(2);
//  M1Bkg->Draw("same");     
//
//  TCanvas *Canvas2 = new TCanvas("Canvas2","Canvas2",11,51,700,500);
//  Canvas2->cd();
//  Canvas2->SetGridx();
//  Canvas2->SetGridy();
//  Canvas2->SetLogx();
//  hist_fit_residual_vsMass->GetYaxis()->SetLimits(-5,5);
//  hist_fit_residual_vsMass->GetYaxis()->SetRangeUser(-5,5);
//  hist_fit_residual_vsMass->GetYaxis()->SetTitle("(data - fit) / #sqrt{data}");
//  hist_fit_residual_vsMass->GetXaxis()->SetRangeUser(minX_mass,maxX_mass);
//  //hist_fit_residual_vsMass->GetXaxis()->SetRangeUser(1000,10000);
//  hist_fit_residual_vsMass->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
//  hist_fit_residual_vsMass->Draw();
//
//  TCanvas *Canvas3 = new TCanvas("Canvas3","Canvas3",11,51,700,500);
//  Canvas3->cd();
//  hist_fit_residual->GetXaxis()->SetTitle("(data - fit) / #sqrt{data}");
//  hist_fit_residual->GetYaxis()->SetTitle("Number of bins");
//  hist_fit_residual->GetYaxis()->SetRangeUser(0,number_of_variableWidth_bins/3);
//  hist_fit_residual->Draw();
//  hist_fit_residual->Fit("gaus","L","",-3,3);
//

//#################### canvas "note" style ####################
 
  int W = 600;
  int H = 650;
  int H_ref = 650; 
  int W_ref = 600; 
  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;



  TCanvas* c = new TCanvas("c","DijetMass cross section with Fit and QCD MC",W,H);
  c->GetWindowHeight();
  c->GetWindowWidth();
  c->SetLogy();
  c->Divide(1,2,0,0,0);


  // ------------ pad 1  ----------------
  c->cd(1);
  p11_1 = (TPad*)c->GetPad(1);
  p11_1->SetPad(0.01,0.23,0.99,0.98);
  p11_1->SetLogy();
  p11_1->SetRightMargin(0.05);
  p11_1->SetTopMargin(0.05);
  p11_1->SetFillColor(0);
  p11_1->SetBorderMode(0);
  p11_1->SetFrameFillStyle(0);
  p11_1->SetFrameBorderMode(0);
//  p11_1->SetLeftMargin( L/W );
//  p11_1->SetRightMargin( R/W );
//  p11_1->SetTopMargin( T/H );
//  p11_1->SetBottomMargin( B/H );
//  p11_1->SetTickx(0);
//  p11_1->SetTicky(0);

  // Pave text
//  TPaveText *pave_fit = new TPaveText(0.1558691,0.30735043,0.3750171,0.4070085,"NDC");
  TPaveText *pave_fit = new TPaveText(0.2058691,0.20735043,0.4750171,0.3670085,"NDC");
  
  //pave_fit->AddText(" #sqrt{s} = 13 TeV");
  pave_fit->AddText("|#eta| < 2.5, |#Delta#eta| < 1.3");
  pave_fit->AddText("M_{jj} > 1.1 TeV");
  pave_fit->AddText("Wide Jets");
  pave_fit->SetFillColor(0);
  pave_fit->SetLineColor(0);
  pave_fit->SetFillStyle(0);
  pave_fit->SetBorderSize(0);
  pave_fit->SetTextFont(42);
  pave_fit->SetTextSize(0.040);
  pave_fit->SetTextAlign(12); 


//  TPaveText *pt1 = new TPaveText(0.1284756,0.9602144,0.3887139,0.9902251,"brNDC");
//  pt1->SetBorderSize(0);
//  pt1->SetFillColor(0);
//  pt1->SetFillStyle(0);
//  pt1->SetLineColor(0);
//  pt1->SetTextAlign(12);
//  pt1->SetTextSize(0.035);
//  TText *text = pt1->AddText("CMS Preliminary");
//
//  TPaveText *pt2 = new TPaveText(0.45,0.96,0.65,0.99,"brNDC");
//  pt2->SetBorderSize(0);
//  pt2->SetFillColor(0);
//  pt2->SetFillStyle(0);
//  pt2->SetLineColor(0);
//  pt2->SetTextAlign(12);
//  pt2->SetTextSize(0.035);
//  TText *text2 = pt2->AddText("#sqrt{s} = 13 TeV");
//
//  TPaveText *pt3 = new TPaveText(0.7687988,0.9602144,0.9297357,0.9902251,"brNDC");
//  pt3->SetBorderSize(0);
//  pt3->SetFillColor(0);
//  pt3->SetFillStyle(0);
//  pt3->SetLineColor(0);
//  pt3->SetTextAlign(12);
//  pt3->SetTextSize(0.035);
//  TText *text3 = pt3->AddText("L= 21.239 pb^{-1}");
//  //TText *text3 = pt3->AddText("L= 1 fb^{-1}");
//
  TH1F *vFrame = p11_1->DrawFrame(minX_mass,0.000005,maxX_mass,5.0);

  vFrame->SetTitle("");
  vFrame->SetXTitle("Dijet Mass (GeV)");
  vFrame->SetYTitle("d#sigma / dm_{jj}   (pb / GeV)");
  vFrame->GetXaxis()->SetTitleSize(0.06);
  vFrame->GetXaxis()->SetTitleOffset(0.95);
  vFrame->GetXaxis()->SetLabelSize(0.05);
  vFrame->GetYaxis()->SetTitleSize(0.06);
  vFrame->GetYaxis()->SetTitleOffset(0.95);
  vFrame->GetYaxis()->SetLabelSize(0.05);

  

//  hist_mass->GetXaxis()->SetRangeUser(minX_mass,maxX_mass);
//  hist_mass->SetTitle("");
//  hist_mass->SetLineColor(1);
//  hist_mass->SetFillColor(1);
//  hist_mass->SetMarkerColor(1);
//  hist_mass->SetMarkerStyle(20);
//  hist_mass->Draw("PE0");

//  gStyle->SetEndErrorSize(0);
  g->SetMarkerSize(0.9);
  g->Draw("PE0");
  M1Bkg->SetLineWidth(2);
  M1Bkg->SetLineStyle(2);
  M1Bkg->SetLineColor(2);
  M1Bkg->Draw("same");
  
  h_w->SetLineColor(kBlue);
  h_w->SetLineWidth(2);
  h_w->SetLineStyle(2);
  //h_w->Draw("same");

  //TLegend *leg = new TLegend(0.5564991,0.4,0.8903575,0.575812);
  TLegend *leg = new TLegend(0.5564991,0.55,0.8903575,0.80);
  leg->SetTextSize(0.03546853);
  leg->SetLineColor(0);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetMargin(0.35);
  leg->AddEntry(hist_mass,"data" ,"PL");
  leg->AddEntry(M1Bkg,"fit to data","L");
  //leg->AddEntry(h_w,"q* (4.5 TeV)","L");
  leg->Draw("same");

//  pt1->Draw("same"); 
//  pt2->Draw("same");
//  pt3->Draw("same");
  pave_fit->Draw("same");

  // writing the lumi information and the CMS "logo"
  CMS_lumi( p11_1, iPeriod, iPos );
  //redraw axis
  p11_1->RedrawAxis();
  p11_1->Update();
  p11_1->GetFrame()->Draw();
  //cout << "MIN: " << p11_1->GetUxmin() << endl;
  //cout << "MAX: " << p11_1->GetUxmax() << endl;

//  //---- Next PAD

  c->cd(2);
  p11_2 = (TPad*)c->GetPad(2);
  p11_2->SetPad(0.01,0.02,0.99,0.24);
  p11_2->SetBottomMargin(0.35);
  p11_2->SetRightMargin(0.05);
  p11_2->SetGridx();
  p11_2->SetGridy();
  //c3_2->SetTickx(50);

  TH1F *vFrame2 = p11_2->DrawFrame(p11_1->GetUxmin(), -3., p11_1->GetUxmax(), 3.);

  vFrame2->SetTitle("");
  vFrame2->SetXTitle("Dijet Mass (GeV)");
  vFrame2->GetXaxis()->SetTitleSize(0.06);
  vFrame2->SetYTitle("(Data-Fit)/#sigma");
  vFrame2->GetYaxis()->SetTitleSize(0.15);
  vFrame2->GetYaxis()->SetTitleOffset(0.40);
  vFrame2->GetYaxis()->SetLabelSize(0.09);
  vFrame2->GetXaxis()->SetTitleSize(0.18);
  vFrame2->GetXaxis()->SetTitleOffset(0.90);
  vFrame2->GetXaxis()->SetLabelSize(0.15);

  hist_fit_residual_vsMass->GetXaxis()->SetRangeUser(minX_mass,maxX_mass);
  hist_fit_residual_vsMass->GetYaxis()->SetRangeUser(-4.,4.);
  hist_fit_residual_vsMass->SetLineWidth(0);
  hist_fit_residual_vsMass->SetFillColor(2);
  hist_fit_residual_vsMass->SetLineColor(1);
  hist_fit_residual_vsMass->Draw("SAMEHIST");

  //histogram signal significance 
  TH1D* hist_sig_significance = new TH1D("hist_sig_significance","hist_sig_significance",nbins_sig,massBinsSig);
  
  //cout << "         Nbins    bin low edge " << endl; 
  for(int i=1; i<number_of_variableWidth_bins+1; i++){
    double fit = M1Bkg->Integral(hist_mass->GetXaxis()->GetBinLowEdge(i), hist_mass->GetXaxis()->GetBinUpEdge(i) );
    double significance = 0;

    if((h_sig->GetBinContent(i)+fit) != 0 && h_sig->GetBinLowEdge(i)<12000){
      significance  = h_sig->GetBinContent(i) / TMath::Sqrt(h_sig->GetBinContent(i) + fit );
    }
    //cout << "sig: " << h_sig->GetBinContent(i) << "  bkg: " <<  hist_mass->GetBinContent(i) << endl; 
    cout <<  "low edge bin: " << h_sig->GetBinLowEdge(i) << "  " << fit << " + " <<  h_sig->GetBinContent(i) << "  sqrt(fit + sig) = " << TMath::Sqrt(h_sig->GetBinContent(i) + fit)  << endl;
    cout << "significance = " << significance << endl;
    hist_sig_significance->SetBinContent(i,significance);

  }
  hist_sig_significance->SetLineColor(kBlue) ;
  hist_sig_significance->SetLineWidth(2) ;
  hist_sig_significance->SetLineStyle(2) ;
//  hist_sig_significance->Draw("same");
  
  
  
  TLine *line = new TLine(minX_mass,0,maxX_mass,0);
  line->Draw("");
  //c->Close();

  //### Output files

  char output_root_file[500];
  sprintf(output_root_file,"dijetFitResults_FuncType%d_nParFit%d_%s.root",FunctionType,nPar,fileNameSuffix); 

  TFile f_output(output_root_file,"RECREATE");
  f_output.cd();
  //Canvas0->Write();
  //Canvas1->Write();
  //Canvas2->Write();
  //Canvas3->Write();
  hist_mass_original->Write();
  hist_binned->Write();
  hist_mass->Write();
  c->Write();
  //r_bin->Write();
  f_output.Close();

  //### Save figures from canvas
  char c0_fileName[200];
  sprintf(c0_fileName,"dijetmass_FuncType%d_nParFit%d_%s.png",FunctionType,nPar,fileNameSuffix);
  char c1_fileName[200];
  sprintf(c1_fileName,"dijetmass_varbin_FuncType%d_nParFit%d_%s.png",FunctionType,nPar,fileNameSuffix);
  char c2_fileName[200];
  sprintf(c2_fileName,"fitresiduals_vs_mass_FuncType%d_nParFit%d_%s.png",FunctionType,nPar,fileNameSuffix);
  char c3_fileName[200];
  sprintf(c3_fileName,"fitresiduals_FuncType%d_nParFit%d_%s.png",FunctionType,nPar,fileNameSuffix);
  char c4_fileName_1[200];
  char c4_fileName_2[200];
  sprintf(c4_fileName_1,"fitAndResiduals_FuncType%d_nParFit%d_%s.png",FunctionType,nPar,fileNameSuffix);
  sprintf(c4_fileName_2,"fitAndResiduals_FuncType%d_nParFit%d_%s.svg",FunctionType,nPar,fileNameSuffix);


  ////Canvas0->SaveAs(c0_fileName);
  ////Canvas1->SaveAs(c1_fileName);
  ////Canvas2->SaveAs(c2_fileName);
  ////Canvas3->SaveAs(c3_fileName);
  c->SaveAs(c4_fileName_1);
  c->SaveAs(c4_fileName_2);

}


