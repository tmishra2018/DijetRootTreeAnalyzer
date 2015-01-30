#include <algorithm>

void SimpleDiJetFitV1()
{

  gROOT->Reset();

  // Set TDR Style	
  gROOT->ProcessLine(".L tdrStyle.C");
  gROOT->ProcessLine("setTDRStyle()");
  gStyle->SetOptFit(1111); 

  //########################
  //##### User Options #####
  //########################
  // MC QCD all  
  char input_root_file[500] = "histo_bkg_mjj.root";
  //char input_directory[500] = "";
  
  char fileNameSuffix[500] = "MC_1fb-1"; //i.e. run period

  char input_1Dhistogram[500] = "hist_allCutsQCD";
  double minX_mass = 1100.;
  double maxX_mass = 6000.;

  //Fit functions
  // 0: DEFAULT (4 par.) - "( [0]*TMath::Power(1-x/8000,[1]) ) / ( TMath::Power(x/8000,[2]+[3]*log(x/8000)) )" 
  //
  // 1: VARIATION 1 (5 par.) - "( [0]*TMath::Power(1-x/8000,[1])*(1+[4]*x/8000) ) / ( TMath::Power(x/8000,[2]+[3]*log(x/8000)) )"
  //
  // 2: VARIATION 2 (6 par.) - "( [0]*TMath::Power(1-x/8000,[1])*(1+[4]*x/8000+[5]*pow(x/8000,2)) ) / ( TMath::Power(x/8000,[2]+[3]*log(x/8000)) )"
  //    --> 2nd order poly extension : inspired by HERA PDF 1.0 [http://arxiv.org/abs/arXiv:0911.0884 , Eq. 4.1]
  //
  // 3: VARIATION 3 (7 par.) - "( [0]*TMath::Power(1-x/8000,[1])*exp([4]*x/8000)*TMath::Power(1+exp([5])*x/8000,[6]) ) / ( TMath::Power(x/8000,[2]+[3]*log(x/8000)) )"
  //    --> "exponential" extension wrt to DEFAULT - inspired by CTEQ 2008 [http://arxiv.org/pdf/hep-ph/0201195v3.pdf , Eq. 4]
  //
  // 4: VARIATION 4 (5 par.) - "( [0]*TMath::Power(1-x/8000,[1]) ) / ( TMath::Power(x/8000,[2]+[3]*log(x/8000)+[4]*TMath::Power(log(x/8000),2)) )" 
  //    --> "log" extension wrt to DEFAULT     
  //
  // 5: VARIATION 5 (6 par.) - "( [0]*TMath::Power(1-x/8000,[1]) ) / ( TMath::Power(x/8000,[2]+[3]*log(x/8000)+[4]*TMath::Power(log(x/8000),2)+[5]*TMath::Power(log(x/8000),3)) )" 
  //    --> "log" extension wrt to DEFAULT     
  //
  // 6: VARIATION 6 (7 par.) - "( [0]*TMath::Power(1-x/8000,[1]) ) / ( TMath::Power(x/8000,[2]+[3]*log(x/8000)+[4]*TMath::Power(log(x/8000),2)+[5]*TMath::Power(log(x/8000),3)+[6]*TMath::Power(log(x/8000),4)) )" 
  //    --> "log" extension wrt to DEFAULT     
  //
  const int FunctionType = 0;

  int number_of_variableWidth_bins = 88 - 1;
  Double_t massBins[88] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649,  693, 740, 788, 838, 890, 944, 1000, 1058, 1100, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7000, 7250,7500,7750,8000}; 

  //================================================================================================================

  int nPar=-1; 

  //### input file and 2D histo
  TFile *file0=TFile::Open( input_root_file );
  //TDirectoryFile* DQMData_Merged_Runs_DataScouting_Run_summary_DiJet = (TDirectoryFile*) file0->Get( input_directory );
  //  DQMData_Merged_Runs_DataScouting_Run_summary_DiJet->cd();
  //TH2D* h_DEta_Mass = (TH2D*) DQMData_Merged_Runs_DataScouting_Run_summary_DiJet->Get( input_2Dhistogram );  
  //h_DEta_Mass->Draw("colz");
  TH1D* hist_mass_original = (TH1D*) file0->Get( input_1Dhistogram );  
  hist_mass_original->Sumw2();


  float bin_width_X =  hist_mass_original->GetBinWidth(1);
  cout << "bin_width_X: " << bin_width_X << endl;
  //   cout << "bin_width_Y:" << bin_width_Y << endl;

  //create dataset 
  // #1 : random dataset 
  //1 Gev binning
  TH1F* hist_mass_1GeV = new TH1F("hist_mass_1GeV", "", 8000, 1., 8000. );
  hist_mass_1GeV->Sumw2(); 
  float bincontent;
  float binerror;
  float binwidth;
  //for (int i=1; i<= 8000){
  //  bincontent=RooRandom::randomGenerator()->Poisson(hist_mass_original->GetBinContent(i));
  //  hist_mass_1GeV->SetBinContent(bincontent);
  //}
  hist_mass_1GeV->FillRandom(hist_mass_original, int(hist_mass_original->Integral()));
 
 //giulia test -> WARNING: not using pseudodataset!!! 
  //TH1F* hist_binned= (TH1F*)hist_mass_original->Rebin(number_of_variableWidth_bins,"hist_binned",massBins);
  hist_binned= (TH1F*)hist_mass_1GeV->Rebin(number_of_variableWidth_bins,"hist_binned",massBins); 
  //variable binning random dataset

  TH1F* hist_mass = (TH1F*)hist_binned->Clone("hist_mass");
  hist_mass->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  for (int i=1; i<= number_of_variableWidth_bins; i++){
    bincontent = hist_binned->GetBinContent(i);
    binwidth = hist_binned->GetBinWidth(i);
    binerror = hist_binned->GetBinError(i);
    hist_mass->SetBinContent(i,bincontent/binwidth);   
    hist_mass->SetBinError(i,binerror/binwidth);   
  }
   
  //hist_mass->Draw();


  //### fit mass histogram with background function
  TF1 *M1Bkg;

  // 0: DEFAULT (4 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )" 
  if( FunctionType==0 )    
    {
      nPar=4;
      M1Bkg = new TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )",minX_mass,maxX_mass);
      // M1Bkg->SetParameter(0,0.05);
      // M1Bkg->SetParameter(1,7.1);
      // M1Bkg->SetParameter(2,5.9);
      // M1Bkg->SetParameter(3,0.2);
      M1Bkg->SetParameter(0,0.098);
      M1Bkg->SetParameter(1,1.49);
      M1Bkg->SetParameter(2,7.75);
      M1Bkg->SetParameter(3,0.46);

      M1Bkg->SetParLimits(0,0,0.3);
      M1Bkg->SetParLimits(1,0,6.5);
      M1Bkg->SetParLimits(2,6.,9.);
      M1Bkg->SetParLimits(3,0,0.8);
    }

  // 1: VARIATION 1 (5 par.) - "( [0]*TMath::Power(1-x/13000,[1])*(1+[4]*x/13000) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
  if( FunctionType==1 )    
    {
      nPar=5;
      M1Bkg = new TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1])*(1+[4]*x/13000) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )",minX_mass,maxX_mass);
      M1Bkg->SetParameter(0,0.005);
      M1Bkg->SetParameter(1,9.3);
      M1Bkg->SetParameter(2,7.2);
      M1Bkg->SetParameter(3,0.4);
      M1Bkg->SetParameter(4,3.1);
      //      M1Bkg->SetParLimits(4,-5,5);
    }

  // 2: VARIATION 2 (6 par.) - "( [0]*TMath::Power(1-x/13000,[1])*(1+[4]*x/13000+[5]*pow(x/13000,2)) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
  //    --> 2nd order poly extension : inspired by HERA PDF 1.0 [http://arxiv.org/abs/arXiv:0911.0884 , Eq. 4.1]
  if( FunctionType==2 )    
    {
      nPar=6;
      M1Bkg = new TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1])*(1+[4]*x/13000+[5]*pow(x/13000,2)) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )",minX_mass,maxX_mass);
      M1Bkg->SetParameter(0,0.005);
      M1Bkg->SetParameter(1,9.3);
      M1Bkg->SetParameter(2,7.2);
      M1Bkg->SetParameter(3,0.4);
      M1Bkg->SetParameter(4,3.1);
      M1Bkg->SetParameter(5,25.6);
      //      M1Bkg->SetParLimits(5,10,50);      
    }

  // 3: VARIATION 3 (7 par.) - "( [0]*TMath::Power(1-x/13000,[1])*exp([4]*x/13000)*TMath::Power(1+exp([5])*x/13000,[6]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
  //    --> "exponential" extension wrt to DEFAULT - inspired by CTEQ 2008 [http://arxiv.org/pdf/hep-ph/0201195v3.pdf , Eq. 4]
  if( FunctionType==3 )    
    {
      nPar=7;
      M1Bkg = new TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1])*exp([4]*x/13000)*TMath::Power(1+exp([5])*x/13000,[6]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )" ,minX_mass,maxX_mass);
      M1Bkg->SetParameter(0,0.005);
      M1Bkg->SetParameter(1,15.1);
      M1Bkg->SetParameter(2,7.2);
      M1Bkg->SetParameter(3,0.4);
      M1Bkg->SetParameter(4,13.0);
      M1Bkg->SetParameter(5,-4.0);
      M1Bkg->SetParameter(6,70.0);
      //      M1Bkg->SetParLimits(4,-1,1);
    }

  // 4: VARIATION 4 (5 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)) )" 
  //    --> "log" extension wrt to DEFAULT     
  if( FunctionType==4 )    
    {
      nPar=5;
      M1Bkg = new TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)) )",minX_mass,maxX_mass);
      M1Bkg->SetParameter(0,0.05);
      M1Bkg->SetParameter(1,7.1);
      M1Bkg->SetParameter(2,5.9);
      M1Bkg->SetParameter(3,0.2);
      M1Bkg->SetParameter(4,0.);
      //M1Bkg->SetParLimits(3,0,0.4);
    }

  // 5: VARIATION 5 (6 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)+[5]*TMath::Power(log(x/13000),3)) )" 
  //    --> "log" extension wrt to DEFAULT     
  if( FunctionType==5 )    
    {
      nPar=6;
      M1Bkg = new TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)+[5]*TMath::Power(log(x/13000),3)) )",minX_mass,maxX_mass);
      M1Bkg->SetParameter(0,0.05);
      M1Bkg->SetParameter(1,7.1);
      M1Bkg->SetParameter(2,5.9);
      M1Bkg->SetParameter(3,0.2);
      M1Bkg->SetParameter(4,0.);
      M1Bkg->SetParameter(5,0.);
      //M1Bkg->SetParLimits(3,0,0.4);
    }

  // 6: VARIATION 6 (7 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)+[5]*TMath::Power(log(x/13000),3)+[6]*TMath::Power(log(x/13000),4)) )" 
  //    --> "log" extension wrt to DEFAULT     
  if( FunctionType==6 )    
    {
      nPar=7;
      M1Bkg = new TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)+[5]*TMath::Power(log(x/13000),3)+[6]*TMath::Power(log(x/13000),4)) )",minX_mass,maxX_mass);
      M1Bkg->SetParameter(0,0.05);
      M1Bkg->SetParameter(1,7.1);
      M1Bkg->SetParameter(2,5.9);
      M1Bkg->SetParameter(3,0.2);
      M1Bkg->SetParameter(4,0.);
      M1Bkg->SetParameter(5,0.);
      M1Bkg->SetParameter(6,0.);
      //M1Bkg->SetParLimits(3,0,0.4);
    }

//  //////////////////////////////////////////////////////////////////////////////////
////  ______ _ _     _     _             _               __    _____       _   _ 
////  |  ___(_) |   | |   (_)           (_)             /  |  |  __ \     | | | |
////  | |_   _| |_  | |__  _ _ __  _ __  _ _ __   __ _  `| |  | |  \/ ___ | | | |
////  |  _| | | __| | '_ \| | '_ \| '_ \| | '_ \ / _` |  | |  | | __ / _ \| | | |
////  | |   | | |_  | |_) | | | | | | | | | | | | (_| | _| |_ | |_\ \  __/\ \_/ /
////  \_|   |_|\__| |_.__/|_|_| |_|_| |_|_|_| |_|\__, | \___/  \____/\___| \___/ 
////                                              __/ |                                               
////                          		         |___/                           
////
// //////////////////////////////////////////////////////////////////////////////////////
// 
//  TFitResultPtr r;
//  int stopProgram=1;
//  for( int loop=0; loop<10; loop++)
//    {
//
//      r = hist_mass_1GeV->Fit("M1Bkg","ILSR","",minX_mass,maxX_mass);      
//      Int_t fitStatus = r;
//      if(fitStatus==0)
//	{
//   	  stopProgram=0;
//	  cout << "###################################" << endl;
//	  cout << "Result fit with binning of 1 GeV : " << endl;
//     	  r->Print("V");  
//     	  break;
//     	}
//    }
//
//  if(stopProgram==1)
//    {
//      cout << "######################" << endl;
//      cout << "######################" << endl;
//      cout << "ERROR : Fit with  1 GeV binning failed!!!!" << endl;
//      cout << "######################" << endl;
//      cout << "######################" << endl;
//      break;
//    }
//
//
//  //### calculate chi2 by hand for 1 GeV bin fit
//  const double alpha = 1 - 0.6827;
//  int NumberOfObservations_FixBin = 0;
//  double chi2_FixBin = 0.;
//  for(int bin=1; bin<hist_mass_1GeV->GetNbinsX()+1; bin++)
//    {
//      if( hist_mass_1GeV->GetXaxis()->GetBinLowEdge(bin)>=minX_mass 
// 	  && hist_mass_1GeV->GetXaxis()->GetBinUpEdge(bin)<=maxX_mass )
// 	{
//	  NumberOfObservations_FixBin++;
//
//	  double data = hist_mass_1GeV->GetBinContent(bin);
//	  double err_data = hist_mass_1GeV->GetBinError(bin);
//	  if( data == 0 )
//	    {
//	      err_data = 1.8;
//	    }
////  	  if(data < 1)
//// 	    {
//// 	      double data_L =  (data==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,data,1.)) ;
//// 	      double data_U =  ROOT::Math::gamma_quantile_c(alpha/2,data+1,1) ;
//// 	      double err_data_l = data - data_L;
//// 	      double err_data_h = data_U - data;
//// 	      err_data = max(err_data_h,err_data_l);	      
//// 	    }
//// 	  else
//// 	    {
//// 	      err_data = sqrt(data);
//// 	    }
//	  //double fit = M1Bkg->Eval(hist_mass->GetBinCenter(bin));
//	  double fit = M1Bkg->Integral(hist_mass_1GeV->GetXaxis()->GetBinLowEdge(bin) 
//				       , hist_mass_1GeV->GetXaxis()->GetBinUpEdge(bin) ); 
//	  fit = fit / ( hist_mass_1GeV->GetBinWidth(bin) ); 
//	  chi2_FixBin += pow( (data - fit) , 2 ) / pow( err_data , 2 );
//	}
//    }
//  int ndf_FixBin = NumberOfObservations_FixBin - nPar  ;
//  cout << "============================" << endl;
//  cout << "NumberOfObservations_FixBin: " << NumberOfObservations_FixBin << endl;
//  cout << "ndf_FixBin: " << ndf_FixBin << endl;
//  cout << "chi2_FixBin: " << chi2_FixBin << endl;
//  cout << "============================" << endl;  
//
											                                                                                         
 
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

  // fit residuals and chi2
  TH1D* hist_fit_residual_vsMass = new TH1D("hist_fit_residual_vsMass","hist_fit_residual_vsMass",number_of_variableWidth_bins,massBins);
  TH1D* hist_fit_residual = new TH1D("hist_fit_residual","hist_fit_residual",10,-5,5);
  int NumberOfObservations_VarBin = 0;
  double chi2_VarBin = 0.;

  //cout << hist_mass_varbin->GetNbinsX() << endl;  
  for(int bin=1; bin<number_of_variableWidth_bins; bin++)
    {
      if( hist_mass->GetXaxis()->GetBinLowEdge(bin)>=minX_mass 
 	  && hist_mass->GetXaxis()->GetBinUpEdge(bin)<=maxX_mass )
 	{	  
	  NumberOfObservations_VarBin++;

 	  //cout << hist_mass_varbin->GetXaxis()->GetBinLowEdge(bin) << endl;
	  double data = hist_mass->GetBinContent(bin);
	  double err_data = hist_mass->GetBinError(bin);
	  if( data == 0 )
	    {
	      err_data = 1.8 / hist_mass->GetBinWidth(bin) ; //giulia : perche`??
	    }
	  //double fit = M1Bkg->Eval(hist_mass->GetBinCenter(bin));
 	  double fit = M1Bkg->Integral(hist_mass->GetXaxis()->GetBinLowEdge(bin) 
 	  			       , hist_mass->GetXaxis()->GetBinUpEdge(bin) ); 
	  fit = fit / ( hist_mass->GetBinWidth(bin) );
	  
	  double err_tot = err_data;	  
	  double fit_residual = (data - fit) / err_tot;
	  double err_fit_residual = 1;
	    	  
	  chi2_VarBin += pow( (data - fit) , 2 ) / pow( err_data , 2 );	 

	  // 	  cout << "data, err_data, fit: " << data << ", " << err_data << ", " << fit << endl;
	  // 	  cout << "bin, fit residual : " << bin << ", " <<fit_residual << endl;	  
	  hist_fit_residual_vsMass->SetBinContent(bin,fit_residual);
	  hist_fit_residual_vsMass->SetBinError(bin,err_fit_residual);
	  hist_fit_residual->Fill(fit_residual);
 	}
    }
  int ndf_VarBin = NumberOfObservations_VarBin - nPar  ;
  cout << "============================" << endl;
  cout << "NumberOfObservations_VarBin: " << NumberOfObservations_VarBin << endl;
  cout << "ndf_VarBin: " << ndf_VarBin << endl;
  cout << "chi2_VarBin: " << chi2_VarBin << endl;
  cout << "============================" << endl;  

  //### Draw plots

  TCanvas *Canvas0 = new TCanvas("Canvas0","Canvas0",11,51,700,500);
  Canvas0->cd();
  Canvas0->SetLogy();
  hist_mass_1GeV->GetYaxis()->SetTitle("Events");
  hist_mass_1GeV->Draw();  
  M1Bkg->SetLineColor(2);
  M1Bkg->Draw("same");     

  TCanvas *Canvas1 = new TCanvas("Canvas1","Canvas1",11,51,700,500);
  Canvas1->cd();
  Canvas1->SetLogy();
  hist_mass->GetYaxis()->SetTitle("Events / bin width");
  hist_mass->Draw();  
  M1Bkg->SetLineColor(2);
  M1Bkg->Draw("same");     

  TCanvas *Canvas2 = new TCanvas("Canvas2","Canvas2",11,51,700,500);
  Canvas2->cd();
  Canvas2->SetGridx();
  Canvas2->SetGridy();
  Canvas2->SetLogx();
  hist_fit_residual_vsMass->GetYaxis()->SetLimits(-5,5);
  hist_fit_residual_vsMass->GetYaxis()->SetRangeUser(-5,5);
  hist_fit_residual_vsMass->GetYaxis()->SetTitle("(data - fit) / #sqrt{data}");
  //hist_fit_residual_vsMass->GetXaxis()->SetRangeUser(minX_mass,maxX_mass);
  hist_fit_residual_vsMass->GetXaxis()->SetRangeUser(1000,10000);
  hist_fit_residual_vsMass->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  hist_fit_residual_vsMass->Draw();

  TCanvas *Canvas3 = new TCanvas("Canvas3","Canvas3",11,51,700,500);
  Canvas3->cd();
  hist_fit_residual->GetXaxis()->SetTitle("(data - fit) / #sqrt{data}");
  hist_fit_residual->GetYaxis()->SetTitle("Number of bins");
  hist_fit_residual->GetYaxis()->SetRangeUser(0,number_of_variableWidth_bins/3);
  hist_fit_residual->Draw();
  hist_fit_residual->Fit("gaus","L","",-3,3);


//#################### canvas "note" style ####################
  
  TCanvas* c = new TCanvas("c","DijetMass cross section with Fit and QCD MC",600,650);
  c->GetWindowHeight();
  c->GetWindowWidth();
  c->SetLogy();
  c->Divide(1,2,0,0,0);
  c->cd(1);
  p11_1 = (TPad*)c->GetPad(1);
  p11_1->SetPad(0.01,0.23,0.99,0.98);
  p11_1->SetLogy();
  p11_1->SetRightMargin(0.05);
  p11_1->SetTopMargin(0.05);

  // Pave text
  TPaveText *pave_fit = new TPaveText(0.18,0.15,0.40,0.27,"NDC");
  pave_fit->AddText(" #sqrt{s} = 13 TeV");
  pave_fit->AddText("|#eta| < 2.5, |#Delta#eta| < 1.3");
  pave_fit->AddText("M_{jj} > 1100 GeV");
  pave_fit->AddText("Wide Jets");
  pave_fit->SetFillColor(0);
  pave_fit->SetLineColor(0);
  pave_fit->SetFillStyle(0);
  pave_fit->SetBorderSize(0);
  pave_fit->SetTextFont(42);
  pave_fit->SetTextSize(0.03);
  pave_fit->SetTextAlign(12); 


  TPaveText *pt1 = new TPaveText(0.1284756,0.9602144,0.3887139,0.9902251,"brNDC");
  pt1->SetBorderSize(0);
  pt1->SetFillColor(0);
  pt1->SetFillStyle(0);
  pt1->SetLineColor(0);
  pt1->SetTextAlign(12);
  pt1->SetTextSize(0.035);
  TText *text = pt1->AddText("CMS");

  TPaveText *pt2 = new TPaveText(0.45,0.96,0.65,0.99,"brNDC");
  pt2->SetBorderSize(0);
  pt2->SetFillColor(0);
  pt2->SetFillStyle(0);
  pt2->SetLineColor(0);
  pt2->SetTextAlign(12);
  pt2->SetTextSize(0.035);
  TText *text2 = pt2->AddText("#sqrt{s} = 13 TeV");

  TPaveText *pt3 = new TPaveText(0.7687988,0.9602144,0.9297357,0.9902251,"brNDC");
  pt3->SetBorderSize(0);
  pt3->SetFillColor(0);
  pt3->SetFillStyle(0);
  pt3->SetLineColor(0);
  pt3->SetTextAlign(12);
  pt3->SetTextSize(0.035);
  TText *text3 = pt3->AddText("L= 1 fb^{-1}");

  TH1F *vFrame = p11_1->DrawFrame(1000.0,0.000000001,8000.0,10.0);

  vFrame->SetTitle("");
  vFrame->SetXTitle("Dijet Mass (GeV)");
  vFrame->GetXaxis()->SetTitleSize(0.06);
  vFrame->SetYTitle("events / bin width");
  
  vFrame->GetYaxis()->SetTitleSize(0.12);
  vFrame->GetYaxis()->SetLabelSize(0.07);
  vFrame->GetYaxis()->SetTitleOffset(0.50);
  vFrame->GetXaxis()->SetTitleOffset(0.90);
  vFrame->GetXaxis()->SetTitleSize(0.18);
  vFrame->GetXaxis()->SetLabelSize(0.1);

  

  hist_mass->GetXaxis()->SetRangeUser(1000,8000.0);
  hist_mass->SetTitle("");
  hist_mass->SetLineColor(1);
  hist_mass->SetFillColor(1);
  hist_mass->SetMarkerColor(1);
  hist_mass->SetMarkerStyle(20);
  hist_mass->Draw("P");
  M1Bkg->SetLineWidth(2);
  M1Bkg->SetLineStyle(2);
  M1Bkg->SetLineColor(2);
  M1Bkg->Draw("same");

  TLegend *leg = new TLegend(0.5564991,0.4,0.8903575,0.575812);
  leg->SetTextSize(0.03146853);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetMargin(0.35);
  leg->AddEntry(hist_mass,"pseudo data" ,"PL");
  leg->AddEntry(M1Bkg,"fit to pseudodata","L");
  leg->Draw("same");

  pt1->Draw("same"); 
  pt2->Draw("same");
  pt3->Draw("same");
  pave_fit->Draw("same");

  //redraw axis
  //p11_1->RedrawAxis();
  //p11_1->Update();
  //cout << "MIN: " << p11_1->GetUxmin() << endl;
  //cout << "MAX: " << p11_1->GetUxmax() << endl;

  //---- Next PAD

  c->cd(2);
  p11_2 = (TPad*)c->GetPad(2);
  p11_2->SetPad(0.01,0.02,0.99,0.24);
  p11_2->SetBottomMargin(0.35);
  p11_2->SetRightMargin(0.05);
  p11_2->SetGridx();
  p11_2->SetGridy();
  //c3_2->SetTickx(50);

  TH1F *vFrame2 = p11_2->DrawFrame(p11_1->GetUxmin(), -4.5, p11_1->GetUxmax(), 4.5);

  vFrame2->SetTitle("");
  vFrame2->SetXTitle("Dijet Mass (GeV)");
  vFrame2->GetXaxis()->SetTitleSize(0.06);
  vFrame2->SetYTitle("(pseudoData-Fit)/#sigma");
  vFrame2->GetYaxis()->SetTitleSize(0.12);
  vFrame2->GetYaxis()->SetLabelSize(0.07);
  vFrame2->GetYaxis()->SetTitleOffset(0.50);
  vFrame2->GetXaxis()->SetTitleOffset(0.90);
  vFrame2->GetXaxis()->SetTitleSize(0.18);
  vFrame2->GetXaxis()->SetLabelSize(0.1);

  hist_fit_residual_vsMass->GetXaxis()->SetRangeUser(1000.,8000.);
  hist_fit_residual_vsMass->GetYaxis()->SetRangeUser(-4.,4.);
  hist_fit_residual_vsMass->SetLineWidth(0);
  hist_fit_residual_vsMass->SetFillColor(2);
  hist_fit_residual_vsMass->SetLineColor(1);
  hist_fit_residual_vsMass->Draw("SAMEHIST");

  TLine *line = new TLine(1000.,0,8000,0);
  line->Draw("");
  //c->Close();

  //### Output files
  char output_root_file[500];
  sprintf(output_root_file,"dijetFitResults_FuncType%d_nParFit%d_%s.root",FunctionType,nPar,fileNameSuffix); 

  TFile f_output(output_root_file,"RECREATE");
  f_output.cd();
  Canvas0->Write();
  Canvas1->Write();
  Canvas2->Write();
  Canvas3->Write();
  hist_mass_1GeV->Write();
  hist_binned->Write();
  hist_mass->Write();
  c->Write();
  r_bin->Write();
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
  char c4_fileName[200];
  sprintf(c4_fileName,"fitAndResiduals_FuncType%d_nParFit%d_%s.png",FunctionType,nPar,fileNameSuffix);

//  char c0_fileName_vec[200];
//  sprintf(c0_fileName,"dijetmass_FuncType%d_nParFit%d_%s.svg",FunctionType,nPar,fileNameSuffix);
//  char c1_fileName_vec[200];
//  sprintf(c1_fileName,"dijetmass_varbin_FuncType%d_nParFit%d_%s.svg",FunctionType,nPar,fileNameSuffix);
//  char c2_fileName_vec[200];
//  sprintf(c2_fileName,"fitresiduals_vs_mass_FuncType%d_nParFit%d_%s.svg",FunctionType,nPar,fileNameSuffix);
//  char c3_fileName_vec[200];
//  sprintf(c3_fileName,"fitresiduals_FuncType%d_nParFit%d_%s.svg",FunctionType,nPar,fileNameSuffix);
//  char c4_filename[200];
//  sprintf(c4_fileName,"fitAndResiduals_FuncType%d_nParFit%d_%s.svg",FunctionType,nPar,fileNameSuffix);
//
//


  Canvas0->SaveAs(c0_fileName);
  Canvas1->SaveAs(c1_fileName);
  Canvas2->SaveAs(c2_fileName);
  Canvas3->SaveAs(c3_fileName);
  c->SaveAs(c4_fileName);

}


  //########################################################################################
  //### OLD FILES WITH BUG (to reproduce the plots presented at the first dijet meeting) ###
  //########################################################################################
  //   // Run2012B - only dijet, razoe, alfaT
  //   char input_root_file[500] = "root://eoscms//eos/cms/store/cmst3/user/santanas/DataScouting/DQM_histograms/DataScouting_V00-01-03_Run2012B_runrange_193752-197044_dijet_alfaT_razor.root";
  //   // Run2012B - all analyses (and 1-2% more events)
  //   //char input_root_file[500] = "root://eoscms//eos/cms/store/cmst3/user/santanas/DataScouting/DQM_histograms/DataScouting_V00-01-05_Run2012B_runrange_193752-197044_dijet_alfaT_razor_dijetpairs_trijetpairs.root";
  //   // Run2012C - all analyses
  //   //char input_root_file[500] = "root://eoscms//eos/cms/store/cmst3/user/santanas/DataScouting/DQM_histograms/DataScouting_V00-01-05_Run2012C_runrange_197885-203755_dijet_alfaT_razor_dijetpairs_trijetpairs.root";
  //   // Run2012B+Run2012C - all analyses
  //   //char input_root_file[500] = "root://eoscms//eos/cms/store/cmst3/user/santanas/DataScouting/DQM_histograms/DataScouting_V00-01-05_Run2012B_Run2012C_runrange_193752-203755_dijet_alfaT_razor_dijetpairs_trijetpairs.root";
  //   char input_directory[500] = "DQMData_Merged Runs_DataScouting_Run summary_DiJet;1";

