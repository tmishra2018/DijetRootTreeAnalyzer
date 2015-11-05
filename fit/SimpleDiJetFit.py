#!/usr/bin/python

from ROOT import *
from array import array
import CMS_lumi, tdrstyle
import subprocess
import os
import imp
import multiprocessing
from itertools import repeat
import math
import optparse

gStyle.SetOptFit(1111) 
#set the tdr style
tdrstyle.setTDRStyle()

lumi=1769
#lumi=809
#lumi=344
#lumi=974.
#lumi = 803.
#lumi = 65.0
#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_7TeV = "4.8 fb^{-1}"
CMS_lumi.lumi_8TeV = "18.3 fb^{-1}"
#CMS_lumi.lumi_13TeV = "974 pb^{-1}"
CMS_lumi.lumi_13TeV = str(int(lumi))+" pb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = ""
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4
 
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("--plotSig",action="store_true",dest="plotSig",default=False)
  
(options, args) = parser.parse_args() 
plotSig = options.plotSig


#Fit function
# 0: DEFAULT (4 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )" 
#
#1: VARIATION 1 (5 par.) - "( [0]*TMath::Power(1-x/13000,[1])*(1+[4]*x/13000) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
#
#2: VARIATION 2 (6 par.) - "( [0]*TMath::Power(1-x/13000,[1])*(1+[4]*x/13000+[5]*pow(x/13000,2)) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
#   --> 2nd order poly extension : inspired by HERA PDF 1.0 [http://arxiv.org/abs/arXiv:0911.0884 , Eq. 4.1]
#
#3: VARIATION 3 (7 par.) - "( [0]*TMath::Power(1-x/13000,[1])*exp([4]*x/13000)*TMath::Power(1+exp([5])*x/13000,[6]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
#   --> "exponential" extension wrt to DEFAULT - inspired by CTEQ 2008 [http://arxiv.org/pdf/hep-ph/0201195v3.pdf , Eq. 4]
#
#4: VARIATION 4 (5 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)) )" 
#   --> "log" extension wrt to DEFAULT     
#
#5: VARIATION 5 (6 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)+[5]*TMath::Power(log(x/13000),3)) )" 
#   --> "log" extension wrt to DEFAULT     
#
#6: VARIATION 6 (7 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)+[5]*TMath::Power(log(x/13000),3)+[6]*TMath::Power(log(x/13000),4)) )" 
#   --> "log" extension wrt to DEFAULT     
#
number_of_variableWidth_bins = 103

massBins =[1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430,10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000];

v_massBins = array("d",massBins)
#sf = 1.07
sf=1.
sigmaQstar4500 = 0.2827E-01
sigmaQstar3000 = 0.8116E+00
sigmaAxigluon3000 =  0.6610
sigmaString5300 = 0.11940
sigmaS82200 = 3.27
sigmaAxigluon5100 = 0.007539 
sigmaString6700 = 0.006634 
sigmaS83200 = 0.137 

label_Qstar3000 = "Excited quark (3 TeV)"
label_Qstar4500 = "Excited quark (4.5 TeV)"
label_String5300 = "String (5.3 TeV)"
label_S82200 = "S8 (2.2 TeV)"
label_Axigluon3000 = "Axigluon (3 TeV)"
label_String6700 = "String (6.7 TeV)"
label_S83200 = "S8 (3.2 TeV)"
label_Axigluon5100 = "Axigluon (5.1 TeV)"

#minX_mass = 2546.
minX_mass = 1181.
#maxX_mass = 3019.
#maxX_mass = 5253.
#maxX_mass = 6564.
maxX_mass = 7589
#FunctionType = -1
FunctionType = 0
#fileNameSuffix = "test_range_1118_3704"
#fileNameSuffix = "test"
#fileNameSuffix = "JEC_L2L3Residuals"
#fileNameSuffix_Qstar4500 = "Qstar4500"
fileNameSuffix_Qstar4500 = "Qstar4500_noSignif"
fileNameSuffix_Qstar3000 = "Qstar3000"
fileNameSuffix_String5300 = "String5300"
fileNameSuffix_S8_2200 = "S8_2200"
fileNameSuffix_Axigluon3000 = "Axigluon3000"
#fileNameSuffix = "RunD_4param_JEC_Summer15_25nsV5_golden_547pb-1_fullPeriod"
#fileNameSuffix = "RunD_4param_JEC_Summer15_25nsV5_golden_period1"
#fileNameSuffix = "RunD_4param_JEC_Summer15_25nsV5_golden_period2"
#fileNameSuffix = "RunD_4param_JEC_Summer15_25nsV5_golden_period3"
#fileNameSuffix = "RunD_4param_JEC_Summer15_25nsV5_golden_periodA"
#fileNameSuffix = "RunD_4param_JEC_Summer15_25nsV5_golden_periodB"
#fileNameSuffix = "RunD_4param_JEC_Summer15_25nsV5_golden_periodC"
#fileNameSuffix = "RunD_4param_JEC_Summer15_25nsV5_DCSonly_974pb-1_fullPeriod"
#fileNameSuffix = "RunD_4param_JEC_Summer15_25nsV5_DCSonly_period1"
#fileNameSuffix = "RunD_4param_JEC_Summer15_25nsV5_DCSonly_period2"
#fileNameSuffix = "RunD_4param_JEC_Summer15_25nsV5_DCSonly_period3"
#fileNameSuffix = "RunD_4param_JEC_Summer15_25nsV5_DCSonly_974pb-1_fitFrom2500"
#fileNameSuffix = "RunD_4param_JEC_Summer15_25nsV5_DCSonly_974pb-1_fitTo3000"
#fileNameSuffix = "RunD_4param_JEC_Summer15_25nsV5_DCSonly_974pb-1"
#fileNameSuffix = "RunD_4param_JEC_Summer15_25nsV5_DCSonly_974pb-1_AK4"
#fileNameSuffix = "RunD_4param_JEC_Summer15_25nsV5_golden_809pb-1"
#fileNameSuffix = "RunD_4param_JEC_Summer15_25nsV5_DCSonly_1769pb-1"
fileNameSuffix = "RunD_4param_JEC_Summer15_25nsV5_DCSonly_1769pb-1_plotSig"

####### INPUT #############
# data 
#input_root_file = "../scripts/plots_data4T_Run2015B_plus_Run2015C_50ns_Spring15_JEC_Summer15_50ns_V4_withSF//histo_data_mjj_fromTree.root"
#input_root_file = "../scripts/Run2015B_plus_Run2015C_50ns_Cert_json_29Aug2015_xsecSpring15_fixedJEC_withSF//histo_data_mjj_fromTree.root"
#input_root_file="../scripts/plots_data4T_Run2015D_DCSonly_390pb-1_JEC_Summer15_25nsV3_withSF/histo_data_mjj_fromTree.root"
#input_root_file="../scripts/plots_data_Run2015D_DCSonly_MC_Spring15_25ns_Summer15_25nsV5_withSF/histo_data_mjj_fromTree.root"
#input_root_file="../scripts/plots_data4T_Run2015D_DCSonly_JEC_Summer15_50nsV5_withSF/histo_data_mjj_fromTree.root"
#input_root_file="../scripts/plots_data4T_Run2015D_golden_JEC_Summer15_25nsV5_205pb-1/histo_data_mjj_fromTree.root"
#input_root_file="../scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF/histo_data_mjj_fromTree.root"
#input_root_file="../scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_period1/histo_data_mjj_fromTree.root"
#input_root_file="../scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_period2/histo_data_mjj_fromTree.root"
#input_root_file="../scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_period3/histo_data_mjj_fromTree.root"
#input_root_file="../scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_periodA/histo_data_mjj_fromTree.root"
#input_root_file="../scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_periodB/histo_data_mjj_fromTree.root"
#input_root_file="../scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_periodC/histo_data_mjj_fromTree.root"
#input_root_file="../scripts/plots_data4T_Run2015D_DCSonly_974pb-1_JEC_Summer15_25nsV5_withSF//histo_data_mjj_fromTree.root"
#input_root_file="../scripts/plots_data4T_Run2015D_DCSonly_974pb-1_JEC_Summer15_25nsV5_withSF_period1/histo_data_mjj_fromTree.root"
#input_root_file="../scripts/plots_data4T_Run2015D_DCSonly_974pb-1_JEC_Summer15_25nsV5_withSF_period2/histo_data_mjj_fromTree.root"
#input_root_file="../scripts/plots_data4T_Run2015D_DCSonly_974pb-1_JEC_Summer15_25nsV5_withSF_period3/histo_data_mjj_fromTree.root"
#input_root_file="../scripts/plots_data4T_Run2015D_DCSonly_974pb-1_JEC_Summer15_25nsV5_withSF//histo_data_Dijet_MassAK4_fromTree.root"
#input_root_file="../scripts/plots_data4T_Run2015D_golden_809pb-1_JEC_Summer15_25nsV5_withSF/histo_data_mjj_fromTree.root"
input_root_file="../scripts/plots_data4T_Run2015D_DCSonly_1769pb-1_JEC_Summer15_25nsV5_withSF/histo_data_mjj_fromTree.root"
#input_root_file="../scripts/plots_data4T_Run2015D_DCSonly_JEC_Summer15_25nsV5_withSF_sample45678_dEta_less_2.6//histo_data_mjj_fromTree.root"

###mc
#input_root_file_mc = "../scripts/plots_data4T_Run2015B_plus_Run2015C_50ns_Spring15_JEC_Summer15_50ns_V4_withSF/histo_data_mjj_fromTree.root"
#input_root_file_mc = "../scripts/Run2015B_plus_Run2015C_50ns_Cert_json_29Aug2015_xsecSpring15_fixedJEC_withSF//histo_data_mjj_fromTree.root"
#input_root_file_mc="../scripts/plots_data4T_Run2015D_DCSonly_390pb-1_JEC_Summer15_25nsV3_withSF/histo_data_mjj_fromTree.root"
#input_root_file_mc="../scripts/plots_data_Run2015D_DCSonly_MC_Spring15_25ns_Summer15_25nsV5_withSF/histo_data_mjj_fromTree.root"
#input_root_file_mc="../scripts/plots_data4T_Run2015D_DCSonly_JEC_Summer15_50nsV5_withSF/histo_data_mjj_fromTree.root"
#input_root_file_mc="../scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_period1/histo_data_mjj_fromTree.root"
#input_root_file_mc="../scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_period2/histo_data_mjj_fromTree.root"
#input_root_file_mc="../scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_period3/histo_data_mjj_fromTree.root"
#input_root_file_mc="../scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_periodA/histo_data_mjj_fromTree.root"
#input_root_file_mc="../scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_periodB/histo_data_mjj_fromTree.root"
#input_root_file_mc="../scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_periodC/histo_data_mjj_fromTree.root"
#input_root_file_mi="../scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF/histo_data_mjj_fromTree.root"
#input_root_file_mc="../scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF/histo_data_mjj_fromTree.root"
#input_root_file_mc="../scripts/plots_data4T_Run2015D_DCSonly_974pb-1_JEC_Summer15_25nsV5_withSF//histo_data_mjj_fromTree.root"
#input_root_file_mc="../scripts/plots_data4T_Run2015D_DCSonly_974pb-1_JEC_Summer15_25nsV5_withSF_period1/histo_data_mjj_fromTree.root"
#input_root_file_mc="../scripts/plots_data4T_Run2015D_DCSonly_974pb-1_JEC_Summer15_25nsV5_withSF_period2/histo_data_mjj_fromTree.root"
#input_root_file_mc="../scripts/plots_data4T_Run2015D_DCSonly_974pb-1_JEC_Summer15_25nsV5_withSF_period3/histo_data_mjj_fromTree.root"
#input_root_file_mc="../scripts/plots_data4T_Run2015D_DCSonly_974pb-1_JEC_Summer15_25nsV5_withSF//histo_data_Dijet_MassAK4_fromTree.root"
#input_root_file_mc="../scripts/plots_data4T_Run2015D_golden_809pb-1_JEC_Summer15_25nsV5_withSF/histo_data_mjj_fromTree.root"
input_root_file_mc="../scripts/plots_data4T_Run2015D_DCSonly_1769pb-1_JEC_Summer15_25nsV5_withSF/histo_data_mjj_fromTree.root"
#input_root_filei_mc="../scripts/plots_data4T_Run2015D_DCSonly_JEC_Summer15_25nsV5_withSF_sample45678_dEta_less_2.6//histo_data_mjj_fromTree.root"

### input file and 1D histo
input_root_file_signal_qq = "ResonanceShapes_qq_13TeV_PU30_Spring15.root"
input_root_file_signal_qg = "ResonanceShapes_qg_13TeV_PU30_Spring15.root"
input_root_file_signal_gg = "ResonanceShapes_gg_13TeV_PU30_Spring15.root"

file0 = TFile.Open( input_root_file )
fileMC = TFile.Open( input_root_file_mc )
file_sig_qq = TFile.Open(input_root_file_signal_qq)
file_sig_qg = TFile.Open(input_root_file_signal_qg)
file_sig_gg = TFile.Open(input_root_file_signal_gg)
input_1Dhistogram = "h_dat"
input_1Dhistogram_mc = "hist_allCutsQCD"
input_sig_Qstar4500 = "h_qg_4500" 
input_sig_Qstar3000 = "h_qg_3000" 
input_sig_Axigluon3000 = "h_qq_3000" 
input_sig_String5300 = "h_qg_5300" 
input_sig_S82200 = "h_gg_2200" 
input_sig_Axigluon5100 = "h_qq_5100" 
input_sig_String6700 = "h_qg_6700" 
input_sig_S83200 = "h_gg_3200" 

hist_mass_original = file0.Get(input_1Dhistogram)
hist_binned = hist_mass_original.Rebin(number_of_variableWidth_bins,"hist_binned",v_massBins)
hist_mass = TH1F("hist_mass","",number_of_variableWidth_bins,v_massBins)
#hist_mass_original.Scale(1/lumi)

hist_mass_original_mc = fileMC.Get(input_1Dhistogram_mc)
hist_mass_original_mc.Scale(sf)
hist_binned_mc = hist_mass_original_mc.Rebin(number_of_variableWidth_bins,"hist_binned_MC",v_massBins)
hist_mass_mc = TH1F("hist_mass_mc","",number_of_variableWidth_bins,v_massBins)

h_sig_Qstar4500 = file_sig_qg.Get(input_sig_Qstar4500)
integral_sig_Qstar4500 = h_sig_Qstar4500.Integral()
h_sig_Qstar4500.Scale( sigmaQstar4500 / integral_sig_Qstar4500)

h_sig_Qstar3000 = file_sig_qg.Get(input_sig_Qstar3000)
integral_sig_Qstar3000 = h_sig_Qstar3000.Integral()
h_sig_Qstar3000.Scale( sigmaQstar3000 / integral_sig_Qstar3000)

h_sig_String5300 = file_sig_qg.Get(input_sig_String5300)
integral_sig_String5300 = h_sig_String5300.Integral()
h_sig_String5300.Scale( sigmaString5300 / integral_sig_String5300)

h_sig_S82200 = file_sig_gg.Get(input_sig_S82200)
integral_sig_S82200 = h_sig_S82200.Integral()
h_sig_S82200.Scale( sigmaS82200 / integral_sig_S82200)

h_sig_Axigluon3000 = file_sig_qq.Get(input_sig_Axigluon3000)
integral_sig_Axigluon3000 = h_sig_Axigluon3000.Integral()
h_sig_Axigluon3000.Scale( sigmaAxigluon3000 / integral_sig_Axigluon3000)

h_sig_String6700 = file_sig_qg.Get(input_sig_String6700)
integral_sig_String6700 = h_sig_String6700.Integral()
h_sig_String6700.Scale( sigmaString6700 / integral_sig_String6700)

h_sig_S83200 = file_sig_gg.Get(input_sig_S83200)
integral_sig_S83200 = h_sig_S83200.Integral()
h_sig_S83200.Scale( sigmaS83200 / integral_sig_S83200)

h_sig_Axigluon5100 = file_sig_qq.Get(input_sig_Axigluon5100)
integral_sig_Axigluon5100 = h_sig_Axigluon5100.Integral()
h_sig_Axigluon5100.Scale( sigmaAxigluon5100 / integral_sig_Axigluon5100)
##########OUTPUT########
#outputDir="fit_Run2015D_DCSonly_json_JEC_Summer15_25nsV5_974pb-1/"
#outputDir="fit_Run2015D_goldenjson_JEC_Summer15_25nsV5_547pb-1/"
#outputDir="fit_Run2015D_studies_vs_time/"
outputDir="fit_Run2015D_JEC_Summer15_25nsV5_sample45678/"
os.system("mkdir -p "+outputDir)

#================================================================================================================
  
def main():
  

  for  i in range (1, number_of_variableWidth_bins):
    #data
    bincontent = hist_binned.GetBinContent(i)
    binwidth = hist_binned.GetBinWidth(i)
    binerror = hist_binned.GetBinError(i)
    hist_mass.SetBinContent(i,bincontent/(binwidth*lumi))   
    #mc
    bincontent_mc = hist_binned_mc.GetBinContent(i)
    binwidth_mc = hist_binned_mc.GetBinWidth(i)
    hist_mass_mc.SetBinContent(i,bincontent_mc/(binwidth_mc*lumi))   

  #hist_mass.Draw()
  #filetest = TFile("filetest.root","recreate")
  #filetest.cd()
  #hist_mass.Write()
  #filetest.Close()
  
  #######################################################
  #data in TGraph format (hist binned)
  alpha = 1 - 0.6827;
  
  x=[]
  y=[]
  exl=[]
  exh=[]
  eyl=[]
  eyh=[]
  x_mc=[]
  y_mc=[]
  
  for i in range(0,number_of_variableWidth_bins):
    n    = hist_binned.GetBinContent(i+1)
    dm   = hist_binned.GetBinWidth(i+1)
    mass = hist_binned.GetBinCenter(i+1)
    xl   = hist_binned.GetBinLowEdge(i+1)
    xh   = xl+dm
    x.append( (xl+xh)/2.)
    exl.append( dm/2.)
    exh.append( dm/2.)
    y.append( n / (dm*lumi))
    l = 0.5*TMath.ChisquareQuantile(alpha/2,2*n)
    h = 0.5*TMath.ChisquareQuantile(1-alpha/2,2*(n+1))
    eyl.append( (n-l)/(lumi*dm) )
    eyh.append( (h-n)/(lumi*dm) )
    #print "%f   %f    %f    %f    %f     %f" % (x[i],y[i],exl[i],exh[i],eyl[i],eyh[i])
    n_mc = hist_binned_mc.GetBinContent(i+1)
    if (i>=41 and i<103):
      x_mc.append((xl+xh)/2.)
      y_mc.append( n_mc / (dm*lumi)) 

  vx = array("f",x)
  vy = array("f",y)
  vexl = array("f",exl)
  vexh = array("f",exh)
  veyl = array("f",eyl)
  veyh = array("f",eyh)
  vx_mc = array("f",x_mc)
  vy_mc = array("f",y_mc)

  #data in TGraph format
  g = TGraphAsymmErrors(number_of_variableWidth_bins,vx,vy,vexl,vexh,veyl,veyh)
  g.SetName("g_data")
  #g.Print()
  #mc
  g_mc = TGraph(number_of_variableWidth_bins-41,vx_mc,vy_mc)
  g_mc.SetName("g_mc")
  #g_mc.Print()
  
  
  nBins_fit = hist_mass.FindBin(maxX_mass)- hist_mass.FindBin(minX_mass) 
  ##count zero bins
  zeroBins = 0
  for i in range(0,nBins_fit):
    if hist_mass.GetBinContent(hist_mass.FindBin(minX_mass)+i)==0: 
      zeroBins +=1
  
  nBins_fit = nBins_fit-zeroBins   
  fitresult = doFitAndChi2(FunctionType,hist_mass,g,hist_mass_original)
  M1Bkg = fitresult[3]
  hist_fit_residual_vsMass = fitresult[4]
  nPar = nBins_fit - fitresult[1]
  #DrawFit(g,g_mc,M1Bkg,hist_fit_residual_vsMass,FunctionType,nPar,h_sig_S82200,h_sig_Axigluon3000,h_sig_String5300,label_S82200,label_Axigluon3000,label_String5300,2200,3000,5300,fileNameSuffix)
  DrawFit(g,g_mc,M1Bkg,hist_fit_residual_vsMass,FunctionType,nPar,h_sig_S83200,h_sig_Axigluon5100,h_sig_String6700,label_S83200,label_Axigluon5100,label_String6700,3200,5100,6700,fileNameSuffix)
  #DrawFit(g,g_mc,M1Bkg,hist_fit_residual_vsMass,FunctionType,nPar,h_sig_String5300,label_String5300,5300,fileNameSuffix_String5300)
  #DrawFit(g,g_mc,M1Bkg,hist_fit_residual_vsMass,FunctionType,nPar,h_sig_Qstar3000,label_Qstar3000,3000,fileNameSuffix_Qstar3000)
  #DrawFit(g,g_mc,M1Bkg,hist_fit_residual_vsMass,FunctionType,nPar,h_sig_Qstar4500,label_Qstar4500,4500,fileNameSuffix_Qstar4500)
  #DrawFit(g,g_mc,M1Bkg,hist_fit_residual_vsMass,FunctionType,nPar,h_sig_S82200,label_S82200,2200,fileNameSuffix_S8_2200)
  #DrawFit(g,g_mc,M1Bkg,hist_fit_residual_vsMass,FunctionType,nPar,h_sig_Axigluon3000,label_Axigluon3000,3000,fileNameSuffix)
  #print "chi2 / dof for f%d = %f / %d" % (FunctionType,fitresult[2],fitresult[1])
  
  result_mc = doChi2MC(g,hist_mass,hist_mass_mc)  
  chi2_mc = result_mc[0] 
  Ndof_mc = result_mc[1] 
  hist_mc_residual_vsMass = result_mc[2] 
  #DrawMC(g,g_mc,hist_mass,hist_mc_residual_vsMass,fileNameSuffix)


def doChi2MC(g,hist_mass,hist_mass_mc):
  hist_mc_residual_vsMass =  TH1D("hist_mc_residual_vsMass","hist_mc_residual_vsMass",number_of_variableWidth_bins,v_massBins)
  NumberOfObservations_VarBin_5entries = 0
  NumberOfObservations_VarBin = 0
  NumberOfVarBins = 0
  chi2_VarBin_5entries = 0
  chi2_VarBin_zeroes = 0
  chi2_VarBin = 0
  for bin in range (1,number_of_variableWidth_bins):
    mc_residual = 0
    hist_mc_residual_vsMass.SetBinContent(bin,mc_residual)
  
    if( hist_mass_mc.GetXaxis().GetBinLowEdge(bin)>=minX_mass and hist_mass_mc.GetXaxis().GetBinUpEdge(bin)<=maxX_mass ):
      NumberOfVarBins += 1
      #print "bin content = " + str(hist_mass.GetBinContent(bin)) + "   graph y = " + str(vy[bin-1]) + "  error y low = " + str(g.GetErrorYlow(bin-1))
      data = hist_mass.GetBinContent(bin)
      err_data_low = g.GetErrorYlow(bin-1) 
      err_data_high= g.GetErrorYhigh(bin-1)
      mc = hist_mass_mc.GetBinContent(bin)
      #print "mc = %f" % mc 
      if(mc > data or mc == 0): err_tot = err_data_high
      else: err_tot = err_data_low
      # if err_tot==0: 
      #   print "!!!! ERR = 0 !!!"
      #   print "mc %f  data %f err_data_low %f    err_data_high %f " % (mc, data,err_data_low,err_data_high)
      mc_residual = (data - mc) / err_tot
      err_mc_residual = 1
      #print "data %f   mc %f   error %f  residual %f" % (data,mc,err_tot,mc_residual)
      chi2_VarBin_zeroes += pow( (data - mc) , 2 ) / pow( err_tot , 2 )

      ##skip bin with zero entries
      if (hist_mass.GetBinContent(bin)>0): 
	NumberOfObservations_VarBin+=1
        chi2_VarBin += pow( (data - mc) , 2 ) / pow( err_tot , 2 )	 
      
      ##skip bin with less than 5 entries
      if (hist_mass.GetBinContent(bin)*hist_mass.GetBinWidth(bin)>=5): 
      #if (hist_mass.GetBinContent(bin)*lumi*hist_mass.GetBinWidth(bin)>=5): 
	NumberOfObservations_VarBin_5entries+=1
        chi2_VarBin_5entries += pow( (data - mc) , 2 ) / pow( err_tot , 2 )	 
      
      
      hist_mc_residual_vsMass.SetBinContent(bin,mc_residual)
    #print "bin : %d   mc_residual : %f" % (bin,mc_residual) 
    
  ndf_VarBin_5entries = NumberOfObservations_VarBin_5entries #-1 
  ndf_VarBin = NumberOfObservations_VarBin #-1 
  ndf_VarBin_withzeroes = NumberOfVarBins #-1 
  print "============ MC ==============" 
  print "NumberOfObservations_VarBin: %d" %  NumberOfObservations_VarBin
  print "chi2 / dof: %f / %d" % (chi2_VarBin,ndf_VarBin)
  print "chi2 / dof (with zeroes) : %f / %d " % (chi2_VarBin_zeroes,ndf_VarBin_withzeroes)
  #print "ndf_VarBin_5entries: %d" % ndf_VarBin_5entries 
  #print "chi2_VarBin_5entries: %f" % chi2_VarBin_5entries
  print "============================"   
  return [chi2_VarBin,ndf_VarBin,hist_mc_residual_vsMass]
  

def doFitAndChi2(FunctionType,hist_mass,g,hist_mass_original):
  ### fit mass histogram with background function
  # -2: VARIATION-1 (3 par.) - " [0] / ( TMath::Power(x/13000,[2]) )" 
  if( FunctionType==-2 ):    
    nPar=2
    M1Bkg = TF1("M1Bkg"," [0] / ( TMath::Power(x/13000,[1]) )",minX_mass,maxX_mass)
    M1Bkg.SetParameter(0,0.08)
    M1Bkg.SetParameter(1,2)
    M1Bkg.SetParLimits(0,0.,10);
    M1Bkg.SetParLimits(1,0.,20);
  
  
  # -1: VARIATION-1 (3 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]) )" 
  if( FunctionType==-1 ):    
    nPar=3
    M1Bkg = TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]) )",minX_mass,maxX_mass)
    M1Bkg.SetParameter(0,0.00001)
    M1Bkg.SetParameter(1,10)
    M1Bkg.SetParameter(2,5)
  
    #1 fb-1
    M1Bkg.SetParLimits(0,0.,1000);
    #10 fb-1
    #M1Bkg->SetParLimits(0,0,1.);
    #M1Bkg.SetParLimits(1,0.,100)
    #M1Bkg.SetParLimits(2,0.,20)
  
  # 0: DEFAULT (4 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )" 
  if( FunctionType==0 ):    
    nPar=4
    M1Bkg = TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )",minX_mass,maxX_mass)
    #M1Bkg.SetParameter(0,0.0014)
    M1Bkg.SetParameter(0,1.4)
    M1Bkg.SetParameter(1,12)
    M1Bkg.SetParameter(2,2)
    M1Bkg.SetParameter(3,-0.5)
  
    #M1Bkg.SetParLimits(0,0.,1000);
    #1 fb-1
    #M1Bkg.SetParLimits(0,0.,10);
    #10 fb-1
    #M1Bkg->SetParLimits(0,0,1.);
   # M1Bkg.SetParLimits(1,0.,10000)
   # M1Bkg.SetParLimits(2,1.,10000)
   # M1Bkg.SetParLimits(1,0,100)
   # M1Bkg.SetParLimits(2,0.,20)
   # M1Bkg.SetParLimits(3,-10.,10)
   # M1Bkg.FixParameter(3,0.)
     
  
  # 1: VARIATION 1 (5 par.) - "( [0]*TMath::Power(1-x/13000,[1])*(1+[4]*x/13000) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
  if( FunctionType==1 ):    
    nPar=5
    M1Bkg = TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1])*(1+[4]*x/13000) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )",minX_mass,maxX_mass)
    M1Bkg.SetParameter(0,0.005)
    M1Bkg.SetParameter(1,9.3)
    M1Bkg.SetParameter(2,7.2)
    M1Bkg.SetParameter(3,0.4)
    M1Bkg.SetParameter(4,3.1)
  #   M1Bkg.SetParLimits(4,-5,5)
  
  # 2: VARIATION 2 (6 par.) - "( [0]*TMath::Power(1-x/13000,[1])*(1+[4]*x/13000+[5]*pow(x/13000,2)) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
    #  --> 2nd order poly extension : inspired by HERA PDF 1.0 [http://arxiv.org/abs/arXiv:0911.0884 , Eq. 4.1]
  if( FunctionType==2 ):    
    nPar=6
    M1Bkg = TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1])*(1+[4]*x/13000+[5]*pow(x/13000,2)) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )",minX_mass,maxX_mass);
    M1Bkg.SetParameter(0,0.005)
    M1Bkg.SetParameter(1,9.3)
    M1Bkg.SetParameter(2,7.2)
    M1Bkg.SetParameter(3,0.4)
    M1Bkg.SetParameter(4,3.1)
    M1Bkg.SetParameter(5,25.6)
    #M1Bkg.SetParLimits(5,10,50)      
      
  
  # 3: VARIATION 3 (7 par.) - "( [0]*TMath::Power(1-x/13000,[1])*exp([4]*x/13000)*TMath::Power(1+exp([5])*x/13000,[6]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
  #    --> "exponential" extension wrt to DEFAULT - inspired by CTEQ 2008 [http://arxiv.org/pdf/hep-ph/0201195v3.pdf , Eq. 4]
  if( FunctionType==3 ) :   
    nPar=7
    M1Bkg =  TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1])*exp([4]*x/13000)*TMath::Power(1+exp([5])*x/13000,[6]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )" ,minX_mass,maxX_mass)
    M1Bkg.SetParameter(0,0.005)
    M1Bkg.SetParameter(1,15.1)
    M1Bkg.SetParameter(2,7.2)
    M1Bkg.SetParameter(3,0.4)
    M1Bkg.SetParameter(4,13.0)
    M1Bkg.SetParameter(5,-4.0)
    M1Bkg.SetParameter(6,70.0)
    # M1Bkg.SetParLimits(4,-1,1)
  
  
  # 4: VARIATION 4 (5 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)) )" 
  #    --> "log" extension wrt to DEFAULT     
  if( FunctionType==4 ) :   
    nPar=5
    M1Bkg = TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)) )",minX_mass,maxX_mass);
    M1Bkg.SetParameter(0,0.08)
    M1Bkg.SetParameter(1,12)
    M1Bkg.SetParameter(2,2)
    M1Bkg.SetParameter(3,-0.5)
    M1Bkg.SetParameter(4,0.)
    #M1Bkg.SetParLimits(3,0,0.4)
    #M1Bkg.FixParameter(3,0.)
    #M1Bkg.FixParameter(4,0.)
    #M1Bkg.SetParLimits(1,0.,10000)
    #M1Bkg.SetParLimits(2,1.,10000)
    #M1Bkg.SetParLimits(0,0.,10);
    #M1Bkg.SetParLimits(1,0.,100)
    #M1Bkg.SetParLimits(2,0.,20)
    #M1Bkg.SetParLimits(3,-10,10)
    #M1Bkg.SetParLimits(4,-10,10)
   
  
  # 5: VARIATION 5 (6 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)+[5]*TMath::Power(log(x/13000),3)) )" 
  #    --> "log" extension wrt to DEFAULT     
  if( FunctionType==5 ):    
    nPar=6
    M1Bkg = TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)+[5]*TMath::Power(log(x/13000),3)) )",minX_mass,maxX_mass)
    M1Bkg.SetParameter(0,0.08)
    M1Bkg.SetParameter(1,12)
    M1Bkg.SetParameter(2,2)
    M1Bkg.SetParameter(3,-0.5)
    M1Bkg.SetParameter(4,0.)
    M1Bkg.SetParameter(5,0.)
    #M1Bkg.SetParLimits(3,0,0.4)
    #M1Bkg.FixParameter(3,0.)
    #M1Bkg.FixParameter(4,0.)
    #M1Bkg.FixParameter(5,0.)
    #M1Bkg.SetParLimits(1,0.,10000)
    #M1Bkg.SetParLimits(2,1.,10000)
    #M1Bkg.SetParLimits(3,0.,10000)
    #M1Bkg.SetParLimits(0,0.,10);
    #M1Bkg.SetParLimits(1,1.,100)
    #M1Bkg.SetParLimits(2,1.,20)
    #M1Bkg.SetParLimits(3,-10,10)
    #M1Bkg.SetParLimits(4,-10,10)
    #M1Bkg.SetParLimits(5,-10,10)
   
  
  # 6: VARIATION 6 (7 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)+[5]*TMath::Power(log(x/13000),3)+[6]*TMath::Power(log(x/13000),4)) )" 
  #  //    --> "log" extension wrt to DEFAULT     
  if( FunctionType==6 ) :   
    nPar=7
    M1Bkg = TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)+[5]*TMath::Power(log(x/13000),3)+[6]*TMath::Power(log(x/13000),4)) )",minX_mass,maxX_mass)
    M1Bkg.SetParameter(0,0.08)
    M1Bkg.SetParameter(1,12)
    M1Bkg.SetParameter(2,2)
    M1Bkg.SetParameter(3,-0.5)
    M1Bkg.SetParameter(4,0.)
    M1Bkg.SetParameter(5,0.)
    M1Bkg.SetParameter(6,0.)
    #M1Bkg.SetParLimits(3,0,0.4)
    #M1Bkg.FixParameter(3,0.)
    #M1Bkg.FixParameter(4,0.)
    #M1Bkg.FixParameter(5,0.)
    #M1Bkg.FixParameter(6,0.)
    #M1Bkg.SetParLimits(1,0.,10000)
    #M1Bkg.SetParLimits(2,1.,10000)
    #M1Bkg.SetParLimits(3,0.,10000)
    #M1Bkg.SetParLimits(0,0.,10);
    #M1Bkg.SetParLimits(1,0.,100)
    #M1Bkg.SetParLimits(2,0.,20)
    #M1Bkg.SetParLimits(3,-10,10)
    #M1Bkg.SetParLimits(4,-10,10)
    #M1Bkg.SetParLimits(5,-10,10)
    #M1Bkg.SetParLimits(6,-10,10)
  
  
  
  #TFitResultPtr r;
  stopProgram=1;
  for loop in range (0,10):
    r = hist_mass_original.Fit("M1Bkg","ELSR","",minX_mass,maxX_mass)      
    #r = hist_mass_original.Fit("M1Bkg","ELLSR","",minX_mass,maxX_mass)      
    #r = hist_mass_original.Fit("M1Bkg","MSR","",minX_mass,maxX_mass)      
    fitStatus = int(r)
    print "fit status : %d" % fitStatus
    if(fitStatus==0):
      stopProgram=0
      r.Print("V")  
      break
   
  
  if(stopProgram==1):
    print "######################" 
    print"######################" 
    print "ERROR : Fit failed!!!!" 
    print "######################" 
    print "######################" 
  
  
  # fit residuals and chi2
  hist_fit_residual_vsMass =  TH1D("hist_fit_residual_vsMass","hist_fit_residual_vsMass",number_of_variableWidth_bins,v_massBins)
  hist_fit_residual = TH1D("hist_fit_residual","hist_fit_residual",10,-5,5)
  NumberOfVarBins = 0
  NumberOfObservations_VarBin = 0
  NumberOfObservations_VarBin_5entries = 0
  chi2_VarBin = 0.
  chi2_VarBin_notNorm = 0.
  chi2_VarBin_5entries = 0.
  chi2_VarBin_zeroes = 0.

  for bin in range (1,number_of_variableWidth_bins):
  
    if( hist_mass.GetXaxis().GetBinLowEdge(bin)>=minX_mass and hist_mass.GetXaxis().GetBinUpEdge(bin)<=maxX_mass ):
      NumberOfVarBins += 1
      #print "bin content = " + str(hist_mass.GetBinContent(bin)) + "   graph y = " + str(vy[bin-1]) + "  error y low = " + str(g.GetErrorYlow(bin-1))
      data = hist_mass.GetBinContent(bin) 
      err_data_low = g.GetErrorYlow(bin-1) 
      err_data_high= g.GetErrorYhigh(bin-1)
      fit = M1Bkg.Integral(hist_mass.GetXaxis().GetBinLowEdge(bin) , hist_mass.GetXaxis().GetBinUpEdge(bin) )
      fit = fit / ( hist_mass.GetBinWidth(bin) * lumi )
      if(fit > data): err_tot = err_data_high
      else: err_tot = err_data_low
      fit_residual = (data - fit) / err_tot
      err_fit_residual = 1
      ##skip bin with zero entries
      #print "res = "+str(pow( (data - fit) , 2 ) / pow( err_tot , 2 ))
      chi2_VarBin_zeroes += pow( (data - fit) , 2 ) / pow( err_tot , 2 )
      if (hist_mass.GetBinContent(bin)>0): 
	NumberOfObservations_VarBin+=1
        chi2_VarBin += pow( (data - fit) , 2 ) / pow( err_tot , 2 )	 
        chi2_VarBin_notNorm += pow( (data - fit) , 2 ) 	 
  
      ##skip bin with less than 5 entries
      if (hist_mass.GetBinContent(bin)*hist_mass.GetBinWidth(bin)*lumi>=5): 
	NumberOfObservations_VarBin_5entries+=1
        chi2_VarBin_5entries += pow( (data - fit) , 2 ) / pow( err_tot , 2 )	 
  
      hist_fit_residual_vsMass.SetBinContent(bin,fit_residual)
      hist_fit_residual_vsMass.SetBinError(bin,err_fit_residual)
      hist_fit_residual.Fill(fit_residual)
    
  ndf_VarBin_5entries = NumberOfObservations_VarBin_5entries - nPar# -1
  ndf_VarBin = NumberOfObservations_VarBin - nPar# -1
  ndf_VarBin_withzeroes = NumberOfVarBins - nPar# -1
  #expected = M1Bkg.Integral(5058,10072) * lumi
  expected = M1Bkg.Integral(5058,10072)
  observed = hist_binned.Integral(hist_mass.FindBin(5058),hist_mass.FindBin(10072)) 
  print "============================" 
  print "NumberOfObservations_VarBin: %d" %  NumberOfObservations_VarBin
  print "chi2 / dof: %f / %d" % (chi2_VarBin,ndf_VarBin)
  print "chi2 / dof (with zeroes) : %f / %d " % (chi2_VarBin_zeroes,ndf_VarBin_withzeroes)
  print "chi2 / dof (> 5 evt) : %f / %d " % (chi2_VarBin_5entries,ndf_VarBin_5entries)
  #print "ndf_VarBin with 5entries: %d" % ndf_VarBin_5entries
  #print "chi2_VarBin_5entries: %f" % chi2_VarBin_5entries
  #print "chi2_VarBin_notNorm: %f" % chi2_VarBin_notNorm
  print "expected events > 5 TeV (fit) : %f" % expected
  print "observed events > 5 TeV  : %f" % observed
  print "============================"   
  return [chi2_VarBin_notNorm,ndf_VarBin,chi2_VarBin,M1Bkg,hist_fit_residual_vsMass]


def DrawMC(g,g_mc,hist_mass,hist_mc_residual_vsMass,fileNameSuffix):
#  //### Draw plots
  W = 600
  H = 650
  H_ref = 650 
  W_ref = 600 
  T = 0.08*H_ref
  B = 0.12*H_ref
  L = 0.12*W_ref
  R = 0.04*W_ref
  
  c = TCanvas("c","DijetMass cross section with QCD MC",W,H)
  c.GetWindowHeight()
  c.GetWindowWidth()
  c.SetLogy()
  c.Divide(1,2,0,0,0)
  
  
  #------------ pad 1  ----------------
  c.cd(1)
  p11_1 = c.GetPad(1)
  p11_1.SetPad(0.01,0.23,0.99,0.98)
  p11_1.SetLogy()
  p11_1.SetRightMargin(0.05)
  p11_1.SetTopMargin(0.05)
  p11_1.SetFillColor(0)
  p11_1.SetBorderMode(0)
  p11_1.SetFrameFillStyle(0)
  p11_1.SetFrameBorderMode(0)
  
  #Pave text
  #pave_fit = TPaveText(0.1558691,0.30735043,0.3750171,0.4070085,"NDC")
  pave_fit = TPaveText(0.2058691,0.1235043,0.4750171,0.2870085,"NDC")
    
  pave_fit.AddText("Wide Jets")
  #pave_fit.AddText("AK4 Jets")
  pave_fit.AddText("M_{jj} > 1.2 TeV")
  pave_fit.AddText("|#eta| < 2.5, |#Delta#eta| < 1.3")
  #pave_fit.AddText("|#eta| < 2.5, |#Delta#eta| < 1.3")
  
  pave_fit.SetFillColor(0)
  pave_fit.SetLineColor(0)
  pave_fit.SetFillStyle(0)
  pave_fit.SetBorderSize(0)
  pave_fit.SetTextFont(42)
  pave_fit.SetTextSize(0.040)
  pave_fit.SetTextAlign(12) 
  
  
  vFrame = p11_1.DrawFrame(minX_mass,0.00000005,maxX_mass,5.0)
  
  vFrame.SetTitle("")
  vFrame.SetXTitle("Dijet Mass [GeV]")
  vFrame.SetYTitle("d#sigma / dm_{jj}   [pb / GeV]")
  vFrame.GetXaxis().SetTitleSize(0.06)
  vFrame.GetXaxis().SetTitleOffset(0.95)
  vFrame.GetXaxis().SetLabelSize(0.05)
  vFrame.GetYaxis().SetTitleSize(0.06)
  vFrame.GetYaxis().SetTitleOffset(0.95)
  vFrame.GetYaxis().SetLabelSize(0.05)
  
  g_mc.SetLineWidth(2)
  g_mc.SetLineColor(kBlue)
  g_mc.Draw("c")
  g.SetMarkerSize(0.9)
  g.SetMarkerStyle(20)
  g.Draw("pe0 same")
    
  leg = TLegend(0.5564991,0.55,0.9203575,0.805812)
  #leg =  TLegend(0.5564991,0.55,0.8903575,0.80)
  leg.SetTextSize(0.03546853)
  leg.SetLineColor(0)
  leg.SetLineStyle(1)
  leg.SetLineWidth(1)
  leg.SetFillColor(0)
  leg.SetFillStyle(0)
  leg.SetMargin(0.35)
  leg.AddEntry(g,"data" ,"PL")
  leg.AddEntry(g_mc,"mc","L")
  leg.Draw("same")
  pave_fit.Draw("same")
  
  # writing the lumi information and the CMS "logo"
  #  CMS_lumi( p11_1, iPeriod, iPos );
  #redraw axis
  p11_1.RedrawAxis()
  p11_1.Update()
  p11_1.GetFrame().Draw()
  #draw the lumi text on the canvas
  CMS_lumi.CMS_lumi(p11_1, iPeriod, iPos)
  
  #--- Next PAD
  
  c.cd(2)
  p11_2 = c.GetPad(2)
  p11_2.SetPad(0.01,0.02,0.99,0.24)
  p11_2.SetBottomMargin(0.35)
  p11_2.SetRightMargin(0.05)
  p11_2.SetGridx()
  p11_2.SetGridy()
  
  vFrame2 = p11_2.DrawFrame(p11_1.GetUxmin(), -4.5, p11_1.GetUxmax(), 4.5)
  
  vFrame2.SetTitle("")
  vFrame2.SetXTitle("Dijet Mass [GeV]")
  vFrame2.GetXaxis().SetTitleSize(0.06)
  vFrame2.SetYTitle("(Data-MC)/#sigma")
  vFrame2.GetYaxis().SetTitleSize(0.15)
  vFrame2.GetYaxis().SetTitleOffset(0.40)
  vFrame2.GetYaxis().SetLabelSize(0.09)
  vFrame2.GetXaxis().SetTitleSize(0.18)
  vFrame2.GetXaxis().SetTitleOffset(0.90)
  vFrame2.GetXaxis().SetLabelSize(0.15)
  
  hist_mc_residual_vsMass.GetXaxis().SetRangeUser(minX_mass,maxX_mass)
  hist_mc_residual_vsMass.GetYaxis().SetRangeUser(-2.5,2.5)
  hist_mc_residual_vsMass.SetLineWidth(0)
  hist_mc_residual_vsMass.SetFillColor(kBlue)
  hist_mc_residual_vsMass.SetLineColor(1)
  hist_mc_residual_vsMass.Draw("SAMEHIST")
  
  line = TLine(minX_mass,0,maxX_mass,0)
  line.Draw("")
  p11_2.RedrawAxis()
  line2=TLine()
  line2.DrawLine(p11_2.GetUxmin(), p11_2.GetUymax(), p11_2.GetUxmax(), p11_2.GetUymax())
  line2.DrawLine(p11_2.GetUxmax(), p11_2.GetUymin(), p11_2.GetUxmax(), p11_2.GetUymax())
  	
  ### Output files
  
  output_root_file = "dijetFitResults_MC_%s.root" % (fileNameSuffix) 
  
  f_output = TFile(output_root_file,"RECREATE")
  f_output.cd()
  g.Write()
  #hist_mass_original.Write()
  #hist_binned.Write()
  hist_mass.Write()
  c.Write()
  #r_bin->Write()
  f_output.Close()
  c_fileName = "MCandResiduals_%s.png" %(fileNameSuffix)
  c.SaveAs(outputDir+"/"+c_fileName)
  c_fileName = "MCandResiduals_%s.pdf" %(fileNameSuffix)
  c.SaveAs(outputDir+"/"+c_fileName)


def DrawFit(g,g_mc,M1Bkg,hist_fit_residual_vsMass,FunctionType,nPar,h_sig,h_sig2,h_sig3,label_sig,label_sig2,label_sig3,mass_sig,mass_sig2,mass_sig3,fileNameSuffix):
#def DrawFit(g,g_mc,M1Bkg,hist_fit_residual_vsMass,FunctionType,nPar,h_sig,label_sig,mass_sig,fileNameSuffix):

  ##rescale fit function
  M1Bkg_xsec = TF1("M1Bkg_xsec","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )",minX_mass,maxX_mass)
  M1Bkg_xsec.SetParameter(0,M1Bkg.GetParameter(0)/lumi) 
  M1Bkg_xsec.SetParameter(1,M1Bkg.GetParameter(1)) 
  M1Bkg_xsec.SetParameter(2,M1Bkg.GetParameter(2)) 
  M1Bkg_xsec.SetParameter(3,M1Bkg.GetParameter(3)) 
  M1Bkg_xsec.SetLineColor(2)
  M1Bkg_xsec.SetLineWidth(2)

  nbins_sig = h_sig.GetNbinsX()
  massBinsSig_list=[]
  for i in range(0,nbins_sig+1):
    massBinsSig_list.append(h_sig.GetXaxis().GetBinLowEdge(i+1))  

  massBinsSig = array("d",massBinsSig_list)
  h_w = TH1D("h_w","", h_sig.GetNbinsX(), massBinsSig)
  h_w2 = TH1D("h_w2","", h_sig2.GetNbinsX(), massBinsSig)
  h_w3 = TH1D("h_w3","", h_sig3.GetNbinsX(), massBinsSig)
  for i in range(1,nbins_sig+1):
    bincontent = h_sig.GetBinContent(i) / h_sig.GetBinWidth(i)
    h_w.SetBinContent(i, bincontent )
    bincontent2 = h_sig2.GetBinContent(i) / h_sig2.GetBinWidth(i)
    h_w2.SetBinContent(i, bincontent2 )
    bincontent3 = h_sig3.GetBinContent(i) / h_sig3.GetBinWidth(i)
    h_w3.SetBinContent(i, bincontent3 )
 
  h_w.GetXaxis().SetRangeUser(mass_sig*0.8, mass_sig*1.2)
  h_w2.GetXaxis().SetRangeUser(mass_sig2*0.8, mass_sig2*1.2)
  h_w3.GetXaxis().SetRangeUser(mass_sig3*0.8, mass_sig3*1.2)
  #h_w.Print();
   
  #  //### Draw plots
  W = 600
  H = 700
  H_ref = 700 
  W_ref = 600 
  T = 0.08*H_ref
  B = 0.12*H_ref
  L = 0.12*W_ref
  R = 0.04*W_ref
  
  c = TCanvas("c","DijetMass cross section with Fit and QCD MC",W,H)
  c.GetWindowHeight()
  c.GetWindowWidth()
  c.SetLogy()
  c.Divide(1,2,0,0,0)
  
  
  #------------ pad 1  ----------------
  c.cd(1)
  p11_1 = c.GetPad(1)
  p11_1.SetPad(0.01,0.26,0.99,0.98)
  p11_1.SetLogy()
  p11_1.SetRightMargin(0.05)
  p11_1.SetTopMargin(0.05)
  p11_1.SetFillColor(0)
  p11_1.SetBorderMode(0)
  p11_1.SetFrameFillStyle(0)
  p11_1.SetFrameBorderMode(0)
  
  #Pave text
  #pave_fit = TPaveText(0.1558691,0.30735043,0.3750171,0.4070085,"NDC")
  pave_fit = TPaveText(0.2358691,0.04035043,0.5050171,0.1870085,"NDC")
    
  pave_fit.AddText("Wide Jets")
  pave_fit.AddText("|#eta| < 2.5, |#Delta#eta_{jj}| < 1.3")
  #pave_fit.AddText("M_{jj} > 1.2 TeV")
  pave_fit.SetFillColor(0)
  pave_fit.SetLineColor(0)
  pave_fit.SetFillStyle(0)
  pave_fit.SetBorderSize(0)
  pave_fit.SetTextFont(42)
  pave_fit.SetTextSize(0.040)
  pave_fit.SetTextAlign(12) 
  
  
  vFrame = p11_1.DrawFrame(minX_mass,0.00000005,maxX_mass,5.0)
  
  vFrame.SetTitle("")
  #vFrame.SetXTitle("Dijet Mass AK4 [GeV]")
  vFrame.SetXTitle("Dijet Mass [GeV]")
  vFrame.SetYTitle("d#sigma / dm_{jj}   [pb / GeV]")
  vFrame.GetXaxis().SetTitleSize(0.06)
  vFrame.GetXaxis().SetTitleOffset(0.95)
  vFrame.GetXaxis().SetLabelSize(0.05)
  vFrame.GetYaxis().SetTitleSize(0.06)
  #vFrame.GetYaxis().SetTitleOffset(1.0)
  vFrame.GetYaxis().SetLabelSize(0.05)
 
  g.GetXaxis().SetNdivisions(405)
  g.SetMarkerSize(0.9)
  g.SetMarkerStyle(20)
  #g.Draw("pe0")
  h_w.SetLineColor(kMagenta-2)
  h_w.SetLineWidth(2)
  h_w.SetLineStyle(4)
  h_w2.SetLineColor(kAzure-7)
  h_w2.SetLineWidth(2)
  h_w2.SetLineStyle(6)
  h_w3.SetLineColor(kTeal+3)
  h_w3.SetLineWidth(2)
  h_w3.SetLineStyle(8)
  if plotSig: 
    h_w.Draw("c same")
    h_w2.Draw("c same")
    h_w3.Draw("c same")
  g_mc.SetLineWidth(2)
  g_mc.SetLineStyle(2)
  g_mc.SetLineColor(kBlue)
  g_mc.Draw("c")
  g.SetMarkerSize(0.9)
  g.SetMarkerStyle(20)
  g.Draw("pe0 same")
  #M1Bkg.SetLineWidth(2)
  #M1Bkg.SetLineStyle(1)
  #M1Bkg.SetLineColor(2)
  #M1Bkg.Draw("same")
  M1Bkg_xsec.Draw("same")
  g.Draw("pe0 same")
    
  #leg = TLegend(0.5564991,0.55,0.8903575,0.705812)
  leg =  TLegend(0.4564991,0.62,0.9303575,0.90)
  leg.SetTextSize(0.038)
  leg.SetLineColor(0)
  leg.SetLineStyle(1)
  leg.SetLineWidth(1)
  leg.SetFillColor(0)
  leg.SetFillStyle(0)
  leg.SetMargin(0.35)
  leg.AddEntry(hist_mass,"Data" ,"PLE");
  leg.AddEntry(M1Bkg_xsec,"Fit","L");
  leg.AddEntry(g_mc,"QCD MC","L");
  if plotSig:
    leg.AddEntry(h_w,label_sig,"L");
    leg.AddEntry(h_w2,label_sig2,"L");
    leg.AddEntry(h_w3,label_sig3,"L");
  leg.Draw("same")
  pave_fit.Draw("same")
  
  # writing the lumi information and the CMS "logo"
  #  CMS_lumi( p11_1, iPeriod, iPos );
  #redraw axis
  p11_1.RedrawAxis()
  p11_1.Update()
  p11_1.GetFrame().Draw()
  #draw the lumi text on the canvas
  CMS_lumi.CMS_lumi(p11_1, iPeriod, iPos)
  
  #--- Next PAD
  
  c.cd(2)
  p11_2 = c.GetPad(2)
  p11_2.SetPad(0.01,0.02,0.99,0.27)
  p11_2.SetBottomMargin(0.35)
  p11_2.SetRightMargin(0.05)
  p11_2.SetGridx()
  p11_2.SetGridy()
  
  vFrame2 = p11_2.DrawFrame(p11_1.GetUxmin(), -2.2, p11_1.GetUxmax(), 2.2)
  
  vFrame2.SetTitle("")
  #vFrame2.SetXTitle("Dijet Mass AK4 [GeV]")
  vFrame2.SetXTitle("Dijet Mass [GeV]")
  vFrame2.GetXaxis().SetTitleSize(0.06)
  vFrame2.SetYTitle("#frac{(Data-Fit)}{#sigma_{Data}}")
  vFrame2.GetYaxis().SetTitleSize(0.15)
  vFrame2.GetYaxis().SetTitleOffset(0.40)
  vFrame2.GetYaxis().SetLabelSize(0.09)
  vFrame2.GetXaxis().SetTitleSize(0.15)
  vFrame2.GetXaxis().SetTitleOffset(0.90)
  vFrame2.GetXaxis().SetLabelSize(0.12)
  
  hist_fit_residual_vsMass.GetXaxis().SetNdivisions(405)
  hist_fit_residual_vsMass.GetXaxis().SetRangeUser(minX_mass,maxX_mass)
  hist_fit_residual_vsMass.GetYaxis().SetRangeUser(-2.5,2.5)
  hist_fit_residual_vsMass.SetLineWidth(0)
  hist_fit_residual_vsMass.SetFillColor(2)
  hist_fit_residual_vsMass.SetLineColor(1)
  hist_fit_residual_vsMass.Draw("SAMEHIST")
 
  hist_sig_significance = histSignificance(M1Bkg,hist_mass,h_sig)
  hist_sig_significance2 = histSignificance(M1Bkg,hist_mass,h_sig2)
  hist_sig_significance3 = histSignificance(M1Bkg,hist_mass,h_sig3)

  hist_sig_significance.GetXaxis().SetRangeUser(mass_sig*0.8,mass_sig*1.2)
  hist_sig_significance.SetLineColor(kMagenta-2)
  hist_sig_significance.SetLineWidth(2)
  hist_sig_significance.SetLineStyle(4)
  hist_sig_significance2.GetXaxis().SetRangeUser(mass_sig2*0.8,mass_sig2*1.2)
  hist_sig_significance2.SetLineColor(kAzure-7)
  hist_sig_significance2.SetLineWidth(2)
  hist_sig_significance2.SetLineStyle(6)
  hist_sig_significance3.GetXaxis().SetRangeUser(mass_sig3*0.8,mass_sig3*1.2)
  hist_sig_significance3.SetLineColor(kTeal+3)
  hist_sig_significance3.SetLineWidth(2)
  hist_sig_significance3.SetLineStyle(8)
  
  if plotSig:
    hist_sig_significance.Draw("c same")
    hist_sig_significance2.Draw("c same")
    hist_sig_significance3.Draw("c same")

  line = TLine(minX_mass,0,maxX_mass,0)
  line.Draw("")
  p11_2.RedrawAxis()
  line2=TLine()
  line2.DrawLine(p11_2.GetUxmin(), p11_2.GetUymax(), p11_2.GetUxmax(), p11_2.GetUymax())
  line2.DrawLine(p11_2.GetUxmax(), p11_2.GetUymin(), p11_2.GetUxmax(), p11_2.GetUymax())
  	
  ### Output files
  
  output_root_file = "dijetFitResults_FuncType%d_nParFit%d_%s.root" % (FunctionType,nPar,fileNameSuffix) 
  
  f_output = TFile(output_root_file,"RECREATE")
  f_output.cd()
  g.Write()
  #hist_mass_original.Write()
  #hist_binned.Write()
  hist_mass.Write()
  c.Write()
  #r_bin->Write()
  f_output.Close()
  c_fileName = "fitAndResiduals_FuncType%d_nParFit%d_%s.png" %(FunctionType,nPar,fileNameSuffix)
  c.SaveAs(outputDir+"/"+c_fileName)
  c_fileName = "fitAndResiduals_FuncType%d_nParFit%d_%s.pdf" %(FunctionType,nPar,fileNameSuffix)
  c.SaveAs(outputDir+"/"+c_fileName)
  c_fileName = "fitAndResiduals_FuncType%d_nParFit%d_%s.root" %(FunctionType,nPar,fileNameSuffix)
  c.SaveAs(outputDir+"/"+c_fileName)
  


def histSignificance(M1Bkg,hist_mass,h_sig):
  nbins_sig = h_sig.GetNbinsX()
  massBinsSig_list=[]
  for i in range(0,nbins_sig+1):
    massBinsSig_list.append(h_sig.GetXaxis().GetBinLowEdge(i+1))  

  massBinsSig = array("d",massBinsSig_list)
  ##histogram signal significance  
  hist_sig_significance = TH1D("hist_sig_significance","hist_sig_significance",nbins_sig,massBinsSig)
  
  for  i in range(1,number_of_variableWidth_bins+1):
    if( hist_mass.GetXaxis().GetBinLowEdge(i)>=minX_mass and hist_mass.GetXaxis().GetBinUpEdge(i)<=maxX_mass ):
      fit = M1Bkg.Integral(hist_mass.GetXaxis().GetBinLowEdge(i), hist_mass.GetXaxis().GetBinUpEdge(i) )
      significance = 0
    
      #if((h_sig.GetBinContent(h_sig.FindBin(hist_mass.GetBinLowEdge(i)))+fit) != 0 and h_sig.GetBinLowEdge(h_sig.FindBin(hist_mass.GetBinLowEdge(i)))<12000):
      #significance  = h_sig.GetBinContent(h_sig.FindBin(hist_mass.GetBinLowEdge(i)))*lumi / TMath.Sqrt(h_sig.GetBinContent(h_sig.FindBin(hist_mass.GetBinLowEdge(i)))*lumi + fit*lumi)
      ## fit does not have to be multiplied by lumi 
      significance  = h_sig.GetBinContent(h_sig.FindBin(hist_mass.GetBinLowEdge(i)))*lumi / TMath.Sqrt(h_sig.GetBinContent(h_sig.FindBin(hist_mass.GetBinLowEdge(i)))*lumi + fit)

      #print  "low edge bin: " +str( h_sig.GetBinLowEdge(i)) +"  bin low edge bkg = "+str(hist_mass.GetBinLowEdge(i))+ "  "+str( fit) + " + " +str( h_sig.GetBinContent(i))+ "  sqrt(fit + sig) = " +str( TMath.Sqrt(h_sig.GetBinContent(i) + fit))
      #print "significance = "+ str(significance)
      hist_sig_significance.SetBinContent(i,significance);

  return hist_sig_significance



#if __name__ == "__main__":

#----- keep the GUI alive ------------
if __name__ == '__main__':
  main() 
#  rep = ''
#  while not rep in ['q','Q']:
#    rep = raw_input('enter "q" to quit: ')
#    if 1 < len(rep):
#      rep = rep[0]
#

