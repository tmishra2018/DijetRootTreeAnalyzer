#!usr/bin/python
import ROOT
from ROOT import *
from setTDRStyle import setTDRStyle
import numpy as np
from array import array 
from ctypes import *

import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("--inputfile",action="store",type="string",dest="inputfile",default='histo.root')
parser.add_option("--inputhisto",action="store",type="string",dest="inputhisto",default='hist_allCutsQCD')
parser.add_option("--nFits",action="store",type="int",dest="nFits",default='1000')
(options, args) = parser.parse_args()

inputfile = options.inputfile
inputhisto = options.inputhisto
nFits = options.nFits

gROOT.Reset()
setTDRStyle()
gROOT.ForceStyle()
gROOT.SetStyle('tdrStyle')
gStyle.SetOptStat(1111)

minX_mass = 1118.
maxX_mass = 6099.

number_of_variableWidth_bins = 93
#number_of_variableWidth_bins =2 
#list_massBins = [0,3854,10000]
list_massBins = [0, 1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687,  1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509,  4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10000]

#massBins = np.asarray(list_massBins)
massBins = array('d')
massBins.fromlist(list_massBins)

print massBins

file0 = TFile(inputfile)
hist_mass_original = file0.Get( inputhisto ) 
#hist_mass_original.Sumw2()

hist_binned= TH1F(hist_mass_original.Rebin(number_of_variableWidth_bins,"hist_binned",massBins))

hist_mass = hist_binned.Clone("hist_mass")
hist_mass.GetXaxis().SetTitle("M_{jj} WideJets [GeV]")
for i in range(0,number_of_variableWidth_bins):
  bincontent = hist_binned.GetBinContent(i)
  binwidth = hist_binned.GetBinWidth(i)
  binerror = hist_binned.GetBinError(i)
  hist_mass.SetBinContent(i,bincontent/binwidth)   
  hist_mass.SetBinError(i,binerror/binwidth)   
  
file_debug = TFile("chi2_studies_1fb-1.root","recreate")
file_debug.cd()
#hist_mass_original.Write()
#hist_binned.Write()
#hist_mass.Write()



nPar=4
M1Bkg =  TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )",minX_mass,maxX_mass)
M1Bkg.SetParameter(0,0.000884)
M1Bkg.SetParameter(1,5.497)
M1Bkg.SetParameter(2,6.63)
M1Bkg.SetParameter(3,0.281)

#M1Bkg.SetParLimits(0,0,0.002)
#M1Bkg.SetParLimits(1,0,10)
#M1Bkg.SetParLimits(2,5.5,8.)
#M1Bkg.SetParLimits(3,0,0.5)
M1Bkg.SetParLimits(0,0,10000)
M1Bkg.SetParLimits(1,0,100)
M1Bkg.SetParLimits(2,0,100)
M1Bkg.SetParLimits(3,0,1)


##Fit the nominal distribution

print("start the fit")
stopProgram = 1
for loop in range(0,10):
           
  r_bin = TFitResultPtr(hist_mass.Fit("M1Bkg","ILSR","",minX_mass,maxX_mass))      
  r_bin.SetName("fit_result_M1Bkg_binned")
  fitStatus = r_bin
  if(fitStatus!=False):
    stopProgram=0
    r_bin.Print("V") 
    break
    
  if(stopProgram==1):
 
    print( "######################" )
    print( "######################" )
    print( "ERROR : Fit failed!!!!" )
    print( "######################" )
    print( "######################" )
    break

fit_results= []  
pseudodatasets = []
numEv = hist_mass_original.Integral()
h_pull_p0 = TH1F("h_pull_p0", "", 100, -5., 5.)
h_pull_p1 = TH1F("h_pull_p1", "", 100, -5., 5.)
h_pull_p2 = TH1F("h_pull_p2", "", 100, -5., 5.)
h_pull_p3 = TH1F("h_pull_p3", "", 100, -5., 5.)
h_chi2 = TH1F("h_chi2","",80,0,80)
h_ll = TH1F("h_ll","",100,0,1)

for i_fit in range(0,nFits):
  pseudodataset_1GeV = TH1F("pseudodataset_1GeV_"+str(i_fit),"",10000,0.,10000.)
  #pseudodataset_1GeV.FillRandom("M1Bkg",int(numEv))
  pseudodataset_1GeV.FillRandom(hist_mass,int(numEv))
  pseudodataset_binned = pseudodataset_1GeV.Rebin(number_of_variableWidth_bins,"pseudodataset_binned_"+str(i_fit),massBins)
  pseudodataset = TH1F("pseudodataset_binned_w_"+str(i_fit), "", number_of_variableWidth_bins, massBins)	
  for i in range(0,number_of_variableWidth_bins): 
    bincontent = pseudodataset_binned.GetBinContent(i)
    binwidth = pseudodataset_binned.GetBinWidth(i)
    binerror = pseudodataset_binned.GetBinError(i)
    pseudodataset.SetBinContent(i, bincontent / binwidth)
    pseudodataset.SetBinError(i, binerror / binwidth)
    if pseudodataset.GetBinContent(i)==0:
      pseudodataset.SetBinError(i, 1.8 / binwidth)

  pseudodatasets.append(pseudodataset) 
  fit_results.append( pseudodataset.Fit("M1Bkg","ILSR","",minX_mass,maxX_mass) )
  fit_results[i_fit].SetName("fit_result"+str(i_fit))
  p0 = fit_results[i_fit].Parameter(0)
  p1 = fit_results[i_fit].Parameter(1)
  p2 = fit_results[i_fit].Parameter(2)
  p3 = fit_results[i_fit].Parameter(3)
  errors = fit_results[i_fit].GetErrors()
  chi2 = fit_results[i_fit].Chi2()
  ll = fit_results[i_fit].MinFcnValue()
  exp_p0 = r_bin.Parameter(0)
  exp_p1 = r_bin.Parameter(1)
  exp_p2 = r_bin.Parameter(2)
  exp_p3 = r_bin.Parameter(3)
  pull_p0 = (p0 - exp_p0) / errors[0] 
  pull_p1 = (p1 - exp_p1) / errors[1] 
  pull_p2 = (p2 - exp_p2) / errors[2] 
  pull_p3 = (p3 - exp_p3) / errors[3]

  #Fill pull histograms
  h_pull_p0.Fill(pull_p0)
  h_pull_p1.Fill(pull_p1)
  h_pull_p2.Fill(pull_p2)
  h_pull_p3.Fill(pull_p3)
  h_chi2.Fill(chi2)
  h_ll.Fill(ll)

#chi2 nominal function
chi2= TF1("chi2","ROOT::Math::chisquared_pdf(x,35,0)",0,100);

#Write histograms
M1Bkg.Write()
h_pull_p0.Write()
h_pull_p1.Write()
h_pull_p2.Write()
h_pull_p3.Write()
h_chi2.Write()
chi2.Write()
h_ll.Write()

c = TCanvas("c","",800,800)
leg_chi2 = TLegend(0.6,0.6,0.86,0.75)                                            
leg_chi2.SetFillColor(0)                                                         
leg_chi2.AddEntry(h_chi2,"toy distribution","l")                                 
leg_chi2.AddEntry(chi2,"#chi^{2} function for 35 ndf","l")                       
leg_chi2.Draw("")

c.cd()
chi2.SetLineColor(2)
chi2.SetLineWidth(2)
h_chi2.DrawNormalized("")
chi2.Draw("same")
leg_chi2.Draw()
c.SaveAs("chi2_study_1fb-1.png")
c.SaveAs("chi2_study_1fb-1.pdf")

c.Write()

file_debug.Close()

