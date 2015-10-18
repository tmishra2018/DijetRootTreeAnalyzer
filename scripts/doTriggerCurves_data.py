#! /usr/bin/env python
import os
import sys
import string
import re
import argparse
import math
from ROOT import *
from array import array
import CMS_lumi
from setTDRStyle import setTDRStyle


gROOT.SetBatch(kTRUE);
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetTitleFont(42, "XYZ")
gStyle.SetTitleSize(0.06, "XYZ")
gStyle.SetLabelFont(42, "XYZ")
gStyle.SetLabelSize(0.05, "XYZ")
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadTopMargin(0.05)
gStyle.SetPadBottomMargin(0.15)
gROOT.ForceStyle()

gROOT.Reset()
setTDRStyle()
gROOT.ForceStyle()
gROOT.SetStyle('tdrStyle')


#################
number_of_variableWidth_bins = 103

massBins =[1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430,10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000];

v_massBins = array("d",massBins)
#######################################################

usage = "usage: "

parser = argparse.ArgumentParser(description='Process options.')

parser.add_argument("-i", "--inputList", type=str, dest="inputList", default="",
    help="input list of files to be merged",
    )
parser.add_argument("-o", "--outputDir", type=str, dest="outputDir", default="./",
    help="output dir",
    )
parser.add_argument("--lumi", type=float, dest="LUMI", default="1000",
    help="luminosity in pb",
    )
parser.add_argument("--mode", type=int, dest="MODE", default="1",
    help="mode for eff asymmetric uncertainties",
    )

args = parser.parse_args()
print args 

inputList = args.inputList
LUMI = args.LUMI
MODE = args.MODE
outputDir = args.outputDir
#####################

minX_mass = 526
maxX_mass = 2231 
minX_fit = 1027 
maxX_fit = 1530 

CMS_lumi.extraText = "Preliminary"
#CMS_lumi.lumi_sqrtS = str(int(LUMI))+" pb^{-1} (13 TeV)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
CMS_lumi.lumi_sqrtS = "2015, 13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 11
iPeriod = 0


#---- read the list -----------------
lines = [line.strip() for line in open(inputList)]

#---- split sample name and xsec
fileNames = []
xsecs = []
ii = 0
for line in lines:
  parts = line.split()
  fileNames.append(parts[0])
  #xsecs.append(parts[1])
  #print ("dataset : %s    xsec : %s" % (fileNames[ii], xsecs[ii]))
  print ("dataset : %s " % fileNames[ii])
  ii+=1

#---- open the files --------------------
h_mjj_HLTpass_noTrig_list = []
h_mjj_HLTpass_PFHT475_list = []
h_mjj_HLTpass_PFHT800_list = []
h_mjj_HLTpass_PFHT650MJJ900_list = []
h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_list = []
h_mjj_HLTpass_PFHT800_noPFHT475_list = []

print fileNames
filelist =[]
for f in fileNames:
  inf = TFile.Open(f)
  filelist.append(inf)
  
i_f = 0

####### mjj histograms matching HLT
for inf in filelist:
 
  Nev = inf.Get('DijetFilter/EventCount/EventCounter').GetBinContent(1)
  print ('processed events: %s' % Nev)
  wt = 1.0

  h_mjj_HLTpass_noTrig = inf.Get("h_mjj_HLTpass_noTrig")
  h_mjj_HLTpass_PFHT475 = inf.Get("h_mjj_HLTpass_PFHT475")
  h_mjj_HLTpass_PFHT800 = inf.Get("h_mjj_HLTpass_PFHT800")
  h_mjj_HLTpass_PFHT650MJJ900 = inf.Get("h_mjj_HLTpass_PFHT650MJJ900")
  h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900 = inf.Get("h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900")
  h_mjj_HLTpass_PFHT800_noPFHT475 = inf.Get("h_mjj_HLTpass_PFHT800_noPFHT475")

  h_mjj_HLTpass_noTrig_list.append(h_mjj_HLTpass_noTrig)
  h_mjj_HLTpass_PFHT475_list.append(h_mjj_HLTpass_PFHT475 )
  h_mjj_HLTpass_PFHT800_list.append(h_mjj_HLTpass_PFHT800)
  h_mjj_HLTpass_PFHT650MJJ900_list.append(h_mjj_HLTpass_PFHT650MJJ900)
  h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_list.append(h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900)
  h_mjj_HLTpass_PFHT800_noPFHT475_list.append(h_mjj_HLTpass_PFHT800_noPFHT475)
  #running on data
  #wt = LUMI*float(xsecs[i_f])/Nev
  print('weight : %f' % wt)
  h_mjj_HLTpass_noTrig_list[i_f].Scale(wt)
  h_mjj_HLTpass_PFHT475_list[i_f].Scale(wt)
  h_mjj_HLTpass_PFHT800_list[i_f].Scale(wt) 
  h_mjj_HLTpass_PFHT650MJJ900_list[i_f].Scale(wt)
  h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_list[i_f].Scale(wt)
  h_mjj_HLTpass_PFHT800_noPFHT475_list[i_f].Scale(wt) 

  h_mjj_HLTpass_PFHT475_list[i_f].Print()
  print str(i_f) 
  h_mjj_HLTpass_PFHT475_list[0].Print()

  i_f += 1

#h_mjj_HLTpass_PFHT475_list[0].Print()
#print h_mjj_HLTpass_PFHT475_list

h_mjj_HLTpass_noTrig_all = h_mjj_HLTpass_PFHT475_list[0].Clone("h_mjj_HLTpass_noTrig_all")
h_mjj_HLTpass_PFHT475_all = h_mjj_HLTpass_PFHT475_list[0].Clone("h_mjj_HLTpass_PFHT475_all")
h_mjj_HLTpass_PFHT800_all = h_mjj_HLTpass_PFHT800_list[0].Clone("h_mjj_HLTpass_PFHT800_all")
h_mjj_HLTpass_PFHT650MJJ900_all = h_mjj_HLTpass_PFHT650MJJ900_list[0].Clone("h_mjj_HLTpass_PFHT650MJJ900_all")
h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_all = h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_list[0].Clone("h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_all")
h_mjj_HLTpass_PFHT800_noPFHT475_all = h_mjj_HLTpass_PFHT800_noPFHT475_list[0].Clone("h_mjj_HLTpass_PFHT800_noPFHT475_all")

if(i_f>1):
  for i in range(1,len(fileNames)):
    h_mjj_HLTpass_noTrig_all.Add(h_mjj_HLTpass_noTrig_list[i])
    h_mjj_HLTpass_PFHT475_all.Add(h_mjj_HLTpass_PFHT475_list[i])
    h_mjj_HLTpass_PFHT800_all.Add(h_mjj_HLTpass_PFHT800_list[i])
    h_mjj_HLTpass_PFHT650MJJ900_all.Add(h_mjj_HLTpass_PFHT650MJJ900_list[i])
    h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_all.Add(h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_list[i])
    h_mjj_HLTpass_PFHT800_noPFHT475_all.Add(h_mjj_HLTpass_PFHT800_noPFHT475_list[i])

h_mjj_HLTpass=[h_mjj_HLTpass_noTrig_all,
    h_mjj_HLTpass_PFHT475_all,
    h_mjj_HLTpass_PFHT800_all, 
    h_mjj_HLTpass_PFHT650MJJ900_all, 
    h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_all,
    h_mjj_HLTpass_PFHT800_noPFHT475_all
    ]

g_eff = []
nMassBins = h_mjj_HLTpass_PFHT475_all.GetNbinsX()
scale=1.

##### efficiency and error for trigger curves
for ii in range(2,5):
  x = []
  y = []
  exl = []
  exh = []
  eyl = []
  eyh = []
  for i in range(0,nMassBins):
    N1 = h_mjj_HLTpass[1].GetBinContent(i+1)
    N2 = h_mjj_HLTpass[ii].GetBinContent(i+1)
    p  = 0
    eU = 0
    eL = 0
    if (N1 > 0):
      p = N2/N1
      n = N1+N2
      w = N2/n
      if (MODE==1): ## Wilson for binomial
	scale = 1.0 # makes sense only for the unprescaled trigger 
        d = math.sqrt(p*(1-p)/N1+0.25/(N1*N1))
        eU = (p+0.5/N1+d)/(1+1/N1)-p
        eL = p-(p+0.5/N1-d)/(1+1/N1)
         
      else: ## Wilson for Poisson ratio
        d  = math.sqrt(w*(1-w)/n+0.25/(n*n))
        UB = (w+0.5/n+d)/(1+1/n)
        LB = (w+0.5/n-d)/(1+1/n)
        eU = UB/(1-UB)-p
        eL = p-LB/(1-LB)     
    print "mass:"+str(h_mjj_HLTpass[0].GetBinCenter(i+1))+"  "+str(N1)+" "+str(N2)+" "+str(p)+" "+str(eL)+" "+str(eU)

    x.append( h_mjj_HLTpass[0].GetBinCenter(i+1))
    y.append( p*scale)
    exl.append( h_mjj_HLTpass[0].GetBinWidth(i+1)/2)
    exh.append( h_mjj_HLTpass[0].GetBinWidth(i+1)/2)      
    eyl.append( eL*scale)
    eyh.append( eU*scale)
       
  vx = array("f",x)
  vy = array("f",y)
  vexl = array("f",exl)
  vexh = array("f",exh)
  veyl = array("f",eyl)
  veyh = array("f",eyh)

  g_eff.append(TGraphAsymmErrors(nMassBins,vx,vy,vexl,vexh,veyl,veyh))  

   
g_eff[0].SetName("g_eff_PFHT800")
g_eff[1].SetName("g_eff_PFHT650MJJ900")
g_eff[2].SetName("g_eff_PFHT800_OR_MJJ900")

####fit efficiency
#turnon = TF1("turnon","(1+TMath::TanH([0]+[1]*x))/2",minX_fit,maxX_fit)
##turnon = TF1("turnon","(TMath::TanH([0]+[1]*x))",900,1500)
#turnon.SetParameter(0,-25)
#turnon.SetParameter(1,0.045)
##turnon.FixParameter(0,-25)
##turnon.FixParameter(1,0.025)
#fermi = TF1("fermi","1/(1+TMath::Exp([0]*(x-[1])))",minX_fit,maxX_fit)
#fermi.SetParameter(1,1100)
#
#stopProgram=1
#for loop in range (0,10):
#  res = g_eff[0].Fit("fermi","SR","",minX_fit,maxX_fit)
#  #res = g_eff[0].Fit("turnon","SR","",minX_fit,maxX_fit)
#  #res = g_eff[0].Fit("turnon","SR","",minX_mass,maxX_mass)
#  fitStatus = int(res)
#  print "fit status : %d" % fitStatus
#  if(fitStatus==0):
#    stopProgram=0
#    res.Print("V")
#    break
#					    
#if(stopProgram==1):
#  print "######################"
#  print"######################"
#  print "ERROR : Fit failed!!!!"
#  print "######################"
#  print "######################"
#
## fit residuals and chi2
#hist_fit_residual_vsMass =  TH1D("hist_fit_residual_vsMass","hist_fit_residual_vsMass",number_of_variableWidth_bins,v_massBins)
#NumberOfVarBins = 0
#NumberOfObservations_VarBin = 0
#chi2_VarBin = 0.
#chi2_VarBin_notNorm = 0.
#chi2_VarBin_zeroes = 0.
#nPar=2
#
#for bin in range (1,number_of_variableWidth_bins):
#
#  if( v_massBins[bin-1]>=minX_fit and v_massBins[bin]<=maxX_fit ):
#    NumberOfVarBins += 1
#    #print "bin content = " + str(hist_mass.GetBinContent(bin)) + "   graph y = " + str(vy[bin-1]) + "  error y low = " + str(g.GetErrorYlow(bin-1))
#    data = g_eff[0].GetY()[bin-1]
#    err_data_low = g_eff[0].GetErrorYlow(bin-1) 
#    err_data_high= g_eff[0].GetErrorYhigh(bin-1)
#    #fit = turnon.Integral(v_massBins[bin-1] , v_massBins[bin] )
#    fit = fermi.Integral(v_massBins[bin-1] , v_massBins[bin] )
#    fit = fit / ( v_massBins[bin]-v_massBins[bin-1] )
#    if(fit > data): err_tot = err_data_high
#    else: err_tot = err_data_low
#    fit_residual = (data - fit) / err_tot
#    err_fit_residual = 1
#    ##skip bin with zero entries
#    chi2_VarBin_zeroes += pow( (data - fit) , 2 ) / pow( err_tot , 2 )
#    if (g_eff[0].GetY()[bin-1]>0): 
#      NumberOfObservations_VarBin+=1
#      chi2_VarBin += pow( (data - fit) , 2 ) / pow( err_tot , 2 )	 
#      chi2_VarBin_notNorm += pow( (data - fit) , 2 ) 	 
#
#    hist_fit_residual_vsMass.SetBinContent(bin,fit_residual)
#    hist_fit_residual_vsMass.SetBinError(bin,err_fit_residual)
#  
#ndf_VarBin = NumberOfObservations_VarBin - nPar# -1
#ndf_VarBin_withzeroes = NumberOfVarBins - nPar# -1
#print "============================" 
#print "NumberOfObservations_VarBin: %d" %  NumberOfObservations_VarBin
#print "ndf_VarBin: %d" % ndf_VarBin 
#print "ndf_VarBin with zeroes: %d" % ndf_VarBin_withzeroes 
#print "chi2_VarBin with zeroes: %f" % chi2_VarBin_zeroes
#print "chi2_VarBin: %f" % chi2_VarBin
#print "chi2_VarBin_notNorm: %f" % chi2_VarBin_notNorm
#print "============================"   
#
#################   

g_ineff = []
### ineff / relative_stat_unc
for ii in range(2,5):
  x = []
  y = []
  exl = []
  exh = []
  eyl = []
  eyh = []
  for i in range(0,nMassBins):
    N1 = h_mjj_HLTpass[1].GetBinContent(i+1)
    N2 = h_mjj_HLTpass[ii].GetBinContent(i+1)
    N = h_mjj_HLTpass[5].GetBinContent(i+1)
    ineff = 1
    p = 0
    rel_stat_unc = 0.000001    
    eU = 0
    eL = 0
    if (N1>0 and N2>0) :
      p = N2/N1
      ineff = 1-p
      rel_stat_unc = 1/math.sqrt(N)
   
      eL = g_eff[ii-2].GetErrorYhigh(i)/rel_stat_unc
      eU = g_eff[ii-2].GetErrorYlow(i)/rel_stat_unc
      print "eL:"+str(eL)+"   eU:"+str(eU)

    print "mass = "+str(h_mjj_HLTpass[0].GetBinCenter(i+1))+"   N2 = "+str(N2)+"  rel_stat_unc = "+str(rel_stat_unc)+"  ineff = "+str(ineff)+"   y = "+str(ineff / rel_stat_unc)+"  eL = "+str(eL)+"  eU = "+str(eU) 
    x.append( h_mjj_HLTpass[0].GetBinCenter(i+1))
    y.append( ineff / rel_stat_unc )
    exl.append(h_mjj_HLTpass[0].GetBinWidth(i+1)/2)
    exh.append(h_mjj_HLTpass[0].GetBinWidth(i+1)/2)
    eyl.append(eL)
    eyh.append(eU)
    
  print "\n"
  vx = array("f",x)
  vy = array("f",y)
  vexl = array("f",exl)
  vexh = array("f",exh)
  veyl = array("f",eyl)
  veyh = array("f",eyh)
  g_ineff.append(TGraphAsymmErrors(nMassBins,vx,vy,vexl,vexh,veyh,veyh) ) 

g_ineff[0].SetName("g_ineff_PFHT800")
g_ineff[1].SetName("g_ineff_PFHT650MJJ900")
g_ineff[2].SetName("g_ineff_PFHT800_OR_MJJ900")


####  Draw plots
gStyle.SetLabelSize(0.04,"XY")
g_ineff[0].GetXaxis().SetRangeUser(1000,2037)
g_ineff[0].GetYaxis().SetRangeUser(0,60)
g_ineff[0].GetXaxis().SetTitle("Dijet Mass [GeV]")
g_ineff[0].GetYaxis().SetTitle("#frac{(1-eff)}{relative stat. uncertainty}")
g_ineff[0].GetYaxis().SetTitleOffset(1.35)
g_ineff[0].GetYaxis().SetTitleSize(0.05)


g_ineff[0].SetLineColor(1)
#g_ineff[1].SetLineColor(2)
#g_ineff[2].SetLineColor(3)
g_ineff[0].SetMarkerColor(1)
#g_ineff[1].SetMarkerColor(2)
#g_ineff[2].SetMarkerColor(3)
g_ineff[0].SetMarkerStyle(20)
#g_ineff[1].SetMarkerStyle(21)
#g_ineff[2].SetMarkerStyle(22)

#l2 = TLegend(0.5,0.4,0.85,0.6)
l3 = TLegend(0.5,0.7,0.85,0.75)
l3.SetFillStyle(0)
l3.AddEntry(g_ineff[0],"PFHT800","PLE")
#l.AddEntry(g_ineff[1],"PFHT650MJJ900","p")
#l.AddEntry(g_ineff[2],"PFHT800_OR_MJJ900","p")

c_ineff = TCanvas("c_ineff","",600,600)
c_ineff.cd()
g_ineff[0].Draw("APL")
#g_ineff[1].Draw("PL SAME")
#g_ineff[2].Draw("PL SAME")
l3.Draw()
#c_ineff.UseCurrentStyle()
c_ineff.RedrawAxis()
#draw the lumi text on the canvas
CMS_lumi.CMS_lumi(c_ineff, iPeriod, iPos)
c_ineff.SaveAs(outputDir+"/trigger_ineff_data.png")
c_ineff.SaveAs(outputDir+"/trigger_ineff_data.pdf")
c_ineff.SaveAs(outputDir+"/trigger_ineff_data.root")

g_ineff[0].GetYaxis().SetRangeUser(0,4)
g_ineff[0].Draw("APL")
c_ineff.RedrawAxis()
line = TLine(c_ineff.GetUxmin(),0.15,c_ineff.GetUxmax(),0.15);
line.SetLineStyle(2)
line.SetLineWidth(2)
line.Draw("")
l3.Draw()
#draw the lumi text on the canvas
CMS_lumi.CMS_lumi(c_ineff, iPeriod, iPos)
c_ineff.SaveAs(outputDir+"/trigger_ineff_data_zoom.png")
c_ineff.SaveAs(outputDir+"/trigger_ineff_data_zoom.pdf")
c_ineff.SaveAs(outputDir+"/trigger_ineff_data_zoom.root")

#turnon.SetLineColor(2)
#turnon.SetLineWidth(2)
g_eff[0].SetMarkerStyle(20)
#g_eff[1].SetMarkerStyle(21)
#g_eff[2].SetMarkerStyle(22)
g_eff[0].SetMarkerColor(1)
#g_eff[1].SetMarkerColor(2)
#g_eff[2].SetMarkerColor(3)
g_eff[0].GetXaxis().SetTitle("Dijet Mass [GeV]")
g_eff[0].GetYaxis().SetTitle("Relative Efficiency")
g_eff[0].GetXaxis().SetRangeUser(minX_mass,maxX_mass)
g_eff[0].GetYaxis().SetRangeUser(0,1.3)

c = TCanvas("c","",600,600)
c.cd()
g_eff[0].Draw("AP")
#fermi.Draw("same")
#turnon.Draw("same")
#g_eff[1].Draw("P SAME")
#g_eff[2].Draw("P SAME")
l = TLegend(0.5,0.5,0.85,0.55)
l.SetFillStyle(0)
l.AddEntry(g_eff[0],"PF H_{T} > 800 GeV","PLE")
#l.AddEntry(fermi,"fit","L")
l.Draw()
#draw the lumi text on the canvas
CMS_lumi.CMS_lumi(c, iPeriod, iPos)

c.SaveAs(outputDir+"/triggerCurves_data.png")
c.SaveAs(outputDir+"/triggerCurves_data.pdf")
c.SaveAs(outputDir+"/triggerCurves_data.gif")

g_eff[0].GetXaxis().SetRangeUser(minX_mass,maxX_mass)
g_eff[0].GetYaxis().SetRangeUser(0.8,1.05)
g_eff[0].Draw("AP")
#g_eff[1].Draw("P SAME")
#g_eff[2].Draw("P SAME")
l = TLegend(0.2,0.50,0.55,0.55)
l.AddEntry(g_eff[0],"PFHT800","PLE")
l.Draw()
#draw the lumi text on the canvas
CMS_lumi.CMS_lumi(c, iPeriod, iPos)
c.SaveAs(outputDir+"/triggerCurves_data_zoom.png")
c.SaveAs(outputDir+"/triggerCurves_data_zoom.pdf")


h_mjj_HLTpass_PFHT475_all.SetLineColor(1)  
h_mjj_HLTpass_PFHT800_all.SetLineColor(2)    
h_mjj_HLTpass_PFHT650MJJ900_all.SetLineColor(3)    
h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_all.SetLineColor(4)    
h_mjj_HLTpass_PFHT475_all.GetXaxis().SetTitle("mjj [GeV]")  
h_mjj_HLTpass_PFHT475_all.GetYaxis().SetTitle("events")  
h_mjj_HLTpass_PFHT475_all.GetXaxis().SetRangeUser(minX_mass,maxX_mass)  
h_mjj_HLTpass_PFHT475_all.GetXaxis().SetNdivisions(404)  
h_mjj_HLTpass_PFHT475_all.UseCurrentStyle()  
c2 = TCanvas("c2","",600,600)
c2.cd()
c2.SetLogy(1)
h_mjj_HLTpass_PFHT475_all.Draw()
h_mjj_HLTpass_PFHT800_all.Draw("same")
#h_mjj_HLTpass_PFHT650MJJ900_all.Draw("same")
#h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_all.Draw("same")

l2 = TLegend(0.5,0.75,0.85,0.85)
#l2 = TLegend(0.5,0.6,0.85,0.85)
l2.AddEntry(h_mjj_HLTpass_PFHT475_all,"PFHT475","l")
l2.AddEntry(h_mjj_HLTpass_PFHT800_all,"PFHT800","l")
#l2.AddEntry(h_mjj_HLTpass_PFHT650MJJ900_all,"PFHT650MJJ900","l")
#l2.AddEntry(h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_all,"PFHT800_OR_PFHT650MJJ900","l")
l2.Draw()
#draw the lumi text on the canvas
CMS_lumi.CMS_lumi(c2, iPeriod, iPos)

c2.SaveAs(outputDir+"/triggerPass_data.png")
c2.SaveAs(outputDir+"/triggerPass_data.pdf")


file_out = TFile(outputDir+"/triggerCurves_data.root","recreate")
g_eff[0].Write()
g_eff[1].Write()
g_eff[2].Write()
g_ineff[0].Write()
g_ineff[1].Write()
g_ineff[2].Write()
h_mjj_HLTpass_PFHT475_all.Write()
h_mjj_HLTpass_PFHT800_all.Write()
h_mjj_HLTpass_PFHT650MJJ900_all.Write()
h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_all.Write()
h_mjj_HLTpass_PFHT800_noPFHT475_all.Write()

c_ineff.Write()
c.Write() 
c2.Write()
