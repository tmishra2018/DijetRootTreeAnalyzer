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

CMS_lumi.extraText = "Simulation Preliminary"
CMS_lumi.lumi_sqrtS = str(LUMI)+" fb^{-1} (13 TeV)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
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
  xsecs.append(parts[1])
  print ("dataset : %s    xsec : %s" % (fileNames[ii], xsecs[ii]))


  ii+=1

#---- open the files --------------------
h_mjj_HLTpass_noTrig_list = []
h_mjj_HLTpass_PFHT475_list = []
h_mjj_HLTpass_PFHT800_list = []
h_mjj_HLTpass_PFHT650MJJ900_list = []
h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_list = []

print fileNames
filelist =[]
for f in fileNames:
  inf = TFile(f)
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

  h_mjj_HLTpass_noTrig_list.append(h_mjj_HLTpass_noTrig)
  h_mjj_HLTpass_PFHT475_list.append(h_mjj_HLTpass_PFHT475 )
  h_mjj_HLTpass_PFHT800_list.append(h_mjj_HLTpass_PFHT800)
  h_mjj_HLTpass_PFHT650MJJ900_list.append(h_mjj_HLTpass_PFHT650MJJ900)
  h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_list.append(h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900)
  #running on data
  #wt = LUMI*float(xsecs[i_f])/Nev
  print('weight : %f' % wt)
  h_mjj_HLTpass_noTrig_list[i_f].Scale(wt)
  h_mjj_HLTpass_PFHT475_list[i_f].Scale(wt)
  h_mjj_HLTpass_PFHT800_list[i_f].Scale(wt) 
  h_mjj_HLTpass_PFHT650MJJ900_list[i_f].Scale(wt)
  h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_list[i_f].Scale(wt)
  h_mjj_HLTpass_PFHT475_list[i_f].Print()
  print str(i_f) 
  #h_mjj_HLTpass_PFHT475_list[0].Print()

  i_f += 1

#h_mjj_HLTpass_PFHT475_list[0].Print()
#print h_mjj_HLTpass_PFHT475_list

h_mjj_HLTpass_noTrig_all = h_mjj_HLTpass_PFHT475_list[0].Clone("h_mjj_HLTpass_noTrig_all")
h_mjj_HLTpass_PFHT475_all = h_mjj_HLTpass_PFHT475_list[0].Clone("h_mjj_HLTpass_PFHT475_all")
h_mjj_HLTpass_PFHT800_all = h_mjj_HLTpass_PFHT800_list[0].Clone("h_mjj_HLTpass_PFHT800_all")
h_mjj_HLTpass_PFHT650MJJ900_all = h_mjj_HLTpass_PFHT650MJJ900_list[0].Clone("h_mjj_HLTpass_PFHT650MJJ900_all")
h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_all = h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_list[0].Clone("h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_all")

if(i_f>1):
  for i in range(1,9):
    h_mjj_HLTpass_noTrig_all.Add(h_mjj_HLTpass_noTrig_list[i])
    h_mjj_HLTpass_PFHT475_all.Add(h_mjj_HLTpass_PFHT475_list[i])
    h_mjj_HLTpass_PFHT800_all.Add(h_mjj_HLTpass_PFHT800_list[i])
    h_mjj_HLTpass_PFHT650MJJ900_all.Add(h_mjj_HLTpass_PFHT650MJJ900_list[i])
    h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_all.Add(h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_list[i])

h_mjj_HLTpass=[h_mjj_HLTpass_noTrig_all,h_mjj_HLTpass_PFHT475_all,h_mjj_HLTpass_PFHT800_all, h_mjj_HLTpass_PFHT650MJJ900_all, h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_all]

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
    N1 = h_mjj_HLTpass[0].GetBinContent(i+1)
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
    print str(N1)+" "+str(N2)+" "+str(p)+" "+str(eL)+" "+str(eU)

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
###############   

g_ineff = []
### ineff / relative_stat_unc
for ii in range(2,5):
  x = []
  y = []
  ex = []
  ey = []
  for i in range(0,nMassBins):
    N1 = h_mjj_HLTpass[0].GetBinContent(i+1)
    N2 = h_mjj_HLTpass[ii].GetBinContent(i+1)
    ineff = 1
    p = 0
    rel_stat_unc = 0.000001    
    if (N1>0 and N2>0) :
      p = N2/N1
      ineff = 1-p
      rel_stat_unc = 1/math.sqrt(N2)
   
    print "N2 = "+str(N2)+"  rel_stat_unc = "+str(rel_stat_unc)+"  ineff = "+str(ineff)+"   y = "+str(ineff / rel_stat_unc) 
    x.append( h_mjj_HLTpass[0].GetBinCenter(i+1))
    y.append( ineff / rel_stat_unc )
    ex.append(h_mjj_HLTpass[0].GetBinWidth(i+1)/2)
    ey.append(0)

  print "\n"
  vx = array("f",x)
  vy = array("f",y)
  vex = array("f",ex)
  vey = array("f",ey)
  g_ineff.append(TGraphErrors(nMassBins,vx,vy,vex,vey) ) 



####  Draw plots
gStyle.SetLabelSize(0.04,"XY")
g_ineff[0].GetXaxis().SetRangeUser(1000,2037)
g_ineff[0].GetYaxis().SetRangeUser(0,4)
g_ineff[0].GetXaxis().SetTitle("mjj [GeV]")
g_ineff[0].GetYaxis().SetTitle("#frac{(1-eff)}{relative stat. uncertainty}")
g_ineff[0].GetYaxis().SetTitleSize(0.02)
g_ineff[0].GetYaxis().SetTitleOffset(0.05)
g_ineff[0].UseCurrentStyle()

g_ineff[0].SetLineColor(1)
g_ineff[1].SetLineColor(2)
g_ineff[2].SetLineColor(3)
g_ineff[0].SetMarkerColor(1)
g_ineff[1].SetMarkerColor(2)
g_ineff[2].SetMarkerColor(3)
g_ineff[0].SetMarkerStyle(20)
g_ineff[1].SetMarkerStyle(21)
g_ineff[2].SetMarkerStyle(22)

l = TLegend(0.5,0.4,0.85,0.6)
l.AddEntry(g_ineff[0],"PFHT800","p")
l.AddEntry(g_ineff[1],"PFHT650MJJ900","p")
l.AddEntry(g_ineff[2],"PFHT800_OR_MJJ900","p")

c_ineff = TCanvas("c_ineff","",600,600)
c_ineff.cd()
g_ineff[0].Draw("APL")
g_ineff[1].Draw("PL SAME")
g_ineff[2].Draw("PL SAME")
l.Draw()
c_ineff.SaveAs(outputDir+"/trigger_ineff_MC.png")

g_eff[0].SetMarkerStyle(20)
g_eff[1].SetMarkerStyle(21)
g_eff[2].SetMarkerStyle(22)
g_eff[0].SetMarkerColor(1)
g_eff[1].SetMarkerColor(2)
g_eff[2].SetMarkerColor(3)
g_eff[0].GetXaxis().SetTitle("mjj [GeV]")
g_eff[0].GetYaxis().SetTitle("efficiency")
g_eff[0].GetXaxis().SetRangeUser(526,6099)

c = TCanvas("c","",600,600)
c.cd()
g_eff[0].Draw("AP")
g_eff[1].Draw("P SAME")
g_eff[2].Draw("P SAME")
l.Draw()

c.SaveAs(outputDir+"/triggerCurves_MC.png")

g_eff[0].GetXaxis().SetRangeUser(526,2037)
g_eff[0].GetYaxis().SetRangeUser(0.8,1.02)
g_eff[0].Draw("AP")
g_eff[1].Draw("P SAME")
g_eff[2].Draw("P SAME")
l.Draw()
c.SaveAs(outputDir+"/triggerCurves_MC_zoom.png")


h_mjj_HLTpass_PFHT475_all.SetLineColor(1)  
h_mjj_HLTpass_PFHT800_all.SetLineColor(2)    
h_mjj_HLTpass_PFHT650MJJ900_all.SetLineColor(3)    
h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_all.SetLineColor(4)    
h_mjj_HLTpass_PFHT475_all.GetXaxis().SetTitle("mjj [GeV]")  
h_mjj_HLTpass_PFHT475_all.GetYaxis().SetTitle("events")  
h_mjj_HLTpass_PFHT475_all.GetXaxis().SetRangeUser(526,7060)  

c2 = TCanvas("c2","",600,600)
c2.cd()
c2.SetLogy(1)
h_mjj_HLTpass_PFHT475_all.Draw()
h_mjj_HLTpass_PFHT800_all.Draw("same")
h_mjj_HLTpass_PFHT650MJJ900_all.Draw("same")
h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_all.Draw("same")

l2 = TLegend(0.5,0.6,0.85,0.85)
l2.AddEntry(h_mjj_HLTpass_PFHT475_all,"PFHT475","l")
l2.AddEntry(h_mjj_HLTpass_PFHT800_all,"PFHT800","l")
l2.AddEntry(h_mjj_HLTpass_PFHT650MJJ900_all,"PFHT650MJJ900","l")
l2.AddEntry(h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_all,"PFHT800_OR_PFHT650MJJ900","l")
l2.Draw()

c2.SaveAs(outputDir+"/triggerPass_MC.png")


file_out = TFile(outputDir+"/triggerCurves_MC.root","recreate")
g_eff[0].Write()
g_eff[1].Write()
g_eff[2].Write()
h_mjj_HLTpass_PFHT475_all.Write()
h_mjj_HLTpass_PFHT800_all.Write()
h_mjj_HLTpass_PFHT650MJJ900_all.Write()
h_mjj_HLTpass_PFHT800_OR_PFHT650MJJ900_all.Write()

c_ineff.Write()
c.Write() 
c2.Write()
