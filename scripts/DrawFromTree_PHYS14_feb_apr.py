#!usr/bin/python

from setTDRStyle import setTDRStyle
import sys, os, subprocess, string, re
from ROOT import *
from array import array
import CMS_lumi

import optparse

CMS_lumi.extraText = "Simulation Preliminary"
CMS_lumi.lumi_sqrtS = "1 fb^{-1} (13 TeV)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 11
iPeriod = 0

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

#######################################################

usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("--var",action="store",type="string",dest="var",default='ptHat')
parser.add_option("--xmin",action="store",type="float",dest="xmin",default=1)
parser.add_option("--xmax",action="store",type="float",dest="xmax",default=1)
parser.add_option("--xtitle",action="store",type="string",dest="xtitle",default='')
parser.add_option("--bins",action="store",type="int",dest="bins",default=11111111111)
parser.add_option("--rebin",action="store",type="int",dest="rebin",default=1)
parser.add_option("--logy",action="store_true",default=False,dest="logy")
parser.add_option("--outputDir",action="store",type="string",default="./",dest="outputDir")
parser.add_option("--inputList",action="store",type="string",default="list.txt",dest="inputList")
parser.add_option("--inputList_PHYS14_feb15",action="store",type="string",default="list.txt",dest="inputList_PHYS14_feb15")
parser.add_option("--lumi",action="store",type="float",default="1000.",dest="lumi")

(options, args) = parser.parse_args()

var = options.var
xmin = options.xmin
xmax = options.xmax
bins = options.bins
xtitle = options.xtitle
rebin = options.rebin
logy = options.logy
outputDir = options.outputDir
inputList = options.inputList
inputList_PHYS14_feb15 = options.inputList_PHYS14_feb15
lumi = options.lumi

gROOT.Reset()
setTDRStyle()
gROOT.ForceStyle()
gROOT.SetStyle('tdrStyle')

#fileNames = ['QCD_Pt-301to470','QCD_Pt-470to600','QCD_Pt-600to800', 'QCD_Pt-800to1000', 'QCD_Pt-1000to1400', 'QCD_Pt-1400to1800', 'QCD_Pt-1800to2400', 'QCD_Pt-2400to3200', 'QCD_Pt-3200']
#xsections = [7475 ,587., 167, 28.25, 8.195, 0.7346, 0.102, 0.00644, 0.000163]
#colorF    = [ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9,ROOT.kBlue-9]
#colorL    = [ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack,ROOT.kBlack]
hist_allCuts      = []
hist_allCuts_PHYS14_feb15 = []

LUMI      = lumi
#PATH      = inputDir

#---- read the list -----------------
lines = [line.strip() for line in open(inputList)]
lines_PHYS14_feb15 = [line_PHYS14_feb15.strip() for line_PHYS14_feb15 in open(inputList_PHYS14_feb15)]

#---- split sample name and xsec
fileNames = []
fileNames_PHYS14_feb15 = []
xsecs = []
ii = 0
for line in lines:
  parts = line.split()
  fileNames.append(parts[0])
  xsecs.append(parts[1])
  print ("dataset : %s    xsec : %s" % (fileNames[ii], xsecs[ii]))
  ii+=1

ii = 0
for line_PHYS14_feb15 in lines_PHYS14_feb15:
  parts_PHYS14_feb15 = line_PHYS14_feb15.split()
  fileNames_PHYS14_feb15.append(parts_PHYS14_feb15[0])
  xsecs.append(parts[1])
  print ("dataset PHYS14_feb15 : %s    xsec : %s" % (fileNames_PHYS14_feb15[ii], xsecs[ii]))


  ii+=1


#---- open the files --------------------
i_f = 0
for f in fileNames:
  inf = TFile.Open(f)
  print inf.GetName()
  
  Nev = inf.Get('DijetFilter/EventCount/EventCounter').GetBinContent(1)
  print ('processed events: %s' % Nev)
  wt = 1.0
  #if i_f < 3:

  h_allCuts = TH1F("h_allCuts", "", bins, xmin, xmax)
  h_allCuts.Sumw2()
  tree = inf.Get('rootTupleTree/tree')
  tree.Project(h_allCuts.GetName(), var,' deltaETAjj < 1.3')
  Npassed = h_allCuts.GetEntries()
  eff = float(Npassed)/Nev
  print('eff : %f' % eff)
  wt = LUMI*float(xsecs[i_f])*eff/Nev
  print('weight : %f' % wt)
  h_allCuts.Scale(wt)
  h_allCuts.Rebin(rebin)
  h_allCuts.SetDirectory(0)
  h_allCuts.SetFillColor(ROOT.kBlue-9)
  h_allCuts.SetLineColor(ROOT.kBlue-9)
  h_allCuts.SetMarkerColor(ROOT.kBlue-9)
  hist_allCuts.append(h_allCuts)
   
  i_f += 1

#------- open the files PHYS14_feb15 -------------
i_f = 0
for f in fileNames_PHYS14_feb15:
  inf = TFile.Open(f)
  print inf.GetName()
  
  Nev = inf.Get('DijetFilter/EventCount/EventCounter').GetBinContent(1)
  print ('processed events PHYS14_feb15 : %s' % Nev)
  wt = 1.0
  #if i_f < 3:

  h_allCuts_PHYS14_feb15 = TH1F("h_allCuts_PHYS14_feb15", "", bins, xmin, xmax)
  h_allCuts_PHYS14_feb15.Sumw2()
  tree = inf.Get('rootTupleTree/tree')
  tree.Project(h_allCuts_PHYS14_feb15.GetName(), var,' deltaETAjj < 1.3')
  Npassed = h_allCuts_PHYS14_feb15.GetEntries()
  eff = float(Npassed)/Nev
  print('eff : %f' % eff)
  wt = LUMI*float(xsecs[i_f])*eff/Nev
  print('weight : %f' % wt)
  h_allCuts_PHYS14_feb15.Scale(wt)
  h_allCuts_PHYS14_feb15.Rebin(rebin)
  h_allCuts_PHYS14_feb15.SetDirectory(0)
  h_allCuts_PHYS14_feb15.SetLineColor(ROOT.kRed)
  h_allCuts_PHYS14_feb15.SetMarkerColor(ROOT.kRed)
  hist_allCuts_PHYS14_feb15.append(h_allCuts_PHYS14_feb15)
   
  i_f += 1


#kFactor = NDAT/NQCD
kFactor = 1.3
print ("kFactor = %f" % kFactor)

                                      
for i in range(0,len(fileNames)) :
  hist_allCuts[i].Scale(kFactor)
  hist_allCuts_PHYS14_feb15[i].Scale(kFactor)


NQCD_allCuts = hist_allCuts[0].Integral()
NQCD_allCuts_PHYS14_feb15 = hist_allCuts_PHYS14_feb15[i].Integral()

for i in range(0,len(fileNames)) :
  NQCD_allCuts += hist_allCuts[i].Integral()
  NQCD_allCuts_PHYS14_feb15 += hist_allCuts_PHYS14_feb15[i].Integral()
  

hist_allCutsQCD = hist_allCuts[0].Clone('hist_allCutsQCD')
hist_allCutsQCD_PHYS14_feb15 = hist_allCuts_PHYS14_feb15[0].Clone('hist_allCutsQCD_PHYS14_feb15')

for i in range(1,len(fileNames)):
  hist_allCutsQCD.Add(hist_allCuts[i])
  hist_allCutsQCD_PHYS14_feb15.Add(hist_allCuts_PHYS14_feb15[i])

#hsQCD_allCuts = THStack('QCD_allCuts','QCD_allCuts')

#for i in range(0,len(fileNames)) :
#  hsQCD_allCuts.Add(hist_allCuts[i])


print ("---- After scaling signal to bkg (if not plotting mjj) -----")
print ("bkg integral all cuts = %f" % NQCD_allCuts)
print ("bkg integral all cuts PHYS14_feb15 = %f" % NQCD_allCuts_PHYS14_feb15)


#----- Drawing  and save on file -----------------------

outFile = TFile(outputDir+"histo_signal_bkg_"+var+".root", "recreate")
outFile.cd()
hist_allCutsQCD.Write()
hist_allCutsQCD_PHYS14_feb15.Write()

can_allCuts = TCanvas('can_allCuts_'+var,'can_allCuts_'+var,900,600)

leg = TLegend(0.6, 0.7, 0.85, 0.85)
leg.SetLineColor(0)
leg.SetFillColor(0)
leg.SetBorderSize(0)
leg.AddEntry(hist_allCutsQCD, "QCD PHYS14", "f")
leg.AddEntry(hist_allCutsQCD_PHYS14_feb15, "QCD PHYS14_feb15", "l")

can_allCuts.cd()

#----- pad 1 -----------
pad1 = TPad("pad1", "pad1",0.01,0.13,0.75,1.)  
pad1.SetRightMargin(0.1)

pad1.SetLogy()
pad1.Draw()
pad1.cd()
pad1.Clear()
     
if logy:
  pad1.SetLogy()
  max=hist_allCutsQCD.GetBinContent(hist_allCutsQCD.GetMaximumBin())
  hist_allCutsQCD.SetMaximum(500*max)

#hist_allCutsQCD.Reset()
hist_allCutsQCD.GetXaxis().SetRangeUser(xmin,xmax)
hist_allCutsQCD.GetXaxis().SetTitle(xtitle)
hist_allCutsQCD.GetXaxis().SetTitleFont(42)
hist_allCutsQCD.GetXaxis().SetTitleSize(0.05)
#hist_allCutsQCD.GetXaxis().SetLabelOffset(999)
hist_allCutsQCD.GetXaxis().SetLabelSize(0)
size = hist_allCutsQCD.GetBinWidth(1)
title_y = "events / "+str(size)+" GeV"
hist_allCutsQCD.GetYaxis().SetTitleFont(42)
hist_allCutsQCD.GetYaxis().SetTitleSize(0.04)
hist_allCutsQCD.GetYaxis().SetTitle(title_y)

#maximumBin = array('f',  [hist_allCutsQCD.GetBinContent(hist_allCutsQCD.GetMaximumBin()), hist_allCutsSig_1000.GetBinContent(hist_allCutsSig_1000.GetMaximumBin()), hist_allCutsSig_5000.GetBinContent(hist_allCutsSig_5000.GetMaximumBin()), hist_allCutsSig_8000.GetBinContent(hist_allCutsSig_8000.GetMaximumBin())])
#max = TMath.MaxElement(4, maximumBin)
#hist_allCutsQCD.SetMaximum(1.2*max)
#hist_allCutsQCD.SetMinimum(0.00001)
#hist_allCutsQCD.Rebin()
#hist_allCutsQCD_PHYS14_feb15.Rebin()
hist_allCutsQCD.Draw("hist")
hist_allCutsQCD_PHYS14_feb15.Draw("hist same")
leg.Draw()
#gPad.RedrawAxis()
#draw the lumi text on the canvas
CMS_lumi.CMS_lumi(pad1, iPeriod, iPos)

can_allCuts.cd()
#-------pad 2------
pad2 = TPad("pad2", "pad2",0.01,0.001,0.75,0.222)
pad2.SetGrid()
	      
pad2.SetTopMargin(0)
pad2.SetBottomMargin(0.4)
pad2.SetRightMargin(0.1)
pad2.Draw()	       
pad2.cd()

ratio = hist_allCutsQCD.Clone("ratio")
ratio.Divide(hist_allCutsQCD_PHYS14_feb15)
ratio.SetFillColor(0)
ratio.SetLineColor(ROOT.kBlack)
ratio.SetMarkerColor(ROOT.kBlack)
ratio.Draw("p")
ratio.GetYaxis().SetRangeUser(0., 2.)
ratio.GetYaxis().SetNdivisions(4, ROOT.kTRUE)
ratio.GetYaxis().SetTitleFont(42)
ratio.GetYaxis().SetTitle("#frac{PHYS14}{PHYS14_feb15}")
ratio.GetXaxis().SetTitleSize(0.2)
ratio.GetXaxis().SetLabelSize(0.16)
ratio.GetYaxis().SetLabelSize(0.16)
ratio.GetYaxis().SetTitleSize(0.15)
#ratio.GetYaxis().SetTitleOffset(0.65)
ratio.GetXaxis().SetTitleOffset(0.8)
ratio.Draw()

can_allCuts.Write()   
can_allCuts.SaveAs(outputDir+var+'_allCuts.png')
can_allCuts.SaveAs(outputDir+var+'_allCuts.svg')
can_allCuts.Close()

outFile.Close()

##----- keep the GUI alive ------------
#if __name__ == '__main__':
#  rep = ''
#  while not rep in ['q','Q']:
#    rep = raw_input('enter "q" to quit: ')
#    if 1 < len(rep):
#      rep = rep[0]
