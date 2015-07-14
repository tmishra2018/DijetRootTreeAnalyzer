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
parser.add_option("--bins",action="store",type="int",dest="bins",default=100)
parser.add_option("--rebin",action="store",type="int",dest="rebin",default=1)
parser.add_option("--logy",action="store_true",default=False,dest="logy")
parser.add_option("--outputDir",action="store",type="string",default="./",dest="outputDir")
parser.add_option("--inputList_1",action="store",type="string",default="list1.txt",dest="inputList_1")
parser.add_option("--inputList_2",action="store",type="string",default="list2.txt",dest="inputList_2")
parser.add_option("--name1",action="store",type="string",default="name1",dest="name1")
parser.add_option("--name2",action="store",type="string",default="name2",dest="name2")
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
inputList_1 = options.inputList_1
inputList_2 = options.inputList_2
name1 = options.name1
name2 = options.name2
lumi = options.lumi

gROOT.Reset()
setTDRStyle()
gROOT.ForceStyle()
gROOT.SetStyle('tdrStyle')
#####################

minX_mass = 1118
maxX_mass = 6099 

massBins_list = [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000]
    
massBins = array("d",massBins_list)


#fileNames = ['QCD_Pt-301to470','QCD_Pt-470to600','QCD_Pt-600to800', 'QCD_Pt-800to1000', 'QCD_Pt-1000to1400', 'QCD_Pt-1400to1800', 'QCD_Pt-1800to2400', 'QCD_Pt-2400to3200', 'QCD_Pt-3200']
#xsections = [7475 ,587., 167, 28.25, 8.195, 0.7346, 0.102, 0.00644, 0.000163]
#colorF    = [ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9,ROOT.kBlue-9]
#colorL    = [ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack,ROOT.kBlack]
hist_allCuts      = []
hist_allCuts_2 = []

LUMI      = lumi
#PATH      = inputDir

#---- read the list -----------------
lines1 = [line1.strip() for line1 in open(inputList_1)]
lines2 = [line2.strip() for line2 in open(inputList_2)]

#---- split sample name and xsec
fileNames = []
fileNames2 = []
xsecs = []
ii = 0
for line in lines1:
  parts = line.split()
  fileNames.append(parts[0])
  xsecs.append(parts[1])
  print ("dataset 1 : %s    xsec : %s" % (fileNames[ii], xsecs[ii]))
  ii+=1

ii = 0
for line in lines2:
  parts2 = line.split()
  fileNames2.append(parts2[0])
  xsecs.append(parts[1])
  print ("dataset 2 : %s    xsec : %s" % (fileNames2[ii], xsecs[ii]))


  ii+=1


#---- open the files --------------------
var1 = ""
if var=="pTWJ_j1" : var1 = "pT_j1"
elif var=="pTWJ_j2" : var1 = "pT_j2"
elif var=="etaWJ_j1" : var1 = "eta_j2"
elif var=="etaWJ_j2" : var1 = "eta_j2"
elif var=="phiWJ_j1" : var1 = "phi_j1"
elif var=="phiWJ_j2" : var1 = "phi_j2"
else : var1 = var

dataset1 = []
dataset2 = []

i_f = 0
for f in fileNames:
  inf = TFile.Open(f)
  print inf.GetName()
  dataset1.append ( os.path.basename(fileNames[i_f]) )
  dataset1[i_f] = dataset1[i_f].split("rootfile_")[1]
  dataset1[i_f] = dataset1[i_f].split("_reduced_skim.root")[0]
  dataset1[i_f] = dataset1[i_f].split("__")[0]

  Nev = inf.Get('DijetFilter/EventCount/EventCounter').GetBinContent(1)
  print ('processed events list1 : %s' % Nev)
  wt = 1.0
  #if i_f < 3:

  h_allCuts = TH1F("h_allCuts", "", bins, xmin, xmax)
  h_allCuts.Sumw2()
  tree = inf.Get('rootTupleTree/tree')
  tree.Project(h_allCuts.GetName(), var1,' deltaETAjj < 1.3')
  Npassed = h_allCuts.GetEntries()
  eff = float(Npassed)/Nev
  print('eff : %f' % eff)
  wt = LUMI*float(xsecs[i_f])/Nev
  print('weight : %f' % wt)
  h_allCuts.Scale(wt)
  h_allCuts.SetDirectory(0)
  h_allCuts.SetFillColor(kBlue-9)
  h_allCuts.SetLineColor(kBlue-9)
  h_allCuts.SetMarkerColor(kBlue-9)
  hist_allCuts.append(h_allCuts)
   
  i_f += 1

#------- open the files in 2nd list -------------

i_f = 0
for f in fileNames2:
  inf = TFile.Open(f)
  print inf.GetName()
  dataset2.append ( os.path.basename(fileNames[i_f]) )
  dataset2[i_f] = dataset2[i_f].split("rootfile_")[1]
  dataset2[i_f] = dataset2[i_f].split("_reduced_skim.root")[0]
  dataset2[i_f] = dataset2[i_f].split("__")[0]
  
  Nev = inf.Get('DijetFilter/EventCount/EventCounter').GetBinContent(1)
  print ('processed events list 2: %s' % Nev)
  wt = 1.0
  #if i_f < 3:

  h_allCuts_2 = TH1F("h_allCuts_2", "", bins, xmin, xmax)
  h_allCuts_2.Sumw2()
  tree = inf.Get('rootTupleTree/tree')
  tree.Project(h_allCuts_2.GetName(), var,' deltaETAjj < 1.3')
  Npassed = h_allCuts_2.GetEntries()
  eff = float(Npassed)/Nev
  print('eff : %f' % eff)
  wt = LUMI*float(xsecs[i_f])/Nev
  print('weight : %f' % wt)
  h_allCuts_2.Scale(wt)
  h_allCuts_2.SetDirectory(0)
  h_allCuts_2.SetLineColor(kRed)
  h_allCuts_2.SetMarkerColor(kRed)
  hist_allCuts_2.append(h_allCuts_2)
   
  i_f += 1


#kFactor = NDAT/NQCD
kFactor = 1.
print ("kFactor = %f" % kFactor)

                                      
for i in range(0,9) :
  hist_allCuts[i].Scale(kFactor)
  hist_allCuts_2[i].Scale(kFactor)


NQCD_allCuts = hist_allCuts[0].Integral()
NQCD_allCuts_2 = hist_allCuts_2[0].Integral()

for i in range(0,9) :
  NQCD_allCuts += hist_allCuts[i].Integral()
  NQCD_allCuts_2 += hist_allCuts_2[i].Integral()
  

hist_allCutsQCD = hist_allCuts[0].Clone(name1)
hist_allCutsQCD_2 = hist_allCuts_2[0].Clone(name2)

for i in range(1,9):
  hist_allCutsQCD.Add(hist_allCuts[i])
  hist_allCutsQCD_2.Add(hist_allCuts_2[i])

#hsQCD_allCuts = THStack('QCD_allCuts','QCD_allCuts')

#for i in range(0,len(fileNames)) :
#  hsQCD_allCuts.Add(hist_allCuts[i])

for i in range(9,len(fileNames)) :
  hist_allCuts[i].Scale(1./hist_allCuts[i].Integral())
  hist_allCuts[i].SetName("h_"+dataset1[i])
  hist_allCuts_2[i].Scale(1./hist_allCuts_2[i].Integral())
  hist_allCuts_2[i].SetName("h_"+dataset2[i])


##print ("---- After scaling signal to bkg (if not plotting mjj) -----")
print ("bkg integral all cuts = %f" % NQCD_allCuts)
print ("bkg integral all cuts 2 = %f" % NQCD_allCuts_2)


#----- Drawing  and save on file -----------------------

outFile = TFile(outputDir+"histo_signal_bkg_"+var+".root", "recreate")
outFile.cd()
hist_allCutsQCD.Write()
hist_allCutsQCD_2.Write()
for i in range(9, len(fileNames)):
  hist_allCuts[i].Write()
  hist_allCuts_2[i].Write()


can_allCuts = TCanvas('can_allCuts_'+var,'can_allCuts_'+var,900,600)

leg = TLegend(0.6, 0.7, 0.85, 0.85)
leg.SetLineColor(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.AddEntry(hist_allCutsQCD, name1, "f")
leg.AddEntry(hist_allCutsQCD_2, name2, "l")

can_allCuts.cd()

#----- pad 1 -----------
pad1 = TPad("pad1", "pad1",0.01,0.13,0.75,1.)  
pad1.SetRightMargin(0.1)

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

###
title_y = "arb. units"
for i in range(9,len(fileNames)):
  hist_allCuts[i].GetXaxis().SetRangeUser(xmin,xmax)
  hist_allCuts[i].GetXaxis().SetTitle(xtitle)
  hist_allCuts[i].GetXaxis().SetTitleFont(42)
  hist_allCuts[i].GetXaxis().SetTitleSize(0.05)
  hist_allCuts[i].GetXaxis().SetLabelSize(0)
  hist_allCuts[i].GetYaxis().SetTitleFont(42)
  hist_allCuts[i].GetYaxis().SetTitleSize(0.04)
  hist_allCuts[i].GetYaxis().SetTitle(title_y)


#maximumBin = array('f',  [hist_allCutsQCD.GetBinContent(hist_allCutsQCD.GetMaximumBin()), hist_allCutsSig_1000.GetBinContent(hist_allCutsSig_1000.GetMaximumBin()), hist_allCutsSig_5000.GetBinContent(hist_allCutsSig_5000.GetMaximumBin()), hist_allCutsSig_8000.GetBinContent(hist_allCutsSig_8000.GetMaximumBin())])
#max = TMath.MaxElement(4, maximumBin)
#hist_allCutsQCD.SetMaximum(1.2*max)
#hist_allCutsQCD.SetMinimum(0.00001)
#hist_allCutsQCD.Rebin()
#hist_allCutsQCD_PHYS14_feb15.Rebin()

hist_allCuts_rebin =[]
hist_allCuts_2_rebin =[]

if var=="mjj":
  hist_allCutsQCD_rebin = hist_allCutsQCD.Rebin(len(massBins_list)-1,name1+"_rebin",massBins)
  hist_allCutsQCD_2_rebin = hist_allCutsQCD_2.Rebin(len(massBins_list)-1,name2+"_rebin",massBins)
  hist_allCutsQCD_rebin.GetXaxis().SetRangeUser(minX_mass,maxX_mass)  
  hist_allCutsQCD_2_rebin.GetXaxis().SetRangeUser(minX_mass,maxX_mass)  
  #hist_allCutsQCD_rebin.GetYaxis().SetRangeUser(0.1,1000000)
  #hist_allCutsQCD_2_rebin.GetYaxis().SetRangeUser(0.1,1000000)
  for i in range(9,len(fileNames)):
    hist_allCuts_rebin.append( hist_allCuts[i].Rebin(len(massBins_list)-1,dataset1[i]+"_"+name1+"_rebin",massBins) )
    hist_allCuts_2_rebin.append( hist_allCuts_2[i].Rebin(len(massBins_list)-1,dataset2[i]+"_"+name2+"_rebin",massBins) )
    hist_allCuts_rebin[i-9].GetXaxis().SetRangeUser(1,10430)
    hist_allCuts_2_rebin[i-9].GetXaxis().SetRangeUser(1,10430)

else :
  hist_allCutsQCD.Rebin(rebin)
  hist_allCutsQCD_2.Rebin(rebin)
  hist_allCutsQCD_rebin = hist_allCutsQCD.Clone(name1+"_rebin")
  hist_allCutsQCD_2_rebin = hist_allCutsQCD_2.Clone(name2+"_rebin")
  for i in range(9,len(fileNames)):
    hist_allCuts[i].Rebin(rebin)
    hist_allCuts_2[i].Rebin(rebin)
    hist_allCuts_rebin.append(hist_allCuts[i].Clone(dataset1[i]+"_"+name1+"_rebin")) 
    hist_allCuts_2_rebin.append(hist_allCuts_2[i].Clone(dataset2[i]+"_"+name2+"_rebin"))

hist_allCutsQCD_rebin.Draw("hist")
hist_allCutsQCD_2_rebin.Draw("hist same")
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

ratio = hist_allCutsQCD_rebin.Clone("ratio")
ratio.Divide(hist_allCutsQCD_2_rebin)
ratio.SetFillColor(0)
ratio.SetLineColor(kBlack)
ratio.SetMarkerColor(kBlack)
ratio.Draw("p")
ratio.GetYaxis().SetRangeUser(0., 2.)
ratio.GetYaxis().SetNdivisions(4, kTRUE)
ratio.GetYaxis().SetTitleFont(42)
ratio.GetYaxis().SetTitle("#frac{"+name1+"}{"+name2+"}")
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

######
# do the same for all signals
####
for i in range(9, len(fileNames)):

  pad1.cd()
  pad1.Clear()
  hist_allCuts_rebin[i-9].Draw("hist")
  hist_allCuts_2_rebin[i-9].Draw("hist same")
  leg.Draw()
  #gPad.RedrawAxis()
  #draw the lumi text on the canvas
  #CMS_lumi.CMS_lumi(pad1, iPeriod, iPos)

  can_allCuts.cd()
  #-------pad 2------
  pad2.cd()
  pad2.Clear()
  ratio = hist_allCuts_rebin[i-9].Clone("ratio")
  ratio.Divide(hist_allCuts_2_rebin[i-9])
  ratio.SetFillColor(0)
  ratio.SetLineColor(kBlack)
  ratio.SetMarkerColor(kBlack)
  ratio.Draw("p")
  ratio.GetYaxis().SetRangeUser(0., 2.)
  ratio.GetYaxis().SetNdivisions(4, kTRUE)
  ratio.GetYaxis().SetTitleFont(42)
  ratio.GetYaxis().SetTitle("#frac{"+name1+"}{"+name2+"}")
  ratio.GetXaxis().SetTitleSize(0.2)
  ratio.GetXaxis().SetLabelSize(0.16)
  ratio.GetYaxis().SetLabelSize(0.16)
  ratio.GetYaxis().SetTitleSize(0.15)
  #ratio.GetYaxis().SetTitleOffset(0.65)
  ratio.GetXaxis().SetTitleOffset(0.8)
  ratio.Draw()

  can_allCuts.SetName("c_"+fileNames[i])   
  can_allCuts.Write()   
  can_allCuts.SaveAs(outputDir+dataset1[i]+"_"+var+'_allCuts.png')
  can_allCuts.SaveAs(outputDir+dataset1[i]+"_"+var+'_allCuts.pdf')


can_allCuts.Close()

outFile.Close()

##----- keep the GUI alive ------------
#if __name__ == '__main__':
#  rep = ''
#  while not rep in ['q','Q']:
#    rep = raw_input('enter "q" to quit: ')
#    if 1 < len(rep):
#      rep = rep[0]
