#!usr/bin/python
import ROOT
from ROOT import TFile, TH1F, THStack, TCanvas, TMath, gROOT, gPad, TLegend, TPad, TString
from setTDRStyle import setTDRStyle
from array import array

import optparse
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
parser.add_option("--inputList_Madgraph",action="store",type="string",default="list.txt",dest="inputList_Madgraph")
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
inputList_Madgraph = options.inputList_Madgraph
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
hist_allCuts_Madgraph = []

LUMI      = lumi
#PATH      = inputDir

#---- read the list -----------------
lines = [line.strip() for line in open(inputList)]
lines_Madgraph = [line_Madgraph.strip() for line_Madgraph in open(inputList_Madgraph)]

#---- split sample name and xsec
fileNames = []
fileNames_Madgraph = []
xsecs = []
xsecs_Madgraph = []
ii = 0
for line in lines:
  parts = line.split()
  fileNames.append(parts[0])
  xsecs.append(parts[1])
  print ("dataset : %s    xsec : %s" % (fileNames[ii], xsecs[ii]))
  ii+=1

ii = 0
for line_Madgraph in lines_Madgraph:
  parts_Madgraph = line_Madgraph.split()
  fileNames_Madgraph.append(parts_Madgraph[0])
  xsecs_Madgraph.append(parts_Madgraph[1])
  print ("dataset Madgraph : %s    xsec : %s" % (fileNames_Madgraph[ii], xsecs_Madgraph[ii]))


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
  print('Not putting efficiency in the weight!!!!')
  #wt = LUMI*float(xsecs[i_f])*eff/Nev
  wt = LUMI*float(xsecs[i_f])/Nev
  print('weight : %f' % wt)
  h_allCuts.Scale(wt)
  h_allCuts.Rebin(rebin)
  h_allCuts.SetDirectory(0)
  h_allCuts.SetFillColor(ROOT.kBlue-9)
  h_allCuts.SetLineColor(ROOT.kBlue-9)
  h_allCuts.SetMarkerColor(ROOT.kBlue-9)
  hist_allCuts.append(h_allCuts)
   
  i_f += 1

#------- open the files Madgraph -------------
i_f = 0
for f in fileNames_Madgraph:
  inf = TFile.Open(f)
  print inf.GetName()
  
  Nev = inf.Get('DijetFilter/EventCount/EventCounter').GetBinContent(1)
  print ('processed events Madgraph : %s' % Nev)
  wt = 1.0
  #if i_f < 3:

  h_allCuts_Madgraph = TH1F("h_allCuts_Madgraph", "", bins, xmin, xmax)
  h_allCuts_Madgraph.Sumw2()
  tree = inf.Get('rootTupleTree/tree')
  tree.Project(h_allCuts_Madgraph.GetName(), var,' deltaETAjj < 1.3')
  Npassed = h_allCuts_Madgraph.GetEntries()
  eff = float(Npassed)/Nev
  print('eff : %f' % eff)
  print('Not putting efficiency in the weight!!!!')
  #wt = LUMI*float(xsecs[i_f])*eff/Nev
  wt = LUMI*float(xsecs_Madgraph[i_f])/Nev
  print(xsecs_Madgraph[i_f]+"*"+str(LUMI)+"/"+str(Nev))
  print('weight Madgraph: %f' % wt)
  h_allCuts_Madgraph.Scale(wt)
  h_allCuts_Madgraph.Rebin(rebin)
  h_allCuts_Madgraph.SetDirectory(0)
  h_allCuts_Madgraph.SetLineColor(ROOT.kRed)
  h_allCuts_Madgraph.SetMarkerColor(ROOT.kRed)
  hist_allCuts_Madgraph.append(h_allCuts_Madgraph)
   
  i_f += 1


#kFactor = NDAT/NQCD
kFactor = 1.3
print ("kFactor = %f" % kFactor)

                                      
for i in range(0,len(fileNames)) :
  hist_allCuts[i].Scale(kFactor)
for i in range(0,len(fileNames_Madgraph)) :
  hist_allCuts_Madgraph[i].Scale(kFactor)


NQCD_allCuts = hist_allCuts[0].Integral()
NQCD_allCuts_Madgraph = hist_allCuts_Madgraph[0].Integral()

for i in range(0,len(fileNames)) :
  NQCD_allCuts += hist_allCuts[i].Integral()
for i in range(0,len(fileNames_Madgraph)) :
  NQCD_allCuts_Madgraph += hist_allCuts_Madgraph[i].Integral()
  

hist_allCutsQCD = hist_allCuts[0].Clone('hist_allCutsQCD')
hist_allCutsQCD_Madgraph = hist_allCuts_Madgraph[0].Clone('hist_allCutsQCD_Madgraph')

for i in range(1,len(fileNames)):
  hist_allCutsQCD.Add(hist_allCuts[i])
for i in range(0,len(fileNames_Madgraph)) :  
  hist_allCutsQCD_Madgraph.Add(hist_allCuts_Madgraph[i])

#hsQCD_allCuts = THStack('QCD_allCuts','QCD_allCuts')

#for i in range(0,len(fileNames)) :
#  hsQCD_allCuts.Add(hist_allCuts[i])


print ("---- After scaling signal to bkg (if not plotting mjj) -----")
print ("bkg integral all cuts = %f" % NQCD_allCuts)
print ("bkg integral all cuts Madgraph = %f" % NQCD_allCuts_Madgraph)


#----- Drawing  and save on file -----------------------

outFile = TFile(outputDir+"histo_signal_bkg_"+var+".root", "recreate")
outFile.cd()
hist_allCutsQCD.Write()
hist_allCutsQCD_Madgraph.Write()

can_allCuts = TCanvas('can_allCuts_'+var,'can_allCuts_'+var,900,600)

leg = TLegend(0.6, 0.7, 0.85, 0.85)
leg.SetLineColor(0)
leg.SetFillColor(0)
leg.SetBorderSize(0)
leg.AddEntry(hist_allCutsQCD, "QCD Pythia", "f")
leg.AddEntry(hist_allCutsQCD_Madgraph, "QCD Madgraph", "l")

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
#hist_allCutsQCD_Madgraph.Rebin()
hist_allCutsQCD.Draw("hist")
hist_allCutsQCD_Madgraph.Draw("hist same")
leg.Draw()
#gPad.RedrawAxis()

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
ratio.Divide(hist_allCutsQCD_Madgraph)
ratio.SetFillColor(0)
ratio.SetLineColor(ROOT.kBlack)
ratio.SetMarkerColor(ROOT.kBlack)
ratio.Draw("hist")
ratio.GetYaxis().SetRangeUser(0., 2.)
ratio.GetYaxis().SetNdivisions(4, ROOT.kTRUE)
ratio.GetYaxis().SetTitleFont(42)
ratio.GetYaxis().SetTitle("#frac{Pythia}{Madgraph}")
ratio.GetXaxis().SetTitleSize(0.2)
ratio.GetXaxis().SetLabelSize(0.16)
ratio.GetYaxis().SetLabelSize(0.16)
ratio.GetYaxis().SetTitleSize(0.15)
#ratio.GetYaxis().SetTitleOffset(0.65)
ratio.GetXaxis().SetTitleOffset(0.8)
ratio.Draw("hist")

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
