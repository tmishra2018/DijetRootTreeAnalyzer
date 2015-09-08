#!usr/bin/python

import sys, os, subprocess, string, re
from ROOT import *
from array import array
import CMS_lumi
import optparse
from setTDRStyle import setTDRStyle

CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "(13 TeV)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 11
iPeriod = 0
#######################
minX_mass = 1000.
maxX_mass = 7060. 

massBins_list = [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000]
    
massBins = array("d",massBins_list)
#############################

gROOT.SetBatch(kTRUE);
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetTitleFont(42, "XYZ")
gStyle.SetTitleSize(0.06, "XYZ")
gStyle.SetLabelFont(42, "XYZ")
gStyle.SetLabelSize(0.05, "XYZ")
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

usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
#parser.add_option("--var",action="store",type="string",dest="var",default='ptHat')
#parser.add_option("--xmin",action="store",type="float",dest="xmin",default=1)
#parser.add_option("--xmax",action="store",type="float",dest="xmax",default=1)
#parser.add_option("--xtitle",action="store",type="string",dest="xtitle",default='')
#parser.add_option("--bins",action="store",type="int",dest="bins",default=11111111111)
#parser.add_option("--rebin",action="store",type="int",dest="rebin",default=1)
#parser.add_option("--logy",action="store_true",default=False,dest="logy")
parser.add_option("--outputDir",action="store",type="string",default="./",dest="outputDir")
parser.add_option("--inputList_minus",action="store",type="string",default="list_minus.txt",dest="inputList_minus")
parser.add_option("--inputList",action="store",type="string",default="list.txt",dest="inputList")
parser.add_option("--inputList_plus",action="store",type="string",default="list_plus.txt",dest="inputList_plus")
parser.add_option("--lumi",action="store",type="float",default="41.8",dest="lumi")

(options, args) = parser.parse_args()

#var = options.var
#xmin = options.xmin
#xmax = options.xmax
#bins = options.bins
#xtitle = options.xtitle
#rebin = options.rebin
#logy = options.logy
outputDir = options.outputDir
inputList = options.inputList
inputList_minus = options.inputList_minus
inputList_plus = options.inputList_plus
lumi = options.lumi

#fileNames = ['QCD_Pt-301to470','QCD_Pt-470to600','QCD_Pt-600to800', 'QCD_Pt-800to1000', 'QCD_Pt-1000to1400', 'QCD_Pt-1400to1800', 'QCD_Pt-1800to2400', 'QCD_Pt-2400to3200', 'QCD_Pt-3200']
#xsections = [7475 ,587., 167, 28.25, 8.195, 0.7346, 0.102, 0.00644, 0.000163]
#colorF    = [ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9,ROOT.kBlue-9]
#colorL    = [ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack,ROOT.kBlack]
hist_allCuts      = []
hist_allCuts_plus      = []
hist_allCuts_minus      = []

LUMI      = lumi
#PATH      = inputDir

#---- read the list -----------------
lines = [line.strip() for line in open(inputList)]
lines_plus = [line.strip() for line in open(inputList_plus)]
lines_minus = [line.strip() for line in open(inputList_minus)]

#---- split sample name and xsec
fileNames = []
fileNames_plus = []
fileNames_minus = []
xsecs = []
hist_JEC_minus = []
hist_JEC_plus = []

os.system("mkdir -p "+options.outputDir)

outFile = TFile(outputDir+"rootfile_JESunc.root", "recreate")

ii = 0
for line in lines:
  parts = line.split()
  fileNames.append(parts[0])
  xsecs.append(parts[1])
  print ("dataset : %s    xsec : %s" % (fileNames[ii], xsecs[ii]))
  ii+=1
ii = 0
for line in lines_plus:
  parts = line.split()
  fileNames_plus.append(parts[0])
  print ("dataset : %s    xsec : %s" % (fileNames_plus[ii], xsecs[ii]))
  ii+=1
ii = 0
for line in lines_minus:
  parts = line.split()
  fileNames_minus.append(parts[0])
  print ("dataset : %s    xsec : %s" % (fileNames_minus[ii], xsecs[ii]))
  ii+=1

#---- open the files --------------------
#i_f = 0
for i_f in range(0,9):
  inf = TFile.Open(fileNames[i_f])
  inf_minus = TFile.Open(fileNames_minus[i_f])
  inf_plus = TFile.Open(fileNames_plus[i_f])
  print inf.GetName()
  
  Nev = inf.Get('DijetFilter/EventCount/EventCounter').GetBinContent(1)
  print ('processed events: %s' % Nev)
  wt = 1.0
  #if i_f < 3:

  h_mass  = TH1F("h_mass_"+str(i_f),"",len(massBins)-1,massBins)
  h_mass_plus  = TH1F("h_mass_plus_"+str(i_f),"",len(massBins)-1,massBins)
  h_mass_minus  = TH1F("h_mass_minus_"+str(i_f),"",len(massBins)-1,massBins)
  h_JES_plus = TH1F("h_JES_plus_"+str(i_f), "",len(massBins)-1,massBins)
  h_JES_minus = TH1F("h_JES_minus_"+str(i_f), "",len(massBins)-1,massBins)
  tree = inf.Get('rootTupleTree/tree')
  tree_plus = inf_plus.Get('rootTupleTree/tree')
  tree_minus = inf_minus.Get('rootTupleTree/tree')
 
  #tree.Project("h_mass"+str(i_f),"mjj","")

  for ev in range(0,100000):
  #for ev in range(0,tree.GetEntries()):
    if (ev % 10000 == 0): 
      print "processing %d" % ev
    
    tree.GetEntry(ev)
    run = tree.run 
    lumi = tree.lumi 
    event = tree.event 
    mjj = tree.mjj
    pTWJ_j1 = tree.pTWJ_j1
    pTWJ_j2 = tree.pTWJ_j2
    etaWJ_j1 = tree.etaWJ_j1
    etaWJ_j2 = tree.etaWJ_j2
    if (mjj>1118):  
      h_mass.Fill(mjj)
    
    tree_minus.GetEntry(ev)
    run_minus = tree_minus.run 
    lumi_minus = tree_minus.lumi 
    event_minus = tree_minus.event 
    pTWJ_j1_minus = tree_minus.pTWJ_j1
    pTWJ_j2_minus = tree_minus.pTWJ_j2
    etaWJ_j1_minus = tree_minus.etaWJ_j1
    etaWJ_j2_minus = tree_minus.etaWJ_j2
    
    tree_plus.GetEntry(ev)
    run_plus = tree_plus.run 
    lumi_plus = tree_plus.lumi
    event_plus = tree_plus.event 
    #print "%f   %f    %f    %f" % (run,lumi,event,mjj)
    #print "%f   %f    %f    %f" % (run_minus,lumi_minus,event_minus,mjj_minus)
    err = 0 

    tree.GetEntry(ev)
    mjj = tree.mjj
    tree_minus.GetEntry(ev)
    mjj_minus = tree_minus.mjj
    if (run_minus == run and lumi_minus==lumi and event_minus==event and mjj_minus>1118):
      h_mass_minus.Fill(mjj_minus)
      err = (mjj_minus/mjj -1) * 100
      if err>0:
	print "%f   %f    %f  %f  %f   %f  %f   %f     %f" % (run,lumi,event,pTWJ_j1,pTWJ_j2,etaWJ_j1,etaWJ_j2,mjj,err)
	print "%f   %f    %f  %f  %f   %f  %f   %f     %f" % (run_minus,lumi_minus,event_minus,pTWJ_j1_minus,pTWJ_j2_minus,etaWJ_j1_minus,etaWJ_j2_minus,mjj_minus,err)
        print ""
      bin = h_JES_minus.FindBin(mjj)
      #h_JES_minus.SetBinContent(bin,h_JES_minus.GetBinContent(bin)+err)
      h_JES_minus.Fill(mjj_minus,err)
      
    tree.GetEntry(ev)
    mjj = tree.mjj
    tree_plus.GetEntry(ev)
    mjj_plus = tree_plus.mjj
    if (run_plus == run and lumi_plus==lumi and event_plus==event and mjj_plus>1118):
      err = 0 
      h_mass_plus.Fill(mjj_plus)
      err = (mjj_plus/mjj -1) * 100
      #print "%f   %f    %f   %f   %f   %f" % (run,lumi,event,mjj,mjj_plus,err)
      bin = h_JES_plus.FindBin(mjj)
      #h_JES_plus.SetBinContent(bin,h_JES_plus.GetBinContent(bin)+err)
      h_JES_plus.Fill(mjj_plus,err)

  #for bin in range(1,len(massBins)):
  #  content_minus = h_JES_minus.GetBinContent(bin)
  #  content_plus = h_JES_plus.GetBinContent(bin)
  #  #print content_minus
  #  #print content_plus
  #  entries = h_mass.GetBinContent(bin) 
  #  entries_plus = h_mass_plus.GetBinContent(bin) 
  #  entries_minus = h_mass_minus.GetBinContent(bin) 
  #  #print"entries %d" % entries 
  #  if entries_plus>0:
  #    h_JES_plus.SetBinContent(bin, content_plus/entries_plus )
  #  if entries_minus>0:
  #    h_JES_minus.SetBinContent(bin, content_minus/entries_minus )
   
 
  h_JES_minus.Print()
  h_JES_plus.Print()
  h_JES_minus.SetLineColor(kRed)
  
  hist_JEC_minus.append(h_JES_minus)
  hist_JEC_plus.append(h_JES_plus)
  outFile.cd()
  h_mass.Write()
  h_mass_plus.Write()
  h_mass_minus.Write()
  h_JES_minus.Write()
  h_JES_plus.Write()
  #i_f+=1

print hist_JEC_minus
print hist_JEC_plus

h_mass_all = outFile.Get("h_mass_0").Clone("h_mass_all")
h_mass_plus_all = outFile.Get("h_mass_plus_0").Clone("h_mass_plus_all")
h_mass_minus_all = outFile.Get("h_mass_minus_0").Clone("h_mass_minus_all")
h_JES_minus_all = outFile.Get("h_JES_minus_0").Clone("h_JES_minus_all")
h_JES_plus_all = outFile.Get("h_JES_plus_0").Clone("h_JES_plus_all")
for i in range(1,9):
  h_mass_all.Add( outFile.Get("h_mass_"+str(i)))
  h_mass_plus_all.Add( outFile.Get("h_mass_plus_"+str(i)))
  h_mass_minus_all.Add( outFile.Get("h_mass_minus_"+str(i)))
  h_JES_minus_all.Add( outFile.Get("h_JES_minus_"+str(i)))
  h_JES_plus_all.Add( outFile.Get("h_JES_plus_"+str(i)))

for bin in range(1,len(massBins)):
  content_minus = h_JES_minus_all.GetBinContent(bin)
  content_plus = h_JES_plus_all.GetBinContent(bin)
  #print content_minus
  #print content_plus
  entries = h_mass_all.GetBinContent(bin) 
  entries_plus = h_mass_plus_all.GetBinContent(bin) 
  entries_minus = h_mass_minus_all.GetBinContent(bin) 
  #print"entries %d" % entries 
  if entries_plus>0:
    h_JES_plus_all.SetBinContent(bin, content_plus/entries_plus )
  if entries_minus>0:
    h_JES_minus_all.SetBinContent(bin, content_minus/entries_minus )
  print h_JES_minus.GetBinContent(bin) 
  print h_JES_plus.GetBinContent(bin) 

outFile.cd()
h_mass_all.Write()
h_JES_minus_all.Write()
h_JES_plus_all.Write()

can = TCanvas('can','',600,600)

leg = TLegend(0.6, 0.7, 0.85, 0.85)
leg.SetLineColor(0)
leg.SetFillColor(0)
leg.AddEntry(h_JES_minus_all, "JES unc. minus", "l")
leg.AddEntry(h_JES_plus_all, "JES unc.plus", "l")

can.cd()
can.SetTitle("JES Uncertainty at 13 TeV")
h_JES_minus_all.SetLineColor(kRed)
h_JES_minus_all.SetFillColor(kYellow)
h_JES_plus_all.SetFillColor(kYellow)
h_JES_minus_all.GetXaxis().SetRangeUser(minX_mass,maxX_mass)
h_JES_minus_all.GetYaxis().SetRangeUser(-5,5)
h_JES_minus_all.GetXaxis().SetTitle("Dijet Mass (GeV)")
h_JES_minus_all.GetYaxis().SetTitle("JES Uncertainty (%)")
h_JES_minus_all.Draw("hist")
h_JES_plus_all.Draw("hist same")
leg.Draw()

#draw the lumi text on the canvas
CMS_lumi.CMS_lumi(can, iPeriod, iPos)
gPad.RedrawAxis()

can.Write()   
can.SaveAs(outputDir+'JESuncertainty.png')
can.SaveAs(outputDir+'JESuncertainty.pdf')
can.Close()

outFile.Close()

##----- keep the GUI alive ------------
#if __name__ == '__main__':
#  rep = ''
#  while not rep in ['q','Q']:
#    rep = raw_input('enter "q" to quit: ')
#    if 1 < len(rep):
#      rep = rep[0]
