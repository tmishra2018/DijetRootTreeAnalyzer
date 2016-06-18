#! /usr/bin/env python

import ROOT as rt
import os.path
import sys, glob, re
from array import *
from optparse import OptionParser
     
def file_key(filename):
    massPoint = re.findall("[0-9]+.000000",filename)
    gluinoMass    = massPoint[0]
    LSPMass  = massPoint[1]
    return float(gluinoMass)
    
def getHybridCLsArrays(directory, model, Box):
    tfile = rt.TFile.Open("%s/xsecUL_Asymptotic_%s_%s.root"%(directory,model,Box))
    xsecTree = tfile.Get("xsecTree")
    
    gluinoMassArray = array('d')
    gluinoMassArray_er = array('d')
    observedLimit = array('d')
    observedLimit_er = array('d')
    expectedLimit = array('d')
    expectedLimit_minus1sigma = array('d')
    expectedLimit_plus1sigma = array('d')
    expectedLimit_minus2sigma = array('d')
    expectedLimit_plus2sigma = array('d')

    
    xsecTree.Draw('>>elist','','entrylist')
    elist = rt.gDirectory.Get('elist')
    entry = -1
    while True:
        entry = elist.Next()
        if entry == -1: break
        xsecTree.GetEntry(entry)

        gluinoMassArray.append(xsecTree.mass)
        gluinoMassArray_er.append(0.0)
        
        exec 'xsecULObs = xsecTree.xsecULObs_%s'%Box
        exec 'xsecULExp = xsecTree.xsecULExp_%s'%Box
        exec 'xsecULExpPlus = xsecTree.xsecULExpPlus_%s'%Box
        exec 'xsecULExpMinus = xsecTree.xsecULExpMinus_%s'%Box
        exec 'xsecULExpPlus2 = xsecTree.xsecULExpPlus2_%s'%Box
        exec 'xsecULExpMinus2 = xsecTree.xsecULExpMinus2_%s'%Box

            
            
        xsecULObs = xsecULObs
        xsecULExp = xsecULExp
        observedLimit.append(xsecULObs)#*crossSections[i])
        observedLimit_er.append(0.0)#*crossSections[i])

        expectedLimit.append(xsecULExp)#*crossSections[i])
        
            

        xsecULExpPlus = max(xsecULExpPlus,xsecULExp)
        xsecULExpMinus = min(xsecULExpMinus,xsecULExp)
        xsecULExpPlus2 = max(xsecULExpPlus2,xsecULExpPlus)
        xsecULExpMinus2 = min(xsecULExpMinus2,xsecULExpMinus)

        expectedLimit_minus1sigma.append(xsecULExp - xsecULExpMinus)#*crossSections[i])
        expectedLimit_plus1sigma.append(xsecULExpPlus - xsecULExp)#*crossSections[i])
        expectedLimit_minus2sigma.append(xsecULExp - xsecULExpMinus2)#*crossSections[i])
        expectedLimit_plus2sigma.append(xsecULExpPlus2 - xsecULExp)#*crossSections[i])
    

    return gluinoMassArray, gluinoMassArray_er, observedLimit, observedLimit_er, expectedLimit, expectedLimit_minus1sigma, expectedLimit_plus1sigma, expectedLimit_minus2sigma, expectedLimit_plus2sigma


def getSignificanceArrays(directory, model, Box):
    tfile = rt.TFile.Open("%s/xsecUL_ProfileLikelihood_%s_%s.root"%(directory,model,Box))
    xsecTree = tfile.Get("xsecTree")
    
    gluinoMassArray = array('d')
    gluinoMassArray_er = array('d')
    observedLimit = array('d')
    observedLimit_er = array('d')
    expectedLimit = array('d')
    expectedLimit_minus1sigma = array('d')
    expectedLimit_plus1sigma = array('d')
    expectedLimit_minus2sigma = array('d')
    expectedLimit_plus2sigma = array('d')

    
    xsecTree.Draw('>>elist','','entrylist')
    elist = rt.gDirectory.Get('elist')
    entry = -1
    while True:
        entry = elist.Next()
        if entry == -1: break
        xsecTree.GetEntry(entry)

        gluinoMassArray.append(xsecTree.mass)
        gluinoMassArray_er.append(0.0)
        
        exec 'xsecULObs = xsecTree.xsecULObs_%s'%Box
        exec 'xsecULExp = xsecTree.xsecULExp_%s'%Box
        exec 'xsecULExpPlus = xsecTree.xsecULExpPlus_%s'%Box
        exec 'xsecULExpMinus = xsecTree.xsecULExpMinus_%s'%Box
        exec 'xsecULExpPlus2 = xsecTree.xsecULExpPlus2_%s'%Box
        exec 'xsecULExpMinus2 = xsecTree.xsecULExpMinus2_%s'%Box

            
            
        xsecULObs = xsecULObs
        xsecULExp = xsecULExp
        observedLimit.append(xsecULObs)#*crossSections[i])
        observedLimit_er.append(0.0)#*crossSections[i])

        expectedLimit.append(xsecULExp)#*crossSections[i])
        
            

        xsecULExpPlus = max(xsecULExpPlus,xsecULExp)
        xsecULExpMinus = min(xsecULExpMinus,xsecULExp)
        xsecULExpPlus2 = max(xsecULExpPlus2,xsecULExpPlus)
        xsecULExpMinus2 = min(xsecULExpMinus2,xsecULExpMinus)

        expectedLimit_minus1sigma.append(xsecULExp - xsecULExpMinus)#*crossSections[i])
        expectedLimit_plus1sigma.append(xsecULExpPlus - xsecULExp)#*crossSections[i])
        expectedLimit_minus2sigma.append(xsecULExp - xsecULExpMinus2)#*crossSections[i])
        expectedLimit_plus2sigma.append(xsecULExpPlus2 - xsecULExp)#*crossSections[i])
    

    return gluinoMassArray, gluinoMassArray_er, observedLimit, observedLimit_er, expectedLimit, expectedLimit_minus1sigma, expectedLimit_plus1sigma, expectedLimit_minus2sigma, expectedLimit_plus2sigma
    
def setstyle():
    # For the canvas:
    rt.gStyle.SetCanvasBorderMode(0)
    rt.gStyle.SetCanvasColor(rt.kWhite)
    rt.gStyle.SetCanvasDefH(400) #Height of canvas
    rt.gStyle.SetCanvasDefW(600) #Width of canvas
    rt.gStyle.SetCanvasDefX(0)   #POsition on screen
    rt.gStyle.SetCanvasDefY(0)
    
    # For the Pad:
    rt.gStyle.SetPadBorderMode(0)
    # rt.gStyle.SetPadBorderSize(Width_t size = 1)
    rt.gStyle.SetPadColor(rt.kWhite)
    rt.gStyle.SetPadGridX(False)
    rt.gStyle.SetPadGridY(False)
    rt.gStyle.SetGridColor(0)
    rt.gStyle.SetGridStyle(3)
    rt.gStyle.SetGridWidth(1)
    
    # For the frame:
    rt.gStyle.SetFrameBorderMode(0)
    rt.gStyle.SetFrameBorderSize(1)
    rt.gStyle.SetFrameFillColor(0)
    rt.gStyle.SetFrameFillStyle(0)
    rt.gStyle.SetFrameLineColor(1)
    rt.gStyle.SetFrameLineStyle(1)
    rt.gStyle.SetFrameLineWidth(1)
    
    # set the paper & margin sizes
    rt.gStyle.SetPaperSize(20,26)
    rt.gStyle.SetPadTopMargin(0.09)
    rt.gStyle.SetPadRightMargin(0.065)
    rt.gStyle.SetPadBottomMargin(0.15)
    rt.gStyle.SetPadLeftMargin(0.17)
    
    # use large Times-Roman fonts
    rt.gStyle.SetTitleFont(42,"xyz")  # set the all 3 axes title font
    rt.gStyle.SetTitleFont(42," ")    # set the pad title font
    rt.gStyle.SetTitleSize(0.06,"xyz") # set the 3 axes title size
    rt.gStyle.SetTitleSize(0.06," ")   # set the pad title size
    rt.gStyle.SetTitleSize(0.052,"y")   # set the pad title size
    rt.gStyle.SetTitleOffset(1.2,"y")   # set the pad title size
    rt.gStyle.SetLabelFont(42,"xyz")
    rt.gStyle.SetLabelSize(0.05,"xyz")
    rt.gStyle.SetLabelColor(1,"xyz")
    rt.gStyle.SetTextFont(42)
    rt.gStyle.SetTextSize(0.08)
    rt.gStyle.SetStatFont(42)
    
    # use bold lines and markers
    rt.gStyle.SetMarkerStyle(8)
    rt.gStyle.SetLineStyleString(2,"[12 12]") # postscript dashes
    
    #..Get rid of X error bars
    rt.gStyle.SetErrorX(0.001)
    
    # do not display any of the standard histogram decorations
    rt.gStyle.SetOptTitle(0)
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetOptFit(11111111)
    
    # put tick marks on top and RHS of plots
    rt.gStyle.SetPadTickX(1)
    rt.gStyle.SetPadTickY(1)
    
    ncontours = 999
    
    stops = [ 0.00, 0.34, 0.61, 0.84, 1.00 ]
    red =   [ 1.0,   0.95,  0.95,  0.65,   0.15 ]
    green = [ 1.0,  0.85, 0.7, 0.5,  0.3 ]
    blue =  [ 0.95, 0.6 , 0.3,  0.45, 0.65 ]
    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)
        
    npoints = len(s)
    rt.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    rt.gStyle.SetNumberContours(ncontours)
   
    rt.gStyle.cd()
        

if __name__ == '__main__':
    
    rt.gROOT.SetBatch()
    parser = OptionParser()
    parser.add_option('-b','--box',dest="box", default="CaloDijet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="gg",type="string",
                  help="signal model name")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Input/Output directory to store output")    
    parser.add_option('-l','--lumi',dest="lumi", default=1.,type="float",
                  help="integrated luminosity in fb^-1")
    parser.add_option('--signif',dest="doSignificance",default=False,action='store_true',
                  help="for significance instead of limit")
    
    (options,args) = parser.parse_args()
    Box = options.box
    model = options.model
    directory      = options.outDir

    box = Box.lower()
    
    thyXsec = {500: 10, 1000: 10, 1450: 10}
    thyXsecErr = {500: 0.1, 1000: 0.1, 1450: 0.1}

    mass_xsec = array('d')
    sig_xsec = array('d')
    sig_xsec_up = array('d')
    sig_xsec_down = array('d')
    sig_xsec_err = array('d')

    for mg in thyXsec:
        mass_xsec.append(mg)
        sig_xsec.append(thyXsec[mg])
        sig_xsec_up.append(thyXsec[mg]*(1.+thyXsecErr[mg]))
        sig_xsec_down.append(thyXsec[mg]*(1.-thyXsecErr[mg]))
        sig_xsec_err.append(thyXsec[mg]*thyXsecErr[mg])
        
    X_xsec = array('d')
    Y_xsec = array('d')
    X_xsec_up = array('d')
    Y_xsec_up = array('d')
    X_xsec_down = array('d')
    Y_xsec_down = array('d')
    X_xsec_err  = array('d')
    Y_xsec_err = array('d')

    N_g_xsec = len(mass_xsec)
    for i in range(0,N_g_xsec):
        X_xsec.append(mass_xsec[i])
        Y_xsec.append(sig_xsec[i])
        X_xsec_up.append(mass_xsec[i])
        Y_xsec_up.append(sig_xsec_up[i])
        X_xsec_down.append(mass_xsec[i])
        Y_xsec_down.append(sig_xsec_down[i])
        X_xsec_err.append(0)
        Y_xsec_err.append(sig_xsec_err[i])

    setstyle()
    c = rt.TCanvas("c","c",500,400)
    #c = rt.TCanvas("c","c",800,800)    
    #c.SetBottomMargin(0.14)
    #c.SetTopMargin(0.06)
    #c.SetLeftMargin(0.15)
    #c.SetRightMargin(0.05)
    if options.doSignificance:
        c.SetLogy(0)
    else:        
        c.SetLogy()
    
    xsec_gr_nom = rt.TGraph(N_g_xsec, X_xsec, Y_xsec)
    xsec_gr_nom.SetMarkerSize(0)
    xsec_gr_nom.SetLineWidth(2)
    xsec_gr_nom.SetLineStyle(1)
    xsec_gr_nom.SetLineColor(rt.kOrange)
    
    xsec_gr = rt.TGraphErrors(N_g_xsec, X_xsec, Y_xsec, X_xsec_err, Y_xsec_err)  
    xsec_gr.SetMarkerStyle(23)
    xsec_gr_nom.SetLineWidth(2)
    xsec_gr.SetMarkerSize(0)
    xsec_gr.SetLineColor(rt.kOrange)
    xsec_gr.SetFillColor(rt.kBlue-7)

    h_limit = rt.TMultiGraph()

    if options.doSignificance:
        h_limit.SetTitle(" ;Resonance Mass m_{X} [GeV];Local Significance n#sigma")
    else:
        h_limit.SetTitle(" ;Resonance Mass m_{X} [GeV];95% C.L. upper limit on cross section [pb]")
    
    if options.doSignificance: 
        gluinoMassArray, gluinoMassArray_er, observedLimit, observedLimit_er, expectedLimit, expectedLimit_minus1sigma, expectedLimit_plus1sigma, expectedLimit_minus2sigma, expectedLimit_plus2sigma = getSignificanceArrays(directory, model, Box)
    else:        
        gluinoMassArray, gluinoMassArray_er, observedLimit, observedLimit_er, expectedLimit, expectedLimit_minus1sigma, expectedLimit_plus1sigma, expectedLimit_minus2sigma, expectedLimit_plus2sigma = getHybridCLsArrays(directory, model, Box)
    
    rt.gStyle.SetOptStat(0)
    
    nPoints = len(observedLimit)
    
    gr_observedLimit = rt.TGraph(nPoints, gluinoMassArray, observedLimit)
    gr_observedLimit.SetMarkerColor(1)
    gr_observedLimit.SetMarkerStyle(22)
    gr_observedLimit.SetMarkerSize(1)
    gr_observedLimit.SetLineWidth(3)
    gr_observedLimit.SetLineColor(rt.kBlack)


    gr_expectedLimit = rt.TGraph(nPoints, gluinoMassArray, expectedLimit)
    gr_expectedLimit.SetLineWidth(3)
    gr_expectedLimit.SetLineStyle(2)
   
    gr_expectedLimit2sigma = rt.TGraphAsymmErrors(nPoints, gluinoMassArray, expectedLimit, gluinoMassArray_er, gluinoMassArray_er, expectedLimit_minus2sigma, expectedLimit_plus2sigma)
    gr_expectedLimit2sigma.SetLineColor(5)
    gr_expectedLimit2sigma.SetFillColor(5)
    gr_expectedLimit2sigma.SetFillStyle(1001)
   
    gr_expectedLimit1sigma = rt.TGraphAsymmErrors(nPoints, gluinoMassArray, expectedLimit, gluinoMassArray_er, gluinoMassArray_er, expectedLimit_minus1sigma, expectedLimit_plus1sigma)
    
    #col1 = rt.gROOT.GetColor(rt.kGreen-7)
    #col1.SetAlpha(0.5)
    gr_expectedLimit1sigma.SetLineColor(rt.kGreen-7)
    gr_expectedLimit1sigma.SetFillColor(rt.kGreen-7)
    #gr_expectedLimit1sigma.SetFillStyle(3001)

    h_limit.Add(gr_expectedLimit2sigma)
    h_limit.Add(gr_expectedLimit1sigma)
    h_limit.Add(gr_observedLimit)
    #h_limit.Add(xsec_gr)
    #h_limit.Add(xsec_gr_nom)

        
    h_limit.Draw("a3")
    #h_limit.GetXaxis().SetLimits(620,1580)
    h_limit.GetXaxis().SetLimits(500,1600)
    if options.doSignificance:
        h_limit.SetMaximum(3)
        h_limit.SetMinimum(0)
    else:
        h_limit.SetMaximum(1000)
        h_limit.SetMinimum(1e-1)

        
    h_limit.Draw("a3")
    
    if options.doSignificance:
        gr_observedLimit.SetMarkerStyle(21)
        gr_observedLimit.SetMarkerSize(0.6)
        gr_observedLimit.SetLineColor(rt.kRed)
        gr_observedLimit.SetMarkerColor(rt.kBlue)
        gr_observedLimit.Draw("lp SAME")
    else:
        gr_expectedLimit.Draw("c same")
        #xsec_gr_nom.Draw("c same")
        gr_observedLimit.Draw("c SAME")

    
    l = rt.TLatex()
    l.SetTextAlign(11)
    l.SetTextSize(0.05)
    l.SetNDC()
    l.SetTextFont(62)
    l.DrawLatex(0.17,0.92,"CMS")
    #l.DrawLatex(0.16,0.95,"CMS")
    l.SetTextFont(52)
    l.DrawLatex(0.26,0.92,"Preliminary")
    #l.DrawLatex(0.27,0.95,"Preliminary")
    l.SetTextFont(42)
    l.DrawLatex(0.65,0.92,"%.0f pb^{-1} (13 TeV)"%(options.lumi*1000))
    #l.DrawLatex(0.59,0.95,"%.0f pb^{-1} (13 TeV)"%(options.lumi*1000))

    if model=="gg":
        l.DrawLatex(0.3,0.8,"gg #rightarrow X #rightarrow jj")
    elif model=="gq":        
        l.DrawLatex(0.3,0.8,"gq #rightarrow X #rightarrow jj")
    elif model=="qq":
        l.DrawLatex(0.3,0.8,"qq #rightarrow X #rightarrow jj")

    if options.doSignificance:
        leg = rt.TLegend(0.55,0.803,0.92,0.87)
    else:        
        leg = rt.TLegend(0.55,0.67,0.92,0.87)
    #leg.AddEntry(xsec_gr, "#sigma_{NLO+NLL} (#tilde{g}#tilde{g}) #pm 1 #sigma (theory)","lf")

    demo = gr_expectedLimit1sigma.Clone()
    demo.SetLineColor(rt.kBlack)
    demo.SetLineStyle(2)
    demo.SetLineWidth(3)
    
    leg.SetTextFont(42)
    leg.SetFillColor(rt.kWhite)
    leg.SetLineColor(rt.kWhite)
    #leg.AddEntry(gr_observedLimit, "observed 0-lep+1-lep","l")
    #leg.AddEntry(gr_expectedLimit, "expected 0-lep+1-lep","l")
    
    if options.doSignificance:
        leg.AddEntry(gr_observedLimit, "observed","lp")
    else:
        leg.AddEntry(gr_observedLimit, "observed","l")
    if not options.doSignificance:
        leg.AddEntry(demo, "expected #pm1#sigma","lf")
    gr_expectedLimit2sigma.SetLineStyle(2)
    gr_expectedLimit2sigma.SetLineWidth(3)
    gr_expectedLimit2sigma.SetLineColor(rt.kBlack)
    
    if not options.doSignificance:
        leg.AddEntry(gr_expectedLimit2sigma, "expected #pm2#sigma","lf")

    leg.Draw("SAME")

    if options.doSignificance:
        c.SaveAs(directory+"/signif_"+model+"_"+box+".pdf")
        c.SaveAs(directory+"/signif_"+model+"_"+box+".C")
    else:
        c.SaveAs(directory+"/limits_"+model+"_"+box+".pdf")
        c.SaveAs(directory+"/limits_"+model+"_"+box+".C")

