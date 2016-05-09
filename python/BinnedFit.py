from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
from array import *
from itertools import *
from operator import *
from WriteDataCard import initializeWorkspace,convertToTh1xHist
import os
import random
import sys
import math

densityCorr = False

def binnedFit(pdf, data, fitRange='Full',useWeight=False):

    if useWeight:
        fr = pdf.fitTo(data,rt.RooFit.Range(fitRange),rt.RooFit.Extended(True),rt.RooFit.SumW2Error(True),rt.RooFit.Save(),rt.RooFit.Minimizer('Minuit2','migrad'),rt.RooFit.Strategy(2))
        migrad_status = fr.status()
        hesse_status = -1        
    else:
        nll = pdf.createNLL(data,rt.RooFit.Range(fitRange),rt.RooFit.Extended(True),rt.RooFit.Offset(False))
        m2 = rt.RooMinimizer(nll)
        m2.setStrategy(2)
        #m2.setEps(1e-5)
        m2.setMaxFunctionCalls(100000)
        m2.setMaxIterations(100000)
        migrad_status = m2.minimize('Minuit2','migrad')
        improve_status = m2.minimize('Minuit2','improve')
        hesse_status = m2.minimize('Minuit2','hesse')
        minos_status = m2.minos()
        fr = m2.save()

    if fr.covQual() != 3:
        print ""
        print "CAUTION: COVARIANCE QUALITY < 3"
        print ""
        
    if migrad_status != 0:
        print ""
        print "CAUTION: MIGRAD STATUS ! = 0"
        print ""

    if hesse_status != 0:
        print ""
        print "CAUTION: HESSE STATUS ! = 0"
        print ""
        
    return fr


def convertSideband(name,w,x):
    if name=="Full":
        return "Full"
    names = name.split(',')
    nBins = (len(x)-1)
    iBinX = -1
    sidebandBins = []
    for ix in range(1,len(x)):
        iBinX+=1
        w.var('mjj').setVal((x[ix]+x[ix-1])/2.)
        inSideband = 0
        for fitname in names:
            inSideband += ( w.var('mjj').inRange(fitname) )
        if inSideband: sidebandBins.append(iBinX)

    sidebandGroups = []
    for k, g in groupby(enumerate(sidebandBins), lambda (i,x):i-x):
        consecutiveBins = map(itemgetter(1), g)
        sidebandGroups.append([consecutiveBins[0],consecutiveBins[-1]+1])
        
    newsidebands = ''
    nameNoComma = name.replace(',','')
        
    for iSideband, sidebandGroup in enumerate(sidebandGroups):
        if not w.var('th1x').hasRange('%s%i'%(nameNoComma,iSideband)):
            w.var('th1x').setRange("%s%i"%(nameNoComma,iSideband),sidebandGroup[0],sidebandGroup[1])
        newsidebands+='%s%i,'%(nameNoComma,iSideband)
    newsidebands = newsidebands[:-1]
    return newsidebands

def convertFunctionToHisto(background_,name_,N_massBins_,massBins_):

    background_hist_ = rt.TH1D(name_,name_,N_massBins_,massBins_)

    for bin in range (0,N_massBins_):
        xbinLow = massBins_[bin]
        xbinHigh = massBins_[bin+1]
        binWidth_current = xbinHigh - xbinLow
        value = background_.Integral(xbinLow , xbinHigh) / binWidth_current
        background_hist_.SetBinContent(bin+1,value)

    return background_hist_

def calculateChi2AndFillResiduals(data_obs_TGraph_,background_hist_,hist_fit_residual_vsMass_,workspace_,prinToScreen_=0):
    
    N_massBins_ = data_obs_TGraph_.GetN()
    MinNumEvents = 10
    nParFit = 4
    if workspace_.var('meff')>0 and workspace_.var('seff')>0 :
        nParFit = 6

    chi2_FullRangeAll = 0
    chi2_PlotRangeAll = 0
    chi2_PlotRangeNonZero = 0
    chi2_PlotRangeMinNumEvents = 0 

    N_FullRangeAll = 0
    N_PlotRangeAll = 0
    N_PlotRangeNonZero = 0
    N_PlotRangeMinNumEvents = 0 

    for bin in range (0,N_massBins_):
        ## Values and errors
        value_data = data_obs_TGraph_.GetY()[bin]
        err_low_data = data_obs_TGraph_.GetEYlow()[bin]
        err_high_data = data_obs_TGraph_.GetEYhigh()[bin]
        xbinCenter = data_obs_TGraph_.GetX()[bin] 
        xbinLow = data_obs_TGraph_.GetX()[bin]-data_obs_TGraph_.GetEXlow()[bin] 
        xbinHigh = data_obs_TGraph_.GetX()[bin]+data_obs_TGraph_.GetEXhigh()[bin]
        binWidth_current = xbinHigh - xbinLow
        #value_fit = background_.Integral(xbinLow , xbinHigh) / binWidth_current
        value_fit = background_hist_.GetBinContent(bin+1)
        
        ## Fit residuals

        err_tot_data = 0
        if (value_fit >= value_data):
            err_tot_data = err_high_data  
        else:
            err_tot_data = err_low_data  
        if plotRegion=='Full' or (plotRegion=='Low,High' and (xbinCenter<workspace_.var('mjj').getMin('Blind') or xbinCenter>workspace_.var('mjj').getMax('Blind') )):   
            fit_residual = (value_data - value_fit) / err_tot_data
            err_fit_residual = 1
        else:
            fit_residual = 0
            err_fit_residual = 0

        ## Fill histo with residuals

        hist_fit_residual_vsMass_.SetBinContent(bin+1,fit_residual)
        hist_fit_residual_vsMass_.SetBinError(bin+1,err_fit_residual)

        ## Chi2

        chi2_FullRangeAll += pow(fit_residual,2)
        N_FullRangeAll += 1
        if plotRegion=='Full' or (plotRegion=='Low,High' and (xbinLow >= workspace_.var('mjj').getMin() and xbinHigh<=workspace_.var('mjj').getMax())):
            #print '%i: obs %.0f, exp %.2f, chi2 %.2f'%(bin, value_data* binWidth_current * lumi, value_fit* binWidth_current * lumi, pow(fit_residual,2))
            chi2_PlotRangeAll += pow(fit_residual,2)
            N_PlotRangeAll += 1
            if (value_data > 0):
                chi2_PlotRangeNonZero += pow(fit_residual,2)
                N_PlotRangeNonZero += 1
                if(value_data * binWidth_current * lumi > MinNumEvents):
                    chi2_PlotRangeMinNumEvents += pow(fit_residual,2)
                    N_PlotRangeMinNumEvents += 1
    
    #==================
    # Calculate chi2/ndf
    #==================

    # ndf
    ndf_FullRangeAll = N_FullRangeAll - nParFit    
    ndf_PlotRangeAll = N_PlotRangeAll - nParFit    
    ndf_PlotRangeNonZero = N_PlotRangeNonZero - nParFit    
    ndf_PlotRangeMinNumEvents = N_PlotRangeMinNumEvents - nParFit    

    chi2_ndf_FullRangeAll = chi2_FullRangeAll / ndf_FullRangeAll
    chi2_ndf_PlotRangeAll = chi2_PlotRangeAll / ndf_PlotRangeAll
    chi2_ndf_PlotRangeNonZero = chi2_PlotRangeNonZero / ndf_PlotRangeNonZero
    chi2_ndf_PlotRangeMinNumEvents = chi2_PlotRangeMinNumEvents / ndf_PlotRangeMinNumEvents

    return [chi2_FullRangeAll, ndf_FullRangeAll, chi2_PlotRangeAll, ndf_PlotRangeAll, chi2_PlotRangeNonZero, ndf_PlotRangeNonZero, chi2_PlotRangeMinNumEvents, ndf_PlotRangeMinNumEvents]


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store cards")
    parser.add_option('-l','--lumi',dest="lumi", default=1.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('-b','--box',dest="box", default="CaloDijet",type="string",
                  help="box name")
    parser.add_option('--no-fit',dest="noFit",default=False,action='store_true',
                  help="Turn off fit (useful for visualizing initial parameters)")
    parser.add_option('--fit-region',dest="fitRegion",default="Full",type="string",
                  help="Fit region")
    parser.add_option('--plot-region',dest="plotRegion",default="Full",type="string",
                  help="Plot region")
    parser.add_option('-i','--input-fit-file',dest="inputFitFile", default=None,type="string",
                  help="input fit file")
    parser.add_option('-w','--weight',dest="useWeight",default=False,action='store_true',
                  help="use weight")
    parser.add_option('-s','--signal',dest="signalFileName", default=None,type="string",
                  help="input dataset file for signal pdf")
    parser.add_option('-m','--model',dest="model", default="gg",type="string",
                  help="signal model")
    parser.add_option('--mass',dest="mass", default=750,type="float",
                  help="mgluino")
    parser.add_option('--xsec',dest="xsec", default=1,type="float",
                  help="cross section in pb")


    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)
    
    box = options.box
    lumi = options.lumi
    noFit = options.noFit
    fitRegion = options.fitRegion
    plotRegion = options.plotRegion
    
    myTH1 = None
    for f in args:
        if f.lower().endswith('.root'):
            rootFile = rt.TFile(f)
            myTH1 = rootFile.Get('h_mjj_HLTpass_HT250_1GeVbin')
    if myTH1 is None:
        print "give a root file as input"
  
    w = rt.RooWorkspace("w"+box)

    paramNames, bkgs = initializeWorkspace(w,cfg,box)

        
    if options.inputFitFile is not None:
        inputRootFile = rt.TFile.Open(options.inputFitFile,"r")
        wIn = inputRootFile.Get("w"+box).Clone("wIn"+box)            
        if wIn.obj("fitresult_extDijetPdf_data_obs") != None:
            frIn = wIn.obj("fitresult_extDijetPdf_data_obs")
        elif wIn.obj("nll_extDijetPdf_data_obs") != None:
            frIn = wIn.obj("nll_extDijetPdf_data_obs")
        elif wIn.obj("fitresult_extDijetPdf_data_obs_with_constr") != None:
            frIn = wIn.obj("fitresult_extDijetPdf_data_obs_with_constr")
        elif wIn.obj("nll_extDijetPdf_data_obs_with_constr") != None:
            frIn = wIn.obj("nll_extDijetPdf_data_obs_with_constr")
                        
        print "restoring parameters from fit"
        frIn.Print("V")
        for p in rootTools.RootIterator.RootIterator(frIn.floatParsFinal()):
            w.var(p.GetName()).setVal(p.getVal())
            w.var(p.GetName()).setError(p.getError())

    
    x = array('d', cfg.getBinning(box)[0]) # mjj binning
    
    th1x = w.var('th1x')
    nBins = (len(x)-1)
    th1x.setBins(nBins)

    # get signal histo if any
    signalHistos = []
    signalHistosOriginal = []
    signalHistosRebin = []
    if options.signalFileName is not None:
        signalFile = rt.TFile.Open(options.signalFileName)
        names = [k.GetName() for k in signalFile.GetListOfKeys()]
        for name in names:
            d = signalFile.Get(name)
            if isinstance(d, rt.TH1):
                if name=='h_%s_%i'%(options.model,options.mass):
                    d.Scale(options.xsec*lumi/d.Integral())
                    d.Rebin(len(x)-1,name+'_rebin',x)
                    d_rebin = rt.gDirectory.Get(name+'_rebin')
                    d_rebin.SetDirectory(0)
    
                    signalHistosOriginal.append(d)
                    signalHistosRebin.append(d_rebin)
    
                    d_th1x = convertToTh1xHist(d_rebin)
                    signalHistos.append(d_th1x)
    
    sideband = convertSideband(fitRegion,w,x)
    plotband = convertSideband(plotRegion,w,x)

    extDijetPdf = w.pdf('extDijetPdf')

    myTH1.Rebin(len(x)-1,'data_obs_rebin',x)
    myRebinnedTH1 = rt.gDirectory.Get('data_obs_rebin')
    myRebinnedTH1.SetDirectory(0)
    
    myRealTH1 = convertToTh1xHist(myRebinnedTH1)        
    
    dataHist = rt.RooDataHist("data_obs","data_obs",rt.RooArgList(th1x), rt.RooFit.Import(myRealTH1))
    dataHist.Print('v')
    
    
    rootTools.Utils.importToWS(w,dataHist)

    if noFit and options.inputFitFile is not None:
        fr = frIn
        fr.Print('v')    
        rootTools.Utils.importToWS(w,fr)
    elif noFit:
        fr = rt.RooFitResult()
    else:
        fr = binnedFit(extDijetPdf,dataHist,sideband,options.useWeight)
        total = extDijetPdf.expectedEvents(rt.RooArgSet(th1x))
            
        fr.Print('v')    
        rootTools.Utils.importToWS(w,fr)
        

    asimov = extDijetPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('central'),rt.RooFit.Asimov())
        
    opt = [rt.RooFit.CutRange(myRange) for myRange in plotband.split(',')]
    asimov_reduce = asimov.reduce(opt[0])
    dataHist_reduce = dataHist.reduce(opt[0])
    for iOpt in range(1,len(opt)):
        asimov_reduce.add(asimov.reduce(opt[iOpt]))
        dataHist_reduce.add(dataHist.reduce(opt[iOpt]))
    rt.TH1D.SetDefaultSumw2()
    
    # start writing output
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetOptTitle(0)
    c = rt.TCanvas('c','c',600,700)
    rootFile = rt.TFile.Open(options.outDir + '/' + 'Plots_%s'%box + '.root','recreate')
    tdirectory = rootFile.GetDirectory(options.outDir)
    if tdirectory==None:
        print "making directory"
        rootFile.mkdir(options.outDir)
        tdirectory = rootFile.GetDirectory(options.outDir)
        tdirectory.Print('v')
        
    h_th1x = asimov_reduce.createHistogram('h_th1x',th1x)
    h_data_th1x = dataHist_reduce.createHistogram('h_data_th1x',th1x)
    
    boxLabel = "%s %s Fit" % (box,fitRegion)
    plotLabel = "%s Projection" % (plotRegion)


    
    background_pdf = w.pdf('%s_bkg_unbin'%box)
    background= background_pdf.asTF(rt.RooArgList(w.var('mjj')),rt.RooArgList(w.var('p0')))
    int_b = background.Integral(w.var('mjj').getMin(),w.var('mjj').getMax())
    p0_b = w.var('Ntot_bkg').getVal() / (int_b * lumi)
    background.SetParameter(0,p0_b)
    
    g_data = rt.TGraphAsymmErrors(myRebinnedTH1)
    
    alpha = 1-0.6827
    for i in range(0,g_data.GetN()):
        N = g_data.GetY()[i]
        binWidth = g_data.GetEXlow()[i] + g_data.GetEXhigh()[i]
        L = 0
        if N!=0:
            L = rt.Math.gamma_quantile(alpha/2,N,1.)
        U = rt.Math.gamma_quantile_c(alpha/2,N+1,1)

        #g_data.SetPointEYlow(i, (N-L))
        #g_data.SetPointEYhigh(i, (U-N))
        #g_data.SetPoint(i, g_data.GetX()[i], N)
        g_data.SetPointEYlow(i, (N-L)/(binWidth * lumi))
        g_data.SetPointEYhigh(i, (U-N)/(binWidth * lumi))
        g_data.SetPoint(i, g_data.GetX()[i], N/(binWidth * lumi))
        
        if plotRegion=='Low,High' and (g_data.GetX()[i]>w.var('mjj').getMin('Blind') and g_data.GetX()[i]<w.var('mjj').getMax('Blind')):
            g_data.SetPointEYlow(i, 0)
            g_data.SetPointEYhigh(i, 0)
            g_data.SetPoint(i, g_data.GetX()[i], 0)

            
    h_background = convertFunctionToHisto(background,"h_background",len(x)-1,x)
    h_fit_residual_vs_mass = rt.TH1D("h_fit_residual_vs_mass","h_fit_residual_vs_mass",len(x)-1,x)
    list_chi2AndNdf_background = calculateChi2AndFillResiduals(g_data,h_background,h_fit_residual_vs_mass,w,0)

    g_data.SetMarkerStyle(20)
    g_data.SetMarkerSize(0.9)
    g_data.SetLineColor(rt.kBlack)
    g_data_clone = g_data.Clone('g_data_clone')
    g_data_clone.SetMarkerSize(0)
    myRebinnedTH1.SetLineColor(rt.kWhite)
    myRebinnedTH1.SetMarkerSize(0)


    c.Divide(1,2,0,0,0)
    
    pad_1 = c.GetPad(1)
    pad_1.SetPad(0.01,0.36,0.99,0.98)
    pad_1.SetLogy()
    pad_1.SetRightMargin(0.05)
    pad_1.SetTopMargin(0.05)
    pad_1.SetLeftMargin(0.175)
    pad_1.SetFillColor(0)
    pad_1.SetBorderMode(0)
    pad_1.SetFrameFillStyle(0)
    pad_1.SetFrameBorderMode(0)
    
    pad_2 = c.GetPad(2)
    pad_2.SetLeftMargin(0.175)
    pad_2.SetPad(0.01,0.02,0.99,0.37)
    pad_2.SetBottomMargin(0.35)
    pad_2.SetRightMargin(0.05)
    pad_2.SetGridx()
    pad_2.SetGridy()

    pad_1.cd()
    
    myRebinnedDensityTH1 = myRebinnedTH1.Clone('data_obs_density')
    for i in range(1,nBins+1):
        myRebinnedDensityTH1.SetBinContent(i, myRebinnedTH1.GetBinContent(i)/ myRebinnedTH1.GetBinWidth(i))
        myRebinnedDensityTH1.SetBinError(i, myRebinnedTH1.GetBinError(i)/ myRebinnedTH1.GetBinWidth(i))        
        if plotRegion=='Low,High' and (myRebinnedDensityTH1.GetXaxis().GetBinCenter(i)>w.var('mjj').getMin('Blind') and myRebinnedDensityTH1.GetXaxis().GetBinCenter(i)<w.var('mjj').getMax('Blind')):
            myRebinnedDensityTH1.SetBinContent(i,0)
            myRebinnedDensityTH1.SetBinError(i,0)
    myRebinnedDensityTH1.GetXaxis().SetRangeUser(w.var('mjj').getMin(),w.var('mjj').getMax())
    myRebinnedDensityTH1.GetYaxis().SetTitle('d#sigma / dm_{jj} [pb / GeV]')
    myRebinnedDensityTH1.GetYaxis().SetTitleOffset(1)
    myRebinnedDensityTH1.GetYaxis().SetTitleSize(0.07)
    myRebinnedDensityTH1.GetYaxis().SetLabelSize(0.05)
    
    myRebinnedDensityTH1.SetLineColor(rt.kWhite)
    myRebinnedDensityTH1.SetMarkerColor(rt.kWhite)
    myRebinnedDensityTH1.SetLineWidth(0)
    myRebinnedDensityTH1.SetMaximum(1e3)
    myRebinnedDensityTH1.SetMinimum(2e-5)
    #myRebinnedDensityTH1.SetMinimum(2e-10)
    myRebinnedDensityTH1.Draw("pe")    
    g_data_clone.Draw("pezsame")
    background.Draw("csame")
    g_data.Draw("pezsame")

    if options.signalFileName is not None:
        sigHist = signalHistosRebin[0]
        
        g_signal = rt.TGraphAsymmErrors(sigHist)
        g_signal.SetLineColor(rt.kBlue)
        g_signal.SetLineWidth(3)
        g_signal.SetLineStyle(2)

        lastX = 0
        lastY = 0
        for i in range(0,g_signal.GetN()):
            N = g_signal.GetY()[i]
            binWidth = g_signal.GetEXlow()[i] + g_signal.GetEXhigh()[i]
            g_signal.SetPoint(i, g_signal.GetX()[i], N/(binWidth * lumi))
            g_signal.SetPointEYlow(i, 0)
            g_signal.SetPointEYhigh(i, 0)
            if g_signal.GetX()[i]>options.mass*1.5:
                g_signal.SetPoint(i,lastX,lastY)
            else:                
                lastX = g_signal.GetX()[i]
                lastY = g_signal.GetY()[i]
        #sys.exit()
        g_signal.Draw("lxsame")

    
    rt.gPad.SetLogy()
    
    l = rt.TLatex()
    l.SetTextAlign(11)
    l.SetTextSize(0.05)
    l.SetTextFont(42)
    l.SetNDC()
    l.DrawLatex(0.7,0.96,"%i pb^{-1} (%i TeV)"%(lumi,w.var('sqrts').getVal()/1000.))
    l.SetTextFont(62)
    l.SetTextSize(0.06)
    l.DrawLatex(0.2,0.96,"CMS")
    l.SetTextFont(52)
    l.SetTextSize(0.05)
    l.DrawLatex(0.3,0.96,"Preliminary")
    if options.signalFileName!=None:
        leg = rt.TLegend(0.6,0.52,0.89,0.88)
    else:        
        leg = rt.TLegend(0.7,0.7,0.89,0.88)
    leg.SetTextFont(42)
    leg.SetFillColor(rt.kWhite)
    leg.SetFillStyle(0)
    leg.SetLineWidth(0)
    leg.SetLineColor(rt.kWhite)
    leg.AddEntry(g_data,"Data","pe")
    leg.AddEntry(background,"Fit","l")
    if options.signalFileName!=None:
        leg.AddEntry(g_signal,"%s (%i GeV)"%(options.model,options.mass),"l")
        leg.AddEntry(None,"%.1f pb"%(options.xsec),"")
    leg.Draw()
    
    pave_sel = rt.TPaveText(0.2,0.03,0.5,0.25,"NDC")
    pave_sel.SetFillColor(0)
    pave_sel.SetBorderSize(0)
    pave_sel.SetFillStyle(0)
    pave_sel.SetTextFont(42)
    pave_sel.SetTextSize(0.045)
    pave_sel.SetTextAlign(11)
    pave_sel.AddText("#chi^{{2}} / ndf = {0:.1f} / {1:d} = {2:.1f}".format(
                          list_chi2AndNdf_background[4], list_chi2AndNdf_background[5],
                          list_chi2AndNdf_background[4]/list_chi2AndNdf_background[5]))
    pave_sel.AddText("Wide Jets")
    pave_sel.AddText("%i < m_{jj} < %i GeV"%(w.var('mjj').getMin(),w.var('mjj').getMax()))
    pave_sel.AddText("|#eta| < 2.5, |#Delta#eta| < 1.3")
    pave_sel.Draw("SAME")
    
    list_parameter = [p0_b, p0_b*(w.var('Ntot_bkg').getErrorHi() - w.var('Ntot_bkg').getErrorLo())/(2.0*w.var('Ntot_bkg').getVal()),                      
                      w.var('p1').getVal(), (w.var('p1').getErrorHi() - w.var('p1').getErrorLo())/2.0,
                      w.var('p2').getVal(), (w.var('p2').getErrorHi() - w.var('p2').getErrorLo())/2.0,
                      w.var('p3').getVal(), (w.var('p3').getErrorHi() - w.var('p3').getErrorLo())/2.0,
                      w.var('meff').getVal(), (w.var('meff').getErrorHi() - w.var('meff').getErrorLo())/2.0,
                      w.var('seff').getVal(), (w.var('seff').getErrorHi() - w.var('seff').getErrorLo())/2.0]


    pave_param = rt.TPaveText(0.55,0.03,0.9,0.25,"NDC")
    pave_param.SetTextFont(42)
    pave_param.SetFillColor(0)
    pave_param.SetBorderSize(0)
    pave_param.SetFillStyle(0)
    pave_param.SetTextAlign(11)
    pave_param.SetTextSize(0.045)
    pave_param.AddText("p_{0}"+" = {0:.2g} #pm {1:.2g}".format(list_parameter[0], list_parameter[1]))
    pave_param.AddText("p_{1}"+" = {0:.2f} #pm {1:.2f}".format(list_parameter[2], list_parameter[3]))
    pave_param.AddText("p_{2}"+" = {0:.2f} #pm {1:.2f}".format(list_parameter[4], list_parameter[5]))
    pave_param.AddText("p_{3}"+" = {0:.2f} #pm {1:.2f}".format(list_parameter[6], list_parameter[7]))
    if w.var('meff').getVal()>0 and w.var('seff').getVal()>0:
        pave_param.AddText("m_{eff}"+" = {0:.2f} #pm {1:.2f}".format(list_parameter[8], list_parameter[9]))
        pave_param.AddText("#sigma_{eff}"+" = {0:.2f} #pm {1:.2f}".format(list_parameter[10], list_parameter[11]))
    pave_param.Draw("SAME")    
    
    pad_1.Update()

    pad_2.cd()
    
    h_fit_residual_vs_mass.GetXaxis().SetRangeUser(w.var('mjj').getMin(),w.var('mjj').getMax())
    h_fit_residual_vs_mass.GetYaxis().SetRangeUser(-3.5,3.5)
    h_fit_residual_vs_mass.GetYaxis().SetNdivisions(210,True)
    h_fit_residual_vs_mass.SetLineWidth(1)
    h_fit_residual_vs_mass.SetFillColor(rt.kRed)
    h_fit_residual_vs_mass.SetLineColor(rt.kBlack)
    
    h_fit_residual_vs_mass.GetYaxis().SetTitleSize(2*0.06)
    h_fit_residual_vs_mass.GetYaxis().SetLabelSize(2*0.05)
    h_fit_residual_vs_mass.GetYaxis().SetTitleOffset(0.5)
    h_fit_residual_vs_mass.GetYaxis().SetTitle('#frac{(Data-Fit)}{#sigma_{Data}}')
        
    h_fit_residual_vs_mass.GetXaxis().SetTitleSize(2*0.06)
    h_fit_residual_vs_mass.GetXaxis().SetLabelSize(2*0.05)
    h_fit_residual_vs_mass.GetXaxis().SetTitle('m_{jj} [GeV]')
    
    
    h_fit_residual_vs_mass.Draw("histsame")
    
        
    if options.signalFileName is not None:
        sigHistResidual = sigHist.Clone(sigHist.GetName()+"_residual")
        sigHistResidual.SetLineColor(rt.kBlue)
        sigHistResidual.SetLineWidth(2)
        sigHistResidual.SetLineStyle(2)
        for bin in range (0,g_data.GetN()):
            value_data = g_data.GetY()[bin]
            err_tot_data = g_data.GetEYhigh()[bin]
            binWidth = g_data.GetEXlow()[i] + g_data.GetEXhigh()[i]
            value_signal = sigHist.GetBinContent(bin+1)/(binWidth*lumi)
        
            ## Signal residuals
            if err_tot_data>0:                
                sig_residual = (value_signal) / err_tot_data
            else:
                sig_residual = 0                                
    
            ## Fill histo with residuals
            sigHistResidual.SetBinContent(bin+1,sig_residual)
        sigHistResidual.Draw("histsame")

    
    c.Print(options.outDir+"/fit_mjj_%s_%s.pdf"%(fitRegion.replace(',','_'),box))
    c.Print(options.outDir+"/fit_mjj_%s_%s.C"%(fitRegion.replace(',','_'),box))
    tdirectory.cd()
    c.Write()
    

    outFileName = "DijetFitResults_%s.root"%(box)
    outFile = rt.TFile.Open(options.outDir+"/"+outFileName,'recreate')
    outFile.cd()
    w.Write()
    outFile.Close()
