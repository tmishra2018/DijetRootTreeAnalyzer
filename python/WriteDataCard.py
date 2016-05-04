from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
from array import *
import os
import sys

def fixPars(w, label, doFix=True, setVal=None):
    parSet = w.allVars()
    for par in rootTools.RootIterator.RootIterator(parSet):
        if label in par.GetName():
            par.setConstant(doFix)
            if setVal is not None: par.setVal(setVal)

def initializeWorkspace(w,cfg,box,scaleFactor=1.,penalty=False,x=None,emptyHist1D=None):
    
    if x is None:
        x = array('d', cfg.getBinning(box)[0]) # mjj binning
    nBins = len(x)-1
    maxBins = nBins
    
    variables = cfg.getVariablesRange(box, "variables",w)
    w.var('th1x').setBins(maxBins)
    parameters = cfg.getVariables(box, "combine_parameters")
    paramNames = []
    for parameter in parameters:
        if penalty and '_norm' in parameter:
            continue
        w.factory(parameter)
        paramName = parameter.split('[')[0]
        if paramName not in ['sqrts','meff','seff','p0']:
            paramNames.append(paramName)
            w.var(paramName).setConstant(False)
            
        
        # float normalization parameters
        fixPars(w,"Ntot",False)
        
        # fix Gaussian constraint parameters
        fixPars(w,"In")
        fixPars(w,"Mean")
        fixPars(w,"Sigma")

        # fix center of mass energy, trigger turn-on, and p0                                                                                   
        for myvar in ['sqrts','meff','seff','p0']:
            fixPars(w,myvar)

        
        
    if emptyHist1D==None:
        emptyHist1D = rt.TH1D("emptyHist1D","emptyHist1D",len(x)-1,x)
        iBinX = -1
        for ix in range(1,len(x)):
            iBinX+=1
            emptyHist1D.SetBinContent(ix,1)
            emptyHist1D.SetBinError(ix,0)
        
    
    commands = cfg.getVariables(box, "combine_pdfs")
    bkgs = []
    for command in commands:
        lower = command.lower()
        if lower.find('sum::')!=-1 or lower.find('prod::')!=-1 or lower.find('expr::')!=-1 or lower.find('roogaussian::')!=-1:
            w.factory(command)
        else:
            myclass = command.split('::')[0]
            remaining = command.split('::')[1]
            name = remaining.split('(')[0]
            altname = '_'.join(reversed(name.split('_')))
            mytuple = remaining.replace(name,'').replace('(','').replace(')','')
            mylist = mytuple.split(',')
            arglist = [name, name]
            for myvar in mylist:
                if w.var(myvar)!=None:
                    arglist.append(w.var(myvar))
                elif w.function(myvar)!=None:
                    arglist.append(w.function(myvar))
                        
            args = tuple(arglist)
            pdf = getattr(rt,myclass)(*args)
            if hasattr(pdf,'setTH1Binning'):
                pdf.setTH1Binning(emptyHist1D)
            rootTools.Utils.importToWS(w,pdf)
            bkg = name.split("_")
            if box in bkg: bkg.remove(box)
            bkgs.append("_".join(bkg))
            
    w.Print('v')
    return paramNames, bkgs


def writeDataCard(box,model,txtfileName,bkgs,paramNames,w,penalty,fixed,shapes=[]):
        obsRate = w.data("data_obs").sumEntries()
        nBkgd = len(bkgs)
        rootFileName = txtfileName.replace('.txt','.root')
        signals = len(model.split('p'))
        if signals>1:
                rates = [w.data("%s_%s"%(box,sig)).sumEntries() for sig in model.split('p')]
                processes = ["%s_%s"%(box,sig) for sig in model.split('p')]
                lumiErrs = [1.027 for sig in model.split('p')]
        else:
                rates = [w.data("%s_%s"%(box,model)).sumEntries()]
                processes = ["%s_%s"%(box,model)]
                lumiErrs = [1.027]
        rates.extend([w.var('Ntot_%s'%(bkg)).getVal() for bkg in bkgs])
        processes.extend(["%s_%s"%(box,bkg) for bkg in bkgs])
        lumiErrs.extend([1.00 for bkg in bkgs])
        divider = "------------------------------------------------------------\n"
        datacard = "imax 1 number of channels\n" + \
                   "jmax %i number of processes minus 1\n"%(nBkgd+signals-1) + \
                   "kmax * number of nuisance parameters\n" + \
                   divider + \
                   "observation	%.3f\n"%obsRate + \
                   divider + \
                   "shapes * * %s w%s:$PROCESS w%s:$PROCESS_$SYSTEMATIC\n"%(rootFileName,box,box) + \
                   divider
        binString = "bin"
        processString = "process"
        processNumberString = "process"
        rateString = "rate"
        lumiString = "lumi\tlnN"
        for i in range(0,len(bkgs)+signals):
            binString +="\t%s"%box
            processString += "\t%s"%processes[i]
            processNumberString += "\t%i"%(i-signals+1)
            rateString += "\t%.3f" %rates[i]
            lumiString += "\t%.3f"%lumiErrs[i]
        binString+="\n"; processString+="\n"; processNumberString+="\n"; rateString +="\n"; lumiString+="\n"
        datacard+=binString+processString+processNumberString+rateString+divider
        # now nuisances
        datacard+=lumiString
        for shape in shapes:
            shapeString = '%s\tshape\t'%shape
            for sig in range(0,signals):
                shapeString += '\t1.0'
            for i in range(0,len(bkgs)):
                shapeString += '\t-'
            shapeString += '\n'
            datacard+=shapeString
        for paramName in paramNames:
            if fixed:
                fixPars(w,paramName)    
            elif 'Mean' in paramName or 'Sigma' in paramName:
                fixPars(w,paramName)           
            elif penalty:                    
                mean = w.var(paramName).getVal()
                sigma = w.var(paramName).getError()                
                if "Ntot" in paramName:                    
                    effectString = ''
                    for sig in range(0,signals):
                        effectString += "\t1.0"           
                    for bkg in bkgs:
                        if bkg in paramName:
                            effectString += "\t%.3f"%(1.0+sigma/mean)                            
                        else:
                            effectString += "\t1.0"                    
                    datacard += "%s\tlnN%s\n"%(paramName.replace("Ntot","Norm"),effectString)
                else:
                    datacard += "%s\tparam\t%e\t%e\n"%(paramName,mean,sigma)
                         
            else:
                if "Ntot" in paramName:
                    continue
                else:
                    datacard += "%s\tflatParam\n"%(paramName)
            
        txtfile = open(txtfileName,"w")
        txtfile.write(datacard)
        txtfile.close()

def convertToTh1xHist(hist):
    
    hist_th1x = rt.TH1D(hist.GetName()+'_th1x',hist.GetName()+'_th1x',hist.GetNbinsX(),0,hist.GetNbinsX())
    for i in range(1,hist.GetNbinsX()+1):
        hist_th1x.SetBinContent(i,hist.GetBinContent(i))
        hist_th1x.SetBinError(i,hist.GetBinError(i))

    return hist_th1x


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store cards")
    parser.add_option('-l','--lumi',dest="lumi", default=1.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('--jes',dest="jesUnc", default=0.05,type="float",
                  help="jes uncertainty, default = 0.05")
    parser.add_option('--jer',dest="jerUnc", default=0.2,type="float",
                  help="jer uncertainty, default = 0.2")
    parser.add_option('-b','--box',dest="box", default="CaloDijet",type="string",
                  help="box name")
    parser.add_option('--asimov',dest="asimov",default=False,action='store_true',
                  help="replace real data with asimov dataset from input fit result")
    parser.add_option('--penalty',dest="penalty",default=False,action='store_true',
                  help="penalty terms on background shape + norm parameters from input fit result")
    parser.add_option('--fixed',dest="fixed",default=False,action='store_true',
                  help="fixed background shape + norm parameters")
    parser.add_option('-i','--input-fit-file',dest="inputFitFile", default=None,type="string",
                  help="input fit file")
    parser.add_option('-m','--model',dest="model", default="gg",type="string",
                  help="signal model name")
    parser.add_option('--mass',dest="mass", default=750,type="float",
                  help="mass of resonance")
    parser.add_option('--xsec',dest="xsec", default=1,type="float",
                  help="xsec of resonance")
    parser.add_option('--no-signal-sys',dest="noSignalSys",default=False,action='store_true',
                  help="no signal shape systematic uncertainties")

    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)

    box = options.box
    lumi = options.lumi
    
    signalXsec = options.xsec

    signalFileName = ''
    model = options.model
    massPoint = options.mass

    jesUnc = options.jesUnc
    jerUnc = options.jerUnc

    myTH1 = None
    for f in args:
        if f.lower().endswith('.root'):
            if f.lower().find('resonanceshapes')!=-1:
                signalFileName = f
            else:
                rootFile = rt.TFile(f)
                myTH1 = rootFile.Get('h_mjj_HLTpass_HT250_1GeVbin')

    w = rt.RooWorkspace("w"+box)
    
    paramNames, bkgs = initializeWorkspace(w,cfg,box,scaleFactor=1,penalty=options.penalty)
    
    
    th1x = w.var('th1x')
    
    if myTH1 is None:
        print "give a background root file as input"        
    
    x = array('d', cfg.getBinning(box)[0]) # mjj binning
        
    myTH1.Rebin(len(x)-1,'data_obs_rebin',x)
    myRebinnedTH1 = rt.gDirectory.Get('data_obs_rebin')
    myRebinnedTH1.SetDirectory(0)

    myRealTH1 = convertToTh1xHist(myRebinnedTH1)
    

    if options.inputFitFile is not None:
        inputRootFile = rt.TFile.Open(options.inputFitFile,"r")
        wIn = inputRootFile.Get("w"+box).Clone("wIn"+box)            
        if wIn.obj("fitresult_extDijetPdf_data_obs") != None:
            frIn = wIn.obj("fitresult_extDijetPdf_data_obs")
        elif wIn.obj("nll_extDijetPdf_data_obs") != None:
            frIn = wIn.obj("nll_extDijetPdf_data_obs")
        elif wIn.obj("fitresult_extDijetPdf_data_obs_with_constr") != None:
            fr = wIn.obj("fitresult_extDijetPdf_data_obs_with_constr")
        elif wIn.obj("nll_extDijetPdf_data_obs_with_constr") != None:
            frIn = wIn.obj("nll_extDijetPdf_data_obs_with_constr")
        print "restoring parameters from fit"
        frIn.Print("V")
        for p in rootTools.RootIterator.RootIterator(frIn.floatParsFinal()):
            w.var(p.GetName()).setVal(p.getVal())
            w.var(p.GetName()).setError(p.getError())
            if "Ntot" in p.GetName():
                normNameList = p.GetName().replace("Ntot_","").split("_")
                normNameList.reverse()
                normNameList.append("norm")
                normName = "_".join(normNameList)
                #w.var(normName).setError(w.var(p.GetName()).getError()/w.var(p.GetName()).getVal())
                
    if options.asimov:
        asimov = w.pdf('extDijetPdf').generateBinned(rt.RooArgSet(th1x),rt.RooFit.Asimov())
        asimov.SetName('data_obs')
        asimov.SetTitle('data_obs')
        dataHist = asimov
    else:
        dataHist = rt.RooDataHist("data_obs", "data_obs", rt.RooArgList(th1x), rt.RooFit.Import(myRealTH1))
        
    rootTools.Utils.importToWS(w,dataHist)
    
    signalHistos = []
    signalFile = rt.TFile.Open(signalFileName)
    names = [k.GetName() for k in signalFile.GetListOfKeys()]
    for name in names:
        d = signalFile.Get(name)
        if isinstance(d, rt.TH1):
            #d.SetDirectory(rt.gROOT)
            if name=='h_%s_%i'%(model,massPoint):
                d.Scale(signalXsec*lumi/d.Integral())
                d.Rebin(len(x)-1,name+'_rebin',x)
                d_rebin = rt.gDirectory.Get(name+'_rebin')
                d_rebin.SetDirectory(0)

                d_th1x = convertToTh1xHist(d_rebin)
                signalHistos.append(d_th1x)
                
                sigDataHist = rt.RooDataHist('%s_%s'%(box,model),'%s_%s'%(box,model), rt.RooArgList(th1x), rt.RooFit.Import(d_th1x))
                rootTools.Utils.importToWS(w,sigDataHist)

    if options.noSignalSys:
        shapes = []
    else:
        shapes = ['jes','jer']


    # JES and JER uncertainties
    hSig = signalHistos[0]
    sigCDF = rt.TGraph(hSig.GetNbinsX()+1)

    sigCDF.SetPoint(0,0.,0.)
    integral = 0.
    for i in range(1, hSig.GetNbinsX()+1):
        x = hSig.GetXaxis().GetBinLowEdge(i+1)
        integral = integral + hSig.GetBinContent(i)
        sigCDF.SetPoint(i,x,integral)
        
    for shape in shapes:
        hUp = hSig.Clone(hSig.GetName()+'_'+shape+'Up')
        hDown = hSig.Clone(hSig.GetName()+'_'+shape+'Down')

        # produce JES/JER signal shapes
        for i in range(1, hSig.GetNbinsX()+1):
            xLow = hSig.GetXaxis().GetBinLowEdge(i)
            xUp = hSig.GetXaxis().GetBinLowEdge(i+1)
            if shape=='jes':
                jes = 1. - jesUnc
                xLowPrime = jes*xLow
                xUpPrime = jes*xUp
            elif shape=='jer':      
                jer = 1. - jerUnc              
                xLowPrime = jer*(xLow-float(massPoint)+float(massPoint))
                xUpPrime = jer*(xUp-float(massPoint)+float(massPoint))
            hUp.SetBinContent(i, sigCDF.Eval(xUpPrime) - sigCDF.Eval(xLowPrime))
            if shape=='jes':
                jes = 1. + jesUnc
                xLowPrime = jes*xLow
                xUpPrime = jes*xUp
            elif shape=='jer':                    
                jer = 1. + jerUnc
                xLowPrime = jer*(xLow-float(massPoint)+float(massPoint))
                xUpPrime = jer*(xUp-float(massPoint)+float(massPoint))
            hDown.SetBinContent(i, sigCDF.Eval(xUpPrime) - sigCDF.Eval(xLowPrime))
                
        hUp_DataHist = rt.RooDataHist('%s_%s_%sUp'%(box,model,shape),'%s_%s_%sUp'%(box,model,shape),rt.RooArgList(th1x),hUp)
        hDown_DataHist = rt.RooDataHist('%s_%s_%sDown'%(box,model,shape),'%s_%s_%sDown'%(box,model,shape),rt.RooArgList(th1x),hDown)
        
        rootTools.Utils.importToWS(w,hUp_DataHist)
        rootTools.Utils.importToWS(w,hDown_DataHist)            

            
    outFile = 'dijet_combine_%s_%i_lumi-%.3f_%s.root'%(model,massPoint,lumi/1000.,box)
    
    outputFile = rt.TFile.Open(options.outDir+"/"+outFile,"recreate")
    writeDataCard(box,model,options.outDir+"/"+outFile.replace(".root",".txt"),bkgs,paramNames,w,options.penalty,options.fixed,shapes=shapes)
    w.Write()
    os.system("cat %s"%options.outDir+"/"+outFile.replace(".root",".txt"))
