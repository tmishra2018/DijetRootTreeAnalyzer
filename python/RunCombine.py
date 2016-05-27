from optparse import OptionParser
import os
import ROOT as rt
from array import *
from framework import Config
import sys
import glob
import rootTools
import time

NSIGMA = 3.0

def massIterable(massList):    
    if len(massList.split(','))==1:
        massIterableList = [options.mass]
    else:
        massIterableList = list(eval(options.mass))
    return massIterableList
        
def exec_me(command,dryRun=True):
    print command
    if not dryRun: os.system(command)
        
def writeBashScript(options,massPoint):
    lumi = float(options.lumi)
    submitDir = options.outDir

    
    signalSys = ''
    if options.noSignalSys:
        signalSys = '--no-signal-sys'
        
    penaltyString = ''
    if options.penalty:
        penaltyString = '--penalty'

    decoString = ''
    if options.deco:
        decoString  ='--deco'
        
    bayesString = ''
    if options.bayes:
        bayesString  ='--bayes'
        
    toyString = ''
    if options.toys>-1:
        toyString  ='--toys %i'%options.toys

    xsecString = '--xsec %f'%options.xsec
        
        
    # prepare the script to run
    outputname = submitDir+"/submit_"+options.model+"_"+massPoint+"_lumi-%.3f_"%(lumi)+options.box+".src"
        
    ffDir = submitDir+"/logs_"+options.model+"_"+massPoint+"_"+options.box
    user = os.environ['USER']
    pwd = os.environ['PWD']
        
    combineDir = "/afs/cern.ch/work/%s/%s/DIJET/Limits/%s/"%(user[0],user,options.model) # directory where combine output files will be copied
    cmsswBase = "/afs/cern.ch/work/%s/%s/DIJET/CMSSW_7_4_14"%(user[0],user) # directory where 'cmsenv' will be run (needs to have combine setup)
    
    script =  '#!/usr/bin/env bash -x\n'
    script += 'mkdir -p %s\n'%combineDir        
    script += 'echo $SHELL\n'
    script += 'pwd\n'
    script += 'cd %s/src/DijetRootTreeAnalyzer \n'%(cmsswBase)
    script += 'pwd\n'
    script += "export SCRAM_ARCH=slc6_amd64_gcc491\n"
    script += "export CMSSW_BASE=%s\n"%(cmsswBase)
    script += 'eval `scramv1 runtime -sh`\n'
    script += 'cd - \n'
    script += "export TWD=${PWD}/%s_%s_lumi-%.3f_%s\n"%(options.model,massPoint,lumi,options.box)
    script += "mkdir -p $TWD\n"
    script += "cd $TWD\n"
    script += 'pwd\n'
    script += 'git clone git@github.com:CMSDIJET/DijetRootTreeAnalyzer\n'
    script += 'cd DijetRootTreeAnalyzer\n'
    script += 'git checkout -b Limits %s\n'%(options.tag)
    script += 'mkdir -p %s\n'%submitDir
    script += 'wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/master/ResonanceShapes_%s_13TeV_CaloScouting_Spring15.root -P inputs/\n'%(options.model)
    for sys in ['JERUP','JERDOWN','JESUP','JESDOWN']:
        script += 'wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/master/ResonanceShapes_%s_13TeV_CaloScouting_Spring15_%s.root -P inputs/\n'%(options.model,sys)
    script += 'python python/RunCombine.py -i %s -m %s --mass %s -c %s --lumi %f -d %s -b %s %s %s --min-tol %e --min-strat %i --rMax %f %s %s %s %s\n'%(options.inputFitFile,
                                                                                                                                                         options.model,
                                                                                                                                                         massPoint,
                                                                                                                                                         options.config,
                                                                                                                                                         lumi,
                                                                                                                                                         submitDir,
                                                                                                                                                         options.box,
                                                                                                                                                         penaltyString,
                                                                                                                                                         signalSys,
                                                                                                                                                         options.min_tol,
                                                                                                                                                         options.min_strat,
                                                                                                                                                         options.rMax,
                                                                                                                                                         decoString,
                                                                                                                                                         bayesString,
                                                                                                                                                         toyString,
                                                                                                                                                         xsecString)
    script += 'cp %s/higgsCombine* %s/\n'%(submitDir,combineDir)
    script += 'cd ../..\n'
    script += 'rm -rf $TWD\n'
        
    outputfile = open(outputname,'w')
    outputfile.write(script)
    outputfile.close
    
    return outputname,ffDir

def submit_jobs(options,args):    
     
    for massPoint in massIterable(options.mass):

        outputname,ffDir = writeBashScript(options,massPoint)
                    
        pwd = os.environ['PWD']
        os.system("mkdir -p "+pwd+"/"+ffDir)
        os.system("echo bsub -q "+options.queue+" -o "+pwd+"/"+ffDir+"/log.log source "+pwd+"/"+outputname)      
        if not options.dryRun:
            time.sleep(3)
            os.system("bsub -q "+options.queue+" -o "+pwd+"/"+ffDir+"/log.log source "+pwd+"/"+outputname)
    
def main(options,args):    
    box = options.box

    signif = options.signif
    
    lumi = float(options.lumi)
    
    model = options.model

    paramDict = {}
    if options.inputFitFile is not None and options.bayes:
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
        elif wIn.obj("simNll") != None:
            frIn = wIn.obj("simNll")
        paramDict = {}
        for p in rootTools.RootIterator.RootIterator(frIn.floatParsFinal()):
            paramDict[p.GetName()] = [p.getVal(), p.getError()]
        print "grabbing parameter ranges +-%gsigma for bayesian"%NSIGMA

    signalSys = ''
    if options.noSignalSys:
        signalSys = '--no-signal-sys'
    else:
        signalSys  =   '--jesUp inputs/ResonanceShapes_%s_13TeV_CaloScouting_Spring15_JESUP.root --jesDown inputs/ResonanceShapes_%s_13TeV_CaloScouting_Spring15_JESDOWN.root'%(model,model)
        signalSys += ' --jerUp inputs/ResonanceShapes_%s_13TeV_CaloScouting_Spring15_JERUP.root --jerDown inputs/ResonanceShapes_%s_13TeV_CaloScouting_Spring15_JERDOWN.root'%(model,model)
    
    penaltyString = ''
    if options.penalty: penaltyString = '--penalty'

    xsecString = '--xsec %f'%(options.xsec)    

    signalDsName = 'inputs/ResonanceShapes_%s_13TeV_CaloScouting_Spring15.root'%model
    backgroundDsName = 'inputs/data_CaloScoutingHT_Run2015D_BiasCorrected_%s.root'%box

    blindString = ''
    if options.blind:
        blindString = '--noFitAsimov --run expected'

    sysString = ''
    if options.noSys:
        sysString = '-S 0 --freezeNuisances=shapeBkg_%s_bkg_%s__norm,p1,p2,p3,jes,jer'%(box,box)
    
    decoString = ''
    if options.deco:
        decoString  ='--deco'

    for massPoint in massIterable(options.mass):
                
        exec_me('python python/WriteDataCard.py -m %s --mass %s -i %s -l %f -c %s -b %s -d %s %s %s %s %s %s %s'%(model, massPoint, options.inputFitFile,1000*lumi,options.config,box,options.outDir,signalDsName,backgroundDsName,penaltyString,signalSys,xsecString,decoString),options.dryRun)

        if options.bayes:
            rRangeString =  '--setPhysicsModelParameterRanges '
            if options.deco:
                rRangeString += 'shapeBkg_%s_bkg_deco_%s__norm=%f,%f'%(box,box,1-NSIGMA*paramDict['Ntot_bkg'][1]/paramDict['Ntot_bkg'][0],1+NSIGMA*paramDict['Ntot_bkg'][1]/paramDict['Ntot_bkg'][0])
                rRangeString += ':deco_eig1=%f,%f'%(-1.0*NSIGMA,NSIGMA)
                rRangeString += ':deco_eig1=%f,%f'%(-1.0*NSIGMA,NSIGMA)
                rRangeString += ':deco_eig1=%f,%f'%(-1.0*NSIGMA,NSIGMA)
            else:
                rRangeString += 'shapeBkg_%s_bkg_%s__norm=%f,%f'%(box,box,1-NSIGMA*paramDict['Ntot_bkg'][1]/paramDict['Ntot_bkg'][0],1+NSIGMA*paramDict['Ntot_bkg'][1]/paramDict['Ntot_bkg'][0])
                rRangeString += ':p1=%f,%f'%(paramDict['p1'][0]-NSIGMA*paramDict['p1'][1],paramDict['p1'][0]+NSIGMA*paramDict['p1'][1])
                rRangeString += ':p2=%f,%f'%(paramDict['p2'][0]-NSIGMA*paramDict['p2'][1],paramDict['p2'][0]+NSIGMA*paramDict['p2'][1])
                rRangeString += ':p3=%f,%f'%(paramDict['p3'][0]-NSIGMA*paramDict['p3'][1],paramDict['p3'][0]+NSIGMA*paramDict['p3'][1])            
            if options.rMax>-1:
                rRangeString += ':r=0,%f'%(options.rMax)
            toyString = ''
            if options.toys>-1:
                toyString = '-t %i'%options.toys
            exec_me('combine -M MarkovChainMC %s/dijet_combine_%s_%s_lumi-%.3f_%s.txt -n %s_%s_lumi-%.3f_%s --tries 20 --proposal ortho --burnInSteps 100 --iteration 20000 %s %s %s %s'%(options.outDir,model,massPoint,lumi,box,model,massPoint,lumi,box,rRangeString,blindString,sysString,toyString),options.dryRun)
            exec_me('mv higgsCombine%s_%s_lumi-%.3f_%s.MarkovChainMC.mH120*root %s/'%(model,massPoint,lumi,box,options.outDir),options.dryRun)  
        else:
            if signif:
                exec_me('combine -M ProfileLikelihood --signif --expectSignal=1 -t -1 %s/dijet_combine_%s_%s_lumi-%.3f_%s.txt -n %s_%s_lumi-%.3f_%s'%(options.outDir,model,massPoint,lumi,box,model,massPoint,lumi,box),options.dryRun)
                exec_me('mv higgsCombine%s_%s_lumi-%.3f_%s.ProfileLikelihood.mH120.root %s/'%(model,massPoint,lumi,box,options.outDir),options.dryRun)
            else:
                rRangeString = ''
                if options.rMax>-1:                
                    rRangeString =  '--setPhysicsModelParameterRanges r=0,%f'%(options.rMax)
                exec_me('combine -M Asymptotic %s/dijet_combine_%s_%s_lumi-%.3f_%s.txt -n %s_%s_lumi-%.3f_%s --minimizerTolerance %f --minimizerStrategy %i %s --saveWorkspace %s %s'%(options.outDir,model,massPoint,lumi,box,model,massPoint,lumi,box,options.min_tol,options.min_strat,rRangeString,blindString,sysString),options.dryRun)
                exec_me('mv higgsCombine%s_%s_lumi-%.3f_%s.Asymptotic.mH120.root %s/'%(model,massPoint,lumi,box,options.outDir),options.dryRun)  

    
if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-b','--box',dest="box", default="CaloDijet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="gg",type="string",
                  help="signal model name")
    parser.add_option('--mass',dest="mass", default='750',type="string",
                  help="mass of resonance")
    parser.add_option('-l','--lumi',dest="lumi", default="1.918",type="string",
                  help="lumi array in fb^-1, e.g.: 1.918")
    parser.add_option('--signif',dest="signif",default=False,action='store_true',
                  help="calculate significance instead of limit")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store cards")
    parser.add_option('--min-tol',dest="min_tol",default=0.001,type="float",
                  help="minimizer tolerance (default = 0.001)")
    parser.add_option('--min-strat',dest="min_strat",default=2,type="int",
                  help="minimizer strategy (default = 2)")
    parser.add_option('--dry-run',dest="dryRun",default=False,action='store_true',
                  help="Just print out commands to run")
    parser.add_option('--penalty',dest="penalty",default=False,action='store_true',
                  help="penalty terms on background parameters")
    parser.add_option('-i','--input-fit-file',dest="inputFitFile", default='FitResults/BinnedFitResults.root',type="string",
                  help="input fit file")
    parser.add_option('--no-signal-sys',dest="noSignalSys",default=False,action='store_true',
                  help="do not create signal shape systematic histograms / uncertainties")
    parser.add_option('--no-sys',dest="noSys",default=False,action='store_true',
                  help="no systematic uncertainties when running combine")
    parser.add_option('--blind',dest="blind",default=False,action='store_true',
                  help="run only blinded expected limits")
    parser.add_option('--rMax',dest="rMax",default=-1,type="float",
                  help="maximum r value (for better precision)")
    parser.add_option('--xsec',dest="xsec",default=1,type="float",
                  help="xsec for signal in pb (r = 1)")
    parser.add_option('--jobs',dest="jobs",default=False,action='store_true',
                  help="submit jobs")
    parser.add_option('--bayes',dest="bayes",default=False,action='store_true',
                  help="bayesian limits")
    parser.add_option('--deco',dest="deco",default=False,action='store_true',
                  help="decorrelate shape parameters")
    parser.add_option('--tag',dest="tag", default='master',type="string",
                  help="tag for repository")
    parser.add_option('-q','--queue',dest="queue",default="1nh",type="string",
                  help="queue: 1nh, 8nh, 1nd, etc.")
    parser.add_option('-t','--toys',dest="toys",default=-1,type="int",
                  help="number of toys (for bayesian expected limits)")


    (options,args) = parser.parse_args()

    if options.jobs:
        submit_jobs(options,args)
    else:            
        main(options,args)
