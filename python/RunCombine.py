from optparse import OptionParser
import os
import ROOT as rt
from array import *
from framework import Config
import sys
import glob

def exec_me(command,dryRun=True):
    print command
    if not dryRun: os.system(command)

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
    parser.add_option('-l','--lumi',dest="lumi", default="1.981",type="string",
                  help="lumi array in fb^-1, e.g.: 1.981")
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
                  help="no signal shape systematic uncertainties")
    parser.add_option('--rMax',dest="rMax",default=-1,type="float",
                  help="maximum r value (for better precision)")
    parser.add_option('--xsec',dest="xsec",default=1,type="float",
                  help="xsec for signal in pb (r = 1)")

    (options,args) = parser.parse_args()

    box = options.box

    signif = options.signif
    
    lumi = float(options.lumi)
    cfg = Config.Config(options.config)
    
    model = options.model

    signalSys = ''
    if options.noSignalSys:
        signalSys = '--no-signal-sys'
    
    penaltyString = ''
    if options.penalty: penaltyString = '--penalty'

    xsecString = '--xsec %f'%(options.xsec)    

    signalDsName = 'inputs/ResonanceShapes_%s_13TeV_Scouting_Spring15.root'%model
    backgroundDsName = 'output/output_JEC.root'

    if len(options.mass.split(','))==1:
        massIterable = [options.mass]
    else:
        massIterable = list(eval(options.mass))
    for massPoint in massIterable:
        
        exec_me('python python/WriteDataCard.py -m %s --mass %s -i %s -l %f -c %s -b %s -d %s %s %s %s %s %s'%(model, massPoint, options.inputFitFile,1000*lumi,options.config,box,options.outDir,signalDsName,backgroundDsName,penaltyString,signalSys,xsecString),options.dryRun)

        
        if signif:
            exec_me('combine -M ProfileLikelihood --signif --expectSignal=1 -t -1 %s/dijet_combine_%s_%s_lumi-%.3f_%s.txt -n %s_%s_lumi-%.3f_%s'%(options.outDir,model,massPoint,lumi,box,model,massPoint,lumi,box),options.dryRun)
            exec_me('mv higgsCombine%s_%s_lumi-%.3f_%s.ProfileLikelihood.mH120.root %s/'%(model,massPoint,lumi,box,options.outDir),options.dryRun)
        else:
            if options.rMax>-1:
                exec_me('combine -M Asymptotic %s/dijet_combine_%s_%s_lumi-%.3f_%s.txt -n %s_%s_lumi-%.3f_%s --minimizerTolerance %f --minimizerStrategy %i --setPhysicsModelParameterRanges r=0,%f --saveWorkspace --noFitAsimov --run expected'%(options.outDir,model,massPoint,lumi,box,model,massPoint,lumi,box,options.min_tol,options.min_strat,options.rMax),options.dryRun)
            else:
                exec_me('combine -M Asymptotic %s/dijet_combine_%s_%s_lumi-%.3f_%s.txt -n %s_%s_lumi-%.3f_%s --minimizerTolerance %f --minimizerStrategy %i --saveWorkspace --noFitAsimov --run expected'%(options.outDir,model,massPoint,lumi,box,model,massPoint,lumi,box,options.min_tol,options.min_strat),options.dryRun)
            exec_me('mv higgsCombine%s_%s_lumi-%.3f_%s.Asymptotic.mH120.root %s/'%(model,massPoint,lumi,box,options.outDir),options.dryRun)  
