#!/usr/bin/python
import argparse, os, tempfile, shutil, sys,math,pickle,itertools
#parent = os.path.dirname(os.getcwd())
#sys.path.append(parent)
from subprocess import call, PIPE, STDOUT, Popen
import argparse
import datetime


parser = argparse.ArgumentParser(description='Submit jot to batch')
parser.add_argument('--config_file',type=str,dest="config_file",default=1, help=" config_file of the era")
parser.add_argument('--outputname',type=str,dest="outputname",default=1, help="Run era")
parser.add_argument('--list',type=str,dest="list",default=1, help="list of file to run")
parser.add_argument('--dirlist',type=str,dest="dirlist",default=1, help="directory where are lists of file to run")
args = parser.parse_args()
config_file = args.config_file #path to the processed lumi JSON file
outputname = args.outputname # which run
listtorun = args.list
dirlist = args.dirlist
outputdir = "/afs/cern.ch/work/h/hlattaud/private/results_final/"
cmd4="./main "+dirlist+listtorun+" "+config_file+" dijets/events "+outputdir+outputname+" "+outputdir+outputname
cmd1="cd /afs/cern.ch/work/h/hlattaud/private/CMSSW_8_0_25/src/CMSDIJET/responsecomputing/DijetRootTreeAnalyzer/"
#cmd2="source /afs/cern.ch/work/h/hlattaud/private/CMSSW_8_0_25/src/CMSDIJET/responsecomputing/DijetRootTreeAnalyzer/set_env.sh"
cmd2="export SCRAM_ARCH=slc6_amd64_gcc530"
cmd3="eval `scramv1 runtime -sh`"
call(cmd1+" && "+cmd2+" && "+cmd3+" && "+cmd4, shell=True )
#call(cmd2, shell=True)
#os.system(cmd1)
#os.system(cmd2)
#os.system(cmd3)






print(cmd4)
#os.system(cmd4)
#call(cmd4, shell=True, env=True)
