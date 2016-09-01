#! /usr/bin/env python

import os
import sys
import optparse
import datetime
import time

usage = "usage: python scripts/splitAndSubmit.py -e python/bTag_fordqm.py -t config/bTag.cfg  -b BTag2016_fordqm -q cmscaf1nh -i lists/myList_ParkingScoutingMonitor_Run2016D.txt  -o histo_mjj --split 1"

parser = optparse.OptionParser(usage)
parser.add_option("-e", "--executable", dest="executable",
    help="name of the executable",
    )

parser.add_option("-t", "--template", dest="template",
    help="name of the template config to be used",
    )

parser.add_option("-b", "--box", dest="box",
    help="analysis box to be used",
    )

parser.add_option('-q', '--queue',       action='store',     dest='queue',       
    help='run in batch in queue specified as option (default -q cmslong)', 
    default='cmsan',
    metavar="QUEUE")

parser.add_option("-i", "--input", dest="input",
    help="path and name of the fileList",
    )

parser.add_option("-o", "--output", dest="output",
    help="the root file outName",
    metavar="OUTDIR")

parser.add_option("--split", dest="filesperjob", type=int,
    help="files to analyze per job ",
    default=10)

parser.add_option('-I', '--interactive',      
    action='store_true',
    dest='interactive',      
    help='run the jobs interactively, 2 jobs at a time',
    default=False)

(opt, args) = parser.parse_args()
################################################


###
pwd = os.environ['PWD']
current_time = datetime.datetime.now()
simpletimeMarker = "_%04d%02d%02d_%02d%02d%02d" % (current_time.year,current_time.month,current_time.day,current_time.hour,current_time.minute,current_time.second) 
timeMarker = "submit_%04d%02d%02d_%02d%02d%02d" % (current_time.year,current_time.month,current_time.day,current_time.hour,current_time.minute,current_time.second) 
workingDir = pwd+"/mybatch/"+timeMarker
executable = pwd+"/"+opt.executable

os.system("mkdir -p "+workingDir)
os.system("cp "+opt.template+" "+workingDir)

head, tail = os.path.split(opt.template)
template = workingDir+"/"+tail

inputlist = []
njobs_list = []

num_lines = sum(1 for line in open(opt.input, "r"))

ins = open(opt.input, "r")
##loop over lists (one for datasets) to create splitted lists
count = 0
jobCount = 0
for line in  ins:
    count = count+1
    line = "root://eoscms.cern.ch//eos/cms"+line.rstrip('\n')
    inputlist.append(line)
    if count%opt.filesperjob == 0 or count==num_lines:
        jobCount = jobCount+1
        os.system("mkdir "+workingDir+"/"+str(jobCount))
        
        with open(template) as fi:
            contents = fi.read()
            replaced_contents = contents.replace('INPUTLIST', str(inputlist))
            replaced_contents = replaced_contents.replace('OUTPUTFILE', "\""+opt.output+"_"+str(jobCount)+".root\"")
        with open(workingDir+"/"+str(jobCount)+"/config.cfg", "w") as fo:
            fo.write(replaced_contents)

        inputlist = []
        
        os.system("echo cd "+pwd+" > launch.sh")
        os.system("echo 'eval `scramv1 runtime -sh`\n' >> launch.sh")
        os.system("echo cd - >> launch.sh")
        os.system("echo python "+executable+" -c "+workingDir+"/"+str(jobCount)+"/config.cfg -b "+opt.box+" >> launch.sh")
        os.system("echo mv "+opt.output+"_"+str(jobCount)+".root "+workingDir+" >> launch.sh")
        os.system("chmod 755 launch.sh")
        os.system("mv launch.sh "+workingDir+"/"+str(jobCount))
        njobs_list.append("bsub -q "+opt.queue+" -o "+workingDir+"/"+str(jobCount)+"/log.out -e "+workingDir+"/"+str(jobCount)+"/log.err "+workingDir+"/"+str(jobCount)+"/launch.sh")
        
for job in njobs_list:
    print job
