#! /usr/bin/env python

import os
import sys
import optparse
import datetime
import time

usage = "usage: To be run from DijetRootTreeAnalyzer/ :  python scripts/submit_batch_EOS_split.py -q 8nh -i test/santanas/list/list_ScoutingPF__17_01_2016/ -o /eos/cms/store/group/phys_exotica/dijet/Dijet13TeVScouting/rootTrees_reduced/TEST/ --split 2 -m ScoutingPFHT --tag ScoutingPF__17_01_2016 -c config/cutFile_mainDijetScoutingSelection.txt"

#parser = optparse.OptionParser("submitAllGJetsID.py")
parser = optparse.OptionParser(usage)
parser.add_option('-q', '--queue',       action='store',     dest='queue',       
    help='run in batch in queue specified as option (default -q cmslong)', 
    default='cmsan',
    metavar="QUEUE")

parser.add_option("-i", "--input", dest="input",
    help="directory containing the lists of samples to be analyzed",
    )

parser.add_option("-o", "--output", dest="output",
    help="the directory OUTDIR contains the output of the program",
    metavar="OUTDIR")

parser.add_option("--timestamp", dest="timestamp",
    help="rerun failed jobs for specific timestamp",
    default="")

parser.add_option("-m", "--match", dest="match",
    help="run only the samples containing this string in the name",
    default="")

parser.add_option("--tag", dest="tag",
    help="useful tag to identify the version",
    default="")

parser.add_option("--split", dest="filesperjob", type=int,
    help="files to analyze per job ",
    default=10)

parser.add_option('-I', '--interactive',      
    action='store_true',
    dest='interactive',      
    help='run the jobs interactively, 2 jobs at a time',
    default=False)

parser.add_option("-c", "--cut", dest="cutfile",
    help="cutfile",                  
    )

(opt, args) = parser.parse_args()
################################################

###
current_time = datetime.datetime.now()
simpletimeMarker = "_%04d%02d%02d_%02d%02d%02d" % (current_time.year,current_time.month,current_time.day,current_time.hour,current_time.minute,current_time.second) 
timeMarker = "mycutFile_%04d%02d%02d_%02d%02d%02d__" % (current_time.year,current_time.month,current_time.day,current_time.hour,current_time.minute,current_time.second) 
cutfileName = timeMarker+os.path.split(opt.cutfile)[1]
print cutfileName

if not opt.timestamp:
  newTag = opt.tag+simpletimeMarker
else:
  newTag = opt.tag+"_"+opt.timestamp

###
os.system("/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select mkdir -p "+opt.output+"/"+newTag)
#os.system("rm -rf batch")
os.system("mkdir -p batch")
pwd = os.environ['PWD']
###
os.system("cp "+opt.cutfile+" batch/"+cutfileName)
###

if not opt.match:
  print opt.input
  os.system("ls "+opt.input+" | grep txt | grep -v \"\~\" > config/lists_to_run.txt")#NEW
else:
  os.system("ls "+opt.input+" | grep txt | grep "+opt.match+" | grep -v \"\~\"  > config/lists_to_run.txt")#NEW

ins = open("config/lists_to_run.txt", "r") 


#lists of lists (in each position there is a list of commads)
inputlists = []
njobs_list = []
#commands = []
#hadd_cmd = []
#filenames_skim = []
splittedDir = "splitted"+"_"+newTag #simpletimeMarker

##create a directory to store logfiles containing the "tag" in the name
os.system("mkdir -p batch/"+newTag)
submitCommandsFile = open("batch/"+newTag+"/bsub_commands.txt","a+")
#os.system("touch batch/"+newTag+"/submitCommands.txt")
os.system("mkdir -p "+opt.input+"/"+splittedDir)

##loop over lists (one for datasets) to create splitted lists
for line in  ins:
  lists_dataset = []
  sample = os.path.splitext(line)[0]
  print ("process %s" % sample)
  line = line.rstrip('\n')
  sample = sample.rstrip('\n')
  # name of final rootfile merged
  #filenames_skim.append(opt.output+simpletimeMarker+"/rootfile_"+sample+"_"+newTag+"_"+str(jj)+"_reduced_skim.root")
  # list of commands for current dataset
  #commands_dataset = []
  
  jf=0  ##counter for num of files
  jj=0  ##counter for num of jobs
  #open list

  total_num = 0
  for file in open(opt.input+"/"+line,"r"):
    total_num +=1
  print total_num

  list = open(opt.input+"/"+line,"r") 
  print ""

  for file in list:
    print "file:%i  filesperjob:%i  job:%i op.modulo:%i  list %s " % (jf, opt.filesperjob,jj,(int(jf+1) % int(opt.filesperjob)), opt.input+"/"+line)
    print file
    modulo = int(jf+1) % int(opt.filesperjob)
    #print "modulo = %i" % modulo
    splittedlist = open(opt.input+"/"+splittedDir+"/"+sample+"_"+newTag+"_"+str(jj)+".txt","a+")
    splittedlist.write(file)
    if ( modulo == 0 or jf+1==total_num):
      lists_dataset.append(opt.input+"/"+splittedDir+"/"+sample+"_"+newTag+"_"+str(jj)+".txt")
      print "==> job "+str(jj)+"   appending "+opt.input+"/"+splittedDir+"/"+sample+"_"+newTag+"_"+str(jj)+".txt"
      print ""
      jj += 1 #increment counter of jobs  
    jf += 1   #increment counter of files         
  #print "job "+str(jj)+"   appending "+opt.input+"/"+splittedDir+"/"+sample+"_"+newTag+"_"+str(jj)+".txt"
  #lists_dataset.append(opt.input+"/"+splittedDir+"/"+sample+"_"+newTag+"_"+str(jj)+".txt")  
  njobs_list.append(jj-1)
  inputlists.append(lists_dataset)

#print njobs_list
print inputlists

i_f = 0
ins = open("config/lists_to_run.txt", "r") 
for line in  ins:
  lists_dataset = []
  sample = os.path.splitext(line)[0]
  line = line.rstrip('\n')
  sample = sample.rstrip('\n')
  splittedlist = inputlists[i_f]
  print splittedlist
  for jj in range(0,njobs_list[i_f]+1):
    print ""
    print sample+"  job "+str(jj)

    logfile = "batch/"+newTag+"/logfile_"+sample+"_"+newTag+"_"+str(jj)+".log"#NEW
    crablogfile = "batch/"+newTag+"/crablogfile_"+sample+"_"+newTag+"_"+str(jj)+".crablog"#NEW

    ntupleDone = opt.output[1:]+"/"+newTag+"/rootfile_"+sample+"_"+newTag+"_"+str(jj)+"_reduced_skim.root"

    if os.path.isfile(ntupleDone)==True:
      print "This job has already been run, output is here: "+ntupleDone+"\n"
      continue

    ###################################
    #command = "./main "+splittedlist[jj]+" config/cutFile_mainDijetSelection.txt dijets/events "+opt.output+simpletimeMarker+"/rootfile_"+sample+"_"+newTag+"_"+str(jj)+" "+opt.output+simpletimeMarker+"/cutEfficiencyFile_"+sample+"_"+newTag+"_"+str(jj)
    #command = "./main "+splittedlist[jj]+" config/cutFile_mainDijetSelection.txt dijets/events $TWD/rootfile_"+sample+"_"+newTag+"_"+str(jj)+" $TWD/cutEfficiencyFile_"+sample+"_"+newTag+"_"+str(jj)
    #command = "./main "+splittedlist[jj]+" batch/"+cutfileName+" dijets/events $TWD/rootfile_"+sample+"_"+newTag+"_"+str(jj)+" $TWD/cutEfficiencyFile_"+sample+"_"+newTag+"_"+str(jj)
    #command = "./main "+splittedlist[jj]+" batch/"+cutfileName+" dijetscouting/events $TWD/rootfile_"+sample+"_"+newTag+"_"+str(jj)+" $TWD/cutEfficiencyFile_"+sample+"_"+newTag+"_"+str(jj)+" >& "+"$TWD/logfile_"+sample+"_"+newTag+"_"+str(jj)+".log"#NEW
    command = "./main "+splittedlist[jj]+" batch/"+cutfileName+" dijetscouting/events $TWD/rootfile_"+sample+"_"+newTag+"_"+str(jj)+" $TWD/cutEfficiencyFile_"+sample+"_"+newTag+"_"+str(jj)+" >& "+pwd+"/"+logfile #NEW
    ###################################

    #print "submit "+command
    outputname = "batch/"+newTag+"/submit_"+sample+"_"+newTag+"_"+str(jj)+".src"
    outputfile = open(outputname,'w')
    #outputfile.write('#!/bin/bash\n')
    outputfile.write('#!/usr/bin/env bash -x\n')
    outputfile.write('export SCRAM_ARCH=slc6_amd64_gcc491\n')
    outputfile.write('pwd\n')
    outputfile.write('export TWD=$PWD\n')
    outputfile.write('cd '+pwd+' \n')
    outputfile.write('pwd\n')
    outputfile.write('eval `scramv1 runtime -sh`\n')
    outputfile.write(command+"\n")
    outputfile.write("ls -lrth $TWD/\n")
    #outputfile.write("ls $TWD/ | grep "+newTag+"\n")
    outputfile.write("xrdcp -f $TWD/rootfile_"+sample+"_"+newTag+"_"+str(jj)+"_reduced_skim.root root://eoscms/"+opt.output+"/"+newTag+"/rootfile_"+sample+"_"+newTag+"_"+str(jj)+"_reduced_skim.root\n")##NEW
    outputfile.write("xrdcp -f $TWD/rootfile_"+sample+"_"+newTag+"_"+str(jj)+".root root://eoscms/"+opt.output+"/"+newTag+"/rootfile_"+sample+"_"+newTag+"_"+str(jj)+".root\n")##NEW
    outputfile.write("xrdcp -f $TWD/cutEfficiencyFile_"+sample+"_"+newTag+"_"+str(jj)+".dat root://eoscms/"+opt.output+"/"+newTag+"/cutEfficiencyFile_"+sample+"_"+newTag+"_"+str(jj)+".dat\n")##NEW
    #outputfile.write("xrdcp -f $TWD/logfile_"+sample+"_"+newTag+"_"+str(jj)+".log "+pwd+"/"+logfile+"\n")##NEW
    outputfile.write("rm $TWD/rootfile_"+sample+"_"+newTag+"_"+str(jj)+"_reduced_skim.root\n")
    outputfile.write("rm $TWD/rootfile_"+sample+"_"+newTag+"_"+str(jj)+".root\n")
    outputfile.write("rm $TWD/cutEfficiencyFile_"+sample+"_"+newTag+"_"+str(jj)+".dat\n")
    #outputfile.write("rm $TWD/logfile_"+sample+"_"+newTag+"_"+str(jj)+".log\n")
    outputfile.write("ls -lrth $TWD/\n")
    outputfile.close()  
    os.system("chmod 755 "+outputname)

    #print outputname 
    if opt.interactive==False:
      bsubCommand = "bsub -q "+opt.queue+" -o "+pwd+"/"+crablogfile+" source "+pwd+"/"+outputname
      #bsubCommand = "bsub -q "+opt.queue+" -o "+pwd+"/"+crablogfile+" < "+pwd+"/"+outputname
      print bsubCommand ##NEW
      submitCommandsFile.write(bsubCommand+"\n")      
      #time.sleep(3)
      os.system(bsubCommand)##NEW
    else:
      print logfile
      if imc==0: os.system(command+" >&! "+logfile+"&")
      else: os.system(command+" >&! "+logfile)
  i_f += 1
     

submitCommandsFile.close() 
