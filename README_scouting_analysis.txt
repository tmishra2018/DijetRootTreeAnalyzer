#######################################
1) Create big root trees
#######################################

* Setup
https://github.com/CMSDIJET/DijetScoutingRootTreeMaker/blob/master/README.md
cd $CMSSW_BASE/src/
git clone https://github.com/CMSDIJET/Utilities.git CMSDIJET/Utilities

* Submit jobs with crab 
https://github.com/CMSDIJET/DijetScoutingRootTreeMaker/blob/master/prod/submitJobsWithCrab3/createAndSubmitCrab.py
more info at: https://twiki.cern.ch/twiki/bin/view/CMS/ExoDijet13TeVSubmitJobsWithCrab2014

* Copy files to your T2 if needed
https://github.com/CMSDIJET/Utilities/blob/master/scripts/copyWithLcg.py

* The big root trees location+info is at 
https://github.com/CMSDIJET/DijetScoutingRootTreeMaker/blob/master/data/json/README.md

#######################################
2) Create reduced root trees
#######################################

* Create lists
https://github.com/CMSDIJET/Utilities/blob/master/scripts/createList_T2.py
note: you can specify more that one folder. All the samples with the dataset name will be merged in the corresponding list.

* Setup analyzer
Main analysis:
analyzer: https://github.com/CMSDIJET/DijetRootTreeAnalyzer/blob/master/src/analysisClass_mainDijetScoutingSelection.C
cutfile: https://github.com/CMSDIJET/DijetRootTreeAnalyzer/blob/master/config/cutFile_mainDijetScoutingSelection.txt
HLT vs RECO comparison:
analyzer: https://github.com/CMSDIJET/DijetRootTreeAnalyzer/blob/master/src/analysisClass_mainDijetScoutingMonitor.C
cutfile: https://github.com/CMSDIJET/DijetRootTreeAnalyzer/blob/master/config/cutFile_mainDijetScoutingMonitor.txt
Example:
./scripts/make_rootNtupleClass.sh -f test/santanas/DijetScouting/testCode/ScoutingPFCommissioning__Run2015D-v1__RAW_74.root -t dijetscouting/events
ln -sf analysisClass_mainDijetScoutingSelection.C src/analysisClass.C
make clean
make
./main test/santanas/DijetScouting/testCode/testList.txt config/cutFile_mainDijetScoutingSelection.txt dijetscouting/events test/santanas/DijetScouting/testCode/output test/santanas/DijetScouting/testCode/output

* Submit jobs in batch
At Rome T2: 
submit: https://github.com/CMSDIJET/DijetRootTreeAnalyzer/blob/master/scripts/submit_batch_T2_split.py
At CERN: 
submit: https://github.com/CMSDIJET/DijetRootTreeAnalyzer/blob/master/scripts/submit_batch_EOS_split.py
check output: https://github.com/CMSDIJET/DijetRootTreeAnalyzer/blob/master/scripts/check_batch_EOS_split.py
NOTE: In both cases, edit this line if needed:
command = "./main "+splittedlist[jj]+" batch/"+cutfileName+" dijetscouting/events /tmp/rootfile_"+sample+"_"+newTag+"_"+str(jj)+" /tmp/cutEfficiencyFile_"+sample+"_"+newTag+"_"+str(jj)
Example:
python scripts/submit_batch_T2_split.py -i $MYCMSSW/src/CMSDIJET/DijetRootTreeAnalyzer/test/santanas/DijetScouting/list/lists__ScoutingPFCommissioning__14_01_2016/ -o /t3/users/santanas/Dijet13TeVScouting/rootTrees_reduced/ScoutingPFCommissioning__14_01_2016 -q cmsan -m ScoutingPFCommissioning --split 10 --tag ScoutingPFCommissioning__14_01_2016 -c config/cutFile_mainDijetScoutingSelection.txt
note: you must specify the folder that contains all the input lists you want to submit 
note: you can use the option match to submit only a set of lists  ( ex. -m ExpressPhysics__Run2015B-Express-v1__FEVT.txt )
note: you can use the option tag to identify the output files better ( ex.  --tag santanas__40pb-1_30_07_2015 )
note: you can use the option split to specify how many files per job ( ex.  --split 5 )

* Merge output files (needed only if you want a single output reduced tree)
At Rome T2: 
https://github.com/CMSDIJET/DijetRootTreeAnalyzer/blob/master/scripts/mergeJobs.sh
At CERN:
https://github.com/CMSDIJET/DijetRootTreeAnalyzer/blob/master/scripts/mergeJobs_EOS.sh
NOTE: 
the files merged will be in '/..output../merged/' folder

* The reduced root trees location+info is at 
https://github.com/CMSDIJET/DijetScoutingRootTreeMaker/blob/master/data/json/README.md

#######################################
3) Run on reduced trees
#######################################

* Main analysis
analyzer: https://github.com/CMSDIJET/DijetRootTreeAnalyzer/blob/master/src/analysisClass_mainDijetScoutingPlots.C
cutfile: https://github.com/CMSDIJET/DijetRootTreeAnalyzer/blob/master/config/cutFile_mainDijetScoutingPlots.txt

* HLT vs RECO comparison
analyzer: https://github.com/CMSDIJET/DijetRootTreeAnalyzer/blob/master/src/analysisClass_HLTRecoComparison.C
cutfile: https://github.com/CMSDIJET/DijetRootTreeAnalyzer/blob/master/config/cutFile_HLTRecoComparison.txt

* Trigger efficiency
https://github.com/CMSDIJET/DijetRootTreeAnalyzer/blob/master/scripts/triggerEfficiency.C
 
* Fit
...
