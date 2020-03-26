# Instructions for producing Gamma + jet step 2 Analyse Nuples

Step 1 is [here](https://github.com/lucastorterotot/DijetRootTreeMaker/blob/master/instructions/GammaJetTree_Instruction.md).

**Target** CMS Gamma + Jet  - JEC calculation

**Authors** Hugues Lattaud and Lucas Torterotot

**Last update** 26 Mar 2020

## Introduction

This package provides a small facility to analyze one or a chain of root ntuples.

A script (`./scripts/make_rootNtupleClass.sh`) is used to generate automatically
(using the root command `RootNtupleMaker->MakeClass`) a class (`include/rootNtupleClass.h`
and `src/rootNtupleClass.C`) with the variable definitions of a given root ntuple
(to be provided by the user).

The class baseClass (`include/baseClass.h` and `src/baseClass.C`) inherits from the
automatically generated `rootNtupleClass`.
`baseClass` provides the methods that are common to all analysis, such as the method
to read a list of root files and form a chain. It will, asap, also provide a method
to read a list of selection cuts.

The class `analysisClass` (`include/analysisClass.h` and `src/analysisClass.C`) inherits
from `baseClass`.
The user's code should be placed in the method `Loop()` of `analysisClass`, which reimplements
the method `Loop()` of `rootNtupleClass`.

The main program (`src/main.C`) receives the configuration parameters (such as the input
chain of root files and a file to provide a cut list) and executes the `analysisClass` code.

## Installation recipe
### CMSSW release
```
cmsrel CMSSW_8_0_31
cd CMSSW_8_0_31/src/
cmsenv
```
### Get the Tree Analyzer code
Get this git repository and choose the branch corresponding to your era of interest by doing
```
git clone -o lucas git@github.com:lucastorterotot/DijetRootTreeAnalyzer.git -b <BRANCH> CMSDIJET/DijetRootTreeAnalyzer
```
where `<BRANCH>` is corresponding to your era of interest:
* 2016: `JEC_JER_2016_master`
* 2017: `JEC_JER_2017_master`
* 2017UL: `JEC_JER_2017UL_CMSSW_8_0_31_master`
* 2018: `JEC_JER_2018_CMSSW_8_0_31_master`
* 2018UL : `JEC_JER_2018UL_CMSSW_8_0_31_master`

### Hot zones maps
Some era use a so-called `EtaPhiCleaning_File` you have to retreive:
* 2017: 
```
cp /afs/cern.ch/work/l/ltortero/public/EtaPhiCleaning_Files/* $CMSSW_BASE/src/CMSDIJET/DijetRootTreeAnalyzer/.
```
* 2017UL: use the one provided in the [JECDatabse repository](https://github.com/cms-jet/JECDatabase) you have cloned while installing the first step of this framework. The following commands, if you followed instructions, should work:
```
cd $CMSSW_BASE/../JECDatabase ; git pull ; cd - ; cp $CMSSW_BASE/../JECDatabase/hotzone_maps/Summer19UL17_V1/hotjets-UL17.root $CMSSW_BASE/src/CMSDIJET/DijetRootTreeAnalyzer/.
```

## Instructions (local)
* Go inside the right directory
```
cd CMSDIJET/DijetRootTreeAnalyzer
```

* Generate the rootNtupleClass with `./scripts/make_rootNtupleClass.sh`. You will be asked for input arguments, for example
```
./scripts/make_rootNtupleClass.sh -f GammajetTree.root -t dijets/events
```
where the rootfile has been produced at the [DijetRootTreeMaker step](https://github.com/lucastorterotot/DijetRootTreeMaker).
Note that there could be minor glitches and then you can fix by hand.
Type yes when asked to overwrite rootNtupleClass.h

* Make a symbolic link to `analysisClass.C` with
```
ln -sf analysisClass_GammaplusJet_New_methods.C src/analysisClass.C 
```
or 
```
ln -sf analysisClass_GammaplusJet_fakejethandlePuppi.C src/analysisClass.C
```
to run on puppi jet.
   
* Compile to test that all is OK so far:
```
make clean && make
```

* Run `./main`. You will be asked for input arguments, for example
```
./main inputList.txt config/cutFile_mainGammaplusJetSelection_Newsmearing.txt dijets/events output/rootFileName output/rootFileName
```
 All the produced files are in the directory `output`.

If everything went fine your root file should contain two tree, one tree with all the information needed for the PUreweighting and the other with all the variables after selection. You are now ready to start the full production.
   
## Instructions (full production)   
To run your analysis on the LXplus batch use the Submit_job_to_batch_short.py script. To use it first make a list of the output of [DijettreeMaker step](https://github.com/lucastorterotot/DijetRootTreeMaker/blob/master/instructions/GammaJetTree_Instruction.txt).

Then move the list in a new directory and split it 10 lines by 10 lines:
```
split -l 10 -d NameOfYourList.txt
```
In the file `Run_Analyser.py`, change the path to the output directory (`outputdir ="path to your output directory"`). For example:
```
mkdir MC_summer
cd MC_Summer

eos ls root://eoscms//eos/cms/store/group/phys_jetmet/hlattaud/GammaJet/GJet-RunIISummer16MiniAODv2/GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_20M/
crab_NEWGJet-RunIISummer16MiniAODv2_Newsmearing > list_MC.txt 
split -l 10 list_MC.txt
   
cd ../
python Submit_job_to_batch_short.py --config_file config/cutFile_mainGammaplusJetSelection_Newsmearing.txt --whichrun MC_summer_datalike --Nlist 4 --dirlist MC_summer/
```
You will be ask for your grid password.
   
You can now go to the next step in [NewGammaJet](https://github.com/lucastorterotot/NewGammaJet) to perform trigger selection and PUreweighting.
   
**NB**
* In the cut file you can chose etheir to apply the JEC on the fly or not.
If you choose to re-apply  the JEC on the fly you have to provide those you want to use in the cut file.
* One can have several analyses in a directory, such as
⋅⋅ src/analysisClass_myCode1.C
⋅⋅ src/analysisClass_myCode2.C
⋅⋅ src/analysisClass_myCode3.C
and move the symbolic link to the one to be used:
```
ln -sf analysisClass_myCode2.C src/analysisClass.C
```
and compile/run as above.
   
## More details

### Example code
The `src/analysisClass_template.C` comes with simple example code. The example code in enclosed by
```
#ifdef USE_EXAMPLE
   ... code ...
#endif //end of USE_EXAMPLE
```
The code is NOT compiled by default. In order to compile it, uncomment the line
```
#FLAGS += -DUSE_EXAMPLE
```
in the Makefile.

### Providing cuts via file
A list of cut variable names and cut limits can be provided through a file (see `config/cutFileExample.txt`).
The variable names in such a file have to be filled with a value calculated by the user `analysisClass` code,
a function `fillVariableWithValue` is provided - see example code.
Once all the cut variables have been filled, the cuts can be evaluated by calling `evaluateCuts` - see
example code. Do not forget to reset the cuts by calling `resetCuts` at each event before filling the
variables - see example code.
The function `evaluateCuts` determines whether the cuts are satisfied or not, stores the pass/failed result
of each cut, calculates cut efficiencies and fills histograms for each cut variable (binning provided by the
cut file, see `config/cutFileExample.txt`).
The user has access to the cut results via a set of functions (see `include/baseClass.h`)
```
  bool baseClass::passedCut(const string& s);
  bool baseClass::passedAllPreviousCuts(const string& s);
  bool baseClass::passedAllOtherCuts(const string& s);
```
where the string to be passed is the cut variable name.
The cuts are evaluated following the order of their apperance in the cut file (`config/cutFileExample.txt`).
One can simply change the sequnce of line in the cut file to have the cuts applied in a different order
and do cut efficiency studies.
Also, the user can assign to each cut a level (0,1,2,3,4 ... n) and use a function
```
  bool baseClass::passedAllOtherSameLevelCuts(const string& s);
```
to have the pass/failed info on all other cuts with the same level.
There is actually also cuts with `level=-1`. These cuts are not actually evaluated, the corresponding lines
in the cut file (`config/cutFileExample.txt`) are used to pass values to the user code (such as fiducial
region limits). The user can access these values (and also those of the cuts with `level >= 0`) by
```
  double baseClass::getCutMinValue1(const string& s);
  double baseClass::getCutMaxValue1(const string& s);
  double baseClass::getCutMinValue2(const string& s);
  double baseClass::getCutMaxValue2(const string& s);
```

### Automatic histograms for cuts
The following histograms are generated for each cut variable with `level >= 0`:
* no cuts applied
* passedAllPreviousCuts
* passedAllOtherSameLevelCuts
* passedAllOtherCuts
* passedAllCut
and by default only the following subset
* no cuts applied
* passedAllPreviousCuts
* passedAllOtherCuts
is saved to the output root file. All histograms can be saved to the output root file by
uncommenting the following line in the Makefile
```
#FLAGS += -DSAVE_ALL_HISTOGRAMS
```

### Automatic cut efficiency
The absolute and relative efficiency is calculated for each cut and stored in an output file
(named `output/cutEfficiencyFile.dat` if the code is executed following the examples).

The user has the option to implement a good run list using a JSON file.  This requires two edits to the cut 
file and one edit to the `analysisClass.C` file.
A line must be inserted at the beginning of the cut file with the word "JSON" first, and then 
the full AFS path of the desiredJSON file. For example:
```
    JSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-163369_7TeV_PromptReco_Collisions11_JSON.txt
```
In addition, the user must define the JSON file selection in the cut file.  This is done in the usual way:
```
    #VariableName                   minValue1(<) maxValue1(>=)      minValue2(<)    maxValue2(>=)   level   histoNbinsMinMax
    #------------                   ------------ -------------      ------------    -------------   -----   ----------------
    PassJSON                        0            1                  -               -               0       2 -0.5 1.5
```
In the `analysisClass.C` file, the user must add the following line within the analysis loop:
```
    fillVariableWithValue ( "PassJSON", passJSON (run, ls, isData));
```

Note that the use of a JSON file (good run list) is optional.  If the user does not list a JSON file in the cut file,
no selection will be made.

This tool was validated on May 18, 2011.  The results of the validation study are [here](https://twiki.cern.ch/twiki/pub/CMS/ExoticaLeptoquarkAnalyzerRootNtupleMakerV2/json_validation.pdf).

### Additional scripts for running on several datasets

See [./doc/howToMakeAnalysisWithRootTuples.txt](https://github.com/lucastorterotot/DijetRootTreeAnalyzer/tree/master/doc/howToMakeAnalysisWithRootTuples.txt).


### Producing an ntuple skim (Dinko Ferencek)

The class `baseClass` provides the ability to produce a skimmed version of the input ntuples. In order to
produce a skim, the following preliminary cut line has to be added to the cut file
```
#VariableName         value1            value2          value3          value4          level
#------------         ------------      -------------   ------------    -------------   -----
produceSkim           1                 -               -               -               -1
```
and call the `fillSkimTree()` method for those events that meet the skimming criteria. One possible example is
```
    if( passedCut("all") ) fillSkimTree();
```
If the above preliminary cut line is not present in the cut file, is commented out or its value1 is set to 0,
the skim creation will be turned off and calling the `fillSkimTree()` method will have no effect.

###JSON parser (Edmund Berry)

See [this hyper news](https://hypernews.cern.ch/HyperNews/CMS/get/exotica-lq/266.html).

### Producing a new ntuple with a subset of cutFile variables and a subset of events (Paolo, Francesco, Edmund)

The class `baseClass` provides the ability to produce a new ntuple with a subset of the variables defined
in the cutFile, and with a subset of events.
In order to do so, the following preliminary cut line has to be added to the cut file
```
#VariableName         value1            value2          value3          value4          level
#------------         ------------      -------------   ------------    -------------   -----
produceReducedSkim              1               -               -               -               -1
```
then each variable that needs to be included in the new tree has to be flagged with SAVE in 
the cutFile at the end of the line where the variabole is defined, as for `pT1stEle` and `pT2ndEle`
below:
```
#VariableName	      minValue1(<) maxValue1(>=)	minValue2(<)	maxValue2(>=)	level	histoNbinsMinMax  OptionalFlag
#------------	      ------------ -------------	------------	-------------	-----	----------------  ------------
nEleFinal	      1		   +inf			-		-		0	11 -0.5 10.5
pT1stEle              85           +inf                 -               -               1       100 0 1000        SAVE
pT2ndEle	      30	   +inf			-	        -	        1	100 0 1000        SAVE
invMass_ee	      0		   80			100	        +inf	        1	120 0 1200
```
(do not put anything for those variables that do not need to be saved, such as for  `nEleFinaland invMass_ee`)

finally, call `fillReducedSkimTree()` in the `analysisClass` for the subset of events that need to be saved, e.g.:
```
    if( passedCut("nEleFinal") ) fillReducedSkimTree();
```
If the above preliminary cut line is not present in the cut file, is commented out or its value1 is set to 0,
the skim creation will be turned off and calling the `fillReducedSkimTree()` method will have no effect.
The new ntuple will be created in a file named as the std output root file with _reduced_skim appended
before the .root and the tree name will be as in the input root file.


