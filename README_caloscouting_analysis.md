
Instructions to run Calo Scouting Dijet Resonance Search from start to finish.

### Set up CMSSW/DijetRootTreeAnalyzer/combine
1. Set up CMSSW/DijetRootTreeAnalyzer/combine
    ```sh
    $ cmsrel CMSSW_7_4_14
    $ cd CMSSW_7_4_14/src
    $ cmsenv
    $ git clone https://github.com/RazorCMS/DijetRootTreeAnalyzer CMSDIJET/DijetRootTreeAnalyzer
    $ git clone -b dijetpdf_74X https://github.com/RazorCMS/HiggsAnalysis-CombinedLimit HiggsAnalysis/CombinedLimit
    $ cd HiggsAnalysis/CombinedLimit
    $ scram b -j 4
    $ cd $CMSSW_BASE/CMSDIJET/DijetRootTreeAnalyzer
    ```
### Getting mjj histogram for fit
2. Add files from EOS together
    ```sh
    $ python python/AddFiles.py -l lists/CaloScoutingHT_JEC_CaloL1L2L3_PFL2L3Residual_20160503_171912_reduced_skim.txt -o inputs/data_CaloScoutingHT_Run2015D_CaloDijet.root
    ```
### Fitting
3. Copy signal resonance files
    ```sh
    $ wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/68e8514f4da8849b99b7dfcf1a7834fa55aeefa6/ResonanceShapes_gg_13TeV_Scouting_Spring15.root -P inputs/
    $ wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/68e8514f4da8849b99b7dfcf1a7834fa55aeefa6/ResonanceShapes_qg_13TeV_Scouting_Spring15.root -P inputs/
    $ wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/68e8514f4da8849b99b7dfcf1a7834fa55aeefa6/ResonanceShapes_qq_13TeV_Scouting_Spring15.root -P inputs/
    ```
1. Perform a binned, blinded, background-only, sideband fit (with signal shape plotted)
    ```sh
    $ mkdir fits_2016_05_07/
    $ mkdir fits_2016_05_07/CaloDijet_Sideband/
    $ python python/BinnedFit.py -c config/dijet.config -b CaloDijet inputs/data_CaloScoutingHT_Run2015D_CaloDijet.root --lumi 1918 --fit-region Low,High --plot-region Low,High -d fits_2016_05_07/CaloDijet_Sideband/ -s inputs/ResonanceShapes_gg_13TeV_Scouting_Spring15.root -m gg --mass 750 --xsec 15
    ```
1. Run and fit 1000 toys and plot GOF
    ```sh
    $ python python/RunToys.py -b CaloDijet --freq -c config/dijet.config --lumi 1918 --fit-region Low,High  -d fits_2016_05_07/CaloDijet_Sideband/ -i fits_2016_05_07/CaloDijet_Sideband/DijetFitResults_CaloDijet.root  -t 1000 -s 0
    $ python python/PlotGOF.py -b CaloDijet -c config/dijet.config -d fits_2016_05_07/CaloDijet_Sideband/ -t fits_2016_05_07/CaloDijet_Sideband/toys_Freq_s0_CaloDijet.root -l 1918 --data
    ```
### Setting limits
5. Run combine to produce blinded, expected limits for gg resonance mass range [50, 1600] in 50 GeV steps
    ```sh
    $ mkdir cards
    $ python python/RunCombine.py -m gg -d cards/ --mass range\(500,1650,50\) -c config/dijet.config -i fits_2016_05_07/CaloDijet_Sideband/DijetFitResults_CaloDijet.root -b CaloDijet --rMax 20 --xsec 10 -l 1.918 --blind
    ```
1. Get combine output files and plot 1D limit
    ```sh
    $ python python/GetCombine.py -d cards/ -m gg --mass range\(500,1650,50\) -b CaloDijet --xsec 10 -l 1.918
    $ python python/Plot1DLimit.py -d cards/ -m gg -b CaloDijet -l 1.918
    ```




