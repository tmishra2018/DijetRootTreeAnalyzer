
Instructions to run Calo Scouting Dijet Resonance Search from start to finish.

### Set up CMSSW/DijetRootTreeAnalyzer/combine
1. Set up CMSSW/DijetRootTreeAnalyzer/combine
    
    ```sh
    $ cmsrel CMSSW_7_4_14
    $ cd CMSSW_7_4_14/src
    $ cmsenv
    $ git clone https://github.com/CMSDIJET/DijetRootTreeAnalyzer CMSDIJET/DijetRootTreeAnalyzer
    $ git clone -b dijetpdf_74X https://github.com/RazorCMS/HiggsAnalysis-CombinedLimit HiggsAnalysis/CombinedLimit
    $ cd HiggsAnalysis/CombinedLimit
    $ scram b -j 4
    $ cd $CMSSW_BASE/CMSDIJET/DijetRootTreeAnalyzer
    ```
    
### Getting mjj histogram for fit
1. Add files from EOS together
    
    ```sh
    $ python python/AddFiles.py -l lists/CaloScoutingHT_2015_JEC_CaloL1L2L3_PFL2L3Residual_BiasCorrected_20160523_161841_reduced_skim.txt -o inputs/data_CaloScoutingHT_Run2015D_BiasCorrected_CaloDijet.root
    ```

### Fitting
1. Copy gg CaloScouting signal resonance files
    
    ```sh
    $ wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/master/ResonanceShapes_gg_13TeV_CaloScouting_Spring15.root -P inputs/
    $ wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/master/ResonanceShapes_gg_13TeV_CaloScouting_Spring15_JERUP.root -P inputs/
    $ wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/master/ResonanceShapes_gg_13TeV_CaloScouting_Spring15_JERDOWN.root -P inputs/
    $ wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/master/ResonanceShapes_gg_13TeV_CaloScouting_Spring15_JESUP.root -P inputs/
    $ wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/master/ResonanceShapes_gg_13TeV_CaloScouting_Spring15_JESDOWN.root -P inputs/
    ```

1. Perform a binned, blinded, background-only, sideband fit (with signal shape plotted)
    
    ```sh
    $ mkdir -p fits_2016_05_23/
    $ mkdir -p fits_2016_05_23/CaloDijet_Sideband/
    $ python python/BinnedFit.py -c config/dijet.config -b CaloDijet
    inputs/data_CaloScoutingHT_Run2015D_BiasCorrected_CaloDijet.root
    --lumi 1918 --fit-region Low,High --plot-region Low,High -d
    fits_2016_05_23/CaloDijet_Sideband/ -s
    inputs/ResonanceShapes_gg_13TeV_CaloScouting_Spring15.root -m gg
    --mass 750 --xsec 15 --fit-spectrum
    ```
	
1. Run and fit 1000 toys and plot GOF for blinded fit
    
    ```sh
    $ python python/RunToys.py -b CaloDijet --freq -c config/dijet.config --lumi 1918 --fit-region Low,High  -d fits_2016_05_23/CaloDijet_Sideband/ -i fits_2016_05_23/CaloDijet_Sideband/DijetFitResults_CaloDijet.root  -t 1000 -s 0
    $ python python/PlotGOF.py -b CaloDijet -c config/dijet.config -d fits_2016_05_23/CaloDijet_Sideband/ -t fits_2016_05_23/CaloDijet_Sideband/toys_Freq_s0_CaloDijet.root -l 1918 --data
    ```

1. Perform an unblinded, background-only, full region fit (with signal shape plotted)
    
    ```sh
    $ mkdir -p fits_2016_05_23/CaloDijet_Full/
    $ python python/BinnedFit.py -c config/dijet.config -b CaloDijet
    inputs/data_CaloScoutingHT_Run2015D_BiasCorrected_CaloDijet.root
    --lumi 1918 --fit-region Full --plot-region Full -d
    fits_2016_05_23/CaloDijet_Sideband/ -s
    inputs/ResonanceShapes_gg_13TeV_CaloScouting_Spring15.root -m gg
    --mass 750 --xsec 15 --fit-spectrum
    ```
	
1. Run and fit 1000 toys and plot GOF for unblinded fit
    
    ```sh
    $ python python/RunToys.py -b CaloDijet --freq -c config/dijet.config --lumi 1918 --fit-region Full  -d fits_2016_05_23/CaloDijet_Full/ -i fits_2016_05_23/CaloDijet_Sideband/DijetFitResults_CaloDijet.root  -t 1000 -s 0
    $ python python/PlotGOF.py -b CaloDijet -c config/dijet.config -d fits_2016_05_23/CaloDijet_Sideband/ -t fits_2016_05_23/CaloDijet_Full/toys_Freq_s0_CaloDijet.root -l 1918 --data
    ```
	
### Setting limits
1. Run combine to produce blinded, expected limits for gg resonance mass range [50, 1600] in 50 GeV steps
    
    ```sh
    $ mkdir cards
    $ python python/RunCombine.py -m gg -d cards/ --mass range\(500,1650,50\) -c config/dijet.config -i fits_2016_05_23/CaloDijet_Sideband/DijetFitResults_CaloDijet.root -b CaloDijet --rMax 20 --xsec 10 -l 1.918 --blind
    ```

1. Get blinded combine output files and plot 1D limit
    ```sh
    $ python python/GetCombine.py -d cards/ -m gg --mass range\(500,1650,50\) -b CaloDijet --xsec 10 -l 1.918
    $ python python/Plot1DLimit.py -d cards/ -m gg -b CaloDijet -l 1.918
    ```
	
1. Run combine to produce unblinded limits for gg resonance mass range [50, 1600] in 50 GeV steps
    
    ```sh
    $ mkdir cards
    $ python python/RunCombine.py -m gg -d cards_unblind/ --mass range\(500,1650,50\) -c config/dijet.config -i fits_2016_05_23/CaloDijet_Full/DijetFitResults_CaloDijet.root -b CaloDijet --rMax 20 --xsec 10 -l 1.918 
    ```

1. Get unblinded combine output files and plot 1D limit
    ```sh
    $ python python/GetCombine.py -d cards_unblind/ -m gg --mass range\(500,1650,50\) -b CaloDijet --xsec 10 -l 1.918
    $ python python/Plot1DLimit.py -d cards_unblind/ -m gg -b CaloDijet -l 1.918
    ```





