
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
    
### Getting mjj histogram for fit (not necessary because histogram is in repo)
2. Add files from EOS together
    
    ```sh
    $ python/AddFiles.py -l lists/CaloScoutingHT_2016_JEC_CaloL1L2L3_PFL2L3Residual_NewBiasCorrected_Golden12910pb_20160723_003113_reduced_skim.txt  -o inputs/data_CaloScoutingHT_Run2016BCD_NewBiasCorrected_Golden12910pb_CaloDijet2016.root
    ```

### Copy CaloScouting signal resonance files
3. Copy CaloScouting signal resonance files
    
    ```sh
    $ wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/master/ResonanceShapes_gg_13TeV_CaloScouting_Spring16.root -P inputs/
    $ wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/master/ResonanceShapes_gg_13TeV_CaloScouting_Spring16_JERUP.root -P inputs/
    $ wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/master/ResonanceShapes_gg_13TeV_CaloScouting_Spring16_JERDOWN.root -P inputs/
    $ wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/master/ResonanceShapes_gg_13TeV_CaloScouting_Spring16_JESUP.root -P inputs/
    $ wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/master/ResonanceShapes_gg_13TeV_CaloScouting_Spring16_JESDOWN.root -P inputs/
    ```

### Fitting (for limits, skip this section)
4. Perform an unblinded, background-only, full region fit (with signal shapes plotted)
    
    ```sh
	$ mkdir -p fits_2016_07_24/
    $ mkdir -p fits_2016_07_24/CaloDijet2016_Full
    $ python python/BinnedFit.py -c config/dijet.config -l 12910 --mass 750_1200_1600 -m gg_qg_qq --xsec 9.5_8.2e-1_2.2e-1 -s inputs/ResonanceShapes_gg_13TeV_CaloScouting_Spring16.root inputs/data_CaloScoutingHT_Run2016BCD_NewBiasCorrectedFlat_Golden12910pb_CaloDijet2016.root -b CaloDijet2016 -d fits_2016_07_24/CaloDijet2016_Full/ --fit-spectrum
    ```
	
5. Run and fit 1000 toys and plot GOF for unblinded fit
    
    ```sh
    $ python python/RunToys.py -b CaloDijet2016 --freq -c config/dijet.config --lumi 12910 --fit-region Full  -d fits_2016_07_24/CaloDijet2016_Full/ -i fits_2016_07_24/CaloDijet2016_Full/DijetFitResults_CaloDijet2016.root  -t 1000 -s 0
    $ python python/PlotGOF.py -b CaloDijet2016 -c config/dijet.config -d fits_2016_07_24/CaloDijet2016_Full/ -t fits_2016_07_24/CaloDijet2016_Full/toys_Freq_s0_CaloDijet2016.root -l 12910 --data
    ```
	
### Setting limits
6. Run combine to produce expected and observed limits for gg resonance mass range [500, 1600] in 50 GeV steps 
    
    ```sh
    $ mkdir -p cards_gg_freq
    $ python python/RunCombine.py -m gg -d cards_gg_freq/ --mass range\(500,1650,50\) -c config/dijet.config -i inputs/DijetFitResults.root -b CaloDijet2016 --rMax 20 --xsec 10 -l 12.910
    ```

7. Convert combine output files and plot 1D limit
    ```sh
	$ python python/GetCombine.py -d cards_gg_freq/ -m gg --mass range\(500,1650,50\) -b CaloDijet2016 --xsec 10 -l 12.910
    $ python python/Plot1DLimit.py -d cards_gg_freq/ -m gg -b CaloDijet2016 -l 12.910 --massMin 600 --massMax 1600 --xsecMin 1e-3 --xsecMax 1e5
    ```

### Bias studies with 4-parameter modified exponential 
8. Run combine-based bias studies with 1000 toys for r = 1 and  using 4-parameter modified exponential (defined in confif/dijet_bias.config)
    ```sh
    	$ mkdir signal_bias/
	$ python python/RunBias.py -c config/dijet_bias.config --mass 750 -m gg -d signal_bias/ -r 1 -l 12.910 --xsec 10 -t 1000 --gen-pdf modexp  --fit-pdf fourparam
    ```
9. Make plots of bias
    ```sh
	$ python python/PlotBias.py -c config/dijet_bias.config --mass 750 -m gg -d signal_bias/ -r 1 -l 12.910 --xsec 10 -t 1000 --gen-pdf modexp  --fit-pdf fourparam
    ```


