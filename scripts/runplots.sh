python DrawFromTree_PHYS14_feb_apr.py --var ptHat --xmin 0 --xmax 5000 --bins 5000  --xtitle "#hat{p_{T}}"   --logy --outputDir ../test/ --inputList ../config/list_forPlots_apr2015.txt --inputList_PHYS14_feb15 ../config/list_forPlots_feb2015.txt  --lumi 1000
python DrawFromTree_PHYS14_feb_apr.py --var mjj --xmin 0 --xmax 10000 --bins 1000  --xtitle "m(jj)"   --logy --outputDir ../test/ --inputList ../config/list_forPlots_apr2015.txt  --inputList_PHYS14_feb15 ../config/list_forPlots_feb2015.txt --lumi 1000
python DrawFromTree_PHYS14_feb_apr.py --var pT_j1 --xmin 0 --xmax 4000 --bins 200  --xtitle "p_{T}(j_{1})"   --logy --outputDir ../test/ --inputList ../config/list_forPlots_apr2015.txt --inputList_PHYS14_feb15 ../config/list_forPlots_feb2015.txt --lumi 1000
python DrawFromTree_PHYS14_feb_apr.py --var pT_j2 --xmin 0 --xmax 4000 --bins 200  --xtitle "p_{T}(j_{2})"   --logy --outputDir ../test/ --inputList ../config/list_forPlots_apr2015.txt --inputList_PHYS14_feb15 ../config/list_forPlots_feb2015.txt  --lumi 1000
python DrawFromTree_PHYS14_feb_apr.py --var eta_j1 --xmin 0 --xmax 1.3 --bins 130  --xtitle "#eta(j_{1})"   --outputDir ../test/ --inputList ../config/list_forPlots_apr2015.txt --inputList_PHYS14_feb15 ../config/list_forPlots_feb2015.txt  --lumi 1000
python DrawFromTree_PHYS14_feb_apr.py --var eta_j2 --xmin 0 --xmax 1.3 --bins 130  --xtitle "#eta(j_{2})"   --outputDir ../test/ --inputList ../config/list_forPlots_apr2015.txt --inputList_PHYS14_feb15 ../config/list_forPlots_feb2015.txt  --lumi 1000
python DrawFromTree_PHYS14_feb_apr.py --var deltaETAjj --xmin 0 --xmax 1.3 --bins 130  --xtitle "Delta#eta(jj)"   --outputDir ../test/ --inputList ../config/list_forPlots_apr2015.txt --inputList_PHYS14_feb15 ../config/list_forPlots_feb2015.txt  --lumi 1000

python DrawFromTree.py --var ptHat --xmin 0 --xmax 5000 --bins 5000  --xtitle "#hat{p_{T}}"   --logy --outputDir ../test/ --inputList ../config/list_forPlots_apr2015.txt   --lumi 1000
python DrawFromTree.py --var mjj --xmin 0 --xmax 10000 --bins 1000  --xtitle "m(jj)"   --logy --outputDir ../test/ --inputList ../config/list_forPlots_apr2015.txt  --lumi 1000
python DrawFromTree.py --var pT_j1 --xmin 0 --xmax 4000 --bins 200  --xtitle "p_{T}(j_{1})"   --logy --outputDir ../test/ --inputList ../config/list_forPlots_apr2015.txt   --lumi 1000
python DrawFromTree.py --var pT_j2 --xmin 0 --xmax 4000 --bins 200  --xtitle "p_{T}(j_{2})"   --logy --outputDir ../test/ --inputList ../config/list_forPlots_apr2015.txt    --lumi 1000
python DrawFromTree.py --var eta_j1 --xmin 0 --xmax 1.3 --bins 130  --xtitle "#eta(j_{1})"   --outputDir ../test/ --inputList ../config/list_forPlots_apr2015.txt    --lumi 1000
python DrawFromTree.py --var eta_j2 --xmin 0 --xmax 1.3 --bins 130  --xtitle "#eta(j_{2})"   --outputDir ../test/ --inputList ../config/list_forPlots_apr2015.txt    --lumi 1000
python DrawFromTree.py --var deltaETAjj --xmin 0 --xmax 1.3 --bins 130  --xtitle "Delta#eta(jj)"   --outputDir ../test/ --inputList ../config/list_forPlots_apr2015.txt   --lumi 1000

