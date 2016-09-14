Fisher tests
============

Code to perform Fisher tests and plot multiple background fits on data.

Fisher_test.C
-------------

Set `parameters` to the results of fits to data, change `f_in` to point to histograms for PF RECO and Calo Scouting, and give `GetObject` the correct histogram names. Run with
```
root -q 'Fisher_test.C(0)'
```
for PF RECO or
```
root -q 'Fisher_test.C(1)'
```
for Calo Scouting.

The macro computes the residual sum squared for each background parametrization excluding bins with zero events. It computes the F statistics for successive parametrizations and then the confidence levels. A function with i parameters with CL(i-1)i > 0.05 do not provide a significant improvement over a function with i-1 parameters.

draw_fits.py
------------

Change `inputHistoFileName` to the histograms for PF RECO and Calo Scouting, modify the histogram names used to get `h_data`, and set the parameters for `background2`, etc. to what was used in `Fisher_test.C`. To run, set `dataset` to `0` for PF RECO or `1` for Calo Scouting and give the command
```
python draw_fits.py
```
