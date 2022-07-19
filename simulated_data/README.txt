This folder contains the data related to the simulated crinoids datasets.

It contains the following folders:
- data_source contains the log files used to calibrate parameters of the simulation
- datasets contains .RData files with the simulated trees, sequences and fossil occurrences for all simulation conditions 
- results contains summaries of the results on the simulated datasets

All files in datasets and results are indexed by the proportion of imprecise-date fossils (prop, 0.1, 0.3 or 0.5), the multiplier used for the age range of precise date fossils (age, 0.1, 0.2 or 0.3) and optionally by the simulation condition.
Simulation conditions are:
 relaxed = relaxed molecular clock
 ss = less fossil samples (subsampling)
 mrphdisc = realistic morphological characters
 no_deposit = without imprecise-date fossils
 burst_deposit = all imprecise-date fossils have the same age
 low_morph = lower morphological clock rate