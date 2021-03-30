source("figures.R")
source("figures_suppl.R")

# To make plots of all results for simulated data (main text)
outfolder = "path/to/results/files/as/RData"
plotfolder = "path/to/folder/to/store/figures"
plot_sim_results(outfolder, plotfolder)

# To make plots of all results for simulated data (additional simulation conditions)
outfolder = "path/to/results/files/as/RData"
plotfolder = "path/to/folder/to/store/figures"
plot_add_sim_results(outfolder, plotfolder)

# To make ridge plots of penguins results
outfolder = "path/to/folder/with/empirical/logfiles"
original_ages = "path/to/taxa_original.tsv"
fossils = "path/to/fossils.RData"
plotfolder = "path/to/folder/to/store/figures"
plot_empirical_results(outfolder, original_ages, fossils, plotfolder)

# To make a ridge plot of an example run
source("read.incomplete.nexus.R")
trees_file = "path/to/example/trees/output"
dataset_file = "path/to/corresponding/simulated/dataset/RData/file"
plotfile = "path/to/file/to/save/figure"
plot_example_ridges(trees_file, dataset_file, plotfile)