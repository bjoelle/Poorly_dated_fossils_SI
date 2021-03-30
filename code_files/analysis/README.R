source("auxf.R")
source("sub_analysis.R")
source("read.incomplete.nexus.R")
source("full_analysis.R")

# Run the full analysis of RevBayes output for simulated datasets
# skip.errors = T will skip runs with low ESS
outputf = "path/to/folders/of/RevBayes/logfiles"
datasetsf = "path/to/folder/of/datasets/stored/as/RData"
run_full_analysis(outputf, datasetsf, skip.errors = T)