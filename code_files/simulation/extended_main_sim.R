# simulate datasets for additional conditions 
# (no deposit, burst deposit, low morphological clock rate & relaxed morphological clock)
run_simulation_extended = function(save.folder) {
  library(FossilSim)
  library(phyclust)
  
  #### Fixed args for all simulations
  args = list(
    # Output parameters
    ntrees = 100,
    save_folder = save.folder,
    seed = 451,
    
    # Tip numbers
    nextant = 25,
    nfossils = 50,
    
    # Molecular parameters
    mol_length = 4500,
    mol_model = "-mHKY -f 0.35 0.16 0.21 0.28 -t 2.33 -a0.35 -g5",
    mol_clock_rate = 5e-3, # from Brunke et al. 2017
    
    # Morphological parameters
    morph_length = 120,
    morph_alpha = 0.55,
    morph_ncats = 5,
    morph_props = c(0.7,0.2,0.1), # these are target proportions, not checked after
    prop_extant_only = 0.05,
    
    # Deposit parameters - 1 = precise-date, 2 = imprecise-date
    prop_undated = 0.1,
    rate1to2 = 0.6,
    rate2to1 = 0.7,
    
    # Fossil parameters
    age_range_mult = 0.1,
    
    # BD parameters
    origin_time = 120,
    spec_rate = 0.05, # from Brunke et al. 2017, range 0.05-0.1
    ext_rate = 0.02, # calibrated for lambda, origin and nextant
    rho = 0.5,
    
    # Imprecise deposit parameters
    undated_min_age = 30,
    undated_max_age = 50 )
  
  for(opt in c("morph_relaxed", "low_morph", "no_deposit", "burst_deposit")) .run.sim(args, opt)
}

.run.sim = function(args, option = c("morph_relaxed", "low_morph", "no_deposit", "burst_deposit")) {
  if(length(option) > 1) option = "morph_relaxed"
  
  set.seed(args$seed)
  full_args = c(args, .dependent_args(opt = option))
  .core_loop(full_args, option)
}

### Arguments dependent on simulation setup
.dependent_args = function(opt) {
  
  if(opt == "low_morph") mcl = 0.01
  else if(opt == "morph_relaxed") mcl = function(n) rexp(n, 10)
  else mcl = 0.1 # from Farrell et al. 2004
  
  args = list(
    morph_clock_rate = mcl,
    
    sampl_rate_up = if(opt == "burst_deposit") 0.2 else 0.04,
    sampl_rate = 0.03, # calibrated for nfossils
    
    name = if(opt == "no_deposit") paste0("_prop_0_age_0.1_", opt) else paste0("_prop_0.1_age_0.1_", opt) 
  )
  
  args
}

### Core simulation
.core_loop = function(args, option = c("morph_relaxed", "low_morph", "no_deposit", "burst_deposit")) {
  
  if(length(option) > 1) option = "morph_relaxed"
  
  name = paste0("DS_seed_", args$seed, args$name)
  dir.create(paste0(args$save_folder, name), showWarnings = F)
  name = paste0(name,"/",name)
  
  # simulating trees and fossils with parameters
  trees = list()
  fossils = list()
  samp_trees = list()
  
  while(length(trees) < args$ntrees) {
    nsim = args$ntrees - length(trees)
    r_phy = r_ext = r_fos = r_prop = 0
    print(paste("Simulations remaining", nsim))
    trees = c(trees, TreeSim::sim.bd.age(args$origin_time, nsim, args$spec_rate, args$ext_rate, complete = T))
    for (i in args$ntrees:(args$ntrees-nsim+1)) {
      if(class(trees[[i]]) != "phylo") { 
        trees = trees[-i] 
        fossils = fossils[-i]
        r_phy = r_phy +1
        next
      }
      
      # filter on number of extant samples
      ext_samples = length(sampled.tree.from.combined(trees[[i]])$tip.label)*args$rho
      if(ext_samples > args$nextant*1.2 || ext_samples < args$nextant*0.8) {
        trees = trees[-i]
        fossils = fossils[-i]
        r_ext = r_ext +1
        next
      }
      
      # add burst deposit and no deposit here
      if(option == "no_deposit") {
        fossils[[i]] = sim.fossils.intervals(tree = trees[[i]], interval.ages = c(0, 130), rates = args$sampl_rate)
        fossils[[i]]$trait = 1
      }
      else if(option == "burst_deposit") {
        start_int = runif(1, 30, 50)
        fssls1 = sim.fossils.intervals(tree = trees[[i]], interval.ages = c(0, 130), rates = args$sampl_rate)
        fssls2 = sim.fossils.intervals(tree = trees[[i]], interval.ages = c(start_int, start_int + 2), rates = args$sampl_rate_up)
        if(length(fssls1$edge) > 0) fssls1$trait = 1
        if(length(fssls2$edge) > 0) fssls2$trait = 2
        fossils[[i]] = rbind(fssls1, fssls2)
      }
      else fossils[[i]] = sim.fossils.intervals(tree = trees[[i]], interval.ages = c(0, args$undated_min_age, args$undated_max_age, 130),
                                                rates = c(args$sampl_rate, args$sampl_rate_up, args$sampl_rate))
      # filter on number of fossils
      if(length(fossils[[i]]$edge) < args$nfossils*0.9 || 
         length(fossils[[i]]$edge) > args$nfossils*1.1) {
        #print(length(fossils[[i]]$edge))
        fossils = fossils[-i]
        trees = trees[-i]
        r_fos = r_fos +1
        next
      }
      
      if(!option %in% c("no_deposit", "burst_deposit")) {
        traits = sim.deposit.values(trees[[i]], c(args$rate1to2, args$rate2to1), args$undated_min_age, args$undated_max_age)
        fossils[[i]] = assign.traits(fossils[[i]], traits)
      }
      
      if(option != "no_deposit") {
        # filter on undated proportion
        undated = which(fossils[[i]]$trait == 2)
        p = length(undated)/length(fossils[[i]]$trait)
        up_tol = 1.1
        low_tol = 0.9
        if(p < args$prop_undated*low_tol || p > args$prop_undated*up_tol) {
          #print(p)
          fossils = fossils[-i]
          trees = trees[-i]
          r_prop = r_prop +1
          next
        }
      }
    }
    print(paste("Rejected for extinction", r_phy, ", for n_extant", r_ext, 
                ", for n_fossils", r_fos, ", for undated prop", r_prop))
  }
  
  mol_seqs = morph_seqs = list()
  for (i in 1:args$ntrees) {
    
    fossils[[i]]$h = (fossils[[i]]$hmin + fossils[[i]]$hmax)/2
    fossils[[i]] = fossils[[i]][order(fossils[[i]]$sp, -fossils[[i]]$h), ]
    
    # adding uncertainty to fossil ages
    undated = which(fossils[[i]]$trait == 2)
    fossils[[i]]$hmax[undated] = args$undated_max_age
    fossils[[i]]$hmin[undated] = args$undated_min_age
    
    if(option != "no_deposit") {
      intervals = sample.intervals(fossils[[i]][-undated,], args$age_range_mult)
      fossils[[i]]$hmax[-undated] = intervals$max
      fossils[[i]]$hmin[-undated] = intervals$min
    }
    else {
      intervals = sample.intervals(fossils[[i]], args$age_range_mult)
      fossils[[i]]$hmax = intervals$max
      fossils[[i]]$hmin = intervals$min
    }
    
    # simulating sequences on the trees
    full = SAtree.from.fossils(trees[[i]],fossils[[i]])
    ftree = full$tree
    fossils[[i]] = full$fossils
    tree = sampled.tree.from.combined(ftree, rho = args$rho)
    samp_trees[[i]] = tree
    
    extant_tips = tree$tip.label[1:(length(tree$tip.label)-length(fossils[[i]]$sp))]
    fossil_tips = tree$tip.label[(length(tree$tip.label)-length(fossils[[i]]$sp)):length(tree$tip.label)]
    
    if(option != "morph_relaxed") {
      morph_seqs[[i]] = sim.morph.seqs(samp_trees[[i]], args$morph_length, args$morph_clock_rate, args$morph_alpha, 
                                     args$morph_ncats, args$morph_props, extant_tips, args$prop_extant_only, 0, 1, fossil_tips[-undated])
    }
    else morph_seqs[[i]] = sim.morph.seqs.relaxed(samp_trees[[i]], args$morph_length, args$morph_clock_rate, args$morph_alpha, 
                                                  args$morph_ncats, args$morph_props, extant_tips, args$prop_extant_only)
    
    mol_seqs[[i]] = sim.mol.seqs(samp_trees[[i]], args$mol_clock_rate, args$mol_length, args$mol_model)
    mol_seqs[[i]] = mol_seqs[[i]][names(mol_seqs[[i]]) %in% extant_tips]
    
    .write.nexus.data(mol_seqs[[i]],file = paste0(args$save_folder, name, "_mol_",i,".nex"))
    .write.nexus.data(morph_seqs[[i]], format = "standard",
                      file = paste0(args$save_folder, name, "_morph_", i, ".nex"))
    write.fossil.ages(samp_trees[[i]], fossils[[i]], file = paste0(args$save_folder, name, "_fossil_ages_", i, ".txt"))
  }
  
  save(trees, fossils, samp_trees, mol_seqs, morph_seqs, file = paste0(args$save_folder, name, ".RData"))
}