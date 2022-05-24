# simulate main datasets & some additional conditions
# (smaller dataset, relaxed molecular clock & complex morphological model)
run_simulation = function(save.folder) {
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
    
    # Morphological parameters
    morph_length = 120,
    morph_clock_rate = 0.1, # from Farrell et al. 2004
    morph_alpha = 0.55,
    morph_ncats = 5,
    morph_props = c(0.7,0.2,0.1), # these are target proportions, not checked after
    
    prop_extant_only = 0.05,
    
    # BD parameters
    origin_time = 120,
    spec_rate = 0.05, # from Brunke et al. 2017, range 0.05-0.1
    ext_rate = 0.02, # calibrated for lambda, origin and nextant
    rho = 0.5,
    
    # Imprecise deposit parameters
    undated_min_age = 30,
    undated_max_age = 50 )
  
  done = rep(F, 18)
  if(file.exists("tmp.RData")) load("tmp.RData")
  k = 1
  for(und_idx in 1:3) {
    for(rge_idx in 1:3) {
      if(!done[k]) .run.sim(args, und_idx, rge_idx)
      done[k] = T
      k = k + 1
      save(done, file = "tmp.RData")
    }
    rge_idx = 2
    if(!done[k]) .run.sim(args, und_idx, rge_idx, subsample = T)
    done[k] = T
    k = k + 1
    save(done, file = "tmp.RData")
    if(!done[k]) .run.sim(args, und_idx, rge_idx, morph_discrep = T)
    done[k] = T
    k = k + 1
    save(done, file = "tmp.RData")
    if(!done[k]) .run.sim(args, und_idx, rge_idx, relaxed = T)
    done[k] = T
    k = k + 1
    save(done, file = "tmp.RData")
  }
}

.run.sim = function(args, und_idx, rge_idx, subsample = F, morph_discrep = F, relaxed = F) {
  set.seed(args$seed)
  full_args = c(args, .dependent_args(und_idx, rge_idx, subsample = subsample, morph_discrep = morph_discrep,
                                      relaxed = relaxed))
  .core_loop(full_args, subsample = subsample, morph_discrep = morph_discrep, relaxed = relaxed)
}

### Arguments dependent on simulation setup
.dependent_args = function(und_prop_idx, age_mult_idx, subsample = F, morph_discrep = F, relaxed = F) {
  
  add = if(subsample) "_ss" else ""
  if(morph_discrep) add = paste0(add, "_mrphdisc")
  if(relaxed) add = paste0(add, "_relaxed")
  
  # Deposit parameters - 1 = precise-date, 2 = imprecise-date
  rate1to2 = c(0.6, 0.8, 1.0)
  rate2to1 = c(0.7, 0.5, 0.4)
  undated_proportions = c(0.1,0.3,0.5)
  
  # Fossil parameters
  sampl_rate = c(0.03, 0.02, 0.01) # calibrated for nfossils
  sampl_rate_up = c(0.04, 0.08, 0.15)
  age_range_mult = c(0.1, 0.2, 0.3)
  
  args = list(
    # Morphological model
    unkn_dated_prop = if(morph_discrep) 0.05 else 0,
    unkn_dated_nchars = if(morph_discrep) 0.1 else 1,
    
    prop_undated = undated_proportions[und_prop_idx],
    rate1to2 = rate1to2[und_prop_idx],
    rate2to1 = rate2to1[und_prop_idx],
    
    sampl_rate = sampl_rate[und_prop_idx],
    sampl_rate_up = sampl_rate_up[und_prop_idx],
    age_range_mult = age_range_mult[age_mult_idx],
    
    # Molecular substitution model
    mol_clock_rate = if(!relaxed) 5e-3 else function(n) rexp(n, 200), # from Brunke et al. 2017
    
    name = paste0("_prop_", undated_proportions[und_prop_idx], "_age_", age_range_mult[age_mult_idx], add) )
  
  if(subsample) {
    args$sampl_rate = args$sampl_rate*0.5
    args$sampl_rate_up = args$sampl_rate_up*0.5
    args$nfossils = round(args$nfossils*0.5)
  }
  
  args
}

### Core simulation
.core_loop = function(args, subsample = F, morph_discrep = F, relaxed = F) {
  
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
      
      fossils[[i]] = sim.fossils.intervals(tree = trees[[i]], interval.ages = c(0, args$undated_min_age, args$undated_max_age, 130),
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
      
      traits = sim.deposit.values(trees[[i]], c(args$rate1to2, args$rate2to1), args$undated_min_age, args$undated_max_age)
      fossils[[i]] = assign.traits(fossils[[i]], traits)
      
      # filter on undated proportion
      undated = which(fossils[[i]]$trait == 2)
      p = length(undated)/length(fossils[[i]]$trait)
      up_tol = if(!subsample) 1.1 else 1.2
      low_tol = if(!subsample) 0.9 else 0.8
      if(p < args$prop_undated*low_tol || p > args$prop_undated*up_tol) {
        #print(p)
        fossils = fossils[-i]
        trees = trees[-i]
        r_prop = r_prop +1
        next
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
    
    intervals = sample.intervals(fossils[[i]][-undated,], args$age_range_mult)
    fossils[[i]]$hmax[-undated] = intervals$max
    fossils[[i]]$hmin[-undated] = intervals$min
    
    # simulating sequences on the trees
    tmp = SAtree.from.fossils(trees[[i]],fossils[[i]])
    ftree = tmp$tree
    fossils[[i]] = tmp$fossils
    tree = sampled.tree.from.combined(ftree, rho = args$rho)
    samp_trees[[i]] = tree
    
    extant_tips = tree$tip.label[1:(length(tree$tip.label)-length(fossils[[i]]$sp))]
    fossil_tips = tree$tip.label[(length(tree$tip.label)-length(fossils[[i]]$sp)):length(tree$tip.label)]
    
    morph_seqs[[i]] = sim.morph.seqs(samp_trees[[i]], args$morph_length, args$morph_clock_rate, args$morph_alpha, 
                                     args$morph_ncats, args$morph_props, extant_tips, args$prop_extant_only, 
                                     args$unkn_dated_prop, args$unkn_dated_nchars, fossil_tips[-undated])
    
    mol_seqs[[i]] = sim.mol.seqs(samp_trees[[i]], args$mol_clock_rate, args$mol_length, args$mol_model, relaxed)
    mol_seqs[[i]] = mol_seqs[[i]][names(mol_seqs[[i]]) %in% extant_tips]
    
    .write.nexus.data(mol_seqs[[i]],file = paste0(args$save_folder, name, "_mol_",i,".nex"))
    .write.nexus.data(morph_seqs[[i]], format = "standard",
                      file = paste0(args$save_folder, name, "_morph_", i, ".nex"))
    write.fossil.ages(samp_trees[[i]], fossils[[i]], file = paste0(args$save_folder, name, "_fossil_ages_", i, ".txt"))
  }
  
  save(trees, fossils, samp_trees, mol_seqs, morph_seqs, file = paste0(args$save_folder, name, ".RData"))
}