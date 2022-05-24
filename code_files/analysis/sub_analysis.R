# run analysis for one dataset
.one_folder_analysis = function(name, outputf, datasetsf, indices = 1:100, seed = 451, skip.errors = F) {
  
  
  load(paste0(datasetsf, "/", name, ".RData"))
  print(paste0("Processing : ", name))
  
  measures = c("mcc_RF", "mean_RF", "fossils_med_rel_error", "fossils_mean_rel_error", "fossils_coverage", 
               "datedf_med_rel_error", "datedf_mean_rel_error", "datedf_coverage", "undatedf_med_rel_error", "undatedf_mean_rel_error", 
               "undatedf_coverage", "undatedf_HPD_width", "datedf_top_error", "undatedf_top_error", "mean_RF_extant", 
               "MRCA_rel_error", "extant_MRCA_rel_error", "MRCA_coverage", "extant_MRCA_coverage")
  
  tmpf = paste0(name, "_partial.RData")
  
  nm = strsplit(name,"_")[[1]]
  x = as.numeric(nm[5])
  y = as.numeric(nm[7])
  z = if(length(nm) > 7) paste(nm[8:length(nm)], collapse = "_") else ""
  
  if(file.exists(paste0(name,"_results.RData"))) {
    load(paste0(name,"_results.RData"))
    final = get(paste0(name,"_results"))
    df = data.frame(prop_undated = x, age_range_mult = y, type = if(z == "") "default" else z)
    for(m in c(measures, "ESS", "nb_NA")) df[[m]] = mean(final[[m]], na.rm = T)
    return(df)
  }
  
  if(file.exists(tmpf)) load(tmpf)
  else {
    results = list()
    ESS = c()
  }
  
  for(i in indices) {
    if(length(results) >= i && !is.na(results[[i]])) next
    
    print(paste0("Trace log ",i))
    path = paste0(outputf, "/", name, "/", name, "_", i)
    res = .combine_log_traces(path)
    if(is.na(res$log)) {
      results[[i]] = NA
      next
    }
    ESS[i] = effectiveSize(as.mcmc(res$log$Posterior))
    
    print(paste0("Tree log ",i))
    path = paste0(outputf, "/", name, "/", name, "_", i)
    trees = .combine_tree_traces(path, res$disc)
    if(is.na(trees)) results[[i]] = NA
    else results[[i]] = .calculate_stats(trees, samp_trees[[i]], fossils[[i]])
    save(results, ESS, file = tmpf)
  }
  
  if(!skip.errors && any(is.na(results))) return(NA)
  
  final = list()
  for(m in measures) {
    final[[m]] = sapply(results, function(r) {
      if(is.na(r)) NA else r[[m]] })
  }
  
  final$nb_NA = sum(is.na(results))
  final$ESS = ESS
  
  df = data.frame(prop_undated = x, age_range_mult = y, type = if(z == "") "default" else z)
  for(m in c(measures, "ESS", "nb_NA")) df[[m]] = mean(final[[m]], na.rm = T)
  
  assign(paste0(name,"_results"), final)
  save(list = c(paste0(name,"_results")), file = paste0(name,"_results.RData"))
  df
}

.calculate_stats = function(trees, true_tree, true_fossils) {
  library(coda)
  
  class(true_tree) = "phylo" # because phangorn's type checking is wonky
  mcct = try(phangorn::mcc(trees), silent = FALSE)
  if(class(mcct) == "try-error") {
    print("Error in tree log, rerun")
    return(NA)
  }
  mccrf = phangorn::RF.dist(mcct, true_tree, normalize = T, rooted = T)
  
  meanrf = mean(sapply(trees, function(t){
    phangorn::RF.dist(t, true_tree, normalize = T)
  }))
  
  fid = true_tree$tip.label[(length(true_tree$tip.label) - length(true_fossils$sp) +1):length(true_tree$tip.label)]
  ages = top_error = list()
  true_fclades = .extant_clades(true_tree, fid)
  
  for(t in trees) {
    tag = ape::node.depth.edgelength(t)
    tag = max(tag) - tag
    tclades = .extant_clades(t, fid)
    for(f in fid) {
      ages[[f]] = c(ages[[f]], tag[which(t$tip.label == f)])
      top_error[[f]] = c(top_error[[f]], setequal(true_fclades[[f]], tclades[[f]])) 
    }
  }
  top_error = sapply(top_error, mean)
  
  mean_error = sapply(1:length(ages), function(i) {
    mean(abs(ages[[i]] - true_fossils$h[i]))/true_fossils$h[i]
  })
  med_error = sapply(1:length(ages), function(i) {
    abs(median(ages[[i]]) - true_fossils$h[i])/true_fossils$h[i]
  })
  
  cover_width = lapply(1:length(ages), function(i) {
    HPD = HPDinterval(as.mcmc(ages[[i]]))
    c((true_fossils$h[i] <= max(HPD) && true_fossils$h[i] >= min(HPD)), max(HPD)-min(HPD))
  })
  
  cover = sapply(cover_width, function(vec) vec[1])
  width = sapply(cover_width, function(vec) vec[2])
  
  undated = which(true_fossils$trait == 2)
  dated = which(true_fossils$trait == 1)
  
  result = list(mcc_RF = mccrf, mean_RF = meanrf, fossils_med_rel_error = mean(med_error), fossils_mean_rel_error = mean(mean_error), 
                fossils_coverage = sum(cover)/length(cover), datedf_med_rel_error = mean(med_error[dated]), datedf_mean_rel_error = mean(mean_error[dated]), 
                datedf_coverage = sum(cover[dated])/length(cover[dated]), undatedf_med_rel_error = mean(med_error[undated]), 
                undatedf_mean_rel_error = mean(mean_error[undated]), undatedf_coverage = sum(cover[undated])/length(cover[undated]),
                undatedf_HPD_width = mean(width[undated]), datedf_top_error = mean(top_error[dated]), undatedf_top_error = mean(top_error[undated]))
  
  mrcas = sapply(trees, function(t) max(ape::node.depth.edgelength(t)))
  true_mrca = max(ape::node.depth.edgelength(true_tree))
  result$MRCA_rel_error = abs(median(mrcas) - true_mrca)/true_mrca
  HPD = HPDinterval(as.mcmc(mrcas))
  result$MRCA_coverage = (true_mrca >= min(HPD)) && (true_mrca <= max(HPD))
  
  extinct_tips = true_tree$tip.label[which(n.ages(true_tree)[1:length(true_tree$tip.label)] > 1e-4)]
  extant_true_tree = ape::drop.tip(true_tree, extinct_tips)
  extant_trees = lapply(trees, function(t) ape::drop.tip(t, extinct_tips))
  
  result$mean_RF_extant = mean(sapply(extant_trees, function(t) phangorn::RF.dist(t, extant_true_tree, normalize = T)))
  
  extant_mrcas = sapply(extant_trees, function(t) max(ape::node.depth.edgelength(t)))
  extant_true_mrca = max(ape::node.depth.edgelength(extant_true_tree))
  result$extant_MRCA_rel_error = abs(median(extant_mrcas) - extant_true_mrca)/extant_true_mrca
  HPD = HPDinterval(as.mcmc(extant_mrcas))
  result$extant_MRCA_coverage = (extant_true_mrca >= min(HPD)) && (true_mrca <= max(HPD))
  
  result
}

.combine_tree_traces = function(x, discard, nr = 2, burnin = 0.25) {
  
  paths = c()
  for(i in 1:nr) {
    if(!i %in% discard) paths = c(paths, paste0(x, "_run_", i))
  }
  if(any(!file.exists(paste0(paths, ".trees")))) return(NA)
  
  trees = lapply(paths, function(p) {
    read.incomplete.nexus(paste0(p, ".trees"))
  })
  
  trees = lapply(trees, function(l) {
    n = length(l)
    l[round(n*burnin):n]
  })
  
  tree = trees[[1]]
  if(nr - length(discard) > 1) {
    for(i in 2:length(trees)) tree = c(tree, trees[[i]])
  }
  tree = lapply(tree, .zero_edge_tree)
  
  tree
}

.combine_log_traces = function(x, nr = 2, burnin = 0.25) {
  library(coda)
  
  paths = c()
  for(i in 1:nr) paths = c(paths, paste0(x, "_run_", i))
  
  logs = try(lapply(paths, function(p) {
    read.table(paste0(p,".log"), header = T, stringsAsFactors = F)
  }))
  if(class(logs) == "try-error") return(list(log = NA))
  logs = lapply(logs, function(l) {
    n = length(l$Posterior)
    l[round(n*burnin):n,]
  })
  
  measures = c("Posterior", "Likelihood", "Prior", "origin_time")
  disc = c()
  for(m in measures) {
    meds = sapply(2:length(logs), function(l) {
      abs(median(logs[[l]][[m]]) - median(logs[[1]][[m]]))/median(logs[[1]][[m]])
    })
    if(any(meds > 0.1)) {
      print(paste("Mismatched chains in", x, "on", m))
      return(list(log = NA))
    }
    ESS = sapply(logs, function(l) {
      effectiveSize(as.mcmc(l[[m]]))
    })
    if(any(ESS < 200)) {
      print(paste("Log file", x, "has low ESS", paste(round(ESS, digits = 2), collapse = " & "), "on", m))
      disc = union(disc, which(ESS < 50))
    }
  }
  
  if(length(disc) == nr) return(list(log = NA))
  if(length(disc) > 0) logs = logs[[-disc]]
  if(nr - length(disc) > 1) {
    log = logs[[1]]
    for(i in 2:length(logs)) log = c(log, logs[[i]])
  }
  else log = logs
  return(list(log = log, disc = disc))
}

n.ages = function(tree) {
  ages = ape::node.depth.edgelength(tree)
  max(ages) - ages
}
