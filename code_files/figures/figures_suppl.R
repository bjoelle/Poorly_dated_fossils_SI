plot_add_sim_results = function(outfolder, plotfolder) {
  library(ggplot2)
  
  df = NULL
  age = 0.2
  props = c(0.1, 0.3, 0.5)
  conds = c("", "_ss", "_relaxed", "_mrphdisc") #, "_large"
  cnames = c("Normal", "Less fossils", "Relaxed clock", "Realistic characters") # , "More fossils"
  
  names = c("mean_RF", "datedf_med_rel_error", "datedf_coverage", "undatedf_med_rel_error", "undatedf_coverage",
            "datedf_top_error","undatedf_top_error", "undatedf_HPD_width")
  plotnames = c("Robinson-Foulds distance", "Absolute relative error of age\n(precise-date fossils)", "Coverage of age (precise-date fossils)",
                "Absolute relative error of age \n(imprecise-date fossils)", "Coverage of age (imprecise-date fossils)", 
                "Proportion of correct positions \n(precise-date fossils)", "Proportion of correct positions \n(imprecise-date fossils)", 
                "HPD width (imprecise-date fossils)")
  
  for(prop in props) {
    for(cnd in conds) {
      file_name = paste0(outfolder, "/DS_seed_451_prop_", prop, "_age_", age, cnd, "_results.RData")
      load(file_name)
      results = get(paste0("DS_seed_451_prop_", prop, "_age_", age, cnd, "_results"))
      for(m in names) {
        df = rbind(df, data.frame(Cond = cnames[which(conds == cnd)], Prop = as.character(prop), mean = mean(results[[m]], na.rm = T),
                                  sd = sd(results[[m]], na.rm = T), Measure = m))
      }
    }
  }
  
  # making sure sd interval for coverage <= 1
  for(mc in c("datedf_coverage", "undatedf_coverage")) {
    rows = which(df$Measure == mc)
    for(r in rows) df$sd[r] = min(df$sd[r], 1 - df$mean[r])
  }
  
  cbPalette <- c("#56B4E9", "#CC79A7", "#009E73", "#F0E442", "#D55E00")
  
  df$lower = sapply(df$mean - df$sd, function(t) max(0,t))
  
  for(i in 1:length(names)) {
    subdf = df[which(df$Measure == names[i]),]
    pl = ggplot(subdf, mapping=aes(x=Prop, y=mean, colour=factor(Cond))) + xlab("Proportion of imprecise-date fossils") + ylab(plotnames[i]) +
      geom_pointrange(aes(ymin=lower, ymax=mean+sd), size = 1, position = position_dodge(width=0.5)) 
    
    if(i == 8) { #HPD width
      pl = pl + geom_hline(yintercept = 20, color = "#D55E00", size = 1)
    }
    
    if(i %in% c(3,5)) { #coverage plots
      pl = pl + ylim(0.2, 1.0)
    }
    
    if(i %in% c(5,8)) { #with legend
      pl = pl + theme(text = element_text(size = 15), legend.position = "bottom", legend.direction = "vertical") + 
        scale_color_manual(values=cbPalette, name = "Simulation condition")
      ggsave(paste0(plotfolder, "/", names[i], ".pdf"), width = 5, height = 6.6)
    }
    else { #no legend
      pl = pl + theme(text = element_text(size = 15), legend.position = "none") + 
        scale_color_manual(values=cbPalette)
      ggsave(paste0(plotfolder, "/", names[i], ".pdf"), width = 5, height = 5)
    }
  }
}

plot_example_ridges = function(trees_file, dataset_file, plotfile = NULL, seed = 451, maxn = 10) {
  library(ggplot2)
  library(ggridges)
  
  filenm = strsplit(trees_file,'/')[[1]]
  filenm = filenm[length(filenm)]
  filenm = substr(filenm, 1, nchar(filenm) -6)
  filenm = strsplit(filenm,"_")[[1]]
  
  prop = filenm[5]
  age = filenm[7]
  cond = if(filenm[8] %in% c("ss", "large", "relaxed", "mrphdisc")) filenm[8] else ""
  idx = as.numeric(if(cond == "") filenm[8] else filenm[9])
  
  name = paste0("DS_seed_", seed, "_prop_", prop, "_age_", age, cond)
  load(paste0(dataset_file, "/", name, ".RData"))
  true_tree = samp_trees[[idx]]
  true_fossils = fossils[[idx]]
  
  trees = read.incomplete.nexus(trees_file)
  n = length(trees)
  trees = trees[round(n*0.25):n]
  
  fid = true_tree$tip.label[(length(true_tree$tip.label) - length(true_fossils$sp) +1):length(true_tree$tip.label)][1:maxn]
  
  est_ages = list()
  for(t in trees) {
    tag = ape::node.depth.edgelength(t)
    tag = max(tag) - tag
    
    for(ff in fid) {
      est_ages[[ff]] = c(est_ages[[ff]], tag[which(t$tip.label == ff)])
    }
  }
  
  min_true = 0.9*min(true_fossils$hmin[1:maxn])
  max_true = 1.1*max(true_fossils$hmax[1:maxn])
  age_range = seq(min_true, max_true, 0.01)
  
  sim_ages = lapply(1:length(fid), function(ff) {
    sapply(age_range, function(a) {
      if(a < true_fossils$hmin[ff] || a > true_fossils$hmax[ff]) 0.01 else 1
    })
  })
  names(sim_ages) = fid
  
  true_ages = true_fossils$h[1:maxn]
  names(true_ages) = fid

  df = data.frame()
  for(ff in fid) {
    df = rbind(df, data.frame(taxa = ff, type = "simulated", age = age_range, height = sim_ages[[ff]]))
    df = rbind(df, data.frame(taxa = ff, type = "estimated", age = est_ages[[ff]], height = NA))
  }
  dftrue = data.frame(taxa = as.factor(fid), true = true_ages, type = "simulated")

  nf = paste0("example file, prop = ", prop, ", age = ", age, ", rep = ", idx)
  
  pl = ggplot(df, aes(x = age, y = taxa, color = type, fill = type)) + 
    geom_density_ridges(data = df[df$type == "simulated",], stat="identity", alpha = 0.5, scale = 0.6, aes(height=height)) +
    geom_density_ridges(data = df[df$type == "estimated",], alpha = 0.5, scale = 0.8) +
    geom_segment(data = dftrue, aes(x = true, xend = true+0.1, y = as.numeric(taxa), yend = as.numeric(taxa) + 0.6), color = "red") +
    scale_color_manual(values = c("#56B4E9", "#CC79A7")) + scale_fill_manual(values = c("#56B4E9", "#CC79A7")) +
    scale_x_continuous(expand = c(0, 0)) + scale_y_discrete(expand = expansion(mult = c(0.01, .1))) + 
    theme(axis.text.y = element_text(vjust = 0, face = "italic", size = rel(1.2)), legend.title = element_blank()) +
    ggtitle(paste0("Estimated and simulated ages,\n", nf)) + 
    labs(y = "Fossil species") + theme(plot.title = element_text(hjust = 0.5))
  
  if(is.null(plotfile)) show(pl)
  else ggsave(paste0(plotfile), height = 7, width = 9)
}