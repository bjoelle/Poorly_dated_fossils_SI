# plot results for additional datasets in supplementary
plot_add_sim_results = function(outfolder, plotfolder, extended = F) {
  library(ggplot2)
  
  df = NULL
  age = if(!extended) 0.2 else 0.1
  props = if(!extended) c(0.1, 0.3, 0.5) else 0.1
  
  if(!extended) {
    conds = c("", "_ss", "_relaxed", "_mrphdisc") #, "_large"
    cnames = c("Normal", "Less fossils", "Relaxed clock", "Realistic characters") # , "More fossils"
  }
  else {
    conds = c("", "_no_deposit", "_burst_deposit", "_low_morph") #, "_morph_relaxed")
    cnames = c("Normal", "No deposit", "Short-time deposit", "Low morphological clock rate") #, "Relaxed morphological clock")
  }
  
  names = c("mean_RF", "datedf_med_rel_error", "datedf_coverage", "undatedf_med_rel_error", "undatedf_coverage",
            "datedf_top_error","undatedf_top_error", "undatedf_HPD_width")
  plotnames = c("Robinson-Foulds distance", "Relative error of age\n(precise-date fossils)", "Coverage of age (precise-date fossils)",
                "Relative error of age \n(imprecise-date fossils)", "Coverage of age (imprecise-date fossils)", 
                "Proportion of correct positions \n(precise-date fossils)", "Proportion of correct positions \n(imprecise-date fossils)", 
                "HPD width (imprecise-date fossils)")
  
  for(cnd in conds) {
    for(prop in props) {
      if(cnd == "_no_deposit") prop = 0
      file_name = paste0(outfolder, "/DS_seed_451_prop_", prop, "_age_", age, cnd, "_results.RData")
      load(file_name)
      results = get(paste0("DS_seed_451_prop_", prop, "_age_", age, cnd, "_results"))
      for(m in names) {
        if(all(is.na(results[[m]]))) df = rbind(df, data.frame(Cond = cnames[which(conds == cnd)], Prop = as.character(prop), values = NA, Measure = m))
        else df = rbind(df, data.frame(Cond = cnames[which(conds == cnd)], Prop = as.character(prop), values = results[[m]][!is.na(results[[m]])], Measure = m))
      }
    }
  }
  
  cbPalette <- c("#56B4E9", "#CC79A7", "#009E73", "#F0E442", "#D55E00")
  
  for(i in 1:length(names)) {
    subdf = df[which(df$Measure == names[i]),]
    subdf$Cond = factor(subdf$Cond, levels = cnames)
    pl = ggplot(subdf, mapping=aes(x=Prop, y=values, colour=factor(Cond))) + xlab("Proportion of imprecise-date fossils") + ylab(plotnames[i]) +
      geom_boxplot(width = 0.5, position = position_dodge(width=0.7))
    
    if(i == 8) { #HPD width
      pl = pl + geom_hline(yintercept = 20, color = "#D55E00", size = 1)
    }
    
    if(i %in% c(3,5)) { #coverage plots
      pl = pl + ylim(0.2, 1.0)
    }
    
    if(i %in% c(6,7)) { #topological error
      pl = pl + ylim(0.0, 1.0)
    }
    
    if(i %in% c(5,8)) { #with legend
      pl = pl + theme(text = element_text(size = 15), legend.position = "bottom", legend.direction = "vertical") + 
        scale_color_manual(values=cbPalette, name = "Simulation condition", drop = F)
      h = 6.6
      ggsave(paste0(plotfolder, "/", names[i], ".pdf"), width = 5, height = h)
    }
    else { #no legend
      pl = pl + theme(text = element_text(size = 15), legend.position = "none") + 
        scale_color_manual(values=cbPalette, drop = F)
      ggsave(paste0(plotfolder, "/", names[i], ".pdf"), width = 5, height = 5)
    }
  }
}

# ridge plot for a simulated run
plot_example_ridges = function(trees_file, dataset_dir, plotfile = NULL, seed = 451, maxn = 10) {
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
  load(paste0(dataset_dir, "/", name, ".RData"))
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
    scale_x_reverse(expand = c(0, 0)) + scale_y_discrete(expand = expansion(mult = c(0.01, .1))) + 
    theme(axis.text.y = element_text(vjust = 0, face = "italic"), legend.title = element_blank(),
          axis.title = element_text(size = 16), axis.text = element_text(size = 14), legend.text = element_text(size = 12)) +
    ggtitle(paste0("Estimated and simulated ages,\n", nf)) + 
    labs(y = "Fossil species") + theme(plot.title = element_text(size = 18, hjust = 0.5))
  
  if(is.null(plotfile)) show(pl)
  else ggsave(paste0(plotfile), height = 7, width = 9)
}

# prior plots for penguins dataset
penguins_prior_plot = function(taxafolder, plotfolder) {
  library(ggplot2)
  library(ggridges)
  
  n = 3
  prop = c(0.5, 1)
  
  for(i in 1:n)  {
    for(p in prop) {
      age_table = read.table(file.path(taxafolder, paste0("taxa_", i, "_", p, ".tsv")), header = T, stringsAsFactors = F)
      age_table = age_table[age_table$min > 1e-3, ]
      
      max_true = 1.1*max(age_table$max)
      age_range = seq(0, max_true, 0.05)
      
      true_ages = lapply(1:length(age_table$taxon), function(id) {
        sapply(age_range, function(a) {
          if(a < age_table$min[id] || a > age_table$max[id]) 0.01 else 1
        })
      })
      names(true_ages) = age_table$taxon
      
      df = data.frame()
      for(ff in 1:length(true_ages)) {
        df = rbind(df, data.frame(taxa = .taxa.name(names(true_ages)[ff]), age = age_range, height = true_ages[[ff]]))
      }
      df$taxa = factor(df$taxa)
      
      int = c("small", "large", "extended")[i]
      nf = paste0(p*100, "% imprecise-date fossils")
      
      pl = ggplot(df, aes(x = age, y = taxa)) + geom_density_ridges(stat="identity", alpha = 0.7, scale = 0.6, aes(height = height), color = "#009E73", fill = "#009E73") +
        scale_x_reverse(expand = c(0, 0)) + scale_y_discrete(expand = expansion(mult = 0.01)) + 
        theme(axis.text.y = element_text(vjust = 0, face = "italic"), legend.title = element_blank(),
              axis.title = element_text(size = 16), axis.text = element_text(size = 14), legend.text = element_text(size = 12)) +
        ggtitle(paste0("Age ranges,\n", int, " interval, ", nf)) + 
        labs(y = "Fossil species") + theme(plot.title = element_text(size = 18, hjust = 0.5))
      
      ggsave(paste0(plotfolder, "/penguins_priors_", i, "_", p, ".pdf"), height = 12, width = 9)
    }
  }
}