plot_sim_results = function(outfolder, plotfolder) {
  library(ggplot2)
  
  df = NULL
  age_ranges = c(0.1, 0.2, 0.3)
  props = c(0.1, 0.3, 0.5)
  names = c("mean_RF", "datedf_med_rel_error", "datedf_coverage", "undatedf_med_rel_error", "undatedf_coverage",
            "datedf_top_error","undatedf_top_error", "undatedf_HPD_width")
  plotnames = c("Robinson-Foulds distance", "Absolute relative error of age \n(precise-date fossils)", "Coverage of age (precise-date fossils)",
                "Absolute relative error of age \n(imprecise-date fossils)", "Coverage of age (imprecise-date fossils)", 
                "Proportion of correct positions \n(precise-date fossils)", "Proportion of correct positions \n(imprecise-date fossils)", 
                "HPD width (imprecise-date fossils)")
  
  for(prop in props) {
    for(age in age_ranges) {
      file_name = paste0(outfolder, "/DS_seed_451_prop_", prop, "_age_", age, "_results.RData")
      load(file_name)
      results = get(paste0("DS_seed_451_prop_", prop, "_age_", age, "_results"))
      for(m in names) {
        df = rbind(df, data.frame(Age = age, Prop = as.character(prop), mean = mean(results[[m]], na.rm = T),
                                  sd = sd(results[[m]], na.rm = T), Measure = m))
      }
    }
  }
  
  cbPalette <- c("#56B4E9", "#CC79A7", "#009E73")
  
  df$lower = sapply(df$mean - df$sd, function(t) max(0,t))
  
  for(i in 1:length(names)) {
    subdf = df[which(df$Measure == names[i]),]
    pl = ggplot(subdf, mapping=aes(x=Prop, y=mean, colour=factor(Age))) + xlab("Proportion of imprecise-date fossils") + ylab(plotnames[i]) +
      geom_pointrange(aes(ymin=lower, ymax=mean+sd), size = 1, position = position_dodge(width=0.5)) 
    
    if(i == 8) { #HPD width
      pl = pl + geom_hline(yintercept = 20, color = "#D55E00", size = 1)
    }
    
    if(i %in% c(3,5)) { #coverage plots
      pl = pl + ylim(0.2, 1.0)
    }
    
    if(i %in% c(5,8)) { #with legend
      pl = pl + theme(text = element_text(size = 15), legend.position = "bottom", legend.direction = "vertical") + 
        scale_color_manual(values=cbPalette, name = "Relative age range of precise-date fossils")
      ggsave(paste0(plotfolder, "/", names[i], ".pdf"), width = 5, height = 6.3)
    }
    else { #no legend
      pl = pl + theme(text = element_text(size = 15), legend.position = "none") + 
        scale_color_manual(values=cbPalette)
      ggsave(paste0(plotfolder, "/", names[i], ".pdf"), width = 5, height = 5)
    }
  }
}

plot_empirical_results = function(outfolder, original_ages, fossils, plotfolder = NULL) {
  library(ggplot2)
  library(ggridges)
  
  n = 2
  prop = c(1, 0.5)
  burnin = 0.1
  orig_table = read.table(original_ages, header = T, stringsAsFactors = F)
  load(fossils)
  
  for(ii in 1:n) {
    intervals = get(paste0("interval_", ii))
    for(p in prop) {
      fl = get(paste0("fossils_", ii, "_", p))
      trees1 = read.incomplete.nexus(file.path(outfolder, paste0("penguins_", ii, "_", p, "_run_1.trees")))
      trees2 = read.incomplete.nexus(file.path(outfolder, paste0("penguins_", ii, "_", p, "_run_2.trees")))
      
      n = length(trees1)
      trees = c(trees1[round(n*burnin):n],trees2[round(n*burnin):n])
      
      est_ages = list()
      for(t in trees) {
        tag = ape::node.depth.edgelength(t)
        tag = max(tag) - tag
        
        for(ff in fl) {
          est_ages[[ff]] = c(est_ages[[ff]], tag[which(t$tip.label == ff)])
        }
      }
      #print(any(est_ages$Palaeospheniscus_patagonicus < 14.6))
      
      min_true = 0.9*intervals$minbound
      max_true = 1.1*intervals$maxbound
      age_range = seq(min_true, max_true, 0.01)
      
      true_ages = lapply(fl, function(ff) {
        id = which(orig_table$taxon == ff)
        sapply(age_range, function(a) {
          if(a < orig_table$min[id] || a > orig_table$max[id]) 0.01 else 1
        })
      })
      names(true_ages) = fl
      
      df = data.frame()
      for(ff in fl) {
        df = rbind(df, data.frame(taxa = gsub("_", " ", ff, fixed = T), type = "observed", age = age_range, height = true_ages[[ff]]))
        df = rbind(df, data.frame(taxa = gsub("_", " ", ff, fixed = T), type = "estimated", age = est_ages[[ff]], height = NA))
      }
      
      df = rbind(df, data.frame(taxa = "Prior", type = "prior", age = age_range, 
                           height = sapply(age_range, function(a) {
                             if(a < intervals$minbound || a > intervals$maxbound) 0.01 else 1
                           })))
      
      int = if(ii == 1) "small" else "large"
      nf = paste("proportion of imprecise-date fossils", p)
      pl = ggplot(df, aes(x = age, y = taxa, color = type, fill = type)) + 
        geom_density_ridges(data = df[df$type == "observed",], stat="identity", alpha = 0.7, scale = 0.6, aes(height=height)) +
        geom_density_ridges(data = df[df$type == "prior",], stat="identity", alpha = 0.7, scale = 0.6, aes(height=height)) +
        geom_density_ridges(data = df[df$type == "estimated",], alpha = 0.7, scale = 0.8, stat = "density", trim = TRUE, aes(height = after_stat(density))) + 
        scale_color_manual(values = c("#56B4E9", "#CC79A7", "#009E73")) + scale_fill_manual(values = c("#56B4E9", "#CC79A7", "#009E73")) +
        scale_x_continuous(expand = c(0, 0)) + scale_y_discrete(expand = expansion(mult = c(0.01, .1))) + 
        theme(axis.text.y = element_text(vjust = 0, face = "italic", size = rel(1.2)), legend.title = element_blank()) +
        ggtitle(paste0("Estimated and observed ages of fossil penguins,\n", int, " interval, ", nf)) + 
        labs(y = "Imprecise-date fossil species") + theme(plot.title = element_text(hjust = 0.5))
      
      if(is.null(plotfolder)) show(pl)
      else ggsave(paste0(plotfolder, "/penguins_", ii, "_", p, ".pdf"), height = 7, width = 9)
    }
  }
}