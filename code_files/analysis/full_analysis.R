# run full postprocessing on RevBayes output
run_full_analysis = function(outputf, datasetsf, indices = 1:100, seed = 451, skip.errors = F) {
  tmpf = "DS_partial.RData"
  if(file.exists(tmpf)) {
    load(tmpf)
  }
  else {
    measures = c("mcc_RF", "mean_RF", "fossils_med_rel_error", "fossils_mean_rel_error", "fossils_coverage", 
                 "datedf_med_rel_error", "datedf_mean_rel_error", "datedf_coverage", "undatedf_med_rel_error", "undatedf_mean_rel_error", 
                 "undatedf_coverage", "datedf_top_error", "undatedf_top_error", "undatedf_HPD_width", "correlation_dtop_undmeanerr", "correlation_dtop_undcov",
                 "ESS", "nb_NA")
    fulldf = data.frame(prop_undated = c(), age_range_mult = c(), type = c())
    for(m in measures) fulldf[[m]] = c()
    done = c(F)
    
    save(fulldf, done, file = tmpf)
  }
  
  undated_proportions = c(0.1,0.3,0.5)
  age_range_mult = c(0.1, 0.2, 0.3)
  conds = c("_large", "_ss", "_relaxed", "_mrphdisc","")
  
  for(x in undated_proportions) {
    if(!is.na(done[paste(as.character(x))]) && done[as.character(x)]) next
    else done[as.character(x)] = F
    
    for(y in age_range_mult) {
      if(!is.na(done[paste(as.character(x),as.character(y))]) && done[paste(as.character(x),as.character(y))]) next
      else done[paste(as.character(x),as.character(y))] = F
      
      for(z in conds) {
        if(y != age_range_mult[2] && z != "") next
        if(!is.na(done[paste(as.character(x),as.character(y),as.character(z))]) && 
           done[paste(as.character(x),as.character(y),as.character(z))]) next
        else done[paste(as.character(x),as.character(y),as.character(z))] = F
        #if(z == "_large") next
        
        fulldf = rbind(fulldf, .one_folder_analysis(x, y, z, outputf, datasetsf, indices, seed, skip.errors))
        done[paste(as.character(x),as.character(y),as.character(z))] = T
        save(fulldf, done, file = tmpf)
      }
      done[paste(as.character(x),as.character(y))] = T
      save(fulldf, done, file = tmpf)
    }
    done[as.character(x)] = T
    save(fulldf, done, file = tmpf)
  }
  
  #file.remove(tmpf)
  fulldf
}