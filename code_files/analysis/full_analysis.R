# run full postprocessing on RevBayes output
run_full_analysis = function(outputf, datasetsf, indices = 1:100, seed = 451, skip.errors = F) {
  tmpf = "DS_partial.RData"
  if(file.exists(tmpf)) {
    load(tmpf)
  }
  else {
    measures = c("mcc_RF", "mean_RF", "fossils_med_rel_error", "fossils_mean_rel_error", "fossils_coverage", 
                 "datedf_med_rel_error", "datedf_mean_rel_error", "datedf_coverage", "undatedf_med_rel_error", "undatedf_mean_rel_error", 
                 "undatedf_coverage", "datedf_top_error", "undatedf_top_error", "undatedf_HPD_width", "mean_RF_extant", 
                 "MRCA_rel_error", "extant_MRCA_rel_error", "MRCA_coverage", "extant_MRCA_coverage",
                 "ESS", "nb_NA")
    fulldf = data.frame(prop_undated = c(), age_range_mult = c(), type = c())
    for(m in measures) fulldf[[m]] = c()
    done = c(F)
    
    save(fulldf, done, file = tmpf)
  }
  
  datasets = list.files(datasetsf, full.names = T)
  for(dataf in datasets) {
    name = basename(tools::file_path_sans_ext(dataf))
    if(name %in% names(done)) next
    
    fulldf = rbind(fulldf, .one_folder_analysis(name, outputf, datasetsf, indices, seed, skip.errors))
    done[name] = T
    save(fulldf, done, file = tmpf)
  }
  
  #file.remove(tmpf)
  fulldf
}