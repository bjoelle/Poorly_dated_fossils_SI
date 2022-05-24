# assign fossils to imprecise/precise categories
assign.traits = function(fossils, traits) {
  aux = function(sp, time) {
    tr = traits[which(traits$sp == sp),]
    if(length(tr$trait) == 1) return(tr$trait[1])
    tr = tr[which(tr$st.time >= time & tr$end.time <= time),]
    if(length(tr$trait) == 1) return(tr$trait[1])
  }
  
  fossils = cbind(fossils, trait = mapply(aux, sp = fossils$edge, time = fossils$hmin))
  fossils
}

# simulate molecular sequences
sim.mol.seqs = function(tree, rate, alignment.length, model, relaxed = F) {
  seqs = list()

  if(!relaxed) new.blens = tree$edge.length * rate
  else {
    rates = rate(length(tree$edge.length))
    new.blens = tree$edge.length * rates
  }
  tree$edge.length = new.blens
  
  seqs = phyclust::seqgen(opts = paste(model, " -l", alignment.length, sep = ""), newick.tree = ape::write.tree(tree)) [-1]
  
  for(i in 1:length(seqs)) {
    tmp = strsplit(seqs[i],split = " ")[[1]]
    names(seqs)[i] = tmp[1]
    seqs[i] = tmp[length(tmp)]
  }
  return(seqs)
}

# simulate character sequences
sim.morph.seqs = function(tree, nchars, rate, alpha, ncats, prop_nstates, extant_tips = c(),
                          ext_only = 0, unkn_prop = 0, unkn_nch = 1, pot_unknw_tips = c()) {
  seq = NULL
  
  while(is.null(seq) || ncol(seq) < nchars) {
    sim_chars = nchars - ifelse(is.null(seq), 0, ncol(seq))
    cat_rates = rgamma(ncats, alpha, rate = alpha/rate)
    rates = cat_rates[sample(ncats, sim_chars, replace = T)]
    
    # all characters assumed symmetric
    mats = lapply(rates, function(x) {
      ns = sample(length(prop_nstates),1, prob = prop_nstates) + 1
      .make.rate.matrix(x, ns)
    })
    seq = cbind(seq, geiger::sim.char(tree, mats, model = "discrete")[,,1])
    
    # filter unchanged chars
    for(l in ncol(seq):1) {
      if(all(seq[,l] == seq[1,l])) seq = seq[,-l]
    }
  }
  
  seq = apply(seq, 2, function(x) { x - 1 })
  
  # characters only sampled for extant tips
  if(ext_only > 0) {
    ext_chars = .sample.prob(1:nchars, ext_only)
    seq[!row.names(seq) %in% extant_tips, ext_chars] = "?"
  }
  
  # fraction of dated fossils with almost no data
  if(unkn_prop > 0) {
    unkn_tips = .sample.prob(pot_unknw_tips, unkn_prop)
    for(tip in unkn_tips) {
      idxes = which(seq[tip, ] != "?")
      kept = sample(idxes, round(nchars*unkn_nch))
      seq[tip, -kept] = "?"
    }
  }
  
  seq = apply(seq, 1, function(x) paste0(x, collapse = ""))
  seq
}

.make.rate.matrix = function(rate, nstates) {
  m = matrix(rate/(nstates-1), nrow = nstates, ncol = nstates)
  for(i in 1:nstates) {
    m[i,i] = -rate
  }
  m
}

# sample age ranges for precise-date fossils
sample.intervals = function(fossils, range_mult) {
  intervals = list(min = c(), max = c())
  for (i in 1:length(fossils$sp)) {
    span = range_mult*fossils$h[i]
    tp = runif(1,0,span)
    intervals$min[i] = max(fossils$h[i] - tp, 0)
    intervals$max[i] = fossils$h[i] + span - tp
  }
  intervals
}

.sample.prob = function(vect, prob) {
  keep = sapply(vect, function(x) if(runif(1) < prob) T else F)
  vect[keep]
}

write.fossil.ages = function(tree, fossils, file) {
  df = data.frame(taxon = tree$tip.label, 
                  min = c(rep(0, length(tree$tip.label) - length(fossils$sp)), fossils$hmin),
                  max = c(rep(0, length(tree$tip.label) - length(fossils$sp)), fossils$hmax))
  write.table(df, file = file, sep = "\t", row.names = F, col.names = T, quote = F)
}

# adapted from ape
.write.nexus.data = function (x, file, format = "dna", datablock = TRUE, interleaved = FALSE, 
                              charsperline = NULL, gap = NULL, missing = NULL) 
{
  format <- match.arg(toupper(format), c("DNA", "PROTEIN", 
                                         "STANDARD", "CONTINUOUS"))
  if (inherits(x, "DNAbin") && format != "DNA") {
    format <- "DNA"
    warning("object 'x' is of class DNAbin: format forced to DNA")
  }
  if (inherits(x, "AAbin") && format != "PROTEIN") {
    format <- "PROTEIN"
    warning("object 'x' is of class AAbin: format forced to PROTEIN")
  }
  indent <- "  "
  maxtax <- 5
  defcharsperline <- 80
  defgap <- "-"
  defmissing <- "?"
  if (is.matrix(x)) {
    if (inherits(x, "DNAbin")) 
      x <- as.list(x)
    else {
      xbak <- x
      x <- vector("list", nrow(xbak))
      for (i in seq_along(x)) x[[i]] <- xbak[i, ]
      names(x) <- rownames(xbak)
      rm(xbak)
    }
  }
  ntax <- length(x)
  nchars <- nchar(x[[1]])
  zz <- file(file, "w")
  if (is.null(names(x))) 
    names(x) <- as.character(1:ntax)
  fcat <- function(..., file = zz) cat(..., file = file, sep = "", 
                                       append = TRUE)
  find.max.length <- function(x) max(nchar(x))
  print.matrix <- function(x, dindent = "    ", collapse = "") {
    Names <- names(x)
    printlength <- find.max.length(Names) + 2
    if (!interleaved) {
      for (i in seq_along(x)) {
        sequence <- paste(x[[i]], collapse = collapse)
        taxon <- Names[i]
        thestring <- sprintf("%-*s%s%s", printlength, 
                             taxon, dindent, sequence)
        fcat(indent, indent, thestring, "\n")
      }
    }
    else {
      ntimes <- ceiling(nchars/charsperline)
      start <- 1
      end <- charsperline
      for (j in seq_len(ntimes)) {
        for (i in seq_along(x)) {
          sequence <- paste(x[[i]][start:end], collapse = collapse)
          taxon <- Names[i]
          thestring <- sprintf("%-*s%s%s", printlength, 
                               taxon, dindent, sequence)
          fcat(indent, indent, thestring, "\n")
        }
        if (j < ntimes) 
          fcat("\n")
        start <- start + charsperline
        end <- end + charsperline
        if (end > nchars) 
          end <- nchars
      }
    }
  }
  if (inherits(x, "DNAbin") || inherits(x, "AAbin")) 
    x <- as.character(x)
  fcat("#NEXUS\n[Data written by write.nexus.data.R, ", date(), 
       "]\n")
  NCHAR <- paste("NCHAR=", nchars, sep = "")
  NTAX <- paste0("NTAX=", ntax)
  DATATYPE <- paste0("DATATYPE=", format)
  if (is.null(charsperline)) {
    if (nchars <= defcharsperline) {
      charsperline <- nchars
      interleaved <- FALSE
    }
    else charsperline <- defcharsperline
  }
  if (is.null(missing)) 
    missing <- defmissing
  MISSING <- paste0("MISSING=", missing)
  if (is.null(gap)) 
    gap <- defgap
  GAP <- paste0("GAP=", gap)
  INTERLEAVE <- if (interleaved) 
    "INTERLEAVE=YES"
  else "INTERLEAVE=NO"
  if (datablock) {
    fcat("BEGIN DATA;\n")
    fcat(indent, "DIMENSIONS ", NTAX, " ", NCHAR, ";\n")
    if (format != "STANDARD") {
      fcat(indent, "FORMAT", " ", DATATYPE, " ", MISSING, 
           " ", GAP, " ", INTERLEAVE, ";\n")
    }
    else {
      fcat(indent, "FORMAT", " ", DATATYPE, " ", MISSING, 
           " ", GAP, " ", INTERLEAVE, " symbols=\"0123456789\";\n")
    }
    fcat(indent, "MATRIX\n")
    if (format != "CONTINUOUS") {
      print.matrix(x)
    }
    else {
      print.matrix(x, collapse = "\t")
    }
    fcat(indent, ";\nEND;\n\n")
  }
  close(zz)
}

# simulate character sequences - relaxed
sim.morph.seqs.relaxed = function(tree, nchars, rate, alpha, ncats, prop_nstates, extant_tips = c(), ext_only = 0) {
  seq = NULL
  
  rates = rate(length(tree$edge.length))
  new.blens = tree$edge.length * rates
  tree$edge.length = new.blens
  
  while(is.null(seq) || ncol(seq) < nchars) {
    sim_chars = nchars - ifelse(is.null(seq), 0, ncol(seq))
    cat_rates = rgamma(ncats, alpha, rate = alpha)
    rates = cat_rates[sample(ncats, sim_chars, replace = T)]
    
    # all characters assumed symmetric
    mats = lapply(rates, function(x) {
      ns = sample(length(prop_nstates),1, prob = prop_nstates) + 1
      .make.rate.matrix(x, ns)
    })
    seq = cbind(seq, geiger::sim.char(tree, mats, model = "discrete")[,,1])
    
    # filter unchanged chars
    for(l in ncol(seq):1) {
      if(all(seq[,l] == seq[1,l])) seq = seq[,-l]
    }
  }
  
  seq = apply(seq, 2, function(x) { x - 1 })
  
  # characters only sampled for extant tips
  if(ext_only > 0) {
    ext_chars = .sample.prob(1:nchars, ext_only)
    seq[!row.names(seq) %in% extant_tips, ext_chars] = "?"
  }
  
  seq = apply(seq, 1, function(x) paste0(x, collapse = ""))
  seq
}