## Adapted from ape to read Nexus files missing the "End;" delimiter
## for incomplete/interrupted runs
read.incomplete.nexus = function (file, tree.names = NULL, force.multi = FALSE) {
  X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
  LEFT <- grep("\\[", X)
  RIGHT <- grep("\\]", X)
  if (length(LEFT)) {
    w <- LEFT == RIGHT
    if (any(w)) {
      s <- LEFT[w]
      X[s] <- gsub("\\[[^]]*\\]", "", X[s])
    }
    w <- !w
    if (any(w)) {
      s <- LEFT[w]
      X[s] <- gsub("\\[.*", "", X[s])
      sb <- RIGHT[w]
      X[sb] <- gsub(".*\\]", "", X[sb])
      if (any(s < sb - 1)) 
        X <- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
    }
  }
  endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
  semico <- grep(";", X)
  i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
  i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
  translation <- if (length(i2) == 1 && i2 > i1) 
    TRUE
  else FALSE
  if (translation) {
    end <- semico[semico > i2][1]
    x <- X[(i2 + 1):end]
    x <- unlist(strsplit(x, "[,; \t]"))
    x <- x[nzchar(x)]
    TRANS <- matrix(x, ncol = 2, byrow = TRUE)
    TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
    n <- dim(TRANS)[1]
  }
  start <- if (translation) 
    semico[semico > i2][1] + 1
  else i1 + 1
  end <- if(length(endblock[endblock > i1]) > 0) endblock[endblock > i1][1] - 1 else length(X)
  tree <- X[start:end]
  rm(X)
  tree <- tree[tree != ""]
  semico <- grep(";", tree)
  Ntree <- length(semico)
  if (Ntree == 1 && length(tree) > 1) 
    STRING <- paste(tree, collapse = "")
  else {
    if (any(diff(semico) != 1)) {
      STRING <- character(Ntree)
      s <- c(1, semico[-Ntree] + 1)
      j <- mapply(":", s, semico)
      if (is.list(j)) {
        for (i in 1:Ntree) STRING[i] <- paste(tree[j[[i]]], 
                                              collapse = "")
      }
      else {
        for (i in 1:Ntree) STRING[i] <- paste(tree[j[, 
                                                     i]], collapse = "")
      }
    }
    else STRING <- tree
  }
  rm(tree)
  STRING <- STRING[grep("^[[:blank:]]*tree.*= *", STRING, 
                        ignore.case = TRUE)]
  Ntree <- length(STRING)
  nms.trees <- sub(" *= *.*", "", STRING)
  nms.trees <- sub("^[[:blank:]]*tree[[:blank:]\\*]*", "", 
                   nms.trees, ignore.case = TRUE)
  STRING <- sub("^.*= *", "", STRING)
  STRING <- gsub(" ", "", STRING)
  colon <- grep(":", STRING)
  if (!length(colon)) {
    trees <- lapply(STRING, ape:::.cladoBuild)
  }
  else if (length(colon) == Ntree) {
    trees <- if (translation) 
      lapply(STRING, ape:::.treeBuildWithTokens)
    else lapply(STRING, ape:::.treeBuild)
  }
  else {
    trees <- vector("list", Ntree)
    trees[colon] <- lapply(STRING[colon], ape:::.treeBuild)
    nocolon <- (1:Ntree)[!1:Ntree %in% colon]
    trees[nocolon] <- lapply(STRING[nocolon], ape:::.cladoBuild)
    if (translation) {
      for (i in 1:Ntree) {
        tr <- trees[[i]]
        for (j in 1:n) {
          ind <- which(tr$tip.label[j] == TRANS[, 1])
          tr$tip.label[j] <- TRANS[ind, 2]
        }
        if (!is.null(tr$node.label)) {
          for (j in 1:length(tr$node.label)) {
            ind <- which(tr$node.label[j] == TRANS[, 
                                                   1])
            tr$node.label[j] <- TRANS[ind, 2]
          }
        }
        trees[[i]] <- tr
      }
      translation <- FALSE
    }
  }
  for (i in 1:Ntree) {
    tr <- trees[[i]]
    if (!translation) 
      n <- length(tr$tip.label)
  }
  if (Ntree == 1 && !force.multi) {
    trees <- trees[[1]]
    if (translation) {
      trees$tip.label <- if (length(colon)) 
        TRANS[, 2]
      else TRANS[, 2][as.numeric(trees$tip.label)]
    }
  }
  else {
    if (!is.null(tree.names)) 
      names(trees) <- tree.names
    if (translation) {
      if (length(colon) == Ntree) 
        attr(trees, "TipLabel") <- TRANS[, 2]
      else {
        for (i in 1:Ntree) trees[[i]]$tip.label <- TRANS[, 
                                                         2][as.numeric(trees[[i]]$tip.label)]
        trees <- ape:::.compressTipLabel(trees)
      }
    }
    class(trees) <- "multiPhylo"
    if (!all(nms.trees == "")) 
      names(trees) <- nms.trees
  }
  trees
}