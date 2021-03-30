# find all unique monophyletic extant clades & their associated fossils
.extant_clades = function(tree, fid) {
  
  auxf = function(node, tmp) {
    desc = which(tree$edge[,1] == node)
    if(length(desc) == 0) {
      if(tree$tip.label[node] %in% fid) tmp$extinct = tree$tip.label[node]
      else tmp$extant = tree$tip.label[node]
      return(tmp)
    }
    tmp1 = auxf(tree$edge[desc[1],2], list(clades = list()))
    tmp2 = auxf(tree$edge[desc[2],2], list(clades = list()))
    
    tmp$clades = c(tmp1$clades, tmp2$clades)
    tmp$extant = c(tmp1$extant, tmp2$extant)
    if(length(tmp1$extant) > 0 && length(tmp2$extant) > 0) {
      for (t in c(tmp1$extinct, tmp2$extinct)) {
        tmp$clades[[t]] = tmp$extant
        tmp$extinct = c()
      }
    }
    else tmp$extinct = c(tmp1$extinct, tmp2$extinct)
    
    return(tmp)
  }
  
  result = list(extant = c(), extinct = c(), clades = list())
  result = auxf(length(tree$tip.label)+1, result)
  for (t in result$extinct) {
    result$clades[[t]] = "stem"
  }
  result$clades
}

# convert tree from SAs as 2-degree nodes to SAs as tips with zero-length edges
.zero_edge_tree = function(tree) {
  newtips = which(tree$node.label != "")
  ntips = length(tree$tip.label)
  tree$edge[which(tree$edge > ntips)] = tree$edge[which(tree$edge > ntips)] + length(newtips)
  for(t in newtips) {
    tree$tip.label = c(tree$tip.label, tree$node.label[t])
    tree$node.label[t] = ""
    tree$edge = rbind(tree$edge, c(t + length(newtips) + ntips, length(tree$tip.label)))
    tree$edge.length = c(tree$edge.length, 0)
  }
  tree
}