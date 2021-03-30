# adapted from FossilSim for interval-specific state simulation
sim.deposit.values = function(tree, change.rates, min_range, max_range, init = 1){
  
  times = n.ages(tree)
  root = length(tree$tip.label) + 1
  traits = data.frame(sp = c(), trait = c(), st.time = c(), end.time = c())
  
  aux = function(sp, tree, trait, traits) {
    
    if(sp != root || !is.null(tree$root.edge)) {
      end = times[sp]
      if(sp == root) start = times[root] + tree$root.edge
      else start = times[tree$edge[which(tree$edge[,2]==sp),1]]
      
      # only proceed with changes if in the zone OR after the zone in the wrong trait
      if(start > min_range || end < max_range || (start < min_range && trait == 2)) {
        mvtime = min(start, max_range) - rexp(1, change.rates[trait])
        while(mvtime > end && 
              (mvtime > min_range || trait == 2)) { # same conditions for proceeding with the change
          if(mvtime < min_range) mvtime = min_range # avoid going out of the zone
          traits = rbind(traits, data.frame(sp = sp, trait = trait, st.time = start, end.time = mvtime))
          trait = (trait %% 2) + 1
          start = mvtime
          mvtime = start - rexp(1, change.rates[trait])
        }
        if(end < min_range && trait == 2) {
          traits = rbind(traits, data.frame(sp = sp, trait = trait, st.time = start, end.time = min_range))
          trait = 1
          start = min_range
        }
        traits = rbind(traits, data.frame(sp = sp, trait = trait, st.time = start, end.time = end))
      }
    }
    
    if(sp >= root) {
      desc = tree$edge[which(tree$edge[,1]==sp),2]
      for(d in desc) traits = aux(d, tree, trait, traits)
    }
    return(traits)
  }
  
  traits = aux(root, tree, init, traits)
  traits
}

n.ages <- function(tree){
  
  depth = ape::node.depth.edgelength(tree)
  node.ages = max(depth) - depth
  names(node.ages) <- 1:(tree$Nnode+length(tree$tip))
  
  # adding possible offset if tree fully extinct
  if(!is.null(tree$root.time)) node.ages = node.ages + tree$root.time - max(node.ages)
  
  return(node.ages)
}