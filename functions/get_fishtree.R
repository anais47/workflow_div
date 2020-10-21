# AIM: Retrieve a phylogenetic tree from the fihstree packages

get_fishtree <- function(species.target=NULL, rank.target=NULL, type.phylo="chronogram"){
  if (is.null(species.target) & is.null(rank.target)){
    
    phylo.tree <- fishtree::fishtree_phylogeny(type=type.phylo)
    
  } else if (!is.null(rank.target)) {
    
    phylo.tree <- fishtree::fishtree_phylogeny(type=type.phylo, rank = rank.target)
    
  } else {
    
    phylo.tree <- fishtree::fishtree_phylogeny(type=type.phylo, species = species.target)
  }
  
  return(phylo.tree)
}
