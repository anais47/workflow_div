# AIMs: Prepare datasets coming from phylogenetic tree, trait table and community matrix:
# 4 functions: 
#### Traits table and community matrix subsetted on each others
#### Phylogenetic tree and community matrix subsetted on each others
#### Phylogenetic tree, traits table and community matrix subsetted on each others
#### A venn diagram showing the shared and lost taxa after the following steps

## Traits table and community matrix subsetted on each others

traits.comm.pruning <- function(otu.comm.in, type.rank = NULL, taxo.table.in , trait.table.in){
  
  library(dplyr);library(tibble)
  
  # If not the case already, change space between genus and species name with an underscore
  
  trait.table.in$taxon <-  gsub("[[:space:]]", "_", trait.table.in$taxon)
  
    
  # Select from the taxo.table.class only the OTUs identified at the species level 
  otu.sp <- taxo.table.in %>%
      filter(rank=="Species" & class==type.rank) %>% select(OTU_ID, valid_name) %>%
      mutate(valid_name=gsub("[[:space:]]", "_", valid_name)) %>%
      filter(valid_name %in% trait.table.in$taxon)
    
  # Make the community data matrix with the Actinopterygii species 
  # (if there is >=2 OTUs per species, the numbers of reads is summed)
  invisible(capture.output(otu.sp.comm <- otu.comm.in %>% rownames_to_column(var="sample") %>%
                               gather("OTU_ID", "count", -sample) %>%
                               inner_join(otu.sp) %>%
                               # to match with phylo tree labels tips,add the "_" between genus and species
                               mutate(valid_name=gsub("[[:space:]]", "_", valid_name)) %>%
                               group_by(sample, valid_name) %>%
                               summarize(count_sp = sum(count)) %>%
                               spread(valid_name, -sample) %>%
                               column_to_rownames(var="sample")))
    
  # keep only taxon present in both community matrix and traits table
  taxa.trait <- trait.table.in %>%
      filter(taxon %in% colnames(otu.sp.comm)) %>%
      # super important to match with outputs of the match.phylo.data
      column_to_rownames(var="taxon")
    
  # make a list of pruned community matrix and pruned trait table
  df_list <- list("traits"=taxa.trait,
                  "comm"=otu.sp.comm)
    
  return(df_list)
}

## Phylogenetic tree and community matrix subsetted on each others

tree.comm.pruning <- function(otu.comm.in, type.rank, taxo.table.in, phylo.tree.in){
  
  library(dplyr);library(tibble);library(fishtree)
  
  # Select from the taxo.table.class only the Actinopterygii OTUs identified at the species level 
  otu.sp <- taxo.table.in %>%
    filter(rank=="Species" & class==type.rank) %>% dplyr::select(OTU_ID, valid_name)
  
  # Make the community data matrix with the Actinopterygii species 
  # (if there is >=2 OTUs per species, the numbers of reads is summed)
  invisible(capture.output(otu.sp.comm <- otu.comm.in %>% rownames_to_column(var="sample") %>%
                             gather("OTU_ID", "count", -sample) %>%
                             inner_join(otu.sp) %>%
                             # to match with phylo tree labels tips,add the "_" between genus and species
                             mutate(valid_name=gsub("[[:space:]]", "_", valid_name)) %>%
                             group_by(sample, valid_name) %>%
                             summarize(count_sp = sum(count)) %>%
                             spread(valid_name, -sample) %>%
                             column_to_rownames(var="sample")))
  
  # Check if there are species from the community data matrix absent in the phylo tree
  missing.sp <- setdiff(colnames(otu.sp.comm), phylo.tree.in$tip.label)
  
  if(length(missing.sp)==0){
    print("All species in the community data matrix are present in the phylogenetic tree --> No pruning required")
  } else {
    print("Not all species in the community data matrix are present in the phylogenetic tree --> Pruning is required")
    
    # Retrieve all Actinopteri species present in fishtree
    tax <- fishtree_taxonomy(rank = "Actinopteri")
    
    # Create objects for loop
    kno.name <- NULL
    unk.name.pbm <- NULL
    phylo.tree.in2 <- phylo.tree.in
    
    for (sp.withUnderscore in missing.sp){
      
      sp.noUnderscore <- gsub("_", " ", sp.withUnderscore) 
      
      # If the missing species is present in the list of taxa of fishtree --> 
      # it means that this species is known but is not placed into the tree
      
      if(sp.noUnderscore %in% tax$Actinopteri$species){
        
        kno.name <- c(kno.name, sp.withUnderscore)
        
      } else {
        
        # If the missing species is absent in the list of taxa of fishtree --> 
        # it means that perhaps there is a mismatch in the taxonomy --> 
        # So we go to look for unaccepted synonyms in WORMS
        
        unk.name.syn <- wm_synonyms(wm_name2id(name=sp.noUnderscore)) %>% pull(scientificname) 
        
        unk.name.syn.tmp <- unk.name.syn[unk.name.syn %in% tax$Actinopteri$species] 
        # If only one of the found synonyms is found in the taxa of fishtree, this name is replaced with the valid name 
        if(length(unk.name.syn.tmp)==1){
          
          unk.name.syn.tmp <- gsub(" ", "_", unk.name.syn.tmp)
          
          phylo.tree.in2$tip.label[phylo.tree.in2$tip.label == unk.name.syn.tmp] <- sp.withUnderscore
          
        } else {
          
          # If 0 or more than one of the found synonyms are found in the taxa of fishtree, 
          # it means that this species is not taken into account
          
          unk.name.pbm <- c(unk.name.pbm, sp.withUnderscore)
          
        }
      }      
    }
    # Pruning and sorting phylo.tree and community data matrix to match another content
    invisible(capture.output(pruned.tree.onComm <- match.phylo.comm(phylo.tree.in2, otu.sp.comm)))
  }       
  
  # Create a list object with new pruned datasets and objects containing not placed taxa in the tree
  df.list <- list("comm" = pruned.tree.onComm$comm, 
                  "species.not.found.inTree"= kno.name,
                  "phylo.tree" = pruned.tree.onComm$phy,
                  "species.with.supto1.synonym"= unk.name.pbm)
  
  
  print(paste(length(kno.name), "species from the", 
              length(colnames(otu.sp.comm)),
              "species of the community data matrix were still absent in the phylogenetic tree after checking synonyms"))
  
  return(df.list)
  
}

## Phylogenetic tree, traits table and community matrix subsetted on each others

tree.traits.comm.pruning <- function(traits.comm, phylo.comm, phylo.tree, traits.table){
  
  # keep only species with both: phylogenetic data and trait data
  sp.comm <- intersect(colnames(traits.comm), colnames(phylo.comm))
  
  # Subsetted community datasets:
  phylo.comm.subs <- phylo.comm %>% select_at(vars(all_of(sp.comm)))
  traits.comm.subs <- traits.comm %>% select_at(vars(all_of(sp.comm)))
  
  # Subset the phylo.tree based on the species present in both dataset
  invisible(capture.output(pruned.tree <- match.phylo.comm(phylo.tree, phylo.comm.subs)))
  
  invisible(capture.output(pruned.traits <- match.phylo.data(pruned.tree$phy, traits.table)))
  
  # Make a list object 
  df.list <- list("comm"=pruned.tree$comm,
                  "phylo.tree"=pruned.tree$phy,
                  "traits"=pruned.traits$data)
  
  return(df.list)
}


## A venn diagram showing the shared and lost taxa after the following steps

vd_3div_plot <- function(taxo.table.in, type.rank, phylo, traits, otu.comm.in){
  
  library(ggVennDiagram);library(dplyr)
  
  # Select from the taxo.table.class only the Actinopterygii OTUs identified at the species level
  comm_ok <- taxo.table.in %>%
    filter(rank=="Species" & class==type.rank) %>% 
    mutate(valid_name = gsub(" ", "_", valid_name)) %>% pull(valid_name) 
  
  
  obj.toComp <- list(comm = comm_ok, 
                     tree.pruned = colnames(phylo),
                     traits.pruned = colnames(traits))
  
  
  ggVennDiagram(obj.toComp, fill="white") # see https://www.rdocumentation.org/packages/ggVennDiagram/versions/0.3
  
}



vd_3div_overlap <- function(taxo.table.in, type.rank, phylo, traits){
  
  library(VennDiagram);library(dplyr)
  
  # Select from the taxo.table.class only the Actinopterygii OTUs identified at the species level
  comm_ok <- taxo.table.in %>%
    filter(rank=="Species" & class==type.rank) %>% 
    mutate(valid_name = gsub(" ", "_", valid_name)) %>% pull(valid_name) 
  
  
  obj.toComp <- list(comm = comm_ok, 
                     tree.pruned = colnames(phylo),
                     traits.pruned = colnames(traits))
  
  overlap <- calculate.overlap(obj.toComp)
  return(overlap)
  
}

## 

tree.phylo.comm.sel <- function(sp.comm.in, trait.table.in, phylo.tree.in){
  
  # Subset the phylo.tree based on the species present in both dataset
  invisible(capture.output(pruned.tree <- match.phylo.comm(phylo.tree.in, sp.comm.in)))
  
  # Subset trait 
  taxa.trait <- trait.table.in %>%
    rownames_to_column(var="taxon") %>%
    filter(taxon %in% colnames(sp.comm.in)) %>%
    # super important to match with outputs of the match.phylo.data
    column_to_rownames(var="taxon")
  
  # Make a list object 
  df.list <- list("comm"=sp.comm.in,
                  "phylo.tree"=pruned.tree$phy,
                  "traits"=taxa.trait)
  
  return(df.list)
} 