# AIM: Get the taxonomic classification  

get_taxo_class <- function(taxo.table.in){
  
  # Select the column with the taxonomic assignment assigned of each OTU/ASV
  taxaList <- taxo.table$user_supplied_name
  
  # Get the unique list to save time
  taxaList.uniq <- unique(taxaList) 
  
  # --- Resolve taxonomic name:   
  # (ex: get the correct spelling when fuzzy spelling, 
  # remove non-scientific name information (spp., subspecies x, ...))
  # This is done by comparing with WORMS and GBIF database :
  # First we look if the taxa is present in WORMS and if not, we go to look in GBIF
  
  # we get the id of worms and gbif databases
  sources <- gnr_datasources()
  #eol <- sources$id[sources$title == 'Encyclopedia of Life' ]
  worms <- sources$id[sources$title == 'World Register of Marine Species']
  gbif <- sources$id[sources$title == 'GBIF Backbone Taxonomy' ]
  
  
  good_name_tmp <- gnr_resolve(taxaList.uniq, data_source_ids=c(worms, gbif), with_canonical_ranks=TRUE)
  # with_canonical_ranks: logical, Returns names with infraspecific ranks, if present. If TRUE, we force canonical=TRUE, otherwise this parameter would have no effect. Default: FALSE
  
  good_name <- good_name_tmp %>% 
    group_by(user_supplied_name) %>%
    spread(data_source_title, matched_name2) %>%
    dplyr::rename(worms='World Register of Marine Species', 
           gbif='GBIF Backbone Taxonomy') %>%
    mutate(query=case_when(!is.na(worms) ~ worms,
                           (is.na(worms) ~ gbif))) %>%
    mutate(check=case_when(!is.na(worms) ~ "matched name found (WORMS)", 
                           is.na(worms) ~ "matched name found (GBIF)")) %>%
    dplyr::select(user_supplied_name, query, check) 
  
  # We add the taxon for which the taxonomic name remains unresolved 
  not_checked <- as.data.frame(attributes(good_name_tmp)$not_known) 
  colnames(not_checked)[1] <- "query"
  not_checked$user_supplied_name <- not_checked$query
  not_checked$check <- "matched name NOT found"
  
  # This is the final object with two columns: 
  # query : name found after resolving taxonomic name or name not resolved
  # check : explanation 
  
  good_name_OK <- rbind(good_name, not_checked)
  
  # --- Retrieve the taxonomic classification for a given taxon
  # We retrieve from the package worrms the taxonomy of each query by putting marine_only=FALSE to   increase the number of species for which the taxonomy is found
  # NB (I used worrms and not taxize with db = 'worms' because I want to retrieve the accepted name and   its taxonomy and it looks like with taxize we only retrieve the taxonomy classification no matter if   the name of the taxa is the accepted one or not
  
  query <- good_name_OK$query
  user_supplied_name <- good_name_OK$user_supplied_name
  
  taxo_worms <- wm_records_names(name = query, marine_only=FALSE)
  
  # We create the object that will contain the worms classification
  taxo_worms_class <- NULL
  
  
  for (taxa in taxo_worms) {
    # if the the tibble from the list returned by wm_records is not with 0 rows --> ie not empty
    if (dim(taxa)[1] != 0){
      
      # we select only the first row 
      # (several taxonomy can be found for a taxa but we select the first row as most likely the more   used)
      r <- as.data.frame(taxa[1,])
      
      # If the taxa name is accepted we select the following columns
      
      if (r$status=="accepted") {
        
        rr <- r %>% select( valid_AphiaID, valid_name, kingdom, phylum, class, order, family, genus,   rank) %>%
          dplyr::rename(valid_ID=valid_AphiaID)
        
      }
      
      # If the taxa name is unaccepted we look for the taxo. classification 
      # of the matching valid name anf we select the following columns
      
      if (r$status=="unaccepted"){
        
        tmp <- r %>% select(valid_name) %>% pull()
        
        tmp.2 <- wm_records_names(name = tmp, marine_only=FALSE)
        rr <- tmp.2[[1]][1,] %>%
          select(valid_AphiaID, valid_name, kingdom, phylum, class, order, family, genus, rank) %>%
          dplyr::rename(valid_ID=valid_AphiaID)
        
      }
      
      rr$source <- "worms"
    }
    
    # if the the tibble from the list returned by wm_records is with 0 rows 
    # --> ie empty, we build an empty dataset to match with the other dataset
    
    if (dim(taxa)[1] == 0){
      
      rr <- data.frame(valid_ID="not_found_in_worms", 
                       valid_name="not_found_in_worms",
                       kingdom="not_found_in_worms",
                       phylum="not_found_in_worms",
                       class="not_found_in_worms",
                       order="not_found_in_worms",
                       family="not_found_in_worms",
                       genus="not_found_in_worms", 
                       rank="not_found_in_worms",
                       source="not_found_in_worms")
    }
    
    taxo_worms_class <- rbind(taxo_worms_class, rr)
    
  }
  
  # we add the query and user_supplied_name columns 
  taxo_worms_class_ok <- cbind(user_supplied_name, query, taxo_worms_class)
  
  
  # We select the query for which worms did not found any taxonomy
  not_fd_worms <- taxo_worms_class_ok %>% filter(source=="not_found_in_worms") 
  
  # We look if we found a taxonomy in gbif
  taxo_gbif_class <- NULL
  
  
  for (i in 1:nrow(not_fd_worms)){
    
    x <- classification(not_fd_worms$query[i],  db = 'gbif')
    r <- x %>% extract2(1)
    
    if(!is.na(r)){
      
      r$source <- "gbif"
      r$query <- not_fd_worms$query[i]
      r$user_supplied_name <- not_fd_worms$user_supplied_name[i]
      
    } else {
      r <- data.frame(name = "not_found_gbif", rank = "not_found_gbif", id = "not_found_gbif", source =   'none',
                      query=not_fd_worms$query[i],   user_supplied_name=not_fd_worms$user_supplied_name[i]) 
    }
    
    taxo_gbif_class <- rbind(taxo_gbif_class, r)
    
  }
  
  # we modify the df to match with the taxo_worms_class_ok
  taxo_gbif_class_ok <- taxo_gbif_class %>%   
    distinct() %>%
    dplyr::select(-id) %>%
    spread(rank, name)
  
  if("species" %in% colnames(taxo_gbif_class_ok)){
    taxo_gbif_class_ok <- taxo_gbif_class_ok %>%
      mutate(valid_name=case_when(!is.na(species) ~species,
                                  is.na(species) ~ genus, 
                                  is.na(genus) ~ family, 
                                  is.na(family) ~ order,
                                  is.na(order) ~ class,
                                  is.na(class) ~ phylum,
                                  is.na(phylum) ~ kingdom)) %>%
      mutate(rank=case_when(!is.na(species) ~ "Species",
                            is.na(species) ~ "Genus", 
                            is.na(genus) ~ "Family", 
                            is.na(family) ~ "Order",
                            is.na(order) ~ "Class",
                            is.na(class) ~ "Phylum",
                            is.na(phylum) ~ "Kingdom")) %>%
      select(-species, -not_found_gbif) 
  } else {
    taxo_gbif_class_ok <- taxo_gbif_class_ok %>%
      mutate(valid_name=NA) %>%
      select(-not_found_gbif) 
  }
  
  taxa_id <- taxo_gbif_class %>%
    distinct(id, name) %>%
    dplyr::rename(valid_name=name,
           valid_ID = id) 
  
  taxo_gbif_class_ok <- taxo_gbif_class_ok %>%
    left_join(., taxa_id, by="valid_name") 
  
  if(!"species" %in% colnames(taxo_gbif_class_ok)){
    taxo_gbif_class_ok$kingdom = NA
    taxo_gbif_class_ok$phylum = NA
    taxo_gbif_class_ok$class = NA
    taxo_gbif_class_ok$order = NA
    taxo_gbif_class_ok$family = NA
    taxo_gbif_class_ok$genus = NA
    taxo_gbif_class_ok$rank = NA
  }
  
  # We select only the taxa identified with worms and we combined them with the others
  tmp <- taxo_worms_class_ok %>% filter(source=="worms")
  taxo_final <- rbind(taxo_gbif_class_ok, tmp)
  
  
  if(dim(taxo_final)[1] == length(taxaList.uniq)){
    print("All taxa from the otu table were searched for taxonomic classification")
  } else {
    print("Something went wrong because not all taxa from the otu table were searched for taxonomic classification")
  } 
  
  # recreate taxo table 
  taxo.table.class <- left_join(taxo.table, taxo_final, by="user_supplied_name")
  
  return(taxo.table.class)
}
