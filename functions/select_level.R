# AIM: Aggregate the samples of the community matrix at the sample level without replicates (ex SAMPLE01_R1 --> SAMPLE01)
#      The number of read per taxa is summed or mean at the defined level

agg_replicates <- function(comm.in, metad, no.replicates, method){
  
  library(dplyr);library(rlang)
  
  no.replicates.list <- syms(no.replicates)
  method <- syms(method)
  
  metad <- metad %>% rownames_to_column(var="sample")
  
  comm.no.replicates <- comm.in %>%
    rownames_to_column(var="sample") %>%
    gather("taxon", "count", -sample) %>%
    left_join(metad, by="sample") %>%
    group_by(!!! no.replicates.list, taxon) %>%
    summarize_at(vars(count), funs(!!!method)) %>%
    ungroup() %>%
    spread(taxon, count) %>%
    column_to_rownames(var=no.replicates)
  
  comm.no.replicates <- comm.no.replicates[, colSums(comm.no.replicates != 0 ) > 0]
  
  return(comm.no.replicates)
}


# AIM: Filter some samples by type of samples (can be done for only one column)

filter_samples <- function(comm.in, metad, col, var){
  
  library(dplyr);library(rlang)
  
  col <- syms(col)
  
  metad2 <- metad %>% rownames_to_column(var="sample")
  
  comm.filt <- comm.in %>%
    rownames_to_column(var="sample") %>%
    gather("taxon", "count", -sample) %>%
    left_join(metad2, by="sample") %>%
    filter((!!!col) %in% var) 
  
  comm.filt <- comm.filt[, !names(comm.filt) %in% colnames(metad)]
  
  comm.filt <- comm.filt %>%
    spread(taxon, count) %>%
    column_to_rownames(var="sample")
  
  comm.filt <- comm.filt[, colSums(comm.filt != 0 ) > 0]
  
  return(comm.filt)
}



# AIM: Aggregate the samples of the community matrix at higher level (ex: at the site, season, ...)  
#      The number of read per taxa is summed or mean at the defined level

select_level <- function(sp.comm.in, metad, levels=NULL, method){
  
  library(dplyr);library(rlang)
  
  levels.list <- syms(levels)
  method <- syms(method)
  
  metad <- metad %>% rownames_to_column(var="sample") 
  
  sp.comm.levels <- sp.comm.in %>%
                     rownames_to_column(var="sample") %>%
                     gather("taxon", "count", -sample) %>%
                     left_join(metad, by="sample") %>%
                     group_by(!!! levels.list, taxon) %>%
                     summarize_at(vars(count), funs(!!!method)) %>%
                     ungroup() %>%
                     spread(taxon, count) %>%
                     unite("sample", levels) %>%
                     column_to_rownames(var="sample")

  sp.comm.levels <- sp.comm.levels[, colSums(sp.comm.levels != 0 ) > 0]
                     
  
  return(sp.comm.levels)
}



## OLD VERSION FOR OTU
#select_level <- function(otu.comm.in, metad, levels=NULL){
  
#  library(dplyr);library(rlang)
  
#  levels.list <- syms(levels)
  
#  metad <- metad %>% rownames_to_column(var="sample") 
  
#  otu.comm.levels <- otu.comm.in %>%
#    rownames_to_column(var="sample") %>%
#   gather("OTU_ID", "count", -sample) %>%
#   left_join(metad, by="sample") %>%
#   group_by(!!! levels.list, OTU_ID) %>%
#   summarize(new_count=sum(count)) %>%
#   spread(OTU_ID, new_count) %>%
#   unite("sample", levels) %>%
#   column_to_rownames(var="sample")
  
# otu.comm.levels <- otu.comm.levels[, colSums(otu.comm.levels != 0 ) > 0]
  
  
# return(otu.comm.levels)
#}


