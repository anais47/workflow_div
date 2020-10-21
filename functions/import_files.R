# AIM: Import files in the workflow

import_files <- function(otu, sample.otu.col, metadata, sample.metaD.col, taxo, trait, phylo){
  
  otu.comm <- read.csv2(file=paste("./", "input/", otu, ".csv", sep=""), row.names = sample.otu.col)
  
  sample.metaD <- read.csv2(file=paste("./", "input/", metadata, ".csv", sep=""), row.names = sample.metaD.col)
  
  if(is.null(sample.metaD$decimalLongitude)){
    sample.metaD$decimalLongitude <- "not_provided_by_user"
  }
  
  if(is.null(sample.metaD$decimalLatitude)){
    sample.metaD$decimalLatitude <- "not_provided_by_user"
  }
  
  if(is.null(sample.metaD$year)){
    sample.metaD$year <- "not_provided_by_user"
  }
  
  if(is.null(sample.metaD$eventDate)){
    sample.metaD$eventDate <- "not_provided_by_user"
  }
  
  taxo.table <- read.csv2(file=paste("./", "input/", taxo, ".csv", sep=""))
  
  if(is.null(trait)){
    trait.table.user <- NULL
  } else {
    trait.table.user <- read.csv2(file=paste("./", "input/", trait, ".csv", sep=""), na.strings=c("","NA")) # na.strings=c("","NA") --> to replace empty case wih NA
  }
  
  if(is.null(phylo)){
    phylo.tree.user <- NULL
  } else {
    phylo.tree.user <- ape::read.tree(file=paste("./", "input/", phylo, ".tre", sep="")) # not tested to see if works
  }
  
  df.list <- list("otu.comm"=otu.comm, 
                  "sample.metaD"=sample.metaD,
                  "taxo.table"=taxo.table,
                  "trait.table.user"=trait.table.user,
                  "phylo.tree.user"=phylo.tree.user)
  
  return(df.list)  
}
