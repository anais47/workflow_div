# AIM: Check files in the provided path exist

check_file_exist <- function(otu, metadata, taxo, trait, phylo){
  
  if(file.exists(paste("./", "input/", otu, ".csv", sep=""))){
    print("OTU file exists in the provided path")
  } else {
    print("ERROR: OTU file does not exist in the provided path")
  }
  
  if(file.exists(paste("./", "input/", metadata, ".csv", sep=""))){
    print("Sample metadata file exists in the provided path")
  } else {
    print("ERROR: Sample metadata file does not exist in the provided path")
  }
  
  if(file.exists(paste("./", "input/", taxo, ".csv", sep=""))){
    print("Taxonomic assignment file exists in the provided path")
  } else {
    print("ERROR: Taxonomic assignment file doest not exist in the provided path")
  }
  
  if(!is.null(paste("./", "input/", trait, ".csv", sep="")) &
     file.exists(paste("./", "input/", trait, ".csv", sep=""))){
    print("Trait table file exists in the provided path")
  } else {
    print("ERROR: Trait table file does not exist in the provided path --> If the user did not provide a trait table file it is OK, if not, there is an issue")
  }
  
  if(!is.null(paste("./", "input/", phylo, ".csv", sep="")) &
     file.exists(paste("./", "input/", phylo, ".csv", sep=""))){
    print("Phylogenetic table file exists in the provided path")
  } else {
    print("ERROR: Phylogenetic tree file does not exist in the provided path --> If the user did not provide a phylogenetic tree it is OK, if not, there is an issue")
  }
  
}

