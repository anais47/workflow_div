# AIM: Check if otu table, taxonomic table and sample metadata table have the same ID (sample ID and OTU ID) 

check_match <- function(otu, metadata, taxo){
  
  if(dim(otu)[2]==dim(taxo)[1]){
    print("Same number of otu in otu table and taxonomic table --> OK !")
  } else {
    print("ERROR: Different number of otu in otu table and taxonomic table --> the user must solve this issue and import again the .csv files !")
  }
  
  if(dim(otu)[1]==dim(metadata)[1]){
    print("Same number of samples in otu table and metadata table --> OK !")
  } else {
    print("ERROR: Different number of samples in otu table and metadata table --> the user must solve this issue and import again the .csv files !")
  }
  
  if(all(colnames(otu)[order(colnames(otu))]==taxo$OTU_ID[order(taxo$OTU_ID)])){
    print("OTU_ID are the same in the otu table and taxonomic table --> OK !")
  } else {
    print("ERROR: OTU_ID are not the same in the otu table and taxonomic table --> the user must solve this issue and import again the .csv files !")
  }
  
  if(all(rownames(otu)[order(rownames(otu))]==rownames(metadata)[order(rownames(metadata))])){
    print("Samples_ID are the same in the otu table and metadata table --> OK !")
  } else {
    print("ERROR: Samples_ID are not the same in the otu table and metadata table --> the user must solve this issue and import again the .csv files !")
  }
  
}