load.packages <- function(){
  
  message("Loading required libraries")

  requiredLib <- c("dplyr", "tibble", "tidyr", "knitr", "magrittr","rlang", "fs", # general packages
                 "taxize", "rfishbase", "worrms", # taxonomy-related packages
                 "robis", "rgbif",  "sp", "ggmap", # occurences-related packages
                 "fishtree", "ape", "picante", # phylogeny packages
                 "FD", # functional diversity package
                 "ggVennDiagram",
                 "RColorBrewer" # palettes/colors package
  ) 

  for (lib in requiredLib) {
    if (!require(lib, character.only = TRUE)) {
      install.packages(lib, )
    }
  require(lib, character.only = TRUE)
  }
}
