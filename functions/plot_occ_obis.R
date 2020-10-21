### AIM: Plot occurences of species with robis

# Three options are possible to gather the occurences data based on a specific region:
# areaid = c("id", "id", ...) # targeting an already defined region (country, LME, ...)
# areaid_chosen <- c("40024", "40022") : ex for Celtic-Biscay Shelf and North Sea	from large marine ecoregion dataset
# areaid = NULL # targeting all occurences over the world
# area_polygon # region defined by the user with GPS coordinates
# ex: areapolyg_chosen <- "POLYGON ((-4 48.5, -4 48.8, -3.7 48.8, -3.7 48.5, -4 48.5))" #https://obis.org/maptool/#"

## To look which area are available, do: 
# areas <- area()
# View(areas)
# The user must choose the targeted area and record the id

# selection of attributes of records, at least,
# we kept GPS coordinates ("decimalLongitude", "decimalLatitude"), 
# depth, year of collection ("year) and abundance ("individualCount").

# FUNCTION 1 ----
# Aim: Retrieve occurences of targeted marine species from obis

get.occ.species.obis <- function(species.name, 
                                 startdate_chosen = NULL, 
                                 type_zone, coord_zone,
                                 otu.comm.target, sample.metaD.target, taxo.table.class.target,
                                 save.occ.csv="no",
                                 pathToCsv, 
                                 nameCsv){
  
  # Change the data for robis package : year to year-month-day
  if(!is.null(startdate_chosen)){
    
    startdate_chosen=as.Date(paste(startdate_chosen, "-01-01", sep=""))
  }
  
  # Check species in present in the database
  
  if(dim(robis::occurrence(species.name))[1] == 0){
    
    stop("Species not found in the robis database")
    
  } else {
    
    # Function which allows the user to retrieve occurences data 
    # from obis by selecting the arguments of its choice (from area or from polygon)
    
    if(type_zone== "area"){
      occ.obis <- robis::occurrence(species.name, 
                                    areaid = coord_zone, 
                                    startdate = startdate_chosen)
      
    } else if(type_zone == "polygon"){
      occ.obis <- robis::occurrence(species.name, 
                                    geometry = coord_zone,
                                    startdate = startdate_chosen) 
      
    } else {
      stop("type_zone must be one of 'area' or 'polygon'",
           call. = FALSE)
    }
    
    # Check if data had been found in the targeted area for the targeted species in robis 
    
    if(dim(occ.obis)[1] == 0){
      
      stop("Species not present in the chosen geographic region")
      
    } else {
      # Retrieval of the presence of species in the dataset of the user
      # site - decimalLongitude - decimalLatitude
      
      metad <- sample.metaD %>% rownames_to_column(var="sample") 
      
      occ.study_user <- otu.comm %>%
        rownames_to_column(var="sample") %>%
        gather("OTU_ID", "counts", -sample) %>%
        # join with the taxonomy table
        left_join(taxo.table.class, by="OTU_ID") %>%
        # join with the sample metadata table containing the GPS coordinates
        left_join(metad, by="sample") %>%
        # summarize number of reads at the sampling site level unit for each taxa
        group_by(Site, valid_name, decimalLongitude, decimalLatitude, year, eventDate) %>%
        summarise(counts_taxa_site=sum(counts)) %>%
        # removal of absence 
        filter(counts_taxa_site>0) %>%
        # add the origin of the occurences
        mutate(origin_occ="Study user data") %>%
        # rename to match the name of the query taxa with the one of robis::occurences()
        dplyr::rename(scientificName=valid_name) %>%
        mutate(recordedBy="Study user data",
               individualCount=counts_taxa_site)
      
      
      # Make a dataframe containing presence data from the used database and from the user study for database
      
      sp2plot.db <- occ.obis %>% 
        #filter(scientificName==sp2plot) %>% 
        select(recordedBy, scientificName,  
               decimalLongitude, decimalLatitude, 
               year, eventDate, individualCount) %>%
        #select(scientificName, decimalLongitude, decimalLatitude) %>%
        # add the origin of the occurences
        mutate(origin_occ="Public data")
      
      ## for user study
      sp2plot.user <- occ.study_user %>%
        filter(scientificName==species.name) %>%
        as.data.frame() %>%
        select(recordedBy, scientificName,  
               decimalLongitude, decimalLatitude, 
               year, eventDate, individualCount, origin_occ) #%>%
      #select(scientificName, decimalLongitude, decimalLatitude, origin_occ)
      
      # combine databaset and user GPS coordinates for the selected species
      sp2plot.db_user <- rbind(sp2plot.db, sp2plot.user)
      
      # add color based on origin
      sp2plot.db_user <- sp2plot.db_user %>%
        mutate(colors=case_when(origin_occ=="Public data" ~ "#ee3300",
                                origin_occ=="Study user data" ~ "#86b300"))
      
      if(save.occ.csv=="yes"){
      write.csv2(sp2plot.db_user, file=paste(pathToCsv, nameCsv, ".csv", sep=""))
      }
      
      return(sp2plot.db_user)
      
    }
  }
}

# FUNCTION 2 ----
plot.occ.species.obis <- function(data.for.plot){
  
  labels <- data.for.plot$origin_occ %>% lapply(htmltools::HTML)
  
  leaflet::leaflet() %>%
    leaflet::addTiles("http://{s}.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}.png") %>%
    leaflet::addCircleMarkers(popup = paste( "<b>Species: </b>",
                                             data.for.plot$scientificName, 
                                             "<br/>",
                                             "<b>Origin: </b>",
                                             data.for.plot$origin_occ,
                                             "<br/>",
                                             "<b>RecordedBy: </b>",
                                             data.for.plot$recordedBy,
                                             "<br/>",
                                             "<b>Year: </b>",
                                             data.for.plot$year, 
                                             "<br/>",
                                             "<b>EventDate: </b>",
                                             data.for.plot$eventDate,
                                             "<br/>",
                                             "<b>Count: </b>",
                                             data.for.plot$individualCount),
                              lat = data.for.plot$decimalLatitude, 
                              lng = data.for.plot$decimalLongitude, 
                              radius = 3.5, 
                              weight = 0, 
                              fillColor = data.for.plot$colors, 
                              fillOpacity = 1, label = labels) 
}     
      
# FUNCTION 3 ----

occ.genus.obis <- function(genus,
                           startdate_chosen = NULL, 
                           type_zone, coord_zone){
  
  # Change the data for robis package : year to year-month-day
  if(!is.null(startdate_chosen)){
    
    startdate_chosen=as.Date(paste(startdate_chosen, "-01-01", sep=""))
  }
  
  if(dim(robis::checklist(genus))[1] == 0){
    
    stop("Genus not found in the robis database")
    
  } else {
    
    # Function which allows the user to retrieve occurences data 
    # from obis by selecting the arguments of its choice (from area or from polygon)
    
    if(type_zone== "area"){
    chek.gen <- robis::checklist(genus, 
                                 areaid = coord_zone, 
                                 startdate = startdate_chosen) 
      
    } else if(type_zone == "polygon"){
    chek.gen <- robis::checklist(genus, 
                                 geometry = coord_zone,
                                 startdate = startdate_chosen) 
    }
  }
  
  if(dim(chek.gen)[1] == 0){
      
    stop("Genus not present in the chosen geographic region")
      
  } else {
    
    chek.gen <- chek.gen %>% select(scientificName)
      
    return(chek.gen)
    
  }
}
  

# FUNCTION 4 ----
plot.occ.genus.obis <- function(genus,
                                startdate_chosen = NULL, 
                                type_zone, coord_zone){
  
  # Change the data for robis package : year to year-month-day
  if(!is.null(startdate_chosen)){
    
    startdate_chosen=as.Date(paste(startdate_chosen, "-01-01", sep=""))
  }
  
  # Check species in present in the database
  
  if(dim(robis::occurrence(genus))[1] == 0){
    
    stop("Genus not found in the robis database")
    
  } else {
    
    # Function which allows the user to retrieve occurences data 
    # from obis by selecting the arguments of its choice (from area or from polygon)
    
    if(type_zone== "area"){
      occ.obis <- robis::occurrence(scientificname=genus, 
                                    areaid = coord_zone, 
                                    startdate = startdate_chosen)
      
    } else if(type_zone == "polygon"){
      occ.obis <- robis::occurrence(scientificname=genus, 
                                    geometry = coord_zone,
                                    startdate = startdate_chosen) 
      
    } else {
      stop("type_zone must be one of 'area' or 'polygon'",
           call. = FALSE)
    }
    
    # Check if data had been found in the targeted area for the targeted species in robis 
    
    if(dim(occ.obis)[1] == 0){
      
      stop("Genus not present in the chosen geographic region")
      
    } else {

## ADD color to each species of the genus

# stop the code if more than 12 species inside the genus are found 
# --> max number of color for the palette
# would be unreadable anyway

if(length(unique(occ.obis$scientificName)) >= 12){
  stop("Number of species to plot exceed the number of color possible (n=12)")
} else {
  
  if(length(unique(occ.obis$scientificName)) < 3 ){
    col <- c("#A6CEE3", "#B2DF8A")
    col <- col[1:length(unique(occ.obis$scientificName))]
    names(col) <- unique(occ.obis$scientificName)
    
    occ.obis$color <- col[occ.obis$scientificName]
    
  }else{
    
    col <- RColorBrewer::brewer.pal(n=length(unique(occ.obis$scientificName)), name="Paired")
    names(col) <- unique(occ.obis$scientificName)
    occ.obis$color <- col[occ.obis$scientificName]
    
  }
}

sp2plot <- occ.obis %>% 
  select(recordedBy, scientificName,  
         decimalLongitude, decimalLatitude, 
         year, eventDate, color) 

labels <- sp2plot$scientificName %>% lapply(htmltools::HTML)

leaflet::leaflet() %>%
  leaflet::addTiles("http://{s}.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}.png") %>%
  leaflet::addCircleMarkers(popup = paste( "<b>Species: </b>",
                                           sp2plot$scientificName, 
                                           "<br/>",
                                           "<b>RecordedBy: </b>",
                                           sp2plot$recordedBy,
                                           "<br/>",
                                           "<b>Year: </b>",
                                           sp2plot$year, 
                                           "<br/>",
                                           "<b>EventDate: </b>",
                                           sp2plot$eventDate),
                            lat = sp2plot$decimalLatitude, 
                            lng = sp2plot$decimalLongitude, 
                            radius = 3.5, 
                            weight = 0, 
                            fillColor = unname(sp2plot$color), 
                            fillOpacity = 1, label = labels) 
    }
  }
}




