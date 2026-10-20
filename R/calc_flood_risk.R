#' Calculate Flood Sensitivity Indicator
#' 
#' Calculates flood sensitivity using current FEMA flood zone maps. Returns total area and % cover of
#' high risk flood zones (those with A or V categories) within water supply buffer. This indicates
#' a 1% annual chance of flooding.
#' 
#' Pulls from FEMA WMS endpoint at: https://hazards.fema.gov/arcgis/rest/services/public/NFHL/MapServer
#' 
#' @param data sf object with water supply locations
#' @param dist distance (in meters) of point buffer
#' @param save logical, whether to save output
#' @param output_dir character, directory path for saving files
#' @param verbose logical, whether to show progress messages
#' 
#' @return tibble with water supply ID, park info, and flood sensitivity indicators. Returns 'NA' where FEMA data is not
#' available, and 0 when it is available but there are no high risk flood zones in the buffer.
#' 
#' @export
calc_flood_risk <-
  function(data,
           dist = 0,
           save = FALSE,
           output_dir = NULL,
           verbose = TRUE) {
    
    # Input validation
    stopifnot(
      "water_supplies must be an sf object" = inherits(data, "sf"),
      "save is TRUE but output_dir not provided" = !save ||
        !is.null(output_dir)
    )
    
    # buffer points to specified distance
    data <- data %>% 
      st_transform(crs = 32616) %>% 
      st_buffer(dist = dist)
    
    # match CRS to that of FEMA flood layer
    if (st_crs(data) != st_crs(4269)) {
      data <- st_transform(data, crs = 4269)
    }
    
    # Start timing
    start_time <- Sys.time()
    
    # Create progress bar
    if (verbose) {
      pb <- progress::progress_bar$new(
        format = "Processing flood sensitivity [:bar] :percent (:current/:total) ETA: :eta",
        total = nrow(data),
        clear = FALSE,
        width = 80
      )
    }
    
    # for each water supply....
    results <-  map_dfr(seq_len(nrow(data)), function(i) {
      if (verbose) pb$tick()
    
      
      boundary <- data[i, ]
      
      ## get bounding box
      bb <- st_bbox(boundary)
      
      ## extract bbox numbers
      bb.ordered <- paste(bb[1], bb[2], bb[3], bb[4], sep = "%2C")
      
      ## construct URL (OG version)
      url <- paste0(
        "https://hazards.fema.gov/arcgis/rest/services/public/NFHL/MapServer/",
        28,
        "/query?",
        "&geometry=",
        bb.ordered,
        "&geometryType=esriGeometryEnvelope",
        "&outFields=*",
        "&returnGeometry=true",
        "&returnZ=false",
        "&returnM=false",
        "&returnExtentOnly=false",
        "&f=geoJSON"
      )
      
      ## Claude version
      # url <- paste0(
      #   "https://hazards.fema.gov/arcgis/rest/services/public/NFHL/MapServer/28/query?",
      #   "geometry=", bb[1], ",", bb[2], ",", bb[3], ",", bb[4],
      #   "&geometryType=esriGeometryEnvelope",
      #   "&inSR=4326",  # Input spatial reference
      #   "&outSR=4326", # Output spatial reference
      #   "&spatialRel=esriSpatialRelIntersects",
      #   "&outFields=*",
      #   "&returnGeometry=true",
      #   "&f=geojson"
      # )
      
      ## read in floodplain
      # floodHaz <- sf::read_sf(url) %>%
      #   st_make_valid()
      
      ## read in floodplain with error handling
      floodHaz <- tryCatch(
        {
          sf::read_sf(url) %>%
            st_make_valid()
        },
        error = function(e) {
          if (verbose) {
            message("\nWarning: Unable to read flood data for ",
                    data[i, ]$wsd_source_id)
          }
          # Return NULL to indicate error
          return(NULL)
        }
      )
      
      # if error occurred during read
      if (is.null(floodHaz)) {
        df <- tibble(
          wsd_source_id = data[i, ]$wsd_source_id,
          park_unit = data[i, ]$park_unit,
          park_name = data[i, ]$park_name,
          raw_flood_risk_area_m2 = NA,
          raw_flood_risk_percent = NA
        )
        
      } else if (nrow(floodHaz) == 0) {
        df <- tibble(
          wsd_source_id = data[i, ]$wsd_source_id,
          park_unit = data[i, ]$park_unit,
          park_name = data[i, ]$park_name,
          raw_flood_risk_area_m2 = NA,
          raw_flood_risk_percent = NA
        )
        
      } else {
        ## filter to zones that hve A or V in them - high-risk flood zones/1% flood prob https://floodpartners.com/flood-zones/
        floodRisk <- floodHaz %>%
          filter(stringr::str_detect(FLD_ZONE, "A|V") &
                   FLD_ZONE != "AREA NOT INCLUDED") # not sure if 'area not included' is in fema raw data
        
        
        # if no high risk flood zones
        if (nrow(floodRisk) == 0) {
          df <- tibble(
            wsd_source_id = data[i, ]$wsd_source_id,
            park_unit = data[i, ]$park_unit,
            park_name = data[i, ]$park_name,
            raw_flood_risk_area_m2 = 0,
            raw_flood_risk_percent = 0
          )
        } else {
          # match CRS
          boundary <- st_transform(boundary, st_crs(floodRisk))
          
          # Set attribute-geometry relationship to suppress warning
          st_agr(floodRisk) <- "constant"
          st_agr(boundary) <- "constant"
          
          # return overlapping areas
          floodArea <- st_intersection(floodRisk, boundary)
          
          # if no intersection
          if (nrow(floodArea) == 0) {
            df <- tibble(
              wsd_source_id = data[i, ]$wsd_source_id,
              park_unit = data[i, ]$park_unit,
              park_name = data[i, ]$park_name,
              raw_flood_risk_area_m2 = 0,
              raw_flood_risk_percent = 0
            )
          } else {
            df <- tibble(
              wsd_source_id = data[i, ]$wsd_source_id,
              park_unit = data[i, ]$park_unit,
              park_name = data[i, ]$park_name,
              raw_flood_risk_area_m2 = as.numeric(sum(st_area(
                floodArea
              ))),
              raw_flood_risk_percent = as.numeric(sum(st_area(
                floodArea
              )) / st_area(boundary) * 100)
            )
            
          }
          return(df)
        }
        
      }
    }
    )
    
    # Average across geometries if a wsd_source_id had multiple entries in water_supplies
    results <- results %>% 
      summarise(
        raw_flood_risk_area_m2  = mean(raw_flood_risk_area_m2, na.rm = TRUE),
        raw_flood_risk_percent  = mean(raw_flood_risk_percent, na.rm = TRUE),
        park_unit = first(park_unit),
        park_name = first(park_name),
        .by = wsd_source_id
      )
    
    # Calculate standardized indicators using helper function
    final_results <- results %>%
      add_standardized_indicators("raw_flood_risk")
    
    # Save results
    if (save) {
      if (!dir.exists(output_dir))
        dir.create(output_dir, recursive = TRUE)
      
      filename <- file.path(output_dir,
                            paste0("flood_sensitivity_", Sys.Date(), ".csv"))
      write_csv(final_results, filename)
      
      if (verbose) {
        message("Results saved to: ", filename)
      }
    }
    
    # Print timing
    if (verbose) {
      end_time <- Sys.time()
      total_time <- end_time - start_time
      message(
        "\nFlood sensitivity processing completed in ",
        round(total_time, 2),
        " ",
        attr(total_time, "units")
      )
    }
    
    return(final_results)
    
  }

