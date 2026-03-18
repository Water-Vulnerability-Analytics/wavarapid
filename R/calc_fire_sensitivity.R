#' Calculate Fire Sensitivity Indicator
#' 
#' Calculates fire sensitivity using Wildfire Hazard Potential (WHP) data
#' from Dillon et al. 2023. Computes area-weighted mean WHP within water
#' supply areas.
#' 
#' @param water_supplies sf object with water supply locations
#' @param whp_raster file path to wildfire hazard potential raster
#' @param whp_type either "raw" for the raw continuous raster or "cat" for the categorized raster
#' @param save logical, whether to save output
#' @param output_dir character, directory path for saving files
#' @param verbose logical, whether to show progress messages
#' 
#' @return tibble with water supply ID, park info, and fire sensitivity indicators
#' 
calc_fire_sensitivity <- function(water_supplies,
                                  whp_raster,
                                  whp_type,
                                  save = FALSE,
                                  output_dir = NULL,
                                  verbose = TRUE) {
  
  # Input validation
  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    #"whp_raster must be a terra raster" = inherits(whp_raster, "SpatRaster"),
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir)
  )
  
  whp_raster <- terra::rast(whp_raster)
  
  # if using the categorical layer, classify here to use binary layer in iteration below
  if(whp_type == "cat") {
    ## classify raster
    binary <- terra::ifel(whp_raster %in% 4:5, 1, 0)
  }
  
  
  # Start timing
  start_time <- Sys.time()
  
  # Create progress bar
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "Processing fire sensitivity [:bar] :percent (:current/:total) ETA: :eta",
      total = nrow(water_supplies),
      clear = FALSE,
      width = 80
    )
  }
  
  # Process each water supply
  results <- map_dfr(seq_len(nrow(water_supplies)), function(i) {
    if (verbose) pb$tick()
    
    # RAW WHP
    
    if(whp_type == "raw") {
      # Calculate area-weighted mean WHP for water supply buffer
      ws_whp <- terra::zonal(
        whp_raster,
        terra::vect(st_transform(water_supplies[i, ], crs(whp_raster))),
        fun = "mean",
        weights = TRUE,
        na.rm = TRUE
      )
      
      # Create results tibble
      tibble(
        wsd_source_id = water_supplies[i, ]$wsd_source_id,
        park_unit = water_supplies[i, ]$park_unit,
        park_name = water_supplies[i, ]$park_name,
        raw_wildfire_hazard_mean = pull(ws_whp)
      )
      
    } else {
      # Calculate area-weighted % high/very high WHP category for water supply buffer
      
      ## pull out park shapefile and transform
      ws <- terra::vect(st_transform(water_supplies[i, ], crs(whp_raster)))
      
      # Extract values with cell overlap fractions
      extracted <- terra::extract(binary, ws, weights = TRUE, exact = TRUE)
      
      # weights column = fraction of cell in polygon
      total_weight <- sum(extracted$weight, na.rm = TRUE)
      coverage_weight <- sum(extracted$weight[extracted$class_desc == 1], na.rm = TRUE)
      
      percent <- (coverage_weight / total_weight) * 100
      
      # Create results tibble
      tibble(
        wsd_source_id = water_supplies[i, ]$wsd_source_id,
        park_unit = water_supplies[i, ]$park_unit,
        park_name = water_supplies[i, ]$park_name,
        raw_wildfire_hazard_mean = percent
      )
      
      
    }
   
    
   
  })
  
  # Average across geometries if a wsd_source_id had multiple entries in water_supplies
  results <- results |>
    summarise(
      raw_wildfire_hazard_mean = mean(raw_wildfire_hazard_mean, na.rm = TRUE),
      park_unit = first(park_unit),
      park_name = first(park_name),
      .by = wsd_source_id
    )
  
  # Calculate standardized indicators using helper function  
  final_results <- results %>%
    add_standardized_indicators("raw_wildfire_hazard")
  
  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(output_dir, paste0("fire_sensitivity_", Sys.Date(), ".csv"))
    write_csv(final_results, filename)
    
    if (verbose) {
      message("Results saved to: ", filename)
    }
  }
  
  # Print timing
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message("\nFire sensitivity processing completed in ", 
            round(total_time, 2), " ", attr(total_time, "units"))
  }
  
  return(final_results)
}