#' @title Fire sensitivity: current wildfire hazard potential
#'
#' @description
#' Summarizes the Wildfire Hazard Potential (WHP, Dillon et al. 2023) inside
#' each water supply area. Use buffered water supplies (the rapid workflow
#' uses 4 km buffers): the value is computed over each polygon. If a
#' `wsd_source_id` has several polygons, their values are averaged.
#'
#' @param water_supplies sf object. The water supply areas, as polygons with
#'   `wsd_source_id`, `park_unit` and `park_name` columns, e.g. points from
#'   [wavacore::get_wsd()] buffered with [sf::st_buffer()]. Required.
#' @param whp_raster String. Path to the WHP raster: the continuous raster
#'   for `whp_type = "raw"`, or the classified one (e.g.
#'   `all/usfs_fire/whp2023_GeoTIF/whp2023_cls_conus.tif` in the shared data
#'   folder) for `whp_type = "cat"`. Required.
#' @param whp_type String. Which raster `whp_raster` is, and so what is
#'   calculated. Must be one of:
#'   * `"raw"`: the continuous raster. The result is the area-weighted mean
#'     WHP inside each polygon.
#'   * `"cat"`: the classified raster. The result is the percentage of each
#'     polygon's area in the High or Very High classes (4 and 5).
#'
#'   Either way the result is in `raw_wildfire_hazard_mean`. Required.
#' @param save Logical. Whether to write the results to a CSV file.
#'   * `TRUE`: writes `fire_sensitivity_<date>.csv`, named with today's date,
#'     to `output_dir`, creating the folder if needed. Stops with an error if
#'     `output_dir` is `NULL`.
#'   * `FALSE`: writes nothing.
#'
#'   Default: `FALSE`.
#' @param output_dir String or `NULL`. The folder for the CSV file.
#'   * A folder path: where the file is written when `save = TRUE`. Ignored
#'     when `save = FALSE`.
#'   * `NULL`: only allowed with `save = FALSE`.
#'
#'   Default: `NULL`.
#' @param verbose Logical. Whether to print progress.
#'   * `TRUE`: shows a progress bar over the water supplies and messages with
#'     the run time and the path of the saved file.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return A tibble, one row per water supply, with `wsd_source_id`,
#'   `park_unit`, `park_name`, `raw_wildfire_hazard_mean` (mean WHP, or % of
#'   area in High or Very High, see `whp_type`), and its `norm_` and `rank_`
#'   versions from [add_standardized_indicators()].
#' @references
#' Dillon, G.K. et al. (2023). Wildfire Hazard Potential for the United
#'   States, version 2023. USDA Forest Service Research Data Archive.
#' @export
#' @seealso [calc_fire_exposure()] for the projected change in fire
#'   probability; [calc_vulnerability_index()].
#' @examplesIf file.exists(wavacore::wava_data_path("all", "usfs_fire", "whp2023_GeoTIF", "whp2023_cls_conus.tif"))
#' # 4 km buffers around the water supplies, as in the rapid workflow
#' ws <- wavacore::get_wsd(park = "ZION")
#' ws_buffer <- sf::st_buffer(sf::st_transform(ws, 32616), dist = 4000)
#' fire_sen <- calc_fire_sensitivity(
#'   ws_buffer,
#'   whp_raster = wavacore::wava_data_path("all", "usfs_fire", "whp2023_GeoTIF",
#'                                         "whp2023_cls_conus.tif"),
#'   whp_type = "cat"
#' )
#' fire_sen

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
    binary <- terra::ifel(terra::`%in%`(whp_raster, 4:5), 1, 0)
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
  results <- purrr::map_dfr(seq_len(nrow(water_supplies)), function(i) {
    if (verbose) pb$tick()
    
    # RAW WHP
    
    if(whp_type == "raw") {
      # Calculate area-weighted mean WHP for water supply buffer
      ws_whp <- terra::zonal(
        whp_raster,
        terra::vect(sf::st_transform(water_supplies[i, ], terra::crs(whp_raster))),
        fun = "mean",
        weights = TRUE,
        na.rm = TRUE
      )
      
      # Create results tibble
      dplyr::tibble(
        wsd_source_id = water_supplies[i, ]$wsd_source_id,
        park_unit = water_supplies[i, ]$park_unit,
        park_name = water_supplies[i, ]$park_name,
        raw_wildfire_hazard_mean = dplyr::pull(ws_whp)
      )
      
    } else {
      # Calculate area-weighted % high/very high WHP category for water supply buffer
      
      ## pull out park shapefile and transform
      ws <- terra::vect(sf::st_transform(water_supplies[i, ], terra::crs(whp_raster)))
      
      # Extract values with cell overlap fractions
      extracted <- terra::extract(binary, ws, weights = TRUE, exact = TRUE)
      
      # weights column = fraction of cell in polygon
      total_weight <- sum(extracted$weight, na.rm = TRUE)
      coverage_weight <- sum(extracted$weight[extracted$class_desc == 1], na.rm = TRUE)
      
      percent <- (coverage_weight / total_weight) * 100
      
      # Create results tibble
      dplyr::tibble(
        wsd_source_id = water_supplies[i, ]$wsd_source_id,
        park_unit = water_supplies[i, ]$park_unit,
        park_name = water_supplies[i, ]$park_name,
        raw_wildfire_hazard_mean = percent
      )
      
      
    }
   
    
   
  })
  
  # Average across geometries if a wsd_source_id had multiple entries in water_supplies
  results <- results %>% 
    dplyr::summarise(
      raw_wildfire_hazard_mean = mean(raw_wildfire_hazard_mean, na.rm = TRUE),
      park_unit = dplyr::first(park_unit),
      park_name = dplyr::first(park_name),
      .by = wsd_source_id
    )
  
  # Calculate standardized indicators using helper function  
  final_results <- results %>%
    add_standardized_indicators("raw_wildfire_hazard")
  
  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(output_dir, paste0("fire_sensitivity_", Sys.Date(), ".csv"))
    readr::write_csv(final_results, filename)
    
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
