#' Calculate Storm Surge Exposure Indicator
#'
#' Assigns each NPS water supply the worst-case-model 10-year return period
#' storm surge change (future 2021-2050 minus reference 1985-2014, meters) from
#' the nearest coastal station in the Muis et al. (2023) CMIP6 dataset,
#' discounted by a joint coastal-proximity / low-elevation weight.
#'
#' - CMCC-CM2-VHR4 is excluded from the worst-case model pool.
#' - Water supplies whose nearest surge station is >50 km away, or that are
#'   more than 50 km from the CUSP coastline, receive NA.
#' - Distance/elevation weight: `sqrt((dist_weight^2 + elev_weight^2) / 2)`,
#'   with `elev_weight` capped at the 10 m Low Elevation Coastal Zone
#'   threshold (McGranahan et al. 2007).
#'
#' See `Explore/explore_updated_storm_surge.Rmd` for the reasoning behind the
#' worst-case selection, the CMCC-CM2-VHR4 exclusion, and the weighting scheme.
#'
#' Source: Muis et al. (2023), "Global Projections of Storm Surges Using
#' High-Resolution CMIP6 Climate Models", Earth's Future.
#'
#' @param water_supplies sf object with water supply point locations
#' @param storm_surge_stations character path to the preprocessed storm
#'   surge GeoPackage produced by `process_storm_surge_data()`, or an sf
#'   object already loaded into memory
#' @param coast_geometry character path to coastal_site_geometry.csv (columns
#'   include `coast_dist_m` and `coast_elev_m`) produced by
#'   `process_coastal_geometry()`
#' @param save logical, whether to save output to a CSV file
#' @param output_dir character, directory path for saving files
#' @param verbose logical, whether to show progress messages
#'
#' @return tibble with raw and standardized storm surge exposure indicators
#'
#' @export
calc_storm_surge_exposure <- function(water_supplies,
                                      storm_surge_stations,
                                      coast_geometry = wava_data_path("rapid",
                                                            "inputs", "coast_comp",
                                                            "coastal_site_geometry.csv"),
                                      save = FALSE,
                                      output_dir = NULL,
                                      verbose = TRUE) {
  
  # Input validation
  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "storm_surge_stations must be a file path or sf object" = is.character(storm_surge_stations) || inherits(storm_surge_stations, "sf"),
    "storm_surge_stations file does not exist \u2014 run process_storm_surge_data() first" =
      !is.character(storm_surge_stations) || file.exists(storm_surge_stations),
    "coast_geometry must be a file path or data frame" = is.character(coast_geometry) || is.data.frame(coast_geometry),
    "coast geometry file does not exist \u2014 run process_coastal_geometry() first" =
      !is.character(coast_geometry) || file.exists(coast_geometry),
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir)
  )

  if (verbose) {
    start_time <- Sys.time()
  }

  # Load storm surge stations if a file path was supplied
  if (is.character(storm_surge_stations)) {
    if (verbose) {message( "Reading storm surge stations from: ", storm_surge_stations)}
    storm_surge_stations <- sf::st_read(storm_surge_stations, quiet = TRUE)
  }
  
  # Prepare the storm surge data for indicator calculation
  prepped_stations <- storm_surge_stations %>%
    # Removing CMCC-CM2-VHR4 model since it has spurious results. 
    filter(!grepl("CMCC-CM2-VHR4", model)) %>% 
    group_by(station_idx) %>%
    mutate(raw_surge_change_max = future_mean - reference_mean,
           raw_surge_pct_change_max = ((future_mean - reference_mean) / reference_mean) * 100) %>%
    slice_max(raw_surge_pct_change_max, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(station_idx, model, raw_surge_change_max, raw_surge_pct_change_max, lat, lon)
  
  # Load coastal site geometries if a file path was supplied
  if (is.character(coast_geometry)) {
    if (verbose) {message( "Reading coastal site geometries from: ", coast_geometry)}
    coast_geometry <- readr::read_csv(coast_geometry, show_col_types = FALSE)
  }

  # Project both layers to NAD83 / CONUS Albers for distance calculations
  if (verbose) {message("Projecting to CONUS Albers (EPSG:5070)...")}

  ws_proj <- sf::st_transform(water_supplies, 5070) %>%
    # Join the coastal geometry information
    dplyr::left_join(coast_geometry, by = c("wsd_source_id", "park_unit", "park_name"))

  stations_proj <- sf::st_transform(prepped_stations, 5070)

  # Restrict the nearest-station lookup to sites within 50 km of the coast
  ws_coastal <- dplyr::filter(ws_proj, !is.na(coast_dist_m))

  if (verbose) {
    message("Computing nearest surge station for ", nrow(ws_coastal), " water supplies...")
  }

  # Assign nearest-station storm surge values
  nn_result <- nngeo::st_nn(ws_coastal, stations_proj, k = 1, returnDist = TRUE, progress = FALSE)
  nn_idx <- unlist(nn_result$nn)
  nn_dist_m <- unlist(nn_result$dist)

  station_matches <- ws_coastal %>%
    sf::st_drop_geometry() %>%
    dplyr::select(wsd_source_id, park_unit, park_name, coast_dist_m, coast_elev_m) %>%
    dplyr::mutate(
      raw_surge_change_max = stations_proj$raw_surge_change_max[nn_idx],
      raw_surge_pct_change_max = stations_proj$raw_surge_pct_change_max[nn_idx],
      nn_dist_m = nn_dist_m
    ) %>%
    # Supplies >50 km from the nearest surge station receive NA
    dplyr::mutate(dplyr::across(c(raw_surge_change_max, raw_surge_pct_change_max),
                                ~ dplyr::if_else(nn_dist_m > 50000, NA_real_, .x))) %>%
    dplyr::select(-nn_dist_m)

  if (verbose) {message("Expanding to all water supplies and applying distance/elevation weighting...")}

  # Expand back to all water supplies; non-coastal / >50km sites get NA
  results <- water_supplies %>%
    sf::st_drop_geometry() %>%
    dplyr::select(wsd_source_id, park_unit, park_name) %>%
    dplyr::left_join(station_matches, by = c("wsd_source_id", "park_unit", "park_name")) %>%
    dplyr::mutate(
      # Coastal-proximity weight: 1 at the shoreline and 0 at or beyond 50km
      dist_weight = pmax(0, (1 - (coast_dist_m / 50000))),
      # Low-elevation weight: 1 at sea level and 0 at the 10 m LECZ cap
      elev_weight = 1 - pmin(coast_elev_m / 10, 1),
      # Joint weight: RMS blend to compensate for potential 0s from LECZ cap.
      joint_weight = sqrt((dist_weight^2 + elev_weight^2) / 2),
      raw_surge_change_max = raw_surge_change_max * joint_weight,
      raw_surge_pct_change_max = raw_surge_pct_change_max * joint_weight) %>%
    dplyr::select(-c(coast_dist_m, dist_weight, elev_weight, joint_weight))
  
  if (verbose) {message("Standardizing indicators to 0-1...")}

  # Standardize all indicators to 0-1
  final_results <- results %>%
    add_standardized_indicators("raw_surge_change_max", invert = FALSE) %>%
    add_standardized_indicators("raw_surge_pct_change_max", invert = FALSE)
  
  # Save results
  if (save) {
    
    if (!dir.exists(output_dir)) {dir.create(output_dir,recursive = TRUE)}
    
    filename <- file.path(output_dir, 
                          paste0("storm_surge_exposure_", Sys.Date(), ".csv"))
    
    readr::write_csv(final_results, filename)
    
    if (verbose) {message("Results saved to: ", filename)}
    
  }
  
  # Report timing
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    
    message("\nStorm surge exposure processing completed in ", 
            round(total_time, 2), " ", attr(total_time, "units")
    )
  }
  
  return(final_results)
}