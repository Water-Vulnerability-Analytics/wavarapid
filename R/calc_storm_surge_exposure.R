#' Calculate Storm Surge Exposure Indicator
#'
#' Assigns each NPS water supply the worst-case model 100-year return period
#' storm surge change (future 2021-2050 minus historical 1985-2014, in meters)
#' from the nearest coastal station in the Muis et al. (2023) CMIP6 dataset.
#' Two separate distance checks are applied:
#' 1. Water supplies whose nearest surge *station* is more than 50 km away
#'    receive NA — they are too far from any measurement point to be assigned
#'    a meaningful surge value.
#' 2. For supplies that pass that threshold, the assigned surge value is divided
#'    by `coast_dist_m` from `process_coastal_weights()` (distance to the NOAA
#'    CUSP shoreline) to weight exposure by coastal proximity. These are
#'    different distances: one screens for data coverage, the other scales the
#'    exposure value.
#'
#' Source: Muis et al. (2023), "Global Projections of Storm Surges Using
#' High-Resolution CMIP6 Climate Models", Earth's Future.
#' SSP5-8.5 only — the only scenario available in this dataset. Uses best-fit
#' (central estimate) only; low- and high-bound confidence intervals exist in
#' the source data but are not used. See `process_storm_surge_data()` for
#' full data provenance.
#'
#' @param water_supplies sf object with water supply point locations
#' @param storm_surge_stations character path to the preprocessed CONUS storm
#'   surge GeoPackage produced by `process_storm_surge_data()`, or an sf
#'   object already loaded into memory
#' @param coast_distances character path to coast_distances.csv produced by
#'   `process_coastal_weights()`. Defaults to the standard rapid inputs path.
#' @param save logical, whether to save output to a CSV file
#' @param output_dir character, directory path for saving files
#' @param verbose logical, whether to show progress messages
#'
#' @return tibble with `wsd_source_id`, `park_unit`, `park_name`, and
#'   standardized storm surge exposure columns for both absolute change
#'   (`raw_surge_change_max`, `norm_surge_change_max`, `rank_surge_change_max`)
#'   and percent change (`raw_surge_pct_change_max`, `norm_surge_pct_change_max`,
#'   `rank_surge_pct_change_max`), each weighted by `coast_dist_m`. Water
#'   supplies >50 km from any surge station have `NA` for all columns.
#'
calc_storm_surge_exposure <- function(water_supplies,
                                      storm_surge_stations,
                                      coast_distances = here("data", "rapid", "inputs", "coast_comp", "coast_distances.csv"),
                                      save = FALSE,
                                      output_dir = NULL,
                                      verbose = TRUE) {
  
  # Input validation
  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "storm_surge_stations must be a file path or sf object" =
      is.character(storm_surge_stations) || inherits(storm_surge_stations, "sf"),
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir),
    "coast distances file does not exist" = file.exists(coast_distances)
  )
  
  if (verbose) start_time <- Sys.time()
  
  # Load storm surge stations if a file path was supplied
  if (is.character(storm_surge_stations)) {
    stopifnot(
      "storm_surge_stations file does not exist — run process_storm_surge_data() first" =
        file.exists(storm_surge_stations)
    )
    if (verbose) message("Reading storm surge stations from: ", storm_surge_stations)
    storm_surge_stations <- sf::st_read(storm_surge_stations, quiet = TRUE)
  }
  
  # Project both layers to NAD83 / CONUS Albers (meters) for distance calculations
  if (verbose) message("Projecting to CONUS Albers (EPSG:5070)...")
  ws_proj <- sf::st_transform(water_supplies, 5070)
  stations_proj <- sf::st_transform(storm_surge_stations, 5070)
  
  # Nearest-neighbor join: k = 1, return distances
  if (verbose) message("Computing nearest surge station for ", nrow(ws_proj), " water supplies...")
  
  nn_result <- nngeo::st_nn(
    ws_proj,
    stations_proj,
    k = 1,
    returnDist = TRUE,
    progress = verbose
  )
  
  nn_idx  <- unlist(nn_result$nn)
  nn_dist <- unlist(nn_result$dist) # meters
  
  # Build results table
  results <- water_supplies %>%
    sf::st_drop_geometry() %>%
    select(wsd_source_id, park_unit, park_name) %>%
    mutate(
      nn_dist_m                 = nn_dist,
      raw_surge_change_max      = stations_proj$max_surge_change_m[nn_idx],
      raw_surge_pct_change_max  = stations_proj$max_surge_pct_change[nn_idx],
      raw_surge_model_agreement = stations_proj$model_agreement[nn_idx]
    ) %>%
    # Supplies > 50 km from the nearest surge station receive NA
    mutate(
      raw_surge_change_max     = if_else(nn_dist_m > 50000, NA_real_, raw_surge_change_max),
      raw_surge_pct_change_max = if_else(nn_dist_m > 50000, NA_real_, raw_surge_pct_change_max)
    ) %>%
    select(-nn_dist_m)
  
  # Divide both metrics by distance to coast to weight exposure by coastal proximity
  coast_distances_df <- read_csv(coast_distances, show_col_types = FALSE)
  
  results <- results %>%
    left_join(coast_distances_df, by = c("wsd_source_id", "park_unit", "park_name")) %>%
    # Apply a bounded coastal-proximity weight so surge exposure decreases with
    # distance from the coast without producing extreme values near zero distance.
    mutate(coast_weight = pmax(0, 1 - pmin(coast_dist_m, 50000) / 50000),
           raw_dist_wght_surge_change_max = raw_surge_change_max * coast_weight,
           raw_dist_wght_surge_pct_change_max = raw_surge_pct_change_max * coast_weight,
           old_dist_surge_change_max = raw_surge_change_max / coast_dist_m,
           old_dist_surge_pct_change_max = raw_surge_pct_change_max / coast_dist_m) %>%
    select(-coast_dist_m, -coast_weight)
  
  # Standardize: invert = FALSE — higher surge change = higher exposure
  final_results <- results %>%
    add_standardized_indicators(
      "raw_dist_wght_surge_change_max",
      invert = FALSE
    ) %>%
    add_standardized_indicators(
      "raw_dist_wght_surge_pct_change_max",
      invert = FALSE
    )
  
  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(
      output_dir,
      paste0("storm_surge_exposure_", Sys.Date(), ".csv")
    )
    
    write_csv(final_results, filename)
    
    if (verbose) message("Results saved to: ", filename)
  }
  
  # Report timing
  if (verbose) {
    end_time   <- Sys.time()
    total_time <- end_time - start_time
    
    message(
      "\nStorm surge exposure processing completed in ",
      round(total_time, 2), " ", attr(total_time, "units")
    )
  }
  
  return(final_results)
}