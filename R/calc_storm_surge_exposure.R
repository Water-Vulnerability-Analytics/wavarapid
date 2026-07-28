#' Calculate Storm Surge Exposure Indicator
#'
#' Assigns each NPS water supply the worst-case model 100-year return period
#' storm surge change (future 2021-2050 minus historical 1985-2014, in meters)
#' from the nearest coastal station in the Muis et al. (2023) CMIP6 dataset.
#'
#' Water supplies whose nearest surge station is more than 50 km away receive
#' NA because they are too far from a measurement point to receive a meaningful
#' surge value. Water supplies who are more than 50 km away from the CUSP coastline 
#' also receiv NA.
#'
#' weighting: Coastal proximity and low elevation are multiplied. A supply must
#'    therefore be both close to the coast and at low elevation to receive a
#'    high weight.
#'
#' Source: Muis et al. (2023), "Global Projections of Storm Surges Using
#' High-Resolution CMIP6 Climate Models", Earth's Future.
#'
#' @param water_supplies sf object with water supply point locations
#' @param storm_surge_stations character path to the preprocessed CONUS storm
#'   surge GeoPackage produced by `process_storm_surge_data()`, or an sf
#'   object already loaded into memory
#' @param coast_distances character path to coast_distances.csv produced by
#'   `process_coastal_weights()`
#' @param site_elevations character path to the intermediary elevation CSV
#'   produced by `process_elevation_data()`
#' @param save logical, whether to save output to a CSV file
#' @param output_dir character, directory path for saving files
#' @param verbose logical, whether to show progress messages
#'
#' @return tibble with raw and standardized storm surge exposure indicators
#'   using RMS and joint distance-elevation weighting
#'
calc_storm_surge_exposure <- function(water_supplies,
                                      storm_surge_stations,
                                      coast_distances = here("data",
                                                             "rapid",
                                                             "inputs",
                                                             "coast_comp",
                                                             "coast_distances.csv"),
                                      site_elevations = here("data",
                                                             "rapid",
                                                             "inputs",
                                                             "water_supply_elevations.csv"),
                                      save = FALSE,
                                      output_dir = NULL,
                                      verbose = TRUE) {
  
  # Input validation
  stopifnot(
    "water_supplies must be an sf object" =
      inherits(water_supplies, "sf"),
    "storm_surge_stations must be a file path or sf object" =
      is.character(storm_surge_stations) ||
      inherits(storm_surge_stations, "sf"),
    "save is TRUE but output_dir not provided" =
      !save || !is.null(output_dir),
    "coast distances file does not exist" =
      file.exists(coast_distances),
    "site elevations file does not exist — run process_elevation_data() first" =
      file.exists(site_elevations)
  )
  
  if (verbose) {
    start_time <- Sys.time()
  }
  
  # Load storm surge stations if a file path was supplied
  if (is.character(storm_surge_stations)) {
    stopifnot(
      "storm_surge_stations file does not exist — run process_storm_surge_data() first" =
        file.exists(storm_surge_stations)
    )
    
    if (verbose) {
      message(
        "Reading storm surge stations from: ",
        storm_surge_stations
      )
    }
    
    storm_surge_stations <- sf::st_read(storm_surge_stations,
                                        quiet = TRUE)
  }
  
  # Project both layers to NAD83 / CONUS Albers for distance calculations
  if (verbose) {
    message("Projecting to CONUS Albers (EPSG:5070)...")
  }
  
  ws_proj <- sf::st_transform(water_supplies,
                              5070)
  
  stations_proj <- sf::st_transform(storm_surge_stations,
                                    5070)
  
  # Find the nearest storm surge station
  if (verbose) {
    message("Computing nearest surge station for ",
            nrow(ws_proj),
            " water supplies...")
  }
  
  nn_result <- nngeo::st_nn(ws_proj,
                            stations_proj,
                            k = 1,
                            returnDist = TRUE,
                            progress = verbose)
  
  nn_idx <- unlist(nn_result$nn)
  nn_dist <- unlist(nn_result$dist)
  
  # Assign nearest-station storm surge values
  results <- water_supplies %>%
    sf::st_drop_geometry() %>%
    dplyr::select(wsd_source_id,
                  park_unit,
                  park_name) %>%
    dplyr::mutate(nn_dist_m = nn_dist,
                  raw_surge_change_max =
                    stations_proj$max_surge_change_m[nn_idx],
                  raw_surge_pct_change_max =
                    stations_proj$max_surge_pct_change[nn_idx],
                  raw_surge_pct_models_increasing =
                    stations_proj$pct_models_increasing[nn_idx]) %>%
    # Supplies >50 km from the nearest surge station receive NA
    dplyr::mutate(raw_surge_change_max = dplyr::if_else(nn_dist_m > 50000,
                                                        NA_real_,
                                                        raw_surge_change_max),
                  raw_surge_pct_change_max = dplyr::if_else(nn_dist_m > 50000,
                                                            NA_real_,
                                                            raw_surge_pct_change_max),
                  raw_surge_pct_models_increasing = dplyr::if_else(nn_dist_m > 50000,
                                                                   NA_real_,
                                                                   raw_surge_pct_models_increasing)) %>%
    dplyr::select(-nn_dist_m)
  
  # Read cached elevation and coastline-distance data
  if (verbose) {
    message("Reading cached elevation and coastline-distance data...")
  }
  
  site_elevations_df <- readr::read_csv(site_elevations,
                                        show_col_types = FALSE)
  
  coast_distances_df <- readr::read_csv(coast_distances,
                                        show_col_types = FALSE)
  
  results <- results %>%
    dplyr::left_join(site_elevations_df,
                     by = c("wsd_source_id",
                            "park_unit",
                            "park_name")) %>%
    dplyr::left_join(coast_distances_df,
                     by = c("wsd_source_id",
                            "park_unit",
                            "park_name")) %>%
    mutate(raw_surge_change_max = dplyr::if_else(coast_dist_m > 50000,
                                                 NA_real_,
                                                 raw_surge_change_max))
  
  # Maximum elevation used to normalize site elevation
  max_elevation_m <- results %>%
    dplyr::filter(!is.na(raw_surge_change_max)) %>%
    dplyr::summarize(max_elevation_m = max(pmax(site_elevation_m, 0), na.rm = TRUE)) %>%
    dplyr::pull(max_elevation_m)
  
  # Protect against all elevations being missing or zero
  if (!is.finite(max_elevation_m) || max_elevation_m <= 0) {
    max_elevation_m <- 1
  }
  
  results <- results %>%
    # Coastal-proximity weight:
    # 1 at the shoreline and 0 at or beyond 50 km
    dplyr::mutate(coast_weight = pmax(0, 1 - pmin(coast_dist_m, 50000) / 50000),
                  
                  # Low-elevation weight:
                  # 1 at sea level and 0 at the maximum observed site elevation
                  elev_weight = 1 - pmin(pmax(site_elevation_m, 0) / max_elevation_m, 1),
                  
                  # # RMS / normalized Euclidean approach:
                  # # one factor can partially compensate for the other
                  # rms_weight = sqrt(
                  #   (coast_weight^2 + elev_weight^2) / 2
                  # ),
                  # 
                  # raw_rms_wght_surge_change_max =
                  #   raw_surge_change_max * rms_weight,
                  # 
                  # raw_rms_wght_surge_pct_change_max =
                  #   raw_surge_pct_change_max * rms_weight,
                  
                  # Product approach:
                  # both proximity and low elevation must be high
                  joint_weight = coast_weight * elev_weight,
                  raw_surge_change_max = raw_surge_change_max * joint_weight,
                  raw_surge_pct_change_max = raw_surge_pct_change_max * joint_weight) %>%
    dplyr::select(-c(coast_dist_m,
                     coast_weight,
                     elev_weight,
                     joint_weight))
  
  # Standardize all indicators to 0-1
  final_results <- results %>%
    add_standardized_indicators("raw_surge_change_max", invert = FALSE) %>%
    add_standardized_indicators("raw_surge_pct_change_max", invert = FALSE)
  
  # Save results
  if (save) {
    if (!dir.exists(output_dir)) {
      dir.create(
        output_dir,
        recursive = TRUE
      )
    }
    
    filename <- file.path(
      output_dir,
      paste0(
        "storm_surge_exposure_",
        Sys.Date(),
        ".csv"
      )
    )
    
    readr::write_csv(
      final_results,
      filename
    )
    
    if (verbose) {
      message("Results saved to: ", filename)
    }
  }
  
  # Report timing
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    
    message(
      "\nStorm surge exposure processing completed in ",
      round(total_time, 2),
      " ",
      attr(total_time, "units")
    )
  }
  
  return(final_results)
}