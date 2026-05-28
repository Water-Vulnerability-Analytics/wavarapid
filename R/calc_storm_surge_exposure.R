#' Calculate Storm Surge Exposure Indicator
#'
#' Assigns each NPS water supply the worst-case model 100-year return period
#' storm surge change (future 2021-2050 minus historical 1985-2014, in meters)
#' from the nearest coastal station in the Muis et al. (2023) CMIP6 dataset.
#' Water supplies whose nearest station is more than 50 km away receive NA (too
#' far from any surge measurement point to be considered in this analysis).
#' Assigned values are then divided by `shore_dist_m` from
#' `process_coastal_weights()` to weight exposure by coastal proximity.
#'
#' Source: Muis et al. (2023), "Global Projections of Storm Surges Using
#' High-Resolution CMIP6 Climate Models", Earth's Future.
#' Note: SSP5-8.5 only — no SSP2-4.5 variant exists in this dataset.
#'
#' @param water_supplies sf object with water supply point locations
#' @param storm_surge_stations character path to the preprocessed CONUS storm
#'   surge GeoPackage produced by `process_storm_surge_data()`, or an sf
#'   object already loaded into memory
#' @param shore_distances character path to shore_distances.csv produced by
#'   `process_coastal_weights()`. Defaults to the standard rapid inputs path.
#' @param save logical, whether to save output to a CSV file
#' @param output_dir character, directory path for saving files
#' @param verbose logical, whether to show progress messages
#'
#' @return tibble with `wsd_source_id`, `park_unit`, `park_name`, and
#'   standardized storm surge exposure columns (`raw_surge_change_max`,
#'   `norm_surge_change_max`, `rank_surge_change_max`). Water supplies >50 km
#'   from any surge station have `NA` for all three columns.
#'
calc_storm_surge_exposure <- function(water_supplies,
                                      storm_surge_stations,
                                      shore_distances = here("data", "rapid", "inputs", "coast_comp", "shore_distances.csv"),
                                      save = FALSE,
                                      output_dir = NULL,
                                      verbose = TRUE) {

  # Input validation
  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "storm_surge_stations must be a file path or sf object" =
      is.character(storm_surge_stations) || inherits(storm_surge_stations, "sf"),
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir),
    "shore distances file does not exist" = file.exists(shore_distances)
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
    k          = 1,
    returnDist = TRUE,
    progress   = verbose
  )

  nn_idx  <- unlist(nn_result$nn)
  nn_dist <- unlist(nn_result$dist)   # meters

  # Build results table
  results <- water_supplies %>%
    sf::st_drop_geometry() %>%
    select(wsd_source_id, park_unit, park_name) %>%
    mutate(
      nn_dist_m          = nn_dist,
      raw_surge_change_max = stations_proj$max_surge_change_m[nn_idx]
    ) %>%
    # Supplies > 50 km from the nearest surge station receive NA
    mutate(
      raw_surge_change_max = if_else(nn_dist_m > 50000, NA_real_, raw_surge_change_max)
    ) %>%
    select(-nn_dist_m)

  # Divide by distance to coast to weight exposure by coastal proximity
  shore_distances_df <- read_csv(shore_distances, show_col_types = FALSE)

  results <- results %>%
    left_join(shore_distances_df, by = c("wsd_source_id", "park_unit", "park_name")) %>%
    mutate(raw_surge_change_max = raw_surge_change_max / shore_dist_m) %>%
    select(-shore_dist_m)

  # Standardize: invert = FALSE — higher surge change = higher exposure
  final_results <- results %>%
    add_standardized_indicators("raw_surge_change_max", invert = FALSE)

  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    filename <- file.path(output_dir, paste0("storm_surge_exposure_", Sys.Date(), ".csv"))
    write_csv(final_results, filename)

    if (verbose) message("Results saved to: ", filename)
  }

  # Report timing
  if (verbose) {
    end_time   <- Sys.time()
    total_time <- end_time - start_time
    message("\nStorm surge exposure processing completed in ",
            round(total_time, 2), " ", attr(total_time, "units"))
  }

  return(final_results)
}
