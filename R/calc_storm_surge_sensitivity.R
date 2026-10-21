#' Calculate Storm Surge Sensitivity Indicator
#'
#' Assigns each NPS water supply the 4-model mean 1985-2014 reference-period
#' 10-year return period storm surge level (meters) from the nearest coastal
#' station in the Muis et al. (2023) CMIP6 dataset, discounted by a joint
#' coastal-proximity / low-elevation weight.
#'
#' - Mean across all 4 models; no model excluded.
#' - Water supplies whose nearest surge station is >50 km away, or that are
#'   more than 50 km from the CUSP coastline, receive NA.
#' - Distance/elevation weight is the same RMS blend as
#'   `calc_storm_surge_exposure()`: `sqrt((dist_weight^2 + elev_weight^2) / 2)`,
#'   with `elev_weight` capped at the 10 m Low Elevation Coastal Zone
#'   threshold (McGranahan et al. 2007).
#' - Muis-only; a SLOSH-based second source is not incorporated.
#'
#' See `Explore/explore_updated_storm_surge.Rmd` for the reasoning behind the
#' mean aggregation, the weighting, and the Muis-only (SLOSH-deferred) scope.
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
#' @return tibble with raw and standardized storm surge sensitivity indicators
#'
#' @export
calc_storm_surge_sensitivity <- function(water_supplies,
                                         storm_surge_stations,
                                         coast_geometry = wavacore::wava_data_path("rapid",
                                                               "inputs",
                                                               "coast_comp",
                                                               "coastal_site_geometry.csv"),
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
    "storm_surge_stations file does not exist \u2014 run process_storm_surge_data() first" =
      !is.character(storm_surge_stations) || file.exists(storm_surge_stations),
    "coast_geometry must be a file path or data frame" =
      is.character(coast_geometry) || is.data.frame(coast_geometry),
    "coast geometry file does not exist \u2014 run process_coastal_geometry() first" =
      !is.character(coast_geometry) || file.exists(coast_geometry),
    "save is TRUE but output_dir not provided" =
      !save || !is.null(output_dir)
  )

  if (verbose) {
    start_time <- Sys.time()
  }

  # Load storm surge stations if a file path was supplied
  if (is.character(storm_surge_stations)) {
    if (verbose) {
      message(
        "Reading storm surge stations from: ",
        storm_surge_stations
      )
    }

    storm_surge_stations <- sf::st_read(storm_surge_stations, quiet = TRUE)
  }

  # Reduce to one row per station: mean reference_mean across all 4 models
  prepped_stations <- storm_surge_stations %>%
    sf::st_drop_geometry() %>%
    dplyr::group_by(station_idx) %>%
    dplyr::summarize(reference_mean = mean(reference_mean, na.rm = TRUE),
             lat = dplyr::first(lat),
             lon = dplyr::first(lon),
             .groups = "drop") %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

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

  # Assign nearest-station reference-period storm surge values
  nn_result <- nngeo::st_nn(ws_coastal, stations_proj, k = 1, returnDist = TRUE, progress = FALSE)
  nn_idx <- unlist(nn_result$nn)
  nn_dist_m <- unlist(nn_result$dist)

  station_matches <- ws_coastal %>%
    sf::st_drop_geometry() %>%
    dplyr::select(wsd_source_id, park_unit, park_name, coast_dist_m, coast_elev_m) %>%
    dplyr::mutate(
      raw_surge_reference_m = stations_proj$reference_mean[nn_idx],
      nn_dist_m = nn_dist_m
    ) %>%
    # Supplies >50 km from the nearest surge station receive NA
    dplyr::mutate(raw_surge_reference_m = dplyr::if_else(nn_dist_m > 50000,
                                                         NA_real_,
                                                         raw_surge_reference_m)) %>%
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
      raw_surge_reference_m = raw_surge_reference_m * joint_weight) %>%
    dplyr::select(-c(coast_dist_m, dist_weight, elev_weight, joint_weight))

  if (verbose) {message("Standardizing indicators to 0-1...")}

  # Standardize indicator
  final_results <- results %>%
    add_standardized_indicators("raw_surge_reference_m", invert = FALSE)

  # Save results
  if (save) {
    if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = TRUE)}

    filename <- file.path(output_dir,
                          paste0("storm_surge_sensitivity_", Sys.Date(), ".csv"))

    readr::write_csv(final_results, filename)

    if (verbose) {message("Results saved to: ", filename)}

  }

  # Report timing
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time

    message(
      "\nStorm surge sensitivity processing completed in ",
      round(total_time, 2),
      " ",
      attr(total_time, "units")
    )
  }

  return(final_results)
}
