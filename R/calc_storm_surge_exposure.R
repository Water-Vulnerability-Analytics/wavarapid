#' @title Storm surge exposure: projected change in the 10-year storm surge
#'
#' @description
#' Gives each coastal water supply the projected change in the 10-year
#' return period storm surge at its nearest station in the Muis et al. (2023)
#' CMIP6 dataset: future (2021-2050) minus reference (1985-2014), in meters
#' and percent. At each station the worst-case model (the largest percent
#' change) is used; CMCC-CM2-VHR4 is left out because its results are
#' spurious.
#'
#' The change is then multiplied by a weight for how close the supply is to
#' the coast and how low it is:
#' `sqrt((dist_weight^2 + elev_weight^2) / 2)`, where `dist_weight` falls
#' from 1 at the shoreline to 0 at 50 km and `elev_weight` from 1 at sea
#' level to 0 at 10 m (the Low Elevation Coastal Zone, McGranahan et al.
#' 2007).
#'
#' Supplies more than 50 km from the coast, or more than 50 km from the
#' nearest station, get `NA`. See `Explore/explore_updated_storm_surge.Rmd`
#' in `nps_water_vulnerability` for the reasoning behind the worst-case
#' choice, the model exclusion and the weights.
#'
#' Source: Muis et al. (2023), "Global Projections of Storm Surges Using
#' High-Resolution CMIP6 Climate Models", Earth's Future.
#'
#' @param water_supplies sf object. The water supplies, as points with
#'   `wsd_source_id`, `park_unit` and `park_name` columns, e.g. from
#'   [wavacore::get_wsd()]. Required.
#' @param storm_surge_stations String or sf object. The stations from
#'   [process_storm_surge_data()].
#'   * String: path to `storm_surge_RP10_muis_2023.gpkg`, e.g. in
#'     `rapid/inputs/storm_surge` in the shared data folder. Stops with an
#'     error if the file doesn't exist.
#'   * sf object: the stations, already read.
#'
#'   Required.
#' @param coast_geometry String or data frame. The distance to the coast
#'   and elevation of each supply (`coast_dist_m`, `coast_elev_m`), from
#'   [process_coastal_geometry()], joined by `wsd_source_id`, `park_unit`
#'   and `park_name`.
#'   * String: path to `coastal_site_geometry.csv`. Stops with an error if
#'     the file doesn't exist.
#'   * Data frame: the same table, already read.
#'
#'   Default: `rapid/inputs/coast_comp/coastal_site_geometry.csv` in the
#'   shared data folder, see [wavacore::wava_data_path()].
#' @param save Logical. Whether to write the results to a CSV file.
#'   * `TRUE`: writes `storm_surge_exposure_<date>.csv`, named with today's date, to
#'     `output_dir`, creating the folder if needed. Stops with an error if
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
#'   * `TRUE`: prints each step, the run time and the path of the saved file.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return A tibble, one row per water supply, with `wsd_source_id`,
#'   `park_unit`, `park_name`, `coast_elev_m` and:
#'   * `raw_surge_change_max`: the weighted change in surge (meters);
#'   * `raw_surge_pct_change_max`: the weighted percent change;
#'   * `norm_*` and `rank_*` versions of both from
#'     [add_standardized_indicators()]. A larger increase gets a higher score.
#' @export
#' @seealso [calc_storm_surge_sensitivity()] for the current surge level;
#'   [process_storm_surge_data()] and [process_coastal_geometry()], which make
#'   the inputs; [calc_inundation_exposure()] and [calc_swi_exposure()], the
#'   other coastal indicators.
#' @examplesIf file.exists(wavacore::wava_data_path("rapid", "inputs", "storm_surge", "storm_surge_RP10_muis_2023.gpkg"))
#' ws <- wavacore::get_wsd(park = "BISC")
#' surge_exp <- calc_storm_surge_exposure(
#'   ws,
#'   storm_surge_stations = wavacore::wava_data_path("rapid", "inputs", "storm_surge",
#'                                                   "storm_surge_RP10_muis_2023.gpkg")
#' )
#' surge_exp

calc_storm_surge_exposure <- function(water_supplies,
                                      storm_surge_stations,
                                      coast_geometry = wavacore::wava_data_path("rapid",
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
    dplyr::filter(!grepl("CMCC-CM2-VHR4", model)) %>% 
    dplyr::group_by(station_idx) %>%
    dplyr::mutate(raw_surge_change_max = future_mean - reference_mean,
           raw_surge_pct_change_max = ((future_mean - reference_mean) / reference_mean) * 100) %>%
    dplyr::slice_max(raw_surge_pct_change_max, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(station_idx, model, raw_surge_change_max, raw_surge_pct_change_max, lat, lon)
  
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
