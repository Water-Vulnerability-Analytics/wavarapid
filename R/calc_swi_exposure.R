#' @title Saltwater intrusion exposure
#'
#' @description
#' Gives each water supply the projected inland migration of the
#' saltwater-freshwater interface by 2100 (`Case_C`, from recharge change and
#' sea level rise combined, Adams et al. 2024) at its nearest saltwater
#' intrusion point, then divides it by the supply's distance to the coast.
#' The result is the fraction of the land between the supply and the coast
#' that is projected to be intruded, capped at 1. Negative values mean the
#' interface moves seaward.
#'
#' Supplies more than 50 km from the nearest point, or more than 50 km from
#' the coast (`coast_dist_m` is `NA`), get `NA`.
#'
#' @param water_supplies sf object. The water supplies, as points with
#'   `wsd_source_id`, `park_unit` and `park_name` columns, e.g. from
#'   [wavacore::get_wsd()]. Stops with an error if a column is missing.
#'   Required.
#' @param swi_pts String or sf object. The saltwater intrusion points from
#'   [process_swi_data()], with a `Case_C` column.
#'   * String: path to `swi_points.gpkg`. Stops with an error if the file
#'     doesn't exist.
#'   * sf object: the points, already read.
#'
#'   Default: `rapid/inputs/swi/swi_points.gpkg` in the shared data folder,
#'   see [wavacore::wava_data_path()].
#' @param coast_distances String. Path to `coastal_site_geometry.csv` from
#'   [process_coastal_geometry()], joined by `wsd_source_id`, `park_unit` and
#'   `park_name`; only `coast_dist_m` is used. Stops with an error if the
#'   file doesn't exist. Default:
#'   `rapid/inputs/coast_comp/coastal_site_geometry.csv` in the shared data
#'   folder.
#' @param save Logical. Whether to write the results to a CSV file.
#'   * `TRUE`: writes `swi_exposure_<date>.csv`, named with today's date, to
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
#'   * `TRUE`: prints each step, how many supplies got a value, the run time
#'     and the path of the saved file.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return A tibble, one row per water supply, with `wsd_source_id`,
#'   `park_unit`, `park_name`, `raw_swi_migration` (`Case_C` divided by the
#'   distance to the coast, at most 1, no units), and its `norm_` and `rank_`
#'   versions from [add_standardized_indicators()].
#' @references
#' Adams, K.H., Reager, J.T., Buzzanga, B., David, C., Sawyer, A., and
#'   Hamlington, B. (2024). Climate-induced saltwater intrusion in 2100:
#'   recharge-driven severity, sea level-driven prevalence. Geophysical
#'   Research Letters.
#' @export
#' @seealso [process_swi_data()] and [process_coastal_geometry()], which make
#'   the inputs; [calc_storm_surge_exposure()] and
#'   [calc_inundation_exposure()], the other coastal indicators.
#' @examplesIf file.exists(wavacore::wava_data_path("rapid", "inputs", "swi", "swi_points.gpkg"))
#' ws <- wavacore::get_wsd(park = "BISC")
#' swi_exp <- calc_swi_exposure(ws)
#' swi_exp

calc_swi_exposure <- function(
    water_supplies,
    swi_pts    = wavacore::wava_data_path("rapid", "inputs", "swi", "swi_points.gpkg"),
    coast_distances = wavacore::wava_data_path("rapid", "inputs", "coast_comp", "coastal_site_geometry.csv"),
    save       = FALSE,
    output_dir = NULL,
    verbose    = TRUE) {

  # Input validation 
  stopifnot(
    "water_supplies must be an sf object"         = inherits(water_supplies, "sf"),
    "save is TRUE but output_dir not provided"    = !save || !is.null(output_dir),
    "water_supplies missing wsd_source_id column" = "wsd_source_id" %in% names(water_supplies),
    "water_supplies missing park_unit column"     = "park_unit"     %in% names(water_supplies),
    "water_supplies missing park_name column"     = "park_name"     %in% names(water_supplies),
    "coast distances file does not exist" = file.exists(coast_distances)
  )

  if (verbose) start_time <- Sys.time()

  # Load SWI points 
  if (!inherits(swi_pts, "sf")) {
    stopifnot("swi_pts path does not exist \u2014 run process_swi_data() first" = file.exists(swi_pts))
    if (verbose) message("Loading SWI points ...")
    swi_pts <- sf::st_read(swi_pts, quiet = !verbose)
  }

  ws_proj  <- sf::st_transform(water_supplies, 5070)  # NAD83 / CONUS Albers (meters)
  swi_proj <- sf::st_transform(swi_pts, 5070)

  if (verbose) message(paste0("Finding nearest SWI point for ", nrow(ws_proj), " water supplies ..."))

  nearest_idx     <- sf::st_nearest_feature(ws_proj, swi_proj)
  dist_to_nearest <- as.numeric(sf::st_distance(ws_proj, swi_proj[nearest_idx, ], by_element = TRUE))

  swi_assigned <- water_supplies %>%
    sf::st_drop_geometry() %>%
    dplyr::select(wsd_source_id, park_unit, park_name) %>%
    dplyr::mutate(
      raw_swi_migration = dplyr::if_else(dist_to_nearest <= 50000,
                                  swi_pts$Case_C[nearest_idx],
                                  NA_real_)
    )

  if (verbose) {
    n_coastal <- sum(!is.na(swi_assigned$raw_swi_migration))
    message(paste0("  Assigned Case_C: ", n_coastal,
                   " | NA (> 50 km): ", nrow(swi_assigned) - n_coastal))
  }

  # Divide by distance to coast: raw_swi_migration / coast_dist_m gives the
  # fraction of the land-to-coast corridor projected to be intruded (dimensionless)
  coast_distances_df <- readr::read_csv(coast_distances, show_col_types = FALSE) %>%
    dplyr::select(wsd_source_id, park_unit, park_name, coast_dist_m)

  swi_weighted <- swi_assigned %>%
    dplyr::left_join(coast_distances_df, by = c("wsd_source_id", "park_unit", "park_name")) %>%
    dplyr::mutate(raw_swi_migration = raw_swi_migration / coast_dist_m,
           # if migration exceeds distance to coast, cap at 1:
           raw_swi_migration = pmin(raw_swi_migration, 1)) %>%
    dplyr::select(-coast_dist_m)
  
  # Standardize 
  results <- swi_weighted %>%
    add_standardized_indicators("raw_swi_migration", invert = FALSE)

  # Save 
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    filename <- file.path(output_dir, paste0("swi_exposure_", Sys.Date(), ".csv"))
    readr::write_csv(results, filename)
    if (verbose) message("Results saved to: ", filename)
  }

  if (verbose) {
    elapsed <- difftime(Sys.time(), start_time)
    message(paste0("\nSWI exposure completed in ", round(elapsed, 2), " ", attr(elapsed, "units")))
  }

  return(results)
}
