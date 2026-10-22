#' @title Distance to the coast and elevation of each water supply
#'
#' @description
#' Combines the NOAA CUSP shoreline tiles into one shoreline layer and finds
#' the distance from each water supply to the nearest shoreline feature. For
#' supplies within `max_coast_dist_m` of the coast it then gets the elevation
#' from the USGS Elevation Point Query Service (EPQS), through
#' [elevatr::get_elev_point()], so it needs internet.
#'
#' Supplies farther than `max_coast_dist_m` get `NA` for both distance and
#' elevation. The coastal indicators that read this table
#' ([calc_storm_surge_exposure()], [calc_storm_surge_sensitivity()],
#' [calc_swi_exposure()] and [calc_inundation_exposure()]) give those supplies
#' `NA`.
#'
#' @param water_supplies sf object. The water supplies, as points with
#'   `wsd_source_id`, `park_unit` and `park_name` columns, e.g. from
#'   [wavacore::get_wsd()]. Stops with an error if a column is missing.
#'   Required.
#' @param cusp_dir_path String. The folder with the CUSP tile subfolders;
#'   every `.shp` file below it is read. Stops with an error if the folder
#'   doesn't exist or has no shapefiles. With `save = TRUE`, the combined
#'   shoreline is also written here. Default: `all/coast_comp_data/NOAA_CUSP`
#'   in the shared data folder, see [wavacore::wava_data_path()].
#' @param ncpu Integer. The number of parallel workers for the elevation
#'   queries. Default: `6`.
#' @param max_coast_dist_m Numeric. The largest distance from the coast, in
#'   meters, at which a supply counts as coastal.
#'   * Supplies within it get `coast_dist_m` and an elevation.
#'   * Supplies beyond it get `NA` for both and are not sent to the elevation
#'     service.
#'
#'   Default: `50000` (50 km).
#' @param save Logical. Whether to write the results.
#'   * `TRUE`: writes `coastal_site_geometry.csv` to `output_dir` (creating
#'     the folder if needed) and the combined shoreline,
#'     `cusp_shorelines.gpkg`, to `cusp_dir_path`. Both replace existing
#'     files. Stops with an error if `output_dir` is `NULL`.
#'   * `FALSE`: writes nothing.
#'
#'   Default: `TRUE`.
#' @param output_dir String or `NULL`. The folder for
#'   `coastal_site_geometry.csv`.
#'   * A folder path: where the file is written when `save = TRUE`. Ignored
#'     when `save = FALSE`.
#'   * `NULL`: only allowed with `save = FALSE`.
#'
#'   Default: `rapid/inputs/coast_comp` in the shared data folder, which is
#'   where the coastal indicator functions look for it.
#' @param verbose Logical. Whether to print progress.
#'   * `TRUE`: prints messages for each step, the run time and the paths of
#'     the saved files.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return A tibble (no geometry), one row per water supply, with columns
#'   `wsd_source_id`, `park_unit`, `park_name`, `coast_dist_m` (meters to the
#'   nearest shoreline) and `coast_elev_m` (elevation in meters).
#' @export
#' @seealso [calc_storm_surge_exposure()], [calc_storm_surge_sensitivity()],
#'   [calc_swi_exposure()] and [calc_inundation_exposure()], which read
#'   `coastal_site_geometry.csv`.
#' @examplesIf interactive() && dir.exists(wavacore::wava_data_path("all", "coast_comp_data", "NOAA_CUSP"))
#' # Reads every CUSP shoreline tile (about 3 GB) and queries the USGS
#' # elevation service, so it needs internet and takes several minutes
#' ws <- wavacore::get_wsd(park = "BISC")
#' geometry <- process_coastal_geometry(ws, save = FALSE)
#' geometry

process_coastal_geometry <- function(water_supplies,
                                     cusp_dir_path = wavacore::wava_data_path("all", "coast_comp_data", "NOAA_CUSP"),
                                     ncpu = 6,
                                     max_coast_dist_m = 50000,
                                     save = TRUE,
                                     output_dir = wavacore::wava_data_path("rapid", "inputs", "coast_comp"),
                                     verbose = TRUE) {
  
  # Input validation
  stopifnot(
    "water_supplies must be an sf object"       = inherits(water_supplies, "sf"),
    "water_supplies must include wsd_source_id" = "wsd_source_id" %in% names(water_supplies),
    "water_supplies must include park_unit"     = "park_unit"     %in% names(water_supplies),
    "water_supplies must include park_name"     = "park_name"     %in% names(water_supplies),
    "cusp_dir_path does not exist"              = dir.exists(cusp_dir_path),
    "save is TRUE but output_dir not provided"  = !save || !is.null(output_dir)
  )
  
  if (verbose) {
    message("Starting coastal geometry processing...")
    start_time <- Sys.time()
  }
  
  # Grab coastal water supply sites' distance to the coast
  
  shp_files <- list.files(cusp_dir_path, pattern = "\\.shp$",
                          full.names = TRUE, recursive = TRUE)
  
  if (length(shp_files) == 0) stop("No shapefiles found in: ", cusp_dir_path)
  
  if (verbose) {
    message("Reading and combining ", length(shp_files), " CUSP shoreline tiles...")
  }
  
  tiles <- purrr::map(shp_files, ~ sf::st_read(.x, quiet = TRUE))
  cusp_shoreline <- do.call(rbind, tiles)
  
  cusp_proj <- sf::st_transform(cusp_shoreline, 5070)
  ws_proj <- sf::st_transform(water_supplies, 5070)
  
  nearest_idx <- sf::st_nearest_feature(ws_proj, cusp_proj)
  dist_to_nearest <- sf::st_distance(ws_proj, cusp_proj[nearest_idx, ], by_element = TRUE)
  
  coast_distances <- ws_proj %>%
    dplyr::mutate(
      dist_to_nearest = as.numeric(dist_to_nearest),
      coast_dist_m    = dplyr::if_else(dist_to_nearest > max_coast_dist_m, NA_real_, dist_to_nearest)
    ) %>%
    dplyr::select(wsd_source_id, park_unit, park_name, coast_dist_m) 
  
  # Grab relevant coastal water supply sites' elevations
  
  if (verbose) {
    message("Extracting elevation for ",  
            nrow(dplyr::filter(coast_distances, !is.na(coast_dist_m))),
            " coastal water supplies...")
  }
  
  coast_elevations <- elevatr::get_elev_point(
    locations = dplyr::filter(coast_distances, !is.na(coast_dist_m)),
    src = "epqs", 
    units = "meters",
    ncpu = ncpu,
    serial = FALSE
  ) %>%
    sf::st_drop_geometry() %>%
    dplyr::select(wsd_source_id, park_unit, park_name, 
                  coast_elev_m = elevation)
  
  coast_geometries <- dplyr::left_join(coast_distances, coast_elevations,
                                by = c("wsd_source_id", "park_unit", "park_name")) %>%
    sf::st_drop_geometry()
  
  if (save) {
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    geometry_path <- file.path(output_dir, "coastal_site_geometry.csv")
    readr::write_csv(coast_geometries, geometry_path)
    
    shoreline_path <- file.path(cusp_dir_path, "cusp_shorelines.gpkg")
    sf::st_write(cusp_shoreline, shoreline_path,
                 delete_dsn = TRUE, quiet = TRUE)
    
    if (verbose) {
      message("Coastal geometry saved to: ", geometry_path)
      message("Combined shoreline cached to: ", shoreline_path)
    }
  }
  
  if (verbose) {
    total_time <- Sys.time() - start_time
    message(sprintf("Coastal geometry processing complete in %.1f %s.",
                    as.numeric(total_time), attr(total_time, "units")))
  }

  return(coast_geometries)
}
