#' Process Coastal Geometry
#'
#' Combines the raw NOAA CUSP shoreline tiles into a single shoreline layer,
#' then computes distance to the nearest CUSP shoreline feature for each NPS
#' water supply. Water supplies farther than `max_coast_dist_m` (50 km) from the
#' coast receive `coast_dist_m = NA` and are dropped from the elevation
#' query entirely. This means those elevations also receive `coast_elev_m = NA`. 
#' Any downstream function that computes a weight from `coast_dist_m` and/or 
#' `coast_elev_m` (e.g. `calc_storm_surge_exposure()`, `calc_storm_surge_sensitivity()`,
#' `calc_swi_exposure()`, `calc_inundation_exposure()`) will therefore also
#' return NA for these non-coastal water supplies.
#'
#' Elevation for the remaining coastal supplies is queried point-by-point
#' from the USGS Elevation Point Query Service (EPQS) via
#' `elevatr::get_elev_point()`.
#'
#' When saved, writes coastal_site_geometry.csv (to `output_dir`) with
#' distance + elevation, and cusp_shorelines.gpkg (the combined shoreline
#' layer, cached alongside the raw tiles in `cusp_dir_path`).
#'
#' @param water_supplies sf object with water supply point locations
#' @param cusp_dir_path character, path to directory containing CUSP tile
#'   subfolders. The combined cusp_shorelines.gpkg is also cached here.
#' @param ncpu numeric, number of parallel workers used by
#'   `elevatr::get_elev_point()` for the elevation query. Defaults to 6.
#' @param max_coast_dist_m numeric, maximum distance (m) from the coast a
#'   water supply can be to receive a non-NA `coast_dist_m` / `coast_elev_m`.
#'   Defaults to 50000 (50 km).
#' @param save logical, whether to save output files. Defaults to TRUE.
#' @param output_dir character, directory path for saving
#'   coastal_site_geometry.csv. Required when `save = TRUE`.
#' @param verbose logical, whether to show progress messages.
#'
#' @return tibble (geometry dropped) with wsd_source_id, park_unit,
#'   park_name, coast_dist_m, coast_elev_m. 
#'
#' @export
process_coastal_geometry <- function(water_supplies,
                                     cusp_dir_path = wava_data_path("all", "coast_comp_data", "NOAA_CUSP"),
                                     ncpu = 6,
                                     max_coast_dist_m = 50000,
                                     save = TRUE,
                                     output_dir = wava_data_path("rapid", "inputs", "coast_comp"),
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
  
  tiles <- map(shp_files, ~ sf::st_read(.x, quiet = TRUE))
  cusp_shoreline <- do.call(rbind, tiles)
  
  cusp_proj <- sf::st_transform(cusp_shoreline, 5070)
  ws_proj <- sf::st_transform(water_supplies, 5070)
  
  nearest_idx <- st_nearest_feature(ws_proj, cusp_proj)
  dist_to_nearest <- st_distance(ws_proj, cusp_proj[nearest_idx, ], by_element = TRUE)
  
  coast_distances <- ws_proj %>%
    mutate(
      dist_to_nearest = as.numeric(dist_to_nearest),
      coast_dist_m    = if_else(dist_to_nearest > max_coast_dist_m, NA_real_, dist_to_nearest)
    ) %>%
    select(wsd_source_id, park_unit, park_name, coast_dist_m) 
  
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
  
  coast_geometries <- left_join(coast_distances, coast_elevations,
                                by = c("wsd_source_id", "park_unit", "park_name")) %>%
    sf::st_drop_geometry()
  
  if (save) {
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    geometry_path <- file.path(output_dir, "coastal_site_geometry.csv")
    write_csv(coast_geometries, geometry_path)
    
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
