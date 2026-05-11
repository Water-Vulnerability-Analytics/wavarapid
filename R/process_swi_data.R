#' Process Adams et al. 2024 Saltwater Intrusion Dataset
#'
#' Reads output_data.txt from the Adams et al. 2024 (GRL) dataset, filters to
#' CONUS extent, and saves one output:
#'   1. `swi_points.gpkg` — raw SWI estimation points with Case_C (meters)
#'
#' Case_C = combined recharge-change + sea-level-rise driven inland migration
#' of the saltwater-freshwater interface projected to 2100, in meters.
#'
#' Note: The published SWI point locations originate from a legacy data pipeline
#' (Sawyer/Zhou et al.) and do not have a one-to-one match with any available
#' HydroSHEDS watershed shapefile version.
#'
#' Citation: Adams, K.H.; J.T. Reager; B. Buzzanga; C. David; A. Sawyer;
#' B. Hamlington, 2024, Climate-induced saltwater intrusion in 2100:
#' recharge-driven severity, sea level-driven prevalence, Geophysical
#' Research Letters.
#'
#' @param data_path character, path to output_data.txt from the Adams et al. 2024 dataset.
#'   Defaults to `here("data", "all", "coast_comp_data", "adams_swi_2024", "data", "output_data.txt")`.
#' @param output_dir character, directory where gpkg file is written.
#'   Defaults to `here("data", "rapid", "inputs", "swi")`.
#'
process_swi_data <- function(
    data_path  = here("data", "all", "coast_comp_data", "adams_swi_2024", "data", "output_data.txt"),
    output_dir = here("data", "rapid", "inputs", "swi")) {

  stopifnot("output_data.txt not found — check data_path" = file.exists(data_path))

  swi_pts <- read.table(data_path, header = TRUE, sep = " ") %>%
    as_tibble() %>%
    rename(watershed_id = WatershedID) %>%
    # Case_C: meters of inland SWI migration by 2100 (combined recharge + SLR)
    select(watershed_id, lat, lon, Case_C) %>%
    filter(lon >= -125, lon <= -66, lat >= 24, lat <= 50) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  pts_path <- file.path(output_dir, "swi_points.gpkg")
  st_write(swi_pts, pts_path, delete_dsn = TRUE, quiet = TRUE)

}
