#' Calculate Saltwater Intrusion Exposure Indicator
#'
#' For each NPS water supply, assigns Case_C saltwater intrusion migration
#' (meters) using nearest-neighbor distance to the Adams et al. 2024 SWI
#' estimation points. Water supplies farther than 50 km from any SWI point
#' are assigned NA (not at coastal SWI risk). Those points that are assigned
#' SWI values are then weighted based on their distance to the coast, using the
#' same dataset which was produced in `process_coastal_weights()`.
#'
#' Case_C = combined recharge-change + sea-level-rise driven inland migration
#' of the saltwater-freshwater interface projected to 2100, in meters. Positive
#' values indicate expected inland intrusion (higher exposure); negative values
#' indicate the interface retreats.
#' 
#' @param water_supplies sf object with water supply point locations. Must
#'   contain columns `wsd_source_id`, `park_unit`, and `park_name`.
#' @param swi_pts Saltwater intrusion points. Character path to swi_points.gpkg produced by
#'   `process_swi_data()`. Defaults to `here("data", "rapid", "inputs", "swi", "swi_points.gpkg")`.
#' @param coastal_weights character path to coastal_weights.csv produced by
#'   `process_coastal_weights()`.
#' @param save logical, whether to save output CSV. Defaults to FALSE.
#' @param output_dir character, directory path for saving files. Required when
#'   `save = TRUE`.
#' @param verbose logical, whether to print progress messages. Defaults to TRUE.
#'
#' @return tibble with one row per water supply and columns:
#'   `wsd_source_id`, `park_unit`, `park_name`,
#'   `raw_swi_migration` (Case_C in meters; NA if > 50 km from any SWI point),
#'   `norm_swi_migration`, `rank_swi_migration`
#'
calc_swi_exposure <- function(
    water_supplies,
    swi_pts    = here("data", "rapid", "inputs", "swi", "swi_points.gpkg"),
    coastal_weights = here("data", "rapid", "inputs", "coast_comp", "coastal_weights.csv"),
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
    "coastal weights file does not exist" = file.exists(coastal_weights)
  )

  if (verbose) start_time <- Sys.time()

  # Load SWI points 
  if (!inherits(swi_pts, "sf")) {
    stopifnot("swi_pts path does not exist — run process_swi_data() first" = file.exists(swi_pts))
    if (verbose) message("Loading SWI points ...")
    swi_pts <- st_read(swi_pts, quiet = !verbose)
  }

  ws_proj  <- st_transform(water_supplies, 5070)  # NAD83 / CONUS Albers (meters)
  swi_proj <- st_transform(swi_pts, 5070)

  if (verbose) message(paste0("Finding nearest SWI point for ", nrow(ws_proj), " water supplies ..."))

  nearest_idx     <- st_nearest_feature(ws_proj, swi_proj)
  dist_to_nearest <- as.numeric(st_distance(ws_proj, swi_proj[nearest_idx, ], by_element = TRUE))

  swi_assigned <- water_supplies %>%
    st_drop_geometry() %>%
    select(wsd_source_id, park_unit, park_name) %>%
    mutate(
      raw_swi_migration = if_else(dist_to_nearest <= 50000,
                                  swi_pts$Case_C[nearest_idx],
                                  NA_real_)
    )

  if (verbose) {
    n_coastal <- sum(!is.na(swi_assigned$raw_swi_migration))
    message(paste0("  Assigned Case_C: ", n_coastal,
                   " | NA (> 50 km): ", nrow(swi_assigned) - n_coastal))
  }

  # Weight
  coastal_weights_df <- read_csv(coastal_weights, show_col_types = FALSE)
  
  swi_weighted <- swi_assigned %>%
    left_join(coastal_weights_df, by = c("wsd_source_id", "park_unit", "park_name")) %>% 
    mutate(raw_swi_migration = raw_swi_migration * coastal_weight)
  
  # Standardize 
  results <- swi_weighted %>%
    add_standardized_indicators("raw_swi_migration", invert = FALSE)

  # Save 
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    filename <- file.path(output_dir, paste0("swi_exposure_", Sys.Date(), ".csv"))
    write_csv(results, filename)
    if (verbose) message("Results saved to: ", filename)
  }

  if (verbose) {
    elapsed <- difftime(Sys.time(), start_time)
    message(paste0("\nSWI exposure completed in ", round(elapsed, 2), " ", attr(elapsed, "units")))
  }

  return(results)
}
