#' Calculate historical Withdrawal-to-Availability sensitivity
#'
#' Returns the 2000-2020 watershed WTA ratio for each water supply. WTA is a
#' dimensionless withdrawal-to-availability ratio; it is not area-normalized.
#' Higher values indicate greater historical water-stress pressure. norm_ and
#' rank_ columns are calculated across water supplies.
#'
#' Requires src/rapid/helper_functions.R (loaded by setup.R).
#'
#' @param water_supplies sf POINT object with wsd_source_id, park_unit, and
#'   park_name.
#' @param data_dir Folder containing wta_ratios_2000_2020.csv, or its parent
#'   water_demand_projections folder.
#' @param watersheds sf polygons, or a WBD geodatabase path.
#' @param huc_level HUC width: 2, 4, 6, or 8 (default).
#' @param max_watershed_distance_km Maximum nearest-polygon fallback distance.
#' @param save Whether to save a dated WTA sensitivity CSV.
#' @param output_dir Output folder; required when save is TRUE.
#' @param verbose Whether to report coverage.
#' @return One row per supply with raw_wta, norm_wta, rank_wta, watershed
#'   matching metadata, and data_status.
calc_wta_sensitivity <- function(
    water_supplies,
    data_dir = file.path("data", "Latest county-level data"),
    watersheds = file.path("data", "all", "WBD_National_GDB.gdb"),
    huc_level = 8L,
    max_watershed_distance_km = 5,
    save = FALSE,
    output_dir = NULL,
    verbose = TRUE) {

  stopifnot(
    "water_supplies must be an sf object with a known CRS" =
      inherits(water_supplies, "sf") && !is.na(sf::st_crs(water_supplies)),
    "water_supplies needs wsd_source_id, park_unit, and park_name columns" =
      all(c("wsd_source_id", "park_unit", "park_name") %in% names(water_supplies)),
    "wsd_source_id must be nonmissing and unique" =
      !anyNA(water_supplies$wsd_source_id) && !anyDuplicated(water_supplies$wsd_source_id),
    "huc_level must be 2, 4, 6, or 8" = length(huc_level) == 1L && huc_level %in% c(2, 4, 6, 8),
    "max_watershed_distance_km must be finite and nonnegative" =
      is.finite(max_watershed_distance_km) && max_watershed_distance_km >= 0,
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir),
    "Source src/rapid/helper_functions.R first, or run setup.R" =
      exists("add_standardized_indicators", mode = "function")
  )
  aggregation <- sprintf("HUC%02d", as.integer(huc_level))  # e.g. "HUC08"
  huc_format <- paste0("%0", huc_level, "d")                # pads leading zeros

  # 1. Read historical WTA ratios --------------------------------------------
  if (!file.exists(file.path(data_dir, "wta_ratios_2000_2020.csv"))) {
    data_dir <- file.path(data_dir, "water_demand_projections")
  }
  historical <- readr::read_csv(file.path(data_dir, "wta_ratios_2000_2020.csv"),
                                col_types = readr::cols(huc = "c", ratio = "d"),
                                show_col_types = FALSE) |>
    dplyr::filter(agg == aggregation) |>
    dplyr::transmute(huc = sprintf(huc_format, as.integer(huc)), raw_wta = ratio)

  # 2. Assign each supply to a watershed -------------------------------------
  if (!inherits(watersheds, "sf")) {
    watersheds <- sf::st_read(watersheds, layer = paste0("WBDHU", huc_level), quiet = TRUE)
  }
  id_column <- match(c(paste0("huc", huc_level), tolower(aggregation), "huc"),
                     tolower(names(watersheds)))
  watershed_huc <- sprintf(huc_format, as.integer(watersheds[[stats::na.omit(id_column)[1]]]))
  # Distances in CONUS Albers metres.
  watersheds <- sf::st_transform(sf::st_make_valid(watersheds), 5070)
  points <- sf::st_transform(water_supplies, 5070)

  huc_index <- rep(NA_integer_, nrow(points))
  match_method <- rep("unmatched", nrow(points))
  huc_distance <- rep(NA_real_, nrow(points))
  has_point <- which(!sf::st_is_empty(points))

  # Containing watershed; a point on a shared boundary takes the lowest HUC.
  hits <- sf::st_intersects(points[has_point, ], watersheds)
  huc_index[has_point] <- vapply(hits, function(i) {
    if (length(i)) i[order(watershed_huc[i])][1L] else NA_integer_
  }, integer(1))
  match_method[!is.na(huc_index)] <- "intersects"
  huc_distance[!is.na(huc_index)] <- 0

  # Points outside every watershed (e.g. coastal) use the nearest within the limit.
  outside <- has_point[is.na(huc_index[has_point])]
  if (length(outside)) {
    nearest <- sf::st_nearest_feature(points[outside, ], watersheds)
    km <- as.numeric(sf::st_distance(points[outside, ], watersheds[nearest, ],
                                     by_element = TRUE)) / 1000
    huc_distance[outside] <- km
    accept <- km <= max_watershed_distance_km
    huc_index[outside[accept]] <- nearest[accept]
    match_method[outside[accept]] <- "nearest"
  }

  # 3. Join ratios to supplies and standardize across supplies ---------------
  results <- sf::st_drop_geometry(water_supplies) |>
    dplyr::select(wsd_source_id, park_unit, park_name) |>
    dplyr::mutate(huc = watershed_huc[huc_index],
                  watershed_match_method = match_method,
                  watershed_distance_km = huc_distance) |>
    dplyr::left_join(historical, by = "huc") |>
    dplyr::mutate(data_status = dplyr::case_when(
      is.na(huc) ~ "unmatched_watershed",
      is.na(raw_wta) ~ "missing_historical_ratio",
      TRUE ~ "ok")) |>
    add_standardized_indicators("raw_wta", invert = FALSE) |>
    dplyr::mutate(aggregation = aggregation,
                  historical_period = "2000-2020",
                  wta_units = "WTA ratio")

  if (save) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    filename <- file.path(output_dir, paste0("wta_sensitivity_", tolower(aggregation), "_", Sys.Date(), ".csv"))
    readr::write_csv(results, filename)
    if (verbose) message("Results saved to: ", filename)
  }
  if (verbose) {
    message("WTA sensitivity: ", sum(results$data_status == "ok"), "/",
            nrow(results), " supplies with historical WTA.")
  }
  results
}
