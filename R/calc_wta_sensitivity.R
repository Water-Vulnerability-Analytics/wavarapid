#' Calculate historical Withdrawal-to-Availability sensitivity
#'
#' Returns the 2000-2020 watershed WTA ratio for each water supply. WTA is a
#' dimensionless withdrawal-to-availability ratio; it is not area-normalized.
#' Higher values indicate greater historical water-stress pressure.
#'
#' @param water_supplies sf POINT object with wsd_source_id and park_unit.
#' @param data_dir Folder containing wta_ratios_2000_2020.csv, or its parent
#'   water_demand_projections folder.
#' @param watersheds sf polygons, or a WBD geodatabase path.
#' @param huc_level HUC width: 2, 4, 6, or 8 (default).
#' @param max_watershed_distance_km Maximum nearest-polygon fallback distance.
#' @param save Whether to save a dated WTA sensitivity CSV.
#' @param output_dir Output folder; required when save is TRUE.
#' @param verbose Whether to report coverage and timing.
#' @return One row per supply with raw_wta, norm_wta, rank_wta, matching
#'   metadata, and historical data status.
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
    inherits(water_supplies, "sf"), nrow(water_supplies) > 0,
    !is.na(sf::st_crs(water_supplies)),
    all(c("wsd_source_id", "park_unit", "park_name") %in% names(water_supplies)),
    !anyNA(water_supplies$wsd_source_id), !anyDuplicated(water_supplies$wsd_source_id),
    all(as.character(sf::st_geometry_type(water_supplies)) == "POINT"),
    huc_level %in% c(2L, 4L, 6L, 8L),
    is.numeric(max_watershed_distance_km), length(max_watershed_distance_km) == 1L,
    is.finite(max_watershed_distance_km), max_watershed_distance_km >= 0,
    !save || !is.null(output_dir))
  if (!exists("add_standardized_indicators", mode = "function")) {
    stop("Source src/rapid/helper_functions.R first, or run setup.R.")
  }
  aggregation <- sprintf("HUC%02d", as.integer(huc_level))
  format_huc <- function(x) {
    x <- trimws(as.character(x))
    if (anyNA(x) || any(!grepl("^[0-9]+$", x) | nchar(x) > huc_level)) {
      stop("Invalid HUC identifiers for ", aggregation, ".")
    }
    sprintf(paste0("%0", huc_level, "d"), as.integer(x))
  }
  filename <- "wta_ratios_2000_2020.csv"
  path <- file.path(data_dir, filename)
  if (!file.exists(path)) path <- file.path(data_dir, "water_demand_projections", filename)
  if (!file.exists(path)) stop("Missing WTA file: ", path)
  historical <- readr::read_csv(path, col_types = readr::cols(huc = readr::col_character(),
                                                                ratio = readr::col_double()),
                                 show_col_types = FALSE)
  if (!all(c("huc", "agg", "ratio") %in% names(historical))) {
    stop(filename, " is missing required columns.")
  }
  historical <- historical |>
    dplyr::filter(agg == aggregation) |>
    dplyr::transmute(huc = format_huc(huc), raw_wta = ratio)
  if (!nrow(historical) || anyDuplicated(historical$huc)) {
    stop(filename, " has no unique ", aggregation, " records.")
  }
  if (any(!is.na(historical$raw_wta) & (!is.finite(historical$raw_wta) | historical$raw_wta < 0))) {
    stop(filename, " ratio must be nonnegative finite numeric values or NA.")
  }

  if (!inherits(watersheds, "sf")) {
    if (length(watersheds) != 1L || !file.exists(watersheds)) stop("WBD geodatabase does not exist.")
    watersheds <- sf::st_read(watersheds, layer = paste0("WBDHU", huc_level), quiet = TRUE)
  }
  id_column <- match(c(paste0("huc", huc_level), tolower(aggregation), "huc"),
                     tolower(names(watersheds)))
  id_column <- id_column[!is.na(id_column)]
  stopifnot(nrow(watersheds) > 0, length(id_column) > 0,
            all(as.character(sf::st_geometry_type(watersheds)) %in% c("POLYGON", "MULTIPOLYGON")))
  huc <- format_huc(watersheds[[id_column[1L]]])
  if (anyDuplicated(huc)) stop("watersheds must have one polygon row per HUC.")
  watersheds <- sf::st_transform(sf::st_make_valid(watersheds), 5070)
  points <- sf::st_transform(water_supplies, 5070)
  index <- rep(NA_integer_, nrow(points)); method <- rep("unmatched", nrow(points)); distance <- rep(NA_real_, nrow(points))
  valid <- which(!sf::st_is_empty(points))
  if (length(valid)) {
    hits <- sf::st_intersects(points[valid, ], watersheds)
    index[valid] <- vapply(hits, function(i) if (length(i)) i[order(huc[i])][1L] else NA_integer_, integer(1))
    matched <- which(!is.na(index)); method[matched] <- "intersects"; distance[matched] <- 0
    missing <- valid[is.na(index[valid])]
    if (length(missing)) {
      nearest <- sf::st_nearest_feature(points[missing, ], watersheds)
      km <- as.numeric(sf::st_distance(points[missing, ], watersheds[nearest, ], by_element = TRUE)) / 1000
      distance[missing] <- km; accept <- is.finite(km) & km <= max_watershed_distance_km
      index[missing[accept]] <- nearest[accept]; method[missing[accept]] <- "nearest"
    }
  }
  locations <- sf::st_drop_geometry(water_supplies) |>
    dplyr::select(wsd_source_id, park_unit, park_name) |>
    dplyr::mutate(huc = huc[index], watershed_match_method = method,
                  watershed_distance_km = distance)
  results <- locations |>
    dplyr::left_join(historical, by = "huc") |>
    dplyr::mutate(data_status = dplyr::case_when(
      is.na(huc) ~ "unmatched_watershed",
      is.na(raw_wta) ~ "missing_historical_ratio",
      TRUE ~ "ok")) |>
    add_standardized_indicators("raw_wta", invert = FALSE) |>
    dplyr::mutate(aggregation = aggregation, historical_period = "2000-2020",
                  wta_units = "WTA ratio")
  if (save) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    out <- file.path(output_dir, paste0("wta_sensitivity_", tolower(aggregation), "_", Sys.Date(), ".csv"))
    readr::write_csv(results, out)
    if (verbose) message("Results saved to: ", out)
  }
  if (verbose) message("WTA sensitivity: ", sum(results$data_status == "ok"), "/", nrow(results), " supplies with historical WTA.")
  results
}
