#' Calculate Withdrawal-to-Availability Exposure
#'
#' Calculates change in watershed water stress from the Blumberg and Warziniack
#' WTA tables: 2040-2070 minus 2000-2020. WTA is withdrawals for public supply,
#' irrigation, and thermoelectric cooling divided by modeled surface and
#' subsurface runoff. It is a physical stress proxy, not a measure of legal
#' water rights or an Aqueduct score.
#'
#' The main indicator is ABSOLUTE change in the dimensionless WTA ratio (0.2 to
#' 0.3 is +0.1). Positive change means increasing exposure. Percent change is
#' also returned; a zero baseline gives NA percent change. Ratios above 1 are
#' valid and are not capped.
#'
#' Each model/scenario pair has equal weight. SSP1-2.6 (ssp126) is always
#' excluded, leaving 7 models x 3 scenarios (21 members). For each p10/p50/p90,
#' the actual member nearest the type-7 quantile is reported; ties use model,
#' then scenario, alphabetical order. Agreement is the percentage of finite
#' members with increasing WTA.
#'
#' HUC8 is the default watershed scale (HUC2, 4, and 6 are also available).
#' Missing watersheds/members remain missing. Naturalized local runoff does not
#' represent managed upstream inflows or interbasin transfers. norm_ and rank_
#' columns are calculated across water supplies.
#'
#' Requires dplyr/stringr and src/rapid/helper_functions.R (loaded by setup.R).
#'
#' @param water_supplies sf POINT object with unique wsd_source_id, park_unit,
#'   and park_name. Empty points are retained with missing indicators.
#' @param data_dir Folder containing water_demand_projections, or that folder
#'   itself. Reads wta_ratios_2000_2020.csv and wta_ratios_2040_2070.csv.
#' @param watersheds sf polygons with a huc8 (or matching level) or huc column,
#'   or a WBD geodatabase path with a WBDHU8 (etc.) layer.
#' @param huc_level HUC width: 2, 4, 6, or 8 (default).
#' @param emission_scenarios Optional scenario names; NULL uses all except
#'   ssp126, which is not available.
#' @param climate_models Optional model names; NULL uses all.
#' @param max_watershed_distance_km Maximum nearest-polygon fallback distance
#'   (default 5 km); 0 requires intersection. Boundary ties use the lowest HUC.
#' @param save Whether to save a dated wta_exposure_huc08 CSV by default.
#' @param output_dir Output folder; required when save is TRUE.
#' @param verbose Whether to report timing and coverage.
#' @return Nonspatial tibble with one row per supply, raw_wta_change_p10/p50/p90,
#'   raw_wta_percent_change_p10/p50/p90, raw_wta_models_showing_increase_pct,
#'   their norm_/rank_ columns, baseline_wta, future_wta_p10/p50/p90, selected
#'   model/scenario names, member counts, and data_status.
calc_wta_exposure <- function(
    water_supplies,
    data_dir = file.path("data", "Latest county-level data"),
    watersheds = file.path("data", "all", "WBD_National_GDB.gdb"),
    huc_level = 8L,
    emission_scenarios = NULL,
    climate_models = NULL,
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
  start_time <- Sys.time()
  aggregation <- sprintf("HUC%02d", as.integer(huc_level))  # e.g. "HUC08"
  huc_format <- paste0("%0", huc_level, "d")                # pads leading zeros
  if (verbose) message("Starting WTA exposure calculation (", aggregation, ")...")

  # 1. Read historical and future WTA ratios ---------------------------------
  if (!file.exists(file.path(data_dir, "wta_ratios_2000_2020.csv"))) {
    data_dir <- file.path(data_dir, "water_demand_projections")
  }
  historical <- readr::read_csv(file.path(data_dir, "wta_ratios_2000_2020.csv"),
                                col_types = readr::cols(huc = "c", ratio = "d"),
                                show_col_types = FALSE) |>
    filter(agg == aggregation) |>
    transmute(huc = sprintf(huc_format, as.integer(huc)),
              baseline_wta = ratio, historical_record = TRUE)

  future <- readr::read_csv(file.path(data_dir, "wta_ratios_2040_2070.csv"),
                            col_types = readr::cols(huc = "c", ratio = "d"),
                            show_col_types = FALSE) |>
    filter(agg == aggregation, emission_scenario != "ssp126") |>
    transmute(huc = sprintf(huc_format, as.integer(huc)),
              climate_model, emission_scenario,
              future_wta = ratio, future_record = TRUE)

  # 2. Choose ensemble members -----------------------------------------------
  if (!is.null(emission_scenarios)) {
    stopifnot("Unknown emission_scenarios (ssp126 is excluded)" =
                all(emission_scenarios %in% future$emission_scenario))
    future <- filter(future, emission_scenario %in% emission_scenarios)
  }
  if (!is.null(climate_models)) {
    stopifnot("Unknown climate_models" = all(climate_models %in% future$climate_model))
    future <- filter(future, climate_model %in% climate_models)
  }
  members <- distinct(future, climate_model, emission_scenario)
  if (!nrow(members)) stop("No model/scenario pairs match the selection.")

  # 3. Assign each supply to a watershed -------------------------------------
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
    if (length(i)) i[order(watershed_huc[i])][1] else NA_integer_
  }, integer(1))
  match_method[!is.na(huc_index)] <- "intersects"
  huc_distance[!is.na(huc_index)] <- 0

  # Points outside every watershed (e.g. coastal) use the nearest within the limit.
  outside <- has_point[is.na(huc_index[has_point])]
  if (length(outside) && max_watershed_distance_km > 0) {
    nearest <- sf::st_nearest_feature(points[outside, ], watersheds)
    km <- as.numeric(sf::st_distance(points[outside, ], watersheds[nearest, ],
                                     by_element = TRUE)) / 1000
    huc_distance[outside] <- km
    accept <- km <= max_watershed_distance_km
    huc_index[outside[accept]] <- nearest[accept]
    match_method[outside[accept]] <- "nearest"
  }

  locations <- water_supplies |>
    sf::st_drop_geometry() |>
    select(wsd_source_id, park_unit, park_name) |>
    mutate(huc = watershed_huc[huc_index],
           watershed_match_method = match_method,
           watershed_distance_km = huc_distance)

  # 4. Change from historical, per watershed and member ----------------------
  # Keep every expected member, including absent projections, to report coverage.
  changes <- tidyr::crossing(huc = unique(stats::na.omit(locations$huc)), members) |>
    left_join(historical, by = "huc") |>
    left_join(future, by = c("huc", "climate_model", "emission_scenario")) |>
    # WTA is already a ratio, so the change is not divided by area.
    mutate(wta_change = future_wta - baseline_wta,
           percent_change = if_else(baseline_wta > 0, 100 * wta_change / baseline_wta, NA_real_)) |>
    arrange(huc, climate_model, emission_scenario)

  # 5. Summarize the ensemble for each watershed -----------------------------
  results <- changes |>
    group_by(huc) |>
    summarise(baseline_wta = first(baseline_wta),
              historical_record = any(historical_record %in% TRUE),
              n_models = n(),
              n_future_records = sum(future_record %in% TRUE),
              n_valid_models = sum(is.finite(wta_change)),
              raw_wta_models_showing_increase_pct = if (n_valid_models > 0) {
                100 * sum(wta_change > 0, na.rm = TRUE) / n_valid_models
              } else NA_real_,
              .groups = "drop") |>
    mutate(data_status = case_when(
      !historical_record ~ "missing_historical_watershed",
      !is.finite(baseline_wta) ~ "missing_historical_ratio",
      n_future_records == 0 ~ "missing_future_watershed",
      n_valid_models == 0 ~ "missing_future_ratio",
      n_valid_models < n_models ~ "partial_ensemble",
      TRUE ~ "ok")) |>
    select(-historical_record)

  # For each percentile, keep the actual member nearest the percentile target.
  # A tiny tolerance stops floating-point roundoff from deciding equal ties.
  for (p in c(10, 50, 90)) {
    nearest_member <- changes |>
      filter(is.finite(wta_change)) |>
      group_by(huc) |>
      mutate(distance = abs(wta_change - stats::quantile(wta_change, p / 100,
                                                         type = 7, names = FALSE)),
             tolerance = 16 * .Machine$double.eps * max(1, abs(wta_change))) |>
      filter(distance <= min(distance) + tolerance) |>
      slice(1) |>
      ungroup() |>
      transmute(huc,
                raw_wta_change = wta_change,
                raw_wta_percent_change = percent_change,
                future_wta, climate_model, emission_scenario) |>
      rename_with(~ paste0(.x, "_p", p), -huc)
    results <- left_join(results, nearest_member, by = "huc")
  }

  # 6. Join to supplies and standardize across supplies ----------------------
  final_results <- locations |>
    left_join(results, by = "huc") |>
    add_standardized_indicators("raw_wta_change", invert = FALSE) |>
    add_standardized_indicators("raw_wta_percent_change", invert = FALSE) |>
    add_standardized_indicators("raw_wta_models_showing_increase_pct", invert = FALSE) |>
    mutate(aggregation = aggregation,
           baseline_period = "2000-2020",
           future_period = "2040-2070",
           wta_change_units = "absolute WTA ratio change",
           selected_emission_scenarios = paste(sort(unique(members$emission_scenario)), collapse = ","),
           selected_climate_models = paste(sort(unique(members$climate_model)), collapse = ","),
           data_status = if_else(is.na(huc), "unmatched_watershed", data_status))

  if (save) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    filename <- file.path(output_dir, paste0("wta_exposure_", tolower(aggregation), "_", Sys.Date(), ".csv"))
    readr::write_csv(final_results, filename)
    if (verbose) message("Results saved to: ", filename)
  }
  if (verbose) {
    message("WTA exposure: ", sum(is.finite(final_results$raw_wta_change_p50)), "/",
            nrow(final_results), " supplies with finite changes; ",
            sum(final_results$data_status == "partial_ensemble"), " with partial ensembles.")
    message("Completed in ", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2), " seconds.")
  }
  final_results
}
