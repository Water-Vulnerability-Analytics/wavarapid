#' Calculate Withdrawal-to-Availability Exposure
#'
#' Calculates change in watershed water stress from the supplied Blumberg and
#' Warziniack WTA tables: 2040-2070 minus 2000-2020. WTA is withdrawals for
#' public supply, irrigation, and thermoelectric cooling divided by modeled
#' surface and subsurface runoff. This is a physical stress proxy, not a measure
#' of legal water rights or an Aqueduct score.
#'
#' The main indicator is ABSOLUTE change in the dimensionless WTA ratio: a
#' change from 0.2 to 0.3 is +0.1; a change from 0.3 to 0.2 is -0.1.
#' Positive change means increasing exposure; negative change means decreasing
#' exposure. Percent change is also returned; a zero baseline gives NA percent
#' change while preserving the actual change. Baseline and future levels retain
#' context for already-stressed watersheds. Ratios above 1 are valid and are not capped.
#'
#' Each model/scenario pair has equal weight. SSP1-2.6 (ssp126) is always
#' excluded, leaving 7 models x 3 scenarios (21 members) in the supplied
#' bundle. Following the existing exposure functions, select the
#' actual member nearest each type-7 p10/p50/p90 quantile, rather than return
#' interpolated changes. Ties use model, then scenario alphabetical order.
#' Future levels refer to those same selected members.
#' Agreement is the percentage of finite paired members with increasing WTA.
#'
#' HUC8 is the default watershed scale. WTA files also contain HUC2,
#' HUC4, and HUC6. HUC codes are padded to the selected width, never joined to
#' county FIPS. Missing watersheds/members remain missing. Coverage counts refer
#' to delivered WTA records; they do not verify completeness of the counties
#' aggregated upstream into those records. Naturalized local runoff does not
#' explicitly represent managed upstream inflows or interbasin transfers.
#'
#' Requires dplyr/stringr and src/rapid/helper_functions.R (loaded by setup.R).
#' Normalization and ranks are calculated across water supplies (one row per
#' supply), so a watershed with several supplies is counted once per supply.
#'
#' @param water_supplies sf POINT object with unique wsd_source_id, park_unit,
#'   and park_name. Empty points are retained with missing indicators.
#' @param data_dir Folder containing water_demand_projections, or that folder
#'   itself. Reads wta_ratios_2000_2020.csv and wta_ratios_2040_2070.csv.
#' @param watersheds sf polygons with a huc4 (or corresponding level) or huc
#'   column, or a WBD geodatabase path with a WBDHU4 (etc.) layer.
#' @param huc_level HUC width: 2, 4, 6, or 8 (default).
#' @param emission_scenarios Optional exact scenario names; NULL uses all except
#'   ssp126, which is not available.
#' @param climate_models Optional exact model names; NULL uses all.
#' @param max_watershed_distance_km Maximum nearest-polygon fallback distance.
#'   Default 5 km; 0 requires intersection. All polygons participate in matching,
#'   including those without WTA data. Boundary ties use the lowest HUC code.
#' @param save Whether to save a dated wta_exposure_huc08 CSV by default.
#' @param output_dir Output folder; required when save is TRUE.
#' @param verbose Whether to report timing and coverage.
#' @return Nonspatial tibble with one row per supply, raw_wta_change_p10/p50/p90,
#'   norm_/rank_ counterparts, raw_wta_models_showing_increase_pct and its
#'   standardized values, baseline_wta, corresponding future_wta_p10/p50/p90,
#'   raw_wta_percent_change_p10/p50/p90 with norm_/rank_ counterparts,
#'   selected member names, counts, and status.
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

  # Input validation
  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "water_supplies missing ID or park columns" =
      all(c("wsd_source_id", "park_unit", "park_name") %in% names(water_supplies)),
    "water_supplies must have rows and a known CRS" =
      nrow(water_supplies) > 0 && !is.na(sf::st_crs(water_supplies)),
    "wsd_source_id must be nonmissing and unique" =
      !anyNA(water_supplies$wsd_source_id) && !anyDuplicated(water_supplies$wsd_source_id),
    "water_supplies must contain POINT geometries" =
      all(as.character(sf::st_geometry_type(water_supplies)) == "POINT"),
    "huc_level must be 2, 4, 6, or 8" =
      length(huc_level) == 1L && huc_level %in% c(2L, 4L, 6L, 8L),
    "max_watershed_distance_km must be finite and nonnegative" =
      is.numeric(max_watershed_distance_km) && length(max_watershed_distance_km) == 1L &&
      is.finite(max_watershed_distance_km) && max_watershed_distance_km >= 0,
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir)
  )
  if (!exists("add_standardized_indicators", mode = "function")) {
    stop("Source src/rapid/helper_functions.R first, or run setup.R.")
  }
  start_time <- Sys.time()
  aggregation <- sprintf("HUC%02d", as.integer(huc_level))
  if (verbose) message("Starting WTA exposure calculation (", aggregation, ")...")

  # Preserve leading zeros in both the supplied tables and watershed boundaries.
  format_huc <- function(x) {
    x <- trimws(as.character(x))
    if (anyNA(x) || any(!grepl("^[0-9]+$", x) | nchar(x) > huc_level)) {
      stop("Invalid HUC identifiers for ", aggregation, ".")
    }
    sprintf(paste0("%0", huc_level, "d"), as.integer(x))
  }
  read_wta <- function(period, future = FALSE) {
    filename <- paste0("wta_ratios_", period, ".csv")
    path <- file.path(data_dir, filename)
    if (!file.exists(path)) path <- file.path(data_dir, "water_demand_projections", filename)
    if (!file.exists(path)) stop("Missing WTA file: ", path)
    x <- readr::read_csv(path, col_types = readr::cols(huc = readr::col_character(),
                                                    ratio = readr::col_double()),
                         show_col_types = FALSE)
    if (nrow(readr::problems(x))) stop(filename, " contains malformed data.")
    keys <- c("huc", if (future) c("climate_model", "emission_scenario"))
    if (!all(c(keys, "agg", "ratio") %in% names(x))) stop(filename, " missing required columns.")
    x <- x %>% filter(agg == aggregation) %>% select(all_of(c(keys, "ratio")))
    if (!nrow(x)) stop(filename, " has no ", aggregation, " records.")
    x$huc <- format_huc(x$huc)
    if (anyNA(x[keys]) || anyDuplicated(x[keys]) ||
        any(vapply(x[keys], function(z) any(!nzchar(trimws(z))), logical(1)))) {
      stop(filename, " contains missing identifiers or duplicate records.")
    }
    if (!is.numeric(x$ratio) || any(!is.na(x$ratio) & (!is.finite(x$ratio) | x$ratio < 0))) {
      stop(filename, " ratio must be nonnegative finite numeric values or NA.")
    }
    x
  }
  historical <- read_wta("2000_2020") %>%
    rename(baseline_wta = ratio) %>% mutate(historical_record = TRUE)
  # SSP1-2.6 is always excluded from the ensemble.
  future <- read_wta("2040_2070", future = TRUE) %>%
    filter(emission_scenario != "ssp126") %>%
    rename(future_wta = ratio) %>% mutate(future_record = TRUE)

  # Define the ensemble before restricting to watersheds containing supplies.
  selections <- list(emission_scenario = emission_scenarios, climate_model = climate_models)
  for (column in names(selections)) {
    requested <- selections[[column]]
    if (!is.null(requested) && (!is.character(requested) || !length(requested) ||
        anyNA(requested) || any(!requested %in% future[[column]]))) {
      stop("Unknown or empty ", column, " selection (ssp126 is excluded).")
    }
  }
  for (column in names(selections)) {
    if (!is.null(selections[[column]])) {
      future <- future %>% filter(.data[[column]] %in% selections[[column]])
    }
  }
  members <- future %>% distinct(climate_model, emission_scenario)
  if (!nrow(members)) stop("No model/scenario pairs match the selection.")

  # Assign the containing watershed, with a bounded fallback for coastal points.
  if (!inherits(watersheds, "sf")) {
    if (length(watersheds) != 1L || !file.exists(watersheds)) stop("WBD geodatabase does not exist.")
    watersheds <- sf::st_read(watersheds, layer = paste0("WBDHU", huc_level), quiet = TRUE)
  }
  stopifnot(
    "watersheds must be nonempty polygons with a known CRS" =
      nrow(watersheds) > 0 && !is.na(sf::st_crs(watersheds)) &&
      !any(sf::st_is_empty(watersheds)) &&
      all(as.character(sf::st_geometry_type(watersheds)) %in% c("POLYGON", "MULTIPOLYGON"))
  )
  id_column <- match(c(paste0("huc", huc_level), tolower(aggregation), "huc"),
                     tolower(names(watersheds)))
  id_column <- id_column[!is.na(id_column)]
  if (!length(id_column)) stop("watersheds needs a HUC column for the selected level.")
  huc <- format_huc(watersheds[[id_column[1]]])
  if (anyDuplicated(huc)) stop("watersheds must have one polygon row per HUC.")
  watersheds <- sf::st_transform(sf::st_make_valid(watersheds), 5070)
  points <- sf::st_transform(water_supplies, 5070)
  valid <- which(!sf::st_is_empty(points))
  index <- rep(NA_integer_, nrow(points))
  method <- rep("unmatched", nrow(points))
  distance <- rep(NA_real_, nrow(points))
  if (length(valid)) {
    hits <- sf::st_intersects(points[valid, ], watersheds)
    index[valid] <- vapply(hits, function(i) {
      if (length(i)) i[order(huc[i])][1] else NA_integer_
    }, integer(1))
    matched <- which(!is.na(index))
    method[matched] <- "intersects"
    distance[matched] <- 0
    unmatched <- valid[is.na(index[valid])]
    if (length(unmatched) && max_watershed_distance_km > 0) {
      nearest <- sf::st_nearest_feature(points[unmatched, ], watersheds)
      km <- as.numeric(sf::st_distance(points[unmatched, ], watersheds[nearest, ],
                                      by_element = TRUE)) / 1000
      distance[unmatched] <- km
      accept <- is.finite(km) & km <= max_watershed_distance_km
      index[unmatched[accept]] <- nearest[accept]
      method[unmatched[accept]] <- "nearest"
    }
  }
  locations <- water_supplies %>% sf::st_drop_geometry() %>%
    select(wsd_source_id, park_unit, park_name) %>%
    mutate(huc = huc[index], watershed_match_method = method, watershed_distance_km = distance)

  # Keep every expected member, including absent projections, to report coverage.
  changes <- tidyr::crossing(huc = unique(stats::na.omit(locations$huc)), members) %>%
    left_join(historical, by = "huc") %>%
    left_join(future, by = c("huc", "climate_model", "emission_scenario")) %>%
    # WTA is already a dimensionless withdrawal-to-availability ratio;
    # preserve the ratio change and do not divide it by watershed area.
    mutate(wta_change = future_wta - baseline_wta,
           percent_change = if_else(baseline_wta > 0, 100 * wta_change / baseline_wta, NA_real_)) %>%
    arrange(huc, climate_model, emission_scenario)

  # Select actual ensemble members, as in calc_runoff_exposure().
  select_percentile <- function(change, probability, value = change) {
    valid <- which(is.finite(change))
    if (!length(valid)) return(value[NA_integer_])
    target <- stats::quantile(change[valid], probability, type = 7, names = FALSE)
    distance <- abs(change[valid] - target)
    # Decimal roundoff must not decide an otherwise equal-distance tie.
    tolerance <- 16 * .Machine$double.eps * max(1, abs(change[valid]))
    nearest <- which(distance <= min(distance) + tolerance)[1]
    value[valid[nearest]]
  }
  coverage <- changes %>% group_by(huc) %>%
    summarise(baseline_wta = first(baseline_wta),
              historical_record = any(historical_record %in% TRUE),
              n_models = n(), n_future_records = sum(future_record %in% TRUE),
              n_valid_models = sum(is.finite(wta_change)),
              raw_wta_models_showing_increase_pct = if (n_valid_models > 0) {
                100 * sum(wta_change > 0, na.rm = TRUE) / n_valid_models
              } else NA_real_, .groups = "drop") %>%
    mutate(data_status = case_when(
      !historical_record ~ "missing_historical_watershed",
      !is.finite(baseline_wta) ~ "missing_historical_ratio",
      n_future_records == 0 ~ "missing_future_watershed",
      n_valid_models == 0 ~ "missing_future_ratio",
      n_valid_models < n_models ~ "partial_ensemble",
      TRUE ~ "ok")) %>% select(-historical_record)

  results <- coverage
  for (p in c(10, 50, 90)) {
    selected <- changes %>% group_by(huc) %>%
      summarise(raw_wta_change = select_percentile(wta_change, p / 100),
                raw_wta_percent_change = select_percentile(wta_change, p / 100, percent_change),
                future_wta = select_percentile(wta_change, p / 100, future_wta),
                climate_model = select_percentile(wta_change, p / 100, climate_model),
                emission_scenario = select_percentile(wta_change, p / 100, emission_scenario),
                .groups = "drop") %>%
      rename_with(~ paste0(.x, "_p", p), -huc)
    results <- results %>% left_join(selected, by = "huc")
  }
  # Standardize after joining so ranks are across water supplies.
  final_results <- locations %>% left_join(results, by = "huc") %>%
    add_standardized_indicators("raw_wta_change", invert = FALSE) %>%
    add_standardized_indicators("raw_wta_percent_change", invert = FALSE) %>%
    add_standardized_indicators("raw_wta_models_showing_increase_pct", invert = FALSE) %>%
    mutate(aggregation = aggregation, baseline_period = "2000-2020", future_period = "2040-2070",
           wta_change_units = "absolute WTA ratio change",
           selected_emission_scenarios = paste(sort(unique(members$emission_scenario)), collapse = ","),
           selected_climate_models = paste(sort(unique(members$climate_model)), collapse = ","),
           data_status = if_else(is.na(huc), "unmatched_watershed", data_status))

  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
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
  return(final_results)
}
