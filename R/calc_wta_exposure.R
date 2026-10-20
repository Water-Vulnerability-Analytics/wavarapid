#' Calculate Withdrawal-to-Availability (WTA) exposure
#'
#' Change in watershed water stress, 2040-2070 minus 2000-2020, from the
#' Blumberg and Warziniack WTA tables. WTA = withdrawals (public supply,
#' irrigation, thermoelectric) / modeled runoff; it is a dimensionless ratio, so
#' the change is not divided by area. Each model/scenario pair is one ensemble
#' member; ssp126 is excluded. For p10/p50/p90, reports the actual member
#' nearest the percentile (ties: model, then scenario, alphabetically).
#'
#' Supplies use the watershed they fall in, or the nearest within
#' max_watershed_distance_km. norm_/rank_ columns are across supplies.
#'
#' @param water_supplies sf points with wsd_source_id, park_unit, park_name.
#' @param data_dir Folder with wta_ratios_2000_2020.csv and wta_ratios_2040_2070.csv.
#' @param watersheds sf polygons, or a WBD geodatabase path.
#' @param huc_level HUC width: 2, 4, 6, or 8 (default).
#' @param emission_scenarios,climate_models Optional ensemble subsets.
#' @param max_watershed_distance_km Nearest-watershed limit; 0 = intersect only.
#' @param save,output_dir Save wta_exposure_hucNN_YYYY-MM-DD.csv to output_dir.
#' @param verbose Report coverage.
#' @return Tibble with one row per supply.
calc_wta_exposure <- function(water_supplies,
                              data_dir = "data/Latest county-level data/water_demand_projections",
                              watersheds = "data/all/WBD_National_GDB.gdb",
                              huc_level = 8L,
                              emission_scenarios = NULL,
                              climate_models = NULL,
                              max_watershed_distance_km = 5,
                              save = FALSE,
                              output_dir = NULL,
                              verbose = TRUE) {

  aggregation <- sprintf("HUC%02d", huc_level)  # e.g. "HUC08"
  huc_format <- paste0("%0", huc_level, "d")    # pads leading zeros

  # 1. Read WTA ratios
  historical <- read_csv(file.path(data_dir, "wta_ratios_2000_2020.csv"),
                         col_types = cols(huc = "c", ratio = "d")) |>
    filter(agg == aggregation) |>
    transmute(huc = sprintf(huc_format, as.integer(huc)), baseline_wta = ratio,
              historical_record = TRUE)

  future <- read_csv(file.path(data_dir, "wta_ratios_2040_2070.csv"),
                     col_types = cols(huc = "c", ratio = "d")) |>
    filter(agg == aggregation, emission_scenario != "ssp126") |>
    transmute(huc = sprintf(huc_format, as.integer(huc)), climate_model,
              emission_scenario, future_wta = ratio, future_record = TRUE)
  if (!is.null(emission_scenarios)) future <- filter(future, emission_scenario %in% emission_scenarios)
  if (!is.null(climate_models)) future <- filter(future, climate_model %in% climate_models)
  members <- distinct(future, climate_model, emission_scenario)

  # 2. Assign each supply to a watershed (lowest HUC on a shared boundary)
  if (!inherits(watersheds, "sf")) {
    watersheds <- st_read(watersheds, layer = paste0("WBDHU", huc_level), quiet = TRUE)
  }
  watersheds <- watersheds |>
    mutate(huc = sprintf(huc_format, as.integer(.data[[paste0("huc", huc_level)]]))) |>
    st_make_valid() |>
    st_transform(5070)
  points <- st_transform(water_supplies, 5070)

  index <- sapply(st_intersects(points, watersheds),
                  function(i) if (length(i)) i[order(watersheds$huc[i])][1] else NA_integer_)
  distance <- ifelse(is.na(index), NA_real_, 0)
  method <- ifelse(is.na(index), "unmatched", "intersects")

  outside <- which(is.na(index) & !st_is_empty(points))
  if (length(outside) && max_watershed_distance_km > 0) {
    nearest <- st_nearest_feature(points[outside, ], watersheds)
    km <- as.numeric(st_distance(points[outside, ], watersheds[nearest, ], by_element = TRUE)) / 1000
    distance[outside] <- km
    index[outside[km <= max_watershed_distance_km]] <- nearest[km <= max_watershed_distance_km]
    method[outside[km <= max_watershed_distance_km]] <- "nearest"
  }

  locations <- st_drop_geometry(water_supplies) |>
    select(wsd_source_id, park_unit, park_name) |>
    mutate(huc = watersheds$huc[index],
           watershed_match_method = method,
           watershed_distance_km = distance)

  # 3. Change per watershed and member (all members kept, to report coverage)
  changes <- crossing(huc = unique(na.omit(locations$huc)), members) |>
    left_join(historical, by = "huc") |>
    left_join(future, by = c("huc", "climate_model", "emission_scenario")) |>
    mutate(wta_change = future_wta - baseline_wta,
           percent_change = if_else(baseline_wta > 0, 100 * wta_change / baseline_wta, NA_real_)) |>
    arrange(huc, climate_model, emission_scenario)

  # 4. Summarize the ensemble per watershed
  results <- changes |>
    group_by(huc) |>
    summarise(baseline_wta = first(baseline_wta),
              historical_record = any(historical_record %in% TRUE),
              n_models = n(),
              n_future_records = sum(future_record %in% TRUE),
              n_valid_models = sum(is.finite(wta_change)),
              raw_wta_models_showing_increase_pct = if (n_valid_models > 0) {
                100 * sum(wta_change > 0, na.rm = TRUE) / n_valid_models
              } else NA_real_) |>
    mutate(data_status = case_when(
      !historical_record ~ "missing_historical_watershed",
      !is.finite(baseline_wta) ~ "missing_historical_ratio",
      n_future_records == 0 ~ "missing_future_watershed",
      n_valid_models == 0 ~ "missing_future_ratio",
      n_valid_models < n_models ~ "partial_ensemble",
      TRUE ~ "ok")) |>
    select(-historical_record)

  # Member nearest each percentile (a tiny tolerance ignores float roundoff in ties)
  for (p in c(10, 50, 90)) {
    nearest_member <- changes |>
      filter(is.finite(wta_change)) |>
      group_by(huc) |>
      mutate(distance = abs(wta_change - quantile(wta_change, p / 100))) |>
      filter(distance <= min(distance) + 16 * .Machine$double.eps * max(1, abs(wta_change))) |>
      slice(1) |>
      ungroup() |>
      transmute(huc, raw_wta_change = wta_change, raw_wta_percent_change = percent_change,
                future_wta, climate_model, emission_scenario) |>
      rename_with(~ paste0(.x, "_p", p), -huc)
    results <- left_join(results, nearest_member, by = "huc")
  }

  # 5. Join to supplies and standardize
  final <- locations |>
    left_join(results, by = "huc") |>
    add_standardized_indicators("raw_wta_change") |>
    add_standardized_indicators("raw_wta_percent_change") |>
    add_standardized_indicators("raw_wta_models_showing_increase_pct") |>
    mutate(aggregation = aggregation,
           baseline_period = "2000-2020",
           future_period = "2040-2070",
           wta_change_units = "absolute WTA ratio change",
           selected_emission_scenarios = paste(sort(unique(members$emission_scenario)), collapse = ","),
           selected_climate_models = paste(sort(unique(members$climate_model)), collapse = ","),
           data_status = if_else(is.na(huc), "unmatched_watershed", data_status))

  if (save) {
    filename <- file.path(output_dir, paste0("wta_exposure_", tolower(aggregation), "_", Sys.Date(), ".csv"))
    write_csv(final, filename)
    if (verbose) message("Results saved to: ", filename)
  }
  if (verbose) {
    message("WTA exposure: ", sum(is.finite(final$raw_wta_change_p50)), "/", nrow(final),
            " supplies with finite changes")
  }
  final
}
