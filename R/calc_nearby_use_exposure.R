#' Calculate nearby-use exposure from county demand projections
#'
#' Compares 2040-2070 with 2000-2020 county withdrawals from the Blumberg and
#' Warziniack monthly climatologies. Withdrawals combine public supply, crop
#' irrigation, and thermoelectric cooling, in millions of gallons/day (MGD).
#'
#' Each climate model/emission scenario pair is one equally weighted ensemble
#' member. SSP1-2.6 (ssp126) is always excluded, leaving 7 models x 3 pathways
#' (21 members). For each p10/p50/p90, the actual member whose absolute change
#' per county area is nearest the type-7 percentile target is reported, with
#' percent change and future withdrawal from that same member.
#'
#' - Absolute change: (future - historical) / county area, in MGD/km2.
#' - Percent change: 100 * (future - historical) / historical MGD.
#' - Model agreement: share of complete members with a positive change.
#'
#' Every selected month and every known power plant is required for a complete
#' county total; missing data are never treated as zero. A known plant absent
#' from the future file makes the future total missing (not assumed retired).
#' Counties with no plants in either file contribute zero thermoelectric use.
#'
#' Each supply uses the 2020 TIGER/Line county containing it, or the nearest
#' county within max_county_distance_km. County areas include water territory.
#' norm_ and rank_ columns are calculated across water supplies.
#'
#' Requires src/rapid/helper_functions.R (loaded by setup.R).
#'
#' @param water_supplies sf POINT locations with unique wsd_source_id, park_unit,
#'   and optionally park_name. Empty points are retained with NA indicators.
#' @param data_dir Directory containing water_demand_projections, or that folder
#'   itself.
#' @param save Whether to save nearby_use_exposure_YYYY-MM-DD.csv.
#' @param output_dir Required when save is TRUE.
#' @param verbose Whether to report the saved filename.
#' @param months Calendar months to average (default all 12).
#' @param emission_scenarios Optional pathway names, e.g. "ssp245".
#' @param climate_models Optional climate model names. NULL uses all.
#' @param max_county_distance_km Maximum distance to the nearest county for
#'   points outside every county polygon (default 5 km). Inf allows any distance.
#' @return Tibble, one row per supply: county fips, area and distance;
#'   raw_use_change_p10/p50/p90 (percent change),
#'   raw_absolute_use_change_p10/p50/p90 (MGD/km2), future_withdrawal_mgd_p10/
#'   p50/p90, raw_models_showing_increase_pct, historic_withdrawal_mgd, member
#'   and month counts, data_status, and norm_/rank_ columns.
calc_nearby_use_exposure <- function(water_supplies,
                                     data_dir = file.path("data", "Latest county-level data"),
                                     save = FALSE,
                                     output_dir = NULL,
                                     verbose = TRUE,
                                     months = 1:12,
                                     emission_scenarios = NULL,
                                     climate_models = NULL,
                                     max_county_distance_km = 5) {

  stopifnot(
    "water_supplies must be an sf object with a known CRS" =
      inherits(water_supplies, "sf") && !is.na(sf::st_crs(water_supplies)),
    "water_supplies needs wsd_source_id and park_unit columns" =
      all(c("wsd_source_id", "park_unit") %in% names(water_supplies)),
    "wsd_source_id must be nonmissing and unique" =
      !anyNA(water_supplies$wsd_source_id) && !anyDuplicated(water_supplies$wsd_source_id),
    "months must be unique values from 1 to 12" =
      length(months) > 0 && all(months %in% 1:12) && !anyDuplicated(months),
    "save is TRUE but output_dir was not provided" = !save || !is.null(output_dir)
  )
  months <- sort(as.integer(months))
  # Climatological month lengths, used to average daily rates across months.
  days_in_month <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  if (!file.exists(file.path(data_dir, "pub_irr_monthly_avg_2000_2020.csv"))) {
    data_dir <- file.path(data_dir, "water_demand_projections")
  }

  # 1. Assign each supply to a county --------------------------------------
  counties <- tigris::counties(cb = FALSE, year = 2020, progress_bar = FALSE) |>
    sf::st_make_valid()
  county_area <- as.numeric(sf::st_area(sf::st_transform(counties, 5070))) / 1e6
  # Work in lon/lat so distances are geodesic metres.
  counties <- sf::st_transform(counties, 4326)
  points <- sf::st_transform(water_supplies, 4326)

  county_index <- rep(NA_integer_, nrow(points))
  county_distance <- rep(NA_real_, nrow(points))
  has_point <- which(!sf::st_is_empty(points))

  # Containing county; a point on a shared boundary takes the lowest FIPS.
  hits <- sf::st_intersects(points[has_point, ], counties, model = "closed")
  county_index[has_point] <- vapply(hits, function(i) {
    if (length(i)) i[order(counties$GEOID[i])][1L] else NA_integer_
  }, integer(1))
  county_distance[!is.na(county_index)] <- 0

  # Points outside every county use the nearest county within the limit.
  outside <- has_point[is.na(county_index[has_point])]
  if (length(outside)) {
    nearest <- sf::st_nearest_feature(points[outside, ], counties)
    km <- as.numeric(sf::st_distance(points[outside, ], counties[nearest, ],
                                     by_element = TRUE)) / 1000
    county_distance[outside] <- km
    accept <- km <= max_county_distance_km
    county_index[outside[accept]] <- nearest[accept]
  }

  locations <- sf::st_drop_geometry(water_supplies) |>
    dplyr::select(dplyr::any_of(c("wsd_source_id", "park_unit", "park_name"))) |>
    dplyr::mutate(fips = counties$GEOID[county_index],
                  county_area_km2 = county_area[county_index],
                  county_distance_km = county_distance,
                  county_source = "2020 Census TIGER/Line counties") |>
    tibble::as_tibble()
  supply_fips <- unique(stats::na.omit(locations$fips))

  # 2. Read withdrawals for the supply counties ------------------------------
  pub_columns <- c("population", "pub_per_cap", "irr")
  model_columns <- c("climate_model", "emission_scenario")

  hist_pub <- data.table::fread(
    file.path(data_dir, "pub_irr_monthly_avg_2000_2020.csv"),
    select = c("FIPS", "Month", pub_columns),
    colClasses = list(character = "FIPS", numeric = pub_columns),
    data.table = FALSE, showProgress = FALSE)
  hist_tep <- data.table::fread(
    file.path(data_dir, "tep_monthly_avg_2000_2020.csv"),
    select = c("FIPS", "plant_id", "Month", "tep"),
    colClasses = list(character = c("FIPS", "plant_id"), numeric = "tep"),
    data.table = FALSE, showProgress = FALSE)
  future_pub <- data.table::fread(
    file.path(data_dir, "pub_irr_monthly_avg_2040_2070.csv"),
    select = c("FIPS", "Month", model_columns, pub_columns),
    colClasses = list(character = c("FIPS", model_columns), numeric = pub_columns),
    data.table = FALSE, showProgress = FALSE)
  future_tep <- data.table::fread(
    file.path(data_dir, "tep_monthly_avg_2040_2070.csv"),
    select = c("FIPS", "plant_id", "Month", model_columns, "tep"),
    colClasses = list(character = c("FIPS", "plant_id", model_columns), numeric = "tep"),
    data.table = FALSE, showProgress = FALSE)

  # Five-character FIPS text, matching the county GEOIDs.
  hist_pub$fips <- sprintf("%05d", as.integer(hist_pub$FIPS))
  hist_tep$fips <- sprintf("%05d", as.integer(hist_tep$FIPS))
  future_pub$fips <- sprintf("%05d", as.integer(future_pub$FIPS))
  future_tep$fips <- sprintf("%05d", as.integer(future_tep$FIPS))

  # 3. Choose ensemble members (always excluding ssp126) ---------------------
  members <- tibble::as_tibble(unique(future_pub[model_columns])) |>
    dplyr::filter(emission_scenario != "ssp126")
  if (!is.null(emission_scenarios)) {
    stopifnot("Unknown emission_scenarios (ssp126 is excluded)" =
                all(emission_scenarios %in% members$emission_scenario))
    members <- dplyr::filter(members, emission_scenario %in% emission_scenarios)
  }
  if (!is.null(climate_models)) {
    stopifnot("Unknown climate_models" = all(climate_models %in% members$climate_model))
    members <- dplyr::filter(members, climate_model %in% climate_models)
  }
  if (!nrow(members)) stop("No model/scenario pairs match the selection.")

  # 4. Historical county withdrawals -----------------------------------------
  # One row per county and month. Public supply MGD = gallons/person/day x
  # population (thousands) / 1000. Thermoelectric is the sum of plants.
  hist_plants <- dplyr::count(dplyr::distinct(hist_tep, fips, plant_id), fips,
                              name = "n_expected_plants")
  hist_plant_totals <- hist_tep |>
    dplyr::group_by(fips, Month) |>
    dplyr::summarise(thermoelectric_mgd = sum(tep), n_plant_records = dplyr::n(),
                     .groups = "drop")

  historical <- tidyr::expand_grid(fips = supply_fips, Month = months) |>
    dplyr::left_join(hist_pub |>
                       dplyr::transmute(fips, Month,
                                        public_supply_mgd = pub_per_cap * population / 1000,
                                        irrigation_mgd = irr, public_record = TRUE),
                     by = c("fips", "Month")) |>
    dplyr::left_join(hist_plant_totals, by = c("fips", "Month")) |>
    dplyr::left_join(hist_plants, by = "fips") |>
    dplyr::mutate(
      public_record = dplyr::coalesce(public_record, FALSE),
      n_expected_plants = dplyr::coalesce(n_expected_plants, 0L),
      n_plant_records = dplyr::coalesce(n_plant_records, 0L),
      # No known plants is a true zero; a missing known plant is missing data.
      thermoelectric_mgd = dplyr::case_when(
        n_expected_plants == 0L ~ 0,
        n_plant_records < n_expected_plants ~ NA_real_,
        TRUE ~ thermoelectric_mgd),
      total_mgd = public_supply_mgd + irrigation_mgd + thermoelectric_mgd,
      days = days_in_month[Month]) |>
    # Average the monthly daily rates, weighted by month length.
    dplyr::group_by(fips) |>
    dplyr::summarise(historic_withdrawal_mgd = stats::weighted.mean(total_mgd, days),
                     historical_complete_months = sum(is.finite(total_mgd)),
                     historical_public_months = sum(public_record),
                     .groups = "drop")

  # 5. Future county withdrawals, per ensemble member ------------------------
  # Plants known in either period are expected in the future file.
  future_plants <- dplyr::bind_rows(dplyr::distinct(future_tep, fips, plant_id),
                                    dplyr::distinct(hist_tep, fips, plant_id)) |>
    dplyr::distinct(fips, plant_id) |>
    dplyr::count(fips, name = "n_expected_plants")
  future_plant_totals <- future_tep |>
    dplyr::group_by(fips, Month, climate_model, emission_scenario) |>
    dplyr::summarise(thermoelectric_mgd = sum(tep), n_plant_records = dplyr::n(),
                     .groups = "drop")
  member_keys <- c("fips", "Month", model_columns)

  projected <- tidyr::expand_grid(fips = supply_fips, Month = months) |>
    dplyr::cross_join(members) |>
    dplyr::left_join(future_pub |>
                       dplyr::transmute(fips, Month, climate_model, emission_scenario,
                                        public_supply_mgd = pub_per_cap * population / 1000,
                                        irrigation_mgd = irr, public_record = TRUE),
                     by = member_keys) |>
    dplyr::left_join(future_plant_totals, by = member_keys) |>
    dplyr::left_join(future_plants, by = "fips") |>
    dplyr::mutate(
      public_record = dplyr::coalesce(public_record, FALSE),
      n_expected_plants = dplyr::coalesce(n_expected_plants, 0L),
      n_plant_records = dplyr::coalesce(n_plant_records, 0L),
      thermoelectric_mgd = dplyr::case_when(
        n_expected_plants == 0L ~ 0,
        n_plant_records < n_expected_plants ~ NA_real_,
        TRUE ~ thermoelectric_mgd),
      total_mgd = public_supply_mgd + irrigation_mgd + thermoelectric_mgd,
      days = days_in_month[Month]) |>
    dplyr::group_by(fips, climate_model, emission_scenario) |>
    dplyr::summarise(withdrawal_mgd = stats::weighted.mean(total_mgd, days),
                     n_public_months = sum(public_record),
                     n_complete_months = sum(is.finite(total_mgd)),
                     missing_plant_records = sum(n_expected_plants - n_plant_records),
                     .groups = "drop")

  # 6. Change from historical, per ensemble member ---------------------------
  projected <- projected |>
    dplyr::left_join(historical, by = "fips") |>
    dplyr::left_join(dplyr::distinct(locations, fips, county_area_km2), by = "fips") |>
    dplyr::mutate(
      absolute_change = (withdrawal_mgd - historic_withdrawal_mgd) / county_area_km2,
      percent_change = dplyr::if_else(
        historic_withdrawal_mgd > 0,
        100 * (withdrawal_mgd - historic_withdrawal_mgd) / historic_withdrawal_mgd,
        NA_real_))

  # 7. Summarize the ensemble for each county --------------------------------
  months_needed <- length(months)
  result <- projected |>
    dplyr::group_by(fips) |>
    dplyr::summarise(
      n_valid_models = sum(is.finite(absolute_change)),
      n_model_members = nrow(members),
      raw_models_showing_increase_pct = if (n_valid_models > 0L) {
        100 * sum(absolute_change[is.finite(absolute_change)] > 0) / n_valid_models
      } else NA_real_,
      missing_future_plant_records = sum(missing_plant_records, na.rm = TRUE),
      future_public_months = max(n_public_months, na.rm = TRUE),
      future_complete_months = max(n_complete_months, na.rm = TRUE),
      .groups = "drop") |>
    dplyr::left_join(historical, by = "fips") |>
    dplyr::mutate(data_status = dplyr::case_when(
      historical_public_months == 0 ~ "missing_historical_county",
      !is.finite(historic_withdrawal_mgd) |
        historical_public_months < months_needed |
        historical_complete_months < months_needed ~ "incomplete_historical_data",
      future_public_months == 0 ~ "missing_future_county",
      missing_future_plant_records > 0 |
        future_public_months < months_needed |
        future_complete_months < months_needed ~ "incomplete_future_data",
      n_valid_models == 0 ~ "no_complete_model_members",
      n_valid_models < n_model_members ~ "partial_ensemble",
      historic_withdrawal_mgd == 0 ~ "zero_historical_withdrawal",
      TRUE ~ "ok"))

  # For each percentile, keep the actual member nearest the percentile target.
  for (p in c(10, 50, 90)) {
    nearest_member <- projected |>
      dplyr::filter(is.finite(absolute_change)) |>
      dplyr::group_by(fips) |>
      dplyr::mutate(distance = abs(absolute_change -
                                     stats::quantile(absolute_change, p / 100,
                                                     names = FALSE, type = 7))) |>
      dplyr::slice_min(distance, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::transmute(fips,
                       !!paste0("raw_use_change_p", p) := percent_change,
                       !!paste0("raw_absolute_use_change_p", p) := absolute_change,
                       !!paste0("future_withdrawal_mgd_p", p) := withdrawal_mgd)
    result <- dplyr::left_join(result, nearest_member, by = "fips")
  }

  # 8. Join to supplies and standardize across supplies ----------------------
  final <- locations |>
    dplyr::left_join(result, by = "fips") |>
    dplyr::mutate(data_status = dplyr::coalesce(data_status, "unmatched_county")) |>
    add_standardized_indicators("raw_", invert = FALSE)

  if (save) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    filename <- file.path(output_dir, paste0("nearby_use_exposure_", Sys.Date(), ".csv"))
    readr::write_csv(final, filename)
    if (verbose) message("Results saved to: ", filename)
  }

  final
}
