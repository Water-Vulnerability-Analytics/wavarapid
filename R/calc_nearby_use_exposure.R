#' Calculate nearby-use exposure from county demand projections
#'
#' Compares 2040-2070 with 2000-2020 county withdrawals (public supply + crop
#' irrigation + thermoelectric, in MGD) from Blumberg and Warziniack. Each
#' model/scenario pair is one ensemble member; ssp126 is excluded. For p10/p50/
#' p90, reports the actual member nearest the percentile of absolute change.
#'
#' - Absolute change: (future - historical) / county area, in MGD/km2.
#' - Percent change: 100 * (future - historical) / historical MGD.
#' - A missing month, or a power plant expected in that period but not
#'   reported, makes a county total NA. Plants retired out of the projections
#'   are not expected in the future period.
#'
#' Supplies use the 2020 TIGER/Line county they fall in, or the nearest county
#' within max_county_distance_km. norm_/rank_ columns are across supplies.
#'
#' @param water_supplies sf points with wsd_source_id, park_unit, park_name.
#' @param data_dir Folder with the water_demand_projections CSVs.
#' @param save,output_dir Save nearby_use_exposure_YYYY-MM-DD.csv to output_dir.
#' @param verbose Report the saved filename.
#' @param months Calendar months to average (default all 12).
#' @param emission_scenarios,climate_models Optional ensemble subsets.
#' @param max_county_distance_km Nearest-county limit for points outside counties.
#' @return Tibble with one row per supply.
calc_nearby_use_exposure <- function(water_supplies,
                                     data_dir = "data/Latest county-level data/water_demand_projections",
                                     save = FALSE,
                                     output_dir = NULL,
                                     verbose = TRUE,
                                     months = 1:12,
                                     emission_scenarios = NULL,
                                     climate_models = NULL,
                                     max_county_distance_km = 5) {

  days_in_month <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

  # 1. Assign each supply to a county
  counties <- tigris::counties(cb = FALSE, year = 2020, progress_bar = FALSE) |>
    st_make_valid()
  counties$area_km2 <- as.numeric(counties$ALAND) / 1e6 # land area; polygons keep water for matching
  counties <- st_transform(counties, 4326)
  points <- st_transform(water_supplies, 4326)

  nearest <- st_nearest_feature(points, counties)
  km <- as.numeric(st_distance(points, counties[nearest, ], by_element = TRUE)) / 1000
  nearest[is.na(km) | km > max_county_distance_km] <- NA

  locations <- st_drop_geometry(water_supplies) |>
    select(any_of(c("wsd_source_id", "park_unit", "park_name"))) |>
    mutate(fips = counties$GEOID[nearest],
           county_area_km2 = counties$area_km2[nearest],
           county_distance_km = km,
           county_source = "2020 Census TIGER/Line counties") |>
    as_tibble()

  # 2. Read historical and future withdrawals ("historical" is its own member)
  pub_cols <- c("FIPS", "Month", "population", "pub_per_cap", "irr")
  tep_cols <- c("FIPS", "plant_id", "Month", "tep")
  model_cols <- c("climate_model", "emission_scenario")

  pub <- bind_rows(
    fread(file.path(data_dir, "pub_irr_monthly_avg_2000_2020.csv"), select = pub_cols) |>
      mutate(climate_model = "historical", emission_scenario = "historical"),
    fread(file.path(data_dir, "pub_irr_monthly_avg_2040_2070.csv"), select = c(pub_cols, model_cols))) |>
    transmute(fips = sprintf("%05d", FIPS), Month, climate_model, emission_scenario,
              pub_irr_mgd = pub_per_cap * population / 1000 + irr,  # population is thousands
              in_pub_file = TRUE)

  tep <- bind_rows(
    fread(file.path(data_dir, "tep_monthly_avg_2000_2020.csv"), select = tep_cols) |>
      mutate(climate_model = "historical", emission_scenario = "historical"),
    fread(file.path(data_dir, "tep_monthly_avg_2040_2070.csv"), select = c(tep_cols, model_cols))) |>
    mutate(fips = sprintf("%05d", FIPS))

  # Plants a county should report, scoped to the period being evaluated.
  # A plant in the historical file but absent from the projections has been
  # retired, not lost: Blumberg and Warziniack set future withdrawals to zero
  # for coal and gas plants in states with 100% clean-energy targets, and for
  # any plant with no withdrawals over 2018-2020. Expecting the projections to
  # report those plants would make the county NA for every member and month.
  hist_plants <- tep |> filter(climate_model == "historical") |> distinct(fips, plant_id) |>
    count(fips, name = "hist_plants")
  future_plants <- tep |> filter(climate_model != "historical") |> distinct(fips, plant_id) |>
    count(fips, name = "future_plants")

  # 3. Ensemble members
  members <- pub |> filter(climate_model != "historical", emission_scenario != "ssp126") |>
    distinct(climate_model, emission_scenario)
  if (!is.null(emission_scenarios)) members <- filter(members, emission_scenario %in% emission_scenarios)
  if (!is.null(climate_models)) members <- filter(members, climate_model %in% climate_models)

  # 4. County withdrawals per member, averaged over months (weighted by days)
  totals <- expand_grid(fips = unique(na.omit(locations$fips)), Month = months) |>
    cross_join(bind_rows(tibble(climate_model = "historical", emission_scenario = "historical"),
                         members)) |>
    left_join(pub, by = c("fips", "Month", model_cols)) |>
    left_join(tep |> group_by(fips, Month, climate_model, emission_scenario) |>
                summarise(tep_mgd = sum(tep), n_plants = n(), .groups = "drop"),
              by = c("fips", "Month", model_cols)) |>
    left_join(hist_plants, by = "fips") |>
    left_join(future_plants, by = "fips") |>
    mutate(n_plants = coalesce(n_plants, 0L),
           n_expected = coalesce(if_else(climate_model == "historical", hist_plants, future_plants), 0L),
           # No plants expected this period = zero; a plant that should report
           # this period but doesn't (a month or member gap) = NA
           tep_mgd = case_when(n_expected == 0 ~ 0, n_plants < n_expected ~ NA, TRUE ~ tep_mgd),
           total_mgd = pub_irr_mgd + tep_mgd) |>
    group_by(fips, climate_model, emission_scenario) |>
    summarise(withdrawal_mgd = weighted.mean(total_mgd, days_in_month[Month]),
              n_public_months = sum(!is.na(in_pub_file)),
              n_complete_months = sum(is.finite(total_mgd)),
              missing_plant_records = sum(n_expected - n_plants),
              .groups = "drop")

  historical <- totals |>
    filter(climate_model == "historical") |>
    transmute(fips, historic_withdrawal_mgd = withdrawal_mgd,
              historical_complete_months = n_complete_months,
              historical_public_months = n_public_months)

  # 5. Change from historical, per member
  projected <- totals |>
    filter(climate_model != "historical") |>
    left_join(historical, by = "fips") |>
    left_join(distinct(locations, fips, county_area_km2), by = "fips") |>
    mutate(absolute_change = (withdrawal_mgd - historic_withdrawal_mgd) / county_area_km2,
           percent_change = if_else(historic_withdrawal_mgd > 0,
                                    100 * (withdrawal_mgd - historic_withdrawal_mgd) / historic_withdrawal_mgd,
                                    NA_real_))

  # 6. Summarize the ensemble per county
  n_months <- length(months)
  result <- projected |>
    group_by(fips) |>
    summarise(n_valid_models = sum(is.finite(absolute_change)),
              n_model_members = nrow(members),
              raw_models_showing_increase_pct = if (n_valid_models > 0) {
                100 * sum(absolute_change[is.finite(absolute_change)] > 0) / n_valid_models
              } else NA_real_,
              missing_future_plant_records = sum(missing_plant_records),
              future_public_months = max(n_public_months),
              future_complete_months = max(n_complete_months),
              .groups = "drop") |>
    left_join(historical, by = "fips") |>
    mutate(data_status = case_when(
      historical_public_months == 0 ~ "missing_historical_county",
      !is.finite(historic_withdrawal_mgd) | historical_public_months < n_months |
        historical_complete_months < n_months ~ "incomplete_historical_data",
      future_public_months == 0 ~ "missing_future_county",
      missing_future_plant_records > 0 | future_public_months < n_months |
        future_complete_months < n_months ~ "incomplete_future_data",
      n_valid_models == 0 ~ "no_complete_model_members",
      n_valid_models < n_model_members ~ "partial_ensemble",
      historic_withdrawal_mgd == 0 ~ "zero_historical_withdrawal",
      TRUE ~ "ok"))

  # Member nearest each percentile of absolute change
  for (p in c(10, 50, 90)) {
    nearest_member <- projected |>
      filter(is.finite(absolute_change)) |>
      group_by(fips) |>
      slice_min(abs(absolute_change - quantile(absolute_change, p / 100)), n = 1, with_ties = FALSE) |>
      ungroup() |>
      transmute(fips,
                !!paste0("raw_use_change_p", p) := percent_change,
                !!paste0("raw_absolute_use_change_p", p) := absolute_change,
                !!paste0("future_withdrawal_mgd_p", p) := withdrawal_mgd)
    result <- left_join(result, nearest_member, by = "fips")
  }

  # 7. Join to supplies and standardize
  final <- locations |>
    left_join(result, by = "fips") |>
    mutate(data_status = coalesce(data_status, "unmatched_county")) |>
    add_standardized_indicators("raw_")

  if (save) {
    filename <- file.path(output_dir, paste0("nearby_use_exposure_", Sys.Date(), ".csv"))
    write_csv(final, filename)
    if (verbose) message("Results saved to: ", filename)
  }
  final
}
