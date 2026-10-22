#' @title Nearby use exposure: projected change in county water withdrawals
#'
#' @description
#' Compares county water withdrawals in 2040-2070 with 2000-2020 (public
#' supply, crop irrigation and thermoelectric, in million gallons per day,
#' MGD), from the Blumberg and Warziniack county demand projections. Each
#' climate model and emission scenario pair is one ensemble member; ssp126 is
#' left out.
#'
#' Each water supply gets the 2020 TIGER/Line county it is in (downloaded
#' with [tigris::counties()], so it needs internet), or the nearest county
#' within `max_county_distance_km`. For each county and member, the monthly
#' withdrawals in `months` are averaged, weighted by days, and compared with
#' the historical average:
#' * absolute change: (future - historical) / county land area, MGD per km2;
#' * percent change: 100 * (future - historical) / historical.
#'
#' For the 10th, 50th and 90th percentiles of the absolute change, the
#' member closest to each is reported. A missing month, or a power plant
#' that should report in that period but doesn't, makes a county total `NA`.
#' Plants retired in the projections (no withdrawals in the future file) are
#' not expected in the future period. Gaps stay `NA` rather than being read
#' as zero use. The `norm_` and `rank_` columns are relative to the water
#' supplies passed in.
#'
#' @param water_supplies sf object. The water supplies, as points with
#'   `wsd_source_id`, `park_unit` and `park_name` columns, e.g. from
#'   [wavacore::get_wsd()]. Required.
#' @param data_dir String. The folder with `pub_irr_monthly_avg_2000_2020.csv`,
#'   `pub_irr_monthly_avg_2040_2070.csv`, `tep_monthly_avg_2000_2020.csv` and
#'   `tep_monthly_avg_2040_2070.csv`. Default:
#'   `Latest county-level data/water_demand_projections` in the shared data
#'   folder, see [wavacore::wava_data_path()].
#' @param save Logical. Whether to write the results to a CSV file.
#'   * `TRUE`: writes `nearby_use_exposure_<date>.csv`, named with today's
#'     date, to `output_dir`. The folder must already exist; `output_dir`
#'     can't be `NULL`.
#'   * `FALSE`: writes nothing.
#'
#'   Default: `FALSE`.
#' @param output_dir String or `NULL`. The folder for the CSV file when
#'   `save = TRUE`; ignored otherwise. Default: `NULL`.
#' @param verbose Logical. Whether to print the path of the saved file.
#'   * `TRUE`: prints it when `save = TRUE`.
#'   * `FALSE`: prints nothing.
#'
#'   Default: `TRUE`.
#' @param months Integer vector. The calendar months (1-12) to average, e.g.
#'   `6:8` for summer. Default: `1:12`.
#' @param emission_scenarios Character vector or `NULL`. The emission
#'   scenarios to use.
#'   * A character vector, e.g. `"ssp245"`: keeps only those members.
#'   * `NULL`: keeps every scenario except ssp126.
#'
#'   Default: `NULL`.
#' @param climate_models Character vector or `NULL`. The climate models to
#'   use, in the same way as `emission_scenarios`. Default: `NULL` (all).
#' @param max_county_distance_km Numeric. How far, in km, a supply outside
#'   every county can be from the nearest county and still use it. Supplies
#'   farther away get `NA` and `data_status = "unmatched_county"`. Default:
#'   `5`.
#'
#' @return A tibble, one row per water supply, with `wsd_source_id`,
#'   `park_unit`, `park_name`, the county (`fips`, `county_area_km2`,
#'   `county_distance_km`, `county_source`) and:
#'   * `raw_use_change_p10`, `_p50`, `_p90`: the percent change of the member
#'     closest to that percentile of absolute change;
#'   * `raw_absolute_use_change_p10`, `_p50`, `_p90`: that member's absolute
#'     change (MGD per km2), and `future_withdrawal_mgd_p10`, `_p50`, `_p90`;
#'   * `raw_models_showing_increase_pct`: the percentage of members with an
#'     increase;
#'   * `historic_withdrawal_mgd` and coverage counts (`n_valid_models`,
#'     `n_model_members`, `missing_future_plant_records`, and the numbers of
#'     public and complete months in each period);
#'   * `data_status`: `"ok"`, `"partial_ensemble"`, `"unmatched_county"`,
#'     `"missing_historical_county"`, `"missing_future_county"`,
#'     `"incomplete_historical_data"`, `"incomplete_future_data"`,
#'     `"no_complete_model_members"` or `"zero_historical_withdrawal"` (the
#'     percent change is `NA`);
#'   * `norm_*` and `rank_*` versions of every `raw_` column from
#'     [add_standardized_indicators()].
#' @export
#' @seealso [calc_nearby_use_sensitivity()] for historical use;
#'   [calc_wta_exposure()] for the change in watershed water stress from the
#'   same projections.
#' @examplesIf interactive() && dir.exists(wavacore::wava_data_path("Latest county-level data", "water_demand_projections"))
#' # Downloads the county boundaries and reads about 140 MB of projections,
#' # so it needs internet and takes a few minutes
#' ws <- wavacore::get_wsd(park = c("ZION", "ROMO"))
#' use_exp <- calc_nearby_use_exposure(ws)
#' use_exp[, c("wsd_source_id", "fips", "raw_use_change_p90", "data_status")]

calc_nearby_use_exposure <- function(water_supplies,
                                     data_dir = wavacore::wava_data_path("Latest county-level data", "water_demand_projections"),
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
    sf::st_make_valid()
  counties$area_km2 <- as.numeric(counties$ALAND) / 1e6 # land area; polygons keep water for matching
  counties <- sf::st_transform(counties, 4326)
  points <- sf::st_transform(water_supplies, 4326)

  nearest <- sf::st_nearest_feature(points, counties)
  km <- as.numeric(sf::st_distance(points, counties[nearest, ], by_element = TRUE)) / 1000
  nearest[is.na(km) | km > max_county_distance_km] <- NA

  locations <- sf::st_drop_geometry(water_supplies) |>
    dplyr::select(dplyr::any_of(c("wsd_source_id", "park_unit", "park_name"))) |>
    dplyr::mutate(fips = counties$GEOID[nearest],
           county_area_km2 = counties$area_km2[nearest],
           county_distance_km = km,
           county_source = "2020 Census TIGER/Line counties") |>
    dplyr::as_tibble()

  # 2. Read historical and future withdrawals ("historical" is its own member)
  pub_cols <- c("FIPS", "Month", "population", "pub_per_cap", "irr")
  tep_cols <- c("FIPS", "plant_id", "Month", "tep")
  model_cols <- c("climate_model", "emission_scenario")

  pub <- dplyr::bind_rows(
    data.table::fread(file.path(data_dir, "pub_irr_monthly_avg_2000_2020.csv"), select = pub_cols) |>
      dplyr::mutate(climate_model = "historical", emission_scenario = "historical"),
    data.table::fread(file.path(data_dir, "pub_irr_monthly_avg_2040_2070.csv"), select = c(pub_cols, model_cols))) |>
    dplyr::transmute(fips = sprintf("%05d", FIPS), Month, climate_model, emission_scenario,
              pub_irr_mgd = pub_per_cap * population / 1000 + irr,  # population is thousands
              in_pub_file = TRUE)

  tep <- dplyr::bind_rows(
    data.table::fread(file.path(data_dir, "tep_monthly_avg_2000_2020.csv"), select = tep_cols) |>
      dplyr::mutate(climate_model = "historical", emission_scenario = "historical"),
    data.table::fread(file.path(data_dir, "tep_monthly_avg_2040_2070.csv"), select = c(tep_cols, model_cols))) |>
    dplyr::mutate(fips = sprintf("%05d", FIPS))

  # Plants a county should report, scoped to the period being evaluated.
  # A plant in the historical file but absent from the projections has been
  # retired, not lost: Blumberg and Warziniack set future withdrawals to zero
  # for coal and gas plants in states with 100% clean-energy targets, and for
  # any plant with no withdrawals over 2018-2020. Expecting the projections to
  # report those plants would make the county NA for every member and month.
  hist_plants <- tep |> dplyr::filter(climate_model == "historical") |> dplyr::distinct(fips, plant_id) |>
    dplyr::count(fips, name = "hist_plants")
  future_plants <- tep |> dplyr::filter(climate_model != "historical") |> dplyr::distinct(fips, plant_id) |>
    dplyr::count(fips, name = "future_plants")

  # 3. Ensemble members
  members <- pub |> dplyr::filter(climate_model != "historical", emission_scenario != "ssp126") |>
    dplyr::distinct(climate_model, emission_scenario)
  if (!is.null(emission_scenarios)) members <- dplyr::filter(members, emission_scenario %in% emission_scenarios)
  if (!is.null(climate_models)) members <- dplyr::filter(members, climate_model %in% climate_models)

  # 4. County withdrawals per member, averaged over months (weighted by days)
  totals <- tidyr::expand_grid(fips = unique(stats::na.omit(locations$fips)), Month = months) |>
    dplyr::cross_join(dplyr::bind_rows(dplyr::tibble(climate_model = "historical", emission_scenario = "historical"),
                         members)) |>
    dplyr::left_join(pub, by = c("fips", "Month", model_cols)) |>
    dplyr::left_join(tep |> dplyr::group_by(fips, Month, climate_model, emission_scenario) |>
                dplyr::summarise(tep_mgd = sum(tep), n_plants = dplyr::n(), .groups = "drop"),
              by = c("fips", "Month", model_cols)) |>
    dplyr::left_join(hist_plants, by = "fips") |>
    dplyr::left_join(future_plants, by = "fips") |>
    dplyr::mutate(n_plants = dplyr::coalesce(n_plants, 0L),
           n_expected = dplyr::coalesce(dplyr::if_else(climate_model == "historical", hist_plants, future_plants), 0L),
           # No plants expected this period = zero; a plant that should report
           # this period but doesn't (a month or member gap) = NA
           tep_mgd = dplyr::case_when(n_expected == 0 ~ 0, n_plants < n_expected ~ NA, TRUE ~ tep_mgd),
           total_mgd = pub_irr_mgd + tep_mgd) |>
    dplyr::group_by(fips, climate_model, emission_scenario) |>
    dplyr::summarise(withdrawal_mgd = stats::weighted.mean(total_mgd, days_in_month[Month]),
              n_public_months = sum(!is.na(in_pub_file)),
              n_complete_months = sum(is.finite(total_mgd)),
              missing_plant_records = sum(n_expected - n_plants),
              .groups = "drop")

  historical <- totals |>
    dplyr::filter(climate_model == "historical") |>
    dplyr::transmute(fips, historic_withdrawal_mgd = withdrawal_mgd,
              historical_complete_months = n_complete_months,
              historical_public_months = n_public_months)

  # 5. Change from historical, per member
  projected <- totals |>
    dplyr::filter(climate_model != "historical") |>
    dplyr::left_join(historical, by = "fips") |>
    dplyr::left_join(dplyr::distinct(locations, fips, county_area_km2), by = "fips") |>
    dplyr::mutate(absolute_change = (withdrawal_mgd - historic_withdrawal_mgd) / county_area_km2,
           percent_change = dplyr::if_else(historic_withdrawal_mgd > 0,
                                    100 * (withdrawal_mgd - historic_withdrawal_mgd) / historic_withdrawal_mgd,
                                    NA_real_))

  # 6. Summarize the ensemble per county
  n_months <- length(months)
  result <- projected |>
    dplyr::group_by(fips) |>
    dplyr::summarise(n_valid_models = sum(is.finite(absolute_change)),
              n_model_members = nrow(members),
              raw_models_showing_increase_pct = if (n_valid_models > 0) {
                100 * sum(absolute_change[is.finite(absolute_change)] > 0) / n_valid_models
              } else NA_real_,
              missing_future_plant_records = sum(missing_plant_records),
              future_public_months = max(n_public_months),
              future_complete_months = max(n_complete_months),
              .groups = "drop") |>
    dplyr::left_join(historical, by = "fips") |>
    dplyr::mutate(data_status = dplyr::case_when(
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
      dplyr::filter(is.finite(absolute_change)) |>
      dplyr::group_by(fips) |>
      dplyr::slice_min(abs(absolute_change - stats::quantile(absolute_change, p / 100)), n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::transmute(fips,
                !!paste0("raw_use_change_p", p) := percent_change,
                !!paste0("raw_absolute_use_change_p", p) := absolute_change,
                !!paste0("future_withdrawal_mgd_p", p) := withdrawal_mgd)
    result <- dplyr::left_join(result, nearest_member, by = "fips")
  }

  # 7. Join to supplies and standardize
  final <- locations |>
    dplyr::left_join(result, by = "fips") |>
    dplyr::mutate(data_status = dplyr::coalesce(data_status, "unmatched_county")) |>
    add_standardized_indicators("raw_")

  if (save) {
    filename <- file.path(output_dir, paste0("nearby_use_exposure_", Sys.Date(), ".csv"))
    readr::write_csv(final, filename)
    if (verbose) message("Results saved to: ", filename)
  }
  final
}
