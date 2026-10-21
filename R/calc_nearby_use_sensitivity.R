#' @title Nearby use sensitivity: historical county water withdrawals
#'
#' @description
#' Gives each water supply the 2000-2020 water withdrawals of its county
#' (public supply, crop irrigation and thermoelectric, in million gallons per
#' day, MGD), from the Blumberg and Warziniack county data, per km2 of county
#' land. The main measure is `raw_avg_nearby_use_per_km2`. The data are 12
#' average months, not a series of years, so no trend is calculated.
#'
#' Each water supply gets the 2020 TIGER/Line county it is in (downloaded
#' with [tigris::counties()], so it needs internet), or the nearest county
#' within `max_county_distance_km`. The months in `months` are averaged,
#' weighted by days. A missing month, or a known power plant that doesn't
#' report, makes the county total `NA`; gaps stay `NA` rather than being read
#' as zero use. The `norm_` and `rank_` columns are relative to the water
#' supplies passed in.
#'
#' @param water_supplies sf object. The water supplies, as points with
#'   `wsd_source_id`, `park_unit` and `park_name` columns, e.g. from
#'   [wavacore::get_wsd()]. Required.
#' @param data_dir String. The folder with `pub_irr_monthly_avg_2000_2020.csv`
#'   and `tep_monthly_avg_2000_2020.csv`. Default:
#'   `Latest county-level data/water_demand_projections` in the shared data
#'   folder, see [wavacore::wava_data_path()].
#' @param save Logical. Whether to write the results to a CSV file.
#'   * `TRUE`: writes `nearby_use_sensitivity_<date>.csv`, named with today's
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
#' @param months Integer vector. The calendar months (1-12) to average.
#'   Default: `1:12`.
#' @param max_county_distance_km Numeric. How far, in km, a supply outside
#'   every county can be from the nearest county and still use it. Supplies
#'   farther away get `NA`. Default: `5`.
#'
#' @return A tibble, one row per water supply, with `wsd_source_id`,
#'   `park_unit`, `park_name`, the county (`fips`, `county_area_km2`,
#'   `county_distance_km`, `county_source`), `raw_avg_nearby_use` (MGD),
#'   `raw_avg_nearby_use_per_km2` (MGD per km2), and their `norm_` and
#'   `rank_` versions from [add_standardized_indicators()]. There is no
#'   `data_status` column; every gap here is historical.
#' @export
#' @seealso [calc_nearby_use_exposure()] for the projected change;
#'   [calc_wta_sensitivity()] for watershed water stress.
#' @examplesIf interactive() && dir.exists(wavacore::wava_data_path("Latest county-level data", "water_demand_projections"))
#' # Downloads the county boundaries, so it needs internet
#' ws <- wavacore::get_wsd(park = c("ZION", "ROMO"))
#' use_sen <- calc_nearby_use_sensitivity(ws)
#' use_sen[, c("wsd_source_id", "fips", "raw_avg_nearby_use_per_km2")]

calc_nearby_use_sensitivity <- function(water_supplies,
                                        data_dir = wavacore::wava_data_path("Latest county-level data", "water_demand_projections"),
                                        save = FALSE,
                                        output_dir = NULL,
                                        verbose = TRUE,
                                        months = 1:12,
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

  # 2. Read historical withdrawals
  pub <- data.table::fread(file.path(data_dir, "pub_irr_monthly_avg_2000_2020.csv"),
               select = c("FIPS", "Month", "population", "pub_per_cap", "irr")) |>
    dplyr::transmute(fips = sprintf("%05d", FIPS), Month,
              pub_irr_mgd = pub_per_cap * population / 1000 + irr)  # population is thousands
  tep <- data.table::fread(file.path(data_dir, "tep_monthly_avg_2000_2020.csv"),
               select = c("FIPS", "plant_id", "Month", "tep")) |>
    dplyr::mutate(fips = sprintf("%05d", FIPS))
  plants <- tep |> dplyr::distinct(fips, plant_id) |> dplyr::count(fips, name = "n_expected")

  # 3. County withdrawals, averaged over months (weighted by days)
  result <- tidyr::expand_grid(fips = unique(stats::na.omit(locations$fips)), Month = months) |>
    dplyr::left_join(pub, by = c("fips", "Month")) |>
    dplyr::left_join(tep |> dplyr::group_by(fips, Month) |>
                dplyr::summarise(tep_mgd = sum(tep), n_plants = dplyr::n(), .groups = "drop"),
              by = c("fips", "Month")) |>
    dplyr::left_join(plants, by = "fips") |>
    dplyr::mutate(n_plants = dplyr::coalesce(n_plants, 0L),
           n_expected = dplyr::coalesce(n_expected, 0L),
           # No known plants = zero; a missing known plant = NA
           tep_mgd = dplyr::case_when(n_expected == 0 ~ 0, n_plants < n_expected ~ NA, TRUE ~ tep_mgd),
           total_mgd = pub_irr_mgd + tep_mgd) |>
    dplyr::group_by(fips) |>
    dplyr::summarise(raw_avg_nearby_use = stats::weighted.mean(total_mgd, days_in_month[Month])) |>
    dplyr::left_join(dplyr::distinct(locations, fips, county_area_km2), by = "fips") |>
    # if NA, it is because the county data is missing or incomplete
    dplyr::transmute(fips, raw_avg_nearby_use,
              raw_avg_nearby_use_per_km2 = raw_avg_nearby_use / county_area_km2)

  # 4. Join to supplies and standardize
  final <- locations |>
    dplyr::left_join(result, by = "fips") |>
    add_standardized_indicators("raw_")

  if (save) {
    filename <- file.path(output_dir, paste0("nearby_use_sensitivity_", Sys.Date(), ".csv"))
    readr::write_csv(final, filename)
    if (verbose) message("Results saved to: ", filename)
  }
  final
}
