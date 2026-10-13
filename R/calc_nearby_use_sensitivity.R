#' Calculate historical nearby-use pressure from county monthly climatologies
#'
#' 2000-2020 county withdrawals (public supply + crop irrigation +
#' thermoelectric, in MGD) from Blumberg and Warziniack, per county area. The
#' main measure is raw_avg_nearby_use_per_km2. The data are 12 climatological
#' months, not a yearly series, so no trend can be estimated. A missing month
#' or known power plant makes a county total NA.
#'
#' Supplies use the 2020 TIGER/Line county they fall in, or the nearest county
#' within max_county_distance_km. norm_/rank_ columns are across supplies.
#'
#' @param water_supplies sf points with wsd_source_id, park_unit, park_name.
#' @param data_dir Folder with the water_demand_projections CSVs.
#' @param save,output_dir Save nearby_use_sensitivity_YYYY-MM-DD.csv to output_dir.
#' @param verbose Report the saved filename.
#' @param months Calendar months to average (default all 12).
#' @param max_county_distance_km Nearest-county limit for points outside counties.
#' @return Tibble with one row per supply. Results are NA when no county falls
#'   within max_county_distance_km, the county is absent from the historical
#'   file, or any month or expected plant is unreported. There is no
#'   data_status column, as every gap here is a historical one. Gaps are left
#'   NA rather than treated as zero use.
calc_nearby_use_sensitivity <- function(water_supplies,
                                        data_dir = "data/Latest county-level data/water_demand_projections",
                                        save = FALSE,
                                        output_dir = NULL,
                                        verbose = TRUE,
                                        months = 1:12,
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

  # 2. Read historical withdrawals
  pub <- fread(file.path(data_dir, "pub_irr_monthly_avg_2000_2020.csv"),
               select = c("FIPS", "Month", "population", "pub_per_cap", "irr")) |>
    transmute(fips = sprintf("%05d", FIPS), Month,
              pub_irr_mgd = pub_per_cap * population / 1000 + irr)  # population is thousands
  tep <- fread(file.path(data_dir, "tep_monthly_avg_2000_2020.csv"),
               select = c("FIPS", "plant_id", "Month", "tep")) |>
    mutate(fips = sprintf("%05d", FIPS))
  plants <- tep |> distinct(fips, plant_id) |> count(fips, name = "n_expected")

  # 3. County withdrawals, averaged over months (weighted by days)
  result <- expand_grid(fips = unique(na.omit(locations$fips)), Month = months) |>
    left_join(pub, by = c("fips", "Month")) |>
    left_join(tep |> group_by(fips, Month) |>
                summarise(tep_mgd = sum(tep), n_plants = n(), .groups = "drop"),
              by = c("fips", "Month")) |>
    left_join(plants, by = "fips") |>
    mutate(n_plants = coalesce(n_plants, 0L),
           n_expected = coalesce(n_expected, 0L),
           # No known plants = zero; a missing known plant = NA
           tep_mgd = case_when(n_expected == 0 ~ 0, n_plants < n_expected ~ NA, TRUE ~ tep_mgd),
           total_mgd = pub_irr_mgd + tep_mgd) |>
    group_by(fips) |>
    summarise(raw_avg_nearby_use = weighted.mean(total_mgd, days_in_month[Month])) |>
    left_join(distinct(locations, fips, county_area_km2), by = "fips") |>
    # if NA, it is because the county data is missing or incomplete
    transmute(fips, raw_avg_nearby_use,
              raw_avg_nearby_use_per_km2 = raw_avg_nearby_use / county_area_km2)

  # 4. Join to supplies and standardize
  final <- locations |>
    left_join(result, by = "fips") |>
    add_standardized_indicators("raw_")

  if (save) {
    filename <- file.path(output_dir, paste0("nearby_use_sensitivity_", Sys.Date(), ".csv"))
    write_csv(final, filename)
    if (verbose) message("Results saved to: ", filename)
  }
  final
}
