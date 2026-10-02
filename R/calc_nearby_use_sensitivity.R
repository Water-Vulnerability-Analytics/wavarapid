#' Calculate historical nearby-use pressure from county monthly climatologies
#'
#' Measures 2000-2020 county withdrawal intensity from the Blumberg and
#' Warziniack public supply, crop irrigation, and thermoelectric cooling data.
#' The main sensitivity measure is raw_avg_nearby_use_per_km2 (MGD/km2); higher
#' values represent greater nearby demand pressure. County-wide withdrawals are
#' a proxy for nearby competition, not a direct measure of a supply's scarcity.
#'
#' The files contain 12 climatological months, NOT a year-by-year series, so
#' no historical trend (Sen slope, Mann-Kendall) can be estimated.
#'
#' Public supply MGD = pub_per_cap * population / 1000 (population is in
#' thousands). Plant withdrawals are summed by county and month. Monthly daily
#' rates are averaged weighted by month length (February = 28.25 days). Every
#' selected month and every known plant is required; missing data are not zero.
#' Counties with no listed plants contribute zero thermoelectric use.
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
#' @param save Whether to save nearby_use_sensitivity_YYYY-MM-DD.csv.
#' @param output_dir Required when save is TRUE.
#' @param verbose Whether to report the saved filename.
#' @param months Calendar months to average (default all 12).
#' @param max_county_distance_km Maximum distance to the nearest county for
#'   points outside every county polygon (default 5 km). Inf allows any distance.
#' @return Nonspatial tibble, one row per supply: wsd_source_id, park_unit,
#'   park_name (if supplied), fips, county_area_km2, county_distance_km,
#'   county_source, raw_avg_nearby_use (MGD), raw_avg_nearby_use_per_km2
#'   (MGD/km2), and their norm_/rank_ columns. Supplies in counties with
#'   incomplete data have NA values.
calc_nearby_use_sensitivity <- function(water_supplies,
                                        data_dir = file.path("data", "Latest county-level data"),
                                        save = FALSE,
                                        output_dir = NULL,
                                        verbose = TRUE,
                                        months = 1:12,
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

  # 2. Read historical withdrawals -------------------------------------------
  pub_columns <- c("population", "pub_per_cap", "irr")
  pub <- data.table::fread(
    file.path(data_dir, "pub_irr_monthly_avg_2000_2020.csv"),
    select = c("FIPS", "Month", pub_columns),
    colClasses = list(character = "FIPS", numeric = pub_columns),
    data.table = FALSE, showProgress = FALSE)
  tep <- data.table::fread(
    file.path(data_dir, "tep_monthly_avg_2000_2020.csv"),
    select = c("FIPS", "plant_id", "Month", "tep"),
    colClasses = list(character = c("FIPS", "plant_id"), numeric = "tep"),
    data.table = FALSE, showProgress = FALSE)

  # Five-character FIPS text, matching the county GEOIDs.
  pub$fips <- sprintf("%05d", as.integer(pub$FIPS))
  tep$fips <- sprintf("%05d", as.integer(tep$FIPS))

  # 3. County withdrawals by month, then averaged ----------------------------
  plants <- dplyr::count(dplyr::distinct(tep, fips, plant_id), fips,
                         name = "n_expected_plants")
  plant_totals <- tep |>
    dplyr::group_by(fips, Month) |>
    dplyr::summarise(thermoelectric_mgd = sum(tep), n_plant_records = dplyr::n(),
                     .groups = "drop")

  result <- tidyr::expand_grid(fips = supply_fips, Month = months) |>
    dplyr::left_join(pub |>
                       dplyr::transmute(fips, Month,
                                        public_supply_mgd = pub_per_cap * population / 1000,
                                        irrigation_mgd = irr),
                     by = c("fips", "Month")) |>
    dplyr::left_join(plant_totals, by = c("fips", "Month")) |>
    dplyr::left_join(plants, by = "fips") |>
    dplyr::mutate(
      n_expected_plants = dplyr::coalesce(n_expected_plants, 0L),
      n_plant_records = dplyr::coalesce(n_plant_records, 0L),
      # No known plants is a true zero; a missing known plant is missing data.
      thermoelectric_mgd = dplyr::case_when(
        n_expected_plants == 0L ~ 0,
        n_plant_records < n_expected_plants ~ NA_real_,
        TRUE ~ thermoelectric_mgd),
      total_mgd = public_supply_mgd + irrigation_mgd + thermoelectric_mgd,
      days = days_in_month[Month]) |>
    dplyr::group_by(fips) |>
    dplyr::summarise(raw_avg_nearby_use = stats::weighted.mean(total_mgd, days),
                     .groups = "drop") |>
    dplyr::left_join(dplyr::distinct(locations, fips, county_area_km2), by = "fips") |>
    dplyr::transmute(fips, raw_avg_nearby_use,
                     raw_avg_nearby_use_per_km2 = raw_avg_nearby_use / county_area_km2)

  # 4. Join to supplies and standardize across supplies ----------------------
  final <- locations |>
    dplyr::left_join(result, by = "fips") |>
    add_standardized_indicators("raw_", invert = FALSE)

  if (save) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    filename <- file.path(output_dir, paste0("nearby_use_sensitivity_", Sys.Date(), ".csv"))
    readr::write_csv(final, filename)
    if (verbose) message("Results saved to: ", filename)
  }

  final
}
