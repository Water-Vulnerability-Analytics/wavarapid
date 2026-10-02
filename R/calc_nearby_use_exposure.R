#' Calculate nearby-use exposure from the new county demand projections
#'
#' Compare 2040-2070 with 2000-2020 monthly climatologies supplied by Blumberg
#' and Warziniack. County withdrawals combine public supply, crop irrigation,
#' and plant-level thermoelectric cooling. Values are millions of gallons/day
#' (MGD); public supply is pub_per_cap * population / 1000. Monthly daily rates
#' are averaged using climatological days per month (February = 28.25 days).
#'
#' Each climate_model/emission_scenario pair is one equally weighted ensemble
#' member. SSP1-2.6 (ssp126) is always excluded, leaving 7 models x 3 pathways
#' (21 members) in the supplied bundle. For each p10/p50/p90, output the actual
#' member whose absolute change per county area is nearest the type-7 percentile
#' target; percent change and future withdrawal come from that same member.
#' Absolute changes are reported in MGD/km2; percent change is
#' 100 * (future - historical) / historical MGD.
#' Model agreement is the share of complete pairs with positive absolute change.
#'
#' Known plants absent from the future file are missing, not assumed retired.
#' Counties with no plants in either file contribute zero thermoelectric use.
#' All selected months and all known plants are required for a complete member.
#' Missing members are excluded from quantiles/agreement and reported in counts.
#' No FIPS substitutions are made for missing counties or changed boundaries.
#'
#' County boundaries always use tigris::counties(cb = FALSE, year = 2020).
#' Areas include county land and water territory. Each supply uses the county
#' containing it, or the nearest county within max_county_distance_km.
#'
#' norm_ (0-1 min-max) and rank_ (0-100 cumulative percentile) columns are
#' calculated across water supplies, so a county with several supplies is
#' counted once per supply.
#'
#' Source this file directly, or use setup.R.
#'
#' @param water_supplies sf POINT locations with unique wsd_source_id, park_unit,
#'   and optionally park_name. Empty points are retained with NA indicators.
#' @param data_dir Directory containing water_demand_projections, or that folder
#'   itself. Requires the pub_irr and tep CSVs for both supplied periods.
#' @param save Whether to save nearby_use_exposure_YYYY-MM-DD.csv.
#' @param output_dir Required when save is TRUE.
#' @param verbose Whether to report the saved filename.
#' @param months Unique calendar months to average (default all 12). Use the
#'   same months for sensitivity when comparing seasonal indicators.
#' @param emission_scenarios Optional exact pathway names, e.g. "ssp245".
#'   ssp126 is not available.
#' @param climate_models Optional exact climate model names. NULL uses all.
#' @param max_county_distance_km Maximum distance to the nearest county for
#'   points outside every county polygon (default 5 km). Inf allows any distance.
#' @return Tibble, one row per supply: county fips (five-character text), area
#'   and distance; raw_use_change_p10/p50/p90 (percent change),
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
  
  # Helpers:
  format_fips <- function(x, label) {
    x <- trimws(as.character(x))
    if (anyNA(x) || any(!grepl("^[0-9]{1,5}$", x))) {
      stop(label, " must contain nonmissing county FIPS codes (1-5 digits).", call. = FALSE)
    }
    sprintf("%05d", as.integer(x))
  }
  
  validate_options <- function(months, save, output_dir, verbose) {
    if (!is.numeric(months) || !length(months) || anyNA(months) ||
        any(!months %in% 1:12) || anyDuplicated(months)) {
      stop("months must be unique integers from 1 through 12.", call. = FALSE)
    }
    for (x in list(save, verbose)) {
      if (!is.logical(x) || length(x) != 1L || is.na(x)) {
        stop("save and verbose must each be TRUE or FALSE.", call. = FALSE)
      }
    }
    if (save && (is.null(output_dir) || !is.character(output_dir) ||
                 length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir))) {
      stop("save is TRUE but output_dir was not provided.", call. = FALSE)
    }
    sort(as.integer(months))
  }
  
  assign_counties <- function(water_supplies, max_county_distance_km) {
    if (!inherits(water_supplies, "sf") || is.na(sf::st_crs(water_supplies))) {
      stop("water_supplies must be an sf object with a known CRS.", call. = FALSE)
    }
    required <- c("wsd_source_id", "park_unit")
    if (!all(required %in% names(water_supplies)) || !nrow(water_supplies)) {
      stop("water_supplies must have rows and columns wsd_source_id and park_unit.", call. = FALSE)
    }
    if (anyNA(water_supplies$wsd_source_id) || anyDuplicated(water_supplies$wsd_source_id)) {
      stop("wsd_source_id must be nonmissing and unique.", call. = FALSE)
    }
    if (any(as.character(sf::st_geometry_type(water_supplies)) != "POINT")) {
      stop("water_supplies must contain POINT geometries.", call. = FALSE)
    }
    if (!is.numeric(max_county_distance_km) || length(max_county_distance_km) != 1L ||
        is.na(max_county_distance_km) || max_county_distance_km < 0) {
      stop("max_county_distance_km must be nonnegative (Inf allows any distance).", call. = FALSE)
    }
    counties <- tigris::counties(cb = FALSE, year = 2020, progress_bar = FALSE)
    county_source <- "2020 Census TIGER/Line counties"
    if (!inherits(counties, "sf") || !nrow(counties) || is.na(sf::st_crs(counties)) ||
        any(!as.character(sf::st_geometry_type(counties)) %in% c("POLYGON", "MULTIPOLYGON")) ||
        any(sf::st_is_empty(counties))) {
      stop("counties must be nonempty sf county polygons with a known CRS.", call. = FALSE)
    }
    id_col <- intersect(c("GEOID", "fips", "FIPS"), names(counties))
    if (length(id_col)) {
      fips <- format_fips(counties[[id_col[1L]]], "counties")
    } else if (all(c("STATEFP", "COUNTYFP") %in% names(counties))) {
      fips <- format_fips(sprintf("%02d%03d", as.integer(counties$STATEFP),
                                  as.integer(counties$COUNTYFP)), "counties")
    } else {
      stop("counties needs GEOID, fips, FIPS, or STATEFP and COUNTYFP.", call. = FALSE)
    }
    if (anyDuplicated(fips)) stop("counties must have one polygon row per FIPS.", call. = FALSE)
    counties <- sf::st_make_valid(counties)
    area <- as.numeric(sf::st_area(sf::st_transform(counties, 5070))) / 1e6
    if (any(!is.finite(area) | area <= 0)) stop("County areas must be positive.", call. = FALSE)
    # Work in lon/lat so distances returned by sf are geodesic metres.
    counties <- sf::st_transform(counties, 4326)
    points <- sf::st_transform(water_supplies, 4326)
    valid <- which(!sf::st_is_empty(points))
    index <- rep(NA_integer_, nrow(points))
    distance <- rep(NA_real_, nrow(points))
    if (length(valid)) {
      hits <- sf::st_intersects(points[valid, ], counties, model = "closed")
      # A point on a shared boundary gets one deterministic county, never two rows.
      index[valid] <- vapply(hits, function(i) {
        if (length(i)) i[order(fips[i])][1L] else NA_integer_
      }, integer(1))
      distance[!is.na(index)] <- 0
      # Points outside every county use the nearest county within the limit.
      missing <- valid[is.na(index[valid])]
      if (length(missing)) {
        nearest <- sf::st_nearest_feature(points[missing, ], counties)
        km <- as.numeric(sf::st_distance(points[missing, ], counties[nearest, ],
                                         by_element = TRUE)) / 1000
        distance[missing] <- km
        accept <- is.finite(km) & km <= max_county_distance_km
        index[missing[accept]] <- nearest[accept]
      }
    }
    result <- sf::st_drop_geometry(water_supplies)[, intersect(
      c("wsd_source_id", "park_unit", "park_name"), names(water_supplies)), drop = FALSE]
    result$fips <- fips[index]
    result$county_area_km2 <- area[index]
    result$county_distance_km <- distance
    result$county_source <- county_source
    tibble::as_tibble(result)
  }
  
  read_withdrawals <- function(data_dir, sector, future, fips) {
    period <- if (future) "2040_2070" else "2000_2020"
    filename <- paste0(sector, "_monthly_avg_", period, ".csv")
    path <- file.path(data_dir, filename)
    if (!file.exists(path)) path <- file.path(data_dir, "water_demand_projections", filename)
    if (!file.exists(path)) stop("Missing data file: ", filename, " in ", data_dir, call. = FALSE)
    keys <- c("FIPS", if (sector == "tep") "plant_id", "Month",
              if (future) c("climate_model", "emission_scenario"))
    values <- if (sector == "tep") "tep" else c("population", "pub_per_cap", "irr")
    header <- names(data.table::fread(path, nrows = 0L, showProgress = FALSE))
    if (!all(c(keys, values) %in% header)) {
      stop(filename, " is missing columns: ", paste(setdiff(c(keys, values), header), collapse = ", "),
           call. = FALSE)
    }
    x <- data.table::fread(path, select = c(keys, values),
                           colClasses = list(character = setdiff(keys, "Month"), numeric = values),
                           data.table = FALSE, showProgress = FALSE)
    x$FIPS <- format_fips(x$FIPS, filename)
    # Capture the ensemble before restricting to supply counties (which may lack data).
    members <- if (future) unique(x[c("climate_model", "emission_scenario")]) else NULL
    if (future && (anyNA(members) || any(!nzchar(as.matrix(members))))) {
      stop(filename, " contains missing model/scenario identifiers.", call. = FALSE)
    }
    x <- x[x$FIPS %in% fips, , drop = FALSE]
    if (anyNA(x[keys]) || any(!x$Month %in% 1:12) || anyDuplicated(x[keys])) {
      stop(filename, " contains invalid months, missing keys, or duplicate records.", call. = FALSE)
    }
    for (column in values) {
      if (!is.numeric(x[[column]]) || any(!is.na(x[[column]]) &
                                          (!is.finite(x[[column]]) | x[[column]] < 0))) {
        stop(filename, ": ", column, " must be nonnegative numeric values or NA.", call. = FALSE)
      }
    }
    names(x)[names(x) == "FIPS"] <- "fips"
    list(data = tibble::as_tibble(x), members = tibble::as_tibble(members))
  }
  
  aggregate_monthly <- function(pub, tep, fips, months, members = NULL,
                                historical_plants = NULL) {
    keys <- c("fips", "Month", if (!is.null(members)) c("climate_model", "emission_scenario"))
    grid <- tidyr::expand_grid(fips = fips, Month = months)
    if (!is.null(members)) grid <- dplyr::cross_join(grid, members)
    # population is THOUSANDS of people; pub_per_cap is gallons/person/day.
    pub <- dplyr::mutate(pub, public_supply_mgd = pub_per_cap * population / 1000,
                         irrigation_mgd = irr, public_record = TRUE)
    inventory <- dplyr::distinct(tep, fips, plant_id)
    if (!is.null(historical_plants)) {
      inventory <- dplyr::distinct(dplyr::bind_rows(inventory, historical_plants), fips, plant_id)
    }
    expected <- dplyr::count(inventory, fips, name = "n_expected_plants")
    plant_totals <- tep |>
      dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
      dplyr::summarise(thermoelectric_mgd = sum(tep),
                       n_plant_records = dplyr::n(), .groups = "drop")
    result <- grid |>
      dplyr::left_join(dplyr::select(pub, dplyr::all_of(keys), public_supply_mgd,
                                     irrigation_mgd, public_record), by = keys) |>
      dplyr::left_join(plant_totals, by = keys) |>
      dplyr::left_join(expected, by = "fips") |>
      dplyr::mutate(
        public_record = dplyr::coalesce(public_record, FALSE),
        n_expected_plants = dplyr::coalesce(n_expected_plants, 0L),
        n_plant_records = dplyr::coalesce(n_plant_records, 0L),
        # No plant in the inventory is a structural zero. An absent known plant
        # is missing data, even if other plants in the county have projections.
        thermoelectric_mgd = dplyr::case_when(
          n_expected_plants == 0L ~ 0,
          n_plant_records < n_expected_plants ~ NA_real_,
          TRUE ~ thermoelectric_mgd),
        missing_plant_records = n_expected_plants - n_plant_records,
        total_mgd = public_supply_mgd + irrigation_mgd + thermoelectric_mgd)
    result
  }
  
  average_period <- function(monthly, future = FALSE) {
    # Use identical climatological month lengths in both periods. Do not sum
    # rates in MGD or interpret these 12 climatological months as yearly data.
    days <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    monthly$days <- days[monthly$Month]
    monthly |>
      dplyr::group_by(dplyr::across(dplyr::all_of(
        c("fips", if (future) c("climate_model", "emission_scenario"))))) |>
      dplyr::summarise(
        withdrawal_mgd = stats::weighted.mean(total_mgd, days),
        public_supply_mgd = stats::weighted.mean(public_supply_mgd, days),
        irrigation_mgd = stats::weighted.mean(irrigation_mgd, days),
        thermoelectric_mgd = stats::weighted.mean(thermoelectric_mgd, days),
        peak_month_withdrawal_mgd = if (length(total_mgd) && all(is.finite(total_mgd))) {
          max(total_mgd)
        } else NA_real_,
        n_public_months = sum(public_record),
        n_complete_months = sum(is.finite(total_mgd)),
        missing_plant_records = sum(missing_plant_records),
        .groups = "drop")
  }
  
  select_percentile_member <- function(data, probability, suffix) {
    data |>
      dplyr::filter(is.finite(absolute_change)) |>
      dplyr::group_by(fips) |>
      dplyr::mutate(target = stats::quantile(absolute_change, probability,
                                              names = FALSE, type = 7),
                    distance = abs(absolute_change - target)) |>
      dplyr::slice_min(distance, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::transmute(fips,
        !!paste0("raw_use_change_", suffix) := percent_change,
        !!paste0("raw_absolute_use_change_", suffix) := absolute_change,
        !!paste0("future_withdrawal_mgd_", suffix) := withdrawal_mgd)
  }
  
  standardize <- function(data, columns) {
    # Match helper_functions.R: min-max normalization and cumulative percentile
    # rank, calculated across water supplies (one row per supply).
    for (column in columns) {
      x <- data[[column]]
      data[[sub("^raw_", "norm_", column)]] <- if (all(is.na(x))) {
        rep(NA_real_, length(x))
      } else {
        scales::rescale(x, to = c(0, 1))
      }
      data[[sub("^raw_", "rank_", column)]] <- dplyr::cume_dist(x) * 100
    }
    data
  }
  
  finish_results <- function(locations, result, save, output_dir, verbose, kind) {
    # A left join preserves supply order, including supplies outside data coverage.
    final <- dplyr::left_join(locations, result, by = "fips") |>
      dplyr::mutate(data_status = dplyr::coalesce(data_status, "unmatched_county"))
    final <- standardize(final, names(final)[startsWith(names(final), "raw_")])
    if (save) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      filename <- file.path(output_dir, paste0("nearby_use_", kind, "_", Sys.Date(), ".csv"))
      readr::write_csv(final, filename)
      if (verbose) message("Results saved to: ", filename)
    }
    final
  }
  
  # Calculate the indicator.
  months <- validate_options(months, save, output_dir, verbose)
  locations <- assign_counties(water_supplies, max_county_distance_km)
  fips <- unique(stats::na.omit(locations$fips))
  hist_pub <- read_withdrawals(data_dir, "pub_irr", FALSE, fips)$data
  hist_tep <- read_withdrawals(data_dir, "tep", FALSE, fips)$data
  future_pub <- read_withdrawals(data_dir, "pub_irr", TRUE, fips)
  future_tep <- read_withdrawals(data_dir, "tep", TRUE, fips)$data
  # SSP1-2.6 is always excluded from the ensemble.
  members <- future_pub$members[future_pub$members$emission_scenario != "ssp126", , drop = FALSE]
  filters <- list(emission_scenario = emission_scenarios, climate_model = climate_models)
  for (column in names(filters)) {
    requested <- filters[[column]]
    if (!is.null(requested)) {
      if (!is.character(requested) || !length(requested) || anyNA(requested) ||
          any(!requested %in% members[[column]])) {
        stop("Unknown or empty ", column, " selection (ssp126 is excluded). Available: ",
             paste(unique(members[[column]]), collapse = ", "), call. = FALSE)
      }
      members <- members[members[[column]] %in% requested, , drop = FALSE]
    }
  }
  if (!nrow(members)) stop("No model/scenario pairs match the selection.", call. = FALSE)
  historical <- aggregate_monthly(hist_pub, hist_tep, fips, months) |>
    average_period() |>
    dplyr::transmute(fips, historic_withdrawal_mgd = withdrawal_mgd,
                     historical_complete_months = n_complete_months,
                     historical_public_months = n_public_months)
  projected <- aggregate_monthly(future_pub$data, future_tep, fips, months, members,
                                 historical_plants = dplyr::distinct(hist_tep, fips, plant_id)) |>
    average_period(future = TRUE) |>
    dplyr::left_join(historical, by = "fips") |>
    dplyr::left_join(dplyr::distinct(dplyr::select(locations, fips, county_area_km2)), by = "fips") |>
    dplyr::mutate(absolute_change = (withdrawal_mgd - historic_withdrawal_mgd) / county_area_km2,
      percent_change = dplyr::if_else(historic_withdrawal_mgd > 0,
        100 * (withdrawal_mgd - historic_withdrawal_mgd) / historic_withdrawal_mgd, NA_real_))

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
    dplyr::mutate(
      data_status = dplyr::case_when(
        # Any missing month or plant record makes the historical mean NA, so
        # check whether the county is absent before checking completeness.
        historical_public_months == 0 ~ "missing_historical_county",
        !is.finite(historic_withdrawal_mgd) |
          historical_public_months < length(months) |
          historical_complete_months < length(months) ~ "incomplete_historical_data",
        future_public_months == 0 ~ "missing_future_county",
        missing_future_plant_records > 0 |
          future_public_months < length(months) |
          future_complete_months < length(months) ~ "incomplete_future_data",
        n_valid_models == 0 ~ "no_complete_model_members",
        n_valid_models < n_model_members ~ "partial_ensemble",
        historic_withdrawal_mgd == 0 ~ "zero_historical_withdrawal",
        TRUE ~ "ok"))
  probabilities <- c(p10 = 0.1, p50 = 0.5, p90 = 0.9)
  for (suffix in names(probabilities)) {
    result <- dplyr::left_join(result,
      select_percentile_member(projected, probabilities[[suffix]], suffix), by = "fips")
  }

  finish_results(locations, result, save, output_dir, verbose, "exposure")
}
