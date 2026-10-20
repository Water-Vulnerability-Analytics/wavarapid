#' Calculate Consecutive Dry-Day Exposure Indicator
#'
#' Calculates projected change in annual maximum consecutive dry-day spell
#' length (CDD) from the processed LOCA2 precipitation-extremes dataset. For
#' each model/SSP/run, exposure is the percent change from historical mean CDD
#' to future mean CDD. Absolute changes in days are retained as supporting
#' values. The closest actual model values to the p10, p50, and p90 of the
#' local ensemble are retained, along with the percentage of models showing
#' longer dry spells.
#'
#' Larger positive percent changes indicate greater exposure to extended dry
#' periods.
#'
#' @param water_supplies An `sf` object containing `wsd_source_id` and
#'   `park_unit`.
#' @param period_data Character. Directory containing parquet files created by
#'   `process_precip_period_data()`.
#' @param historical Character vector of length 2 giving the historical
#'   start/end years.
#' @param future Character vector of length 2 giving the future start/end years.
#' @param save Logical. If `TRUE`, writes results to disk.
#' @param output_dir Character. Output directory when `save = TRUE`.
#' @param verbose Logical. If `TRUE`, prints progress messages.
#'
#' @return A tibble with one row per water supply and dry-period exposure
#'   indicators.
#' @export
calc_dry_period_exposure <- function(
    water_supplies,
    period_data = wava_data_path("rapid", "inputs", "loca_cmip6_precip_periods_water_supplies/"),
    historical,
    future,
    save = FALSE,
    output_dir = NULL,
    verbose = TRUE) {

  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "water_supplies must include wsd_source_id" =
      "wsd_source_id" %in% names(water_supplies),
    "water_supplies must include park_unit" =
      "park_unit" %in% names(water_supplies),
    "historical must be length 2" = length(historical) == 2,
    "future must be length 2" = length(future) == 2,
    "save is TRUE but output_dir not provided" =
      !save || !is.null(output_dir)
  )

  if (!dir.exists(period_data)) {
    stop("Processed precipitation-period directory not found: ", period_data)
  }

  hist_start_year <- as.integer(substr(historical[1], 1, 4))
  hist_end_year <- as.integer(substr(historical[2], 1, 4))
  future_start_year <- as.integer(substr(future[1], 1, 4))
  future_end_year <- as.integer(substr(future[2], 1, 4))

  if (anyNA(c(
    hist_start_year,
    hist_end_year,
    future_start_year,
    future_end_year
  ))) {
    stop("Could not parse years from historical/future inputs.")
  }

  if (
    hist_start_year > hist_end_year ||
      future_start_year > future_end_year
  ) {
    stop("Period start years must be less than or equal to end years.")
  }

  if (verbose) {
    message("Starting consecutive dry-day exposure calculation...")
    start_time <- Sys.time()
  }

  period_ds <- arrow::open_dataset(period_data, format = "parquet")

  required_cols <- c(
    "wsd_source_id",
    "park_unit",
    "year",
    "model",
    "ssp",
    "run",
    "dry_spell_max_length"
  )

  missing_cols <- setdiff(required_cols, names(period_ds))

  if (length(missing_cols)) {
    stop(
      "period_data is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  ws_ids <- unique(water_supplies$wsd_source_id)

  annual_data <- period_ds %>%
    dplyr::select(dplyr::all_of(required_cols)) %>%
    dplyr::filter(
      wsd_source_id %in% ws_ids,
      (
        ssp == "historical" &
          year >= hist_start_year &
          year <= hist_end_year
      ) |
        (
          ssp != "historical" &
            year >= future_start_year &
            year <= future_end_year
        )
    ) %>%
    dplyr::collect() %>%
    data.table::as.data.table()

  if (!nrow(annual_data)) {
    stop("No dry-period records remained after filtering.")
  }

  annual_data[
    !is.finite(dry_spell_max_length),
    dry_spell_max_length := NA_real_
  ]

  historical_means <- annual_data[
    ssp == "historical",
    .(
      cdd_historical = mean(dry_spell_max_length, na.rm = TRUE),
      hist_years_used = sum(!is.na(dry_spell_max_length))
    ),
    by = .(wsd_source_id, park_unit, model, run)
  ]

  future_means <- annual_data[
    ssp != "historical",
    .(
      cdd_future = mean(dry_spell_max_length, na.rm = TRUE),
      future_years_used = sum(!is.na(dry_spell_max_length))
    ),
    by = .(wsd_source_id, park_unit, model, ssp, run)
  ]

  dry_change <- merge(
    future_means,
    historical_means,
    by = c("wsd_source_id", "park_unit", "model", "run"),
    all = FALSE
  )

  dry_change[
    ,
    `:=`(
      cdd_change = cdd_future - cdd_historical,
      cdd_percent_change = data.table::fifelse(
        is.finite(cdd_historical) & cdd_historical > 0,
        100 * (cdd_future - cdd_historical) / cdd_historical,
        NA_real_
      )
    )
  ]

  dry_change <- dry_change[
    hist_years_used > 0L &
      future_years_used > 0L &
      is.finite(cdd_change) &
      is.finite(cdd_percent_change)
  ]

  if (!nrow(dry_change)) {
    stop(
      "No paired historical/future dry-period percent changes could be calculated."
    )
  }

  summary_dt <- dry_change[
    ,
    {
      change_quantiles <- stats::quantile(
        cdd_change,
        probs = c(0.10, 0.50, 0.90),
        na.rm = TRUE
      )
      percent_quantiles <- stats::quantile(
        cdd_percent_change,
        probs = c(0.10, 0.50, 0.90),
        na.rm = TRUE
      )

      closest_value <- function(values, target) {
        values[which.min(abs(values - target))]
      }

      .(
        raw_cdd_change_p10 =
          closest_value(cdd_change, change_quantiles[[1]]),
        raw_cdd_change_p50 =
          closest_value(cdd_change, change_quantiles[[2]]),
        raw_cdd_change_p90 =
          closest_value(cdd_change, change_quantiles[[3]]),
        raw_cdd_percent_change_p10 =
          closest_value(cdd_percent_change, percent_quantiles[[1]]),
        raw_cdd_percent_change_p50 =
          closest_value(cdd_percent_change, percent_quantiles[[2]]),
        raw_cdd_percent_change_p90 =
          closest_value(cdd_percent_change, percent_quantiles[[3]]),
        raw_cdd_change_models_showing_increase =
          100 * sum(cdd_percent_change > 0) / .N,
        raw_cdd_change_model_count = .N
      )
    },
    by = .(wsd_source_id, park_unit)
  ]

  final_results <- summary_dt %>%
    # Longer projected dry spells indicate greater exposure. Standardize only
    # percent change so the parallel absolute-day fields are not double-counted.
    add_standardized_indicators("raw_cdd_percent_change_p") %>%
    add_standardized_indicators(
      "raw_cdd_change_models_showing_increase"
    ) %>%
    dplyr::select(wsd_source_id, park_unit, dplyr::everything())

  if (save) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    output_path <- file.path(
      output_dir,
      paste0("dry_period_exposure_", Sys.Date(), ".csv")
    )

    readr::write_csv(final_results, output_path)

    if (verbose) message("Results saved to: ", output_path)
  }

  if (verbose) {
    elapsed <- Sys.time() - start_time
    message(
      "Consecutive dry-day exposure completed in ",
      round(elapsed, 2), " ", attr(elapsed, "units")
    )
  }

  dplyr::as_tibble(final_results)
}
