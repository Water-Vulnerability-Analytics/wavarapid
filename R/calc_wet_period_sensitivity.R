#' Calculate Historical Consecutive Wet-Day Sensitivity Indicator
#'
#' Calculates historical trends in annual maximum consecutive wet-day spell
#' length (CWD) from the processed LOCA2 precipitation-extremes dataset. Sen's
#' slope and the Mann-Kendall test are calculated for each historical climate
#' model at each water supply. The slope is divided by that model's historical
#' mean CWD and multiplied by 100 to express sensitivity as percent change per
#' year. Absolute slopes in days per year are retained as supporting values.
#' The closest actual model values to the p10, p50, and p90 of the historical
#' ensemble are retained.
#'
#' More-negative percent slopes and a larger percentage of models showing
#' decreasing CWD indicate greater sensitivity to historically shortening wet
#' periods.
#'
#' @param water_supplies An `sf` object containing `wsd_source_id` and
#'   `park_unit`.
#' @param period_data Character. Directory containing parquet files created by
#'   `process_precip_period_data()`.
#' @param historical Character vector of length 2 giving the historical
#'   start/end years.
#' @param save Logical. If `TRUE`, writes results to disk.
#' @param output_dir Character. Output directory when `save = TRUE`.
#' @param verbose Logical. If `TRUE`, prints progress messages.
#'
#' @return A tibble with one row per water supply and historical wet-period
#'   sensitivity indicators.
calc_wet_period_sensitivity <- function(
    water_supplies,
    period_data = "data/rapid/inputs/loca_cmip6_precip_periods_water_supplies/",
    historical,
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
    "save is TRUE but output_dir not provided" =
      !save || !is.null(output_dir)
  )

  if (!dir.exists(period_data)) {
    stop("Processed precipitation-period directory not found: ", period_data)
  }

  hist_start_year <- as.integer(substr(historical[1], 1, 4))
  hist_end_year <- as.integer(substr(historical[2], 1, 4))

  if (anyNA(c(hist_start_year, hist_end_year))) {
    stop("Could not parse years from historical input.")
  }

  if (hist_start_year > hist_end_year) {
    stop("Historical start year must be less than or equal to end year.")
  }

  if (verbose) {
    message("Starting historical consecutive wet-day sensitivity calculation...")
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
    "wet_spell_max_length"
  )

  missing_cols <- setdiff(required_cols, names(period_ds))

  if (length(missing_cols)) {
    stop(
      "period_data is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  ws_ids <- unique(water_supplies$wsd_source_id)

  historical_data <- period_ds %>%
    dplyr::select(dplyr::all_of(required_cols)) %>%
    dplyr::filter(
      wsd_source_id %in% ws_ids,
      ssp == "historical",
      year >= hist_start_year,
      year <= hist_end_year
    ) %>%
    dplyr::collect() %>%
    dplyr::mutate(
      wet_spell_max_length = dplyr::if_else(
        is.finite(wet_spell_max_length),
        wet_spell_max_length,
        NA_real_
      )
    )

  if (!nrow(historical_data)) {
    stop("No historical wet-period records remained after filtering.")
  }

  model_trends <- historical_data %>%
    dplyr::arrange(wsd_source_id, park_unit, model, run, year) %>%
    dplyr::group_by(wsd_source_id, park_unit, model, run) %>%
    dplyr::summarize(
      cwd_historical_mean = mean(wet_spell_max_length, na.rm = TRUE),
      data = list(dplyr::pick(year, wet_spell_max_length)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      n_years = purrr::map_int(
        data,
        ~ sum(!is.na(.x$wet_spell_max_length))
      ),
      sen = purrr::map(data, ~ {
        values <- .x$wet_spell_max_length
        if (sum(!is.na(values)) >= 3L) {
          trend::sens.slope(values)
        } else {
          NULL
        }
      }),
      mk = purrr::map(data, ~ {
        values <- .x$wet_spell_max_length
        if (sum(!is.na(values)) >= 3L) {
          trend::mk.test(values)
        } else {
          NULL
        }
      }),
      cwd_slope = purrr::map_dbl(
        sen,
        ~ if (is.null(.x)) NA_real_ else unname(.x$estimates[[1]])
      ),
      cwd_percent_slope = dplyr::if_else(
        is.finite(cwd_historical_mean) & cwd_historical_mean > 0,
        100 * cwd_slope / cwd_historical_mean,
        NA_real_
      ),
      mk_p_value = purrr::map_dbl(
        mk,
        ~ if (is.null(.x)) NA_real_ else .x$p.value
      )
    ) %>%
    dplyr::select(-data, -sen, -mk) %>%
    dplyr::filter(
      is.finite(cwd_slope),
      is.finite(cwd_percent_slope)
    )

  if (!nrow(model_trends)) {
    stop("No historical wet-period percent Sen slopes could be calculated.")
  }

  summary_dt <- data.table::as.data.table(model_trends)[
    ,
    {
      slope_quantiles <- stats::quantile(
        cwd_slope,
        probs = c(0.10, 0.50, 0.90),
        na.rm = TRUE
      )
      percent_slope_quantiles <- stats::quantile(
        cwd_percent_slope,
        probs = c(0.10, 0.50, 0.90),
        na.rm = TRUE
      )

      closest_value <- function(values, target) {
        values[which.min(abs(values - target))]
      }

      .(
        raw_cwd_slope_p10 =
          closest_value(cwd_slope, slope_quantiles[[1]]),
        raw_cwd_slope_p50 =
          closest_value(cwd_slope, slope_quantiles[[2]]),
        raw_cwd_slope_p90 =
          closest_value(cwd_slope, slope_quantiles[[3]]),
        raw_cwd_percent_slope_p10 =
          closest_value(cwd_percent_slope, percent_slope_quantiles[[1]]),
        raw_cwd_percent_slope_p50 =
          closest_value(cwd_percent_slope, percent_slope_quantiles[[2]]),
        raw_cwd_percent_slope_p90 =
          closest_value(cwd_percent_slope, percent_slope_quantiles[[3]]),
        raw_cwd_trend_models_showing_decrease =
          100 * sum(cwd_percent_slope < 0) / .N,
        raw_cwd_trend_models_significant_decrease =
          100 *
            sum(
              cwd_percent_slope < 0 & mk_p_value < 0.05,
              na.rm = TRUE
            ) /
            .N,
        raw_cwd_trend_model_count = .N,
        cwd_median_years_used = stats::median(n_years)
      )
    },
    by = .(wsd_source_id, park_unit)
  ]

  final_results <- summary_dt %>%
    # Decreasing historical CWD indicates greater sensitivity. Standardize
    # only percent slope so absolute days/year are not double-counted.
    add_standardized_indicators(
      "raw_cwd_percent_slope_p",
      invert = TRUE
    ) %>%
    add_standardized_indicators(
      "raw_cwd_trend_models_showing_decrease"
    ) %>%
    dplyr::select(wsd_source_id, park_unit, dplyr::everything())

  if (save) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    output_path <- file.path(
      output_dir,
      paste0("wet_period_sensitivity_", Sys.Date(), ".csv")
    )

    readr::write_csv(final_results, output_path)

    if (verbose) message("Results saved to: ", output_path)
  }

  if (verbose) {
    elapsed <- Sys.time() - start_time
    message(
      "Historical consecutive wet-day sensitivity completed in ",
      round(elapsed, 2), " ", attr(elapsed, "units")
    )
  }

  dplyr::as_tibble(final_results)
}
