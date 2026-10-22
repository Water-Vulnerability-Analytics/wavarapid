#' @title Drought sensitivity: historical drought trend, frequency and severity
#'
#' @description
#' Calculates historical drought indicators at each water supply from the Jay
#' Alder MWBM training run. Monthly precipitation and PET are averaged over
#' each site's grid cells and joined by month, and the climatic water balance
#' (precipitation minus PET) is turned into SPEI with [SPEI::spei()], using
#' the historical period as the reference. `Inf` and `-Inf` SPEI values are
#' capped at `spei_inf_cap` and `-spei_inf_cap`. Sites need SPEI for at least
#' 90% of the historical months.
#'
#' From the monthly SPEI it calculates:
#' * the trend in annual mean SPEI: Sen's slope ([trend::sens.slope()]),
#'   tested with Mann-Kendall ([trend::mk.test()]), for sites with at least 3
#'   years;
#' * drought frequency: the percentage of months with
#'   `SPEI <= drought_threshold`;
#' * drought severity: the mean of `abs(SPEI)` over those months, so larger
#'   means more severe.
#'
#' @param water_supplies sf object. The water supplies, with `wsd_source_id`
#'   and `park_unit` columns, e.g. from [wavacore::get_wsd()]. Only their IDs
#'   are used, to pick rows from the files. Required.
#' @param precip_data String. The folder with the processed training
#'   precipitation file,
#'   `jay_alder_pr_MWBM_LOCA2_training_195001-201412.nc.parquet`; stops with
#'   an error if it isn't there. Default:
#'   `rapid/inputs/loca2_precipitation_water_supplies` in the shared data
#'   folder, see [wavacore::wava_data_path()].
#' @param pet_data String. The folder with the processed training PET file,
#'   `jay_alder_pet_MWBM_LOCA2_training_195001-201412.nc.parquet`; stops with
#'   an error if it isn't there. The name is matched exactly, so on a
#'   case-sensitive file system a file named `Jay_Alder_pet_...` is not found.
#'   Default: `rapid/inputs/loca2_pet_water_supplies` in the shared data
#'   folder.
#' @param spei_scale Integer. The SPEI accumulation scale in months, passed
#'   to [SPEI::spei()]. Default: `12`.
#' @param spei_inf_cap Numeric. The absolute value that replaces `Inf` and
#'   `-Inf` SPEI values. Must be positive. Default: `4`.
#' @param drought_threshold Numeric. The SPEI value at or below which a month
#'   counts as a drought month. Default: `-1` (moderate drought or worse).
#' @param historical Character or numeric vector of 2 years, e.g.
#'   `c(1980, 2010)`. The first and last year of the historical period. The
#'   training run covers 1950-2014. Required.
#' @param water_year Logical. How months are grouped into years for the
#'   annual means and the trend.
#'   * `TRUE`: water years, October to September, labelled by the year they
#'     end in.
#'   * `FALSE`: calendar years.
#'
#'   Default: `TRUE`.
#' @param save Logical. Whether to write the results to a CSV file.
#'   * `TRUE`: writes `drought_sensitivity_<date>.csv`, named with today's
#'     date, to `output_dir`, creating the folder if needed. Stops with an
#'     error if `output_dir` is `NULL`.
#'   * `FALSE`: writes nothing.
#'
#'   Default: `FALSE`.
#' @param output_dir String or `NULL`. The folder for the CSV file.
#'   * A folder path: where the file is written when `save = TRUE`. Ignored
#'     when `save = FALSE`.
#'   * `NULL`: only allowed with `save = FALSE`.
#'
#'   Default: `NULL`.
#' @param verbose Logical. Whether to print progress.
#'   * `TRUE`: prints each step, the run time and the path of the saved file.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return A tibble, one row per water supply, with `wsd_source_id`,
#'   `park_unit` and:
#'   * `raw_slope_spei`: Sen's slope of annual mean SPEI (per year), with
#'     `n_years`, `mk_significant` (p < 0.05) and `trend_direction`
#'     (`"becoming wetter"`, `"becoming drier"` or `"no change"`);
#'   * `raw_avg_spei`, `raw_drought_month_frequency` (%) and
#'     `raw_avg_drought_severity`;
#'   * `raw_spei_*inf*` and `raw_spei_months`: counts and the percentage of
#'     capped `Inf` and `-Inf` months;
#'   * `norm_*` and `rank_*` versions of the slope (inverted, so drying gets
#'     the highest score), the frequency and the severity, from
#'     [add_standardized_indicators()].
#' @export
#' @seealso [calc_drought_exposure()] for the projected change;
#'   [calc_precip_sensitivity()], which reads the same precipitation file.
#' @examplesIf file.exists(wavacore::wava_data_path("rapid", "inputs", "loca2_pet_water_supplies", "jay_alder_pet_MWBM_LOCA2_training_195001-201412.nc.parquet"))
#' ws <- wavacore::get_wsd(park = "ZION")
#' drought_sen <- calc_drought_sensitivity(ws, historical = c(1980, 2010))
#' drought_sen[, c("wsd_source_id", "raw_slope_spei", "trend_direction",
#'                 "raw_drought_month_frequency", "raw_avg_drought_severity")]

calc_drought_sensitivity <- function(water_supplies,
                                     precip_data = wavacore::wava_data_path("rapid", "inputs", "loca2_precipitation_water_supplies/"),
                                     pet_data = wavacore::wava_data_path("rapid", "inputs", "loca2_pet_water_supplies/"),
                                     spei_scale = 12,
                                     spei_inf_cap = 4,
                                     drought_threshold = -1,
                                     historical,
                                     water_year = TRUE,
                                     save = FALSE,
                                     output_dir = NULL,
                                     verbose = TRUE) {
  
  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "water_supplies must include wsd_source_id" = "wsd_source_id" %in% names(water_supplies),
    "water_supplies must include park_unit" = "park_unit" %in% names(water_supplies),
    "spei_scale must be a positive integer" = length(spei_scale) == 1 &&
      is.finite(spei_scale) &&
      spei_scale >= 1 &&
      spei_scale == as.integer(spei_scale),
    "spei_inf_cap must be a positive finite number" = length(spei_inf_cap) == 1 &&
      is.finite(spei_inf_cap) &&
      spei_inf_cap > 0,
    "drought_threshold must be a finite number" = length(drought_threshold) == 1 &&
      is.finite(drought_threshold),
    "historical must be length 2" = length(historical) == 2,
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir)
  )
  
  start_time <- Sys.time()
  
  hist_start_year <- as.integer(historical[1])
  hist_end_year <- as.integer(historical[2])
  
  if (anyNA(c(hist_start_year, hist_end_year))) {
    stop("Could not parse years from historical input.")
  }
  
  if (hist_start_year > hist_end_year) {
    stop("historical start year must be less than or equal to the end year.")
  }
  
  precip_file <- file.path(
    precip_data,
    "jay_alder_pr_MWBM_LOCA2_training_195001-201412.nc.parquet"
  )
  
  pet_file <- file.path(
    pet_data,
    "jay_alder_pet_MWBM_LOCA2_training_195001-201412.nc.parquet"
  )
  
  if (!file.exists(precip_file)) {
    stop("Historical precipitation file not found: ", precip_file)
  }
  
  if (!file.exists(pet_file)) {
    stop("Historical PET file not found: ", pet_file)
  }
  
  if (water_year) {
    hist_report_start_month <- data.table::as.IDate(
      paste0(hist_start_year - 1L, "-10-01")
    )
    hist_report_end_month <- data.table::as.IDate(
      paste0(hist_end_year, "-09-01")
    )
  } else {
    hist_report_start_month <- data.table::as.IDate(
      paste0(hist_start_year, "-01-01")
    )
    hist_report_end_month <- data.table::as.IDate(
      paste0(hist_end_year, "-12-01")
    )
  }
  
  pad_months <- as.integer(spei_scale) - 1L
  
  hist_input_start_month <- seq(
    hist_report_start_month,
    by = "-1 month",
    length.out = pad_months + 1L
  )[pad_months + 1L]
  
  expected_hist_months <- length(
    seq(hist_report_start_month, hist_report_end_month, by = "month")
  )
  
  min_hist_months <- ceiling(expected_hist_months * 0.90)
  ws_ids <- unique(water_supplies$wsd_source_id)
  
  if (verbose) {
    message("Starting historical drought sensitivity calculation...")
    message("Historical precipitation file: ", basename(precip_file))
    message("Historical PET file: ", basename(pet_file))
    message("Historical reporting window: ", hist_report_start_month, " to ", hist_report_end_month)
    message("Historical input window: ", hist_input_start_month, " to ", hist_report_end_month)
    message("Minimum historical SPEI months required: ", min_hist_months, " of ", expected_hist_months)
    message("Infinite SPEI cap: +/-", spei_inf_cap)
  }
  
  # ---------------------------------------------------------------------------
  # Read and aggregate historical precipitation and PET
  # ---------------------------------------------------------------------------
  
  process_historical_var <- function(file, value_col) {
    dt <- arrow::open_dataset(file, format = "parquet") %>%
      dplyr::select(wsd_source_id, park_unit, date, value) %>%
      dplyr::filter(
        wsd_source_id %in% ws_ids,
        date >= as.character(hist_input_start_month),
        date <= as.character(hist_report_end_month)
      ) %>%
      dplyr::collect() %>%
      data.table::as.data.table()
    
    if (!nrow(dt)) {
      return(data.table::data.table())
    }
    
    dt[, date := data.table::as.IDate(date)]
    dt[, ym := date - (data.table::mday(date) - 1L)]
    
    dt <- dt[
      ym >= hist_input_start_month & ym <= hist_report_end_month
    ][
      , .(value = mean(value, na.rm = TRUE)),
      by = .(wsd_source_id, park_unit, ym)
    ]
    
    dt[!is.finite(value), value := NA_real_]
    data.table::setnames(dt, "value", value_col)
    
    dt[]
  }
  
  if (verbose) message("Processing historical precipitation...")
  precip_monthly <- process_historical_var(
    precip_file,
    "monthly_precipitation"
  )
  
  if (verbose) message("Processing historical PET...")
  pet_monthly <- process_historical_var(
    pet_file,
    "monthly_pet"
  )
  
  if (!nrow(precip_monthly)) {
    stop("No historical precipitation records remained after filtering.")
  }
  
  if (!nrow(pet_monthly)) {
    stop("No historical PET records remained after filtering.")
  }
  
  climate_monthly <- merge(
    precip_monthly,
    pet_monthly,
    by = c("wsd_source_id", "park_unit", "ym"),
    all = TRUE
  )
  
  climate_monthly[
    , monthly_balance := monthly_precipitation - monthly_pet
  ]
  
  climate_monthly[
    !is.finite(monthly_balance),
    monthly_balance := NA_real_
  ]
  
  if (!nrow(climate_monthly)) {
    stop("No joinable historical precipitation and PET records were found.")
  }
  
  # ---------------------------------------------------------------------------
  # Calculate historical SPEI for each water supply
  # ---------------------------------------------------------------------------
  
  site_keys <- unique(climate_monthly[, .(wsd_source_id, park_unit)])
  spei_list <- vector("list", nrow(site_keys))
  
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "Calculating historical SPEI [:bar] :percent (:current/:total) ETA: :eta",
      total = nrow(site_keys),
      clear = FALSE,
      width = 80
    )
  }
  
  for (i in seq_len(nrow(site_keys))) {
    if (verbose) pb$tick()
    
    this_wsd <- site_keys$wsd_source_id[i]
    this_pk <- site_keys$park_unit[i]
    
    site_data <- climate_monthly[
      wsd_source_id == this_wsd & park_unit == this_pk,
      .(ym, monthly_balance)
    ]
    
    site_full <- data.table::data.table(
      ym = seq(hist_input_start_month, hist_report_end_month, by = "month")
    )
    
    site_full <- merge(
      site_full,
      site_data,
      by = "ym",
      all.x = TRUE,
      sort = TRUE
    )
    
    report_idx <- site_full$ym >= hist_report_start_month &
      site_full$ym <= hist_report_end_month
    
    ref_dates <- site_full$ym[report_idx & !is.na(site_full$monthly_balance)]
    
    if (!length(ref_dates)) {
      next
    }
    
    balance_ts <- stats::ts(
      site_full$monthly_balance,
      start = c(
        data.table::year(min(site_full$ym)),
        data.table::month(min(site_full$ym))
      ),
      frequency = 12
    )
    
    raw_spei <- tryCatch(
      as.numeric(
        SPEI::spei(
          data = balance_ts,
          scale = spei_scale,
          ref.start = c(
            data.table::year(min(ref_dates)),
            data.table::month(min(ref_dates))
          ),
          ref.end = c(
            data.table::year(max(ref_dates)),
            data.table::month(max(ref_dates))
          ),
          na.rm = TRUE,
          verbose = FALSE
        )[["fitted"]]
      ),
      error = function(e) rep(NA_real_, nrow(site_full))
    )
    
    spei <- raw_spei
    spei[raw_spei == Inf] <- spei_inf_cap
    spei[raw_spei == -Inf] <- -spei_inf_cap
    spei[is.nan(spei)] <- NA_real_
    
    site_report <- data.table::data.table(
      wsd_source_id = this_wsd,
      park_unit = this_pk,
      ym = site_full$ym[report_idx],
      spei = spei[report_idx],
      spei_raw = raw_spei[report_idx]
    )
    
    site_report[, year := if (water_year) {
      data.table::year(ym) + data.table::fifelse(data.table::month(ym) >= 10L, 1L, 0L)
    } else {
      data.table::year(ym)
    }]
    
    if (sum(!is.na(site_report$spei)) < min_hist_months) {
      next
    }
    
    spei_list[[i]] <- site_report
  }
  
  monthly_spei <- data.table::rbindlist(
    spei_list,
    use.names = TRUE,
    fill = TRUE
  )
  
  if (!nrow(monthly_spei)) {
    stop("No historical SPEI values could be calculated.")
  }
  
  # ---------------------------------------------------------------------------
  # Summarize monthly and annual drought conditions
  # ---------------------------------------------------------------------------
  
  monthly_summary <- monthly_spei[
    ,
    .(
      raw_avg_spei = if (any(!is.na(spei))) {
        mean(spei, na.rm = TRUE)
      } else {
        NA_real_
      },
      raw_drought_month_frequency = if (sum(!is.na(spei)) > 0L) {
        100 * sum(spei <= drought_threshold, na.rm = TRUE) / sum(!is.na(spei))
      } else {
        NA_real_
      },
      raw_avg_drought_severity = if (any(spei <= drought_threshold, na.rm = TRUE)) {
        mean(abs(spei[spei <= drought_threshold]), na.rm = TRUE)
      } else {
        0
      },
      raw_spei_inf_months = sum(is.infinite(spei_raw), na.rm = TRUE),
      raw_spei_pos_inf_months = sum(spei_raw == Inf, na.rm = TRUE),
      raw_spei_neg_inf_months = sum(spei_raw == -Inf, na.rm = TRUE),
      raw_spei_months = sum(!is.na(spei)),
      raw_spei_pct_inf_months = if (.N > 0L) {
        100 * sum(is.infinite(spei_raw), na.rm = TRUE) / .N
      } else {
        NA_real_
      }
    ),
    by = .(wsd_source_id, park_unit)
  ]
  
  annual_spei <- monthly_spei[
    ,
    .(
      annual_mean_spei = if (sum(!is.na(spei)) >= 9L) {
        mean(spei, na.rm = TRUE)
      } else {
        NA_real_
      }
    ),
    by = .(wsd_source_id, park_unit, year)
  ]
  
  trend_results <- annual_spei %>%
    dplyr::as_tibble() %>%
    dplyr::arrange(wsd_source_id, year) %>%
    dplyr::group_by(wsd_source_id, park_unit) %>%
    dplyr::summarize(
      data = list(dplyr::pick(dplyr::everything())),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      n_years = purrr::map_int(data, ~ sum(!is.na(.x$annual_mean_spei))),
      sen = purrr::map(data, ~ {
        x <- .x$annual_mean_spei
        if (sum(!is.na(x)) >= 3) trend::sens.slope(x) else NULL
      }),
      mk = purrr::map(data, ~ {
        x <- .x$annual_mean_spei
        if (sum(!is.na(x)) >= 3) trend::mk.test(x) else NULL
      }),
      raw_slope_spei = purrr::map_dbl(
        sen,
        ~ if (is.null(.x)) NA_real_ else unname(.x$estimates[[1]])
      ),
      mk_tau = purrr::map_dbl(
        mk,
        ~ if (is.null(.x)) NA_real_ else unname(.x$estimates[[1]])
      ),
      mk_p_value = purrr::map_dbl(
        mk,
        ~ if (is.null(.x)) NA_real_ else .x$p.value
      ),
      mk_significant = mk_p_value < 0.05
    ) %>%
    dplyr::select(-data, -sen, -mk, -mk_tau, -mk_p_value) %>%
    dplyr::mutate(
      trend_direction = dplyr::case_when(
        is.na(raw_slope_spei) ~ NA_character_,
        raw_slope_spei > 0 ~ "becoming wetter",
        raw_slope_spei < 0 ~ "becoming drier",
        TRUE ~ "no change"
      )
    )
  
  final_results <- trend_results %>%
    dplyr::left_join(
      dplyr::as_tibble(monthly_summary),
      by = c("wsd_source_id", "park_unit")
    ) %>%
    # Lower SPEI slopes indicate increasing drought sensitivity.
    add_standardized_indicators("raw_slope_spei", invert = TRUE) %>%
    # Higher drought frequency and severity indicate greater sensitivity.
    add_standardized_indicators("raw_drought_month_frequency") %>%
    add_standardized_indicators("raw_avg_drought_severity") %>%
    dplyr::select(wsd_source_id, park_unit, dplyr::everything())
  
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(
      output_dir,
      paste0("drought_sensitivity_", Sys.Date(), ".csv")
    )
    
    readr::write_csv(final_results, filename)
    
    if (verbose) {
      message("Results saved to: ", filename)
    }
  }
  
  if (verbose) {
    total_time <- Sys.time() - start_time
    message(
      "\nHistorical drought sensitivity processing completed in ",
      round(total_time, 2), " ", attr(total_time, "units")
    )
  }
  
  return(final_results)
}
