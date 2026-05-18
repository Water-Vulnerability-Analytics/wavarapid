#' Calculate Drought Exposure Indicator
#'
#' Computes projected change in drought conditions for NPS water supplies using
#' the joined monthly drought output from `process_drought_data()`.
#'
#' Drought is characterized by the Standardized Precipitation-Evapotranspiration
#' Index (SPEI), calculated at a user-specified accumulation scale from monthly
#' climatic water balance: precipitation minus PET.
#'
#' This function assumes historical/future assembly has already been handled
#' upstream by `process_drought_data()`. Therefore, `hindcast` is intentionally
#' not an argument here. Use `process_drought_data(hindcast = TRUE)` or
#' `process_drought_data(hindcast = FALSE)` to choose the baseline construction.
#'
#' Infinite SPEI values are directionally meaningful but numerically unbounded.
#' Therefore, `-Inf` values are capped at `-spei_inf_cap` and `Inf` values are
#' capped at `spei_inf_cap` before calculating mean historical/future SPEI.
#' Counts of capped positive and negative infinite values are retained as
#' diagnostics in the final output.
#'
#' @param water_supplies An `sf` object containing water supply locations.
#'   Must include `wsd_source_id` and `park_unit`.
#' @param drought_data Character. Path to joined monthly drought parquet created
#'   by `process_drought_data()`.
#' @param spei_scale Positive integer. Accumulation scale in months passed to
#'   `SPEI::spei()`. Default is `12`.
#' @param spei_inf_cap Numeric. Absolute cap used to replace `Inf` and `-Inf`
#'   SPEI values before calculating mean SPEI. Default is `4`.
#' @param historical Character vector of length 2 giving the historical period,
#'   e.g. `c("1995", "2024")`.
#' @param future Character vector of length 2 giving the future period,
#'   e.g. `c("2035", "2065")`.
#' @param water_year Logical. If `TRUE`, reporting windows run Oct-Sep.
#' @param save Logical. If `TRUE`, writes results to disk.
#' @param output_dir Character. Output directory when `save = TRUE`.
#' @param verbose Logical. If `TRUE`, prints progress messages.
#'
#' @return A tibble with drought exposure indicators by water supply, including
#'   SPEI fit diagnostics for capped `Inf` and `-Inf` monthly values.
calc_drought_exposure <- function(water_supplies,
                                  drought_data = "data/rapid/inputs/loca2_drought_monthly_1980-2010_2035-2065_hindcast-TRUE.parquet",
                                  spei_scale = 12,
                                  spei_inf_cap = 4,
                                  historical,
                                  future,
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
    "historical must be length 2" = length(historical) == 2,
    "future must be length 2" = length(future) == 2,
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir)
  )
  
  data.table::setDTthreads(1L)
  Sys.setenv(ARROW_NUM_THREADS = 1L)
  
  if (verbose) {
    message("Starting drought exposure calculation...")
    start_time <- Sys.time()
  }
  
  hist_start_year   <- as.integer(substr(historical[1], 1, 4))
  hist_end_year     <- as.integer(substr(historical[2], 1, 4))
  future_start_year <- as.integer(substr(future[1], 1, 4))
  future_end_year   <- as.integer(substr(future[2], 1, 4))
  
  if (anyNA(c(hist_start_year, hist_end_year, future_start_year, future_end_year))) {
    stop("Could not parse years from historical/future inputs.")
  }
  
  if (water_year) {
    hist_report_start_month   <- data.table::as.IDate(paste0(hist_start_year - 1L, "-10-01"))
    hist_report_end_month     <- data.table::as.IDate(paste0(hist_end_year, "-09-01"))
    future_report_start_month <- data.table::as.IDate(paste0(future_start_year - 1L, "-10-01"))
    future_report_end_month   <- data.table::as.IDate(paste0(future_end_year, "-09-01"))
  } else {
    hist_report_start_month   <- data.table::as.IDate(paste0(hist_start_year, "-01-01"))
    hist_report_end_month     <- data.table::as.IDate(paste0(hist_end_year, "-12-01"))
    future_report_start_month <- data.table::as.IDate(paste0(future_start_year, "-01-01"))
    future_report_end_month   <- data.table::as.IDate(paste0(future_end_year, "-12-01"))
  }
  
  pad_months <- as.integer(spei_scale) - 1L
  
  hist_input_start_month <- seq(
    hist_report_start_month,
    by = "-1 month",
    length.out = pad_months + 1L
  )[pad_months + 1L]
  
  hist_input_end_month <- hist_report_end_month
  
  future_input_start_month <- seq(
    future_report_start_month,
    by = "-1 month",
    length.out = pad_months + 1L
  )[pad_months + 1L]
  
  future_input_end_month <- future_report_end_month
  
  expected_hist_months <- length(seq(hist_report_start_month, hist_report_end_month, by = "month"))
  expected_future_months <- length(seq(future_report_start_month, future_report_end_month, by = "month"))
  
  min_hist_months <- ceiling(expected_hist_months * 0.90)
  min_future_months <- ceiling(expected_future_months * 0.90)
  
  if (verbose) {
    message("Historical reporting window: ", hist_report_start_month, " to ", hist_report_end_month)
    message("Future reporting window: ", future_report_start_month, " to ", future_report_end_month)
    message("Historical input window: ", hist_input_start_month, " to ", hist_input_end_month)
    message("Future input window: ", future_input_start_month, " to ", future_input_end_month)
    message("Minimum historical SPEI months required: ", min_hist_months, " of ", expected_hist_months)
    message("Minimum future SPEI months required: ", min_future_months, " of ", expected_future_months)
    message("Infinite SPEI cap: +/-", spei_inf_cap)
  }
  
  # ---------------------------------------------------------------------------
  # Read processed monthly drought data
  # ---------------------------------------------------------------------------
  
  ds <- arrow::open_dataset(drought_data, format = "parquet")
  
  required_cols <- c(
    "wsd_source_id",
    "park_unit",
    "ym",
    "year",
    "month",
    "model",
    "ssp",
    "run",
    "drought_period",
    "monthly_precipitation",
    "monthly_pet"
  )
  
  optional_cols <- c(
    "climatic_water_balance",
    "source_priority",
    "source_type",
    "input_window",
    "reporting_period"
  )
  
  available_cols <- names(ds)
  missing_cols <- setdiff(required_cols, available_cols)
  
  if (length(missing_cols)) {
    stop(
      "drought_data is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "\nThis function expects the joined monthly output from process_drought_data()."
    )
  }
  
  select_cols <- unique(c(required_cols, intersect(optional_cols, available_cols)))
  
  climate_monthly <- ds |>
    dplyr::select(dplyr::all_of(select_cols)) |>
    dplyr::filter(
      wsd_source_id %in% water_supplies$wsd_source_id,
      park_unit %in% water_supplies$park_unit
    ) |>
    dplyr::collect() |>
    data.table::as.data.table()
  
  if (!nrow(climate_monthly)) {
    stop("No drought records found for the selected water supplies.")
  }
  
  climate_monthly[, ym := data.table::as.IDate(ym)]
  
  if ("climatic_water_balance" %in% names(climate_monthly)) {
    climate_monthly[, monthly_balance := climatic_water_balance]
  } else {
    climate_monthly[, monthly_balance := monthly_precipitation - monthly_pet]
  }
  
  climate_monthly[!is.finite(monthly_balance), monthly_balance := NA_real_]
  
  if (!"source_priority" %in% names(climate_monthly)) {
    climate_monthly[, source_priority := 1L]
  }
  
  climate_monthly <- climate_monthly[
    (
      drought_period == "historical" &
        ym >= hist_input_start_month &
        ym <= hist_input_end_month
    ) |
      (
        drought_period == "future" &
          ym >= future_input_start_month &
          ym <= future_input_end_month
      )
  ]
  
  if (!nrow(climate_monthly)) {
    stop("No drought records remained after filtering to SPEI input windows.")
  }
  
  # Dedupe accidental duplicate source rows. The processed data should already
  # be deduped, but this keeps the exposure function safe.
  data.table::setorder(
    climate_monthly,
    wsd_source_id,
    park_unit,
    model,
    ssp,
    run,
    drought_period,
    ym,
    source_priority
  )
  
  climate_monthly <- climate_monthly[
    ,
    .SD[1L],
    by = .(wsd_source_id, park_unit, model, ssp, run, drought_period, ym)
  ]
  
  # ---------------------------------------------------------------------------
  # SPEI helper with capped Inf / -Inf diagnostics
  # ---------------------------------------------------------------------------
  
  run_spei <- function(series_dt, date_col, ref_idx, report_idx = NULL) {
    
    if (is.null(report_idx)) {
      report_idx <- rep(TRUE, nrow(series_dt))
    }
    
    empty_diag <- list(
      spei = rep(NA_real_, nrow(series_dt)),
      spei_raw = rep(NA_real_, nrow(series_dt)),
      n_inf = 0L,
      n_pos_inf = 0L,
      n_neg_inf = 0L,
      n_report_inf = 0L,
      n_report_pos_inf = 0L,
      n_report_neg_inf = 0L
    )
    
    if (!nrow(series_dt)) {
      return(empty_diag)
    }
    
    if (!any(ref_idx & !is.na(series_dt$monthly_balance))) {
      return(empty_diag)
    }
    
    date_vals <- series_dt[[date_col]]
    ref_dates <- date_vals[ref_idx & !is.na(series_dt$monthly_balance)]
    
    ts_start <- c(
      as.integer(format(min(date_vals), "%Y")),
      as.integer(format(min(date_vals), "%m"))
    )
    
    ref_start <- c(
      as.integer(format(min(ref_dates), "%Y")),
      as.integer(format(min(ref_dates), "%m"))
    )
    
    ref_end <- c(
      as.integer(format(max(ref_dates), "%Y")),
      as.integer(format(max(ref_dates), "%m"))
    )
    
    balance_ts <- stats::ts(
      series_dt$monthly_balance,
      start = ts_start,
      frequency = 12
    )
    
    raw_out <- tryCatch(
      as.numeric(
        SPEI::spei(
          data = balance_ts,
          scale = spei_scale,
          ref.start = ref_start,
          ref.end = ref_end,
          na.rm = TRUE,
          verbose = FALSE
        )[["fitted"]]
      ),
      error = function(e) rep(NA_real_, nrow(series_dt))
    )
    
    n_pos_inf <- sum(raw_out == Inf, na.rm = TRUE)
    n_neg_inf <- sum(raw_out == -Inf, na.rm = TRUE)
    n_inf <- n_pos_inf + n_neg_inf
    
    n_report_pos_inf <- sum(report_idx & raw_out == Inf, na.rm = TRUE)
    n_report_neg_inf <- sum(report_idx & raw_out == -Inf, na.rm = TRUE)
    n_report_inf <- n_report_pos_inf + n_report_neg_inf
    
    capped_out <- raw_out
    capped_out[raw_out == Inf] <- spei_inf_cap
    capped_out[raw_out == -Inf] <- -spei_inf_cap
    capped_out[is.nan(capped_out)] <- NA_real_
    
    list(
      spei = capped_out,
      spei_raw = raw_out,
      n_inf = n_inf,
      n_pos_inf = n_pos_inf,
      n_neg_inf = n_neg_inf,
      n_report_inf = n_report_inf,
      n_report_pos_inf = n_report_pos_inf,
      n_report_neg_inf = n_report_neg_inf
    )
  }
  
  # ---------------------------------------------------------------------------
  # Build SPEI change for one model / SSP / run
  # ---------------------------------------------------------------------------
  
  build_spei_change <- function(historical_data, future_data, model_id, ssp_id, run_id) {
    
    if (!nrow(historical_data) || !nrow(future_data)) {
      return(data.table::data.table())
    }
    
    site_keys <- unique(
      rbind(
        historical_data[, .(wsd_source_id, park_unit)],
        future_data[, .(wsd_source_id, park_unit)]
      )
    )
    
    out <- vector("list", nrow(site_keys))
    
    for (i in seq_len(nrow(site_keys))) {
      this_wsd <- site_keys$wsd_source_id[i]
      this_pk  <- site_keys$park_unit[i]
      
      hist_site <- historical_data[
        wsd_source_id == this_wsd &
          park_unit == this_pk
      ]
      
      fut_site <- future_data[
        wsd_source_id == this_wsd &
          park_unit == this_pk
      ]
      
      if (!nrow(hist_site) || !nrow(fut_site)) {
        out[[i]] <- data.table::data.table(
          wsd_source_id = this_wsd,
          park_unit = this_pk,
          spei_historical = NA_real_,
          spei_future = NA_real_,
          hist_months_used = 0L,
          future_months_used = 0L,
          spei_inf_months = 0L,
          spei_pos_inf_months = 0L,
          spei_neg_inf_months = 0L,
          spei_report_inf_months = 0L,
          spei_report_pos_inf_months = 0L,
          spei_report_neg_inf_months = 0L,
          spei_hist_report_inf_months = 0L,
          spei_hist_report_pos_inf_months = 0L,
          spei_hist_report_neg_inf_months = 0L,
          spei_future_report_inf_months = 0L,
          spei_future_report_pos_inf_months = 0L,
          spei_future_report_neg_inf_months = 0L,
          spei_total_fit_months = 0L,
          spei_report_months = 0L,
          spei_pct_inf_months = NA_real_,
          spei_pct_report_inf_months = NA_real_,
          spei_future_report_pct_neg_inf_months = NA_real_
        )
        next
      }
      
      hist_full <- data.table::data.table(
        ym = seq(hist_input_start_month, hist_input_end_month, by = "month")
      )
      
      hist_full <- merge(
        hist_full,
        hist_site[, .(ym, monthly_balance)],
        by = "ym",
        all.x = TRUE,
        sort = TRUE
      )
      
      hist_idx <- hist_full$ym >= hist_report_start_month &
        hist_full$ym <= hist_report_end_month
      
      hist_spei_result <- run_spei(
        series_dt = hist_full,
        date_col = "ym",
        ref_idx = hist_idx,
        report_idx = hist_idx
      )
      
      hist_spei <- hist_spei_result$spei
      
      fut_full <- data.table::data.table(
        ym = seq(future_input_start_month, future_input_end_month, by = "month")
      )
      
      fut_full <- merge(
        fut_full,
        fut_site[, .(ym, monthly_balance)],
        by = "ym",
        all.x = TRUE,
        sort = TRUE
      )
      
      # Historical months define the reference distribution.
      # Future months are scored against that historical reference by gluing
      # historical + future into one compact sequence for SPEI::spei().
      future_fit <- data.table::rbindlist(
        list(
          hist_full[, .(
            original_ym = ym,
            period = "historical",
            monthly_balance = monthly_balance
          )],
          fut_full[, .(
            original_ym = ym,
            period = "future",
            monthly_balance = monthly_balance
          )]
        ),
        use.names = TRUE
      )
      
      future_fit[
        ,
        spei_ym := seq(
          from = hist_input_start_month,
          by = "month",
          length.out = .N
        )
      ]
      
      future_ref_idx <- future_fit$period == "historical" &
        future_fit$original_ym >= hist_report_start_month &
        future_fit$original_ym <= hist_report_end_month
      
      future_idx <- future_fit$period == "future" &
        future_fit$original_ym >= future_report_start_month &
        future_fit$original_ym <= future_report_end_month
      
      future_spei_result <- run_spei(
        series_dt = future_fit,
        date_col = "spei_ym",
        ref_idx = future_ref_idx,
        report_idx = future_idx
      )
      
      future_spei <- future_spei_result$spei
      
      hist_n <- sum(!is.na(hist_spei[hist_idx]))
      future_n <- sum(!is.na(future_spei[future_idx]))
      
      total_inf <- hist_spei_result$n_inf + future_spei_result$n_inf
      total_pos_inf <- hist_spei_result$n_pos_inf + future_spei_result$n_pos_inf
      total_neg_inf <- hist_spei_result$n_neg_inf + future_spei_result$n_neg_inf
      
      hist_report_inf <- hist_spei_result$n_report_inf
      hist_report_pos_inf <- hist_spei_result$n_report_pos_inf
      hist_report_neg_inf <- hist_spei_result$n_report_neg_inf
      
      future_report_inf <- future_spei_result$n_report_inf
      future_report_pos_inf <- future_spei_result$n_report_pos_inf
      future_report_neg_inf <- future_spei_result$n_report_neg_inf
      
      report_inf <- hist_report_inf + future_report_inf
      report_pos_inf <- hist_report_pos_inf + future_report_pos_inf
      report_neg_inf <- hist_report_neg_inf + future_report_neg_inf
      
      hist_report_months <- length(hist_spei[hist_idx])
      future_report_months <- length(future_spei[future_idx])
      
      total_fit_months <- length(hist_spei) + length(future_spei)
      total_report_months <- hist_report_months + future_report_months
      
      out[[i]] <- data.table::data.table(
        wsd_source_id = this_wsd,
        park_unit = this_pk,
        
        spei_historical = if (hist_n > 0L) mean(hist_spei[hist_idx], na.rm = TRUE) else NA_real_,
        spei_future = if (future_n > 0L) mean(future_spei[future_idx], na.rm = TRUE) else NA_real_,
        
        hist_months_used = hist_n,
        future_months_used = future_n,
        
        # All fitted months, including padding/calibration/future-fit sequence.
        spei_inf_months = total_inf,
        spei_pos_inf_months = total_pos_inf,
        spei_neg_inf_months = total_neg_inf,
        
        # Reporting-window infinities only.
        spei_report_inf_months = report_inf,
        spei_report_pos_inf_months = report_pos_inf,
        spei_report_neg_inf_months = report_neg_inf,
        
        # Historical reporting-window infinities only.
        spei_hist_report_inf_months = hist_report_inf,
        spei_hist_report_pos_inf_months = hist_report_pos_inf,
        spei_hist_report_neg_inf_months = hist_report_neg_inf,
        
        # Future reporting-window infinities only.
        spei_future_report_inf_months = future_report_inf,
        spei_future_report_pos_inf_months = future_report_pos_inf,
        spei_future_report_neg_inf_months = future_report_neg_inf,
        
        spei_total_fit_months = total_fit_months,
        spei_report_months = total_report_months,
        
        spei_pct_inf_months = if (total_fit_months > 0L) {
          100 * total_inf / total_fit_months
        } else {
          NA_real_
        },
        
        spei_pct_report_inf_months = if (total_report_months > 0L) {
          100 * report_inf / total_report_months
        } else {
          NA_real_
        },
        
        spei_future_report_pct_neg_inf_months = if (future_report_months > 0L) {
          100 * future_report_neg_inf / future_report_months
        } else {
          NA_real_
        }
      )
    }
    
    spei_change <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
    
    spei_change[
      ,
      `:=`(
        spei_change = spei_future - spei_historical,
        model = model_id,
        ssp = ssp_id,
        run = run_id
      )
    ]
    
    spei_change[
      hist_months_used >= min_hist_months &
        future_months_used >= min_future_months &
        is.finite(spei_change)
    ]
  }
  
  # ---------------------------------------------------------------------------
  # Run SPEI by model / SSP / run combo
  # ---------------------------------------------------------------------------
  
  future_ssps <- unique(
    climate_monthly[
      drought_period == "future" &
        !is.na(model) &
        !is.na(ssp) &
        !is.na(run),
      .(model, ssp, run)
    ]
  )
  
  if (!nrow(future_ssps)) {
    stop("No future model/ssp/run combinations were found in drought_data.")
  }
  
  run_one <- function(i) {
    ssp_sub <- future_ssps[i]
    
    if (verbose) {
      message(
        "  Running SPEI for ", ssp_sub$model,
        " / ", ssp_sub$ssp,
        " / ", ssp_sub$run,
        " (", i, "/", nrow(future_ssps), ")"
      )
    }
    
    historical_series <- climate_monthly[
      model == ssp_sub$model &
        ssp == ssp_sub$ssp &
        run == ssp_sub$run &
        drought_period == "historical" &
        ym >= hist_input_start_month &
        ym <= hist_input_end_month
    ]
    
    future_series <- climate_monthly[
      model == ssp_sub$model &
        ssp == ssp_sub$ssp &
        run == ssp_sub$run &
        drought_period == "future" &
        ym >= future_input_start_month &
        ym <= future_input_end_month
    ]
    
    build_spei_change(
      historical_data = historical_series,
      future_data = future_series,
      model_id = ssp_sub$model,
      ssp_id = ssp_sub$ssp,
      run_id = ssp_sub$run
    )
  }
  
  drought_change_list <- lapply(seq_len(nrow(future_ssps)), run_one)
  
  drought_change <- data.table::rbindlist(
    drought_change_list,
    use.names = TRUE,
    fill = TRUE
  )
  
  rm(drought_change_list)
  invisible(gc())
  
  if (!nrow(drought_change)) {
    stop(
      "No drought change values could be calculated from the available climate series. ",
      "Check that each model/ssp/run has enough historical and future months after padding. ",
      "Minimum historical months required: ", min_hist_months, "; ",
      "minimum future months required: ", min_future_months, "."
    )
  }
  
  # ---------------------------------------------------------------------------
  # Summarize ensemble
  # ---------------------------------------------------------------------------
  
  summary_dt <- drought_change[
    is.finite(spei_change),
    {
      q <- stats::quantile(spei_change, c(0.10, 0.50, 0.90), na.rm = TRUE)
      
      rk <- data.table::frank(abs(spei_change - q["10%"]), ties.method = "first")
      p10 <- spei_change[rk == 1L]
      
      rk <- data.table::frank(abs(spei_change - q["50%"]), ties.method = "first")
      p50 <- spei_change[rk == 1L]
      
      rk <- data.table::frank(abs(spei_change - q["90%"]), ties.method = "first")
      p90 <- spei_change[rk == 1L]
      
      total_fit <- sum(spei_total_fit_months, na.rm = TRUE)
      total_report <- sum(spei_report_months, na.rm = TRUE)
      
      total_inf <- sum(spei_inf_months, na.rm = TRUE)
      total_pos_inf <- sum(spei_pos_inf_months, na.rm = TRUE)
      total_neg_inf <- sum(spei_neg_inf_months, na.rm = TRUE)
      
      report_inf <- sum(spei_report_inf_months, na.rm = TRUE)
      report_pos_inf <- sum(spei_report_pos_inf_months, na.rm = TRUE)
      report_neg_inf <- sum(spei_report_neg_inf_months, na.rm = TRUE)
      
      hist_report_inf <- sum(spei_hist_report_inf_months, na.rm = TRUE)
      hist_report_pos_inf <- sum(spei_hist_report_pos_inf_months, na.rm = TRUE)
      hist_report_neg_inf <- sum(spei_hist_report_neg_inf_months, na.rm = TRUE)
      
      future_report_inf <- sum(spei_future_report_inf_months, na.rm = TRUE)
      future_report_pos_inf <- sum(spei_future_report_pos_inf_months, na.rm = TRUE)
      future_report_neg_inf <- sum(spei_future_report_neg_inf_months, na.rm = TRUE)
      
      future_report_months <- sum(future_months_used, na.rm = TRUE)
      
      .(
        raw_spei_change_p10 = p10,
        raw_spei_change_p50 = p50,
        raw_spei_change_p90 = p90,
        raw_spei_models_showing_decrease = 100 * sum(spei_change < 0) / .N,
        raw_spei_model_count = .N,
        
        # SPEI infinite-value diagnostics. These are QC/severity fields.
        raw_spei_inf_months = total_inf,
        raw_spei_pos_inf_months = total_pos_inf,
        raw_spei_neg_inf_months = total_neg_inf,
        raw_spei_total_fit_months = total_fit,
        raw_spei_pct_inf_months = if (total_fit > 0) {
          100 * total_inf / total_fit
        } else {
          NA_real_
        },
        
        # Reporting-window infinite diagnostics.
        raw_spei_report_inf_months = report_inf,
        raw_spei_report_pos_inf_months = report_pos_inf,
        raw_spei_report_neg_inf_months = report_neg_inf,
        raw_spei_report_months = total_report,
        raw_spei_pct_report_inf_months = if (total_report > 0) {
          100 * report_inf / total_report
        } else {
          NA_real_
        },
        
        # Historical reporting-window diagnostics.
        raw_spei_hist_report_inf_months = hist_report_inf,
        raw_spei_hist_report_pos_inf_months = hist_report_pos_inf,
        raw_spei_hist_report_neg_inf_months = hist_report_neg_inf,
        
        # Future reporting-window diagnostics.
        raw_spei_future_report_inf_months = future_report_inf,
        raw_spei_future_report_pos_inf_months = future_report_pos_inf,
        raw_spei_future_report_neg_inf_months = future_report_neg_inf,
        raw_spei_future_report_pct_neg_inf_months = if (future_report_months > 0) {
          100 * future_report_neg_inf / future_report_months
        } else {
          NA_real_
        },
        
        # Model-level infinite diagnostics.
        raw_spei_models_with_inf = 100 * sum(spei_inf_months > 0, na.rm = TRUE) / .N,
        raw_spei_models_with_report_inf = 100 * sum(spei_report_inf_months > 0, na.rm = TRUE) / .N,
        raw_spei_models_with_future_report_neg_inf = 100 *
          sum(spei_future_report_neg_inf_months > 0, na.rm = TRUE) / .N
      )
    },
    by = .(park_unit, wsd_source_id)
  ]
  
  final_results <- summary_dt |>
    dplyr::as_tibble() |>
    add_standardized_indicators("raw_spei_change", invert = TRUE) |>
    add_standardized_indicators("raw_spei_models_showing_decrease") |>
    dplyr::select(wsd_source_id, park_unit, dplyr::everything())
  
  if (save) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    filename <- file.path(output_dir, paste0("drought_exposure_", Sys.Date(), ".csv"))
    readr::write_csv(final_results, filename)
    
    if (verbose) {
      message("Results saved to: ", filename)
    }
  }
  
  if (verbose) {
    elapsed <- Sys.time() - start_time
    message("\nCompleted in ", round(elapsed, 2), " ", attr(elapsed, "units"))
  }
  
  final_results
}
