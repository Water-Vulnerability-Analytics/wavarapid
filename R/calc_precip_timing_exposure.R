#' Calculate Precipitation Timing Exposure Indicator (Center-of-Timing)
#'
#' Computes a change-in-precipitation-timing exposure metric for NPS water
#' supplies using LOCA2-driven MWBM monthly precipitation outputs. Timing is
#' measured with a **circular center of timing (CT)**, the precipitation-weighted
#' temporal centroid of the annual precipitation cycle. The circular calculation
#' is the precipitation analogue of the streamflow center-of-timing concept of
#' Stewart, Cayan & Dettinger (2005) and avoids an artificial discontinuity at
#' the September--October water-year boundary.
#'
#' For each water year and each GCM/SSP/run, monthly precipitation totals are
#' placed at their month-midpoint and combined as weighted vectors. CT is the
#' angle of the resulting vector, expressed as an approximate day of water year
#' (WY starts Oct 1). The vector length R (0--1) is retained as a description of
#' how strongly precipitation is concentrated within the annual cycle (no R
#' cutoff is used).
#'
#' The exposure metric is the **magnitude of CT shift** between periods:
#' Larger abs(Delta CT) indicates greater departure from historical timing.
#' `raw_precip_timing_earlier_frac` reports the % of realizations projecting
#' earlier timing.
#'
#' File-selection and filename-parsing logic mirrors \code{calc_precip_exposure()}
#' (single-year and multi-year files, Livneh observed record, ssp back-fill), so
#' this function is a drop-in timing companion to the precipitation-volume
#' exposure function.
#'
#' @param water_supplies An \code{sf} object with a \code{wsd_source_id} column.
#' @param precip_data Character. Directory of processed LOCA2 MWBM monthly
#'   precipitation parquet files.
#' @param hindcast Logical. If \code{TRUE}, uses hindcast climate data; if
#'   \code{FALSE}, uses LOCA-downscaled observed (Livneh) plus modeled data.
#' @param historical,future Character length-2 vectors of start/end years ("YYYY").
#' @param water_year Logical. If \code{TRUE}, uses water years (Oct-Sep). CT is
#'   defined on the water year, so this should normally be TRUE.
#' @param save Logical. If `TRUE`, results are saved to a file in `output_dir`.
#'   Defaults to `FALSE`.
#' @param output_dir Optional character string specifying the directory where results
#'   will be saved if `save = TRUE`. Ignored otherwise.
#' @param verbose Logical. If `TRUE`, progress messages are printed to the console.
#'   Defaults to `TRUE`.
#'
#' @return A tibble with one row per water supply / park unit: raw p10/p50/p90 of
#'   the CT change, an ensemble-agreement fraction, and standardized indicators.
#'
#' @references
#' Stewart, I.T., Cayan, D.R., & Dettinger, M.D. (2005). Changes toward earlier
#'   streamflow timing across western North America. J. Climate, 18, 1136-1155.
#' Dettinger, M.D., & Cayan, D.R. (1995). Large-scale atmospheric forcing of
#'   recent trends toward early snowmelt runoff in California. J. Climate, 8,
#'   606-623.
calc_precip_timing_exposure <- function(water_supplies,
                                        precip_data = "data/rapid/inputs/loca2_precipitation_water_supplies/",
                                        hindcast = TRUE,
                                        historical,
                                        future,
                                        water_year = TRUE,
                                        save = FALSE,
                                        output_dir = NULL,
                                        verbose = TRUE) {

  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir)
  )

  if (verbose) {
    message("Starting precipitation timing (center-of-timing) exposure calculation...")
    start_time <- Sys.time()
  }

  hist_start_year   <- as.integer(historical[1])
  hist_end_year     <- as.integer(historical[2])
  future_start_year <- as.integer(future[1])
  future_end_year   <- as.integer(future[2])

  # Water year starts Oct 1. Fixed (non-leap) month lengths are adequate: CT is
  # a within-year centroid and the sub-day leap error cancels in the future-
  # minus-historical difference.
  .wy_month_order <- c(10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  .mlen           <- c(31, 30, 31, 31, 28, 31, 30, 31, 30, 31, 31, 30)
  dowy_lookup     <- stats::setNames(cumsum(c(0, head(.mlen, -1))) + .mlen / 2,
                                      .wy_month_order)
  angle_lookup    <- 2 * pi * dowy_lookup / 365

  # ---- File selection (mirrors calc_precip_exposure) --------------------
  if (hindcast == FALSE) {
    files <- list.files(precip_data, full.names = TRUE) %>%
      grep("2015|2045|2075|livneh_unsplit", ., value = TRUE)
  } else {
    files <- list.files(precip_data, full.names = TRUE) %>%
      grep("training|livneh_unsplit", ., invert = TRUE, value = TRUE)
  }

  files_dt <- data.table::data.table(file = files)
  files_dt[, `:=`(
    year_single = as.integer(stringr::str_extract(file, "(?<!\\d)(\\d{4})(?=\\.nc)")),
    year_start  = as.integer(stringr::str_extract(file, "(\\d{4})(?=-\\d{4})")),
    year_end    = as.integer(stringr::str_extract(file, "(?<=-)(\\d{4})(?=\\.)"))
  )]
  files_dt[!is.na(year_single), `:=`(year_start = year_single, year_end = year_single)]

  if (water_year == FALSE) {
    files_subset <- files_dt[
      (year_end >= hist_start_year & year_start <= hist_end_year) |
        (year_end >= future_start_year & year_start <= future_end_year)
    ]
  } else {
    files_subset <- files_dt[
      (year_end >= hist_start_year - 1 & year_start <= hist_end_year) |
        (year_end >= future_start_year - 1 & year_start <= future_end_year)
    ]
  }
  files <- files_subset$file
  if (length(files) == 0) stop("No processed files found. Check input data and paths.")

  if (verbose) {
    message("Processing ", length(files), " precipitation files")
    pb <- progress::progress_bar$new(
      format = "Processing precipitation files [:bar] :percent (:current/:total) ETA: :eta",
      total = length(files), clear = FALSE, width = 80)
  }

  process_raw_grid_data <- function(f) {
    raw <- arrow::read_parquet(f) %>%
      dplyr::filter(wsd_source_id %in% water_supplies$wsd_source_id)

    dt <- data.table::as.data.table(raw)[
      , `:=`(
        year = if (water_year) {
          lubridate::year(date) + data.table::fifelse(lubridate::month(date) >= 10, 1L, 0L)
        } else {
          lubridate::year(date)
        },
        month = lubridate::month(date)
      )
    ]

    # when using observed data, drop modeled overlap with the observed record
    if (!grepl("livneh_unsplit", f) && hindcast == FALSE) {
      dt <- dt[year >= 2019]
    }

    monthly_points <- dt[
      (year >= hist_start_year & year <= hist_end_year) |
        (year >= future_start_year & year <= future_end_year)
    ][
      , .(grid_precip = if (all(is.na(value))) NA_real_ else sum(value, na.rm = TRUE)),
      by = .(wsd_source_id, year, month, park_unit, nc_file, cell_lon, cell_lat)
    ][
      , .(monthly_precip = mean(grid_precip, na.rm = TRUE)),
      by = .(year, month, wsd_source_id, park_unit, nc_file)
    ]

    # Circular center of timing within each year; require all 12 months.
    monthly_points[
      ,
      {
        complete <- length(unique(month)) == 12L && all(!is.na(monthly_precip))
        tot <- if (complete) sum(monthly_precip) else NA_real_
        x <- if (!is.na(tot) && tot > 0) {
          sum(cos(angle_lookup[as.character(month)]) * monthly_precip) / tot
        } else NA_real_
        y <- if (!is.na(tot) && tot > 0) {
          sum(sin(angle_lookup[as.character(month)]) * monthly_precip) / tot
        } else NA_real_
        .(
          timing_angle = if (is.na(x) || is.na(y)) NA_real_ else atan2(y, x) %% (2 * pi),
          timing_x = x,
          timing_y = y,
          timing_r = if (is.na(x) || is.na(y)) NA_real_ else sqrt(x^2 + y^2),
          n_months = if (complete) 12L else sum(!is.na(monthly_precip))
        )
      },
      by = .(year, wsd_source_id, park_unit, nc_file)
    ][n_months == 12L]
  }

  precip_data_list <- purrr::map(files, ~{
    if (verbose) pb$tick()
    process_raw_grid_data(.x)
  })

  combined_data <- dplyr::bind_rows(precip_data_list)
  data.table::setDT(combined_data)

  # ---- Filename parsing (mirrors calc_precip_exposure) ------------------
  suppressMessages(suppressWarnings({
    combined_data[
      nc_file == "Jay_Alder_pr_MWBM_LOCA2_training_195001-201412.nc",
      nc_file := "Jay_Alder_pr_MWBM_LOCA2.historical.training.1950-2100.nc"
    ]
    combined_data[
      , c("var", "model", "ssp", "run", "nc", "misc1", "misc2", "misc3") :=
        data.table::tstrsplit(nc_file, "\\.")
    ]
    combined_data[
      !ssp %in% c("ssp245", "ssp370", "ssp585", "historical"),
      ssp := dplyr::case_when(
        grepl("livneh_unsplit", nc_file, ignore.case = TRUE) ~ "historical",
        TRUE ~ stringr::str_extract(nc_file, "ssp\\d+")
      )
    ]
    combined_data[!is.na(run), c("run", "time") := data.table::tstrsplit(run, "_", fixed = TRUE)]
  }))
  combined_data[, c("var", "nc") := NULL]

  # Circular mean CT over years. Averaging annual x/y components gives each
  # water year equal influence while retaining its seasonal concentration.
  summarize_ct <- function(df) {
    df %>%
      dplyr::group_by(wsd_source_id, park_unit) %>%
      dplyr::summarize(timing_x = mean(timing_x, na.rm = TRUE),
                       timing_y = mean(timing_y, na.rm = TRUE),
                       timing_angle = atan2(timing_y, timing_x) %% (2 * pi),
                       timing_r = sqrt(timing_x^2 + timing_y^2),
                       range = paste0(min(year), "-", max(year)),
                       .groups = "drop")
  }

  if (hindcast == TRUE) {
    ssps <- combined_data %>%
      dplyr::filter(!ssp %in% c("1950-2014", "historical")) %>%
      dplyr::distinct(model, ssp, run)
    precip_change_sub <- vector("list", nrow(ssps))

    for (i in seq_len(nrow(ssps))) {
      ssp_sub <- ssps[i, ]
      hist_hc <- combined_data %>%
        dplyr::filter(model == ssp_sub$model, ssp == "historical", run == ssp_sub$run,
                      year >= hist_start_year, year <= hist_end_year)
      fut_fill <- combined_data %>%
        dplyr::filter(model == ssp_sub$model, ssp == ssp_sub$ssp, run == ssp_sub$run,
                      year >= hist_start_year, year <= hist_end_year,
                      !year %in% hist_hc$year)
      historical_summary <- summarize_ct(dplyr::bind_rows(hist_hc, fut_fill)) %>%
        dplyr::rename(timing_angle_historical = timing_angle,
                      seasonality_historical = timing_r,
                      range_historical = range) %>%
        dplyr::select(-timing_x, -timing_y)
      future_summary <- combined_data %>%
        dplyr::filter(model == ssp_sub$model, ssp == ssp_sub$ssp, run == ssp_sub$run,
                      year >= future_start_year, year <= future_end_year) %>%
        summarize_ct() %>%
        dplyr::rename(timing_angle_future = timing_angle,
                      seasonality_future = timing_r,
                      range_future = range) %>%
        dplyr::select(-timing_x, -timing_y)

      precip_change_sub[[i]] <- future_summary %>%
        dplyr::inner_join(historical_summary, by = c("wsd_source_id", "park_unit")) %>%
        dplyr::mutate(model = ssp_sub$model, ssp = ssp_sub$ssp, run = ssp_sub$run,
                      timing_change = atan2(sin(timing_angle_future - timing_angle_historical),
                                            cos(timing_angle_future - timing_angle_historical)) * 365 / (2 * pi),
                      timing_change_mag = abs(timing_change))
    }
    precip_change <- dplyr::bind_rows(precip_change_sub)

  } else {
    ssps <- combined_data %>%
      dplyr::filter(ssp != "training", ssp != "historical") %>%
      dplyr::distinct(model, ssp, run)

    # observed/training historical CT (shared across future realizations)
    training_hist_ct <- combined_data %>%
      dplyr::filter(ssp == "training" | ssp == "historical",
                    year >= hist_start_year, year <= hist_end_year) %>%
      summarize_ct() %>%
      dplyr::rename(timing_angle_historical = timing_angle,
                    seasonality_historical = timing_r,
                    range_historical = range) %>%
      dplyr::select(-timing_x, -timing_y)

    precip_change_sub <- vector("list", nrow(ssps))
    for (i in seq_len(nrow(ssps))) {
      suppressMessages(suppressWarnings({
        ssp_sub <- ssps[i, ]
        future_summary <- combined_data %>%
          dplyr::filter(model == ssp_sub$model, ssp == ssp_sub$ssp, run == ssp_sub$run,
                        year >= future_start_year, year <= future_end_year) %>%
          summarize_ct() %>%
          dplyr::rename(timing_angle_future = timing_angle,
                        seasonality_future = timing_r,
                        range_future = range) %>%
          dplyr::select(-timing_x, -timing_y)

        precip_change_sub[[i]] <- future_summary %>%
          dplyr::inner_join(training_hist_ct, by = c("wsd_source_id", "park_unit")) %>%
          dplyr::mutate(model = ssp_sub$model, ssp = ssp_sub$ssp, run = ssp_sub$run,
                        timing_change = atan2(sin(timing_angle_future - timing_angle_historical),
                                              cos(timing_angle_future - timing_angle_historical)) * 365 / (2 * pi),
                        timing_change_mag = abs(timing_change))
      }))
    }
    precip_change <- dplyr::bind_rows(precip_change_sub)
  }

  pick_pctl <- function(df, p, outname) {
    df %>%
      dplyr::group_by(park_unit, wsd_source_id) %>%
      dplyr::mutate(target = stats::quantile(timing_change_mag, p, na.rm = TRUE),
                    d = abs(timing_change_mag - target)) %>%
      dplyr::filter(d == min(d, na.rm = TRUE)) %>%
      dplyr::slice(1) %>% dplyr::ungroup() %>%
      dplyr::select(wsd_source_id, park_unit, !!outname := timing_change_mag)
  }

  suppressMessages(suppressWarnings({
    p10 <- pick_pctl(precip_change, 0.10, "raw_precip_timing_change_magnitude_p10")
    p50 <- pick_pctl(precip_change, 0.50, "raw_precip_timing_change_magnitude_p50")
    p90 <- pick_pctl(precip_change, 0.90, "raw_precip_timing_change_magnitude_p90")
    signed_and_r <- precip_change %>%
      dplyr::group_by(park_unit, wsd_source_id) %>%
      dplyr::summarize(
        raw_precip_timing_change_p10 = stats::quantile(timing_change, 0.10, na.rm = TRUE),
        raw_precip_timing_change_p50 = stats::quantile(timing_change, 0.50, na.rm = TRUE),
        raw_precip_timing_change_p90 = stats::quantile(timing_change, 0.90, na.rm = TRUE),
        raw_precip_seasonality_historical_p10 = stats::quantile(seasonality_historical, 0.10, na.rm = TRUE),
        raw_precip_seasonality_historical_p50 = stats::quantile(seasonality_historical, 0.50, na.rm = TRUE),
        raw_precip_seasonality_historical_p90 = stats::quantile(seasonality_historical, 0.90, na.rm = TRUE),
        raw_precip_seasonality_future_p10 = stats::quantile(seasonality_future, 0.10, na.rm = TRUE),
        raw_precip_seasonality_future_p50 = stats::quantile(seasonality_future, 0.50, na.rm = TRUE),
        raw_precip_seasonality_future_p90 = stats::quantile(seasonality_future, 0.90, na.rm = TRUE),
        .groups = "drop"
      )
    model_frac <- precip_change %>%
      dplyr::group_by(park_unit, wsd_source_id) %>%
      dplyr::summarize(raw_precip_timing_earlier_frac =
                         100 * sum(timing_change < 0, na.rm = TRUE) / sum(!is.na(timing_change)),
                       .groups = "drop")
  }))

  final_results <- p10 %>%
    dplyr::left_join(p50, by = c("wsd_source_id", "park_unit")) %>%
    dplyr::left_join(p90, by = c("wsd_source_id", "park_unit")) %>%
    dplyr::left_join(signed_and_r, by = c("wsd_source_id", "park_unit")) %>%
    dplyr::left_join(model_frac, by = c("wsd_source_id", "park_unit")) %>%
    # Larger magnitude = greater departure from historical = more exposed
    add_standardized_indicators("raw_precip_timing_change_magnitude", invert = FALSE) %>%
    add_standardized_indicators("raw_precip_timing_earlier_frac") %>%
    dplyr::select(wsd_source_id, park_unit, dplyr::everything())

  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    filename <- file.path(output_dir, paste0("precip_timing_exposure_", Sys.Date(), ".csv"))
    readr::write_csv(final_results, filename)
    if (verbose) message("Results saved to: ", filename)
  }

  if (verbose) {
    total_time <- Sys.time() - start_time
    message("\nPrecipitation timing (CT) exposure completed in ",
            round(total_time, 2), " ", attr(total_time, "units"))
  }

  final_results
}
