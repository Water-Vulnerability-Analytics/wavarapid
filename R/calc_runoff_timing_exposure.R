#' Calculate Runoff Timing Exposure Indicator (Center-of-Timing)
#'
#' Computes a change-in-runoff-timing exposure metric for NPS water supplies
#' using LOCA2-driven MWBM monthly runoff outputs. Timing is measured with a
#' **circular center of timing (CT)**, the flow-weighted temporal centroid of
#' the annual hydrograph. The circular calculation follows the center-of-timing
#' concept of Stewart, Cayan & Dettinger (2005, J. Climate) while avoiding an
#' artificial discontinuity at the September--October water-year boundary. CT is
#' the standard runoff-timing metric in the western-US snowmelt literature
#' (Stewart et al. 2005; McCabe & Clark 2005; Regonda et al. 2005; Hidalgo et
#' al. 2009; Dudley et al. 2017).
#'
#' For each water year and each GCM/SSP/run, monthly runoff totals are placed at
#' their month-midpoint angles and combined as weighted vectors. CT is the angle
#' of the resulting vector, expressed as an approximate day of water year (WY
#' starts Oct 1).
#'
#' For each site, the direction supported by the most model realizations is
#' selected. The signed p10 shift is used when earlier timing is dominant and
#' the signed p90 shift is used when later timing is dominant. The absolute
#' value of that selected shift is the exposure magnitude. The output also
#' reports the selected direction and the fraction of valid realizations that
#' agree with it.
#'
#' @param water_supplies An \code{sf} object with a \code{wsd_source_id} column.
#' @param runoff_data Character. Directory of processed LOCA2 MWBM monthly
#' runoff parquet files.
#' @param hindcast Logical. If \code{TRUE}, blends model historical runs with
#' overlapping scenario data for the historical period; if \code{FALSE}, uses
#' training data supplemented with modeled data.
#' @param historical,future Character length-2 vectors of start/end years ("YYYY").
#' @param water_year Logical. If \code{TRUE}, uses water years (Oct--Sep). CT is
#' defined on the water year, so this should normally be TRUE.
#' @param save,output_dir,verbose As in the sibling exposure functions.
#'
#' @return A tibble with one row per water supply / park unit containing the
#' selected timing-change magnitude, dominant direction, model-agreement
#' fraction, and standardized magnitude indicators.
#'
#' @references
#' Stewart, I.T., Cayan, D.R., & Dettinger, M.D. (2005). Changes toward earlier
#' streamflow timing across western North America. J. Climate, 18, 1136-1155.
#' Dudley, R.W., Hodgkins, G.A., McHale, M.R., Kolian, M.J., & Renard, B. (2017).
#' Trends in snowmelt-related streamflow timing in the conterminous United
#' States. J. Hydrology, 547, 208-221.
calc_runoff_timing_exposure <- function(water_supplies,
                                        runoff_data = "data/rapid/inputs/loca2_runoff_water_supplies/",
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
    message("Starting runoff timing (center-of-timing) exposure calculation...")
    start_time <- Sys.time()
  }
  
  hist_start_year <- as.integer(historical[1])
  hist_end_year <- as.integer(historical[2])
  future_start_year <- as.integer(future[1])
  future_end_year <- as.integer(future[2])
  
  # Day-of-water-year at each calendar month's midpoint 
  # Water year starts Oct 1. Fixed (non-leap) month lengths are adequate: CT is
  # a within-year centroid and the sub-day leap error cancels in the future-
  # minus-historical difference.
  .wy_month_order <- c(10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  .mlen <- c(31, 30, 31, 31, 28, 31, 30, 31, 30, 31, 31, 30)
  dowy_lookup <- stats::setNames(cumsum(c(0, head(.mlen, -1))) + .mlen / 2,
                                 .wy_month_order)
  angle_lookup <- 2 * pi * dowy_lookup / 365
  
  if (hindcast == FALSE) {
    files <- list.files(runoff_data, full.names = TRUE) %>%
      grep("2015-2100|training", ., value = TRUE)
  } else {
    files <- list.files(runoff_data, full.names = TRUE) %>%
      grep("training", ., invert = TRUE, value = TRUE)
  }
  if (length(files) == 0) stop("No processed files found. Check input data and paths.")
  
  if (verbose) {
    message("Processing ", length(files), " runoff files")
    pb <- progress::progress_bar$new(
      format = "Processing runoff files [:bar] :percent (:current/:total) ETA: :eta",
      total = length(files), clear = FALSE, width = 80)
  }
  
  process_raw_grid_data <- function(f) {
    raw <- arrow::read_parquet(f) %>%
      dplyr::filter(wsd_source_id %in% water_supplies$wsd_source_id)
    
    # monthly total per grid cell, then average across grid-cell neighborhood
    monthly_points <- data.table::as.data.table(raw)[
      , `:=`(
        year = if (water_year) {
          lubridate::year(date) + data.table::fifelse(lubridate::month(date) >= 10, 1L, 0L)
        } else {
          lubridate::year(date)
        },
        month = lubridate::month(date)
      )
    ][
      (year >= hist_start_year & year <= hist_end_year) |
        (year >= future_start_year & year <= future_end_year)
    ][
      , .(grid_monthly_runoff = if (all(is.na(value))) NA_real_ else sum(value, na.rm = TRUE)),
      by = .(wsd_source_id, year, month, park_unit, nc_file, cell_lon, cell_lat)
    ][
      , .(monthly_runoff = mean(grid_monthly_runoff, na.rm = TRUE)),
      by = .(year, month, wsd_source_id, park_unit, nc_file)
    ]
    
    # Circular center of timing within each year; require all 12 months.
    monthly_points[
      ,
      {
        complete <- length(unique(month)) == 12L && all(!is.na(monthly_runoff))
        tot <- if (complete) sum(monthly_runoff) else NA_real_
        x <- if (!is.na(tot) && tot > 0) {
          sum(cos(angle_lookup[as.character(month)]) * monthly_runoff) / tot
        } else NA_real_
        y <- if (!is.na(tot) && tot > 0) {
          sum(sin(angle_lookup[as.character(month)]) * monthly_runoff) / tot
        } else NA_real_
        .(
          timing_angle = if (is.na(x) || is.na(y)) NA_real_ else atan2(y, x) %% (2 * pi),
          timing_x = x,
          timing_y = y,
          n_months = if (complete) 12L else sum(!is.na(monthly_runoff))
        )
      },
      by = .(year, wsd_source_id, park_unit, nc_file)
    ][n_months == 12L]
  }
  
  runoff_data_list <- purrr::map(files, ~{
    if (verbose) pb$tick()
    process_raw_grid_data(.x)
  })
  
  combined_data <- dplyr::bind_rows(runoff_data_list)
  data.table::setDT(combined_data)
  
  suppressMessages(suppressWarnings({
    combined_data[
      nc_file == "Jay_Alder_runoff_MWBM_LOCA2_training_195001-201412.nc",
      nc_file := "Jay_Alder_runoff_MWBM_LOCA2.historical.training.1950-2100.nc"
    ]
    combined_data[, c("var", "model", "ssp", "run", "nc") := data.table::tstrsplit(nc_file, "\\.")]
    combined_data[!is.na(run), c("run", "time") := data.table::tstrsplit(run, "_", fixed = TRUE)]
  }))
  combined_data[, c("var", "nc") := NULL]
  
  # Circular mean CT over years. Averaging annual x/y components gives each
  # water year equal influence.
  summarize_ct <- function(df) {
    df %>%
      dplyr::group_by(wsd_source_id, park_unit) %>%
      dplyr::summarize(timing_x = mean(timing_x, na.rm = TRUE),
                       timing_y = mean(timing_y, na.rm = TRUE),
                       timing_angle = atan2(timing_y, timing_x) %% (2 * pi),
                       range = paste0(min(year), "-", max(year)),
                       .groups = "drop")
  }
  
  if (hindcast == TRUE) {
    ssps <- combined_data %>% dplyr::filter(ssp != "historical") %>%
      dplyr::distinct(model, ssp, run)
    runoff_change_sub <- vector("list", nrow(ssps))
    
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
                      range_historical = range) %>%
        dplyr::select(-timing_x, -timing_y)
      future_summary <- combined_data %>%
        dplyr::filter(model == ssp_sub$model, ssp == ssp_sub$ssp, run == ssp_sub$run,
                      year >= future_start_year, year <= future_end_year) %>%
        summarize_ct() %>%
        dplyr::rename(timing_angle_future = timing_angle,
                      range_future = range) %>%
        dplyr::select(-timing_x, -timing_y)
      
      runoff_change_sub[[i]] <- future_summary %>%
        dplyr::inner_join(historical_summary, by = c("wsd_source_id", "park_unit")) %>%
        dplyr::mutate(model = ssp_sub$model, ssp = ssp_sub$ssp, run = ssp_sub$run,
                      timing_change = atan2(sin(timing_angle_future - timing_angle_historical),
                                            cos(timing_angle_future - timing_angle_historical)) * 365 / (2 * pi),
                      timing_change_magnitude = abs(timing_change))
    }
    runoff_change <- dplyr::bind_rows(runoff_change_sub)
    
  } else {
    ssps <- combined_data %>% dplyr::filter(ssp != "training") %>%
      dplyr::distinct(model, ssp, run)
    training_hist <- combined_data %>%
      dplyr::filter(ssp == "training", year >= hist_start_year, year <= hist_end_year)
    runoff_change_sub <- vector("list", nrow(ssps))
    
    for (i in seq_len(nrow(ssps))) {
      suppressMessages(suppressWarnings({
        ssp_sub <- ssps[i, ]
        model_hist <- combined_data %>%
          dplyr::filter(model == ssp_sub$model, ssp == ssp_sub$ssp, run == ssp_sub$run,
                        year >= hist_start_year, year <= hist_end_year,
                        !year %in% training_hist$year)
        historical_summary <- summarize_ct(dplyr::bind_rows(training_hist, model_hist)) %>%
          dplyr::rename(timing_angle_historical = timing_angle,
                        range_historical = range) %>%
          dplyr::select(-timing_x, -timing_y)
        future_summary <- combined_data %>%
          dplyr::filter(model == ssp_sub$model, ssp == ssp_sub$ssp, run == ssp_sub$run,
                        year >= future_start_year, year <= future_end_year) %>%
          summarize_ct() %>%
          dplyr::rename(timing_angle_future = timing_angle,
                        range_future = range) %>%
          dplyr::select(-timing_x, -timing_y)
        
        runoff_change_sub[[i]] <- future_summary %>%
          dplyr::inner_join(historical_summary, by = c("wsd_source_id", "park_unit")) %>%
          dplyr::mutate(model = ssp_sub$model, ssp = ssp_sub$ssp, run = ssp_sub$run,
                        timing_change = atan2(sin(timing_angle_future - timing_angle_historical),
                                              cos(timing_angle_future - timing_angle_historical)) * 365 / (2 * pi),
                        timing_change_magnitude = abs(timing_change))
      }))
    }
    runoff_change <- dplyr::bind_rows(runoff_change_sub)
  }
  
  suppressMessages(suppressWarnings({
    timing_summary <- runoff_change %>%
      dplyr::group_by(park_unit, wsd_source_id) %>%
      dplyr::summarize(n_valid = sum(!is.na(timing_change)),
                       n_earlier = sum(timing_change < 0, na.rm = TRUE),
                       n_later = sum(timing_change > 0, na.rm = TRUE),
                       runoff_timing_change_direction = dplyr::case_when(
                         n_valid == 0 ~ NA_character_,
                         n_earlier > n_later ~ "earlier",
                         n_later > n_earlier ~ "later",
                         TRUE ~ "tie"
                       ),
                       selected_timing_change = dplyr::case_when(
                         runoff_timing_change_direction == "earlier" ~
                           stats::quantile(timing_change, 0.10, na.rm = TRUE),
                         runoff_timing_change_direction == "later" ~
                           stats::quantile(timing_change, 0.90, na.rm = TRUE),
                         TRUE ~ NA_real_
                       ),
                       raw_runoff_timing_change_magnitude = abs(selected_timing_change),
                       raw_runoff_timing_direction_agreement_fraction = dplyr::if_else(
                         n_valid > 0, max(n_earlier, n_later) / n_valid, NA_real_
                       ),
                       .groups = "drop"
      ) %>%
      dplyr::select(-n_valid, -n_earlier, -n_later, -selected_timing_change)
  }))
  
  final_results <- timing_summary %>%
    # Larger magnitude = greater departure from historical = more exposed
    add_standardized_indicators("raw_runoff_timing_change_magnitude", invert = FALSE) %>%
    dplyr::select(wsd_source_id, park_unit, dplyr::everything())
  
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    filename <- file.path(output_dir, paste0("runoff_timing_exposure_", Sys.Date(), ".csv"))
    readr::write_csv(final_results, filename)
    if (verbose) message("Results saved to: ", filename)
  }
  
  if (verbose) {
    total_time <- Sys.time() - start_time
    message("\nRunoff timing (CT) exposure completed in ",
            round(total_time, 2), " ", attr(total_time, "units"))
  }
  
  final_results
}
