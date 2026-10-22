#' @title Runoff timing exposure: projected shift in the center of timing
#'
#' @description
#' Calculates how much the timing of annual runoff at each water supply
#' shifts between a historical and a future period, from LOCA2-driven MWBM
#' monthly runoff. Timing is measured with the circular center of timing
#' (CT), the flow-weighted centroid of the year's hydrograph. CT is the
#' standard runoff timing metric in the western US snowmelt literature
#' (Stewart et al. 2005; McCabe & Clark 2005; Regonda et al. 2005; Hidalgo et
#' al. 2009; Dudley et al. 2017); the circular form avoids a jump at the
#' September-October water-year boundary.
#'
#' For each year and model run, each month's runoff total is placed at the
#' angle of the month's midpoint and the months are combined as weighted
#' vectors; CT is the angle of the result, as a day of the water year. Only
#' years with all 12 months are used. The historical and future CT are
#' circular means over their years, and the shift is future minus historical
#' in days (negative = earlier).
#'
#' For each site, the direction (earlier or later) that most model runs show
#' is selected, and the shift of the run closest to the 10th percentile
#' (earlier) or 90th percentile (later) is used. On an even split, the tail
#' with the larger absolute shift is used. The size of that shift is the
#' exposure. A second set (`alt_timing`) uses the run closest to the 90th
#' percentile of the absolute shift, whatever its direction.
#'
#' @param water_supplies sf object. The water supplies, with a
#'   `wsd_source_id` column, e.g. from [wavacore::get_wsd()]. Only their IDs
#'   are used, to pick rows from the processed files. Required.
#' @param runoff_data String. The folder of processed runoff parquet files
#'   from [process_loca2_mwbm()]. Stops with an error if no files are left
#'   after the `hindcast` filter. Default:
#'   `rapid/inputs/loca2_runoff_water_supplies` in the shared data folder,
#'   see [wavacore::wava_data_path()].
#' @param hindcast Logical. Where the historical period comes from.
#'   * `TRUE`: each model's own historical run, so the shift is measured
#'     within each model. Historical years after the historical run ends
#'     (2014) come from the same model's scenario run. The training file is
#'     not read.
#'   * `FALSE`: the Jay Alder MWBM training run, the same baseline for every
#'     model. Historical years after it ends come from each model's scenario
#'     run. Only the training file and the `2015-2100` scenario files are
#'     read.
#'
#'   Default: `TRUE`.
#' @param historical Character or numeric vector of 2 years, e.g.
#'   `c(1980, 2010)`. The first and last year of the historical period.
#'   Required.
#' @param future Character or numeric vector of 2 years, e.g.
#'   `c(2035, 2065)`. The first and last year of the future period. Required.
#' @param water_year Logical. How months are grouped into years. CT is
#'   defined on the water year, so this should normally be `TRUE`.
#'   * `TRUE`: water years, October to September, labelled by the year they
#'     end in.
#'   * `FALSE`: calendar years. The month angles still start in October.
#'
#'   Default: `TRUE`.
#' @param save Logical. Whether to write the results to a CSV file.
#'   * `TRUE`: writes `runoff_timing_exposure_<date>.csv`, named with today's
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
#'   * `TRUE`: shows a progress bar over the files and messages with the run
#'     time and the path of the saved file.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return A tibble, one row per water supply, with `wsd_source_id`,
#'   `park_unit` and:
#'   * `runoff_timing_change_direction`: `"earlier"` or `"later"`;
#'   * `raw_runoff_timing_change_magnitude`: the size of the selected shift,
#'     in days;
#'   * `raw_runoff_timing_direction_agreement_fraction`: the fraction (0-1)
#'     of model runs that shift in that direction;
#'   * `runoff_alt_timing_change_direction`,
#'     `raw_runoff_alt_timing_change_magnitude` and
#'     `raw_runoff_alt_timing_direction_agreement_fraction`: the same for the
#'     run closest to the 90th percentile of the absolute shift;
#'   * `norm_*` and `rank_*` versions of the two magnitudes from
#'     [add_standardized_indicators()] (a larger shift gets a higher score).
#' @references
#' Stewart, I.T., Cayan, D.R., & Dettinger, M.D. (2005). Changes toward
#'   earlier streamflow timing across western North America. J. Climate, 18,
#'   1136-1155.
#'
#' Dudley, R.W., Hodgkins, G.A., McHale, M.R., Kolian, M.J., & Renard, B.
#'   (2017). Trends in snowmelt-related streamflow timing in the conterminous
#'   United States. J. Hydrology, 547, 208-221.
#' @export
#' @seealso [calc_runoff_exposure()] for the change in amount;
#'   [calc_precip_timing_exposure()] for precipitation timing;
#'   [process_loca2_mwbm()], which makes the input files.
#' @examplesIf interactive() && dir.exists(wavacore::wava_data_path("rapid", "inputs", "loca2_runoff_water_supplies"))
#' # Reads every processed runoff file (about 4 GB), so it takes several
#' # minutes even for one park
#' ws <- wavacore::get_wsd(park = "ROMO")
#' timing_exp <- calc_runoff_timing_exposure(ws,
#'                                           historical = c(1980, 2010),
#'                                           future = c(2035, 2065))
#' timing_exp

calc_runoff_timing_exposure <- function(water_supplies,
                                        runoff_data = wavacore::wava_data_path("rapid", "inputs", "loca2_runoff_water_supplies/"),
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
  
  hist_start_year   <- as.integer(historical[1])
  hist_end_year     <- as.integer(historical[2])
  future_start_year <- as.integer(future[1])
  future_end_year   <- as.integer(future[2])
  
  # Water year starts Oct 1. Fixed (non-leap) month lengths are adequate: CT is
  # a within-year centroid and the sub-day leap error cancels in the future-
  # minus-historical difference.
  .wy_month_order <- c(10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  .mlen           <- c(31, 30, 31, 31, 28, 31, 30, 31, 30, 31, 31, 30)
  dowy_lookup     <- stats::setNames(cumsum(c(0, utils::head(.mlen, -1))) + .mlen / 2,
                                     .wy_month_order)
  angle_lookup    <- 2 * pi * dowy_lookup / 365
  
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
      # First average all available neighborhood cells for each time step.
      , .(site_runoff = if (any(is.finite(value))) {
        mean(value[is.finite(value)])
      } else {
        NA_real_
      }),
      by = .(wsd_source_id, date, year, month, park_unit, nc_file)
    ][
      # Then aggregate the available site values within the month.
      , .(monthly_runoff = if (any(is.finite(site_runoff))) {
        sum(site_runoff[is.finite(site_runoff)])
      } else {
        NA_real_
      }),
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
    nearest_quantile <- function(x, probability) {
      x <- x[is.finite(x)]
      if (!length(x)) return(NA_real_)
      target <- stats::quantile(x, probability, names = FALSE, type = 7)
      x[which.min(abs(x - target))]
    }
    # Primary timing indicators: choose the actual member nearest each tail,
    # then choose the direction supported by more models.
    # For an exact 50/50 split, choose the direction with the larger absolute
    # tail change (p10 earlier versus p90 later). Agreement remains 0.5.
    timing_summary <- runoff_change %>%
      dplyr::group_by(park_unit, wsd_source_id) %>%
      dplyr::summarize(n_valid = sum(!is.na(timing_change)),
                       n_earlier = sum(timing_change < 0, na.rm = TRUE),
                       n_later = sum(timing_change > 0, na.rm = TRUE),
                       p10_timing_change = nearest_quantile(timing_change, 0.10),
                       p90_timing_change = nearest_quantile(timing_change, 0.90),
                       runoff_timing_change_direction = dplyr::case_when(n_valid == 0 ~ NA_character_,
                                                                         n_earlier > n_later ~ "earlier",
                                                                         n_later > n_earlier ~ "later",
                                                                         abs(p10_timing_change) >= abs(p90_timing_change) ~ "earlier",
                                                                         TRUE ~ "later"),
                       selected_timing_change = dplyr::case_when(runoff_timing_change_direction == "earlier" ~ p10_timing_change,
                                                                 runoff_timing_change_direction == "later" ~ p90_timing_change,
                                                                 TRUE ~ NA_real_),
                       raw_runoff_timing_change_magnitude = abs(selected_timing_change),
                       raw_runoff_timing_direction_agreement_fraction = dplyr::if_else(n_valid > 0,
                                                                                       max(n_earlier, n_later) / n_valid,
                                                                                       NA_real_),
                       .groups = "drop") %>%
      dplyr::select(-c(n_valid,
                       n_earlier,
                       n_later,
                       p10_timing_change,
                       p90_timing_change,
                       selected_timing_change))
    
    # Alternative timing indicators (alt_timing): rank model realizations by
    # absolute timing change, select the actual realization nearest the p90,
    # retain its earlier/later direction, and calculate the fraction of valid
    # model realizations projecting that same direction.
    alt_timing_summary <- runoff_change %>%
      dplyr::filter(!is.na(timing_change)) %>%
      dplyr::group_by(park_unit, wsd_source_id) %>%
      dplyr::mutate(abs_timing_change = abs(timing_change),
                    alt_timing_p90 = stats::quantile(abs_timing_change, 0.90, na.rm = TRUE),
                    diff_from_alt_timing_p90 = abs(abs_timing_change - alt_timing_p90)) %>%
      dplyr::filter(diff_from_alt_timing_p90 == min(diff_from_alt_timing_p90, na.rm = TRUE)) %>%
      dplyr::slice(1) %>%
      dplyr::mutate(runoff_alt_timing_change_direction = dplyr::case_when(
        timing_change < 0 ~ "earlier",
        timing_change > 0 ~ "later",
        TRUE ~ "no change"),
        raw_runoff_alt_timing_change_magnitude = abs_timing_change) %>%
      dplyr::select(wsd_source_id,
                    park_unit,
                    runoff_alt_timing_change_direction,
                    raw_runoff_alt_timing_change_magnitude) %>%
      dplyr::left_join(runoff_change %>%
                         dplyr::filter(!is.na(timing_change)) %>%
                         dplyr::group_by(park_unit, wsd_source_id) %>%
                         dplyr::summarize(n_valid = dplyr::n(),
                                          n_earlier = sum(timing_change < 0, na.rm = TRUE),
                                          n_later = sum(timing_change > 0, na.rm = TRUE),
                                          n_no_change = sum(timing_change == 0, na.rm = TRUE),
                                          .groups = "drop"),
                       by = c("wsd_source_id", "park_unit")) %>%
      dplyr::mutate(raw_runoff_alt_timing_direction_agreement_fraction = dplyr::case_when(
        runoff_alt_timing_change_direction == "earlier" ~ n_earlier / n_valid,
        runoff_alt_timing_change_direction == "later" ~ n_later / n_valid,
        runoff_alt_timing_change_direction == "no change" ~ n_no_change / n_valid,
        TRUE ~ NA_real_)) %>%
      dplyr::select(-n_valid, -n_earlier, -n_later, -n_no_change)
  }))
  
  final_results <- timing_summary %>%
    dplyr::left_join(alt_timing_summary, by = c("wsd_source_id", "park_unit")) %>%
    # Original timing exposure indicator
    add_standardized_indicators("raw_runoff_timing_change_magnitude", invert = FALSE) %>%
    # Alternative p90 absolute-magnitude timing exposure indicator
    add_standardized_indicators("raw_runoff_alt_timing_change_magnitude", invert = FALSE) %>%
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
