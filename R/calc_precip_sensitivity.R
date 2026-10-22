#' @title Precipitation sensitivity: historical trend in annual precipitation
#'
#' @description
#' Calculates the historical trend in annual precipitation at each water
#' supply from the Jay Alder MWBM training run (the file in `precip_data`
#' with "jay_alder" in its name, lowercase). Monthly precipitation is
#' averaged over each site's grid cells and summed to annual totals. The
#' trend is Sen's slope ([trend::sens.slope()]), tested with Mann-Kendall
#' ([trend::mk.test()]); sites with fewer than 3 years get `NA`.
#'
#' @param water_supplies sf object. The water supplies, with a
#'   `wsd_source_id` column, e.g. from [wavacore::get_wsd()]. Only their IDs
#'   are used, to pick rows from the processed file. Required.
#' @param precip_data String. The folder of processed precipitation parquet
#'   files from [process_loca2_mwbm()]. Only the Jay Alder training file is
#'   read; stops with an error if there isn't one. Default:
#'   `rapid/inputs/loca2_precipitation_water_supplies` in the shared data
#'   folder,
#'   see [wavacore::wava_data_path()].
#' @param historical Character or numeric vector of 2 years, e.g.
#'   `c(1980, 2010)`. The first and last year of the trend period. The
#'   training run covers 1950-2014. Required.
#' @param water_year Logical. How months are grouped into years.
#'   * `TRUE`: water years, October to September, labelled by the year they
#'     end in.
#'   * `FALSE`: calendar years.
#'
#'   Default: `TRUE`.
#' @param save Logical. Whether to write the results to a CSV file.
#'   * `TRUE`: writes `precip_sensitivity_<date>.csv`, named with today's
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
#'   * `TRUE`: shows a progress bar and messages with the run time and the
#'     path of the saved file.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return A tibble, one row per water supply, with `wsd_source_id`,
#'   `park_unit` and:
#'   * `raw_avg_precip`: mean annual precipitation (mm);
#'   * `raw_slope_precip`: Sen's slope (mm per year);
#'   * `raw_percent_slope_precip`: the slope as a percentage of the mean
#'     (% per year);
#'   * `n_years`, `mk_significant` (Mann-Kendall p < 0.05) and
#'     `trend_direction` (`"increasing"`, `"decreasing"` or `"no change"`);
#'   * `norm_*` and `rank_*` versions of the three `raw_*` columns from
#'     [add_standardized_indicators()]. All three are inverted, so lower or
#'     more steeply declining precipitation gets a higher score.
#' @export
#' @seealso [calc_precip_exposure()] for the projected change;
#'   [calc_drought_sensitivity()], which uses the same training file;
#'   [process_loca2_mwbm()], which makes the input files.
#' @examplesIf file.exists(wavacore::wava_data_path("rapid", "inputs", "loca2_precipitation_water_supplies", "jay_alder_pr_MWBM_LOCA2_training_195001-201412.nc.parquet"))
#' ws <- wavacore::get_wsd(park = "ZION")
#' precip_sen <- calc_precip_sensitivity(ws, historical = c(1980, 2010))
#' precip_sen[, c("wsd_source_id", "raw_avg_precip", "raw_percent_slope_precip",
#'                "trend_direction")]

calc_precip_sensitivity <- function(water_supplies,
                                    precip_data = wavacore::wava_data_path("rapid", "inputs", "loca2_precipitation_water_supplies/"),
                                    historical,
                                    water_year = TRUE,
                                    save = FALSE,
                                    output_dir = NULL,
                                    verbose = TRUE) {
  
  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir)
  )
  
  start_time <- Sys.time()
  
  hist_start_year <- as.integer(historical[1])
  hist_end_year <- as.integer(historical[2])
  
  if (verbose) {
    message("Starting historical precip sensitivity calculation...")
  }
  
  files <- grep("jay_alder", list.files(precip_data, full.names = TRUE, pattern = "pr"), value = TRUE, invert = FALSE)
  
  if (length(files) == 0) {
    stop("No processed precip files found. Check your input data and file paths.")
  }
  
  if (verbose) {
    message("Processing ", length(files), " precip files")
    pb <- progress::progress_bar$new(
      format = "Processing precip files [:bar] :percent (:current/:total) ETA: :eta",
      total = length(files),
      clear = FALSE,
      width = 80
    )
  }
  
  process_raw_grid_data <- function(f) {
    raw <- arrow::read_parquet(f) %>%
      dplyr::filter(wsd_source_id %in% water_supplies$wsd_source_id)
    
    data.table::as.data.table(raw)[
      ,  year := if (water_year) {
        lubridate::year(date) + data.table::fifelse(lubridate::month(date) >= 10, 1L, 0L)
      } else {
        lubridate::year(date)
      }
    ][
      year >= hist_start_year & year <= hist_end_year
    ][
      , .(site_precip = if (any(is.finite(value))) {
        mean(value[is.finite(value)])
      } else {
        NA_real_
      }),
      by = .(wsd_source_id, date, year, park_unit)
    ][
      , .(total_precip = if (any(is.finite(site_precip))) {
        sum(site_precip[is.finite(site_precip)])
      } else {
        NA_real_
      }),
      by = .(wsd_source_id, park_unit, year)
    ]
  }
  
  annual_data <- purrr::map(files, ~{
    if (verbose) pb$tick()
    process_raw_grid_data(.x)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(wsd_source_id, park_unit, year) %>%
    dplyr::summarize(total_precip = mean(total_precip, na.rm = TRUE), .groups = "drop")
  
  if (nrow(annual_data) == 0) {
    stop("No historical precip records found after filtering.")
  }
  
  final_results <- annual_data %>%
    dplyr::arrange(wsd_source_id, year) %>%
    dplyr::group_by(wsd_source_id, park_unit) %>%
    dplyr::summarize(raw_avg_precip = mean(total_precip, na.rm = TRUE),
                     data = list(dplyr::pick(dplyr::everything())),
                     .groups = "drop") %>%
    dplyr::mutate(n_years = purrr::map_int(data, ~ sum(!is.na(.x$total_precip))),
           sen = purrr::map(data, ~ {
             x <- .x$total_precip
             if (sum(!is.na(x)) >= 3) trend::sens.slope(x) else NULL
           }),
           mk = purrr::map(data, ~ {
             x <- .x$total_precip
             if (sum(!is.na(x)) >= 3) trend::mk.test(x) else NULL
           }),
           raw_slope_precip = purrr::map_dbl(sen, ~ if (is.null(.x)) NA_real_ else .x$estimates[[1]]),
           mk_tau = purrr::map_dbl(mk, ~ if (is.null(.x)) NA_real_ else .x$estimates[[1]]),
           mk_p_value = purrr::map_dbl(mk, ~ if (is.null(.x)) NA_real_ else .x$p.value),
           mk_significant = mk_p_value < 0.05) %>%
    dplyr::select(-data, -sen, -mk, -mk_p_value, -mk_tau) %>%
    dplyr::mutate(trend_direction = dplyr::case_when(is.na(raw_slope_precip) ~ NA_character_,
                                       raw_slope_precip > 0 ~ "increasing",
                                       raw_slope_precip < 0 ~ "decreasing",
                                       TRUE ~ "no change"),
           raw_percent_slope_precip = ifelse(is.na(raw_slope_precip) | is.na(raw_avg_precip) | raw_avg_precip == 0,
                                             NA_real_,
                                             (raw_slope_precip / raw_avg_precip) * 100)) %>%
    # add standardized indicators using helper function
    add_standardized_indicators("raw_slope_precip", invert = TRUE) %>%
    add_standardized_indicators("raw_avg_precip", invert = TRUE) %>%
    add_standardized_indicators("raw_percent_slope_precip", invert = TRUE) %>%
    # reorder columns following standardized pattern
    dplyr::select(wsd_source_id, park_unit, dplyr::everything())
  
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(output_dir, paste0("precip_sensitivity_", Sys.Date(), ".csv"))
    readr::write_csv(final_results, filename)
    
    if (verbose) {
      message("Results saved to: ", filename)
    }
  }
  
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message(
      "\nHistorical precip sensitivity processing completed in ",
      round(total_time, 2), " ", attr(total_time, "units")
    )
  }
  
  return(final_results)
}
