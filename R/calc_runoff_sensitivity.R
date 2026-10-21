#' Calculate Historical Runoff Sensitivity Indicator
#'
#' Calculates a historical runoff sensitivity indicator from processed LOCA2 MWBM
#' runoff data using only the true historical runoff time series. Annual runoff is
#' aggregated from grid cells and used to estimate historical runoff trends with
#' Sen's slope and the Mann-Kendall test for each water supply.
#'
#' @param water_supplies sf object with water supply locations.
#' @param runoff_data character, directory containing processed LOCA2 MWBM runoff files.
#'   Defaults to `wava_data_path("rapid", "inputs", "loca2_runoff_water_supplies/")`.
#' @param historical Character vector of length 2 giving the start and end years
#'   (in `"YYYY"` format) that define the historical period.
#' @param water_year T/F, if TRUE, uses water-year bounds (Oct 1 - Sep 30) before annual aggregation.
#' @param save logical, whether to save output. Default is `FALSE`.
#' @param output_dir character, directory path for saving files. Required if
#'   `save = TRUE`.
#' @param verbose logical, whether to show progress messages. Default is `TRUE`.
#'
#' @return tibble with one row per water supply and historical runoff sensitivity
#'   indicators.
#'
#' @export
calc_runoff_sensitivity <- function(water_supplies,
                                    runoff_data = wavacore::wava_data_path("rapid", "inputs", "loca2_runoff_water_supplies/"),
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
    message("Starting historical runoff sensitivity calculation...")
  }
  
  files <- grep("Jay", list.files(runoff_data, full.names = TRUE, pattern = "runoff"), value = TRUE, invert = FALSE)
  
  if (length(files) == 0) {
    stop("No processed runoff files found. Check your input data and file paths.")
  }
  
  if (verbose) {
    message("Processing ", length(files), " runoff files")
    pb <- progress::progress_bar$new(
      format = "Processing runoff files [:bar] :percent (:current/:total) ETA: :eta",
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
      , .(site_runoff = if (any(is.finite(value))) {
        mean(value[is.finite(value)])
      } else {
        NA_real_
      }),
      by = .(wsd_source_id, date, year, park_unit)
    ][
      , .(total_runoff = if (any(is.finite(site_runoff))) {
        sum(site_runoff[is.finite(site_runoff)])
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
    dplyr::summarize(total_runoff = mean(total_runoff, na.rm = TRUE), .groups = "drop")
  
  if (nrow(annual_data) == 0) {
    stop("No historical runoff records found after filtering.")
  }
  
  final_results <- annual_data %>%
    dplyr::arrange(wsd_source_id, year) %>%
    dplyr::group_by(wsd_source_id, park_unit) %>%
    dplyr::summarize(raw_avg_runoff = mean(total_runoff, na.rm = TRUE),
                     data = list(dplyr::pick(dplyr::everything())),
                     .groups = "drop") %>%
    dplyr::mutate(n_years = purrr::map_int(data, ~ sum(!is.na(.x$total_runoff))),
           sen = purrr::map(data, ~ {
             x <- .x$total_runoff
             if (sum(!is.na(x)) >= 3) trend::sens.slope(x) else NULL
           }),
           mk = purrr::map(data, ~ {
             x <- .x$total_runoff
             if (sum(!is.na(x)) >= 3) trend::mk.test(x) else NULL
           }),
           raw_slope_runoff = purrr::map_dbl(sen, ~ if (is.null(.x)) NA_real_ else .x$estimates[[1]]),
           mk_tau = purrr::map_dbl(mk, ~ if (is.null(.x)) NA_real_ else .x$estimates[[1]]),
           mk_p_value = purrr::map_dbl(mk, ~ if (is.null(.x)) NA_real_ else .x$p.value),
           mk_significant = mk_p_value < 0.05) %>%
    dplyr::select(-data, -sen, -mk, -mk_p_value, -mk_tau) %>%
    dplyr::mutate(trend_direction = dplyr::case_when(is.na(raw_slope_runoff) ~ NA_character_,
                                       raw_slope_runoff > 0 ~ "increasing",
                                       raw_slope_runoff < 0 ~ "decreasing",
                                       TRUE ~ "no change"),
           raw_percent_slope_runoff = ifelse(is.na(raw_slope_runoff) | is.na(raw_avg_runoff) | raw_avg_runoff == 0,
                                             NA_real_,
                                             (raw_slope_runoff / raw_avg_runoff) * 100)) %>%
    # add standardized indicators using helper function
    add_standardized_indicators("raw_slope_runoff", invert = TRUE) %>%
    add_standardized_indicators("raw_avg_runoff", invert = TRUE) %>%
    add_standardized_indicators("raw_percent_slope_runoff", invert = TRUE) %>%
    # reorder columns following standardized pattern
    dplyr::select(wsd_source_id, park_unit, dplyr::everything())
  
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(output_dir, paste0("runoff_sensitivity_", Sys.Date(), ".csv"))
    readr::write_csv(final_results, filename)
    
    if (verbose) {
      message("Results saved to: ", filename)
    }
  }
  
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message(
      "\nHistorical runoff sensitivity processing completed in ",
      round(total_time, 2), " ", attr(total_time, "units")
    )
  }
  
  return(final_results)
}
