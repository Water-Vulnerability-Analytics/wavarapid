#' Calculate Historical SWE Sensitivity Indicator
#'
#' Calculates a historical SWE sensitivity indicator from processed LOCA2 MWBM
#' SWE data using only the true historical SWE time series. Annual SWE is
#' aggregated from grid cells and used to estimate historical SWE trends with
#' Sen's slope and the Mann-Kendall test for each water supply.
#'
#' @param water_supplies sf object with water supply locations.
#' @param swe_data character, directory containing processed LOCA2 MWBM SWE files.
#'   Defaults to `wava_data_path("rapid", "inputs", "loca2_swe_water_supplies/")`.
#' @param historical Character vector of length 2 giving the start and end years
#'   (in `"YYYY"` format) that define the historical period.
#' @param water_year T/F, if TRUE, uses water-year bounds (Oct 1 - Sep 30) before annual aggregation.
#' @param save logical, whether to save output. Default is `FALSE`.
#' @param output_dir character, directory path for saving files. Required if
#'   `save = TRUE`.
#' @param verbose logical, whether to show progress messages. Default is `TRUE`.
#'
#' @return tibble with one row per water supply and historical SWE sensitivity
#'   indicators.
#'
#' @export
calc_swe_sensitivity <- function(water_supplies,
                                 swe_data = wava_data_path("rapid", "inputs", "loca2_swe_water_supplies/"),
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
    message("Starting historical SWE sensitivity calculation...")
  }
  
  files <- grep("Jay", list.files(swe_data, full.names = TRUE, pattern = "snow"), value = TRUE, invert = FALSE)
  
  if (length(files) == 0) {
    stop("No processed SWE files found. Check your input data and file paths.")
  }
  
  if (verbose) {
    message("Processing ", length(files), " SWE files")
    pb <- progress::progress_bar$new(
      format = "Processing SWE files [:bar] :percent (:current/:total) ETA: :eta",
      total = length(files),
      clear = FALSE,
      width = 80
    )
  }
  
  process_raw_grid_data <- function(f) {
    raw <- arrow::read_parquet(f) %>%
      dplyr::filter(wsd_source_id %in% water_supplies$wsd_source_id)
    
    data.table::as.data.table(raw)[
      , year := if (water_year) {
        year(date) + fifelse(month(date) >= 10, 1L, 0L)
      } else {
        year(date)
      }
    ][
      year >= hist_start_year & year <= hist_end_year
    ][
      # mean SWE across grid cells per month
      , .(mean_swe = if (any(is.finite(value))) {
        mean(value[is.finite(value)])
      } else {
        NA_real_
      }),
      by = .(wsd_source_id, date, year, park_unit, nc_file)
    ][
      # average annual max per water source/model
      , .(max_swe = if (any(is.finite(mean_swe))) {
        max(mean_swe[is.finite(mean_swe)])
      } else {
        NA_real_
      }),
      by = .(year, wsd_source_id, park_unit, nc_file)
    ]
    
  }
  
  annual_data <- purrr::map(files, ~{
    if (verbose) pb$tick()
    process_raw_grid_data(.x)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::select(wsd_source_id, park_unit, year, max_swe)
  
  if (nrow(annual_data) == 0) {
    stop("No historical SWE records found after filtering.")
  }
  
  final_results <- annual_data %>%
    dplyr::arrange(wsd_source_id, year) %>%
    dplyr::group_by(wsd_source_id, park_unit) %>%
    dplyr::summarize(raw_avg_swe = mean(max_swe, na.rm = TRUE),
                     data = list(dplyr::pick(dplyr::everything())),
                     .groups = "drop") %>%
    mutate(n_years = map_int(data, ~ sum(!is.na(.x$max_swe))),
           sen = map(data, ~ {
             x <- .x$max_swe
             if (sum(!is.na(x)) >= 3) sens.slope(x) else NULL
           }),
           mk = map(data, ~ {
             x <- .x$max_swe
             if (sum(!is.na(x)) >= 3) mk.test(x) else NULL
           }),
           raw_slope_swe = map_dbl(sen, ~ if (is.null(.x)) NA_real_ else .x$estimates[[1]]),
           mk_tau = map_dbl(mk, ~ if (is.null(.x)) NA_real_ else .x$estimates[[1]]),
           mk_p_value = map_dbl(mk, ~ if (is.null(.x)) NA_real_ else .x$p.value),
           mk_significant = mk_p_value < 0.05) %>%
    select(-data, -sen, -mk, -mk_p_value, -mk_tau) %>%
    mutate(trend_direction = case_when(is.na(raw_slope_swe) ~ "no trend",
                                       raw_slope_swe > 0 ~ "increasing",
                                       raw_slope_swe < 0 ~ "decreasing",
           TRUE ~ "no change"),
           raw_percent_slope_swe = ifelse(is.na(raw_slope_swe) | is.na(raw_avg_swe) | raw_avg_swe == 0,
                                          NA_real_,
                                          (raw_slope_swe / raw_avg_swe) * 100)) %>%
    # add standardized indicators using helper function
    add_standardized_indicators("raw_slope_swe", invert = TRUE) %>%
    add_standardized_indicators("raw_avg_swe", invert = TRUE) %>%
    add_standardized_indicators("raw_percent_slope_swe", invert = TRUE) %>%
    # reorder columns following standardized pattern
    select(wsd_source_id, park_unit, everything())
  
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(output_dir, paste0("swe_sensitivity_", Sys.Date(), ".csv"))
    readr::write_csv(final_results, filename)
    
    if (verbose) {
      message("Results saved to: ", filename)
    }
  }
  
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message(
      "\nHistorical SWE sensitivity processing completed in ",
      round(total_time, 2), " ", attr(total_time, "units")
    )
  }
  
  return(final_results)
}


