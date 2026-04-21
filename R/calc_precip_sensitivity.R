#' Calculate Historical Precipitation Sensitivity Indicator
#'
#' Calculates a historical precip sensitivity indicator from processed LOCA2 MWBM
#' precip data using only the true historical precip time series. Annual precip is
#' aggregated from grid cells and used to estimate historical precip trends with
#' Sen's slope and the Mann-Kendall test for each water supply.
#'
#' @param water_supplies sf object with water supply locations.
#' @param precip_data character, directory containing processed LOCA2 MWBM precip files.
#'   Defaults to `"data/rapid/inputs/loca2_precipitation_water_supplies/"`.
#' @param historical Character vector of length 2 giving the start and end years
#'   (in `"YYYY"` format) that define the historical period.
#' @param water_year T/F, if TRUE, uses water-year bounds (Oct 1 - Sep 30) before annual aggregation.
#' @param save logical, whether to save output. Default is `FALSE`.
#' @param output_dir character, directory path for saving files. Required if
#'   `save = TRUE`.
#' @param verbose logical, whether to show progress messages. Default is `TRUE`.
#'
#' @return tibble with one row per water supply and historical precip sensitivity
#'   indicators.
#'
calc_precip_sensitivity <- function(water_supplies,
                                    precip_data = "data/rapid/inputs/loca2_precipitation_water_supplies/",
                                    historical,
                                    water_year,
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
  
  # use raw precip data (denoted with "livneh")
  files <- grep("livneh", list.files(precip_data, full.names = TRUE, pattern = "precip"), value = TRUE, invert = FALSE)
  
  files_dt <- data.table(file = files)
  
  # extract years
  files_dt[, `:=`(
    year_single = as.integer(str_extract(file, "(?<!\\d)(\\d{4})(?=\\.nc)")),
    year_start = as.integer(str_extract(file, "(\\d{4})(?=-\\d{4})")),
    year_end = as.integer(str_extract(file, "(?<=-)(\\d{4})(?=\\.)"))
  )]
  
  # fill single-year files into start/end
  files_dt[!is.na(year_single), `:=`(
    year_start = year_single,
    year_end = year_single
  )]
  
  if(water_year == FALSE){
    # keep files overlapping either period
    files_subset <- files_dt[
      (year_end >= hist_start_year & year_start <= hist_end_year) 
    ]
  } else {
    # keep files overlapping either period (with wyear accounted for)
    files_subset <- files_dt[
      (year_end >= hist_start_year - 1 & year_start <= hist_end_year) 
    ]
  }
  
  files <- files_subset$file
  
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
      ,       year := if (water_year) {
        year(date) + fifelse(month(date) >= 10, 1L, 0L)
      } else {
        year(date)
      }
    ][
      year >= hist_start_year & year <= hist_end_year
    ][
      , .(grid_precip = sum(value, na.rm = TRUE)),
      by = .(wsd_source_id, year, park_unit, cell_lon, cell_lat)
    ][
      , .(total_precip = mean(grid_precip, na.rm = TRUE)),
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
                     n_years = sum(!is.na(total_precip)),
                     data = list(dplyr::pick(dplyr::everything())),
                     .groups = "drop") %>%
    dplyr::mutate(n_years = map_int(data, ~ sum(!is.na(.x$total_precip))),
                  sen = map(data, ~ {
                    x <- .x$total_precip
                    if (sum(!is.na(x)) >= 3) sens.slope(x) else NULL
                  }),
                  mk = map(data, ~ {
                    x <- .x$total_precip
                    if (sum(!is.na(x)) >= 3) mk.test(x) else NULL
                  }),
                  raw_slope_precip = map_dbl(sen, ~ if (is.null(.x)) NA_real_ else .x$estimates[[1]]),
                  mk_tau = map_dbl(mk, ~ if (is.null(.x)) NA_real_ else .x$estimates[[1]]),
                  mk_p_value = map_dbl(mk, ~ if (is.null(.x)) NA_real_ else .x$p.value),
                  mk_significant = mk_p_value < 0.05) %>%
    dplyr::select(-data, -sen, -mk, -mk_p_value, -mk_tau) %>%
    dplyr::mutate(trend_direction = case_when(is.na(raw_slope_precip) ~ NA_character_,
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
    select(wsd_source_id, park_unit, everything())
  
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
