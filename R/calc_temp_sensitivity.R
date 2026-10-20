#' Calculate Historical Temperature Sensitivity Indicator
#'
#' Calculates a historical temperature sensitivity indicator from processed
#' LOCA2/MWBM temperature data. This follows the precipitation sensitivity helper
#' structure, but uses annual mean temperature and raw temperature slopes in
#' degrees C per year rather than percent-scaled change.
#'
#' @param water_supplies sf object with water supply locations.
#' @param temp_data character, directory containing processed LOCA2/MWBM temperature files.
#'   Defaults to `wava_data_path("rapid", "inputs", "loca2_temperature_water_supplies/")`.
#' @param historical Character vector of length 2 giving the start and end years
#'   (in `"YYYY"` format) that define the historical period.
#' @param water_year T/F, if TRUE, uses water-year bounds (Oct 1 - Sep 30) before annual aggregation.
#' @param save logical, whether to save output. Default is `FALSE`.
#' @param output_dir character, directory path for saving files. Required if
#'   `save = TRUE`.
#' @param verbose logical, whether to show progress messages. Default is `TRUE`.
#'
#' @return tibble with one row per water supply and historical temperature
#'   sensitivity indicators.
#'
#' @export
calc_temp_sensitivity <- function(water_supplies,
                                  temp_data = wava_data_path("rapid", "inputs", "loca2_temperature_water_supplies/"),
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
    message("Starting historical temperature sensitivity calculation...")
  }
  
  all_files <- list.files(temp_data, full.names = TRUE, pattern = "\\.parquet$")
  temp_files <- grep("tas|temp|temperature|Tmin|Tmax", all_files, value = TRUE, ignore.case = TRUE)
  
  if (length(temp_files) == 0) {
    temp_files <- all_files
  }
  
  # Prefer the historical/training products when they are available, matching
  # the precip helper's use of the true historical series.
  files <- grep("Jay|lusu|livneh|training", temp_files, value = TRUE, ignore.case = TRUE)
  
  if (length(files) == 0) {
    files <- grep("historical", temp_files, value = TRUE, ignore.case = TRUE)
  }
  
  if (length(files) == 0) {
    files <- temp_files
  }
  
  if (length(files) == 0) {
    stop("No processed temperature files found. Check your input data and file paths.")
  }
  
  if (verbose) {
    message("Processing ", length(files), " temperature files")
    pb <- progress::progress_bar$new(
      format = "Processing temperature files [:bar] :percent (:current/:total) ETA: :eta",
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
        year(date) + fifelse(month(date) >= 10, 1L, 0L)
      } else {
        year(date)
      }
    ][
      year >= hist_start_year & year <= hist_end_year
    ][
      , .(site_temperature = if (any(is.finite(value))) {
        mean(value[is.finite(value)])
      } else {
        NA_real_
      }),
      by = .(wsd_source_id, date, year, park_unit)
    ][
      , .(mean_temperature = if (any(is.finite(site_temperature))) {
        mean(site_temperature[is.finite(site_temperature)])
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
    dplyr::summarize(mean_temperature = mean(mean_temperature, na.rm = TRUE), .groups = "drop")
  
  if (nrow(annual_data) == 0) {
    stop("No historical temperature records found after filtering.")
  }
  
  final_results <- annual_data %>%
    dplyr::arrange(wsd_source_id, year) %>%
    dplyr::group_by(wsd_source_id, park_unit) %>%
    dplyr::summarize(raw_avg_temperature = mean(mean_temperature, na.rm = TRUE),
                     data = list(dplyr::pick(dplyr::everything())),
                     .groups = "drop") %>%
    mutate(n_years = map_int(data, ~ sum(!is.na(.x$mean_temperature))),
           sen = map(data, ~ {
             x <- .x$mean_temperature
             if (sum(!is.na(x)) >= 3) sens.slope(x) else NULL
           }),
           mk = map(data, ~ {
             x <- .x$mean_temperature
             if (sum(!is.na(x)) >= 3) mk.test(x) else NULL
           }),
           raw_slope_temperature = map_dbl(sen, ~ if (is.null(.x)) NA_real_ else .x$estimates[[1]]),
           mk_tau = map_dbl(mk, ~ if (is.null(.x)) NA_real_ else .x$estimates[[1]]),
           mk_p_value = map_dbl(mk, ~ if (is.null(.x)) NA_real_ else .x$p.value),
           mk_significant = mk_p_value < 0.05) %>%
    select(-data, -sen, -mk, -mk_p_value, -mk_tau) %>%
    mutate(trend_direction = case_when(is.na(raw_slope_temperature) ~ NA_character_,
                                       raw_slope_temperature > 0 ~ "increasing",
                                       raw_slope_temperature < 0 ~ "decreasing",
                                       TRUE ~ "no change")) %>%
    # add standardized indicators using helper function
    add_standardized_indicators("raw_slope_temperature", invert = FALSE) %>%
    add_standardized_indicators("raw_avg_temperature", invert = FALSE) %>%
    # reorder columns following standardized pattern
    select(wsd_source_id, park_unit, everything())
  
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(output_dir, paste0("temperature_sensitivity_", Sys.Date(), ".csv"))
    readr::write_csv(final_results, filename)
    
    if (verbose) {
      message("Results saved to: ", filename)
    }
  }
  
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message(
      "\nHistorical temperature sensitivity processing completed in ",
      round(total_time, 2), " ", attr(total_time, "units")
    )
  }
  
  return(final_results)
}
