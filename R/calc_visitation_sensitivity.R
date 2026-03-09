#' Calculate Historical Visitation Sensitivity Indicator
#'
#' Calculates a historical visitation sensitivity indicator for park water supplies
#' using annual NPS visitation totals. For each park unit, annual visitation is
#' summarized over the specified historical period and used to estimate historical
#' visitation trends with Sen's slope and the Mann-Kendall test.
#'
#' Results are merged with the provided `water_supplies` dataset.
#'
#' @param water_supplies A data frame containing at least a `park_unit` column,
#'   representing unique park identifiers for which visitation sensitivity
#'   indicators will be calculated.
#' @param start_year Numeric. First year of visitation data to request.
#'   Defaults to `2010`.
#' @param end_year Numeric. Last year of visitation data to request.
#'   Defaults to `2024`.
#' @param save Logical. If `TRUE`, results are saved to a file in `output_dir`.
#'   Defaults to `FALSE`.
#' @param output_dir Optional character string specifying the directory where
#'   results will be saved if `save = TRUE`.
#' @param verbose Logical. If `TRUE`, progress messages are printed to the console.
#'   Defaults to `TRUE`.
#'
#' @return
#' A tibble combining the original `water_supplies` data with historical visitation
#' sensitivity indicators.
#'
calc_visitation_sensitivity <- function(water_supplies,
                                        start_year = 1990,
                                        end_year = 2024,
                                        save = FALSE,
                                        output_dir = NULL,
                                        verbose = TRUE) {
  
  stopifnot(
    "water_supplies must contain park_unit" = "park_unit" %in% names(water_supplies),
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir)
  )
  
  start_time <- Sys.time()
  
  parks <- unique(water_supplies$park_unit)
  
  if (verbose) {
    message("Starting historical visitation sensitivity calculation...")
    message("Processing ", length(parks), " park units")
  }
  
  visitation_indicators <- function(park_unit) {
    
    visitation <- getUnitVisitation(units = park_unit,
      startYear = start_year,
      endYear = end_year)
    
    annual_visitation <- visitation %>%
      dplyr::group_by(Year) %>%
      dplyr::summarize(total_visitation = sum(RecreationVisitors, na.rm = TRUE)) %>%
      ungroup() %>%
      dplyr::arrange(Year)
    
    n_years <- sum(!is.na(annual_visitation$total_visitation))
    
    sen <- if (n_years >= 3) trend::sens.slope(annual_visitation$total_visitation) else NULL
    mk  <- if (n_years >= 3) trend::mk.test(annual_visitation$total_visitation) else NULL
    
    raw_sen_slope <- if (is.null(sen)) NA_real_ else sen$estimates[[1]]
    mk_tau <- if (is.null(mk)) NA_real_ else mk$estimates[[1]]
    mk_p_value <- if (is.null(mk)) NA_real_ else mk$p.value
    mk_significant <- mk_p_value < 0.05
    
    trend_direction <- dplyr::case_when(
      is.na(raw_sen_slope) ~ NA_character_,
      raw_sen_slope > 0 ~ "increasing",
      raw_sen_slope < 0 ~ "decreasing",
      TRUE ~ "no change"
    )
    
    tibble::tibble(park_unit = park_unit,
                   raw_avg_visitation = mean(annual_visitation$total_visitation, na.rm = TRUE),
                   n_years = n_years,
                   time_range = paste(min(annual_visitation$Year, na.rm = TRUE),
                                      max(annual_visitation$Year, na.rm = TRUE),
                                      sep = "-"),
                   raw_sen_slope = raw_sen_slope,
                   mk_p_value = mk_p_value,
                   mk_significant = mk_significant,
                   trend_direction = trend_direction,
                   peak_year = annual_visitation %>%
                     dplyr::filter(total_visitation == max(total_visitation, na.rm = TRUE)) %>%
                     dplyr::pull(Year) %>%
                     paste(collapse = ", "))
  }
  
  indicators <- parks %>%
    purrr::map(purrr::possibly(~ visitation_indicators(.x), otherwise = NULL)) %>%
    purrr::keep(~ !is.null(.x)) %>%
    dplyr::bind_rows()
  
  final_results <- indicators %>%
    add_standardized_indicators("raw_slope_visitation") %>%
    add_standardized_indicators("raw_avg_visitation") %>%
    dplyr::left_join(water_supplies, ., by = "park_unit")
  
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(output_dir, paste0("visitation_sensitivity_", Sys.Date(), ".csv"))
    readr::write_csv(final_results, filename)
    
    if (verbose) {
      message("Results saved to: ", filename)
    }
  }
  
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message(
      "\nHistorical visitation sensitivity processing completed in ",
      round(total_time, 2), " ", attr(total_time, "units")
    )
  }
  
  return(final_results)
}