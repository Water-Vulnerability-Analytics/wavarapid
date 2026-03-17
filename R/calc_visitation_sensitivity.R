#' Calculate Historical Visitation Sensitivity Indicator
#'
#' Calculates historical visitation sensitivity indicators for park water supplies
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
#'   Defaults to `1990`.
#' @param end_year Numeric. Last year of visitation data to request.
#'   Defaults to `2024`.
#' @param save Logical. If `TRUE`, results are saved to a file in `output_dir`.
#' @param output_dir Optional character string specifying the directory where
#'   results will be saved if `save = TRUE`.
#' @param verbose Logical. If `TRUE`, progress messages are printed to the console.
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
  
  parks <- water_supplies %>%
    dplyr::distinct(park_unit) %>%
    dplyr::pull(park_unit)
  
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
      dplyr::summarize(total_visitation = sum(RecreationVisitors, na.rm = TRUE),
                       .groups = "drop") %>%
      dplyr::arrange(Year)
    
    n_years <- sum(!is.na(annual_visitation$total_visitation))
    raw_avg_visitation <- mean(annual_visitation$total_visitation, na.rm = TRUE)
    
    sen <- if (n_years >= 3) {
      trend::sens.slope(annual_visitation$total_visitation)
    } else {
      NULL
    }
    
    mk <- if (n_years >= 3) {
      trend::mk.test(annual_visitation$total_visitation)
    } else {
      NULL
    }
    
    raw_sen_slope <- if (is.null(sen)) NA_real_ else as.numeric(sen$estimates[[1]])
    mk_tau <- if (is.null(mk)) NA_real_ else as.numeric(mk$estimates[[1]])
    mk_p_value <- if (is.null(mk)) NA_real_ else mk$p.value
    mk_significant <- ifelse(is.na(mk_p_value), NA, mk_p_value < 0.05)
    
    # annual percent change in visitation relative to a park’s average visitation over the analysis period
    # calling "scaled_trend"
    percent_sen_slope <- ifelse(is.na(raw_sen_slope) | is.na(raw_avg_visitation) | raw_avg_visitation == 0,
                                NA_real_,
                                (raw_sen_slope / raw_avg_visitation) * 100)
    
    trend_direction <- dplyr::case_when(is.na(raw_sen_slope) ~ NA_character_,
                                        raw_sen_slope > 0 ~ "increasing",
                                        raw_sen_slope < 0 ~ "decreasing",
                                        TRUE ~ "no change")
    
    peak_year <- annual_visitation %>%
      dplyr::filter(total_visitation == max(total_visitation, na.rm = TRUE)) %>%
      dplyr::pull(Year) %>%
      paste(collapse = ", ")
    
    tibble::tibble(park_unit = park_unit,
                   raw_avg_visitation = raw_avg_visitation,
                   raw_visitation_slope = raw_sen_slope,
                   raw_scaled_trend_visitation = percent_sen_slope,
                   mk_tau = mk_tau,
                   mk_p_value = mk_p_value,
                   mk_significant = mk_significant,
                   trend_direction = trend_direction,
                   n_years = n_years,
                   time_range = paste(min(annual_visitation$Year, na.rm = TRUE),
                                      max(annual_visitation$Year, na.rm = TRUE),
                                      sep = "-"),
                   peak_year = peak_year)
    
  }
  
  indicators <- parks %>%
    purrr::map(purrr::possibly(
      ~ visitation_indicators(.x), otherwise = NULL)) %>%
    purrr::keep(~ !is.null(.x)) %>%
    dplyr::bind_rows()
  
  final_results <- water_supplies %>%
    dplyr::select(park_unit, wsd_source_id) %>%
    dplyr::left_join(indicators, by = "park_unit") %>%
    add_standardized_indicators("raw_visitation_slope") %>%
    # annual percent change in visitation relative to a park’s average visitation 
    # over the analysis period
    add_standardized_indicators("raw_scaled_trend_visitation") %>%
    add_standardized_indicators("raw_avg_visitation") %>% 
    st_drop_geometry() %>%
    select(wsd_source_id, park_unit, everything()) 
  
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(
      output_dir,
      paste0("visitation_sensitivity_", Sys.Date(), ".csv")
    )
    
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