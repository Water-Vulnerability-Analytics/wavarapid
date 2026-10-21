#' @title Visitation sensitivity: historical trend in park visitation
#'
#' @description
#' Calculates the trend in annual recreation visits for each park unit in
#' `water_supplies`, from the NPS visitation statistics API through
#' [wavacore::getUnitVisitation()], so it needs internet. The trend is Sen's
#' slope ([trend::sens.slope()]), tested with Mann-Kendall
#' ([trend::mk.test()]), for parks with at least 3 years. A few WSD codes are
#' looked up under a different visitation code (ANAC and NACE as NCPE, NAMA
#' as NCPC, VALR as PERL); the output keeps the WSD code. Parks whose request
#' or calculation fails are dropped without a message and get `NA`.
#'
#' Every water supply gets the values of its park.
#'
#' @param water_supplies Data frame or sf object. The water supplies, with
#'   `wsd_source_id` and `park_unit` columns, e.g. from
#'   [wavacore::get_wsd()]. One request is made per park unit. Stops with an
#'   error if `park_unit` is missing. Required.
#' @param start_year Integer. The first year of visitation to request.
#'   Default: `1990`.
#' @param end_year Integer. The last year of visitation to request.
#'   Default: `2024`.
#' @param save Logical. Whether to write the results to a CSV file.
#'   * `TRUE`: writes `visitation_sensitivity_<date>.csv`, named with today's date, to
#'     `output_dir`, creating the folder if needed. Stops with an error if
#'     `output_dir` is `NULL`.
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
#'   * `TRUE`: prints the number of parks, the run time and the path of the
#'     saved file.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return A tibble (no geometry), one row per water supply, with
#'   `wsd_source_id`, `park_unit` and:
#'   * `raw_avg_visitation`: mean annual visits;
#'   * `raw_visitation_slope`: Sen's slope (visits per year);
#'   * `raw_scaled_trend_visitation`: the slope as a percentage of the mean
#'     (% per year);
#'   * `mk_tau`, `mk_p_value`, `mk_significant` (p < 0.05),
#'     `trend_direction`, `n_years`, `time_range` and `peak_year`;
#'   * `norm_*` and `rank_*` versions of the three `raw_*` columns from
#'     [add_standardized_indicators()]. Faster growth and more visits get a
#'     higher score.
#' @export
#' @seealso [calc_visitation_exposure()] for the projected change;
#'   [wavacore::getUnitVisitation()], which downloads the visitation data.
#' @examplesIf interactive()
#' # Downloads visitation statistics from the NPS API, so it needs internet
#' ws <- wavacore::get_wsd(park = c("ZION", "ROMO"))
#' visitation_sen <- calc_visitation_sensitivity(ws, start_year = 1990, end_year = 2024)
#' visitation_sen[, c("wsd_source_id", "park_unit", "raw_avg_visitation",
#'                    "raw_scaled_trend_visitation")]

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

  # WSD codes that differ from the codes used by visitation statistics.
  visitation_code_crosswalk <- c(
    "ANAC" = "NCPE",
    "NACE" = "NCPE",
    "NAMA" = "NCPC",
    "VALR" = "PERL"
  )
  
  parks <- water_supplies %>%
    dplyr::distinct(park_unit) %>%
    dplyr::pull(park_unit)
  
  if (verbose) {
    message("Starting historical visitation sensitivity calculation...")
    message("Processing ", length(parks), " park units")
  }
  
  visitation_indicators <- function(park_unit) {

    # Translate only the API lookup code; retain the WSD code in the output.
    park_unit_clean <- toupper(trimws(as.character(park_unit)))
    visitation_code <- unname(visitation_code_crosswalk[park_unit_clean])
    if (is.na(visitation_code)) visitation_code <- park_unit_clean

    visitation <- wavacore::getUnitVisitation(units = visitation_code,
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
    sf::st_drop_geometry() %>%
    dplyr::select(wsd_source_id, park_unit, dplyr::everything()) 
  
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
