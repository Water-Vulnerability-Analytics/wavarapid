#' @title Visitation exposure: projected change in park visitation to 2050
#'
#' @description
#' Projects annual recreation visits and the seasonality of visits for each
#' park unit in `water_supplies` to 2050, from the NPS visitation statistics
#' API through [wavacore::getUnitVisitation()], so it needs internet. A few
#' WSD codes are looked up under a different visitation code (ANAC and NACE
#' as NCPE, NAMA as NCPC, VALR as PERL); the output keeps the WSD code.
#'
#' For each park:
#' 1. A linear trend is fitted to annual visits from `start_year` to 2024.
#'    The 2024 value of that trend is the starting point (`obs_2024`; it is
#'    the fitted value, not the observed count).
#' 2. The park uses a baseline slope, `baseline_growth_rate` times the 2024
#'    value, instead of its own trend when any of these is true: the trend is
#'    negative; the park has timed entry (ROMO and YOSE from 2020; GLAC, ACAD
#'    and HALE from 2021; ZION, ARCH and SHEN from 2022); it has fewer than
#'    `((2024 - start_year) - 1) * 12` monthly records (with the default
#'    `start_year`, 33 years); or its trend is below the baseline slope.
#' 3. Visits are projected to 2050 (not below 0), and the percent change from
#'    2024 is the visitation indicator.
#' 4. Seasonality is the coefficient of variation (CV) of the monthly shares
#'    of each year's visits, leaving out 2020 and 2021. Its linear trend is
#'    projected to 2050; for parks on the baseline slope the CV trend is
#'    multiplied by `baseline_seasonality_shift`. The seasonality indicator is
#'    `100 * (CV 2024 - CV 2050) / CV 2024`, so positive means visits spread
#'    over more of the year.
#'
#' Parks whose request or calculation fails are dropped without a message
#' and get `NA`.
#'
#' @param water_supplies Data frame or sf object. The water supplies, with a
#'   `park_unit` column, e.g. from [wavacore::get_wsd()]. One request is made
#'   per park unit, and the results are joined back onto these rows.
#'   Required.
#' @param start_year Integer. The first year of visitation to request and
#'   fit. It also sets the short-record test (see Description). Default:
#'   `1990`.
#' @param baseline_growth_rate Numeric. The yearly growth, as a fraction of
#'   the 2024 value, used for parks that get the baseline slope. Default:
#'   `0.0087` (0.87% a year, the long-term NPS system-wide trend).
#' @param baseline_seasonality_shift Numeric from 0 to 1. How much of the
#'   seasonality trend parks on the baseline slope keep.
#'   * `1`: they keep their observed seasonality trend (no change).
#'   * Less than 1, e.g. `0.25`: their seasonality trend is cut to that
#'     fraction, assuming parks with limited growth see little change in
#'     seasonality.
#'
#'   Still being worked out. Default: `1`.
#' @param save Logical. Not used yet: the function never writes a file,
#'   whatever this is. Default: `FALSE`.
#' @param output_dir String or `NULL`. Not used yet, as for `save`.
#'   Default: `NULL`.
#' @param verbose Logical. Whether to print progress.
#'   * `TRUE`: prints how many parks were calculated and which parks got the
#'     baseline slope and why.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return `water_supplies` (with its class and columns) with these columns
#'   joined on by `park_unit`: `pred_2050`, `obs_2024`, `peak_year`,
#'   `raw_visitation_pct_change`, `raw_visitation_cv_change`,
#'   `n_months_analyzed`, `time_range`, `timed_entry_year`,
#'   `has_timed_entry`, `original_slope`, `projection_slope`,
#'   `slope_adjustment` (`"none"`, `"baseline_insufficient_data"`,
#'   `"baseline_negative"`, `"baseline_timed_entry"` or
#'   `"trend_below_baseline"`), `baseline_slope`, `negative_trend`, and the
#'   `norm_` and `rank_` versions of the two `raw_` columns from
#'   [add_standardized_indicators()].
#' @export
#' @seealso [calc_visitation_sensitivity()] for the historical trend;
#'   [wavacore::getUnitVisitation()], which downloads the visitation data.
#' @examplesIf interactive()
#' # Downloads visitation statistics from the NPS API, so it needs internet
#' ws <- wavacore::get_wsd(park = c("ZION", "GRSA"))
#' visitation_exp <- calc_visitation_exposure(ws)
#' unique(sf::st_drop_geometry(visitation_exp)[, c("park_unit", "slope_adjustment",
#'                                                 "raw_visitation_pct_change")])

calc_visitation_exposure <- function(water_supplies,
                                     start_year = 1990,
                                     baseline_growth_rate = 0.0087,
                                     baseline_seasonality_shift = 1,
                                     save = FALSE,
                                     output_dir = NULL,
                                     verbose = TRUE) {

  # WSD codes that differ from the codes used by visitation statistics.
  visitation_code_crosswalk <- c(
    "ANAC" = "NCPE",
    "NACE" = "NCPE",
    "NAMA" = "NCPC",
    "VALR" = "PERL"
  )
  
  # timed entry start years
  timed_entry_parks <- dplyr::tribble(~park_unit, ~timed_entry_year,
                               "ROMO", 2020,
                               "YOSE", 2020,
                               "GLAC", 2021,
                               "ACAD", 2021,
                               "HALE", 2021,
                               "ZION", 2022,
                               "ARCH", 2022,
                               "SHEN", 2022)
  
  parks <- unique(water_supplies$park_unit)
  
  visitation_indicators <- function(park_unit){

    # Translate only the API lookup code; retain the WSD code in the output.
    park_unit_clean <- toupper(trimws(as.character(park_unit)))
    visitation_code <- unname(visitation_code_crosswalk[park_unit_clean])
    if (is.na(visitation_code)) visitation_code <- park_unit_clean

    visitation <- wavacore::getUnitVisitation(units = visitation_code,
                                    startYear = start_year,
                                    endYear = 2024)
    
    # check if park has timed entry
    timed_entry_info <- timed_entry_parks %>% dplyr::filter(.data$park_unit == !!park_unit)
    has_timed_entry <- nrow(timed_entry_info) > 0
    timed_entry_year <- if(has_timed_entry) timed_entry_info$timed_entry_year else NA
    
    # calculate annual totals
    annual_visitation <- visitation %>%
      dplyr::group_by(Year) %>%
      dplyr::summarize(total = sum(RecreationVisitors), .groups = "drop")
    
    seasonal_visitation <- visitation %>%
      dplyr::mutate(date = lubridate::ym(paste0(Year, "-", Month)),
             season = dplyr::case_when(Month %in% c(12, 1, 2) ~ "Winter",
                                Month %in% c(3, 4, 5) ~ "Spring",
                                Month %in% c(6, 7, 8) ~ "Summer",
                                Month %in% c(9, 10, 11) ~ "Fall"),
             month_name = month.abb[Month],
             season_month = factor(paste0(season, "_", month_name),
                                   levels = c("Winter_Jan","Winter_Feb","Spring_Mar","Spring_Apr","Spring_May",
                                              "Summer_Jun","Summer_Jul","Summer_Aug",
                                              "Fall_Sep","Fall_Oct","Fall_Nov",
                                              "Winter_Dec")))
    
    # fit model on full data (start_year-2024)
    lm_total <- stats::lm(total ~ Year, data = annual_visitation)
    original_slope <- stats::coef(lm_total)[2]
    intercept <- stats::coef(lm_total)[1]
    
    # calculate baseline growth in absolute visitors per year
    # baseline_slope <- pred_2024 * baseline_growth_rate
    # 
    # pred_2024 <- annual_visitation %>% filter(Year == 2024) %>% pull(total)
    
    pred_2024 <- stats::predict(lm_total, newdata = data.frame(Year = 2024))
    baseline_slope <- pred_2024 * baseline_growth_rate
    
    # check if baseline growth adjustment is needed
    negative_trend <- original_slope < 0
    insufficient_data <- nrow(visitation) < (((2024-start_year) - 1) * 12)
    below_baseline <- original_slope < baseline_slope
    use_baseline <- negative_trend || has_timed_entry || insufficient_data || below_baseline
    
    # grab final slope to use for projection
    # depending on the park's circumstances
    if(use_baseline) {
      projection_slope <- baseline_slope
      slope_adjustment <- dplyr::case_when(
        insufficient_data ~ "baseline_insufficient_data",
        negative_trend    ~ "baseline_negative",
        has_timed_entry   ~ "baseline_timed_entry",
        below_baseline    ~ "trend_below_baseline",
        TRUE              ~ "baseline_other"
      )
    } else {
      projection_slope <- original_slope
      slope_adjustment <- "none"
    }
    
    # project to 2050
    years_forward <- 2050 - 2024
    pred_2050 <- max(0, pred_2024 + projection_slope * years_forward)
    pct_change <- 100*((pred_2050 - pred_2024) / pred_2024)
    
    # other fun facts:
    peak_year <- annual_visitation %>%
      dplyr::filter(total == max(total, na.rm = TRUE)) %>%
      dplyr::pull(Year) %>%
      paste(collapse = ", ")
    
    # projected change in seasonality (exclude COVID years):
    monthly_percentages <- visitation %>%
      dplyr::filter(!Year %in% c(2020, 2021)) %>%
      dplyr::group_by(Year) %>%
      dplyr::mutate(annual_total = sum(RecreationVisitors),
             monthly_percentage = RecreationVisitors / annual_total * 100) %>%
      dplyr::ungroup()
    
    seasonal_variation <- monthly_percentages %>%
      dplyr::group_by(Year) %>%
      dplyr::summarize(cv_monthly = stats::sd(monthly_percentage) / mean(monthly_percentage) * 100,
                .groups = "drop")
    
    lm_cv <- stats::lm(cv_monthly ~ Year, data = seasonal_variation)
    original_cv_slope <- stats::coef(lm_cv)[2]
    
    cv_2024 <- stats::predict(lm_cv, newdata = data.frame(Year = 2024))
    
    # apply appropriate adjustment to seasonality trend (? maybe - or should it always be what's represented in data?)
    if(use_baseline) {
      # assume stable or very modest seasonality change for baseline growth parks
      projection_cv_slope <- original_cv_slope * baseline_seasonality_shift
    } else {
      projection_cv_slope <- original_cv_slope
    }

    cv_2050 <- cv_2024 + (projection_cv_slope * years_forward)
    cv_change <- 100 * ((cv_2024 - cv_2050) / cv_2024)  # positive = decreasing CV (longer visitation window)
    
    # Return raw indicators (scaling happens later)
    dplyr::tibble(
      park_unit = park_unit,
      pred_2050 = pred_2050,
      obs_2024 = pred_2024,
      peak_year = peak_year,
      raw_visitation_pct_change = pct_change,
      raw_visitation_cv_change = cv_change,
      n_months_analyzed = nrow(visitation),
      time_range = paste(min(visitation$Year), max(visitation$Year), sep = "-"),
      timed_entry_year = timed_entry_year,
      has_timed_entry = has_timed_entry,
      original_slope = original_slope,
      projection_slope = projection_slope,
      slope_adjustment = slope_adjustment,
      baseline_slope = baseline_slope,
      negative_trend = negative_trend
    )
  }
  
  indicators <- parks %>%
    purrr::map(purrr::possibly(~ visitation_indicators(.x), otherwise = NULL)) %>%
    purrr::keep(~!is.null(.)) %>%
    dplyr::bind_rows()
  
  norm_indicators <- indicators %>%
    add_standardized_indicators(., prefix = "raw_visitation") %>%
    dplyr::left_join(water_supplies, ., by = "park_unit")
  
  if(verbose) {
    message("Visitation exposure calculated for ", nrow(indicators), " parks")
    
    timed_entry_baseline <- indicators %>% dplyr::filter(slope_adjustment == "baseline_timed_entry")
    if(nrow(timed_entry_baseline) > 0) {
      message("  Applied baseline growth (", baseline_growth_rate*100, "% annually) to ",
              nrow(timed_entry_baseline), " parks with timed entry:")
      message("    ", paste(timed_entry_baseline$park_unit, collapse = ", "))
    }
    
    negative_baseline <- indicators %>% dplyr::filter(slope_adjustment == "baseline_negative")
    if(nrow(negative_baseline) > 0) {
      message("  Applied baseline growth (", baseline_growth_rate*100, "% annually) to ",
              nrow(negative_baseline), " parks with negative trends:")
      message("    ", paste(negative_baseline$park_unit, collapse = ", "))
    }
    
    insufficient_baseline <- indicators %>% dplyr::filter(slope_adjustment == "baseline_insufficient_data")
    if(nrow(insufficient_baseline) > 0) {
      message("  Applied baseline growth (", baseline_growth_rate*100, "% annually) to ",
              nrow(insufficient_baseline), " parks with insufficient data (< 14 years):")
      message("    ", paste(insufficient_baseline$park_unit, collapse = ", "))
    }
    below_baseline <- indicators %>% dplyr::filter(slope_adjustment == "trend_below_baseline")
    if(nrow(below_baseline) > 0) {
      message("  Applied baseline growth (", baseline_growth_rate*100, "% annually) to ",
              nrow(below_baseline), " parks with less than the baseline growth rate:")
      message("    ", paste(below_baseline$park_unit, collapse = ", "))
    }
    
    
    unadjusted <- indicators %>% dplyr::filter(slope_adjustment == "none")
    message("  ", nrow(unadjusted), " parks using original observed trends")
  }
  
  return(norm_indicators)
}
