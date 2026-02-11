#' Calculate visitation-based exposure indicators for park water supplies
#'
#' This function computes visitation-related exposure indicators for a set of park units
#' based on annual and seasonal visitation trends. It uses NPS visitation data to estimate
#' historical visitation patterns, trends, and projected changes through 2050.
#'
#' For parks with timed entry systems, negative trends, or insufficient data (< 14 years of
#' monthly records available), a minimum baseline growth rate is applied based on long-term
#' NPS system-wide trends. This reflects that these parks are unlikely to experience
#' unconstrained historical growth patterns or do not have enough information to support
#' stable trend estimation.
#'
#' Results are merged with the provided `water_supplies` dataset for further analysis.
#'
#' @param water_supplies A data frame containing at least a `park_unit` column, representing
#'   unique park identifiers for which visitation exposure indicators will be calculated.
#' @param baseline_growth_rate Default growth rate to set to if actual growth
#'  rate per park is either below it, negative, not enough data, or for parks with timed entry.
#' @param baseline_seasonality_shift Numeric. A dampening factor (0-1) applied to the
#'   coefficient of variation (CV) trend for parks using baseline growth adjustments.
#'   For example, if 0.25, this means that the seasonality trend is reduced to 25% of the observed rate.
#'   Lower values assume parks with constrained growth (timed entry, negative trends, or
#'   insufficient data) will experience minimal changes in visitation seasonality patterns.
#'   This is still a WIP!
#' @param save Logical. If `TRUE`, results are saved to a file in `output_dir`.
#'   Defaults to `FALSE`.
#' @param output_dir Optional character string specifying the directory where results
#'   will be saved if `save = TRUE`. Ignored otherwise.
#' @param verbose Logical. If `TRUE`, progress messages are printed to the console.
#'   Defaults to `TRUE`.
#'
#' @return
#' A tibble combining the original `water_supplies` data with visitation exposure indicators.
calc_visitation_exposure <- function(water_supplies,
                                     baseline_growth_rate = 0.0087,
                                     baseline_seasonality_shift = 1,
                                     save = FALSE,
                                     output_dir = NULL,
                                     verbose = TRUE) {
  
  # timed entry start years
  timed_entry_parks <- tribble(~park_unit, ~timed_entry_year,
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
    
    visitation <- getUnitVisitation(units = park_unit, startYear = 2010, endYear = 2024)
    
    # check if park has timed entry
    timed_entry_info <- timed_entry_parks %>% filter(.data$park_unit == !!park_unit)
    has_timed_entry <- nrow(timed_entry_info) > 0
    timed_entry_year <- if(has_timed_entry) timed_entry_info$timed_entry_year else NA
    
    # calculate annual totals
    annual_visitation <- visitation %>%
      group_by(Year) %>%
      summarize(total = sum(RecreationVisitors), .groups = "drop")
    
    seasonal_visitation <- visitation %>%
      mutate(date = ym(paste0(Year, "-", Month)),
             season = case_when(Month %in% c(12, 1, 2) ~ "Winter",
                                Month %in% c(3, 4, 5) ~ "Spring",
                                Month %in% c(6, 7, 8) ~ "Summer",
                                Month %in% c(9, 10, 11) ~ "Fall"),
             month_name = month.abb[Month],
             season_month = factor(paste0(season, "_", month_name),
                                   levels = c("Winter_Jan","Winter_Feb","Spring_Mar","Spring_Apr","Spring_May",
                                              "Summer_Jun","Summer_Jul","Summer_Aug",
                                              "Fall_Sep","Fall_Oct","Fall_Nov",
                                              "Winter_Dec")))
    
    # fit model on full data (2010-2024)
    lm_total <- lm(total ~ Year, data = annual_visitation)
    original_slope <- coef(lm_total)[2]
    intercept <- coef(lm_total)[1]
    
    # calculate baseline growth in absolute visitors per year
    # baseline_slope <- pred_2024 * baseline_growth_rate
    # 
    # pred_2024 <- annual_visitation %>% filter(Year == 2024) %>% pull(total)
    
    pred_2024 <- predict(lm_total, newdata = data.frame(Year = 2024))
    baseline_slope <- pred_2024 * baseline_growth_rate
    
    # check if baseline growth adjustment is needed
    negative_trend <- original_slope < 0
    insufficient_data <- nrow(visitation) < (14 * 12)
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
      filter(total == max(total, na.rm = TRUE)) %>%
      pull(Year) %>%
      paste(collapse = ", ")
    
    # projected change in seasonality (exclude COVID years):
    monthly_percentages <- visitation %>%
      filter(!Year %in% c(2020, 2021)) %>%
      group_by(Year) %>%
      mutate(annual_total = sum(RecreationVisitors),
             monthly_percentage = RecreationVisitors / annual_total * 100) %>%
      ungroup()
    
    seasonal_variation <- monthly_percentages %>%
      group_by(Year) %>%
      summarize(cv_monthly = sd(monthly_percentage) / mean(monthly_percentage) * 100,
                .groups = "drop")
    
    lm_cv <- lm(cv_monthly ~ Year, data = seasonal_variation)
    original_cv_slope <- coef(lm_cv)[2]
    
    cv_2024 <- predict(lm_cv, newdata = data.frame(Year = 2024))
    
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
    tibble(
      park_unit = park_unit,
      pred_2050 = pred_2050,
      obs_2024 = pred_2024,
      peak_year = peak_year,
      visit_raw_pct_change = pct_change,
      visit_raw_cv_change = cv_change,
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
    map(possibly(~ visitation_indicators(.x), otherwise = NULL)) %>%
    keep(~!is.null(.)) %>%
    bind_rows()
  
  norm_indicators <- indicators %>%
    add_standardized_indicators(., prefix = "visit") %>%
    left_join(water_supplies, ., by = "park_unit")
  
  if(verbose) {
    message("Visitation exposure calculated for ", nrow(indicators), " parks")
    
    timed_entry_baseline <- indicators %>% filter(slope_adjustment == "baseline_timed_entry")
    if(nrow(timed_entry_baseline) > 0) {
      message("  Applied baseline growth (", baseline_growth_rate*100, "% annually) to ",
              nrow(timed_entry_baseline), " parks with timed entry:")
      message("    ", paste(timed_entry_baseline$park_unit, collapse = ", "))
    }
    
    negative_baseline <- indicators %>% filter(slope_adjustment == "baseline_negative")
    if(nrow(negative_baseline) > 0) {
      message("  Applied baseline growth (", baseline_growth_rate*100, "% annually) to ",
              nrow(negative_baseline), " parks with negative trends:")
      message("    ", paste(negative_baseline$park_unit, collapse = ", "))
    }
    
    insufficient_baseline <- indicators %>% filter(slope_adjustment == "baseline_insufficient_data")
    if(nrow(insufficient_baseline) > 0) {
      message("  Applied baseline growth (", baseline_growth_rate*100, "% annually) to ",
              nrow(insufficient_baseline), " parks with insufficient data (< 14 years):")
      message("    ", paste(insufficient_baseline$park_unit, collapse = ", "))
    }
    below_baseline <- indicators %>% filter(slope_adjustment == "trend_below_baseline")
    if(nrow(below_baseline) > 0) {
      message("  Applied baseline growth (", baseline_growth_rate*100, "% annually) to ",
              nrow(below_baseline), " parks with less than the baseline growth rate:")
      message("    ", paste(below_baseline$park_unit, collapse = ", "))
    }
    
    
    unadjusted <- indicators %>% filter(slope_adjustment == "none")
    message("  ", nrow(unadjusted), " parks using original observed trends")
  }
  
  return(norm_indicators)
}
