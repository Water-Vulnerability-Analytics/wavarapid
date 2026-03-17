#' Calculate Nearby Water Use Sensitivity Indicator
#'
#' Calculates a county-level nearby water use sensitivity indicator for water supplies
#' based on historical total water withdrawals from USGS data. Water supplies are
#' assigned to the nearest county, historical withdrawals are aggregated across
#' major use sectors, and county-level withdrawal trends are estimated using
#' Sen's slope and Mann-Kendall trend tests. Outputs include raw average
#' withdrawals, raw trend estimates, trend direction, and standardized indicators.
#'
#' @param water_supplies sf object with water supply locations.
#' @param observed_use data frame of historical county-level water withdrawals.
#'   Defaults to `fread("data/all/RPAWaterDemand/1_BaseData/Trash 2022/usgswithdrawals.csv")`.
#'   Must contain columns for year, county/state identifiers, FIPS, and sectoral
#'   withdrawal variables used in the function.
#' @param save logical, whether to save output. Default is `FALSE`.
#' @param output_dir character, directory path for saving files. Required if
#'   `save = TRUE`.
#' @param verbose logical, whether to print progress and timing messages. Default is `TRUE`.

calc_nearby_use_sensitivity <- function(water_supplies,
                                        observed_use = fread("data/all/RPAWaterDemand/1_BaseData/Trash 2022/usgswithdrawals.csv"), #precompiled data from USFS
                                        save = FALSE,
                                        output_dir = NULL,
                                        verbose = TRUE) {
  
  # Input validation
  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir)
  )
  
  start_time <- Sys.time()
  
  suppressMessages({suppressWarnings({
    
    counties <- counties(cb = TRUE, year = 2019)
    
    counties <- st_transform(counties, st_crs(water_supplies)) %>%
      mutate(fips = as.numeric(paste0(STATEFP, COUNTYFP)),
             county_area_km2 = as.numeric(st_area(st_transform(geometry, 5070))) / 1000000)
    
    points_with_counties <- st_join(water_supplies, counties, join = st_nearest_feature) %>%
      mutate(fips = as.numeric(paste0(STATEFP, COUNTYFP))) %>%
      select(park_unit, wsd_source_id, fips, county_area_km2)
    
    # aq = aquaculture; dp = domestic; ic = industrial + mining; ir = irrigated agriculture; ls = livestock; th = thermoelectric
    #1985-2015 observed data:
    use_historic <- observed_use %>%
      select(year, state_name, county_nm, fips = FIPS, PStot,
             aq = AquaSelfTot,
             dp = DomTot,
             ic = IndTot,
             ir = IrrigSelfTot,
             ls = StockSelfTot,
             th = ThermSelfTot,
             min = MinSelfTot) %>%
      mutate(across(c(aq, dp, ic, ir, ls, th, min), ~as.numeric(na_if(as.character(.x), "-"))),
             ic = ic + min,
             fips = as.numeric(fips)) %>%
      select(-min) %>%
      # subset to only counties with water supplies
      filter(fips %in% unique(points_with_counties$fips)) %>%
      pivot_longer(c(aq, dp, ic, ir, ls, th), names_to = "sector", values_to = "withdrawal") %>%
      group_by(state_name, county_nm, fips, year) %>%
      summarize(total_withdrawal = sum(withdrawal, na.rm = TRUE), .groups = "drop") %>%
      left_join(points_with_counties %>% st_drop_geometry() %>% distinct(fips, county_area_km2),
                by = "fips") %>%
      mutate(total_withdrawal_per_km2 = total_withdrawal / county_area_km2)
    
    county_trends <- use_historic %>%
      group_by(state_name, county_nm, fips) %>%
      arrange(year, .by_group = TRUE) %>%
      summarize(raw_avg_nearby_use = mean(total_withdrawal, na.rm = TRUE),
                raw_avg_nearby_use_per_km2 = mean(total_withdrawal_per_km2, na.rm = TRUE),
                data = list(dplyr::pick(dplyr::everything())),
                .groups = "drop") %>%
      mutate(n_years = map_int(data, ~ sum(!is.na(.x$total_withdrawal))),
             sen = map(data, ~ {
               x <- .x$total_withdrawal
               if (sum(!is.na(x)) >= 3) sens.slope(x) else NULL
             }),
             sen_per_km2 = map(data, ~ {
               x <- .x$total_withdrawal_per_km2
               if (sum(!is.na(x)) >= 3) sens.slope(x) else NULL
             }),
             mk = map(data, ~ {
               x <- .x$total_withdrawal
               if (sum(!is.na(x)) >= 3) mk.test(x) else NULL
             }),
             raw_slope_nearby_use = map_dbl(sen, ~ if (is.null(.x)) NA_real_ else .x$estimates[[1]]),
             raw_slope_nearby_use_per_km2 = map_dbl(sen_per_km2, ~ if (is.null(.x)) NA_real_ else .x$estimates[[1]]),
             mk_tau = map_dbl(mk, ~ if (is.null(.x)) NA_real_ else .x$estimates[[1]]),
             mk_p_value = map_dbl(mk, ~ if (is.null(.x)) NA_real_ else .x$p.value),
             mk_significant = mk_p_value < 0.05) %>%
      select(-data, -sen, -sen_per_km2, -mk, -mk_p_value, -mk_tau) %>%
      mutate(trend_direction = case_when(is.na(raw_slope_nearby_use) ~ NA_character_,
                                         raw_slope_nearby_use > 0 ~ "increasing",
                                         raw_slope_nearby_use < 0 ~ "decreasing",
                                         TRUE ~ "no change"),
             trend_direction_per_km2 = case_when(is.na(raw_slope_nearby_use_per_km2) ~ NA_character_,
                                                 raw_slope_nearby_use_per_km2 > 0 ~ "increasing",
                                                 raw_slope_nearby_use_per_km2 < 0 ~ "decreasing",
                                                 TRUE ~ "no change"),
             raw_percent_slope_nearby_use = ifelse(is.na(raw_slope_nearby_use) | is.na(raw_avg_nearby_use) | raw_avg_nearby_use == 0,
                                                   NA_real_,
                                                   (raw_slope_nearby_use / raw_avg_nearby_use) * 100),
             raw_percent_slope_nearby_use_per_km2 = ifelse(is.na(raw_slope_nearby_use_per_km2) | is.na(raw_avg_nearby_use_per_km2) | raw_avg_nearby_use_per_km2 == 0,
                                                           NA_real_,
                                                           (raw_slope_nearby_use_per_km2 / raw_avg_nearby_use_per_km2) * 100))
    
    
    # Combine all indicators using standardized joins
    final_results <- county_trends %>%
      # add standardized indicators using helper function
      add_standardized_indicators("raw_slope_nearby_use") %>%
      add_standardized_indicators("raw_avg_nearby_use") %>%
      add_standardized_indicators("raw_percent_slope_nearby_use") %>%
      add_standardized_indicators("raw_slope_nearby_use_per_km2") %>%
      add_standardized_indicators("raw_avg_nearby_use_per_km2") %>%
      add_standardized_indicators("raw_percent_slope_nearby_use_per_km2") %>%
      right_join(., points_with_counties %>% st_drop_geometry(), by = "fips") %>%
      # reorder columns following standardized pattern
      select(wsd_source_id, park_unit, everything())
    
  })})
  
  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(output_dir, paste0("nearby_use_sensitivity_", Sys.Date(), ".csv"))
    write_csv(final_results, filename)
    
    if (verbose) {
      message("Results saved to: ", filename)
    }
  }
  
  # Print timing
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message("\nNearby water withdrawal sensitivity processing completed in ",
            round(total_time, 2), " ", attr(total_time, "units"))
  }
  
  return(final_results)
  
}