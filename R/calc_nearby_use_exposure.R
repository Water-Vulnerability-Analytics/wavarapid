#' Calculate Water Withdrawal Exposure Indicator
#'
#' Calculates changes in water withdrawals for counties containing water supplies
#' using historical data from USGS and future projections from Warziniack et al. 2022.
#' Computes percentile changes (p10, p50, p90) and the fraction of climate scenarios
#' showing increases in water use across multiple SSP-RCP combinations.
#'
#' @param water_supplies sf object with water supply locations
#' @param observed_use Warziniack et al. 2022 observed water use
#' @param future_use Warziniack et al. 2022 future water use
#' @param save logical, whether to save output
#' @param output_dir character, directory path for saving files
#' @param verbose logical, whether to show progress messages

calc_nearby_use_exposure <- function(water_supplies,
                                     observed_use, # = fread("data/all/RPAWaterDemand/1_BaseData/Trash 2022/usgswithdrawals.csv"),
                                     future_use, # = fread("data/all/RPAWaterDemand/2_Demand_Results/cons_and_withdrawal.csv"),
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
    
    counties <- st_transform(counties, st_crs(water_supplies))
    
    points_with_counties <- st_join(water_supplies, counties, join = st_nearest_feature) %>%
      mutate(fips = as.numeric(paste0(STATEFP, COUNTYFP))) %>%
      select(park_unit, park_name, wsd_source_id, fips)
    
    # aq = aquaculture; dp = domestic; ic = industrial + mining; ir = irrigated agriculture; ls = livestock; th = thermoelectric
    use_future <- future_use %>%
      mutate(fips = as.numeric(fips)) %>%
      filter(fips %in% unique(points_with_counties$fips),
             (yearr >= 2035 & yearr <= 2065)) %>%
      rename(year = yearr) %>%
      group_by(fips, year, scenario) %>%
      summarize(total_withdrawal = sum(withdrawal, na.rm = TRUE)) %>%
      # per county and scenario, grab mean annual withdrawal from 2035-2065
      group_by(fips, scenario) %>%
      summarize(future_withdrawal = mean(total_withdrawal, na.rm = TRUE)) %>%
      mutate(future_withdrawal = ifelse(future_withdrawal == "NaN", 0, future_withdrawal))
    
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
      summarize(total_withdrawal = sum(withdrawal, na.rm = TRUE)) %>%
      group_by(fips) %>%
      # get avg. total use across 1985-2015
      summarize(historic_withdrawal = mean(total_withdrawal, na.rm = TRUE)) %>%
      mutate(historic_withdrawal = ifelse(historic_withdrawal == "NaN", 0, historic_withdrawal))
    
    use_change <- use_future %>%
      inner_join(., use_historic, by = c("fips")) %>%
      select(fips, scenario, future_withdrawal, historic_withdrawal) %>%
      mutate(use_change = future_withdrawal - historic_withdrawal,
             percent_change = (use_change / historic_withdrawal) * 100)
    
    p10_values <- use_change %>%
      group_by(fips) %>%
      mutate(p10 = quantile(percent_change, 0.10, na.rm = TRUE),
             diff_from_p10 = abs(percent_change - p10)) %>%
      filter(diff_from_p10 == min(diff_from_p10, na.rm = TRUE)) %>%
      slice(1) %>%  # take first if ties
      ungroup() %>%
      select(fips,
             raw_use_change_p10 = percent_change)
    
    p50_values <- use_change %>%
      group_by(fips) %>%
      mutate(p50 = median(percent_change, na.rm = TRUE),
             diff_from_p50 = abs(percent_change - p50)) %>%
      filter(diff_from_p50 == min(diff_from_p50, na.rm = TRUE)) %>%
      slice(1) %>%  # take first if ties
      ungroup() %>%
      select(fips,
             raw_use_change_p50 = percent_change)
    
    p90_values <- use_change %>%
      group_by(fips) %>%
      mutate(p90 = quantile(percent_change, 0.90, na.rm = TRUE),
             diff_from_p90 = abs(percent_change - p90)) %>%
      filter(diff_from_p90 == min(diff_from_p90, na.rm = TRUE)) %>%
      slice(1) %>%  # take first if ties
      ungroup() %>%
      select(fips,
             raw_use_change_p90 = percent_change)
    
    # calculate fraction of models showing increase in use
    model_fractions <- use_change %>%
      group_by(fips) %>%
      summarize(raw_models_showing_increase_pct = 100 * sum(percent_change > 0, na.rm = TRUE) / n()) %>%
      ungroup()
    
    # Combine all indicators using standardized joins
    final_results <- p10_values %>%
      left_join(p50_values, by = c("fips")) %>%
      left_join(p90_values, by = c("fips")) %>%
      left_join(model_fractions, by = c("fips")) %>%
      # add standardized indicators using helper function
      add_standardized_indicators("raw_use_change", invert = TRUE) %>%
      add_standardized_indicators("raw_models_showing_increase") %>%
      right_join(., points_with_counties, by = "fips") %>%
      # reorder columns following standardized pattern
      select(wsd_source_id, park_unit, park_name, everything())
    
    # Something weird is going on with some of the FIPS codes... 
    # some RPA ones don't match up with the tigris FIPS. perhaps they've
    # changed?
    # missing <- final_results %>%
    #   filter(is.na(raw_use_change_p10)) %>%
    #   st_as_sf() %>%
    #   st_join(., counties)
    
  })})
  
  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(output_dir, paste0("nearby_use_exposure_", Sys.Date(), ".csv"))
    write_csv(final_results, filename)
    
    if (verbose) {
      message("Results saved to: ", filename)
    }
  }
  
  # Print timing
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message("\nWater withdrawal exposure processing completed in ",
            round(total_time, 2), " ", attr(total_time, "units"))
  }
  
  return(final_results)
  
}