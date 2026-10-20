#' Calculate Temperature Exposure Indicator (Standardized Version)
#'
#' Calculates temperature exposure using processed LOCA2/MWBM temperature data.
#' This follows the precipitation exposure helper structure, but uses annual mean
#' temperature and raw temperature change in degrees C rather than percent change.
#'
#' @param water_supplies sf object with water supply locations
#' @param temp_data character, directory containing processed LOCA2/MWBM temperature files
#' @param hindcast logical, if TRUE uses hindcast climate data, if FALSE uses
#'   training/measured downscaled climate data
#' @param historical Character vector of length 2 giving the start and end years
#'   (in `"YYYY"` format) that defines the "historical" period.
#' @param future Character vector of length 2 giving the start and end years
#'   (in `"YYYY"` format) that defines the "future" period.
#' @param water_year T/F, if TRUE, uses water-year bounds (Oct 1 - Sep 30) before annual aggregation.
#' @param save logical, whether to save output
#' @param output_dir character, directory path for saving files
#' @param verbose logical, whether to show progress messages
#'
#' @return tibble with water supply ID, park info, and temperature exposure indicators
#'
#' @export
calc_temp_exposure <- function(water_supplies,
                               temp_data = wava_data_path("rapid", "inputs", "loca2_temperature_water_supplies/"),
                               hindcast = TRUE,
                               historical,
                               future,
                               water_year = TRUE,
                               save = FALSE,
                               output_dir = NULL,
                               verbose = TRUE) {
  
  # Input validation
  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir)
  )
  
  # Start timing
  if (verbose) {
    message("Starting temperature exposure calculation...")
    start_time <- Sys.time()
  }
  
  # Extract year boundaries from date arguments
  hist_start_year <- as.integer(historical[1])
  hist_end_year <- as.integer(historical[2])
  future_start_year <- as.integer(future[1])
  future_end_year <- as.integer(future[2])
  
  all_files <- list.files(temp_data, full.names = TRUE, pattern = "\\.parquet$")
  
  if(hindcast == FALSE){
    # if we want to use training data... select associated files
    files <- all_files %>%
      grep("2015|2045|2075|livneh|lusu|training", ., value = TRUE, ignore.case = TRUE)
    
  } else if(hindcast == TRUE){
    # if we want to use hindcast climate data...
    files <- all_files %>%
      grep("training|livneh_unsplit", .,  invert = TRUE, value = TRUE, ignore.case = TRUE)
  }
  
  # If file names do not carry the same patterns as runoff/precip, assume the
  # provided folder already contains the intended temperature files.
  if (length(files) == 0) {
    files <- all_files
  }
  
  if (length(files) == 0) {
    stop("No processed temperature files found. Check your input data and file paths.")
  }
  
  files_dt <- data.table(file = files)
  
  # extract years
  files_dt[, `:=`(
    year_single = as.integer(str_extract(file, "(?<!\\d)(\\d{4})(?=\\.nc)")),
    year_start = as.integer(str_extract(file, "(\\d{4})(?=-\\d{4})")),
    year_end = as.integer(str_extract(file, "(?<=-)(\\d{4})(?=\\.)"))
  )]
  
  # fill single-year files into start/end only when a range was not found
  files_dt[!is.na(year_single) & is.na(year_start) & is.na(year_end), `:=`(
    year_start = year_single,
    year_end = year_single
  )]
  
  if(water_year == FALSE){
    # keep files overlapping either period
    files_subset <- files_dt[
      (year_end >= hist_start_year & year_start <= hist_end_year) |
        (year_end >= future_start_year  & year_start <= future_end_year) |
        (is.na(year_start) & is.na(year_end))
    ]
  } else {
    # keep files overlapping either period (with wyear accounted for)
    files_subset <- files_dt[
      (year_end >= hist_start_year - 1 & year_start <= hist_end_year) |
        (year_end >= future_start_year - 1  & year_start <= future_end_year) |
        (is.na(year_start) & is.na(year_end))
    ]
  }
  
  files <- files_subset$file
  
  # Check if we have any files to process
  if (length(files) == 0) {
    stop("No processed temperature files overlap the historical or future period.")
  }
  
  if (verbose) {
    message("Processing ", length(files), " temperature files")
    
    # Create progress bar for file processing
    pb <- progress::progress_bar$new(
      format = "Processing temperature files [:bar] :percent (:current/:total) ETA: :eta",
      total = length(files),
      clear = FALSE,
      width = 80)
  }
  
  # function that employs data.table to get annual averages from the raw tabular
  # data for every water supply, calculates the mean temperature across neighborhood
  # grids and times
  process_raw_grid_data <- function(f) {
    
    raw <- read_parquet(f) %>%
      dplyr::filter(wsd_source_id %in% water_supplies$wsd_source_id)
    
    annual_points <- as.data.table(raw)[
      ,
      year := if (water_year) {
        year(date) + fifelse(month(date) >= 10, 1L, 0L)
      } else {
        year(date)
      }
    ]
    
    # when using observed data, get rid of the modeled data
    # that overlaps with observed record
    if (!grepl("livneh|lusu", f, ignore.case = TRUE) && hindcast == FALSE) {
      annual_points <- annual_points[year >= 2019]
    }
    
    annual_points <- annual_points[
      # Use the extracted year boundaries from function arguments
      (year >= hist_start_year & year <= hist_end_year) |
        (year >= future_start_year & year <= future_end_year)
    ][
      # Average the available neighborhood cells independently at each time step.
      , .(site_temperature = if (any(is.finite(value))) {
        mean(value[is.finite(value)])
      } else {
        NA_real_
      }),
      by = .(wsd_source_id, date, year, park_unit, nc_file)
    ][
      # Then give each available time step equal weight in the annual mean.
      , .(mean_temperature = if (any(is.finite(site_temperature))) {
        mean(site_temperature[is.finite(site_temperature)])
      } else {
        NA_real_
      }),
      by = .(year, wsd_source_id, park_unit, nc_file)
    ]
    
    annual_points
  }
  
  # Process the extracted data files with progress bar
  mean_data <- map(files, ~{
    if (verbose) pb$tick(tokens = list(what = "Processing files"))
    process_raw_grid_data(.x)
  })
  
  combined_data <- mean_data %>% bind_rows()
  
  if (nrow(combined_data) == 0) {
    stop("No temperature records found after filtering.")
  }
  
  setDT(combined_data)
  
  suppressMessages({suppressWarnings({
    # fix common Jay Alder training filename patterns, if using
    combined_data[
      grepl("^Jay_Alder_.*_MWBM_LOCA2_training_[0-9]{6}-[0-9]{6}\\.nc$", nc_file),
      nc_file := sub(
        "^Jay_Alder_(.*)_MWBM_LOCA2_training_[0-9]{6}-[0-9]{6}\\.nc$",
        "Jay_Alder_\\1_MWBM_LOCA2.historical.training.1950-2100.nc",
        nc_file
      )
    ]
    
    # parse filename components
    nc_parts <- tstrsplit(combined_data$nc_file, "\\.", fill = NA)
    if (length(nc_parts) < 8) {
      nc_parts <- c(nc_parts, rep(list(NA_character_), 8 - length(nc_parts)))
    }
    
    combined_data[
      , c("var", "model", "ssp", "run", "nc", "misc1", "misc2", "misc3" ) :=
        nc_parts[1:8]
    ]
    
    # ssp populate with "historical" for hindcast = FALSE (as in, no ssp)
    combined_data[
      !ssp %in% c("ssp245", "ssp370", "ssp585", "historical"),
      ssp := dplyr::case_when(
        grepl("livneh|lusu", nc_file, ignore.case = TRUE) ~ "historical",
        TRUE ~ stringr::str_extract(nc_file, "ssp\\d+")
      )
    ]
    
    # Split run/time if present
    run_parts <- tstrsplit(combined_data$run, "_", fixed = TRUE, fill = NA)
    if (length(run_parts) < 2) {
      run_parts <- c(run_parts, list(NA_character_))
    }
    
    combined_data[
      , c("run", "time") := run_parts[1:2]
    ]
  })})
  
  # Drop unnecessary columns
  combined_data[, c("var", "nc") := NULL]
  
  if(hindcast == TRUE){
    
    # distinct (model, scenario) pairs that exist in the past and future
    ssps <- combined_data %>%
      filter(!ssp %in% c("1950-2014", "historical")) %>%
      distinct(model, ssp, run)
    
    temperature_change_sub <- vector("list", length = nrow(ssps))
    
    for(i in seq_len(nrow(ssps))){
      
      ssp_sub <- ssps[i,]
      
      hindcast_historical <- combined_data %>%
        filter(model == ssp_sub$model,
               ssp == "historical",
               run == ssp_sub$run) %>%
        filter(year >= hist_start_year & year <= hist_end_year)
      
      # if our historical record is actually considered "future" data...
      hindcast_future <- combined_data %>%
        filter(model == ssp_sub$model,
               ssp == ssp_sub$ssp,
               run == ssp_sub$run) %>%
        filter(year >= hist_start_year & year <= hist_end_year,
               !year %in% unique(hindcast_historical$year))
      
      historical_complete <- hindcast_historical %>%
        bind_rows(hindcast_future) %>%
        group_by(wsd_source_id, park_unit) %>%
        summarize(temperature_historical = mean(mean_temperature, na.rm = TRUE),
                  range_historical = paste0(min(year), "-", max(year)),
                  .groups = "drop")
      
      future_complete <- combined_data %>%
        filter(model == ssp_sub$model,
               ssp == ssp_sub$ssp,
               run == ssp_sub$run) %>%
        filter(year >= future_start_year & year <= future_end_year) %>%
        group_by(wsd_source_id, park_unit) %>%
        summarize(temperature_future = mean(mean_temperature, na.rm = TRUE),
                  period_future = paste0(min(year), "-", max(year)),
                  .groups = "drop")
      
      temperature_change_sub[[i]] <- future_complete %>%
        inner_join(., historical_complete, by = c("wsd_source_id", "park_unit")) %>%
        mutate(model = ssp_sub$model,
               ssp = ssp_sub$ssp,
               run = ssp_sub$run,
               temperature_change = temperature_future - temperature_historical)
    }
    
    temperature_change <- bind_rows(temperature_change_sub)
    
  } else if(hindcast == FALSE){
    
    # distinct (model, scenario) pairs that exist in the past and future
    ssps <- combined_data %>%
      filter(ssp != "training" & ssp != "historical") %>%
      distinct(model, ssp, run)
    
    training_historical <- combined_data %>%
      filter(ssp == "training" | ssp == "historical") %>%
      group_by(wsd_source_id, park_unit) %>%
      filter(year >= hist_start_year & year <= hist_end_year) %>%
      summarize(temperature_historical = mean(mean_temperature, na.rm = TRUE),
                range_historical = paste0(min(year), "-", max(year)),
                .groups = "drop")
    
    temperature_change_sub <- vector("list", length = nrow(ssps))
    
    for(i in seq_len(nrow(ssps))){
      suppressMessages({suppressWarnings({
        ssp_sub <- ssps[i,]
        
        future_complete <- combined_data %>%
          filter(model == ssp_sub$model,
                 ssp == ssp_sub$ssp,
                 run == ssp_sub$run) %>%
          filter(year >= future_start_year & year <= future_end_year) %>%
          group_by(wsd_source_id, park_unit) %>%
          summarize(temperature_future = mean(mean_temperature, na.rm = TRUE),
                    period_future = paste0(min(year), "-", max(year)),
                    .groups = "drop")
        
        temperature_change_sub[[i]] <- future_complete %>%
          inner_join(., training_historical, by = c("wsd_source_id", "park_unit")) %>%
          mutate(model = ssp_sub$model,
                 ssp = ssp_sub$ssp,
                 run = ssp_sub$run,
                 temperature_change = temperature_future - temperature_historical)
        
      })})
    }
    
    temperature_change <- bind_rows(temperature_change_sub)
    
  }
  
  if (nrow(temperature_change) == 0) {
    stop("No temperature model combinations found for the requested historical/future periods.")
  }
  
  # calculate percentile indicators:
  suppressMessages({suppressWarnings({
    
    p10_values <- temperature_change %>%
      group_by(park_unit, wsd_source_id) %>%
      mutate(p10 = quantile(temperature_change, 0.10, na.rm = TRUE),
             diff_from_p10 = abs(temperature_change - p10)) %>%
      filter(diff_from_p10 == min(diff_from_p10, na.rm = TRUE)) %>%
      slice(1) %>%  # take first if ties
      ungroup() %>%
      select(wsd_source_id,
             park_unit,
             raw_temperature_change_p10 = temperature_change)
    
    p50_values <- temperature_change %>%
      group_by(park_unit, wsd_source_id) %>%
      mutate(p50 = median(temperature_change, na.rm = TRUE),
             diff_from_p50 = abs(temperature_change - p50)) %>%
      filter(diff_from_p50 == min(diff_from_p50, na.rm = TRUE)) %>%
      slice(1) %>%  # take first if ties
      ungroup() %>%
      select(wsd_source_id,
             park_unit,
             raw_temperature_change_p50 = temperature_change)
    
    p90_values <- temperature_change %>%
      group_by(park_unit, wsd_source_id) %>%
      mutate(p90 = quantile(temperature_change, 0.90, na.rm = TRUE),
             diff_from_p90 = abs(temperature_change - p90)) %>%
      filter(diff_from_p90 == min(diff_from_p90, na.rm = TRUE)) %>%
      slice(1) %>%  # take first if ties
      ungroup() %>%
      select(wsd_source_id,
             park_unit,
             raw_temperature_change_p90 = temperature_change)
    
    # calculate fraction of models showing warming
    model_fractions <- temperature_change %>%
      group_by(park_unit, wsd_source_id) %>%
      summarize(raw_temp_models_showing_increase = 100 * sum(temperature_change > 0, na.rm = TRUE) / n(),
                .groups = "drop")
    
  })})
  
  # Combine all indicators using standardized joins
  final_results <- p10_values %>%
    left_join(p50_values, by = c("wsd_source_id", "park_unit")) %>%
    left_join(p90_values, by = c("wsd_source_id", "park_unit")) %>%
    left_join(model_fractions, by = c("wsd_source_id", "park_unit")) %>%
    # add standardized indicators using helper function
    add_standardized_indicators("raw_temperature_change", invert = FALSE) %>%
    add_standardized_indicators("raw_temp_models_showing_increase") %>%
    # reorder columns following standardized pattern
    select(wsd_source_id, park_unit, everything())
  
  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(output_dir, paste0("temperature_exposure_", Sys.Date(), ".csv"))
    write_csv(final_results, filename)
    
    if (verbose) {
      message("Results saved to: ", filename)
    }
  }
  
  # Print timing
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message("\ntemperature exposure processing completed in ",
            round(total_time, 2), " ", attr(total_time, "units"))
  }
  
  return(final_results)
}
