# HELPER FUNCTIONS =========================================================

#' Add standardized normalized and ranked columns to indicator data
#' 
#' @param data tibble with raw indicator columns
#' @param prefix character, prefix pattern to match columns (e.g., "raw_fire_prob")
#' @param invert logical, whether to invert the direction of scaling so that
#'   higher raw values correspond to lower normalized and ranked scores.
#' @return tibble with added norm_ and rank_ columns
#' 
add_standardized_indicators <- function(data, prefix, invert = FALSE) {
  # Get columns that match the prefix
  cols_to_transform <- names(data)[str_starts(names(data), prefix)]
  
  # Create normalized columns
  for (col in cols_to_transform) {
    norm_col <- str_replace(col, "raw_", "norm_")
    rank_col <- str_replace(col, "raw_", "rank_")
    
    data[[norm_col]] <- scales::rescale(data[[col]], to = c(0, 1))
    data[[rank_col]] <- cume_dist(data[[col]]) * 100
    
    if (invert) {
      data[[norm_col]] <- 1 - data[[norm_col]]
      data[[rank_col]] <- 100 - data[[rank_col]]
    }
    
  }
  
  return(data)
}

# RUNOFF PROCESSING HELPER FUNCTIONS =====================================

#' Process raw grid data to calculate mean runoff statistics
#' 
#' @param f file path to processed runoff CSV file
#' @param water_supplies sf object with water supply locations
#' 
process_raw_grid_data <- function(f, water_supplies) {
  
  raw <- fread(f)
  
  #filter to water supplies of interest
  raw <- raw %>% filter(wsd_source_id %in% water_supplies$wsd_source_id)
  
  # 1) annual totals per point/model
  annual_point <- raw[
    , .(grid_runoff = sum(value, na.rm = TRUE)),
    by = .(year = as.integer(format(date, "%Y")), wsd_source_id, park_unit, nc_file, cell_lon, cell_lat)
  ]
  
  # 2) average annual totals per site/model
  annual <- annual_point[
    , .(total_runoff = mean(grid_runoff, na.rm = TRUE)),
    by = .(year, wsd_source_id, park_unit, nc_file)
  ]
  
  # 3) mean of annual totals per site/model
  mean_annual <- annual[
    , .(mean_runoff = mean(total_runoff, na.rm = TRUE)),
    by = .(wsd_source_id, park_unit, nc_file)
  ]
  
  # 4) carry time_units (first), min_date, max_date from raw
  meta <- raw[
    , .(
      time_units = first(time_units),
      min_date   = min(date, na.rm = TRUE),
      max_date   = max(date, na.rm = TRUE)
    ),
    by = .(wsd_source_id, park_unit, nc_file)
  ]
  
  # 4) combine, then fix filename and parse components
  out <- mean_annual[meta, on = .(wsd_source_id, park_unit, nc_file)]
  
  # fix specific filename pattern
  out[
    nc_file == "Jay_Alder_runoff_MWBM_LOCA2_training_195001-201412.nc",
    nc_file := "Jay_Alder_runoff_MWBM_LOCA2.historical.training.1950-2100.nc"
  ]
  
  # parse filename components
  out[
    , c("var", "model", "period", "run", "nc") := tstrsplit(nc_file, "\\.")
  ]
  
  # Split run/time if present
  out[
    !is.na(run), c("run", "time") := tstrsplit(run, "_", fixed = TRUE)
  ]
  
  # Drop unnecessary columns
  out[, c("var", "nc") := NULL]
  
  # Reorder to match original-style output
  setcolorder(out, c(
    "wsd_source_id", "park_unit", "nc_file",
    "mean_runoff", "time_units", "min_date", "max_date",
    "model", "period", "run", "time"
  ))
  
  return(out[])
}