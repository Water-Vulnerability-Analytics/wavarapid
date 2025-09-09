# HELPER FUNCTIONS =========================================================

#' Add standardized normalized and ranked columns to indicator data
#' 
#' @param data tibble with raw indicator columns
#' @param prefix character, prefix pattern to match columns (e.g., "raw_fire_prob")
#' 
#' @return tibble with added norm_ and rank_ columns
#' 
add_standardized_indicators <- function(data, prefix) {
  # Get columns that match the prefix
  cols_to_transform <- names(data)[str_starts(names(data), prefix)]
  
  # Create normalized columns
  for (col in cols_to_transform) {
    norm_col <- str_replace(col, "raw_", "norm_")
    rank_col <- str_replace(col, "raw_", "rank_")
    
    data[[norm_col]] <- scales::rescale(data[[col]], to = c(0, 1))
    data[[rank_col]] <- cume_dist(data[[col]]) * 100
  }
  
  return(data)
}




# RUNOFF PROCESSING HELPER FUNCTIONS =====================================

#' Process NetCDF files for runoff data extraction
#' 
#' @param nc_path path to NetCDF file
#' @param points sf object with point locations  
#' @param var character, variable name to extract (default "runoff")
#' 
loca2_mwbm_processor <- function(nc_path, output_dir, points, var = "runoff") {
  
  message("Starting vectorized extraction...")
  
  start_time <- Sys.time()
  
  # open using ncdf4 package: 
  nc <- nc_open(nc_path)
  
  # grab layer of user-selected variable
  v <- nc$var[[var]]
  nc  <- nc_open(nc_path)
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  
  # If longitudes are 0–360, convert to -180–180
  # Alder training historical data is in 0-360 whereas others are not...
  if (any(lon > 180)) lon <- ifelse(lon > 180, lon - 360, lon)
  
  # transform points to nc coordinate system if not already
  if (st_crs(points) != st_crs(rast(nc_path))) {
    points <- st_transform(points, st_crs(rast(nc_path)))
  }
  
  grid_sf <- expand.grid(lon = lon, lat = lat) %>%
    st_as_sf(coords = c("lon", "lat"), crs = st_crs(rast(nc_path)))
  
  # get time units 
  time_seq <-  ncvar_get(nc, "time")
  time_units <- nc$dim[["time"]]$units
  
  # Filter time indices based on units
  time_indices_to_keep <- seq_along(time_seq)  # Default: keep all
  
  if (!is.na(time_units) && grepl("days since 1950-01-01", time_units, ignore.case = TRUE)) {
    message("Found 'days since 1950-01-01' units - filtering to 1980-2010 and 2035-2065 periods")
    
    # Convert days since 1950-01-01 to actual years
    start_1981 <- 31 * 365.25    
    end_2010 <- 61 * 365.25      
    start_2035 <- 85 * 365.25    
    end_2065 <- 116 * 365.25     
    
    # grab grids within our desired periods
    period_1981_2010 <- which(time_seq >= start_1981 & time_seq <= end_2010)
    period_2035_2065 <- which(time_seq >= start_2035 & time_seq <= end_2065)
    
    time_indices_to_keep <- c(period_1981_2010, period_2035_2065)
    time_indices_to_keep <- sort(time_indices_to_keep)  # Keep in chronological order
    
    message(sprintf("Filtered from %d to %d time steps", length(time_seq), length(time_indices_to_keep)))
    message(sprintf("1980-2010 period: %d steps, 2035-2065 period: %d steps", 
                    length(period_1981_2010), length(period_2035_2065)))
  }
  
  # update time_seq and nt to reflect filtered data
  time_seq <- time_seq[time_indices_to_keep]
  nt <- length(time_seq)
  
  # convert sf points to terra vect object and extract coordinates
  v_pts <- vect(points)
  xy <- crds(v_pts)
  
  # convert geographic coordinates to array indices by finding nearest nc grid cell center
  lon_index <- map_int(xy[, 1], ~ which.min(abs(lon - .x)))
  # closest latitude index for each point
  lat_index <- map_int(xy[, 2], ~ which.min(abs(lat - .x)))
  # combine longitude and latitude indices into matrix of closest nc grid cell:
  center_indices  <- cbind(lon_index = lon_index, lat_index = lat_index)
  
  # get surrounding (3x3) neighborhoods around each point 
  # (i.e., -1, 0, +1 in both lon and lat directions)
  offsets <- expand.grid(lon_offset = -1:1, lat_offset = -1:1)
  
  # loop through each point to calculate 3x3 neighborhood grid
  neighbor_list <- vector("list", nrow(center_indices))
  
  for (i in 1:nrow(center_indices)) {
    # Get center longitude and latitude indices for current point
    lon_index <- center_indices[i, "lon_index"]
    lat_index <- center_indices[i, "lat_index"]
    
    # apply "offsets" 
    nbr_lons <- pmax(1, pmin(length(lon), lon_index + offsets$lon_offset))
    nbr_lats <- pmax(1, pmin(length(lat), lat_index + offsets$lat_offset))
    
    # store
    neighbor_list[[i]] <- data.frame(
      point_id = i,
      lon_index = nbr_lons,
      lat_index = nbr_lats
    )
  }
  
  neighbor_list <- bind_rows(neighbor_list) %>%
    # Add actual geographic coordinates for each neighborhood cell
    mutate(cell_lon = lon[lon_index],
           cell_lat = lat[lat_index])
  
  # extract dimension names and sizes from the variable
  dnames <- map_chr(v$dim, "name")
  dsizes <- map_int(v$dim, "len")
  
  # identify time dimension by looking for "time" in dimension names (case insensitive)
  time_dim <- which(grepl("time", dnames, ignore.case = TRUE))
  
  # make matrix to store info in
  n_extractions <- nrow(neighbor_list)
  result_matrix <- matrix(NA_real_, nrow = n_extractions, ncol = nt)
  
  message(sprintf("Extracting %d locations across %d time steps...", n_extractions, nt))
  
  # If we have filtered time indices, we need a different approach
  if (length(time_indices_to_keep) < length(ncvar_get(nc, "time"))) {
    # Read all filtered time steps at once by reading the full array and subsetting
    message("Reading filtered time data...")
    
    # Read the entire variable
    full_arr <- ncvar_get(nc, var)
    
    # Subset to only the time indices we want
    if (length(time_dim) == 1) {
      if (time_dim == 1) {
        arr <- full_arr[time_indices_to_keep, , ]
        # Reorder to (lon, lat, time)
        arr <- aperm(arr, c(2, 3, 1))
      } else if (time_dim == 2) {
        arr <- full_arr[, time_indices_to_keep, ]
        # Reorder to (lon, lat, time)
        arr <- aperm(arr, c(1, 3, 2))
      } else {
        # time_dim == 3, already (lon, lat, time)
        arr <- full_arr[, , time_indices_to_keep]
      }
    }
    
    # extract values for all locations and times
    for (i in seq_len(n_extractions)) {
      lon_index <- neighbor_list$lon_index[i]
      lat_index <- neighbor_list$lat_index[i]
      result_matrix[i, ] <- arr[lon_index, lat_index, ]
    }
    
    message("Completed filtered extraction")
    
  } else {
    # use original chunking approach for non-filtered data
    chunk_size <- min(50, nt)  
    n_chunks <- ceiling(nt / chunk_size)
    
    for (chunk in 1:n_chunks) {
      start_t <- (chunk - 1) * chunk_size + 1
      end_t <- min(chunk * chunk_size, nt)
      chunk_nt <- end_t - start_t + 1
      
      # start and count vectors for NetCDF data reading (start at index 1 for all dimensions)
      start <- rep(1L, length(dsizes))
      count <- dsizes
      
      # modify start and count for time dimension if it exists
      if (length(time_dim) == 1) {
        start[time_dim] <- start_t
        count[time_dim] <- chunk_nt
      }
      
      # read the data chunk from nc file
      arr <- ncvar_get(nc, var, start = start, count = count, collapse_degen = FALSE)
      
      # reorder dimensions to ensure (lon, lat, time) order regardless of original dimension order
      if (length(time_dim) == 1) {
        if (time_dim == 1) {
          arr <- aperm(arr, c(2, 3, 1))
        } else if (time_dim == 2) {
          arr <- aperm(arr, c(1, 3, 2))
        }
        # If time_dim == 3, already in (lon, lat, time) order
      }
      
      # get values for all spatial locations and times in this chunk
      for (i in seq_len(n_extractions)) {
        # get lat-long indices for current extraction point
        lon_index <- neighbor_list$lon_index[i]
        lat_index <- neighbor_list$lat_index[i]
        
        # handle single time step vs multiple time steps differently
        if (chunk_nt == 1) {
          # single time step: extract scalar value
          result_matrix[i, start_t] <- arr[lon_index, lat_index]
        } else {
          # multiple time steps: extract time series
          result_matrix[i, start_t:end_t] <- arr[lon_index, lat_index, ]
        }
      }
      
      # progress bar
      if (chunk %% max(1, n_chunks %/% 10) == 0 || chunk == n_chunks) {
        message(sprintf("Processed chunk %d / %d (time steps %d-%d)", chunk, n_chunks, start_t, end_t))
      }
    }
  }
  
  # convert matrix results to tidy data frame
  results <- expand.grid(extraction_id = seq_len(n_extractions),
                         time_index = seq_len(nt),
                         stringsAsFactors = FALSE) %>%
    # add extracted values as vector and corresponding time values
    mutate(value = as.vector(result_matrix),
           time = if (!is.null(time_seq)) time_seq[time_index] else NA) %>%
    # join with neighborhood indices to add spatial information
    left_join(neighbor_list %>% mutate(extraction_id = row_number()), by = "extraction_id") %>%
    select(-extraction_id)
  
  end_time <- Sys.time()
  
  message(sprintf("Extraction completed in %.2f seconds", as.numeric(difftime(end_time, start_time, units = "secs"))))
  
  # extract point attributes (non-spatial data) and add point ID column
  pts_attr <- st_drop_geometry(points) %>% mutate(.pt_id = row_number())
  
  # convert results to data.table and join with point attributes
  results_dt <- as.data.table(results)[
    as.data.table(pts_attr)[, .pt_id := .I], 
    on = c("point_id" = ".pt_id")
  ][, `:=`(
    layer_index = time_index,
    time = as.character(time),
    time_units = time_units,           
    nc_file = basename(nc_path),
    cell_lon = as.character(cell_lon),
    cell_lat = as.character(cell_lat),
    date = as.Date("1950-01-01") + as.numeric(time),  # Convert days since 1950-01-01 to date
    time_index = NULL 
  )]
  
  # Create output directory if it doesn't exist
 #output_dir <- "data/rapid/loca2_mwbm_runoff/"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  write_csv(results_dt, paste0(output_dir, basename(nc_path), ".csv"))
  
  nc_close(nc)
  gc()
}

#' Process raw grid data to calculate mean runoff statistics
#' 
#' @param f file path to processed runoff CSV file
#' 
process_raw_grid_data <- function(f) {
  
  raw <- fread(f)
  
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