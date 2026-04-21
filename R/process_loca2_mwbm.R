#' Process NetCDF files for LOCA2 water balance data extraction
#' 
#' @param nc_path path to NetCDF files - found in "data/all/loca2_mwbm/"
#' @param water_supplies sf object with point locations  
#' @param var character, variable name to extract (default "runoff")
#' @param historical Character vector of length 2 giving the start and end dates
#'   (in `"YYYY-MM-DD"` format) that defines the "historical" period.  
#'   Defaults to `c("1980-01-01", "2024-12-31")`.
#' @param future Character vector of length 2 giving the start and end dates
#'   (in `"YYYY-MM-DD"` format) that defines the "future" period.  
#'   Defaults to `c("2035-01-01", "2065-12-31")`.
#' @param output_dir Directory to save the extracted CSV output files.
#' 
process_loca2_mwbm <- function(nc_path,
                               water_supplies,
                               var, # can be "runoff", "pr", "swe", "pet"
                               historical = c("1980-01-01", "2024-12-31"),
                               future = c("2035-01-01", "2065-12-31"),
                               output_dir){ 
  
  
  message(paste0("Starting vectorized extraction of ", nc_path))
  
  start_time <- Sys.time()
  
  # open using ncdf4 package: 
  nc <- nc_open(nc_path)
  
  # some renaming required
  if(var == "pr" & !grepl("alder", nc_path, ignore.case = TRUE)){
    var <- "pr_tavg"
  } 
  
  if(var == "pr_tavg" & grepl("livneh_unsplit_precip", nc_path, ignore.case = TRUE)){
    var <- "PRCP"
  } 
  
  if(var == "tasmin_tavg" & grepl("lusu", nc_path, ignore.case = TRUE)){
    var <- "Tmin"
  } 
  
  if(var == "tasmax_tavg" & grepl("lusu", nc_path, ignore.case = TRUE)){
    var <- "Tmax"
  } 
  
  if(var == "swe"){
    var <- "snow"
  } 
  
  # grab layer of user-selected variable
  v <- nc$var[[var]]
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  
  # If longitudes are 0–360, convert to -180–180
  # Alder training historical data is in 0-360 whereas others are not...
  if (any(lon > 180)) lon <- ifelse(lon > 180, lon - 360, lon)
  
  # transform water_supplies to nc coordinate system if not already
  if (st_crs(water_supplies) != st_crs(rast(nc_path))) {
    water_supplies <- st_transform(water_supplies, st_crs(rast(nc_path)))
  }
  
  # grid_sf <- expand.grid(lon = lon, lat = lat) %>%
  #   st_as_sf(coords = c("lon", "lat"), crs = st_crs(rast(nc_path)))
  
  # get time units 
  time_seq <-  if(var == "PRCP"){ncvar_get(nc, "Time")} else {ncvar_get(nc, "time") }
  time_units <- if(var == "PRCP"){nc$dim[["Time"]]$units} else {nc$dim[["time"]]$units}
  start_time_date <- lubridate::ymd(sub(" .*", "", sub(".*since\\s+", "", time_units)))
  
  # get mwbm units
  mwbm_units <- ncatt_get(nc, var, "units")$value
  
  # Filter time indices based on units
  time_indices_to_keep <- seq_along(time_seq)  
  
  if (!is.na(time_units) && grepl("days since", time_units, ignore.case = TRUE)) {
    # message(sprintf("Filtering to %s to %s and %s to %s periods", 
    #                 historical[1], historical[2], future[1], future[2]))
    
    # convert date boundaries to days since start_time_date
    # (how these cdf files record times)
    origin_date <- as.Date(start_time_date)
    
    # Use the historical and future arguments
    start_historical <- as.numeric(as.Date(historical[1]) - origin_date)
    end_historical   <- as.numeric(as.Date(historical[2]) - origin_date)
    start_future     <- as.numeric(as.Date(future[1]) - origin_date)
    end_future       <- as.numeric(as.Date(future[2]) - origin_date)
    
    # grab grids within our desired periods
    period_historical <- which(time_seq >= start_historical & time_seq <= end_historical)
    period_future     <- which(time_seq >= start_future & time_seq <= end_future)
    
    time_indices_to_keep <- c(period_historical, period_future)
    time_indices_to_keep <- sort(time_indices_to_keep)  # Keep in chronological order
    
    # message(sprintf("Filtered from %d to %d time steps", length(time_seq), length(time_indices_to_keep)))
    # message(sprintf("Historical period (%s to %s): %d steps, Future period (%s to %s): %d steps", 
    #                 historical[1], historical[2], length(period_historical),
    #                 future[1], future[2], length(period_future)))
  }
  
  # update time_seq and nt to reflect filtered data
  time_seq <- time_seq[time_indices_to_keep]
  nt <- length(time_seq)
  
  # convert sf water_supplies to terra vect object and extract coordinates
  v_pts <- vect(water_supplies)
  xy <- crds(v_pts)
  
  # convert geographic coordinates to array indices by finding nearest nc grid cell center
  lon_index <- map_int(xy[, 1], ~ which.min(abs(lon - .x)))
  # closest latitude index for each point
  lat_index <- map_int(xy[, 2], ~ which.min(abs(lat - .x)))
  # combine longitude and latitude indices into matrix for "closest nc grid cell":
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
  
  # message(sprintf("Extracting %d locations across %d time steps...", n_extractions, nt))
  
  # If we have filtered time indices, we need a different approach
  if (length(time_indices_to_keep) < length(time_seq)) {
    # Read all filtered time steps at once by reading the full array and subsetting
    # message("Reading filtered time data...")
    
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
    # if no filter, use original chunking approach (KW OG)
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
        # message(sprintf("Processed chunk %d / %d (time steps %d-%d)", chunk, n_chunks, start_t, end_t))
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
 
  # extract point attributes (non-spatial data) and add point ID column
  pts_attr <- st_drop_geometry(water_supplies) %>% mutate(.pt_id = row_number())
  
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
    date = as.Date(start_time_date) + as.numeric(time),  # Convert days since start_time_date
    time_index = NULL 
  )]
  
  
  # daily precip is sometimes in kg m-2 s-1 units. convert to mm/day
  if(mwbm_units != "mm"){
    message(mwbm_units)
    results_dt <- results_dt %>%
      mutate(value = value * 86400)
  }
  
  
  # precip is daily, must further aggregate to monthly
  if(var %in% c("pr_tavg", "PRCP", "pr")){
    
    results_dt <- results_dt[
      , date := as.Date(format(date, "%Y-%m-01"))
    ][
      , .(value = sum(value, na.rm = TRUE)),
      by = .(wsd_source_id, park_unit, time_units, nc_file, cell_lon, cell_lat, date)
    ]
    
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  if(var %in% c("Tmin", "Tmax")){
    write_parquet(results_dt, paste0(output_dir, basename(nc_path), "_", var, ".parquet"))
  } else {
    write_parquet(results_dt, paste0(output_dir, basename(nc_path), ".parquet"))
  }
  
  end_time <- Sys.time()
  
  message(sprintf("Extraction completed in %.2f seconds", as.numeric(difftime(end_time, start_time, units = "secs"))))
  
  nc_close(nc)
  
  gc()
}
