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
  
  start_clock_time <- Sys.time()
  
  # open using ncdf4 package and ensure it is always closed on exit
  nc <- nc_open(nc_path)
  on.exit(nc_close(nc), add = TRUE)
  
  # Keep original var label for output path decisions;
  # nc_var is what gets passed to ncvar_get and may be renamed per-file.
  nc_var <- var
  
  if(nc_var == "pr" & !grepl("alder", nc_path, ignore.case = TRUE)){
    nc_var <- "pr_tavg"
  }
  
  if(nc_var == "pr_tavg" & grepl("livneh_unsplit_precip", nc_path, ignore.case = TRUE)){
    nc_var <- "PRCP"
  }
  
  if(nc_var == "tasmin_tavg" & grepl("lusu", nc_path, ignore.case = TRUE)){
    nc_var <- "Tmin"
  }
  
  if(nc_var == "tasmax_tavg" & grepl("lusu", nc_path, ignore.case = TRUE)){
    nc_var <- "Tmax"
  }
  
  if(nc_var == "swe"){
    nc_var <- "snow"
  }
  
  # Guard: check that nc_var actually exists in this file
  available_vars <- names(nc$var)
  if(!nc_var %in% available_vars){
    stop(sprintf(
      "Variable '%s' (mapped from '%s') not found in %s.\nAvailable variables: %s",
      nc_var, var, basename(nc_path), paste(available_vars, collapse = ", ")
    ))
  }
  
  # grab layer of user-selected variable
  v <- nc$var[[nc_var]]
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  
  # If longitudes are 0-360, convert to -180-180
  if (any(lon > 180)) lon <- ifelse(lon > 180, lon - 360, lon)
  
  # transform water_supplies to nc coordinate system if not already
  if (st_crs(water_supplies) != st_crs(rast(nc_path))) {
    water_supplies <- st_transform(water_supplies, st_crs(rast(nc_path)))
  }
  
  # get time units
  time_seq <- if(nc_var == "PRCP") {
    ncvar_get(nc, "Time")
  } else {
    ncvar_get(nc, "time")
  }
  
  time_units <- if(nc_var == "PRCP") {
    nc$dim[["Time"]]$units
  } else {
    nc$dim[["time"]]$units
  }
  
  start_time_date <- lubridate::ymd(
    sub(" .*", "", sub(".*since\\s+", "", time_units))
  )
  
  # get mwbm units
  mwbm_units <- ncatt_get(nc, nc_var, "units")$value
  
  # Filter time indices based on units
  time_indices_to_keep <- seq_along(time_seq)
  
  if (!is.na(time_units) && grepl("days since", time_units, ignore.case = TRUE)) {
    
    origin_date <- as.Date(start_time_date)
    
    start_historical <- as.numeric(as.Date(historical[1]) - origin_date)
    end_historical   <- as.numeric(as.Date(historical[2]) - origin_date)
    start_future     <- as.numeric(as.Date(future[1]) - origin_date)
    end_future       <- as.numeric(as.Date(future[2]) - origin_date)
    
    period_historical <- which(time_seq >= start_historical & time_seq <= end_historical)
    period_future     <- which(time_seq >= start_future & time_seq <= end_future)
    
    time_indices_to_keep <- sort(c(period_historical, period_future))
  }
  
  # update time_seq and nt to reflect filtered data
  time_seq <- time_seq[time_indices_to_keep]
  nt <- length(time_seq)
  
  # convert sf water_supplies to terra vect object and extract coordinates
  v_pts <- vect(water_supplies)
  xy <- crds(v_pts)
  
  lon_index <- map_int(xy[, 1], ~ which.min(abs(lon - .x)))
  lat_index <- map_int(xy[, 2], ~ which.min(abs(lat - .x)))
  center_indices <- cbind(lon_index = lon_index, lat_index = lat_index)
  
  offsets <- expand.grid(lon_offset = -1:1, lat_offset = -1:1)
  
  neighbor_list <- vector("list", nrow(center_indices))
  
  for (i in 1:nrow(center_indices)) {
    lon_index <- center_indices[i, "lon_index"]
    lat_index <- center_indices[i, "lat_index"]
    
    nbr_lons <- pmax(1, pmin(length(lon), lon_index + offsets$lon_offset))
    nbr_lats <- pmax(1, pmin(length(lat), lat_index + offsets$lat_offset))
    
    neighbor_list[[i]] <- data.frame(
      point_id = i,
      lon_index = nbr_lons,
      lat_index = nbr_lats
    )
  }
  
  neighbor_list <- bind_rows(neighbor_list) %>%
    mutate(
      cell_lon = lon[lon_index],
      cell_lat = lat[lat_index]
    )
  
  dnames <- map_chr(v$dim, "name")
  dsizes <- map_int(v$dim, "len")
  
  time_dim <- which(grepl("time", dnames, ignore.case = TRUE))
  
  n_extractions <- nrow(neighbor_list)
  result_matrix <- matrix(NA_real_, nrow = n_extractions, ncol = nt)
  
  # ============================================================================
  # OPTIMIZED EXTRACTION: Use segment-based reading
  # Segments are contiguous runs of time_indices_to_keep. Each segment is then
  # split into sub-chunks of at most max_chunk_size steps to cap memory usage.
  # ============================================================================
  
  if (length(time_indices_to_keep) > 1) {
    time_diffs <- diff(time_indices_to_keep)
    breaks <- which(time_diffs > 1)
  } else {
    breaks <- integer(0)
  }
  
  segment_starts <- c(1, breaks + 1)
  segment_ends <- c(breaks, length(time_indices_to_keep))
  n_segments <- length(segment_starts)
  
  max_chunk_size <- 600
  result_position <- 1
  
  for (seg in 1:n_segments) {
    seg_start_pos <- segment_starts[seg]
    seg_end_pos <- segment_ends[seg]
    seg_length <- seg_end_pos - seg_start_pos + 1
    
    n_subchunks <- ceiling(seg_length / max_chunk_size)
    
    for (sub in 1:n_subchunks) {
      sub_start  <- seg_start_pos + (sub - 1) * max_chunk_size
      sub_end    <- min(seg_start_pos + sub * max_chunk_size - 1, seg_end_pos)
      sub_length <- sub_end - sub_start + 1
      
      # The actual file-level time indices for this sub-chunk
      file_indices <- time_indices_to_keep[sub_start:sub_end]
      file_start   <- file_indices[1]
      file_end     <- file_indices[length(file_indices)]
      file_count   <- file_end - file_start + 1
      
      # Relative positions within the array to extract.
      # Handles any gaps in file_indices so we never mis-align values
      # into the wrong output time slots.
      relative_indices <- file_indices - file_start + 1
      
      start <- rep(1L, length(dsizes))
      count <- dsizes
      
      if (length(time_dim) == 1) {
        start[time_dim] <- file_start
        count[time_dim] <- file_count
      }
      
      arr <- ncvar_get(
        nc,
        nc_var,
        start = start,
        count = count,
        collapse_degen = FALSE
      )
      
      # Reorder dimensions to (lon, lat, time)
      if (length(time_dim) == 1) {
        if (time_dim == 1) {
          arr <- aperm(arr, c(2, 3, 1))
        } else if (time_dim == 2) {
          arr <- aperm(arr, c(1, 3, 2))
        } else if (time_dim == 3) {
          # already in correct order
        } else {
          stop(paste("Unexpected time dimension position:", time_dim))
        }
      } else {
        stop("No time dimension found in NetCDF file")
      }
      
      # Extract values for all neighborhood cells from this chunk.
      # Using relative_indices ensures we only pull the time steps we want,
      # even if file_count > sub_length due to gaps within the range.
      for (i in seq_len(n_extractions)) {
        lon_idx <- neighbor_list$lon_index[i]
        lat_idx <- neighbor_list$lat_index[i]
        
        result_matrix[
          i,
          result_position:(result_position + sub_length - 1)
        ] <- arr[lon_idx, lat_idx, relative_indices]
      }
      
      result_position <- result_position + sub_length
    }
  }
  
  # ============================================================================
  # Faster results_dt creation
  # rep.int() does not support each =, so use rep() for time_index.
  # result_matrix is unrolled column-major, which matches rep.int ordering
  # (all extraction IDs for time step 1, then time step 2, etc.).
  # ============================================================================
  
  neighbor_dt <- as.data.table(neighbor_list)
  neighbor_dt[, extraction_id := .I]
  
  pts_attr <- as.data.table(sf::st_drop_geometry(water_supplies))
  pts_attr[, point_id := .I]
  
  results_dt <- data.table(
    extraction_id = rep.int(seq_len(n_extractions), times = nt),
    time_index = rep(seq_len(nt), each = n_extractions),
    value = as.vector(result_matrix)
  )
  
  results_dt[
    neighbor_dt,
    on = "extraction_id",
    `:=`(
      point_id = i.point_id,
      cell_lon = i.cell_lon,
      cell_lat = i.cell_lat
    )
  ]
  
  results_dt[
    pts_attr,
    on = "point_id",
    `:=`(
      wsd_source_id = i.wsd_source_id,
      park_unit = i.park_unit
    )
  ]
  
  results_dt[, `:=`(
    time_units = time_units,
    nc_file = basename(nc_path),
    cell_lon = as.character(cell_lon),
    cell_lat = as.character(cell_lat),
    date = as.Date(start_time_date) + as.numeric(time_seq[time_index]),
    extraction_id = NULL,
    time_index = NULL,
    point_id = NULL
  )]
  
  if(mwbm_units != "mm" & grepl("pr", var)){
    message(mwbm_units)
    
    results_dt[, ym := as.Date(format(date, "%Y-%m-01"))]
    results_dt[, days_in_month := lubridate::days_in_month(ym)]
    
    results_dt <- results_dt %>%
      mutate(value = value * 86400 * days_in_month)
    
  } else if(mwbm_units == "K"){
    
    message(mwbm_units)
    
    results_dt <- results_dt %>%
      mutate(value = value - 273.15)
    
  }
  
  # Use original var (not nc_var) for aggregation check so the branch is
  # reliable regardless of what nc_var was renamed to
  if(var %in% c("pr") & grepl("livneh", nc_path)){
    
    results_dt <- results_dt[
      , date := as.Date(format(date, "%Y-%m-01"))
    ][
      , .(value = sum(value, na.rm = TRUE)),
      by = .(wsd_source_id, park_unit, time_units, nc_file, cell_lon, cell_lat, date)
    ]
  }
  
  results_dt <- results_dt[
    ,
    .(wsd_source_id, park_unit, time_units, nc_file, cell_lat, cell_lon, date, value)
  ]
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Use nc_var for output naming since it reflects the actual variable extracted
  if(nc_var %in% c("Tmin", "Tmax")){
    output_path <- file.path(output_dir, paste0(basename(nc_path), "_", nc_var, ".parquet"))
  } else {
    output_path <- file.path(output_dir, paste0(basename(nc_path), ".parquet"))
  }
  
  write_parquet(results_dt, output_path)
  
  end_clock_time <- Sys.time()
  
  message(sprintf(
    "Extraction completed in %.2f seconds",
    as.numeric(difftime(end_clock_time, start_clock_time, units = "secs"))
  ))
  
  # nc_close() is now handled by on.exit() at the top of the function
}