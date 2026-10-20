#' Process NetCDF files for LOCA2 water balance data extraction
#'
#' @param nc_path path to NetCDF files - found in wava_data_path("all", "loca2_mwbm/")
#' @param water_supplies sf object with point locations
#' @param var character, variable name to extract. Supported values include
#'   `"aet"`, `"deficit"`, `"pet"`, `"runoff"`, `"snow"`/`"swe"`,
#'   `"stor"`, `"pr"`, `"tas"`, `"tasmin_tavg"`, and `"tasmax_tavg"`.
#' @param historical Character vector of length 2 giving the start and end dates
#'   (in `"YYYY-MM-DD"` format) that defines the "historical" period.
#'   Defaults to `c("1980-01-01", "2024-12-31")`.
#' @param future Character vector of length 2 giving the start and end dates
#'   (in `"YYYY-MM-DD"` format) that defines the "future" period.
#'   Defaults to `c("2035-01-01", "2065-12-31")`.
#' @param output_dir Directory to save the extracted CSV output files.
#' @param grab_nearest Logical. If `TRUE`, sites whose original 3 x 3
#'   neighborhood contains no finite values are recentered on the nearest grid
#'   cell with finite data within 50 km. Defaults to `FALSE`.
#'
#' @export
process_loca2_mwbm <- function(nc_path,
                               water_supplies,
                               var,
                               historical = c("1980-01-01", "2024-12-31"),
                               future = c("2035-01-01", "2065-12-31"),
                               output_dir,
                               grab_nearest = FALSE){

  supported_vars <- c(
    "aet",
    "deficit",
    "pet",
    "runoff",
    "snow",
    "swe",
    "stor",
    "pr",
    "pr_tavg",
    "tas",
    "tasmin_tavg",
    "tasmax_tavg"
  )

  if (!var %in% supported_vars) {
    stop(
      "Unsupported var '", var, "'. Supported values are: ",
      paste(supported_vars, collapse = ", "), "."
    )
  }
  
  
  message(paste0("Starting vectorized extraction of ", nc_path))
  
  start_clock_time <- Sys.time()
  
  # open using ncdf4 package and ensure it is always closed on exit
  nc <- nc_open(nc_path)
  on.exit(nc_close(nc), add = TRUE)
  
  # Some files use different internal NetCDF variable names for the same
  # user-facing variable. Keep the original var label for later output and
  # aggregation decisions; nc_var is what gets passed to ncvar_get.
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
  
  # Guard: check that nc_var actually exists in this file before extraction.
  # This makes variable-name mismatches fail early with a useful message.
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
  # Alder training historical data is in 0-360 whereas others are not...
  if (any(lon > 180)) lon <- ifelse(lon > 180, lon - 360, lon)
  
  # transform water_supplies to nc coordinate system if not already
  if (st_crs(water_supplies) != st_crs(rast(nc_path))) {
    water_supplies <- st_transform(water_supplies, st_crs(rast(nc_path)))
  }
  
  # get time units
  # Some Livneh precipitation files use capitalized Time instead of time.
  time_seq <- if(nc_var == "PRCP") {
    ncvar_get(nc, "Time")
  } else {
    ncvar_get(nc, "time")
  }
  
  # Pull the matching time-unit metadata so numeric time values can be
  # converted back to dates below.
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
    
    # Convert date boundaries to the same "days since origin" units used by
    # the NetCDF time coordinate.
    origin_date <- as.Date(start_time_date)
    
    # Use the user-specified historical and future windows.
    start_historical <- as.numeric(as.Date(historical[1]) - origin_date)
    end_historical   <- as.numeric(as.Date(historical[2]) - origin_date)
    start_future     <- as.numeric(as.Date(future[1]) - origin_date)
    end_future       <- as.numeric(as.Date(future[2]) - origin_date)
    
    # Grab only layers that fall within the requested reporting periods.
    period_historical <- which(time_seq >= start_historical & time_seq <= end_historical)
    period_future     <- which(time_seq >= start_future & time_seq <= end_future)
    
    time_indices_to_keep <- sort(c(period_historical, period_future))
  }
  
  if (length(time_indices_to_keep) == 0) {
    message("No dates in ", basename(nc_path), " overlap the requested historical/future windows; skipping.")
    return(invisible(NULL))
  }
  
  # update time_seq and nt to reflect filtered data
  time_seq <- time_seq[time_indices_to_keep]
  nt <- length(time_seq)
  
  # convert sf water_supplies to terra vect object and extract coordinates
  v_pts <- vect(water_supplies)
  xy <- crds(v_pts)
  
  # Convert point coordinates to array indices by finding the nearest NetCDF
  # grid-cell center for each water supply.
  lon_index <- map_int(xy[, 1], ~ which.min(abs(lon - .x)))
  lat_index <- map_int(xy[, 2], ~ which.min(abs(lat - .x)))
  center_indices <- cbind(lon_index = lon_index, lat_index = lat_index)
  
  # Build a 3x3 neighborhood around each nearest grid cell
  # (-1, 0, +1 in both lon and lat directions).
  offsets <- expand.grid(lon_offset = -1:1, lat_offset = -1:1)
  
  neighbor_list <- vector("list", nrow(center_indices))
  
  for (i in 1:nrow(center_indices)) {
    # Get center longitude and latitude indices for the current point.
    lon_index <- center_indices[i, "lon_index"]
    lat_index <- center_indices[i, "lat_index"]
    
    # Apply offsets while clamping to valid grid boundaries.
    nbr_lons <- pmax(1, pmin(length(lon), lon_index + offsets$lon_offset))
    nbr_lats <- pmax(1, pmin(length(lat), lat_index + offsets$lat_offset))
    
    neighbor_list[[i]] <- data.frame(
      point_id = i,
      lon_index = nbr_lons,
      lat_index = nbr_lats
    )
  }
  
  neighbor_list <- bind_rows(neighbor_list) %>%
    # Clamping at raster edges can map multiple offsets to the same cell.
    # Keep each available grid cell once rather than weighting edge cells twice.
    distinct(point_id, lon_index, lat_index, .keep_all = TRUE) %>%
    # Add actual geographic coordinates for each neighborhood cell.
    mutate(
      cell_lon = lon[lon_index],
      cell_lat = lat[lat_index]
    )
  
  # Extract dimension names and sizes from the selected NetCDF variable.
  dnames <- map_chr(v$dim, "name")
  dsizes <- map_int(v$dim, "len")
  
  # Identify time dimension by looking for "time" in dimension names.
  time_dim <- which(grepl("time", dnames, ignore.case = TRUE))

  # Pre-allocate matrix for extracted values:
  # rows = point-neighborhood cells; columns = selected time steps.
  extract_cells <- function(cell_lookup) {
    n_extractions <- nrow(cell_lookup)
    result_matrix <- matrix(NA_real_, nrow = n_extractions, ncol = nt)

    # Use segment-based reading. Segments are contiguous runs of selected time
    # indices, split into chunks of at most 600 steps to cap memory use.
    if (length(time_indices_to_keep) > 1) {
      time_diffs <- diff(time_indices_to_keep)
      breaks <- which(time_diffs > 1)
    } else {
      breaks <- integer(0)
    }

    segment_starts <- c(1, breaks + 1)
    segment_ends <- c(breaks, length(time_indices_to_keep))
    result_position <- 1

    for (seg in seq_along(segment_starts)) {
      seg_start_pos <- segment_starts[seg]
      seg_end_pos <- segment_ends[seg]
      seg_length <- seg_end_pos - seg_start_pos + 1
      n_subchunks <- ceiling(seg_length / 600)

      for (sub in seq_len(n_subchunks)) {
        sub_start <- seg_start_pos + (sub - 1) * 600
        sub_end <- min(seg_start_pos + sub * 600 - 1, seg_end_pos)
        sub_length <- sub_end - sub_start + 1
        file_indices <- time_indices_to_keep[sub_start:sub_end]
        file_start <- file_indices[1]
        file_end <- file_indices[length(file_indices)]
        file_count <- file_end - file_start + 1
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

        if (length(time_dim) == 1) {
          if (time_dim == 1) {
            arr <- aperm(arr, c(2, 3, 1))
          } else if (time_dim == 2) {
            arr <- aperm(arr, c(1, 3, 2))
          } else if (time_dim != 3) {
            stop(paste("Unexpected time dimension position:", time_dim))
          }
        } else {
          stop("No time dimension found in NetCDF file")
        }

        for (i in seq_len(n_extractions)) {
          result_matrix[
            i,
            result_position:(result_position + sub_length - 1)
          ] <- arr[
            cell_lookup$lon_index[i],
            cell_lookup$lat_index[i],
            relative_indices
          ]
        }

        result_position <- result_position + sub_length
      }
    }

    result_matrix
  }

  result_matrix <- extract_cells(neighbor_list)

  if (grab_nearest) {
    point_has_data <- rowsum(
      as.integer(rowSums(is.finite(result_matrix)) > 0),
      group = neighbor_list$point_id,
      reorder = FALSE
    )[, 1] > 0
    missing_point_ids <- which(!point_has_data)

    if (length(missing_point_ids) > 0) {
      search_cells <- bind_rows(lapply(missing_point_ids, function(point_id) {
        point_lon <- xy[point_id, 1]
        point_lat <- xy[point_id, 2]
        lat_buffer <- 50 / 111.32
        lon_buffer <- 50 / (111.32 * max(abs(cos(point_lat * pi / 180)), 1e-6))
        candidate_lons <- which(abs(lon - point_lon) <= lon_buffer)
        candidate_lats <- which(abs(lat - point_lat) <= lat_buffer)

        if (!length(candidate_lons) || !length(candidate_lats)) return(NULL)

        candidates <- expand.grid(
          lon_index = candidate_lons,
          lat_index = candidate_lats
        )
        candidates$distance_m <- geosphere::distHaversine(
          cbind(lon[candidates$lon_index], lat[candidates$lat_index]),
          c(point_lon, point_lat)
        )
        candidates <- candidates[candidates$distance_m <= 50000, , drop = FALSE]
        if (!nrow(candidates)) return(NULL)
        candidates$point_id <- point_id
        candidates
      }))

      if (nrow(search_cells) > 0) {
        search_values <- extract_cells(search_cells)
        search_cells$has_data <- rowSums(is.finite(search_values)) > 0
        nearest_centers <- search_cells %>%
          filter(has_data) %>%
          group_by(point_id) %>%
          slice_min(distance_m, n = 1, with_ties = FALSE) %>%
          ungroup()

        if (nrow(nearest_centers) > 0) {
          fallback_neighbors <- bind_rows(lapply(seq_len(nrow(nearest_centers)), function(i) {
            point_id <- nearest_centers$point_id[i]
            center_lon <- nearest_centers$lon_index[i]
            center_lat <- nearest_centers$lat_index[i]
            data.frame(
              point_id = point_id,
              lon_index = pmax(1, pmin(length(lon), center_lon + offsets$lon_offset)),
              lat_index = pmax(1, pmin(length(lat), center_lat + offsets$lat_offset))
            )
          })) %>%
            distinct(point_id, lon_index, lat_index, .keep_all = TRUE) %>%
            mutate(
              cell_lon = lon[lon_index],
              cell_lat = lat[lat_index]
            )

          fallback_values <- extract_cells(fallback_neighbors)
          replace_point_ids <- unique(fallback_neighbors$point_id)
          keep_original <- !neighbor_list$point_id %in% replace_point_ids
          neighbor_list <- bind_rows(
            neighbor_list[keep_original, , drop = FALSE],
            fallback_neighbors
          )
          result_matrix <- rbind(
            result_matrix[keep_original, , drop = FALSE],
            fallback_values
          )
        }
      }
    }
  }

  n_extractions <- nrow(neighbor_list)

  # ============================================================================
  # Faster results_dt creation
  # rep.int() does not support each =, so use rep() for time_index.
  # result_matrix is unrolled column-major, which matches rep.int ordering
  # (all extraction IDs for time step 1, then time step 2, etc.).
  # ============================================================================
  
  # Convert the neighborhood lookup to data.table and add extraction IDs for
  # joining extracted values back to grid-cell metadata.
  neighbor_dt <- as.data.table(neighbor_list)
  neighbor_dt[, extraction_id := .I]
  
  # Extract point attributes (non-spatial data) for joining back to results.
  pts_attr <- as.data.table(sf::st_drop_geometry(water_supplies))
  pts_attr[, point_id := .I]
  
  # Convert matrix results to long format. as.vector(result_matrix) unrolls the
  # matrix by column, so extraction IDs repeat for each time step.
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
  
  # Precipitation can be stored in flux units (kg m-2 s-1, equivalent to
  # mm/sec for water). Convert to monthly depth using seconds/day and days/month.
  if(mwbm_units != "mm" & grepl("pr", var)){
    message(mwbm_units)
    
    results_dt[, ym := as.Date(format(date, "%Y-%m-01"))]
    results_dt[, days_in_month := lubridate::days_in_month(ym)]
    
    results_dt <- results_dt %>%
      mutate(value = value * 86400 * days_in_month)
    
  } else if(mwbm_units == "K"){
    
    # Temperature variables are stored in Kelvin; convert to Celsius.
    message(mwbm_units)
    
    results_dt <- results_dt %>%
      mutate(value = value - 273.15)
    
  }
  
  # Use original var (not nc_var) for aggregation check so the branch is
  # reliable regardless of what nc_var was renamed to
  # Livneh precipitation is daily, so aggregate daily values to monthly totals.
  # Use original var here because nc_var may have been renamed to PRCP.
  if(var %in% c("pr") & grepl("livneh", nc_path)){
    
    results_dt <- results_dt[
      , date := as.Date(format(date, "%Y-%m-01"))
    ][
      , .(value = if (any(is.finite(value))) {
        sum(value[is.finite(value)])
      } else {
        NA_real_
      }),
      by = .(wsd_source_id, park_unit, time_units, nc_file, cell_lon, cell_lat, date)
    ]
  }
  
  results_dt <- results_dt[
    ,
    .(wsd_source_id, park_unit, time_units, nc_file, cell_lat, cell_lon, date, value)
  ]
  
  # Create output directory if it does not exist.
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
