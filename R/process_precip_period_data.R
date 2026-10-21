#' @title Extract dry and wet spell lengths at water supplies
#'
#' @description
#' Reads one NetCDF file of the LOCA2 precipitation-extremes collection and
#' writes the yearly maximum spell lengths at each water supply to a parquet
#' file. Run it once per file; [calc_dry_period_exposure()],
#' [calc_wet_period_exposure()], [calc_dry_period_sensitivity()] and
#' [calc_wet_period_sensitivity()] read the folder.
#'
#' Each value is the mean over the 3 x 3 block of grid cells around the
#' supply's nearest cell (the same block as [process_loca2_mwbm()]), ignoring
#' missing cells. The two variables use a 2.54 mm/day threshold:
#' - `dry_spell_max_length`: the year's longest run of dry days (CDD);
#' - `wet_spell_max_length`: the year's longest run of wet days (CWD).
#'
#' The model, scenario and run are read from the file name, which must follow
#' `MODEL.SCENARIO.RUN.YEARS.LOCA_16thdeg_pr_metrics.nc`.
#'
#' The parquet file is `<nc file name>.parquet` in `output_dir`, with one row
#' per water supply and year and columns `wsd_source_id`, `park_unit`,
#' `year`, `model`, `ssp`, `run`, `nc_file`, `dry_spell_max_length` and
#' `wet_spell_max_length`.
#'
#' @param nc_path String. Path to one precipitation-metrics NetCDF file, from
#'   `all/loca_cmip6_precip_products` in the shared data folder. Stops with an
#'   error if the file doesn't exist, lacks either variable, or has a name
#'   that can't be parsed. Required.
#' @param water_supplies sf object. The water supplies, as `POINT` geometries
#'   with `wsd_source_id` and `park_unit` columns, e.g. from
#'   [wavacore::get_wsd()]. Stops with an error if any point is outside the
#'   LOCA2 grid. Required.
#' @param output_dir String. The folder for the parquet file, created if it
#'   doesn't exist. Default:
#'   `rapid/inputs/loca_cmip6_precip_periods_water_supplies` in the shared
#'   data folder, see [wavacore::wava_data_path()].
#' @param overwrite Logical. What to do when the output file already exists.
#'   * `TRUE`: processes the file again and replaces it.
#'   * `FALSE`: skips the file and returns the existing path.
#'
#'   Default: `FALSE`.
#' @param verbose Logical. Whether to print progress.
#'   * `TRUE`: prints the file being read, each variable, and the run time.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return The path of the parquet file, invisibly.
#' @export
#' @seealso [calc_dry_period_exposure()], [calc_wet_period_exposure()],
#'   [calc_dry_period_sensitivity()] and [calc_wet_period_sensitivity()],
#'   which read the parquet files; [process_loca2_mwbm()] for the other LOCA2
#'   variables.
#' @examplesIf interactive() && dir.exists(wavacore::wava_data_path("all", "loca_cmip6_precip_products"))
#' # One file is about 0.7 GB, so this takes a few minutes
#' ws <- wavacore::get_wsd(park = "ZION")
#' nc_name <- "ACCESS-CM2.historical.r1i1p1f1.1950-2014.LOCA_16thdeg_pr_metrics.nc"
#' nc <- wavacore::wava_data_path("all", "loca_cmip6_precip_products", "historical", nc_name)
#' out_file <- process_precip_period_data(nc, ws, output_dir = tempdir())
#' head(arrow::read_parquet(out_file))

process_precip_period_data <- function(
    nc_path,
    water_supplies,
    output_dir = wavacore::wava_data_path("rapid", "inputs", "loca_cmip6_precip_periods_water_supplies/"),
    overwrite = FALSE,
    verbose = TRUE) {

  required_variables <- c("dry_spell_max_length", "wet_spell_max_length")

  stopifnot(
    "nc_path must identify one existing file" = length(nc_path) == 1 &&
      file.exists(nc_path),
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "water_supplies must include wsd_source_id" =
      "wsd_source_id" %in% names(water_supplies),
    "water_supplies must include park_unit" =
      "park_unit" %in% names(water_supplies)
  )

  if (!all(sf::st_geometry_type(water_supplies) %in% "POINT")) {
    stop("water_supplies must contain POINT geometries.")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  output_path <- file.path(output_dir, paste0(basename(nc_path), ".parquet"))

  if (file.exists(output_path) && !overwrite) {
    if (verbose) message("Output already exists; skipping: ", output_path)
    return(invisible(output_path))
  }

  if (verbose) {
    message("Processing precipitation-period metrics from ", basename(nc_path))
    start_time <- Sys.time()
  }

  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  missing_variables <- setdiff(required_variables, names(nc$var))

  if (length(missing_variables)) {
    stop("NetCDF is missing required variables: ",
      paste(missing_variables, collapse = ", "))
  }

  # Parse model metadata from the fixed dataset filename convention:
  # MODEL.SCENARIO.RUN.YEARS.LOCA_16thdeg_pr_metrics.nc
  filename_parts <- strsplit(basename(nc_path), "\\.")[[1]]

  if (length(filename_parts) < 3L) {
    stop("Could not parse model, scenario, and run from: ", basename(nc_path))
  }

  model_id <- filename_parts[1]
  ssp_id <- filename_parts[2]
  run_id <- filename_parts[3]

  # Read the coordinate axes directly from NetCDF, following the faster array
  # extraction approach used by process_loca2_mwbm().
  lon <- ncdf4::ncvar_get(nc, "lon")
  lat <- ncdf4::ncvar_get(nc, "lat")

  template_variable <- nc$var[[required_variables[1]]]
  template_dim_names <- vapply(
    template_variable$dim,
    function(x) x$name,
    character(1)
  )
  template_dim_sizes <- vapply(
    template_variable$dim,
    function(x) x$len,
    integer(1)
  )

  lon_dim <- which(tolower(template_dim_names) %in% c("lon", "longitude"))
  lat_dim <- which(tolower(template_dim_names) %in% c("lat", "latitude"))
  time_dim <- which(grepl("time", template_dim_names, ignore.case = TRUE))

  if (length(lon_dim) != 1L || length(lat_dim) != 1L || length(time_dim) != 1L) {
    stop(
      "Could not identify one longitude, latitude, and time dimension for ",
      required_variables[1], "."
    )
  }

  time_name <- template_dim_names[time_dim]
  time_values <- ncdf4::ncvar_get(nc, time_name)
  time_units <- nc$dim[[time_name]]$units

  if (
    is.null(time_units) ||
      length(time_units) != 1L ||
      is.na(time_units) ||
      !grepl("^days since\\s+", time_units, ignore.case = TRUE)
  ) {
    stop("Expected the time coordinate to use 'days since ...' units.")
  }

  time_origin_text <- trimws(
    sub("^.*since\\s+", "", time_units, ignore.case = TRUE)
  )
  time_origin <- as.Date(sub("[ T].*$", "", time_origin_text))
  layer_dates <- time_origin + as.numeric(time_values)

  if (
    length(layer_dates) != template_dim_sizes[time_dim] ||
      anyNA(layer_dates)
  ) {
    stop("Could not recover one annual date for every NetCDF time step.")
  }

  # The NetCDF longitude axis is 0-360. Convert supply-point longitudes to the
  # same convention before locating their nearest grid-cell centers.
  supply_points <- sf::st_transform(water_supplies, 4326)
  point_xy <- sf::st_coordinates(supply_points)

  if (any(lon > 180)) {
    point_xy[, 1] <- ifelse(
      point_xy[, 1] < 0,
      point_xy[, 1] + 360,
      point_xy[, 1]
    )
  }

  lon_half_step <- stats::median(abs(diff(lon))) / 2
  lat_half_step <- stats::median(abs(diff(lat))) / 2
  outside_grid <-
    point_xy[, 1] < min(lon) - lon_half_step |
    point_xy[, 1] > max(lon) + lon_half_step |
    point_xy[, 2] < min(lat) - lat_half_step |
    point_xy[, 2] > max(lat) + lat_half_step

  if (any(outside_grid)) {
    outside_ids <- water_supplies$wsd_source_id[outside_grid]
    stop(
      "Water supplies fall outside the LOCA2 grid: ",
      paste(outside_ids, collapse = ", ")
    )
  }

  # At an exact cell boundary, use the eastern and southern cell centers to
  # retain the prior Terra extraction's boundary behavior.
  center_indices <- cbind(
    lon_index = vapply(
      point_xy[, 1],
      function(x) {
        distance <- abs(lon - x)
        candidates <- which(distance == min(distance))
        candidates[which.max(lon[candidates])]
      },
      integer(1)
    ),
    lat_index = vapply(
      point_xy[, 2],
      function(x) {
        distance <- abs(lat - x)
        candidates <- which(distance == min(distance))
        candidates[which.min(lat[candidates])]
      },
      integer(1)
    )
  )

  # Build the same 3 x 3 index neighborhood used by process_loca2_mwbm().
  offsets <- expand.grid(lon_offset = -1:1, lat_offset = -1:1)
  neighbor_lookup <- data.table::rbindlist(
    lapply(seq_len(nrow(center_indices)), function(point_id) {
      data.table::data.table(
        point_id = point_id,
        lon_index = pmax(
          1L,
          pmin(
            length(lon),
            center_indices[point_id, "lon_index"] + offsets$lon_offset
          )
        ),
        lat_index = pmax(
          1L,
          pmin(
            length(lat),
            center_indices[point_id, "lat_index"] + offsets$lat_offset
          )
        )
      )
    })
  )
  neighbor_lookup <- unique(
    neighbor_lookup,
    by = c("point_id", "lon_index", "lat_index")
  )

  n_points <- nrow(water_supplies)
  n_layers <- length(time_values)

  processed <- data.table::data.table(
    point_id = rep(seq_len(n_points), times = n_layers),
    year = rep(as.integer(format(layer_dates, "%Y")), each = n_points)
  )

  for (variable in required_variables) {
    if (verbose) message("  Extracting ", variable)

    metric_variable <- nc$var[[variable]]
    dim_names <- vapply(
      metric_variable$dim,
      function(x) x$name,
      character(1)
    )
    dim_sizes <- vapply(
      metric_variable$dim,
      function(x) x$len,
      integer(1)
    )
    variable_lon_dim <- which(tolower(dim_names) %in% c("lon", "longitude"))
    variable_lat_dim <- which(tolower(dim_names) %in% c("lat", "latitude"))
    variable_time_dim <- which(grepl("time", dim_names, ignore.case = TRUE))
    extra_dims <- setdiff(
      seq_along(dim_names),
      c(variable_lon_dim, variable_lat_dim, variable_time_dim)
    )

    if (
      length(variable_lon_dim) != 1L ||
        length(variable_lat_dim) != 1L ||
        length(variable_time_dim) != 1L ||
        dim_sizes[variable_lon_dim] != length(lon) ||
        dim_sizes[variable_lat_dim] != length(lat) ||
        dim_sizes[variable_time_dim] != n_layers ||
        any(dim_sizes[extra_dims] != 1L)
    ) {
      stop(
        variable,
        " does not have the expected lon, lat, time, and singleton extra dimensions."
      )
    }

    n_extractions <- nrow(neighbor_lookup)
    neighbor_values <- matrix(
      NA_real_,
      nrow = n_extractions,
      ncol = n_layers
    )

    # Read contiguous time chunks directly into an array, reorder to
    # lon-lat-time, and extract all neighborhood indices. Annual files fit in
    # one chunk, while the cap retains the MWBM processor's memory safeguard.
    max_chunk_size <- 600L
    n_chunks <- ceiling(n_layers / max_chunk_size)

    for (chunk_id in seq_len(n_chunks)) {
      chunk_start <- (chunk_id - 1L) * max_chunk_size + 1L
      chunk_end <- min(chunk_id * max_chunk_size, n_layers)
      chunk_length <- chunk_end - chunk_start + 1L

      start <- rep(1L, length(dim_sizes))
      count <- dim_sizes
      start[variable_time_dim] <- chunk_start
      count[variable_time_dim] <- chunk_length

      metric_array <- ncdf4::ncvar_get(
        nc,
        variable,
        start = start,
        count = count,
        collapse_degen = FALSE
      )

      dimension_order <- c(
        variable_lon_dim,
        variable_lat_dim,
        variable_time_dim,
        extra_dims
      )
      metric_array <- aperm(metric_array, dimension_order)
      dim(metric_array) <- c(length(lon), length(lat), chunk_length)

      for (extraction_id in seq_len(n_extractions)) {
        neighbor_values[
          extraction_id,
          chunk_start:chunk_end
        ] <- metric_array[
          neighbor_lookup$lon_index[extraction_id],
          neighbor_lookup$lat_index[extraction_id],
          seq_len(chunk_length)
        ]
      }
    }

    neighbor_values[!is.finite(neighbor_values)] <- NA_real_

    neighbor_sums <- rowsum(replace(neighbor_values, is.na(neighbor_values), 0),
      group = neighbor_lookup$point_id,
      reorder = FALSE)

    neighbor_counts <- rowsum(1L * !is.na(neighbor_values),
      group = neighbor_lookup$point_id,
      reorder = FALSE)

    point_means <- neighbor_sums / neighbor_counts
    point_means[neighbor_counts == 0] <- NA_real_

    data.table::set(processed,
      j = variable,
      value = as.vector(point_means))

    rm(metric_array, neighbor_values)
    invisible(gc())
  }

  point_attributes <- data.table::as.data.table(
    sf::st_drop_geometry(water_supplies)[
      ,
      c("wsd_source_id", "park_unit"),
      drop = FALSE
    ]
  )
  point_attributes[, point_id := .I]

  processed[
    point_attributes,
    on = "point_id",
    `:=`(
      wsd_source_id = i.wsd_source_id,
      park_unit = i.park_unit
    )
  ]

  processed[
    ,
    `:=`(
      nc_file = basename(nc_path),
      model = model_id,
      ssp = ssp_id,
      run = run_id,
      point_id = NULL
    )
  ]

  data.table::setcolorder(processed, c("wsd_source_id", "park_unit", "year", "model", "ssp", "run", "nc_file", "dry_spell_max_length", "wet_spell_max_length"))

  arrow::write_parquet(processed, output_path)

  if (verbose) {
    elapsed <- Sys.time() - start_time
    message(
      "Precipitation-period extraction completed in ",
      round(elapsed, 2), " ", attr(elapsed, "units"),
      ": ", output_path
    )
  }

  invisible(output_path)
}
