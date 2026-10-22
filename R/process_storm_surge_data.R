#' @title Prepare the storm surge stations (Muis et al. 2023)
#'
#' @description
#' Reads every 10-year return period storm surge NetCDF file in `data_dir`
#' (4 CMIP6 models: CMCC-CM2-VHR4, EC-Earth3P-HR, HadGEM3-GC31-HM and
#' HadGEM3-GC31-HM-SST; 3 periods; 3 confidence-interval variants) and reshapes
#' them into one row per station and model. It doesn't calculate a change or
#' pick a model; [calc_storm_surge_exposure()] and
#' [calc_storm_surge_sensitivity()] do that from its output.
#'
#' The files come from the Copernicus Climate Data Store
#' (<https://cds.climate.copernicus.eu/datasets/sis-water-level-change-indicators-cmip6>):
#' variable surge level, derived variable absolute value, product type single
#' model, statistic 10-year return period, scenario SSP5-8.5 (the only one
#' available), all 3 confidence intervals (best fit, low bound, high bound)
#' and all 3 periods: historical (1951-1980), reference (1985-2014, the
#' baseline for the future change) and future (2021-2050). The function
#' doesn't check that the reference files are there.
#'
#' Source: Muis et al. (2023), "Global Projections of Storm Surges Using
#' High-Resolution CMIP6 Climate Models", Earth's Future.
#'
#' @param data_dir String. The folder with the RP10 `surge` NetCDF files.
#'   Every `.nc` file in it is read, so it must not also hold the `waterlevel`
#'   variant of the download (same models and periods, a different quantity).
#'   The model, period and confidence interval are read from each file name.
#'   Stops with an error if the folder doesn't exist or has no `.nc` files.
#'   Required.
#' @param save Logical. Whether to write the stations to a GeoPackage.
#'   * `TRUE`: writes `storm_surge_RP10_muis_2023.gpkg` to `output_dir`,
#'     creating the folder if needed and replacing an existing file. Stops
#'     with an error if `output_dir` is `NULL`.
#'   * `FALSE`: writes nothing.
#'
#'   Default: `TRUE`.
#' @param output_dir String or `NULL`. The folder for the GeoPackage.
#'   * A folder path: where the file is written when `save = TRUE`. Ignored
#'     when `save = FALSE`.
#'   * `NULL`: only allowed with `save = FALSE`.
#'
#'   Default: `rapid/inputs/storm_surge` in the shared data folder, see
#'   [wavacore::wava_data_path()].
#' @param verbose Logical. Whether to print progress.
#'   * `TRUE`: prints messages for each step, the run time and the path of the
#'     saved file.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return An sf object of points (CRS 4326), one row per station and model,
#'   with columns `station_idx`, `model`, `historical_p05`,
#'   `historical_mean`, `historical_p95` (1951-1980), `reference_p05`,
#'   `reference_mean`, `reference_p95` (1985-2014), `future_p05`,
#'   `future_mean`, `future_p95` (2021-2050), `lat` and `lon`. Surge levels
#'   are in meters; `_mean` is the best fit and `_p05`/`_p95` the low and
#'   high bounds.
#' @export
#' @seealso [calc_storm_surge_exposure()] and
#'   [calc_storm_surge_sensitivity()], which read the GeoPackage.
#' @examplesIf dir.exists(wavacore::wava_data_path("all", "coast_comp_data", "storm_surge_muis_2023_RP10"))
#' # Reshape the NetCDF files without writing the GeoPackage
#' stations <- process_storm_surge_data(
#'   data_dir = wavacore::wava_data_path("all", "coast_comp_data", "storm_surge_muis_2023_RP10"),
#'   save = FALSE
#' )
#' table(stations$model)

process_storm_surge_data <- function(data_dir,
                                     save = TRUE,
                                     output_dir = wavacore::wava_data_path("rapid", "inputs", "storm_surge"),
                                     verbose = TRUE) {
  
  # Input validation
  stopifnot(
    "data_dir does not exist"                   = dir.exists(data_dir),
    "save is TRUE but output_dir not provided"  = !save || !is.null(output_dir)
  )

  if (verbose) {
    message("Starting storm surge data processing...")
    start_time <- Sys.time()
  }

  # List all RP10 surge NetCDF files (directory must hold only the surge variant)
  all_files <- list.files(data_dir, pattern = "\\.nc$", full.names = TRUE)

  if (length(all_files) == 0) {
    stop("No RP10 surge NetCDF files found in: ", data_dir)
  }

  if (verbose) {
    message("Found ", length(all_files), " RP10 surge NetCDF files in: ", data_dir)
    message("Generating station index and coordinate key...")
  }

  # Generate station index and coordinate key.
  sample_nc <- ncdf4::nc_open(all_files[1])
  
  station_idx_coords <- tibble::tibble(
    station_idx = as.integer(ncdf4::ncvar_get(sample_nc, "stations")),
    lat = ncdf4::ncvar_get(sample_nc, "station_y_coordinate"),
    lon = ncdf4::ncvar_get(sample_nc, "station_x_coordinate")
  )
  
  ncdf4::nc_close(sample_nc)

  if (verbose) {
    message("Reading and reshaping ", 
            length(all_files), 
            " files across models, periods, and confidence intervals...")
  }

  # Generate a sf for each stations' model values across the confidence intervals.
  stations_sf <- purrr::map(all_files, function(f){

    # Use file name to get metadata for df
    f_metadata <- stringr::str_split_1(basename(f), "_")

    # Grab data from nc file for df
    nc_data <- ncdf4::nc_open(f)

    surge_vec <- switch(
      EXPR = f_metadata[7],
      "best-fit" = ncdf4::ncvar_get(nc_data, "return_mean_surge_level"),
      "high-bound-confidence-interval" = ncdf4::ncvar_get(nc_data, "return_95_surge_level"),
      "low-bound-confidence-interval" = ncdf4::ncvar_get(nc_data, "return_5_surge_level")
    )

    ncdf4::nc_close(nc_data)

    # Generate df for surge data
    tibble::tibble(
      station_idx    = station_idx_coords$station_idx,
      model          = f_metadata[1],
      period_variant = paste(f_metadata[5], f_metadata[7], sep = "_"),
      value          = surge_vec
    )
    
  }) |>
    purrr::compact() |>
    dplyr::bind_rows() |>
    tidyr::pivot_wider(names_from = period_variant, values_from = value) |>
    dplyr::left_join(station_idx_coords, by = "station_idx") |>
    # Order and rename the columns
    dplyr::select(
      station_idx, model,
      historical_p05  = `1951-1980_low-bound-confidence-interval`,
      historical_mean = `1951-1980_best-fit`,
      historical_p95  = `1951-1980_high-bound-confidence-interval`,
      reference_p05   = `1985-2014_low-bound-confidence-interval`,
      reference_mean  = `1985-2014_best-fit`,
      reference_p95   = `1985-2014_high-bound-confidence-interval`,
      future_p05      = `2021-2050_low-bound-confidence-interval`,
      future_mean     = `2021-2050_best-fit`,
      future_p95      = `2021-2050_high-bound-confidence-interval`,
      lat, lon
    ) |>
    # Convert to sf
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

  if (verbose) {
    message("Reshaped to ", nrow(stations_sf), " station-model rows...")
  }

  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    out_path <- file.path(output_dir, "storm_surge_RP10_muis_2023.gpkg")
    sf::st_write(stations_sf, out_path, delete_dsn = TRUE, quiet = TRUE)

    if (verbose) {
      message("Storm surge stations saved to: ", out_path)
    }
  }

  if (verbose) {
    total_time <- Sys.time() - start_time
    message(sprintf("Storm surge data processing complete in %.1f %s.",
                    as.numeric(total_time), attr(total_time, "units")))
  }

  return(stations_sf)
}
