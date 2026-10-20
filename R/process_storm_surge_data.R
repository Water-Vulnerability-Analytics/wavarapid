#' Process Storm Surge Data (Muis 2023, 10-year return period)
#'
#' Reads every RP10 storm surge NetCDF file in `data_dir` for 4 CMIP6
#' models (CMCC-CM2-VHR4, EC-Earth3P-HR, HadGEM3-GC31-HM, HadGEM3-GC31-HM-SST),
#' across all available periods and all 3 confidence-interval variants, and
#' reshapes them into one row per station per model. Every `.nc` file in
#' `data_dir` is read, so the directory must contain only the `surge` variant
#' of the CDS download — the `waterlevel` variant (same models/periods/CI
#' variants, different physical quantity) must not be placed here.
#' This function does not compute a surge change or select a worst-case model
#' across the ensemble. Those calculations happen in `calc_storm_surge_exposure()`,
#' which calls this functions output (GeoPackage) as its input.
#'
#' Source: Muis et al. (2023), "Global Projections of Storm Surges Using
#' High-Resolution CMIP6 Climate Models", Earth's Future.
#' Data downloaded from Copernicus Climate Data Store:
#' https://cds.climate.copernicus.eu/datasets/sis-water-level-change-indicators-cmip6
#'   - Variable: surge level
#'   - Derived variable: absolute value
#'   - Product type: single model
#'   - Statistic: 10-year return period
#'   - Confidence interval: all 3 variants are read and kept — best-fit
#'     (central estimate, `_mean`), low-bound (`_p05`), and high-bound (`_p95`)
#'   - Periods: `historical` (1951-1980), `reference` (1985-2014, the
#'     baseline used for computing future change), and `future` (2021-2050) --
#'     all 3 must be present in `data_dir` for `calc_storm_surge_exposure()`/
#'     `calc_storm_surge_sensitivity()` downstream; this function does not
#'     currently validate that `reference` files were actually downloaded
#'   - Scenario: SSP5-8.5 (the only scenario available in this dataset)
#'   - Models: all 4 HighResMIP ensemble members, kept separate (not reduced
#'     to a worst-case model here)
#'
#' @param data_dir character, path to a directory containing only the RP10
#'   `surge` NetCDF files, all periods and CI variants (e.g.
#'   `data/all/coast_comp_data/storm_surge_muis_2023_RP10/`). Every `.nc` file
#'   in this directory is read.
#' @param save logical, whether to save output GeoPackage. Defaults to TRUE.
#' @param output_dir character, directory path where the GeoPackage will be saved.
#'   Required when `save = TRUE`. Defaults to `wava_data_path("rapid", "inputs", "storm_surge")`.
#' @param verbose logical, whether to show progress messages. Defaults to TRUE.
#'
#' @return sf object (EPSG:4326) of storm surge stations, one row per
#'   `station_idx` per `model`, with columns `historical_p05`/`historical_mean`/
#'   `historical_p95` (1951-1980), `reference_p05`/`reference_mean`/`reference_p95`
#'   (1985-2014), and `future_p05`/`future_mean`/`future_p95` (2021-2050), all
#'   in meters
#'
#' @export
process_storm_surge_data <- function(data_dir,
                                     save = TRUE,
                                     output_dir = wava_data_path("rapid", "inputs", "storm_surge"),
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
