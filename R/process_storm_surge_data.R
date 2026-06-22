#' Process Storm Surge Data (Muis 2023, 100-year return period)
#'
#' Reads the 100-year return period (RP100) storm surge NetCDF files for 4 CMIP6
#' models (CMCC-CM2-VHR4, EC-Earth3P-HR, HadGEM3-GC31-HM, HadGEM3-GC31-HM-SST),
#' computes the change in RP100 absolute surge level (future 2021-2050 minus
#' historical 1985-2014) for each station and each model, filters to CONUS
#' stations, and retains the model with the highest change per station
#' (worst-case across the HighResMIP ensemble).
#' Saves a GeoPackage for use by `calc_storm_surge_exposure()`.
#'
#' Source: Muis et al. (2023), "Global Projections of Storm Surges Using
#' High-Resolution CMIP6 Climate Models", Earth's Future.
#' Data downloaded from Copernicus Climate Data Store:
#' https://cds.climate.copernicus.eu/datasets/sis-water-level-change-indicators-cmip6
#'   - Variable: surge level
#'   - Derived variable: absolute value
#'   - Product type: single model
#'   - Statistic: 100-year return period
#'   - Confidence interval: best-fit (central estimate); low- and high-bound
#'     files are also available in the data directory but not used here
#'   - Experiment: historical (1985-2014) and future (2021-2050)
#'   - Scenario: SSP5-8.5 (the only scenario available in this dataset)
#'   - Models: all 4 HighResMIP ensemble members
#'
#' @param data_dir character, path to directory containing the RP100 best-fit
#'   NetCDF files (e.g. `data/all/coast_comp_data/storm_surge_muis_2023/`)
#' @param save logical, whether to save output GeoPackage. Defaults to TRUE.
#' @param output_dir character, directory path where the GeoPackage will be saved.
#'   Required when `save = TRUE`. Defaults to `here("data", "rapid", "inputs", "storm_surge")`.
#'
#' @return sf object (EPSG:4326) of CONUS storm surge stations with worst-case
#'   model RP100 absolute surge change (`max_surge_change_m`, meters), percent
#'   change relative to the historical baseline (`max_surge_pct_change`), and
#'   the model name (`max_surge_model`)
#'
process_storm_surge_data <- function(data_dir,
                                     save       = TRUE,
                                     output_dir = here("data", "rapid", "inputs", "storm_surge")) {

  # Input validation
  stopifnot(
    "data_dir does not exist"                   = dir.exists(data_dir),
    "save is TRUE but output_dir not provided"  = !save || !is.null(output_dir)
  )

  # List best-fit RP100 files only
  all_files <- list.files(data_dir, pattern = "rp100_best-fit.*\\.nc$", full.names = TRUE)

  if (length(all_files) == 0) {
    stop("No rp100_best-fit NetCDF files found in: ", data_dir)
  }

  # Parse model names from filenames: everything before _historical or _future
  model_names <- unique(
    sub("_(historical|future)_.*", "", basename(all_files))
  )

  # For each model, read both periods and compute surge change
  model_results <- map(model_names, function(model) {

    hist_file <- grep(paste0(model, "_historical.*1985-2014"), all_files, value = TRUE)
    fut_file  <- grep(paste0(model, "_future"),                all_files, value = TRUE)

    if (length(hist_file) != 1 || length(fut_file) != 1) {
      warning("Expected exactly one historical (1985-2014) and one future file for model: ", model,
              " — skipping.")
      return(NULL)
    }

    # Read historical
    nc_hist    <- ncdf4::nc_open(hist_file)
    lat_vec    <- ncdf4::ncvar_get(nc_hist, "station_y_coordinate")
    lon_vec    <- ncdf4::ncvar_get(nc_hist, "station_x_coordinate")
    surge_hist <- ncdf4::ncvar_get(nc_hist, "return_mean_surge_level")
    ncdf4::nc_close(nc_hist)

    # Read future
    nc_fut    <- ncdf4::nc_open(fut_file)
    surge_fut <- ncdf4::ncvar_get(nc_fut, "return_mean_surge_level")
    ncdf4::nc_close(nc_fut)

    tibble(
      station_idx  = seq_along(lat_vec),
      lon          = lon_vec,
      lat          = lat_vec,
      model        = model,
      surge_hist   = surge_hist,
      surge_change = surge_fut - surge_hist
    )
  })

  # Drop any NULLs (skipped models)
  model_results <- compact(model_results)

  if (length(model_results) == 0) stop("No valid model pairs found.")

  all_stations <- bind_rows(model_results)

  # Filter to CONUS bounding box: lat 24-50 N, lon -125 to -65 W
  conus_stations <- all_stations %>%
    filter(lat >= 24, lat <= 50,
           lon >= -125, lon <= -65)

  # For each station, keep the model with the highest surge change
  worst_case <- conus_stations %>%
    group_by(station_idx, lon, lat) %>%
    slice_max(surge_change, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(max_surge_pct_change = (surge_change / surge_hist) * 100) %>%
    rename(
      max_surge_change_m = surge_change,
      max_surge_model    = model
    ) %>%
    select(-surge_hist)

  # Convert to sf
  stations_sf <- worst_case %>%
    select(lon, lat, max_surge_change_m, max_surge_pct_change, max_surge_model) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    out_path <- file.path(output_dir, "storm_surge_exposure_conus.gpkg")
    sf::st_write(stations_sf, out_path, delete_dsn = TRUE, quiet = TRUE)
  }

  return(stations_sf)
}
