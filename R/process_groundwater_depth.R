#' @title Modeled depth to the water table at each water supply
#'
#' @description
#' Reads the modeled water table depth (WTD) at each water supply point from
#' the CONUS random forest WTD model. Points with no value (a NoData cell, or
#' just outside the raster) get the value of the nearest valid cell. The
#' search starts at `start_buffer_m` and doubles until it finds one, up to
#' `max_buffer_m`; points with no valid cell within that stay `NA`. The search
#' runs in parallel on all but two CPU cores (with
#' [future::multisession]), and the future plan is set back to sequential
#' afterwards.
#'
#' The points are read in raster-row order, because the raster is stored in
#' full-width strips rather than tiles; the results are put back in the
#' original order.
#'
#' Reading the value at the point, rather than a mean over a buffer, matched
#' WSD sites with recorded static water levels best (see
#' `Explore/explore_groundwater_depth.Rmd` in `nps_water_vulnerability`).
#'
#' Citation: Zenodo record 18504963,
#' <https://ngwa.onlinelibrary.wiley.com/doi/10.1111/gwat.13362>.
#'
#' @param water_supplies sf object. The water supplies, as points with
#'   `wsd_source_id`, `park_unit` and `park_name` columns, e.g. from
#'   [wavacore::get_wsd()]. Required.
#' @param wtd_path String. Path to the CONUS WTD raster (meters). Can also be
#'   a `/vsicurl/` URL, so GDAL reads only the parts it needs. Default: the
#'   copy downloaded from Zenodo, in `all/groundwater_sensitivity` in the
#'   shared data folder (about 37 GB), see [wavacore::wava_data_path()].
#' @param start_buffer_m Numeric. The first search radius, in meters, for
#'   points with no value. Default: `500`.
#' @param max_buffer_m Numeric. The largest search radius, in meters. Points
#'   with no valid cell within it get `NA`. Default: `50000` (50 km).
#' @param save Logical. Whether to write the results to a CSV file.
#'   * `TRUE`: writes `wtd_all_sites.csv` to `output_dir`, creating the folder
#'     if needed and replacing an existing file. Stops with an error if
#'     `output_dir` is `NULL`.
#'   * `FALSE`: writes nothing.
#'
#'   Default: `FALSE`.
#' @param output_dir String or `NULL`. The folder for `wtd_all_sites.csv`.
#'   * A folder path: where the file is written when `save = TRUE`. Ignored
#'     when `save = FALSE`.
#'   * `NULL`: only allowed with `save = FALSE`.
#'
#'   Default: `rapid/inputs/groundwater_depth` in the shared data folder.
#' @param verbose Logical. Whether to print progress.
#'   * `TRUE`: prints how many points are being filled, each point left as
#'     `NA`, and the run time.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return A data frame, one row per water supply, with columns
#'   `wsd_source_id`, `park_unit`, `park_name` and `wtd_point_m` (depth to
#'   the water table in meters).
#' @export
#' @seealso [calc_groundwater_depth_sensitivity()], which uses
#'   `wtd_all_sites.csv` where a supply has no recorded well depth.
#' @examplesIf interactive() && file.exists(wavacore::wava_data_path("all", "groundwater_sensitivity", "wtd_mean_estimate_RF_additional_inputs_dummy_drop0LP_1s_CONUS2_m_v_20240813.tif"))
#' # The default raster is about 37 GB, so the first read takes a long time
#' ws <- wavacore::get_wsd(park = "ZION")
#' wtd <- process_groundwater_depth(ws)
#' wtd

process_groundwater_depth <- function(
    water_supplies,
    wtd_path       = wavacore::wava_data_path("all", "groundwater_sensitivity",
                           "wtd_mean_estimate_RF_additional_inputs_dummy_drop0LP_1s_CONUS2_m_v_20240813.tif"),
    start_buffer_m = 500,
    max_buffer_m   = 50000,
    save           = FALSE,
    output_dir     = wavacore::wava_data_path("rapid", "inputs", "groundwater_depth"),
    verbose        = TRUE) {

  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "water_supplies must include wsd_source_id" = "wsd_source_id" %in% names(water_supplies),
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir)
  )
  
  if (verbose) {
    message("Starting groundwater depth processing...")
    start_time <- Sys.time()
  }
  
  wtd <- terra::rast(wtd_path)
  
  points_native <- water_supplies %>%
    sf::st_transform(terra::crs(wtd, proj = TRUE))
  
  pts_vect <- terra::vect(points_native)
  
  # Raster is stored as full-width scanline strips, not tiles, so extracting
  # in row order means each strip is read/decompressed once instead of being
  # revisited out of sequence.
  row_idx <- terra::rowFromY(wtd, terra::geom(pts_vect)[, "y"])
  row_order <- order(row_idx)
  
  extracted <- terra::extract(wtd, pts_vect[row_order])[, 2]
  wtd_point_m <- numeric(length(row_order))
  wtd_point_m[row_order] <- extracted
  
  # Fill sites with no raster overlap (NoData cell or just outside the
  # extent) using the nearest valid cell within a growing search buffer
  na_idx <- which(is.na(wtd_point_m))
  
  if (length(na_idx) > 0) {
    
    if (verbose) {
      message(sprintf("Filling %d site(s) with no direct raster overlap using nearest valid cell...",
                      length(na_idx)))
    }
    
    # sf/SpatVector objects hold external C++ pointers that don't survive
    # serialization to future workers, so pass plain coordinates + crs WKT
    # instead of subsetting an sf object inside the parallel closure
    xy_native <- terra::geom(pts_vect)[, c("x", "y")]
    wtd_crs <- terra::crs(wtd)
    
    fill_nearest_valid <- function(i) {
      
      # Each parallel worker is a separate process, so the SpatRaster's
      # underlying pointer has to be re-opened here rather than shared
      wtd_worker <- terra::rast(wtd_path)
      
      pt_xy <- xy_native[i, ]
      pt_vect <- terra::vect(matrix(pt_xy, nrow = 1), type = "points", crs = wtd_crs)
      buffer_m <- start_buffer_m

      repeat {

        # crop() alone only restricts to buffer()'s rectangular bounding
        # box, not the circular buffer itself, which would let a valid cell
        # near a box corner (up to buffer_m * sqrt(2) away) get accepted as
        # if it were within buffer_m. mask() to the circular polygon so only
        # cells genuinely within the search radius are ever considered.
        pt_buf <- terra::buffer(pt_vect, width = buffer_m)
        pt_crop <- terra::crop(wtd_worker, pt_buf)
        pt_mask <- terra::mask(pt_crop, pt_buf)
        valid_pts <- terra::as.points(pt_mask, values = TRUE, na.rm = TRUE)

        if (length(valid_pts) > 0) {
          # terra::nearest is CRS-aware (planar distance for projected rasters,
          # geodesic for geographic ones), unlike a hand-rolled xy distance calc
          nn <- terra::nearest(pt_vect, valid_pts)
          return(terra::values(valid_pts)[nn$to_id, 1])
        }

        if (buffer_m >= max_buffer_m) {
          if (verbose) {
            message(sprintf("  No valid cell found within %.0f m for site %s; leaving as NA.",
                            max_buffer_m, water_supplies$wsd_source_id[i]))
          }
          return(NA_real_)
        }

        # Clamp to max_buffer_m so the final attempt lands exactly at the
        # documented cap instead of doubling past it (e.g. 32km -> 64km)
        buffer_m <- min(buffer_m * 2, max_buffer_m)
      }
    }
    
    n_workers <- max(1, future::availableCores() - 2)
    future::plan(future::multisession, workers = n_workers)
    wtd_point_m[na_idx] <- furrr::future_map_dbl(na_idx, fill_nearest_valid)
    future::plan(future::sequential)
  }
  
  wtd_all_sites <- water_supplies %>%
    sf::st_drop_geometry() %>%
    dplyr::mutate(wtd_point_m = wtd_point_m) %>%
    dplyr::select(wsd_source_id, park_unit, park_name, wtd_point_m)
  
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    utils::write.csv(wtd_all_sites, file.path(output_dir, "wtd_all_sites.csv"), row.names = FALSE)
  }
  
  if (verbose) {
    message(sprintf("Water Table Depth processing complete in %.1f minutes.",
                    as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
  }
  
  return(wtd_all_sites)

}
