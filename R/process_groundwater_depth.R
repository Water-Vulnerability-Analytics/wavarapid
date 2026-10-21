#' Process Modeled Distance to Groundwater (Water Table) Depth for Water Supplies
#'
#' Extracts modeled water table depth (WTD) at each water supply point from
#' the CONUS-wide random forest WTD model. For those points that do not overlap 
#' with the raster dataset and therefore are NA we assign them the closest 
#' available raster value within a 50km buffer.
#' 
#' The source raster is stored as full-width LZW-compressed scanline strips
#' not tiles — so points are sorted by raster row before extraction to avoid 
#' re-decompressing the same strip out of sequence, then restored to their 
#' original order.
#' 
#' Saves one output:
#'   1. `wtd_all_sites.csv` — wsd_source_id, park_unit, park_name, wtd_point_m
#'
#' Point extraction (rather than a buffer mean) was validated against WSD
#' sites with recorded static water levels in
#' `Explore/explore_groundwater_depth.Rmd` and found to best match the
#' published model's own performance metrics.
#'
#' Citation: Zenodo record 18504963,
#' https://ngwa.onlinelibrary.wiley.com/doi/10.1111/gwat.13362
#'
#' @param water_supplies sf object with wsd_source_id, park_unit, park_name.
#' @param wtd_path character, path to the CONUS WTD raster. Defaults to the
#'   local copy downloaded from Zenodo; can also be a `/vsicurl/` URL.
#' @param start_buffer_m numeric, initial search radius (m) used to fill
#'   sites with no raster overlap. Defaults to 500.
#' @param max_buffer_m numeric, largest search radius (m) tried before giving
#'   up on a site. Defaults to 50000.
#' @param save logical, whether to save output CSV. Defaults to TRUE.
#' @param output_dir character, directory where CSV is written.
#'   Defaults to `wava_data_path("rapid", "inputs", "groundwater_depth")`.
#' @param verbose logical, print progress messages. Defaults to TRUE.
#'
#' @return data.frame of wsd_source_id, park_unit, park_name, wtd_point_m.

#' @export
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
