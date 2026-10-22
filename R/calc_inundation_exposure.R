#' @title Sea level rise exposure: change in inundated area
#'
#' @description
#' Measures how much of each coastal water supply's 4 km buffer is flooded
#' under the NOAA Sea Level Rise Viewer surfaces, now (0 ft) and under the
#' projected 2050 sea level rise picked by [process_inundation_data()]. For
#' each buffer it counts the flooded raster cells in both scenarios and
#' calculates the change in cells and in percentage points, and whether the
#' water supply point itself becomes inundated. Buffers within one raster
#' tile are processed together with exactextractr; buffers that span
#' several tiles are mosaicked and processed one at a time.
#'
#' Supplies with no scenario (no `proj_file`) and supplies more than 50 km
#' from the coast (`coast_dist_m` is `NA` in `coast_distances`) get `NA` for
#' every indicator. They are set to `NA` before standardizing, so they don't
#' affect the scores of the coastal supplies.
#'
#' @param ws_scenarios sf object. Every water supply, with the buffered
#'   scenarios from [process_inundation_data()] joined on: coastal rows have
#'   the buffer polygon and the scenario columns, other rows have empty
#'   geometry and `NA`s (see the example). Needs `wsd_source_id`,
#'   `park_unit`, `buffer_km`, `inun_proj_mm`, `inun_region`,
#'   `inun_proj_closest_mm`, `depth_ft`, `depth_m`, `inun_surf_url`,
#'   `proj_file`, `ref_file` and a geometry column named `geom` (as read
#'   from `buff_ws_inun_regions.gpkg`). Required.
#' @param inundation_raster_data_dir String. The folder with the downloaded
#'   rasters, in one subfolder per state (`<state>/<file>`). Default:
#'   `rapid/inputs/coast_comp/Extent_Rasters` in the shared data folder,
#'   which is where [process_inundation_data()] downloads them, see
#'   [wavacore::wava_data_path()].
#' @param coast_distances String. Path to `coastal_site_geometry.csv` from
#'   [process_coastal_geometry()]; only `coast_dist_m` is used. Stops with an
#'   error if the file doesn't exist. Default:
#'   `rapid/inputs/coast_comp/coastal_site_geometry.csv` in the shared data
#'   folder.
#' @param save Logical. Whether to write the results to a CSV file.
#'   * `TRUE`: writes `inundation_exposure_<date>.csv`, named with today's date, to
#'     `output_dir`, creating the folder if needed. Stops with an error if
#'     `output_dir` is `NULL`.
#'   * `FALSE`: writes nothing.
#'
#'   Default: `FALSE`.
#' @param output_dir String or `NULL`. The folder for the CSV file.
#'   * A folder path: where the file is written when `save = TRUE`. Ignored
#'     when `save = FALSE`.
#'   * `NULL`: only allowed with `save = FALSE`.
#'
#'   Default: `NULL`.
#' @param verbose Logical. Whether to print progress.
#'   * `TRUE`: shows a progress bar over the raster files and messages with
#'     the run time and the path of the saved file.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return A tibble, one row per water supply, with `wsd_source_id`,
#'   `park_unit`, `buffer_km`, `inun_proj_mm`, `inun_proj_closest_mm`,
#'   `depth_ft` and:
#'   * `raw_coverage_pct_ref` and `raw_coverage_pct_proj`: the percentage of
#'     the buffer flooded now and under the projection;
#'   * `raw_cell_change`: the change in flooded cells (projection minus now);
#'   * `raw_pct_point_change`: the change in percentage points;
#'   * `norm_*` and `rank_*` versions of those four, from
#'     [add_standardized_indicators()];
#'   * `inun_change`: `"no_to_no"`, `"no_to_some"`, `"some_to_some"` or
#'     `"some_to_no"`, for the buffer;
#'   * `inun_centroid_change`: `"no_to_no"`, `"no_to_yes"`, `"yes_to_yes"` or
#'     `"yes_to_no"`, for the water supply point.
#' @export
#' @seealso [process_inundation_data()], which makes the scenarios and
#'   downloads the rasters; [process_coastal_geometry()];
#'   [calc_storm_surge_exposure()] and [calc_swi_exposure()], the other
#'   coastal indicators.
#' @examplesIf interactive() && file.exists(wavacore::wava_data_path("rapid", "inputs", "coast_comp", "buff_ws_inun_regions.gpkg"))
#' # The scenarios, joined onto the supplies of one coastal park, as in the
#' # rapid workflow
#' ws <- wavacore::get_wsd(park = "BISC")
#' scenarios <- sf::st_read(wavacore::wava_data_path("rapid", "inputs", "coast_comp",
#'                                                   "buff_ws_inun_regions.gpkg"),
#'                          quiet = TRUE)
#' scenarios <- scenarios[!is.na(scenarios$proj_file), names(scenarios) != "park_unit"]
#' ws_ids <- sf::st_drop_geometry(ws)[, c("wsd_source_id", "park_unit")]
#' ws_scenarios <- dplyr::right_join(scenarios, ws_ids, by = "wsd_source_id")
#'
#' slr_exp <- calc_inundation_exposure(ws_scenarios)
#' slr_exp[, c("wsd_source_id", "raw_coverage_pct_ref", "raw_coverage_pct_proj",
#'             "inun_change", "inun_centroid_change")]

calc_inundation_exposure <- function(ws_scenarios,
                                     inundation_raster_data_dir = wavacore::wava_data_path("rapid", "inputs", "coast_comp", "Extent_Rasters"),
                                     coast_distances = wavacore::wava_data_path("rapid", "inputs", "coast_comp", "coastal_site_geometry.csv"),
                                     save = FALSE,
                                     output_dir = NULL,
                                     verbose = TRUE) {
  # Input validation
  required_columns <- c("inun_region", "ref_file", "proj_file")

  stopifnot("ws_scenarios must be an sf object" = inherits(ws_scenarios, "sf"),
            "save is TRUE but output_dir not provided" = !save || !is.null(output_dir),
            "ws_scenarios missing required columns" = all(required_columns %in% names(ws_scenarios)),
            "coast distances file does not exist" = file.exists(coast_distances))
  
  if (verbose) start_time <- Sys.time()

  # Separate non-analysis sites (no SLR scenario info) — recombined at the end
  ws_not_in_analysis <- ws_scenarios %>% dplyr::filter(is.na(proj_file))
  ws_in_analysis     <- ws_scenarios %>% dplyr::filter(!is.na(proj_file))

  ws_scenarios_mod <- ws_in_analysis %>%
    dplyr::ungroup() %>% 
    dplyr::distinct() %>%
    dplyr::mutate(identity = paste0(wsd_source_id, "-", buffer_km)) %>%
    dplyr::select(c(park_unit, identity, wsd_source_id, buffer_km, inun_proj_mm,
             inun_region, inun_proj_closest_mm, depth_ft, depth_m, inun_surf_url,
             proj_file, ref_file, geom))
  
  region_list <- ws_scenarios_mod %>% 
    split(f = .$proj_file) 
  
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "Processing inundation exposure [:bar] :percent (:current/:total) ETA: :eta",
      total  = length(region_list),
      clear  = FALSE,
      width  = 80
    )
  }
 
  # Main loop over regions ====
  ws_inun_stats <- purrr::map_dfr(region_list, function(region_sf) {
    
    if (verbose) pb$tick()
    
    # tiles needed for this region 
    locs_needed <- ws_scenarios_mod %>%
      dplyr::filter(identity %in% region_sf$identity) %>%
      dplyr::distinct(proj_file, ref_file) %>%
      dplyr::mutate(state = substr(proj_file, 1, 2))
    
    # identify IDs with overlapping inundation layers
    overlap_ids <-  ws_scenarios_mod %>%
      dplyr::filter(identity %in% region_sf$identity) %>%
      dplyr::group_by(identity) %>%
      dplyr::mutate(count = dplyr::n()) %>%
      dplyr::filter(count > 1) %>%
      dplyr::pull(identity) %>%
      unique()

    overlap_files <- ws_scenarios_mod %>%
      dplyr::filter(identity %in% region_sf$identity) %>%
      dplyr::group_by(identity) %>%
      dplyr::mutate(count = dplyr::n()) %>%
      dplyr::filter(count > 1)

    # build full paths to all tiles (using per-row state to avoid recycling mismatches)
    ref_tiles  <- sort(unique(file.path(inundation_raster_data_dir, locs_needed$state, locs_needed$ref_file)))
    proj_tiles <- sort(unique(file.path(inundation_raster_data_dir, locs_needed$state, locs_needed$proj_file)))
    
    stopifnot(length(ref_tiles)  >= 1, all(file.exists(ref_tiles)))
    stopifnot(length(proj_tiles) >= 1, all(file.exists(proj_tiles)))
    
    # split buffer indices into batch (all sites contained within single inundation
    # layer) vs. one-by-one (those with overlapping inundation layers)
    idx_batch <- region_sf %>% 
      dplyr::filter(!identity %in% overlap_ids)
    
    idx_one_by_one <- ws_scenarios_mod %>% 
      dplyr::filter(identity %in% overlap_ids)
    
    # BATCH: non-overlapping sites on a single tile.
    # Use {exactextractr} for super fast polygon-raster extraction
    if (nrow(idx_batch) > 0) {
      
      batch_sf <- sf::st_as_sf(idx_batch)
      
      counts_list <- list()
      centroid_list <- list()

      for (scenario in c("ref", "proj")) {

        file_col  <- if (scenario == "ref") "ref_file" else "proj_file"

        tile_file <- unique(region_sf[[file_col]])

        stopifnot("batch path expects exactly one tile per scenario \u2014 proj_file split should guarantee this" = length(tile_file) == 1)

        tile_path <- file.path(inundation_raster_data_dir, substr(tile_file, 1, 2), tile_file)

        stopifnot(file.exists(tile_path))

        rast <- terra::rast(tile_path)

        # exactextractr::exact_extract returns one row per raster cell that
        # intersects the polygon, regardless of whether the cell value is NA.
        # This is correct for NOAA inundation rasters where NA = non-inundated:
        # nrow(df) gives the true total cell count (denominator) and
        # sum(!is.na(df$value)) counts only flooded cells (numerator).
        # The batch path therefore handles all-NA (non-inundated) tiles correctly
        # without any special-casing — unlike the one-by-one mosaic path (see
        # prep_tile above), which previously discarded all-NA tiles and lost
        # those cells from the total denominator.
        # Explicitly reproject batch_sf to match the raster.
        extracted <- suppressWarnings(
          exactextractr::exact_extract(rast, sf::st_transform(batch_sf, terra::crs(rast)), progress = FALSE)
        )

        counts_list[[scenario]] <- data.frame(identity = batch_sf$identity,
                                              buffer_km = batch_sf$buffer_km,
                                              total = vapply(extracted, \(df) nrow(df), numeric(1)),
                                              flooded = vapply(extracted, \(df) sum(!is.na(df$value)), numeric(1)),
                                              scenario = scenario)

        # Extract at centroids using already-loaded raster (terra::extract used here
        # since exactextractr does not support point geometries)
        batch_centroids <- terra::vect(sf::st_transform(suppressWarnings(sf::st_centroid(batch_sf)), terra::crs(rast)))
        centroid_list[[scenario]] <- !is.na(terra::extract(rast, batch_centroids)[, 2])
      }

      batch_counts <- dplyr::bind_rows(counts_list) %>%
        tidyr::pivot_wider(names_from  = scenario,
                    values_from = c(total, flooded),
                    names_glue  = "raw_{.value}_cells_{scenario}") %>%
        dplyr::mutate(inun_centroid_change = dplyr::case_when(
          !centroid_list[["ref"]] & !centroid_list[["proj"]] ~ "no_to_no",
          !centroid_list[["ref"]] &  centroid_list[["proj"]] ~ "no_to_yes",
           centroid_list[["ref"]] &  centroid_list[["proj"]] ~ "yes_to_yes",
           centroid_list[["ref"]] & !centroid_list[["proj"]] ~ "yes_to_no",
          TRUE ~ NA_character_
        ))
    } else {
    batch_counts <- dplyr::tibble()
    }
    
    # ONE-BY-ONE: overlapping sites may span multiple tiles,
    # so build a mosaic per buffer then use exactextractr
    if (nrow(idx_one_by_one) > 0) {
      
      one_by_one_list <- list()
      centroid_obo_list <- list()

      combos <- unique(idx_one_by_one$identity)

      for (i in seq_len(length(combos))) {

        b <- idx_one_by_one %>%
          dplyr::filter(identity == combos[i]) %>%
          dplyr::distinct(identity, .keep_all = TRUE)

        b_sf <- sf::st_as_sf(b)
        centroid_vals_scenarios <- list(ref = FALSE, proj = FALSE)

        for (scenario in c("ref", "proj")) {

          tiles <- if (scenario == "ref") ref_tiles else proj_tiles

          tiles_words <- overlap_files %>%
            dplyr::filter(identity == b_sf$identity) %>%
            dplyr::pull(if (scenario == "ref") ref_file else proj_file)

          tiles <- tiles[basename(tiles) %in% tiles_words]

          if (length(tiles) == 0) {
            warning("No matching tiles for rowid ", b_sf$identity, ", scenario ", scenario)
            next
          }

          out <- NULL

          # fxn to crop and mask a single raster tile to a buffer polygon.
          # NOAA inundation rasters use NA to represent non-inundated cells, so an
          # all-NA tile is valid: it means the buffer overlaps the tile extent but
          # no flooding occurs there. We must NOT return NULL for all-NA tiles —
          # doing so would exclude those cells from the total count, inflating the
          # coverage percentage (or producing NA totals) for multi-tile water supplies.
          prep_tile <- function(r_path, b_sf) {

            r0 <- terra::rast(r_path)
            b_projected <- sf::st_transform(b_sf, terra::crs(r0))

            r1 <- tryCatch(terra::crop(r0, b_projected), error = function(e) NULL)
            if (is.null(r1) || !inherits(r1, "SpatRaster")) return(NULL)
            if (terra::ncell(r1) == 0) return(NULL)  # truly outside tile extent

            r2 <- tryCatch(terra::mask(r1, terra::vect(b_projected)), error = function(e) NULL)
            if (is.null(r2) || !inherits(r2, "SpatRaster")) return(NULL)
            # Do NOT filter all-NA r2: an all-NA masked tile means zero flooding
            # in this tile, but its cells must still count toward the total denominator.

            return(r2)

          }

         tile_list <- list()

          for (r in tiles) {
            t <- tryCatch(prep_tile(r, b_sf), error = function(e) NULL)
            if (!is.null(t) && inherits(t, "SpatRaster")) tile_list[[length(tile_list) + 1]] <- t
          }
          if (length(tile_list) > 0) {
            # Normalize resolutions before mosaicking: NOAA tiles from adjacent
            # regions can have different pixel sizes (across-region differences
            # of up to ~6e-6 degrees, plus floating-point jitter of ~1e-7 within
            # the same nominal family) that cause terra::mosaic to error with
            # "resolution does not match". Resample all tiles to the coarsest
            # resolution among the set — this avoids fabricating sub-pixel
            # detail that doesn't exist in the coarser source. Nearest-neighbour
            # is used to preserve NA = non-inundated cell semantics.
            if (length(tile_list) > 1) {
              coarsest_idx <- which.max(sapply(tile_list, function(t) terra::res(t)[1]))
              template     <- tile_list[[coarsest_idx]]
              tile_list    <- lapply(tile_list, function(t) {
                if (isTRUE(all.equal(terra::res(t), terra::res(template)))) return(t)
                terra::resample(t, template, method = "near")
              })
            }
            out <- tryCatch(
              terra::mosaic(terra::sprc(tile_list), fun = "max"),
              error = function(e) {
                message("mosaic failed: ", e$message)
                NULL
              }
            )
          }

          if (is.null(out)) {

            one_by_one_list[[paste(b_sf$identity, scenario, sep = "_")]] <- data.frame(identity = b_sf$identity,
                                                                                       buffer_km = b_sf$buffer_km,
                                                                                       total = NA,
                                                                                       flooded = 0,
                                                                                       scenario = scenario)
            centroid_vals_scenarios[[scenario]] <- FALSE

          } else {

            extracted <- suppressWarnings(
              exactextractr::exact_extract(out, sf::st_transform(b_sf, terra::crs(out)), progress = FALSE)
            )[[1]]

            one_by_one_list[[paste(b_sf$identity, scenario, sep = "_")]] <- data.frame(identity = b_sf$identity,
                                                                                       buffer_km = b_sf$buffer_km,
                                                                                       total = nrow(extracted),
                                                                                       flooded = sum(!is.na(extracted$value)),
                                                                                       scenario = scenario)

            # Extract at centroid using already-built mosaic (terra::extract used here
            # since exactextractr does not support point geometries)
            b_centroid_vect <- terra::vect(sf::st_transform(suppressWarnings(sf::st_centroid(b_sf)), terra::crs(out)))
            centroid_vals_scenarios[[scenario]] <- !is.na(terra::extract(out, b_centroid_vect)[, 2])
          }
        }

        centroid_obo_list[[combos[i]]] <- data.frame(
          identity  = combos[i],
          inun_centroid_change = dplyr::case_when(
            !centroid_vals_scenarios[["ref"]] & !centroid_vals_scenarios[["proj"]] ~ "no_to_no",
            !centroid_vals_scenarios[["ref"]] &  centroid_vals_scenarios[["proj"]] ~ "no_to_yes",
             centroid_vals_scenarios[["ref"]] &  centroid_vals_scenarios[["proj"]] ~ "yes_to_yes",
             centroid_vals_scenarios[["ref"]] & !centroid_vals_scenarios[["proj"]] ~ "yes_to_no",
            TRUE ~ NA_character_
          )
        )
      }

      one_by_one_counts <- dplyr::bind_rows(one_by_one_list) %>%
        tidyr::pivot_wider(names_from = scenario,
                    values_from = c(total, flooded),
                    names_glue = "raw_{.value}_cells_{scenario}") %>%
        dplyr::mutate(raw_total_cells_ref = dplyr::if_else(is.na(raw_total_cells_ref), raw_total_cells_proj, raw_total_cells_ref)) %>%
        dplyr::left_join(dplyr::bind_rows(centroid_obo_list), by = "identity")
    } else {
      one_by_one_counts <- dplyr::tibble()
    }
    
    region_results <- dplyr::bind_rows(one_by_one_counts, batch_counts) %>%
      dplyr::left_join(region_sf, by = c("buffer_km", "identity")) %>%
      dplyr::mutate(raw_coverage_pct_ref = dplyr::if_else(raw_total_cells_ref > 0,
                                            (raw_flooded_cells_ref / raw_total_cells_ref) * 100,
                                            NA_real_),
             raw_coverage_pct_proj = dplyr::if_else(raw_total_cells_proj > 0,
                                             (raw_flooded_cells_proj / raw_total_cells_proj) * 100,
                                             NA_real_),
             raw_cell_change = raw_flooded_cells_proj - raw_flooded_cells_ref,
             raw_pct_point_change = raw_coverage_pct_proj - raw_coverage_pct_ref,
             inun_change = dplyr::case_when(raw_coverage_pct_ref == 0 & raw_coverage_pct_proj == 0 ~ "no_to_no",
                                     raw_coverage_pct_ref == 0 & raw_coverage_pct_proj > 0  ~ "no_to_some",
                                     raw_coverage_pct_ref > 0  & raw_coverage_pct_proj == 0 ~ "some_to_no",
                                     raw_coverage_pct_ref > 0  & raw_coverage_pct_proj > 0  ~ "some_to_some",
                                     TRUE ~ NA_character_)) %>%
      dplyr::select(identity, wsd_source_id, buffer_km,
             raw_total_cells_ref,    # Total raster cells intersecting polygon in reference scenario
             raw_flooded_cells_ref,  # Count of flooded cells in reference scenario
             raw_coverage_pct_ref,   # Reference flood coverage as proportion of intersecting area
             raw_total_cells_proj,   # Total raster cells intersecting polygon in projection scenario
             raw_flooded_cells_proj, # Count of flooded cells in projection scenario
             raw_coverage_pct_proj,  # Projection flood coverage as proportion of intersecting area
             raw_cell_change,        # Net change in flooded cell count (proj - ref)
             raw_pct_point_change,   # Percentage point difference in coverage (proj - ref)
             inun_change,            # Categorical indicator for inundation change patterns
             inun_centroid_change)   # Whether the water supply point itself becomes inundated
    
    return(region_results)
    
  })
  
  # Clean up final results

  # SLR metadata per water supply (distinct handles multi-region sites)
  ws_geom_meta <- ws_in_analysis %>%
    sf::st_drop_geometry() %>%
    dplyr::select(wsd_source_id, inun_proj_mm, inun_proj_closest_mm, depth_ft) %>%
    dplyr::distinct()

  # Processed coastal results with metadata
  coastal_results <- ws_inun_stats %>%
    dplyr::left_join(ws_geom_meta, by = "wsd_source_id") %>%
    dplyr::select(wsd_source_id, buffer_km,
           inun_proj_mm, inun_proj_closest_mm, depth_ft,
           raw_total_cells_ref, raw_flooded_cells_ref, raw_coverage_pct_ref,
           raw_total_cells_proj, raw_flooded_cells_proj, raw_coverage_pct_proj,
           raw_cell_change, raw_pct_point_change, inun_change, inun_centroid_change) %>%
    dplyr::mutate(dplyr::across(c(raw_total_cells_ref, raw_flooded_cells_ref, raw_coverage_pct_ref,
                    raw_total_cells_proj, raw_flooded_cells_proj, raw_coverage_pct_proj,
                    raw_cell_change, raw_pct_point_change),
                  ~ ifelse(inun_change == "no_to_no", 0, .x))) %>%
    dplyr::distinct()

  # All unique water supplies (base for expansion)
  all_ws_base <- ws_scenarios %>%
    sf::st_drop_geometry() %>%
    dplyr::select(wsd_source_id, park_unit) %>%
    dplyr::distinct()

  # Water supplies >50 km from the coast get every indicator column set to NA.
  # coast_dist_m is already NA'd upstream (process_coastal_geometry()) beyond
  # 50 km, so we only need to check is.na(coast_dist_m) here — no new distance
  # math. This screening must happen BEFORE add_standardized_indicators(),
  # matching the convention in calc_swi_exposure()/calc_storm_surge_exposure():
  # proj_file/ref_file assignment (inundation-region bbox intersection) is a
  # different spatial test than the 50 km CUSP shoreline distance, so
  # coastal_results can contain sites that are >50 km by coast_dist_m. Masking
  # their raw values to NA first (rather than after standardizing) keeps them
  # out of the range()/cume_dist() pool that norm_*/rank_* are computed from
  # for genuinely in-range sites. Metadata columns (buffer_km, inun_proj_mm,
  # inun_proj_closest_mm, depth_ft) are left untouched; only the risk/indicator
  # columns are screened.
  coast_distances_df <- readr::read_csv(coast_distances, show_col_types = FALSE) %>%
    dplyr::select(wsd_source_id, coast_dist_m) %>%
    dplyr::distinct()

  indicator_cols <- c(
    "raw_coverage_pct_ref", "raw_coverage_pct_proj",
    "raw_cell_change", "raw_pct_point_change",
    "inun_change", "inun_centroid_change"
  )

  # Expand to all water supplies; non-coastal rows get NA for all indicator columns.
  # group_by(buffer_km): non-coastal rows have buffer_km = NA → their own group;
  final_results <- all_ws_base %>%
    dplyr::left_join(coastal_results, by = "wsd_source_id") %>%
    dplyr::left_join(coast_distances_df, by = "wsd_source_id") %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(indicator_cols), ~ replace(., is.na(coast_dist_m), NA))) %>%
    dplyr::select(-coast_dist_m) %>%
    dplyr::group_by(buffer_km) %>%
    add_standardized_indicators("raw_coverage_pct_ref") %>%
    add_standardized_indicators("raw_coverage_pct_proj") %>%
    add_standardized_indicators("raw_cell_change") %>%
    add_standardized_indicators("raw_pct_point_change") %>%
    dplyr::ungroup() %>%
    dplyr::select(wsd_source_id, park_unit, buffer_km,
           inun_proj_mm, inun_proj_closest_mm, depth_ft,
           raw_coverage_pct_ref,  norm_coverage_pct_ref,  rank_coverage_pct_ref,
           raw_coverage_pct_proj, norm_coverage_pct_proj, rank_coverage_pct_proj,
           raw_cell_change,       norm_cell_change,       rank_cell_change,
           raw_pct_point_change,  norm_pct_point_change,  rank_pct_point_change,
           inun_change, inun_centroid_change)

  # Save if requested
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    filename <- file.path(output_dir, paste0("inundation_exposure_", Sys.Date(), ".csv"))
    readr::write_csv(final_results, filename)
    if (verbose) message("Results saved to: ", filename)
  }
  
  if (verbose) {
    end_time   <- Sys.time()
    total_time <- end_time - start_time
    message("\nInundation exposure processing completed in ",
            round(total_time, 2), " ", attr(total_time, "units"))
  }
  
  return(final_results)
  
}
