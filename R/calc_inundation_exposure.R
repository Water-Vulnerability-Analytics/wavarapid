#' Calculate Sea Level Rise Inundation Exposure Indicator
#'
#' Calculates inundation exposure for buffered water supply polygons using NOAA
#' Sea Level Rise Viewer inundation surfaces. For each polygon, the function
#' extracts flooded cell counts and percent coverage under both a reference
#' (current) and a projected SLR scenario, then computes absolute cell change,
#' percentage point change, and whether the water supply centroid itself becomes
#' inundated. Sites whose buffer overlaps a single inundation tile are processed
#' in batch via exactextractr; sites spanning multiple tiles are mosaicked and
#' processed individually.
#'
#' @param ws_scenarios sf object of ALL water supplies right-joined onto buffered
#'   coastal SLR scenarios. Coastal rows have buffer polygon geometry and SLR
#'   scenario columns (e.g. `inun_region`, `ref_file`, `proj_file`,
#'   `buffer_km`); non-coastal rows have NA geometry and NA SLR columns. Must
#'   include columns `wsd_source_id`, `park_unit`, `inun_region`,
#'   `ref_file`, and `proj_file` (NAs acceptable for non-coastal rows).
#' @param inundation_raster_data_dir character, path to inundation surface raster
#'   state subdirectories
#' @param save logical, whether to save output
#' @param output_dir character, directory path for saving files
#' @param verbose logical, whether to show progress messages
#'
#' @return tibble with water supply ID, buffer size, SLR projection metadata,
#'   raw/normalized/ranked coverage indicators for reference and projected
#'   scenarios (`raw/norm/rank_coverage_pct_ref`, `raw/norm/rank_coverage_pct_proj`),
#'   raw/normalized/ranked cell change and percentage point change
#'   (`raw/norm/rank_cell_change`, `raw/norm/rank_pct_point_change`),
#'   an `inun_change` category (`no_to_no`, `no_to_some`, `some_to_some`,
#'   `some_to_no`), and `inun_centroid_change` (categorical: `no_to_no`,
#'   `no_to_yes`, `yes_to_yes`, `yes_to_no` — whether the water supply centroid
#'   point was inundated in the reference and/or projected scenario).

calc_inundation_exposure <- function(ws_scenarios,
                                     inundation_raster_data_dir = here("data", "rapid", "inputs", "coast_comp", "Extent_Rasters"),
                                     save = FALSE,
                                     output_dir = NULL,
                                     verbose = TRUE) {
  # Input validation
  required_columns <- c("inun_region", "ref_file", "proj_file")
   
  stopifnot("ws_scenarios must be an sf object" = inherits(ws_scenarios, "sf"),
            "save is TRUE but output_dir not provided" = !save || !is.null(output_dir),
            "ws_scenarios missing required columns" = all(required_columns %in% names(ws_scenarios)))
  
  if (verbose) start_time <- Sys.time()

  # Separate non-analysis sites (no SLR scenario info) — recombined at the end
  ws_not_in_analysis <- ws_scenarios %>% filter(is.na(proj_file))
  ws_in_analysis     <- ws_scenarios %>% filter(!is.na(proj_file))

  ws_scenarios_mod <- ws_in_analysis %>%
    ungroup() %>% 
    distinct() %>%
    mutate(identity = paste0(wsd_source_id, "-", buffer_km)) %>%
    select(c(park_unit, identity, wsd_source_id, buffer_km, inun_proj_mm,
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
  ws_inun_stats <- map_dfr(region_list, function(region_sf) {
    
    if (verbose) pb$tick()
    
    # tiles needed for this region 
    locs_needed <- ws_scenarios_mod %>%
      filter(identity %in% region_sf$identity) %>%
      distinct(proj_file, ref_file) %>%
      mutate(state = substr(proj_file, 1, 2))
    
    # identify IDs with overlapping inundation layers
    overlap_ids <-  ws_scenarios_mod %>%
      filter(identity %in% region_sf$identity) %>%
      group_by(identity) %>%
      mutate(count = n()) %>%
      filter(count > 1) %>%
      pull(identity) %>%
      unique()

    overlap_files <- ws_scenarios_mod %>%
      filter(identity %in% region_sf$identity) %>%
      group_by(identity) %>%
      mutate(count = n()) %>%
      filter(count > 1)

    # build full paths to all tiles (using per-row state to avoid recycling mismatches)
    ref_tiles  <- sort(unique(here(inundation_raster_data_dir, locs_needed$state, locs_needed$ref_file)))
    proj_tiles <- sort(unique(here(inundation_raster_data_dir, locs_needed$state, locs_needed$proj_file)))
    
    stopifnot(length(ref_tiles)  >= 1, all(file.exists(ref_tiles)))
    stopifnot(length(proj_tiles) >= 1, all(file.exists(proj_tiles)))
    
    # split buffer indices into batch (all sites contained within single inundation
    # layer) vs. one-by-one (those with overlapping inundation layers)
    idx_batch <- region_sf %>% 
      filter(!identity %in% overlap_ids)
    
    idx_one_by_one <- ws_scenarios_mod %>% 
      filter(identity %in% overlap_ids)
    
    # BATCH: non-overlapping sites on a single tile.
    # Use {exactextractr} for super fast polygon-raster extraction
    if (nrow(idx_batch) > 0) {
      
      batch_sf <- sf::st_as_sf(idx_batch)
      
      counts_list <- list()
      centroid_list <- list()

      for (scenario in c("ref", "proj")) {

        file_col  <- if (scenario == "ref") "ref_file" else "proj_file"

        tile_file <- unique(region_sf[[file_col]])

        stopifnot(length(tile_file) == 1)

        tile_path <- here(inundation_raster_data_dir, substr(tile_file, 1, 2), tile_file)

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
        extracted <- suppressWarnings(
          exactextractr::exact_extract(rast, batch_sf, progress = FALSE)
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

      batch_counts <- bind_rows(counts_list) %>%
        pivot_wider(names_from  = scenario,
                    values_from = c(total, flooded),
                    names_glue  = "raw_{.value}_cells_{scenario}") %>%
        mutate(inun_centroid_change = case_when(
          !centroid_list[["ref"]] & !centroid_list[["proj"]] ~ "no_to_no",
          !centroid_list[["ref"]] &  centroid_list[["proj"]] ~ "no_to_yes",
           centroid_list[["ref"]] &  centroid_list[["proj"]] ~ "yes_to_yes",
           centroid_list[["ref"]] & !centroid_list[["proj"]] ~ "yes_to_no",
          TRUE ~ NA_character_
        ))
    } else {
    batch_counts <- tibble()
    }
    
    # ONE-BY-ONE: overlapping sites may span multiple tiles,
    # so build a mosaic per buffer then use exactextractr
    if (nrow(idx_one_by_one) > 0) {
      
      one_by_one_list <- list()
      centroid_obo_list <- list()

      combos <- unique(idx_one_by_one$identity)

      for (i in seq_len(length(combos))) {

        b <- idx_one_by_one %>%
          filter(identity == combos[i]) %>%
          distinct(identity, .keep_all = TRUE)

        b_sf <- sf::st_as_sf(b)
        centroid_vals_scenarios <- list(ref = FALSE, proj = FALSE)

        for (scenario in c("ref", "proj")) {

          tiles <- if (scenario == "ref") ref_tiles else proj_tiles

          tiles_words <- overlap_files %>%
            filter(identity == b_sf$identity) %>%
            pull(if (scenario == "ref") ref_file else proj_file)

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

      one_by_one_counts <- bind_rows(one_by_one_list) %>%
        pivot_wider(names_from = scenario,
                    values_from = c(total, flooded),
                    names_glue = "raw_{.value}_cells_{scenario}") %>%
        mutate(raw_total_cells_ref = if_else(is.na(raw_total_cells_ref), raw_total_cells_proj, raw_total_cells_ref)) %>%
        left_join(bind_rows(centroid_obo_list), by = "identity")
    } else {
      one_by_one_counts <- tibble()
    }
    
    region_results <- bind_rows(one_by_one_counts, batch_counts) %>%
      left_join(region_sf, by = c("buffer_km", "identity")) %>%
      mutate(raw_coverage_pct_ref = if_else(raw_total_cells_ref > 0,
                                            (raw_flooded_cells_ref / raw_total_cells_ref) * 100,
                                            NA_real_),
             raw_coverage_pct_proj = if_else(raw_total_cells_proj > 0,
                                             (raw_flooded_cells_proj / raw_total_cells_proj) * 100,
                                             NA_real_),
             raw_cell_change = raw_flooded_cells_proj - raw_flooded_cells_ref,
             raw_pct_point_change = raw_coverage_pct_proj - raw_coverage_pct_ref,
             inun_change = case_when(raw_coverage_pct_ref == 0 & raw_coverage_pct_proj == 0 ~ "no_to_no",
                                     raw_coverage_pct_ref == 0 & raw_coverage_pct_proj > 0  ~ "no_to_some",
                                     raw_coverage_pct_ref > 0  & raw_coverage_pct_proj == 0 ~ "some_to_no",
                                     raw_coverage_pct_ref > 0  & raw_coverage_pct_proj > 0  ~ "some_to_some",
                                     TRUE ~ NA_character_)) %>%
      select(identity, wsd_source_id, buffer_km,
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
    st_drop_geometry() %>%
    select(wsd_source_id, inun_proj_mm, inun_proj_closest_mm, depth_ft) %>%
    distinct()

  # Processed coastal results with metadata
  coastal_results <- ws_inun_stats %>%
    left_join(ws_geom_meta, by = "wsd_source_id") %>%
    select(wsd_source_id, buffer_km,
           inun_proj_mm, inun_proj_closest_mm, depth_ft,
           raw_total_cells_ref, raw_flooded_cells_ref, raw_coverage_pct_ref,
           raw_total_cells_proj, raw_flooded_cells_proj, raw_coverage_pct_proj,
           raw_cell_change, raw_pct_point_change, inun_change, inun_centroid_change) %>%
    mutate(across(c(raw_total_cells_ref, raw_flooded_cells_ref,
                    raw_total_cells_proj, raw_flooded_cells_proj, raw_coverage_pct_proj,
                    raw_cell_change, raw_pct_point_change),
                  ~ ifelse(inun_change == "no_to_no", 0, .x))) %>%
    distinct()

  # All unique water supplies (base for expansion)
  all_ws_base <- ws_scenarios %>%
    st_drop_geometry() %>%
    select(wsd_source_id, park_unit) %>%
    distinct()

  # Expand to all water supplies; non-coastal rows get NA for all indicator columns.
  # group_by(buffer_km): non-coastal rows have buffer_km = NA → their own group;
  final_results <- all_ws_base %>%
    left_join(coastal_results, by = "wsd_source_id") %>%
    group_by(buffer_km) %>%
    add_standardized_indicators("raw_coverage_pct_ref") %>%
    add_standardized_indicators("raw_coverage_pct_proj") %>%
    add_standardized_indicators("raw_cell_change") %>%
    add_standardized_indicators("raw_pct_point_change") %>%
    ungroup() %>%
    select(wsd_source_id, park_unit, buffer_km,
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
    write_csv(final_results, filename)
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
