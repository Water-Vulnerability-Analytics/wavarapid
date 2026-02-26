#' Calculate inundation exposure metrics for water-supply polygons
#'
#' Computes inundation exposure statistics for each water-supply polygon by
#' extracting values from a reference inundation raster and a projected inundation
#' raster, grouped by inundation surface region (e.g., state/region file pairs).
#' For each polygon, the function calculates total intersecting raster cells,
#' flooded cell counts, percent coverage in each scenario, absolute and percentage
#' point changes, and a relative percent change (with special handling when the
#' reference coverage is zero).
#'
#' Non-overlapping sites (count <= 2) fall on a single raster tile and are
#' extracted all at once using exactextractr for speed. Sites with overlapping
#' geometries (count > 2) are processed individually with mosaicking across
#' tiles as needed.
#'
#' @param ws_scenarios A modified water_supplies sf feature that must include
#'   at minimum the following columns:
#'   \describe{
#'     \item{\code{inun_region}}{Grouping variable used to split polygons into regions.}
#'     \item{\code{state}}{State subdirectory name used to locate rasters.}
#'     \item{\code{ref_file}}{Filename of the reference inundation raster.}
#'     \item{\code{proj_file}}{Filename of the projected inundation raster.}
#'   }
#' @param inundation_raster_data_dir path to inundation surface raster state subdirectories
#' @param save logical, whether to save output
#' @param output_dir character, directory path for saving files
#' @param verbose logical, whether to show progress messages

calc_inundation_exposure <- function(ws_scenarios = sf::st_read(here("data", "rapid", "inputs", "coast_comp", 
                                                                     "buff_ws_inun_regions.gpkg")) %>%
                                       # Filter for water supplies that have relevant inundation surface information
                                       filter(!is.na(proj_file)),
                                     inundation_raster_data_dir = here("data", "rapid", "inputs", "coast_comp", "Extent_Rasters"),
                                     save = FALSE,
                                     output_dir = NULL,
                                     verbose = TRUE) {
  
  # for testing:
  # inun_region_bboxs <- st_read(here("data/all/coast_comp_data/", "NOAA_SLRV_2022", "inundation_data", "inun_region_bboxs.gpkg"))
  
  # Input validation
  required_columns <- c("inun_region", "state", "ref_file", "proj_file")
   
  stopifnot("ws_scenarios must be an sf object" = inherits(ws_scenarios, "sf"),
            "save is TRUE but output_dir not provided" = !save || !is.null(output_dir),
            "ws_scenarios missing required columns" = all(required_columns %in% names(ws_scenarios)))
  
  if (verbose) start_time <- Sys.time()
  
  ws_scenarios_mod <- ws_scenarios %>% 
    ungroup() %>% 
    distinct() %>%
    # for testing -
    # filter(wsd_source_id == "CALO_WS_007_01") %>%
    mutate(identity = paste0(wsd_source_id, "-", buffer_km)) %>%
    select(c(park_unit, state, identity, wsd_source_id, buffer_km, inun_proj_mm,
             inun_region ,inun_proj_closest_mm, depth_ft, depth_m, inun_surf_url, 
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
    
    # for testing - 
    # region_sf <- region_list[[1]]

    if (verbose) pb$tick()
    
    # tiles needed for this region 
    locs_needed <- ws_scenarios_mod %>%
      filter(identity %in% region_sf$identity) %>%
      distinct(proj_file, ref_file) %>%
      mutate(state = substr(proj_file, 1, 2))
    
    state <- unique(locs_needed$state)
    
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
    
    # build full paths to all tiles
    ref_tiles  <- sort(unique(here(inundation_raster_data_dir, state, locs_needed$ref_file)))
    proj_tiles <- sort(unique(here(inundation_raster_data_dir, state, locs_needed$proj_file)))
    
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
      
      for (scenario in c("ref", "proj")) {
        
        file_col  <- if (scenario == "ref") "ref_file" else "proj_file"
        
        tile_file <- unique(region_sf[[file_col]])
        
        stopifnot(length(tile_file) == 1)
        
        tile_path <- here(inundation_raster_data_dir, state, tile_file)
        
        stopifnot(file.exists(tile_path))
        
        rast <- terra::rast(tile_path)
        
        extracted <- exactextractr::exact_extract(rast, batch_sf)
        
        counts_list[[scenario]] <- data.frame(identity = batch_sf$identity,
                                              buffer_km = batch_sf$buffer_km,
                                              total = vapply(extracted, \(df) nrow(df), numeric(1)),
                                              flooded = vapply(extracted, \(df) sum(!is.na(df$value)), numeric(1)),
                                              scenario = scenario)
      }
      
      batch_counts <- bind_rows(counts_list) %>%
        pivot_wider(names_from  = scenario,
                    values_from = c(total, flooded),
                    names_glue  = "raw_{.value}_cells_{scenario}")
    } else {
    batch_counts <- tibble()
    }
    
    
    # ONE-BY-ONE: overlapping sites may span multiple tiles,
    # so build a mosaic per buffer then use exactextractr
    if (nrow(idx_one_by_one) > 0) {
      
      one_by_one_list <- list()
      
      combos <- unique(idx_one_by_one$identity)
      
      for (i in seq_len(length(combos))) {
        
        b <- idx_one_by_one %>% 
          filter(identity == combos[i]) %>% 
          distinct(identity, .keep_all = TRUE)
        
        b_sf <- sf::st_as_sf(b)
        
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
          
          # fxn to crop and mask a single raster tile to a buffer polygon
          prep_tile <- function(r_path, b_sf) {
            
            r0 <- terra::rast(r_path)
            b_projected <- sf::st_transform(b_sf, terra::crs(r0))
            
            r1 <- tryCatch(terra::crop(r0, b_projected), error = function(e) NULL)
            if (is.null(r1) || !inherits(r1, "SpatRaster")) return(NULL)
            if (terra::ncell(r1) == 0 || all(is.na(terra::values(r1)))) return(NULL)
            
            r2 <- tryCatch(terra::mask(r1, terra::vect(b_projected)), error = function(e) NULL)
            if (is.null(r2) || !inherits(r2, "SpatRaster")) return(NULL)
            if (all(is.na(terra::values(r2)))) return(NULL)
            
            return(r2)
            
          }
          
         tile_list <- list()
         
          for (r in tiles) {
            t <- tryCatch(prep_tile(r, b_sf), error = function(e) NULL)
            if (!is.null(t) && inherits(t, "SpatRaster")) tile_list[[length(tile_list) + 1]] <- t
          }
          if (length(tile_list) > 0) {
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
            
          } else {
            
            extracted <- exactextractr::exact_extract(out, 
                                                      sf::st_transform(b_sf, terra::crs(out)))[[1]]
            
            one_by_one_list[[paste(b_sf$identity, scenario, sep = "_")]] <- data.frame(identity = b_sf$identity,
                                                                                       buffer_km = b_sf$buffer_km,
                                                                                       total = nrow(extracted),
                                                                                       flooded = sum(!is.na(extracted$value)),
                                                                                       scenario = scenario)
          }
        }
      }
      
      one_by_one_counts <- bind_rows(one_by_one_list) %>%
        pivot_wider(names_from  = scenario,
                    values_from = c(total, flooded),
                    names_glue  = "raw_{.value}_cells_{scenario}") %>%
        mutate(raw_total_cells_ref = if_else(is.na(raw_total_cells_ref), raw_total_cells_proj, raw_total_cells_ref))
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
             raw_pct_change = case_when(is.na(raw_coverage_pct_ref) | is.na(raw_coverage_pct_proj) ~ NA_real_,
                                        raw_coverage_pct_ref == 0 & raw_coverage_pct_proj == 0 ~ 0,
                                        raw_coverage_pct_ref == 0 & raw_coverage_pct_proj > 0 ~ Inf,
                                        TRUE ~ (raw_pct_point_change / raw_coverage_pct_ref) * 100),
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
             raw_pct_change,         # Relative percentage change from reference baseline
             inun_change)            # Categorical indicator for inundation change patterns
    
    return(region_results)
    
  })
  
  # Clean up final results 
  final_results <- ws_inun_stats %>%
    left_join(ws_scenarios %>%
                select(park_unit, wsd_source_id, inun_proj_mm, inun_proj_closest_mm, depth_ft) %>%
                distinct(),
              by = "wsd_source_id") %>%
    select(wsd_source_id,
           buffer_km, # NOTE: We do not pivot based on buffer_km here. Data has 1km and 4km stacked.
           inun_proj_mm, inun_proj_closest_mm, depth_ft,
           raw_total_cells_ref, raw_flooded_cells_ref, raw_coverage_pct_ref,
           raw_total_cells_proj, raw_flooded_cells_proj, raw_coverage_pct_proj,
           raw_cell_change, raw_pct_point_change, raw_pct_change, inun_change) %>%
    mutate(across(c(raw_total_cells_ref, raw_flooded_cells_ref, raw_coverage_pct_ref,
                    raw_total_cells_proj, raw_flooded_cells_proj, raw_coverage_pct_proj,
                    raw_cell_change, raw_pct_point_change, raw_pct_change),
                  ~ ifelse(inun_change == "no_to_no", NA, .x))) %>%
    distinct() %>%
    add_standardized_indicators("raw_coverage_pct_ref") %>%
    add_standardized_indicators("raw_coverage_pct_proj") %>%
    add_standardized_indicators("raw_cell_change") %>%
    add_standardized_indicators("raw_pct_point_change") %>%
    # Percent change special case: Inf arises from no_to_some transitions.
    # Temporarily replace Inf with NA for ranking, then restore afterwards.
    mutate(raw_pct_change = ifelse(is.infinite(raw_pct_change), NA, raw_pct_change)) %>%
    add_standardized_indicators(prefix = "raw_pct_change") %>%
    mutate(raw_pct_change = ifelse(is.na(raw_pct_change) & inun_change == "no_to_some", Inf, raw_pct_change))
  
  # Save if requested ====
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
