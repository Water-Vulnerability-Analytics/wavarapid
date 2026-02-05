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
#' @param water_supplies A modified water_supplies sf feature that must include
#'   at minimum the following columns:
#'   \describe{
#'     \item{\code{inun_region}}{Grouping variable used to split polygons into regions.}
#'     \item{\code{state}}{State subdirectory name used to locate rasters.}
#'     \item{\code{ref_file}}{Filename of the reference inundation raster.}
#'     \item{\code{proj_file}}{Filename of the projected inundation raster.}
#'   }
#' @param save logical, whether to save output
#' @param output_dir character, directory path for saving files
#' @param verbose logical, whether to show progress messages
calc_inundation_exposure <- function(water_supplies,
                                     inundation_raster_data_dir,
                                     save = FALSE,
                                     output_dir = NULL,
                                     verbose = TRUE) {
  
  # Input validation
  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir)
  )
  
  # If verbose start time 
  if (verbose) {start_time <- Sys.time()} 
  
  region_list <- water_supplies %>%
    split(f = .$inun_region) 
  
  # If verbose start progress bar for each region
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "Processing inundation exposure [:bar] :percent (:current/:total) ETA: :eta",
      total = length(region_list),
      clear = FALSE,
      width = 80
    )
  } 
  
  # Process each water supply by inundation surface region 
  ws_inun_stats <- map_dfr(region_list, function(region_sf) {
    
    if (verbose) pb$tick()
    
    # Setup ====
    crs_val <- "+proj=longlat +datum=WGS84 +no_defs"
    
    ## Build the file paths for the inundation surface rasters
    state <- unique(region_sf$state)
    ref_file <- unique(region_sf$ref_file)
    proj_file <- unique(region_sf$proj_file)
    
    ref_f <- here(inundation_raster_data_dir, state, ref_file)
    proj_f <- here(inundation_raster_data_dir, state, proj_file)
    
    ## Convert the sf object to terra
    ws_terra <- terra::vect(region_sf)
    
    # Load and extract raster information ====
    extract_lists <- map(
      list(ref_f, proj_f),
      function(rast_f) {
        rast <- terra::rast(rast_f)
        terra::crs(rast) <- crs_val
        extracts <- terra::extract(rast, ws_terra)
        rm(rast)
        gc()
        Sys.sleep(1)
        return(extracts)
      }
    )
    
    # Calculate statistics ====
    ref_extracts <- extract_lists[[1]]  
    proj_extracts <- extract_lists[[2]] 
    
    ref_extract_stats <- ref_extracts %>%
      group_by(ID) %>% 
      summarise(
        total_cells = n(), # All cells touching polygon
        flooded_cells = sum(.data[[names(ref_extracts)[2]]], na.rm = T), # Flooded cells
        .groups = "drop"
      ) %>%
      ungroup()
    
    proj_extract_stats <- proj_extracts %>%
      group_by(ID) %>% 
      summarise(
        total_cells = n(), # All cells touching polygon
        flooded_cells = sum(.data[[names(proj_extracts)[2]]], na.rm = T), # Flooded cells
        .groups = "drop"
      ) %>%
      ungroup()
    
    all_extract_stats <- left_join(
      ref_extract_stats, proj_extract_stats,
      by = "ID", suffix = c("_ref", "_proj")
    ) %>%
      mutate(
        coverage_pct_ref = (flooded_cells_ref/total_cells_ref) * 100,
        coverage_pct_proj = (flooded_cells_proj/total_cells_proj) * 100,
        cell_change = flooded_cells_proj - flooded_cells_ref, 
        pct_point_change = coverage_pct_proj - coverage_pct_ref,
        pct_change = ifelse(coverage_pct_ref == 0, 
                            ifelse(coverage_pct_proj == 0, 0, Inf), 
                            pct_point_change/coverage_pct_ref * 100) 
      ) %>%
      # organize columns and define them
      select(
        ID,                    # Water supply polygon identifier from extract function
        total_cells_ref,       # Total raster cells intersecting polygon in reference scenario
        flooded_cells_ref,     # Count of flooded cells in reference scenario
        coverage_pct_ref,      # Reference flood coverage as proportion of intersecting area
        total_cells_proj,      # Total raster cells intersecting polygon in projection scenario  
        flooded_cells_proj,    # Count of flooded cells in projection scenario
        coverage_pct_proj,     # Projection flood coverage as proportion of intersecting area
        cell_change,           # Net change in flooded cell count (proj - ref)
        pct_point_change,      # Percentage point difference in coverage (proj - ref)
        pct_change             # Relative percentage change from reference baseline
      )
    
    # Join all extract statistics with the water supply region data
    result <- region_sf %>%
      mutate(row_id = row_number()) %>%
      left_join(all_extract_stats, by = c("row_id" = "ID")) %>%
      select(-row_id) %>%
      sf::st_drop_geometry()
    
    return(result)
  })
  
  # Clean up final results ====
  final_results <- ws_inun_stats %>%
    select(ID, wsd_source_id, 
           buffer_km, # NOTE: We do not pivot based on the buffer km here.
           inun_proj_mm, inun_proj_closest_mm, depth_ft, 
           inun_region, inun_surf_url, proj_file, ref_file,
           total_cells_ref, flooded_cells_ref, coverage_pct_ref,
           total_cells_proj, flooded_cells_proj, coverage_pct_proj,
           cell_change, pct_point_change, pct_change) %>%
    mutate(
      across(c(total_cells_ref, flooded_cells_ref, coverage_pct_ref,
               total_cells_proj, flooded_cells_proj, coverage_pct_proj,
               cell_change, pct_point_change, pct_change),
             ~ ifelse(pct_change == 0, NA, .x))
    )
  
  # If save, save results ====
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(output_dir, paste0("inundation_exposure_", Sys.Date(), ".csv"))
    write_csv(final_results, filename)
    
    if (verbose) {
      message("Results saved to: ", filename)
    }
  }
  
  # If verbose end time and progress bar ====
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message("\nInundation exposure processing completed in ", 
            round(total_time, 2), " ", attr(total_time, "units"))
  }
  
  return(final_results)
}
