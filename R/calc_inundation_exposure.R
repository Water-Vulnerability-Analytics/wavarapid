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
#' @param inundation_raster_data_dir path to inundation surface raster state subdirectories
#' @param save logical, whether to save output
#' @param output_dir character, directory path for saving files
#' @param verbose logical, whether to show progress messages

calc_inundation_exposure <- function(water_supplies,
                                     inundation_raster_data_dir = here("data", "rapid", "inputs", "coast_comp", "Extent_Rasters"),
                                     save = FALSE, 
                                     output_dir = NULL,
                                     verbose = TRUE) {
  
  # Input validation
  required_columns <- c("inun_region", "state", "ref_file", "proj_file")
  
  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir),
    "water_supplies missing required columns" = all(required_columns %in% names(water_supplies))
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
        rast <- terra::crop(rast, terra::project(ws_terra, rast))
        extracts <- terra::extract(rast, ws_terra)
        rm(rast)
        gc()
        return(extracts)
      }
    )
    
    # Calculate statistics ====
    ref_extracts <- extract_lists[[1]]
    proj_extracts <- extract_lists[[2]]
    rm(extract_lists)
    gc()

    ref_extract_stats <- ref_extracts %>%
      group_by(ID) %>%
      summarise(
        raw_total_cells = n(), # All cells touching polygon
        raw_flooded_cells = sum(.data[[names(ref_extracts)[2]]], na.rm = T), # Flooded cells
        .groups = "drop"
      )

    proj_extract_stats <- proj_extracts %>%
      group_by(ID) %>%
      summarise(
        raw_total_cells = n(), # All cells touching polygon
        raw_flooded_cells = sum(.data[[names(proj_extracts)[2]]], na.rm = T), # Flooded cells
        .groups = "drop"
      )

    all_extract_stats <- left_join(
      ref_extract_stats, proj_extract_stats,
      by = "ID", suffix = c("_ref", "_proj")
    ) %>%
      mutate(
        raw_coverage_pct_ref = (raw_flooded_cells_ref / raw_total_cells_ref) * 100,
        raw_coverage_pct_proj = (raw_flooded_cells_proj / raw_total_cells_proj) * 100,
        raw_cell_change = raw_flooded_cells_proj - raw_flooded_cells_ref,
        raw_pct_point_change = raw_coverage_pct_proj - raw_coverage_pct_ref,
        raw_pct_change = ifelse(raw_coverage_pct_ref == 0,
                                ifelse(raw_coverage_pct_proj == 0, 0, Inf),
                                raw_pct_point_change / raw_coverage_pct_ref * 100),
        inun_change = case_when(
          raw_coverage_pct_ref == 0 & raw_coverage_pct_proj == 0 ~ "no_to_no",   # No inundation to no inundation
          raw_coverage_pct_ref == 0 & raw_coverage_pct_proj > 0 ~ "no_to_some",  # No inundation to some inundation
          raw_coverage_pct_ref > 0 & raw_coverage_pct_proj == 0 ~ "some_to_no",  # Some inundation to no inundation
          raw_coverage_pct_ref > 0 & raw_coverage_pct_proj > 0  ~ "some_to_some",# Some inundation to some inundation
          TRUE ~ NA_character_                                                     # Handle any edge cases
        )
      ) %>%
      # organize columns and define them
      select(
        ID,                     # Water supply polygon identifier from extract function
        raw_total_cells_ref,    # Total raster cells intersecting polygon in reference scenario
        raw_flooded_cells_ref,  # Count of flooded cells in reference scenario
        raw_coverage_pct_ref,   # Reference flood coverage as proportion of intersecting area
        raw_total_cells_proj,   # Total raster cells intersecting polygon in projection scenario  
        raw_flooded_cells_proj, # Count of flooded cells in projection scenario
        raw_coverage_pct_proj,  # Projection flood coverage as proportion of intersecting area
        raw_cell_change,        # Net change in flooded cell count (proj - ref)
        raw_pct_point_change,   # Percentage point difference in coverage (proj - ref)
        raw_pct_change,         # Relative percentage change from reference baseline
        inun_change             # Categorical indicator for inundation change patterns
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
           buffer_km, # NOTE: We do not pivot based on the buffer km here. This data has 1km and 4km stacked.
           inun_proj_mm, inun_proj_closest_mm, depth_ft, 
           inun_region, inun_surf_url, proj_file, ref_file,
           raw_total_cells_ref, raw_flooded_cells_ref, raw_coverage_pct_ref,
           raw_total_cells_proj, raw_flooded_cells_proj, raw_coverage_pct_proj,
           raw_cell_change, raw_pct_point_change, raw_pct_change, inun_change) %>%
    mutate(
      across(c(raw_total_cells_ref, raw_flooded_cells_ref, raw_coverage_pct_ref,
               raw_total_cells_proj, raw_flooded_cells_proj, raw_coverage_pct_proj,
               raw_cell_change, raw_pct_point_change, raw_pct_change),
             ~ ifelse(inun_change == "no_to_no", NA, .x))
    ) %>%
    add_standardized_indicators("raw_coverage_pct_ref") %>%
    add_standardized_indicators("raw_coverage_pct_proj") %>%
    add_standardized_indicators("raw_cell_change") %>%
    add_standardized_indicators("raw_pct_point_change") %>%
    # Percent change is a special case where there could be infinite values in it 
    # due to the no inundation to some inundation situation.
    mutate(
      # Temporarily replace Inf with NA for ranking calculations
      raw_pct_change = ifelse(is.infinite(raw_pct_change), NA, raw_pct_change)
    ) %>%
    add_standardized_indicators(prefix = "raw_pct_change") %>%
    mutate(
      # Put the Inf values back in the original column after ranking
      raw_pct_change = ifelse(is.na(raw_pct_change) & inun_change == "no_to_some", Inf, raw_pct_change)
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
