#' Process coastal inundation inputs for water-supply risk analysis
#'
#' Buffers water supply point locations (1 km and 4 km), extracts Interagency Task
#' Force (ITF) sea-level-rise projections for 2050 (83rd percentile), builds a master
#' NOAA Sea Level Rise Viewer (SLRV) inundation-extent URL table from a compiled text list,
#' determines the relevant NOAA inundation region per water supply using a two-step process:
#' first a same-state spatial join to keep all candidate regions, then distance-based
#' disambiguation for sites with multiple overlapping candidates — selects
#' the nearest available NOAA inundation extent depth that matches/exceeds the ITF projection,
#' downloads the needed NOAA extent rasters, and writes intermediate and final geopackages/rasters.
#' Uses ITF 83rd percentile projections for 2050. Excludes AK, HI, Guam, CNMI, USVI.
#' Matches ITF projections to nearest available NOAA inundation depths.
#' Downloads reference (0 ft SLR) and scenario rasters as needed.
#'
#' Data sources:
#' - ITF projections: https://zenodo.org/records/5951626
#' - Technical report: https://earth.gov/sealevel/us/internal_resources/756/noaa-nos-techrpt01-global-regional-SLR-scenarios-US.pdf
#' - NOAA SLRV data: https://coast.noaa.gov/slrdata/
#' - Interactive tool: https://sealevel.nasa.gov/task-force-scenario-tool
#'
#' @param coast_dir_path Character path to the coastal comparison data directory that contains
#'   ITF NetCDF inputs and NOAA SLRV URL list inputs, and where outputs will be written.
#'   Defaults to `"data/all/coast_comp_data"`.
#' @param water_supplies An `sf` object of point geometries representing water supply locations.
#'   Must include a `state` column (2-letter state code) used to group downloads/processing.
#' @param output_dir Directory to save the datasets that will be used as inputs for the
#'   calculation of inundation exposure.

process_inundation_data <- function(coast_dir_path = here("data", "all", "coast_comp_data"),
                                    water_supplies,
                                    output_dir = here("data", "rapid", "inputs", "coast_comp")) {
  
  # 1) Buffer water supply points at 1 km and 4 km
  
  transformed_ws <- water_supplies %>%
    st_transform(5070) %>% # NAD83 / Conus Albers (meters)
    rowid_to_column(var = "ID")
  
  buffer_1km_ws <- transformed_ws %>%
    st_buffer(1000) %>%
    st_transform(4326) %>%
    mutate(buffer_km = 1)
  
  buffer_4km_ws <- transformed_ws %>%
    st_buffer(4000) %>%
    st_transform(4326) %>%
    mutate(buffer_km = 4)
  
  all_buffer_ws <- bind_rows(buffer_1km_ws, buffer_4km_ws)
  
  sf::st_write(all_buffer_ws,
               here(coast_dir_path, "altered_water_supplies", "buff_ws.gpkg"),
               append = FALSE)
  
  rm(buffer_1km_ws, buffer_4km_ws)
  gc()
  
  # 2) Load and process ITF SLR projections (83rd pct, 2050)
  
  inundation_projections_path <- here(coast_dir_path, "interagency_report",
                                      "Results", "TR_gridded_projections.nc")
  
  projections_raster <- terra::rast(inundation_projections_path)
  terra::crs(projections_raster) <- "+proj=longlat +datum=WGS84 +no_defs"
  
  projection_layer <- projections_raster[["rsl_total_IntHigh_percentiles=83_years=2050"]]
  proj_layer_rotated <- terra::rotate(projection_layer)
  
  # Values <= -20000 represent no-change and are set to NA
  final_itf_slr_proj <- terra::ifel(proj_layer_rotated <= -20000, NA, proj_layer_rotated)
  
  terra::writeRaster(final_itf_slr_proj,
                     here(coast_dir_path, "interagency_report", "generated_products",
                          "itf_slr_projections_2050_83rd.tif"),
                     overwrite = TRUE)
  
  rm(projections_raster, projection_layer, proj_layer_rotated)
  gc()
  
  # ---- Coastal subset for visualisation only (not used in analysis) ----
  
  conus_bbox <- sf::st_as_sfc(sf::st_bbox(c(xmin = -125, ymin = 20, xmax = -66, ymax = 50)))
  sf::st_crs(conus_bbox) <- 4269
  
  us_coastline <- tigris::coastline(year = 2020)
  conus_coastline <- sf::st_crop(us_coastline, conus_bbox) %>%
    sf::st_transform(crs = 4326)
  
  projection_coast_mask <- terra::mask(final_itf_slr_proj, conus_coastline, touches = TRUE)
  data_mask <- !is.na(projection_coast_mask)
  expanded_mask <- terra::focal(data_mask, w = 3, fun = max, na.rm = TRUE)
  projections_expanded <- terra::mask(final_itf_slr_proj, expanded_mask, maskvalue = 0)
  relevant_projections <- terra::trim(projections_expanded)
  
  terra::writeRaster(relevant_projections,
                     here(coast_dir_path, "interagency_report", "generated_products",
                          "itf_slr_proj_coast_subset.tif"),
                     overwrite = TRUE)
  
  rm(us_coastline, conus_bbox, conus_coastline, projection_coast_mask,
     data_mask, expanded_mask, projections_expanded, relevant_projections)
  gc()
  
  # 3) Extract ITF projections for each buffered water supply
  
  water_supply_inundation_values <- terra::extract(final_itf_slr_proj, transformed_ws) 
  
  ws_itf_slr_proj <- left_join(all_buffer_ws, water_supply_inundation_values, by = "ID") %>%
    rename(inun_proj_mm = `rsl_total_IntHigh_percentiles=83_years=2050`)
  
  sf::st_write(ws_itf_slr_proj,
               here(coast_dir_path, "altered_water_supplies", "buff_ws_itf_slr_proj.gpkg"),
               append = FALSE)
  
  rm(water_supply_inundation_values)
  gc()
  
  # 4) Build master NOAA SLRV URL table
  
  tif_txt_path <- here(coast_dir_path, "NOAA_SLRV_2022", "ALL_extent_raster_URLlist.txt")
  tif_vec <- read_lines(tif_txt_path)
  
  tif_df <- tibble(url = tif_vec) %>%
    mutate(raw_file = basename(url),
           file = str_replace(raw_file, "_Redacted", ""),
           info = str_split(str_remove(file, ".tif"), "_"),
           state = map_chr(info, ~ .x[1]),
           depth_ft = map(info, ~ {
             n <- length(.x)
             if (n >= 2) {
               last_two   <- .x[(n-1):n]
               is_numeric <- str_detect(last_two, "^\\d+(\\.\\d+)?$")
               if (all(is_numeric))  return(as.numeric(str_flatten(last_two, collapse = ".")))
               if (is_numeric[2])    return(as.numeric(.x[n]))
             } else if (n >= 1 && str_detect(.x[n], "^\\d+(\\.\\d+)?$")) {
               return(as.numeric(.x[n]))
             }
             return(NA_real_)
           })
    ) %>%
    select(-info) %>%
    filter(!is.na(depth_ft)) %>%
    mutate(depth_ft  = as.numeric(depth_ft),
           depth_m   = depth_ft * 0.3048,
           depth_mm  = round(depth_m * 1000, 1),
           inun_region = str_remove(file, "_[0-9_]+\\.tif$")) %>%
    select(url, file = raw_file, state, inun_region, depth_ft, depth_m, depth_mm) %>%
    filter(!state %in% c("AK", "AS", "HI", "MP", "PR", "VI", "Guam", "CNMI", "USVI"))
  
  write_csv(tif_df, here(coast_dir_path, "NOAA_SLRV_2022", "ALL_extent_raster_URLlist.csv"))
  
  # 5) Download reference (0 ft SLR) inundation rasters
  
  #' @noRd
  download_NOAA_inun_slr_tif_data <- function(tif_url_df,
                                              save_dir = here(output_dir, "Extent_Rasters")) {
    tif_url_df %>%
      pwalk(function(url, file, state, ...) {
        dest_dir  <- here(save_dir, state)
        dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
        dest_file <- here(save_dir, state, file)
        
        for (i in 1:3) {
          tryCatch({
            download.file(url = url, destfile = dest_file, mode = "wb")
            cat("\nDownload succeeded (attempt:", i, ") for:", file, "\n")
            break
          }, error = function(e) {
            cat("\nDownload failed (attempt:", i, ") for:", file, "\n")
            Sys.sleep(i)
            if (i == 3) cat("\n⚠️ All 3 attempts failed for:", file, "\n")
          })
        }
        
        gc()
        Sys.sleep(1)
      })
  }
  
  reference_tif_df <- tif_df %>% filter(depth_mm == 0.0)
  
  download_NOAA_inun_slr_tif_data(tif_url_df = reference_tif_df,
                                  save_dir = here(output_dir, "Extent_Rasters"))
  
  # 6) Build inundation region bounding boxes from reference rasters
  
  raster_f <- list.files(path = here(output_dir, "Extent_Rasters"),
                         pattern = "_0(_0)?(_Redacted)?\\.tif$",
                         full.names = TRUE,
                         recursive = TRUE)
  
  inun_region_bboxs <- map(raster_f, function(f) {
    region_name <- tools::file_path_sans_ext(basename(f))
    region_name <- sub("connectRaster.*", "connectRaster", region_name)
    
    rast    <- terra::rast(f)
    bbox_sf <- sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(rast))) %>%
      sf::st_set_crs(terra::crs(rast)) %>%
      sf::st_transform(4326) %>%
      mutate(inun_region = region_name)
    
    rm(rast)
    Sys.sleep(1)
    gc()
    
    bbox_sf
  }) %>%
    bind_rows() %>%
    mutate(state = str_extract(inun_region, "^[A-Z]{2}")) %>%
    distinct()
  
  sf::st_write(inun_region_bboxs,
               here(coast_dir_path, "NOAA_SLRV_2022", "inundation_data", "inun_region_bboxs.gpkg"),
               delete_dsn = TRUE)
  
  
  # 7) Spatial join: water supplies -> inundation regions
  
  rel_state_codes <- unique(inun_region_bboxs$state)
  
  inun_rel_ws <- ws_itf_slr_proj %>%
    filter(state %in% rel_state_codes, !is.na(inun_proj_mm))
  
  # spatial join — may produce multiple rows per site where bboxes overlap
  ws_inun_joined <- sf::st_join(inun_rel_ws,
                                inun_region_bboxs %>%
                                  rename(inun_state = state) %>%
                                  select(inun_region, inun_state),
                                join = sf::st_intersects,
                                left = TRUE)
  
  ws_inun_regions <- ws_inun_joined %>%
    # Nullify projection if no region matched
    mutate(inun_proj_mm = ifelse(is.na(inun_region), NA, inun_proj_mm),
           # Map ITF projection to nearest available NOAA inundation depth increment
           inun_proj_closest_mm = case_when(inun_proj_mm <= 0.0    ~ 0.0,
                                            inun_proj_mm <= 152.4  ~ 152.4,
                                            inun_proj_mm <= 304.8  ~ 304.8,
                                            inun_proj_mm <= 457.2  ~ 457.2,
                                            inun_proj_mm <= 609.6  ~ 609.6,
                                            inun_proj_mm <= 762.0  ~ 762.0,
                                            inun_proj_mm <= 914.4  ~ 914.4,
                                            inun_proj_mm <= 1066.8 ~ 1066.8,
                                            inun_proj_mm <= 1219.2 ~ 1219.2,
                                            inun_proj_mm <= 1371.6 ~ 1371.6,
                                            inun_proj_mm <= 3048.0 ~ 3048.0,
                                            TRUE ~ NA_real_)) %>%
    relocate(inun_proj_closest_mm, .after = inun_proj_mm)
  
  rm(ws_inun_joined, inun_rel_ws)
  gc()
  
  
  # 8) Download scenario inundation rasters
  
  relevant_ws_inun_proj <- ws_inun_regions %>%
    drop_na(inun_proj_closest_mm) %>%
    select(depth_mm = inun_proj_closest_mm, inun_region) %>%
    # CA EKA region only available in 1 ft increments
    mutate(depth_mm = ifelse(inun_region == "CA_EKA_connectRaster", 609.6, depth_mm),
           depth_mm = round(depth_mm, 1)) %>%
    distinct(depth_mm, inun_region)
  
  download_tif_df <- relevant_ws_inun_proj %>%
    left_join(tif_df, by = c("inun_region", "depth_mm"))
  
  download_NOAA_inun_slr_tif_data(tif_url_df = download_tif_df,
                                  save_dir = here(output_dir, "Extent_Rasters"))
  
  # 9) Assemble final water supply dataset with all inundation columns
  
  # Columns added:
  # buffer_km             - buffer size in km
  # inun_proj_mm          - ITF projected SLR in mm
  # inun_proj_closest_mm  - nearest available NOAA depth increment >= ITF projection
  # depth_ft / depth_m    - NOAA inundation depth in original units
  # inun_region           - NOAA inundation region identifier
  # inun_surf_url         - URL to the projected scenario raster
  # proj_file             - local filename of the projected scenario raster
  # ref_file              - local filename of the reference (0 ft SLR) raster
  
  updated_ws_inun_proj <- ws_itf_slr_proj %>%
    inner_join(ws_inun_regions %>% st_drop_geometry() %>% select(ID, buffer_km, inun_region, inun_proj_closest_mm),
              by = c("ID", "buffer_km")) %>%
    left_join(download_tif_df,
              by = c("state", "inun_region", "inun_proj_closest_mm" = "depth_mm")) %>%
    left_join(reference_tif_df %>% select(state, inun_region, reference_file = file),
              by = c("inun_region", "state")) %>%
    rename(inun_surf_url = url, proj_file = file, ref_file = reference_file) %>%
    relocate(depth_ft, depth_m, .after = inun_proj_closest_mm) %>%
    relocate(geometry, .after = last_col()) %>%
    filter(!is.na(ref_file))
  
  sf::st_write(updated_ws_inun_proj,
               here(output_dir, "buff_ws_inun_regions.gpkg"),
               append = FALSE)
  
  rm(ws_inun_regions, download_tif_df, reference_tif_df, relevant_ws_inun_proj,
     inun_region_bboxs, raster_f, rel_state_codes, tif_df, all_buffer_ws,
     final_itf_slr_proj, ws_itf_slr_proj)
  gc()
  
  return(updated_ws_inun_proj)
}