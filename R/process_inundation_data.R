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
#'   Must include a `source_state` column (2-letter state code) used to group downloads/processing.
#' @param output_dir Directory to save the datasets that will be used as inputs for the
#'   calculation of inundation exposure.
process_inundation_data <- function(coast_dir_path = here("data", "all", "coast_comp_data"),
                                    water_supplies,
                                    output_dir = here("data", "rapid", "inputs", "coast_comp")) {
  
  # Generate the buffered water supply dataset 
  
  # For this dataset we will buffer the water supply points by 1 and 4 km. We will
  # use this dataset to do our inundation risk analysis.
  
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
  
  # Save the buffered water supplies
  sf::st_write(all_buffer_ws,
               here(coast_dir_path, "altered_water_supplies", "buff_ws.gpkg"), # Buffered water supplies
               append = FALSE)
  
  rm(buffer_1km_ws, buffer_4km_ws)
  gc()
  
  # Load and Process the Interagency Task Force (ITF) Sea Level Rise (SLR) Projections 
  
  # These data were downloaded from:
  #   [Global and Regional Sea Level Rise Scenarios for the United States: Updated Mean Projections and Extreme Water Level Probabilities Along U.S. Coastlines](https://zenodo.org/records/5951626)
  # 
  # The report that this data came from can be found here:
  #   [2022 Technical Report](https://earth.gov/sealevel/us/internal_resources/756/noaa-nos-techrpt01-global-regional-SLR-scenarios-US.pdf)
  # 
  # An interactive tool displaying these data can be found here:
  #   [Interagency Sea Level Rise Scenario Tool](https://sealevel.nasa.gov/task-force-scenario-tool)
  # 
  # These data were made by the [Interagency Sea Level Task Force](https://earth.gov/sealevel/us/about/interagency-sea-level-task-force/). This working group
  # consists of members from NOAA, NASA, JPL, USGS, EPA, FEMA, DoD, and USACE.
  # 
  # This data includes a 1km gridded raster of the 83rd percentile of ITF Sea 
  # Level Rise Projections for the year 2050 found in the "Results" section at the 
  # download link above. The full dataset includes U.S. Sea Level Rise Scenarios 
  # to 2150 by decade that incorporate estimates of vertical land motion. Data are 
  # available at 1-degree grids along the U.S. coastline and downscaled 
  # specifically at NOAA tide-gauge locations. Global mean sea level rise 
  # scenarios are downscaled with output directly from the United Nations 
  # Intergovernmental Panel on Climate Change (IPCC) Sixth Assessment Report 
  # (AR6; IPCC, 2021a) and regionalized for the U.S. coastline.
  
  ## Generate the complete gridded ITF SLR projections ----
  
  inundation_projections_path <- here(coast_dir_path, "interagency_report",
                                      "Results", "TR_gridded_projections.nc")
  
  # read in nc file as a raster with terra
  projections_raster <- terra::rast(inundation_projections_path)
  
  # assign projection manually
  terra::crs(projections_raster) <- "+proj=longlat +datum=WGS84 +no_defs"
  
  # pick the projection that we want from this raster
  # 83rd percentile projections out to the year 2050
  projection_layer <- projections_raster[["rsl_total_IntHigh_percentiles=83_years=2050"]]
  
  # rotate raster from 0–360 to -180–180
  proj_layer_rotated <- terra::rotate(projection_layer)
  
  # Convert data that essentially represent no change to NA values
  # This is the final manipulation before the projections get used for the analysis.
  final_itf_slr_proj <- terra::ifel(proj_layer_rotated <= -20000, NA, proj_layer_rotated)
  
  # Save the gridded projections as a tif
  terra::writeRaster(final_itf_slr_proj,
                     here(coast_dir_path, "interagency_report", "generated_products",
                          "itf_slr_projections_2050_83rd.tif"),
                     overwrite = TRUE)
  
  rm(projections_raster, projection_layer, proj_layer_rotated)
  gc()
  
  # Generate the water supply data with ITF SLR projections 
  
  # Now that we have the buffered water supply points and the ITF SLR gridded projections
  # we can join them together to find the SLR projections for each of our water supply points.
  
  # Extract the inundation projections specific to each water supply
  water_supply_inundation_values <- terra::extract(final_itf_slr_proj, transformed_ws)
  
  # Join water supply inundation projections to the water supply buffer data frame
  ws_itf_slr_proj <- left_join(all_buffer_ws, water_supply_inundation_values, by = "ID") %>%
    rename(inun_proj_mm = `rsl_total_IntHigh_percentiles=83_years=2050`)
  
  # Save the dataset
  sf::st_write(ws_itf_slr_proj,
               here(coast_dir_path, "altered_water_supplies", "buff_ws_itf_slr_proj.gpkg"), # Buffered water supplies with ITF SLR projections
               append = FALSE)
  
  rm(water_supply_inundation_values)
  gc()
  
  # Download the NOAA inundation surface datasets 
  
  # We will be using the data used in [NOAA's Sea Level Rise Viewer](https://coast.noaa.gov/slr/#/layer/slr/0/-11581024.663779823/5095888.569004184/4/satellite/none/0.8/2050/interHigh/noAccretion/NOS_Minor)
  # as our inundation surface datasets. These data are not stored in a way that ties them
  # to any specific scenario of SLR, instead they are the inundation surfaces given
  # some amount of SLR. This means that we will need to use the information from the
  # ITF SLR projections in order to find the correct inundation surface given the amount
  # of SLR the ITF projected for some region specific to each watershed.
  
  # The process for generating the NOAA inundation surface datasets is not simple.
  # This is because the inundation surface datasets are large and stored in
  # [Sea Level Rise Viewer Data Download](https://coast.noaa.gov/slrdata/).
  
  # There is no clear way to easily programmatically download data from this source.
  # They provide text data with URL links for the raster data, but that text data is also not
  # readily accessible via an API. Here is an [example](https://coast.noaa.gov/slrdata/Extent_Rasters/CA/URLlist_CA.txt)
  # of what one of these text files looks like. The process to programmatically download this data is,
  # broadly: compile these text files into one master text file, convert this master
  # text file into a form that is easily stored as a dataframe, use the information
  # in this dataframe to download the appropriate data via the provided URLs.
  
  # We will use this section to download the appropriate NOAA SLRV inundation surface datasets.
  
  ## Generate the master inundation dataset URL link dataframe ----
  
  # Here we convert the master inundation dataset URL text file into the master
  # inundation dataset URL dataframe. The inundation surface data is stored as feet
  # so we will need to convert it to millimeters to join it to the ITF SLR projections,
  # which are themselves in mm.
  
  # Note that the creation of the master inundation dataset URL text file was made outside
  # of R. I went into each state's download repository to copy the URL text file
  # for each state into the master text file.
  
  # File path for the compiled master inundation dataset URL text file
  tif_txt_path <- here(coast_dir_path, "NOAA_SLRV_2022", "ALL_extent_raster_URLlist.txt")
  
  # Read in the text file as a vector
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
           inun_region = str_remove(file, "_[0-9_]+\\.tif$"),
           inun_region = str_remove(inun_region, "_connectRaster|_connnectRaster")) %>%
    select(url, file = raw_file, state, inun_region, depth_ft, depth_m, depth_mm) %>%
    filter(!state %in% c("AK", "AS", "HI", "MP", "PR", "VI", "Guam", "CNMI", "USVI"))
  
  write_csv(tif_df, here(coast_dir_path, "NOAA_SLRV_2022", "ALL_extent_raster_URLlist.csv"))
  
  ## Download the reference inundation surface datasets ----
  
  # To determine the difference between the present day sea level conditions and
  # the conditions in our scenarios we need to download the surface inundation datasets
  # which predict 0 feet of SLR. We will use these datasets as our reference datasets.
  # We will download all of these since we do not yet know which inundation surface regions
  # will be relevant to us. We will use these reference datasets to figure that out.
  # Noted here is the fact that the Washington Pacific Coast Inundation Region is
  # tagged as redacted in the NOAA SLRV Data Download URL, but it is still unclear
  # why that is the case.
  
  # To do this we will use `download_NOAA_inun_slr_tif_data()`:
  
  ### DOWNLOADING NOAA SURFACE INUNDATION TIF DATASETS HELPER FUNCTION ----
  
  #' Download NOAA SLRV inundation surface tif data from master URL df
  #'
  #' @param tif_url_df df which holds the URLs for inundation surface raster tifs from NOAA SLRV
  #' @param save_dir Character path to directory where files should be saved
  
  download_NOAA_inun_slr_tif_data <- function(tif_url_df,
                                              save_dir = here(output_dir, "Extent_Rasters")) {
    tif_url_df %>%
      pwalk(function(url, file, state, ...) {
        
        # Create state directory if it doesn't exist
        dest_dir  <- here(save_dir, state)
        dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
        
        # Check if data exists, if it does, skip the download
        dest_file <- here(save_dir, state, file)
        
        # Try to download data 3 times
        for (i in 1:3) {
          tryCatch({
            download.file(url = url, destfile = dest_file, mode = "wb")
            cat("\nDownload succeeded (attempt:", i, ") for:", file, "\n")
            break # exit the for loop on success
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
  
  reference_tif_df <- tif_df %>% 
    group_by(inun_region) %>%
    slice_min(depth_ft)
  
  download_NOAA_inun_slr_tif_data(tif_url_df = reference_tif_df,
                                  save_dir = here(output_dir, "Extent_Rasters"))
  
  ## Download the scenario inundation surface datasets ----
  
  # First, we will use the reference inundation surface datasets to determine
  # the specific inundation region and SLR scenario for each of our water supplies.
  
  # Since it takes so long to pull in the data and RStudio cannot hold all of the
  # inundation surface rasters in memory at the same time we will do this part of the
  # analysis in two parts. First, we will create and save a shapefile of the boundaries
  # of the inundation regions. Then, we will spatially join these shapes with our
  # buffered watersheds. Finally, we will use this data to download the inundation
  # surface scenarios.
  
  # Make a list of the file paths of the reference inundation surfaces
  raster_f <- list.files(path = here(output_dir, "Extent_Rasters"),
                         full.names = TRUE,
                         recursive = TRUE)
  
  raster_f <- raster_f[basename(raster_f) %in% reference_tif_df$file]
  
  # Generate inundation region bounding boxes from the downloaded reference rasters
  inun_region_bboxs <- map(raster_f, function(f) {
    region_name <- tools::file_path_sans_ext(basename(f))
    region_name <- sub("_connectRaster.*", "", region_name)
    
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
  
  # Spatially join buffered water supplies against inundation region bounding boxes.
  # Sites with multiple candidates exist.
  ws_inun_joined <- sf::st_join(ws_itf_slr_proj,
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
  
  rm(ws_inun_joined)
  gc()
  
  # Second, we will use the region and SLR scenario information to download only the
  # relevant inundation surface datasets from NOAA.
  
  relevant_ws_inun_proj <- ws_inun_regions %>%
    drop_na(inun_proj_closest_mm) %>%
    select(depth_mm = inun_proj_closest_mm, inun_region) %>%
    mutate(depth_mm = round(depth_mm, 1)) %>%
    distinct(depth_mm, inun_region)
  
  download_tif_df <- relevant_ws_inun_proj %>%
    left_join(tif_df, by = c("inun_region", "depth_mm")) %>%
    drop_na(url)
  
  download_NOAA_inun_slr_tif_data(tif_url_df = download_tif_df,
                                  save_dir = here(output_dir, "Extent_Rasters"))
  
  ## Generate the water supply data with the inundation region information ----
  
  # Now, we will save the water supply data with the information that we gathered
  # from this exercise. Those new columns are:
  
  # `buffer_km`: The size of the buffer in km around the water supply.
  # `inun_proj_mm`: Projected SLR in millimeters from ITF projections.
  # `inun_proj_closest_mm`: Nearest available NOAA inundation extent scenario in millimeters that matches or exceeds the ITF projection.
  # `depth_ft`: Inundation depth in feet from the original NOAA dataset naming convention.
  # `depth_m`: Inundation depth converted to meters.
  # `inun_region`: NOAA inundation region identifier that spatially contains each water supply location.
  # `inun_surf_url`: Direct URL link to download the corresponding NOAA inundation surface raster file.
  # `proj_file`: Local base file name of the downloaded inundation surface raster for the projected scenario.
  # `ref_file`: Local base file name of the reference (0 feet SLR) inundation surface raster for baseline conditions.
  
  updated_ws_inun_proj <- ws_itf_slr_proj %>%
    inner_join(ws_inun_regions %>% st_drop_geometry() %>% select(ID, buffer_km, inun_region, inun_proj_closest_mm),
               by = c("ID", "buffer_km")) %>%
    left_join(download_tif_df %>% select(-state),
              by = c("inun_region", "inun_proj_closest_mm" = "depth_mm")) %>%
    left_join(reference_tif_df %>% select(inun_region, reference_file = file),
              by = "inun_region") %>%
    rename(inun_surf_url = url, proj_file = file, ref_file = reference_file) %>%
    relocate(depth_ft, depth_m, .after = inun_proj_closest_mm) %>%
    relocate(geometry, .after = last_col()) %>%
    filter(!is.na(ref_file), !is.na(proj_file))
  
  sf::st_write(updated_ws_inun_proj,
               here(output_dir, "buff_ws_inun_regions.gpkg"), # Buffered water supplies with inundation regions
               append = FALSE)
  
  # Environment clean up
  rm(ws_inun_regions, download_tif_df, reference_tif_df, relevant_ws_inun_proj,
     inun_region_bboxs, raster_f, tif_df, all_buffer_ws, final_itf_slr_proj, 
     ws_itf_slr_proj)
  gc()
  
  # Now we have the appropriate data to determine what the difference in inundation risk will be.
  
  return(updated_ws_inun_proj)
}