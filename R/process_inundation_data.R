#' Process coastal inundation inputs for water-supply risk analysis
#'
#' Buffers water supply point locations (1 km and 4 km), extracts Interagency Task
#' Force (ITF) sea-level-rise projections for 2050 (83rd percentile), builds a master
#' NOAA Sea Level Rise Viewer (SLRV) inundation-extent URL table from a compiled text list,
#' determines the relevant NOAA inundation region per water supply (by state), selects
#' the nearest available NOAA inundation extent depth that matches/exceeds the ITF projection,
#' downloads the needed NOAA extent rasters, and writes intermediate and final geopackages/rasters.
#'
#' @param water_supplies An `sf` object of point geometries representing water supply locations.
#'   Must include a `state` column (2-letter state code) used to group downloads/processing.
#' @param coast_dir_path Character path to the coastal comparison data directory that contains
#'   ITF NetCDF inputs and NOAA SLRV URL list inputs, and where outputs will be written.
#'   Defaults to `"data/all/coast_comp_data"`.

process_inundation_data <- function(water_supplies,
                                    coast_dir_path = "data/all/coast_comp_data"){
  
  ## Generate the buffered water supply dataset
  
  # For this dataset we will buffer the water supply points by 1 and 4 km. We will
  # use this dataset to do our inundation risk analysis.
  
  transformed_ws <- water_supplies %>%
    st_transform(5070) # NAD83 / Conus Albers (meters)
  
  buffer_1km_ws <- transformed_ws %>%
    st_buffer(1000) %>% # 1 km
    st_transform(4326) %>% # back to 4326
    rowid_to_column(var = "ID") %>%
    mutate(buffer_km = 1)
  
  buffer_4km_ws <- transformed_ws %>%
    st_buffer(4000) %>% # 4 km
    st_transform(4326) %>% # back to 4326
    rowid_to_column(var = "ID") %>%
    mutate(buffer_km = 4)
  
  all_buffer_ws <- bind_rows(buffer_1km_ws, buffer_4km_ws)
  
  # Save the buffered water supplies
  sf::st_write(all_buffer_ws, 
               here(coast_dir_path, "altered_water_supplies", 
                    "buff_ws.gpkg"), # Buffered water supplies
               append = FALSE)
  
  rm(transformed_ws, buffer_1km_ws, buffer_4km_ws)
  
  gc()
  
  ## Generate the Interagency Task Force (ITF) Sea Level Rise (SLR) Projections
  
  # These data were downloaded from:
  #   [Global and Regional Sea Level Rise Scenarios for the United States: Updated Mean Projections and Extreme Water Level Probabilities Along U.S. Coastlines](https://zenodo.org/records/5951626)
  # 
  # The report that this data came from can be found here:
  #   [2022 Technical Report](https://earth.gov/sealevel/us/resources/2022-sea-level-rise-technical-report/)
  # 
  # An interactive tool displaying these data can be found here:
  #   [Interagency Sea Level Rise Scenario Tool](https://sealevel.nasa.gov/task-force-scenario-tool)
  # 
  # These data were made by the [Interagency Sea Level Task Force](https://earth.gov/sealevel/us/about/interagency-sea-level-task-force/). This working group
  # consists of members from NOAA, NASA, JPL, USGS, EPA, FEMA, DoD, and USACE.
  # 
  # These data are the projections (scenarios) of SLR from the interagency report.
  
  ### Generate the complete gridded ITF SLR projections
  
  inundation_projections_path <- here(coast_dir_path, "interagency_report", 
                                      "Results", "TR_gridded_projections.nc")
  
  # read in nc file as a raster with terra
  projections_raster <- terra::rast(inundation_projections_path)
  
  # assign projection manually 
  terra::crs(projections_raster) <- "+proj=longlat +datum=WGS84 +no_defs"
  
  # pick the projection that we want from this raster
  ## 83rd percentile projections out to the year 2050
  projection_layer <- projections_raster[["rsl_total_IntHigh_percentiles=83_years=2050"]]
  
  # rotate raster from 0–360 to -180–180
  proj_layer_rotated <- terra::rotate(projection_layer)
  
  # Convert data that essentially represent no change to NA values
  # This is the final manipulation before the projections get used for the analysis.
  final_itf_slr_proj <- terra::ifel(proj_layer_rotated <= -20000, NA, proj_layer_rotated)
  
  # Save the gridded projections as a tif
  terra::writeRaster(final_itf_slr_proj, 
                     here(coast_dir_path, "interagency_report", 
                          "generated_products", "itf_slr_projections_2050_83rd.tif"),
                     overwrite = TRUE)
  
  rm(inundation_projections_path, projections_raster, projection_layer, 
     proj_layer_rotated)
  
  gc()
  
  
  # ### Generate US Coast subset of ITF SLR projections
  # 
  # This version of the ITF SLR projection dataset will not be used for the analysis.
  # This version of the projection dataset is for visualization purposes.
  # 
  # Excludes Alaska.
  
  
  # Define CONUS bounding box
  conus_bbox <- sf::st_bbox(c(xmin = -125, ymin = 20, xmax = -66, ymax = 50))
  conus_bbox <- sf::st_as_sfc(conus_bbox)
  sf::st_crs(conus_bbox) <- 4269  # Match the coastline CRS
  
  # Get coastline and crop with the bounding box
  us_coastline <- tigris::coastline(year = 2020)
  conus_coastline <- sf::st_crop(us_coastline, conus_bbox)
  conus_coastline <- sf::st_transform(conus_coastline, crs = 4326)
  
  # Mask the us coastline on the raster data. Masking on the point buffers was inconsistent
  projection_coast_mask <- terra::mask(final_itf_slr_proj, conus_coastline, touches = TRUE)
  
  # Create a binary mask of where data exists
  data_mask <- !is.na(projection_coast_mask)
  
  # Expand the mask using focal 
  expanded_mask <- terra::focal(data_mask, w = 3, fun = max, na.rm = TRUE)
  
  # Apply the expanded mask to the original data
  projections_expanded <- terra::mask(final_itf_slr_proj, expanded_mask, maskvalue = 0, inverse = FALSE)
  
  relevant_projections <- terra::trim(projections_expanded)
  
  # Save the coastal gridded projections to a tif
  terra::writeRaster(relevant_projections,
                     here(coast_dir_path,
                          "interagency_report",  "generated_products", 
                          "itf_slr_proj_coast_subset.tif"),
                     overwrite = TRUE)  
  
  rm(us_coastline, conus_bbox, projection_coast_mask,
     data_mask, expanded_mask, projections_expanded, relevant_projections,
     conus_coastline)
  
  gc()
  
  
  ## Generate the water supply data with ITF SLR projections
  # 
  # Now that we have the buffered water supply points and the ITF SLR gridded projections
  # we can join them together to find the SLR projections for each of our water supply points.
  
  
  # Extract the inundation projections specific to each water supply
  water_supply_inundation_values <- terra::extract(final_itf_slr_proj, all_buffer_ws) %>%
    group_by(ID) %>%
    summarise(
      across(everything(), ~ if(all(is.na(.x))) NA_real_ else max(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  # Join water supply inundation projections to the water supply buffer data frame
  ws_itf_slr_proj <- left_join(all_buffer_ws, 
                               water_supply_inundation_values, 
                               by = "ID") %>%
    rename(inun_proj_mm = `rsl_total_IntHigh_percentiles=83_years=2050`)
  
  # Save the dataset
  sf::st_write(ws_itf_slr_proj,
               here(coast_dir_path, "altered_water_supplies",
                    "buff_ws_itf_slr_proj.gpkg"), # Buffered water supplies with ITF SLR projections
               append = FALSE)
  
  rm(water_supply_inundation_values)
  
  gc()
  
  
  ## Download the NOAA inundation surface datasets
  
  # We will be using the data used in [NOAA's Sea Level Rise Viewer](https://coast.noaa.gov/slr/#/layer/slr/0/-11581024.663779823/5095888.569004184/4/satellite/none/0.8/2050/interHigh/noAccretion/NOS_Minor)
  # as our inundation surface datasets. These data are not stored in a way that ties them
  # to any specific scenario of SLR, instead they are the inundation surfaces given
  # some amount of SLR. This means that we will need to use the information from the
  # ITF SLR projections in order to find the correct inundation surface given the amount
  # of SLR the ITF projected for some region specific to each watershed. 
  # 
  # The process for generating the NOAA inundation surface datasets is not simple.
  # This is because the inundation surface datasets are large and stored in
  # [Sea Level Rise Viewer Data Download](https://coast.noaa.gov/slrdata/). 
  # 
  # There is no clear way to easily programmatically download data from this source.
  # They provide text data with URL links for the raster data, but that text data is also not
  # readily accessible via an API. Here is an [example](https://coast.noaa.gov/slrdata/Extent_Rasters/CA/URLlist_CA.txt)
  # of what one of these text files looks like. The process to programmatically download this data is,
  # broadly: compile these text files into one master text file, convert this master
  # text file into a form that is easily stored as a dataframe, use the information
  # in this dataframe to download the appropriate data via the provided URLs.
  # 
  # We will use this section to download the appropriate NOAA SLRV inundation surface datasets.
  # 
  # ### Generate the master inundation dataset URL link dataframe
  # 
  # Here we convert the master inundation dataset URL text file into the master 
  # inundation dataset URL dataframe. The inundation surface data is stored as feet
  # so we will need to convert it to millimeters to join it to the ITF SLR projections, 
  # which are themselves in mm.
  # 
  # Note that the creation of the master inundation dataset URL text file was made outside
  # of R. I went into each state's download repository to copy the URL text file 
  # for each state into the master text file. 
  
  
  # File path for the compiled master inundation dataset URL text file
  tif_txt_path <- here(coast_dir_path, "NOAA_SLRV_2022", "ALL_extent_raster_URLlist.txt")
  
  # Read in the text file as a vector
  tif_vec <- read_lines(tif_txt_path)
  
  tif_df <- tibble(url = tif_vec) %>% # Initialize the master URL df with the URL list
    mutate(
      file = basename(url), 
      info = str_split(str_remove(file, ".tif"), "_"),
      state = map_chr(info, ~ .x[1]),
      depth_ft_in = map(info, ~ {
        n <- length(.x)
        if (n >= 2) {
          last_two <- .x[(n-1):n]
          is_numeric <- str_detect(last_two, "^\\d+(\\.\\d+)?$")
          
          if (all(is_numeric)) {
            # Both are numeric - return both
            # Combine them with a period in the middle
            depth <- str_flatten(last_two, collapse = ".")
            return(as.numeric(depth))
          } else if (is_numeric[2]) {
            # Only the last one is numeric - return just that
            return(as.numeric(.x[n]))
          }
        } else if (n >= 1 && str_detect(.x[n], "^\\d+(\\.\\d+)?$")) {
          # Only one element and it's numeric
          return(as.numeric(.x[n]))
        }
        return(NA_real_)
      })
    ) %>% 
    select(-info) %>% 
    filter(!is.na(depth_ft_in)) %>% 
    mutate(
      depth_ft_in = as.numeric(depth_ft_in),
      depth_m = depth_ft_in * 0.3048,
      depth_mm = depth_m * 1000,
      depth_mm = round(depth_mm, 1),
      inun_region = str_remove(file, "_[0-9_]+\\.tif$")
    ) %>%
    select(
      url, file, 
      state, inun_region,
      depth_ft_in, depth_m, depth_mm
    ) %>%
    filter(!state %in% c("AK", "Guam", "HI", "CNMI", "USVI"))
  
  write_csv(tif_df, here(coast_dir_path, "NOAA_SLRV_2022", "ALL_extent_raster_URLlist.csv"))
  
  
  ### Download the reference inundation surface datasets
  
  # To determine the difference between the present day sea level conditions and
  # the conditions in our scenarios we need to download the surface inundation datasets
  # which predict 0 feet of SLR. We will use these datasets as our reference datasets.
  # We will download all of these since we do not yet know which inundation surface regions 
  # will be relevant to us. We will use these reference datasets to figure that out.
  # 
  # To do this we will use `download_NOAA_inun_surf_tif_data()`:
  
  # DOWNLOADING NOAA SURFACE INUNDATION TIF DATASETS HELPER FUNCTION =============
  
  #' Download NOAA SLRV inundation surface tif data from master URL df
  #' 
  #' @param slr_tif_url_df df which holds the URLs for inundation surface raster tifs from NOAA SLRV
  
  download_NOAA_inun_surf_tif_data <- function(tif_url_df) {
    tif_url_df %>% 
      pwalk(
        function(url, file, state, ...){
          
          # Create state directory if it doesn't exist
          dest_dir <- here("data", "all", "coast_comp_data", "NOAA_SLRV_2022", "Extent_Rasters", state)
          dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
          
          # Try to download data 3 times
          for (i in 1:3){
            tryCatch({
              download.file(
                url = url,
                destfile = here("data", "all", "coast_comp_data", "NOAA_SLRV_2022", "Extent_Rasters", state, file), 
                mode = "wb"
              ) 
              cat("\nDownload succeeded (attempt:", i,") for:", file, "\n")
              break # exit the for loop on success
            }, error = function(e) {
              cat("\nDownload failed (attempt:", i,") for:", file, "\n")
              Sys.sleep(i)
              if (i == 3) cat("\n⚠️ All 3 attempts failed for:", file, "\n")
            })
          }
          
          gc()
          Sys.sleep(1)
          
          # All 3 download attempts failed
          return()
          
        }
      )
  }
  
  
  
  reference_tif_df <- tif_df %>%
    filter(depth_mm == 0.0)
  
  download_NOAA_inun_surf_tif_data(reference_tif_df)
  
  
  ### Download the scenario inundation surface datasets
  
  # First, we will use the reference inundation surface datasets to determine 
  # the specific inundation region and SLR scenario for each of our water supplies.
  # 
  # Since it takes so long to pull in the data and RStudio cannot hold all of the
  # inundation surface rasters in memory at the same time we will do this part of the
  # analysis by state. We will load in the relevant data by state, determine the 
  # appropriate inundation region for each water supply, and then remove the raster
  # data from memory before moving on to the next state.
  
  
  # Make a list of relevant state codes
  rel_state_codes <- list.dirs(
    path = here(coast_dir_path, "NOAA_SLRV_2022", "Extent_Rasters"),
    full.names = FALSE,
    recursive = FALSE
  )
  
  # Water supplies that are relevant to this analysis split by state
  inun_rel_ws_list <- ws_itf_slr_proj %>%
    filter(
      state %in% rel_state_codes,
      !is.na(inun_proj_mm)
    ) %>%
    split(f = .$state)
  
  ws_inun_region_list <- imap(
    inun_rel_ws_list,
    function(df, state) {
      
      # Find the relevant state's tif files
      state_tif_files_dir <- here(coast_dir_path, 
                                  "NOAA_SLRV_2022", "Extent_Rasters", 
                                  state)
      
      state_tif_file_paths <- list.files(path = state_tif_files_dir,
                                         # pattern excludes non-zero digits in the mm section of the file name
                                         pattern = ".*[A-Za-z0_-]{4}\\.tif$",
                                         full.names = TRUE)
      
      # Read in the tif files and make a dataframe that contains each regions' bbox
      inun_rast_bboxs <- map_dfr(
        state_tif_file_paths,
        function(tif_f) {
          
          # Get the name of the region
          region_name <- tools::file_path_sans_ext(basename(tif_f))
          region_name <- sub("connectRaster.*", "connectRaster", region_name)
          
          # Read in the raster
          inun_raster <- terra::rast(here(tif_f))
          crs(inun_raster) <- "+proj=longlat +datum=WGS84 +no_defs"
          
          # Convert raster to bbox to mask with water supply inundation projection data
          bbox <- terra::ext(inun_raster)
          bbox_poly <- terra::as.polygons(bbox)
          terra::crs(bbox_poly) <- "+proj=longlat +datum=WGS84 +no_defs"
          bbox_sf <- sf::st_as_sf(bbox_poly)
          sf::st_crs(bbox_sf) <- 4326
          
          # Add the region to the bbox
          bbox_sf <- bbox_sf %>%
            mutate(inun_region = region_name)
          
          return(bbox_sf)
        }
      ) 
      
      # Spatial join to assign each water supply to a region
      water_supply_region_mapping <- sf::st_join(
        df, inun_rast_bboxs,
        join = st_intersects,
        largest = TRUE # If multiple matches, take largest overlap
      ) %>%
        # If the water supply does not fall in an inundation projection region,
        # it will have no inundation risk
        mutate(
          inun_proj_mm = ifelse(is.na(inun_region), NA, inun_proj_mm),
          # Find the closest greatest value available in the NOAA surface projections
          inun_proj_closest_mm = case_when(
            inun_proj_mm <= 0.0 ~ 0.0,
            inun_proj_mm <= 152.4 ~ 152.4,
            inun_proj_mm <= 304.8 ~ 304.8,
            inun_proj_mm <= 457.2 ~ 457.2,
            inun_proj_mm <= 609.6 ~ 609.6,
            inun_proj_mm <= 762.0 ~ 762.0,
            inun_proj_mm <= 914.4 ~ 914.4,
            inun_proj_mm <= 1066.8 ~ 1066.8,
            inun_proj_mm <= 1219.2 ~ 1219.2,
            inun_proj_mm <= 1371.6 ~ 1371.6,
            inun_proj_mm <= 3048.0 ~ 3048.0,
            TRUE ~ NA_real_
          )
        ) %>%
        relocate(inun_proj_closest_mm, .after = inun_proj_mm) %>%
        sf::st_drop_geometry(.)
      
      return(water_supply_region_mapping)
      
    })
  
  
  # Second, we will use the region and SLR scenario information to download only the
  # relevant inundation surface datasets from NOAA. 
  
  
  relevant_ws_inun_proj <- bind_rows(ws_inun_region_list) %>%
    group_by(state, inun_region) %>%
    slice_max(inun_proj_closest_mm, n = 1, with_ties = FALSE) %>%
    drop_na(inun_proj_closest_mm) %>%
    ungroup() %>%
    select(depth_mm = inun_proj_closest_mm, inun_region) %>%
    mutate(
      # CA EKA region is only available in 1 foot increments, so we specify that here
      depth_mm = ifelse(inun_region == "CA_EKA_connectRaster", 609.6, depth_mm), 
      depth_mm = round(depth_mm, 1)
    )
  
  download_tif_df <- relevant_ws_inun_proj %>%
    left_join(tif_df, by = c("inun_region", "depth_mm")) 
  
  download_NOAA_inun_surf_tif_data(download_tif_df) 
  
  
  
  ### Generate the water supply data with the inundation region information
  
  # Now, we will save the water supply data with the information that we gathered
  # from this exercise. Those new columns are:
  #   
  # - `buffer_km`: The size of the buffer in km around the water supply.
  # - `inun_proj_mm`: Projected SLR in millimeters from ITF projections.
  # - `inun_proj_closest_mm`: Nearest available NOAA inundation extent scenario in millimeters that matches or exceeds the ITF projection.
  # - `depth_ft`: Inundation depth in feet from the original NOAA dataset naming convention.
  # - `depth_m`: Inundation depth converted to meters.
  # - `inun_region`: NOAA inundation region identifier that spatially contains each water supply location.
  # - `inun_surf_url`: Direct URL link to download the corresponding NOAA inundation surface raster file.
  # - `proj_file`: Local base file name of the downloaded inundation surface raster for the projected scenario.
  # - `ref_file`: Local base file name of the reference (0 feet SLR) inundation surface raster for baseline conditions.
  
  
  updated_ws_inun_proj <- ws_itf_slr_proj %>%
    left_join(bind_rows(ws_inun_region_list)) %>%
    left_join(download_tif_df, 
              by = c("state", "inun_region", "inun_proj_closest_mm" = "depth_mm")) %>%
    left_join(select(reference_tif_df, 
                     state, reference_file = file, inun_region), 
              by = c("inun_region", "state")) %>%
    relocate(geometry, .after = last_col()) %>%
    relocate(depth_ft_in, depth_m, .after = inun_proj_closest_mm) %>%
    rename(
      "depth_ft" = "depth_ft_in",
      "inun_surf_url" = "url", 
      "proj_file" = "file", 
      "ref_file" = "reference_file"
    )
  
  sf::st_write(updated_ws_inun_proj,
               here(coast_dir_path, "NOAA_SLRV_2022", "inundation_data", 
                    "ws_inun_regions.gpkg"), # Water supplies with inundation regions
               append = FALSE)
  
  
  
  # Environment clean up
  rm(ws_inun_region_list, download_tif_df, reference_tif_df, relevant_ws_inun_proj,
     inun_rel_ws_list, rel_state_codes, tif_df, all_buffer_ws, final_itf_slr_proj, 
     ws_itf_slr_proj)
  
  gc()
  # Now we have the appropriate data to determine what the difference in inundation
  # risk will be.  
  return(updated_ws_inun_proj)
  
}
