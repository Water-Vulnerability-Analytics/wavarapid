#' @title Flood sensitivity: area in high-risk FEMA flood zones
#'
#' @description
#' Calculates how much of the area around each water supply is in a
#' high-risk flood zone (FEMA zones starting with A or V, the 1% annual
#' chance floodplain). Each water supply is buffered by `dist` meters, and
#' the FEMA National Flood Hazard Layer is queried over the buffer's
#' bounding box, one request per supply, from
#' <https://hazards.fema.gov/arcgis/rest/services/public/NFHL/MapServer>. It
#' needs internet, and a national run takes hours.
#'
#' A supply gets:
#' * `NA` when the request fails or FEMA returns no flood mapping there;
#' * `0` when there is mapping but no high-risk zone in the buffer;
#' * otherwise, the area and percentage of the buffer in high-risk zones.
#'
#' If a `wsd_source_id` has several rows, their values are averaged.
#'
#' @param data sf object. The water supplies, as points with `wsd_source_id`,
#'   `park_unit` and `park_name` columns, e.g. from [wavacore::get_wsd()].
#'   Required.
#' @param dist Numeric. The buffer radius in meters, applied in UTM zone 16N
#'   (EPSG 32616).
#'   * A positive distance, e.g. `4000` as in the rapid workflow: the
#'     flood zones inside that circle are measured.
#'   * `0`: buffers the points to empty polygons, so the queries can't work
#'     and every supply gets `NA`.
#'
#'   Default: `0`.
#' @param save Logical. Whether to write the results to a CSV file.
#'   * `TRUE`: writes `flood_sensitivity_<date>.csv`, named with today's
#'     date, to `output_dir`, creating the folder if needed. Stops with an
#'     error if `output_dir` is `NULL`.
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
#'   * `TRUE`: shows a progress bar over the water supplies, a message for
#'     each failed request, and the run time and the path of the saved file.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return A tibble, one row per water supply, with `wsd_source_id`,
#'   `raw_flood_risk_area_m2` (area of high-risk zones in the buffer, square
#'   meters), `raw_flood_risk_percent` (percentage of the buffer),
#'   `park_unit`, `park_name`, and the `norm_` and `rank_` versions of the
#'   two `raw_` columns from [add_standardized_indicators()].
#' @export
#' @seealso [calc_vulnerability_index()], which uses the percentage as
#'   `sen_flood_risk`.
#' @examplesIf interactive()
#' # Queries the FEMA flood map service, so it needs internet
#' ws <- wavacore::get_wsd(park = "ZION")
#' flood_sen <- calc_flood_risk(ws, dist = 4000)
#' flood_sen

calc_flood_risk <-
  function(data,
           dist = 0,
           save = FALSE,
           output_dir = NULL,
           verbose = TRUE) {
    
    # Input validation
    stopifnot(
      "water_supplies must be an sf object" = inherits(data, "sf"),
      "save is TRUE but output_dir not provided" = !save ||
        !is.null(output_dir)
    )
    
    # buffer points to specified distance
    data <- data %>% 
      sf::st_transform(crs = 32616) %>% 
      sf::st_buffer(dist = dist)
    
    # match CRS to that of FEMA flood layer
    if (sf::st_crs(data) != sf::st_crs(4269)) {
      data <- sf::st_transform(data, crs = 4269)
    }
    
    # Start timing
    start_time <- Sys.time()
    
    # Create progress bar
    if (verbose) {
      pb <- progress::progress_bar$new(
        format = "Processing flood sensitivity [:bar] :percent (:current/:total) ETA: :eta",
        total = nrow(data),
        clear = FALSE,
        width = 80
      )
    }
    
    # for each water supply....
    results <-  purrr::map_dfr(seq_len(nrow(data)), function(i) {
      if (verbose) pb$tick()
    
      
      boundary <- data[i, ]
      
      ## get bounding box
      bb <- sf::st_bbox(boundary)
      
      ## extract bbox numbers
      bb.ordered <- paste(bb[1], bb[2], bb[3], bb[4], sep = "%2C")
      
      ## construct URL (OG version)
      url <- paste0(
        "https://hazards.fema.gov/arcgis/rest/services/public/NFHL/MapServer/",
        28,
        "/query?",
        "&geometry=",
        bb.ordered,
        "&geometryType=esriGeometryEnvelope",
        "&outFields=*",
        "&returnGeometry=true",
        "&returnZ=false",
        "&returnM=false",
        "&returnExtentOnly=false",
        "&f=geoJSON"
      )
      
      ## Claude version
      # url <- paste0(
      #   "https://hazards.fema.gov/arcgis/rest/services/public/NFHL/MapServer/28/query?",
      #   "geometry=", bb[1], ",", bb[2], ",", bb[3], ",", bb[4],
      #   "&geometryType=esriGeometryEnvelope",
      #   "&inSR=4326",  # Input spatial reference
      #   "&outSR=4326", # Output spatial reference
      #   "&spatialRel=esriSpatialRelIntersects",
      #   "&outFields=*",
      #   "&returnGeometry=true",
      #   "&f=geojson"
      # )
      
      ## read in floodplain
      # floodHaz <- sf::read_sf(url) %>%
      #   st_make_valid()
      
      ## read in floodplain with error handling
      floodHaz <- tryCatch(
        {
          sf::read_sf(url) %>%
            sf::st_make_valid()
        },
        error = function(e) {
          if (verbose) {
            message("\nWarning: Unable to read flood data for ",
                    data[i, ]$wsd_source_id)
          }
          # Return NULL to indicate error
          return(NULL)
        }
      )
      
      # if error occurred during read
      if (is.null(floodHaz)) {
        df <- dplyr::tibble(
          wsd_source_id = data[i, ]$wsd_source_id,
          park_unit = data[i, ]$park_unit,
          park_name = data[i, ]$park_name,
          raw_flood_risk_area_m2 = NA,
          raw_flood_risk_percent = NA
        )
        
      } else if (nrow(floodHaz) == 0) {
        df <- dplyr::tibble(
          wsd_source_id = data[i, ]$wsd_source_id,
          park_unit = data[i, ]$park_unit,
          park_name = data[i, ]$park_name,
          raw_flood_risk_area_m2 = NA,
          raw_flood_risk_percent = NA
        )
        
      } else {
        ## filter to zones that hve A or V in them - high-risk flood zones/1% flood prob https://floodpartners.com/flood-zones/
        floodRisk <- floodHaz %>%
          dplyr::filter(stringr::str_detect(FLD_ZONE, "A|V") &
                   FLD_ZONE != "AREA NOT INCLUDED") # not sure if 'area not included' is in fema raw data
        
        
        # if no high risk flood zones
        if (nrow(floodRisk) == 0) {
          df <- dplyr::tibble(
            wsd_source_id = data[i, ]$wsd_source_id,
            park_unit = data[i, ]$park_unit,
            park_name = data[i, ]$park_name,
            raw_flood_risk_area_m2 = 0,
            raw_flood_risk_percent = 0
          )
        } else {
          # match CRS
          boundary <- sf::st_transform(boundary, sf::st_crs(floodRisk))
          
          # Set attribute-geometry relationship to suppress warning
          sf::st_agr(floodRisk) <- "constant"
          sf::st_agr(boundary) <- "constant"
          
          # return overlapping areas
          floodArea <- sf::st_intersection(floodRisk, boundary)
          
          # if no intersection
          if (nrow(floodArea) == 0) {
            df <- dplyr::tibble(
              wsd_source_id = data[i, ]$wsd_source_id,
              park_unit = data[i, ]$park_unit,
              park_name = data[i, ]$park_name,
              raw_flood_risk_area_m2 = 0,
              raw_flood_risk_percent = 0
            )
          } else {
            df <- dplyr::tibble(
              wsd_source_id = data[i, ]$wsd_source_id,
              park_unit = data[i, ]$park_unit,
              park_name = data[i, ]$park_name,
              raw_flood_risk_area_m2 = as.numeric(sum(sf::st_area(
                floodArea
              ))),
              raw_flood_risk_percent = as.numeric(sum(sf::st_area(
                floodArea
              )) / sf::st_area(boundary) * 100)
            )
            
          }
          return(df)
        }
        
      }
    }
    )
    
    # Average across geometries if a wsd_source_id had multiple entries in water_supplies
    results <- results %>% 
      dplyr::summarise(
        raw_flood_risk_area_m2  = mean(raw_flood_risk_area_m2, na.rm = TRUE),
        raw_flood_risk_percent  = mean(raw_flood_risk_percent, na.rm = TRUE),
        park_unit = dplyr::first(park_unit),
        park_name = dplyr::first(park_name),
        .by = wsd_source_id
      )
    
    # Calculate standardized indicators using helper function
    final_results <- results %>%
      add_standardized_indicators("raw_flood_risk")
    
    # Save results
    if (save) {
      if (!dir.exists(output_dir))
        dir.create(output_dir, recursive = TRUE)
      
      filename <- file.path(output_dir,
                            paste0("flood_sensitivity_", Sys.Date(), ".csv"))
      readr::write_csv(final_results, filename)
      
      if (verbose) {
        message("Results saved to: ", filename)
      }
    }
    
    # Print timing
    if (verbose) {
      end_time <- Sys.time()
      total_time <- end_time - start_time
      message(
        "\nFlood sensitivity processing completed in ",
        round(total_time, 2),
        " ",
        attr(total_time, "units")
      )
    }
    
    return(final_results)
    
  }

