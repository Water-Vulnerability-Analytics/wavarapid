#' Calculate Groundwater Depth Sensitivity Indicator
#'
#' Calculates depth-to-water sensitivity for all water supplies. For
#' groundwater/spring `source_type` supplies, depth-to-water is taken from the
#' WSD's recorded `well_depth`, assumed to sit 20% shallower than the well
#' bottom (i.e. `well_depth * 0.8`), falling back to modeled water table depth
#' (WTD) from `process_groundwater_depth()` (CONUS-wide RF WTD model) where
#' `well_depth` is missing. All other (non-groundwater/spring) source types
#' are hardcoded to a depth-to-water of 0. A shallower water table means higher
#' sensitivity, so `invert = TRUE` is used when standardizing the raw
#' (continuous) depth.
#'
#' The raw depth is also binned into a DRASTIC-style depth-to-water sensitivity
#' score, using thresholds of 1.5/4.5/9.1/15.2/22.9/30.5 m (5/15/30/50/75/100 ft)
#' mapped to scores of 1/.9/.7/.5/.3/.2/.1. Unlike the raw depth, the bin score
#' is not inverted when standardized. A depth of 0 (non-groundwater/spring
#' source types) naturally falls into the top bin (1) under these thresholds.
#'
#' Citation: Zenodo record 18504963,
#' https://ngwa.onlinelibrary.wiley.com/doi/10.1111/gwat.13362
#'
#' Citation: Aller, Linda, Jay H. Lehr, Rebecca Petty, and Truman Bennett.
#' "DRASTIC: A Standardized System to Evaluate Groundwater Pollution Potential
#' Using Hydrogeologic Setting." Journal Geological Society of India 29, no. 1
#' (1987): 23-37. https://doi.org/10.17491/jgsi/1987/290112.
#'
#' @param water_supplies sf object with water supply locations, including
#'   `wsd_source_id`, `park_unit`, `park_name`, `source_type`, `well_depth`
#'   (recorded well depth in feet).
#' @param groundwater_data character path to `wtd_all_sites.csv` (output of
#'   `process_groundwater_depth()`), or a data.frame with the same columns.
#' @param save logical, whether to save output
#' @param output_dir character, directory path for saving files
#' @param verbose logical, whether to show progress messages
#'
#' @return tibble with water supply ID, park info, `raw_groundwater_depth_source`
#'   ("well_depth", "modeled_wtd", or "surface_zero", indicating which
#'   source was used), and raw/norm/rank columns for `raw_groundwater_depth_m`
#'   (continuous depth, inverted) and `raw_groundwater_depth_bin` (DRASTIC-style
#'   bin score, not inverted)

#' @export
calc_groundwater_depth_sensitivity <- function(water_supplies,
                                               groundwater_data,
                                               save = FALSE,
                                               output_dir = NULL,
                                               verbose = TRUE) {

  # Input validation
  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "water_supplies must include source_type" = "source_type" %in% names(water_supplies),
    "water_supplies must include well_depth" = "well_depth" %in% names(water_supplies),
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir)
  )

  # Start timing
  if (verbose) {
    message("Starting groundwater depth sensitivity calculation...")
    start_time <- Sys.time()
  }

  # Load pre-processed modeled WTD if a file path was given
  if (is.character(groundwater_data)) {
    groundwater_data <- readr::read_csv(groundwater_data, show_col_types = FALSE)
  }

  # Score all water supplies. Groundwater/spring source types get a depth
  # from the well-depth/modeled-WTD priority chain below; all other
  # (non-groundwater/spring) source types are hardcoded to a depth of 0.
  
  all_supplies <- water_supplies %>%
    sf::st_drop_geometry() %>%
    dplyr::select(wsd_source_id, park_unit, park_name, source_type,
                  well_depth) %>%
    dplyr::mutate(well_depth = as.numeric(well_depth)) %>%
    dplyr::distinct(wsd_source_id, .keep_all = TRUE)

  # Join modeled water table depth, then prefer the WSD's recorded well depth
  # (ground truth), assumed to sit 20% shallower than the well bottom, over
  # the modeled WTD.

  results <- all_supplies %>%
    dplyr::left_join(
      groundwater_data %>% dplyr::select(wsd_source_id, wtd_point_m),
      by = "wsd_source_id"
    ) %>%
    dplyr::mutate(
      is_gw_spring = grepl("groundwater|spring", source_type),
      well_depth_m = dplyr::if_else((well_depth * 0.8 * 0.3048) >= 200, 200, (well_depth * 0.8 * 0.3048)),
      raw_groundwater_depth_m = dplyr::case_when(
        !is_gw_spring              ~ 0,
        !is.na(well_depth_m)       ~ well_depth_m,
        !is.na(wtd_point_m)        ~ wtd_point_m,
        .default = NA_real_
      ),
      raw_groundwater_depth_source = dplyr::case_when(
        !is_gw_spring              ~ "surface_zero",
        !is.na(well_depth_m)       ~ "well_depth",
        !is.na(wtd_point_m)        ~ "modeled_wtd",
        .default = NA_character_
      )
    ) %>%
    dplyr::select(-source_type, -well_depth, -is_gw_spring,
                  -well_depth_m, -wtd_point_m)

  # Bin the raw depth into a DRASTIC-style sensitivity score, then standardize
  # both the raw depth and the bin score separately
  final_results <- results %>%
    dplyr::mutate(raw_groundwater_depth_bin = dplyr::case_when(
      dplyr::between(raw_groundwater_depth_m, 0, 1.5) ~ 1,
      dplyr::between(raw_groundwater_depth_m, 1.5, 4.5) ~ .9,
      dplyr::between(raw_groundwater_depth_m, 4.5, 9.1) ~ .7,
      dplyr::between(raw_groundwater_depth_m, 9.1, 15.2) ~ .5,
      dplyr::between(raw_groundwater_depth_m, 15.2, 22.9) ~ .3,
      dplyr::between(raw_groundwater_depth_m, 22.9, 30.5) ~ .2,
      raw_groundwater_depth_m >= 30.5 ~ .1,
      .default = NA
    )) %>%
    add_standardized_indicators("raw_groundwater_depth_m", invert = TRUE) %>%
    add_standardized_indicators("raw_groundwater_depth_bin", invert = FALSE)

  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    filename <- file.path(output_dir, paste0("groundwater_depth_sensitivity_", Sys.Date(), ".csv"))
    readr::write_csv(final_results, filename)

    if (verbose) {
      message("Results saved to: ", filename)
    }
  }

  # Print timing
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message("\nGroundwater depth sensitivity processing completed in ",
            round(total_time, 2), " ", attr(total_time, "units"))
  }

  return(final_results)
}
