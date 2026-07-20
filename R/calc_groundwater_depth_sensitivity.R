#' Calculate Groundwater Depth Sensitivity Indicator
#'
#' Calculates groundwater depth sensitivity for groundwater-sourced water supplies.
#' The WSD's recorded `well_depth` is treated as ground-truth and used as the
#' primary depth-to-water estimate, assuming a 10 ft well screen (i.e. the water
#' table sits ~10 ft above the bottom of the well: `well_depth_ft - 10`). Where
#' `well_depth` is missing, this falls back to modeled water table depth (WTD)
#' from `process_groundwater_depth()` (CONUS-wide RF WTD model). A shallower
#' water table means higher sensitivity, so `invert = TRUE` is used when
#' standardizing the raw (continuous) depth.
#'
#' The raw depth is also binned into a DRASTIC-style depth-to-water sensitivity
#' score, using thresholds of 1.5/4.5/9.1/15.2/22.9/30.5 m (5/15/30/50/75/100 ft)
#' mapped to scores of 1/.9/.7/.5/.3/.2/.1. Unlike the raw depth, the bin score
#' is not inverted when standardized.
#'
#' Only applies to water supplies with a groundwater `source_type`
#' ("shallow groundwater", "deep groundwater", "unknown groundwater");
#' all other source types are excluded from the output.
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
#'   ("well_depth" or "modeled_wtd", indicating which source was used), and
#'   raw/norm/rank columns for `raw_groundwater_depth_m` (continuous depth,
#'   inverted) and `raw_groundwater_depth_bin` (DRASTIC-style bin score, not
#'   inverted)

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
    groundwater_data <- read_csv(groundwater_data, show_col_types = FALSE)
  }

  # Restrict to groundwater-sourced water supplies
  groundwater_supplies <- water_supplies %>%
    sf::st_drop_geometry() %>%
    dplyr::filter(grepl("groundwater", source_type)) %>%
    dplyr::mutate(well_depth = as.numeric(well_depth)) %>%
    dplyr::select(wsd_source_id, park_unit, park_name, well_depth) %>%
    dplyr::distinct(wsd_source_id, .keep_all = TRUE)

  # Join modeled water table depth, then prefer the WSD's recorded well depth
  # (ground truth) over the modeled WTD, assuming a 10 ft well screen so the
  # water table sits ~10 ft above the bottom of the well. Depths that would go
  # negative under that assumption (well shallower than 10 ft) are clipped to 0.
  results <- groundwater_supplies %>%
    dplyr::left_join(
      groundwater_data %>% dplyr::select(wsd_source_id, wtd_point_m),
      by = "wsd_source_id"
    ) %>%
    dplyr::mutate(
      well_depth_m = pmax(0, (well_depth - 10) * 0.3048),
      raw_groundwater_depth_m = dplyr::coalesce(well_depth_m, wtd_point_m),
      raw_groundwater_depth_source = dplyr::case_when(
        !is.na(well_depth_m) ~ "well_depth",
        !is.na(wtd_point_m)  ~ "modeled_wtd",
        .default = NA_character_
      )
    ) %>%
    dplyr::select(-well_depth, -well_depth_m, -wtd_point_m)

  # Bin the raw depth into a DRASTIC-style sensitivity score, then standardize
  # both the raw depth and the bin score separately
  final_results <- results %>%
    mutate(raw_groundwater_depth_bin = case_when(
      between(raw_groundwater_depth_m, 0, 1.5) ~ 1,
      between(raw_groundwater_depth_m, 1.5, 4.5) ~ .9,
      between(raw_groundwater_depth_m, 4.5, 9.1) ~ .7,
      between(raw_groundwater_depth_m, 9.1, 15.2) ~ .5,
      between(raw_groundwater_depth_m, 15.2, 22.9) ~ .3,
      between(raw_groundwater_depth_m, 22.9, 30.5) ~ .2,
      raw_groundwater_depth_m >= 30.5 ~ .1,
      .default = NA
    )) %>%
    add_standardized_indicators("raw_groundwater_depth_m", invert = TRUE) %>%
    add_standardized_indicators("raw_groundwater_depth_bin", invert = FALSE)

  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    filename <- file.path(output_dir, paste0("groundwater_depth_sensitivity_", Sys.Date(), ".csv"))
    write_csv(final_results, filename)

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
