#' @title Groundwater sensitivity: depth to water
#'
#' @description
#' Gives each water supply a depth to water, in meters:
#' * groundwater and spring sources (`source_type` containing
#'   "groundwater" or "spring", lowercase): the WSD's recorded `well_depth`
#'   (feet) times 0.8, assuming the water sits 20% above the well bottom,
#'   capped at 200 m; where there is no well depth, the modeled water table
#'   depth from [process_groundwater_depth()];
#' * every other source type: 0.
#'
#' A shallower water table means higher sensitivity, so the depth is
#' inverted when standardized. The depth is also binned into a DRASTIC-style
#' score: 0-1.5 m = 1, 1.5-4.5 m = 0.9, 4.5-9.1 m = 0.7, 9.1-15.2 m = 0.5,
#' 15.2-22.9 m = 0.3, 22.9-30.5 m = 0.2 and over 30.5 m = 0.1 (5, 15, 30, 50,
#' 75 and 100 ft). The bin score is not inverted. Sources with a depth of 0
#' fall in the top bin.
#'
#' Citations: Zenodo record 18504963,
#' <https://ngwa.onlinelibrary.wiley.com/doi/10.1111/gwat.13362>; Aller, L.,
#' Lehr, J.H., Petty, R., and Bennett, T. (1987). DRASTIC: A standardized
#' system to evaluate groundwater pollution potential using hydrogeologic
#' setting. Journal Geological Society of India, 29(1), 23-37,
#' <https://doi.org/10.17491/jgsi/1987/290112>.
#'
#' @param water_supplies sf object. The water supplies, with
#'   `wsd_source_id`, `park_unit`, `park_name`, `source_type` and
#'   `well_depth` (feet) columns, e.g. from [wavacore::get_wsd()]. Stops with
#'   an error if `source_type` or `well_depth` is missing. Repeated
#'   `wsd_source_id`s are dropped. Required.
#' @param groundwater_data String or data frame. The modeled water table
#'   depths, with `wsd_source_id` and `wtd_point_m` columns.
#'   * String: path to `wtd_all_sites.csv` from [process_groundwater_depth()],
#'     e.g. `rapid/inputs/groundwater_depth/wtd_all_sites.csv` in the shared
#'     data folder.
#'   * Data frame: the same table, already read.
#'
#'   Required.
#' @param save Logical. Whether to write the results to a CSV file.
#'   * `TRUE`: writes `groundwater_depth_sensitivity_<date>.csv`, named with
#'     today's date, to `output_dir`, creating the folder if needed. Stops
#'     with an error if `output_dir` is `NULL`.
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
#'   * `TRUE`: prints messages with the run time and the path of the saved
#'     file.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return A tibble, one row per water supply, with `wsd_source_id`,
#'   `park_unit`, `park_name`, `raw_groundwater_depth_m`,
#'   `raw_groundwater_depth_source` (`"well_depth"`, `"modeled_wtd"` or
#'   `"surface_zero"`), `raw_groundwater_depth_bin`, and the `norm_` and
#'   `rank_` versions of the depth (inverted) and the bin score (not
#'   inverted) from [add_standardized_indicators()].
#' @export
#' @seealso [process_groundwater_depth()], which makes `wtd_all_sites.csv`;
#'   [calc_vulnerability_index()], which uses the depth as
#'   `sen_source_type`.
#' @examplesIf file.exists(wavacore::wava_data_path("rapid", "inputs", "groundwater_depth", "wtd_all_sites.csv"))
#' ws <- wavacore::get_wsd(park = "ZION")
#' gw_sen <- calc_groundwater_depth_sensitivity(
#'   ws,
#'   groundwater_data = wavacore::wava_data_path("rapid", "inputs", "groundwater_depth",
#'                                               "wtd_all_sites.csv")
#' )
#' gw_sen[, c("wsd_source_id", "raw_groundwater_depth_m",
#'            "raw_groundwater_depth_source", "raw_groundwater_depth_bin")]

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
