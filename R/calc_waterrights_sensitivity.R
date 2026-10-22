#' @title Water rights sensitivity: current water stress (Aqueduct 4.0)
#'
#' @description
#' Gives each water supply the baseline (1979-2019) Aqueduct 4.0 water stress
#' score of the HydroBASINS level 6 sub-basin it is in, from the table made
#' by [process_waterrightsindicator_data()].
#'
#' Water stress is total withdrawals (domestic, industrial, irrigation,
#' livestock) divided by the renewable water available, so a high score
#' means many users already compete for the same water. That is when water
#' rights administration starts to bind: where water is plentiful a park's
#' right is rarely tested, and where competition is high the state's doctrine
#' decides whether it can divert its full entitlement.
#'
#' The score (0-5) is used rather than the raw ratio: it is a log-normalized
#' version of the same signal, has published category thresholds that apply
#' to both baseline and future, and is better distributed than the raw ratio
#' (unbounded, right-skewed, and 9999 for arid basins). A score of 5 covers
#' both "Extremely High (>80%)" and "Arid and Low Water Use", which Aqueduct
#' treats as the top priority.
#'
#' Citation: Kuzma, S. et al. (2023) Aqueduct 4.0: Updated decision-relevant
#' global water risk indicators. World Resources Institute.
#'
#' @param water_supplies sf object. The water supplies, with a
#'   `wsd_source_id` column. Only used as a check: the function stops with an
#'   error if the `wsd_source_id`s in `stress_data` are not exactly the same
#'   set, which catches a stress table made from a different or older list
#'   of supplies. Required.
#' @param stress_data String or data frame. The long water stress table from
#'   [process_waterrightsindicator_data()], with columns `wsd_source_id`,
#'   `park_unit`, `park_name`, `period`, `stress_score` and `stress_label`.
#'   Only the `period == "baseline"` rows are used.
#'   * String: path to `waterrights_stress_long.csv`. Stops with an error if
#'     the file doesn't exist.
#'   * Data frame: the same table, already read.
#'
#'   Default: `rapid/inputs/WaterRightsIndicator/waterrights_stress_long.csv`
#'   in the shared data folder, see [wavacore::wava_data_path()].
#' @param save Logical. Whether to write the results to a CSV file.
#'   * `TRUE`: writes `waterrights_sensitivity_<date>.csv`, named with today's date, to
#'     `output_dir`, creating the folder if needed. Stops with an error if
#'     `output_dir` is `NULL`.
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
#'   * `TRUE`: prints each step, the run time and the path of the saved file.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return A tibble, one row per water supply, with `wsd_source_id`,
#'   `park_unit`, `park_name`, `water_stress_label` (the Aqueduct category),
#'   `raw_water_stress` (baseline score, 0-5) and its `norm_` and `rank_`
#'   versions from [add_standardized_indicators()]. Higher stress gets a
#'   higher score.
#' @export
#' @seealso [calc_waterrights_exposure()] for the projected stress;
#'   [process_waterrightsindicator_data()], which makes `stress_data`.
#' @examplesIf file.exists(wavacore::wava_data_path("rapid", "inputs", "WaterRightsIndicator", "waterrights_stress_long.csv"))
#' # The stress table and water_supplies must cover the same supplies, so
#' # keep the same ones in both
#' stress <- readr::read_csv(wavacore::wava_data_path("rapid", "inputs",
#'                                                    "WaterRightsIndicator",
#'                                                    "waterrights_stress_long.csv"),
#'                           show_col_types = FALSE)
#' ws <- wavacore::get_wsd(park = "ZION")
#' ws <- ws[ws$wsd_source_id %in% stress$wsd_source_id, ]
#' stress <- stress[stress$wsd_source_id %in% ws$wsd_source_id, ]
#'
#' wr_sen <- calc_waterrights_sensitivity(ws, stress_data = stress)
#' wr_sen

calc_waterrights_sensitivity <- function(
    water_supplies,
    stress_data = wavacore::wava_data_path("rapid",
                       "inputs",
                       "WaterRightsIndicator",
                       "waterrights_stress_long.csv"),
    save = FALSE,
    output_dir = NULL,
    verbose = TRUE) {

  # Input validation
  stopifnot(
    "water_supplies must be an sf object"         = inherits(water_supplies, "sf"),
    "stress_data must be a file path or data.frame" =
      is.character(stress_data) || is.data.frame(stress_data),
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir)
  )
  
  # Start timing
  if (verbose) {
    message("Starting water rights sensitivity calculation...")
    start_time <- Sys.time()
  }

  # Load the long stress table if a file path was given
  if (is.character(stress_data)) {
    stopifnot(
      "stress_data path does not exist \u2014 run process_waterrightsindicator_data() first" =
        file.exists(stress_data)
    )
    stress_data <- readr::read_csv(stress_data, show_col_types = FALSE)
  }

  required_columns <- c("wsd_source_id", "park_unit", "park_name",
                        "period", "stress_score", "stress_label")

  stopifnot(
    "stress_data missing required columns" = all(required_columns %in% names(stress_data))
  )
  
  # stress_data must cover exactly the same water supplies that were passed in
  ws_ids     <- unique(water_supplies$wsd_source_id)
  stress_ids <- unique(stress_data$wsd_source_id)

  stopifnot(
    "wsd_source_id in stress_data does not match water_supplies \u2014 rerun process_waterrightsindicator_data() with this water_supplies object" =
      setequal(ws_ids, stress_ids)
  )
  
  if (verbose) {
    message("  ", length(stress_ids), " water supplies matched between stress_data and water_supplies")
  }
  
  # Sensitivity is the baseline period only — one row per water supply
  results <- stress_data %>%
    dplyr::filter(period == "baseline") %>%
    dplyr::distinct(wsd_source_id, .keep_all = TRUE) %>%
    dplyr::select(wsd_source_id, park_unit, park_name,
           raw_water_stress   = stress_score,
           water_stress_label = stress_label)

  if (nrow(results) == 0) {
    stop("No baseline rows found in stress_data. Check that period == 'baseline' exists.")
  }

  if (verbose) {
    message("  ", nrow(results), " water supplies with a baseline water stress score")
  }

  # invert = FALSE: higher water stress = higher sensitivity
  final_results <- results %>%
    add_standardized_indicators("raw_water_stress", invert = FALSE) %>%
    # reorder columns following standardized pattern
    dplyr::select(wsd_source_id, park_unit, park_name, dplyr::everything())

  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    filename <- file.path(output_dir,
                          paste0("waterrights_sensitivity_", Sys.Date(), ".csv"))
    readr::write_csv(final_results, filename)

    if (verbose) {
      message("Results saved to: ", filename)
    }
  }

  # Print timing
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message("\nWater rights sensitivity processing completed in ",
            round(total_time, 2), " ", attr(total_time, "units"))
  }

  return(final_results)
}
