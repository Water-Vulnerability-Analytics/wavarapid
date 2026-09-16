#' Calculate Water Rights Sensitivity Indicator
#'
#' Calculates water rights sensitivity from Aqueduct 4.0 baseline (1979-2019)
#' water stress for the HydroBASINS Level 6 sub-basin containing each NPS water
#' supply. Requires the long stress table produced by
#' `process_waterrightsindicator_data()`.
#'
#' Sensitivity is how constrained a water system already is. Water stress is the
#' ratio of total withdrawals (domestic, industrial, irrigation, livestock) to
#' the renewable blue water available in the sub-basin, so a high score means
#' many users are already competing for the same water. Competition is the
#' condition under which water rights administration becomes binding: where
#' supply is ample a park's legal right is rarely tested, and where competition
#' is high the governing doctrine determines whether it can divert its full
#' entitlement.
#'
#' The Aqueduct score (0-5) is used rather than the raw withdrawal ratio. The
#' score is a log-normalized rescaling of the same trend signal, carries
#' published category thresholds that apply consistently across baseline and
#' future scenarios, and is better distributed than the raw ratio, which is
#' unbounded, heavily right-skewed, and uses a 9999 sentinel for arid basins.
#'
#' Higher stress = higher sensitivity, so `invert = FALSE` is used when
#' standardizing. A score of 5 covers both "Extremely High (>80%)" and "Arid and
#' Low Water Use"; Aqueduct treats both as top priority.
#'
#' `water_supplies` is used only as a check: the function stops if the distinct
#' `wsd_source_id` values in `stress_data` are not the same set as those in
#' `water_supplies`, which catches a stress table built from a different (or
#' stale) supply list.
#'
#' Citation: Kuzma, S. et al. (2023) Aqueduct 4.0: Updated decision-relevant
#' global water risk indicators. World Resources Institute.
#'
#' @param water_supplies sf object with water supply point locations, used to
#'   verify that `stress_data` covers the same set of `wsd_source_id` values.
#' @param stress_data character path to `waterrights_stress_long.csv` (output of
#'   `process_waterrightsindicator_data()`), or a data.frame with the same
#'   columns. Must include `wsd_source_id`, `park_unit`, `park_name`, `period`,
#'   `stress_score`, and `stress_label`.
#' @param save logical, whether to save output. Default is `FALSE`.
#' @param output_dir character, directory path for saving files. Required if
#'   `save = TRUE`.
#' @param verbose logical, whether to show progress messages. Default is `TRUE`.
#'
#' @return tibble with one row per water supply: `wsd_source_id`, `park_unit`,
#'   `park_name`, `water_stress_label` (Aqueduct category), and raw/norm/rank
#'   columns for `raw_water_stress` (baseline score 0-5, not inverted)
#'
calc_waterrights_sensitivity <- function(
    water_supplies,
    stress_data = here("data",
                       "rapid",
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
      "stress_data path does not exist — run process_waterrightsindicator_data() first" =
        file.exists(stress_data)
    )
    stress_data <- read_csv(stress_data, show_col_types = FALSE)
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
    "wsd_source_id in stress_data does not match water_supplies — rerun process_waterrightsindicator_data() with this water_supplies object" =
      setequal(ws_ids, stress_ids)
  )
  
  if (verbose) {
    message("  ", length(stress_ids), " water supplies matched between stress_data and water_supplies")
  }
  
  # Sensitivity is the baseline period only — one row per water supply
  results <- stress_data %>%
    filter(period == "baseline") %>%
    distinct(wsd_source_id, .keep_all = TRUE) %>%
    select(wsd_source_id, park_unit, park_name,
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
    select(wsd_source_id, park_unit, park_name, everything())

  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    filename <- file.path(output_dir,
                          paste0("waterrights_sensitivity_", Sys.Date(), ".csv"))
    write_csv(final_results, filename)

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
