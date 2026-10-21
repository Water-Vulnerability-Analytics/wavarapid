#' @title Water rights exposure: projected water stress (Aqueduct 4.0)
#'
#' @description
#' Gives each water supply the Aqueduct 4.0 water stress score of its
#' HydroBASINS level 6 sub-basin at one planning horizon (2030, 2050 or
#' 2080), under the worst of the three CMIP6 pathways (RCP2.6, RCP7.0,
#' RCP8.5), from the table made by [process_waterrightsindicator_data()].
#'
#' Exposure is an end state, not a change: the stress a system faces at that
#' horizon, whichever pathway happens. On ties, the highest forcing is
#' reported as `worst_rcp` and all tied pathways are listed in `tied_rcps`;
#' `"rcp26, rcp70, rcp85"` means stress at that horizon doesn't depend on the
#' pathway.
#'
#' The change columns are not the indicator. They show how the supply got
#' there from the baseline: two supplies can both score 5 while one rose
#' from "Low" and the other was always "Extremely High".
#' `water_stress_change` is future minus baseline, so positive means worse
#' (from -5 to 5). "Arid and Low Water Use" and "Extremely High (>80%)" share
#' a rank, so a move between them shows no rank change;
#' `water_stress_label_changed` catches those relabels.
#'
#' Standardizing uses only the selected horizon.
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
#'   `park_unit`, `park_name`, `rcp`, `period`, `stress_score` and
#'   `stress_label`.
#'   * String: path to `waterrights_stress_long.csv`. Stops with an error if
#'     the file doesn't exist.
#'   * Data frame: the same table, already read.
#'
#'   Default: `rapid/inputs/WaterRightsIndicator/waterrights_stress_long.csv`
#'   in the shared data folder, see [wavacore::wava_data_path()].
#' @param period String. The one planning horizon to report. The horizon is
#'   added to the indicator column names and the file name.
#'   * `"2030"`: early century.
#'   * `"2050"`: mid century. Matches NPS planning horizons and comes before
#'     Aqueduct holds irrigation demand fixed.
#'   * `"2080"`: late century. Aqueduct's crop extent projections end in 2050
#'     and are held constant after it, and livestock demand stays at 2014
#'     levels throughout, so report 2080 results with that caveat.
#'
#'   Anything else stops with an error. Default: `"2050"`.
#' @param save Logical. Whether to write the results to a CSV file.
#'   * `TRUE`: writes `waterrights_exposure_<period>_<date>.csv`, named with
#'     the period and today's date, to `output_dir`, creating the folder if
#'     needed. Stops with an error if
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
#'   * `TRUE`: prints each step, how many supplies got a value, the run time
#'     and the path of the saved file.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return A tibble, one row per water supply, with `wsd_source_id`,
#'   `park_unit`, `park_name` and:
#'   * `worst_rcp` and `tied_rcps`: the pathway behind the score;
#'   * `raw_water_stress_<period>` (score 0-5) with its `norm_` and `rank_`
#'     versions from [add_standardized_indicators()] (higher stress gets a
#'     higher score), and `water_stress_label_<period>`, the Aqueduct
#'     category;
#'   * the path from the baseline: `baseline_water_stress`,
#'     `baseline_water_stress_label`, `water_stress_change`,
#'     `water_stress_change_direction`, `water_stress_rank_change`,
#'     `water_stress_label_changed` and `water_stress_label_transition`.
#' @export
#' @seealso [calc_waterrights_sensitivity()] for the baseline stress;
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
#' # Mid century (the default) and late century
#' calc_waterrights_exposure(ws, stress_data = stress)
#' calc_waterrights_exposure(ws, stress_data = stress, period = "2080")

calc_waterrights_exposure <- function(
    water_supplies,
    stress_data = wavacore::wava_data_path("rapid",
                       "inputs",
                       "WaterRightsIndicator",
                       "waterrights_stress_long.csv"),
    period = "2050",   #options: "2030", "2050", or "2080"
    save = FALSE,
    output_dir = NULL,
    verbose = TRUE) {
  
  # Input validation
  stopifnot(
    "water_supplies must be an sf object"         = inherits(water_supplies, "sf"),
    "stress_data must be a file path or data.frame" =
      is.character(stress_data) || is.data.frame(stress_data),
    "period must be a single value" = length(period) == 1,
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir)
  )
  
  # Exactly one horizon; errors listing the valid options if not matched
  period <- match.arg(as.character(period), c("2030", "2050", "2080"))
  
  # Start timing
  if (verbose) {
    message("Starting water rights exposure calculation (", period, ")...")
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
                        "rcp", "period", "stress_score", "stress_label")
  
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
  
  # Horizon-suffixed indicator names, e.g. raw_water_stress_2050
  raw_col   <- paste0("raw_water_stress_", period)
  norm_col  <- paste0("norm_water_stress_", period)
  rank_col  <- paste0("rank_water_stress_", period)
  label_col <- paste0("water_stress_label_", period)
  
  # Severity rank: Arid ties with Extremely High (both score 5), so a move
  # between those two categories correctly reads as no rank change
  stress_rank <- c("Low (<10%)" = 1, "Low - Medium (10-20%)" = 2,
                   "Medium - High (20-40%)" = 3, "High (40-80%)" = 4,
                   "Extremely High (>80%)" = 5, "Arid and Low Water Use" = 5)
  
  # Forcing order, so the highest pathway wins a tie on the worst score
  rcp_levels <- c("rcp26", "rcp70", "rcp85")
  
  # Baseline (Sensitivity) values, carried as trajectory context
  baseline_stress <- stress_data %>%
    dplyr::filter(period == "baseline") %>%
    dplyr::distinct(wsd_source_id, .keep_all = TRUE) %>%
    dplyr::select(wsd_source_id,
           baseline_water_stress       = stress_score,
           baseline_water_stress_label = stress_label)
  
  # Worst-case RCP at the selected horizon — one row per water supply
  worst_future <- stress_data %>%
    dplyr::filter(period == !!period, !is.na(stress_score)) %>%
    dplyr::mutate(rcp = factor(as.character(rcp), levels = rcp_levels)) %>%
    dplyr::group_by(wsd_source_id) %>%
    dplyr::mutate(
      is_worst  = stress_score == max(stress_score),                       # flag the max-stress rcp(s)
      tied_rcps = paste(sort(unique(as.character(rcp[is_worst]))),         # always list low -> high
                        collapse = ", ")
    ) %>%
    dplyr::filter(is_worst) %>%
    dplyr::arrange(dplyr::desc(rcp), .by_group = TRUE) %>%                               # highest forcing first
    dplyr::slice(1) %>%                                                           # keep 1 row per water supply
    dplyr::ungroup() %>%
    dplyr::select(wsd_source_id, park_unit, park_name,
           worst_rcp = rcp, tied_rcps,
           raw_water_stress   = stress_score,
           water_stress_label = stress_label)
  
  if (nrow(worst_future) == 0) {
    stop("No rows found in stress_data for period == '", period, "'.")
  }
  
  # Trajectory from baseline to the selected end-state
  results <- worst_future %>%
    dplyr::left_join(baseline_stress, by = "wsd_source_id") %>%
    dplyr::mutate(
      worst_rcp = as.character(worst_rcp),
      
      # numeric change: future - baseline (positive = worsening)
      water_stress_change = round(raw_water_stress - baseline_water_stress, 3),
      
      water_stress_change_direction = factor(
        dplyr::case_when(water_stress_change > 0 ~ "Worsened",
                  water_stress_change < 0 ~ "Improved",
                  TRUE                    ~ "Unchanged"),
        levels = c("Improved", "Unchanged", "Worsened")),
      
      # category movement (Arid ties with Extremely High -> reads 0)
      water_stress_rank_change =
        unname(stress_rank[as.character(water_stress_label)]) -
        unname(stress_rank[as.character(baseline_water_stress_label)]),
      
      # literal label move - catches Arid <-> Extremely High relabels
      water_stress_label_changed =
        as.character(water_stress_label) != as.character(baseline_water_stress_label),
      
      water_stress_label_transition = paste(baseline_water_stress_label, "->",
                                            water_stress_label)
    ) %>%
    # append the horizon to the indicator columns
    dplyr::rename_with(~ paste0(.x, "_", period),
                .cols = c(raw_water_stress, water_stress_label))
  
  if (verbose) {
    message("  ", period, " worst-case water stress assigned: ",
            sum(!is.na(results[[raw_col]])),
            " | tied across all three RCPs: ",
            sum(results$tied_rcps == "rcp26, rcp70, rcp85", na.rm = TRUE))
  }
  
  # invert = FALSE: higher water stress = higher exposure
  final_results <- results %>%
    add_standardized_indicators(raw_col, invert = FALSE) %>%
    # reorder columns following standardized pattern
    dplyr::select(wsd_source_id, park_unit, park_name,
           worst_rcp, tied_rcps,
           dplyr::all_of(c(raw_col, norm_col, rank_col, label_col)),
           baseline_water_stress, baseline_water_stress_label,
           water_stress_change, water_stress_change_direction,
           water_stress_rank_change, water_stress_label_changed,
           water_stress_label_transition) %>%
    dplyr::arrange(wsd_source_id)
  
  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(output_dir,
                          paste0("waterrights_exposure_", period, "_", Sys.Date(), ".csv"))
    readr::write_csv(final_results, filename)
    
    if (verbose) {
      message("Results saved to: ", filename)
    }
  }
  
  # Print timing
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message("\nWater rights exposure processing completed in ",
            round(total_time, 2), " ", attr(total_time, "units"))
  }
  
  return(final_results)
}
