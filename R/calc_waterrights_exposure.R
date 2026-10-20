#' Calculate Water Rights Exposure Indicator
#'
#' Calculates water rights exposure from Aqueduct 4.0 projected water stress for
#' the HydroBASINS Level 6 sub-basin containing each NPS water supply, at one of
#' three planning horizons: early-century (2030), mid-century (2050, default) or
#' late-century (2080). Requires the long stress table produced by
#' `process_waterrightsindicator_data()`.
#'
#' Exposure is the water stress a system faces at the selected horizon: the
#' Aqueduct score (0-5) under the worst of the three CMIP6 pathways (RCP2.6,
#' RCP7.0, RCP8.5). It is an *end-state*, not a delta — what the system is up
#' against at that horizon, regardless of which pathway materializes. On ties,
#' the highest forcing is reported as `worst_rcp` and all tied pathways are
#' listed in `tied_rcps`; `"rcp26, rcp70, rcp85"` means stress at that horizon is
#' insensitive to the emissions pathway.
#'
#' Mid-century is the default because it matches NPS planning horizons and
#' precedes the point where Aqueduct holds irrigation demand static. Aqueduct's
#' crop extent projections end in 2050 and are held constant beyond it, and
#' livestock demand is held at 2014 levels throughout, so a `period = "2080"`
#' run rests on materially stronger demand-side assumptions than 2030 or 2050
#' and should be reported with that caveat.
#'
#' The change fields are **not** the metric. They describe the trajectory from
#' baseline to that end-state, and exist to distinguish systems that arrive at a
#' high score by deteriorating from ones that were already there. Two systems can
#' share a score of 5 while one escalated from "Low" and the other never left
#' "Extremely High" — the same exposure, very different stories.
#' `water_stress_change` is future minus baseline, so positive values indicate
#' worsening stress (range -5 to 5).
#'
#' Because "Arid and Low Water Use" and "Extremely High (>80%)" share a severity
#' rank, a move between them reads as no rank change; `water_stress_label_changed`
#' catches those relabels that the score and rank both miss.
#'
#' Higher stress = higher exposure, so `invert = FALSE` is used when
#' standardizing. Standardization pools only within the selected horizon.
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
#'   columns. Must include `wsd_source_id`, `park_unit`, `park_name`, `rcp`,
#'   `period`, `stress_score`, and `stress_label`.
#' @param period character, the single planning horizon to report: `"2030"`
#'   (early-century), `"2050"` (mid-century, the default) or `"2080"`
#'   (late-century). Only one may be selected; the horizon is appended to the
#'   indicator column names and the output filename.
#' @param save logical, whether to save output. Default is `FALSE`.
#' @param output_dir character, directory path for saving files. Required if
#'   `save = TRUE`.
#' @param verbose logical, whether to show progress messages. Default is `TRUE`.
#'
#' @return tibble with one row per water supply: `wsd_source_id`, `park_unit`,
#'   `park_name`, the pathway driving the metric (`worst_rcp`, `tied_rcps`),
#'   raw/norm/rank columns for `raw_water_stress_<period>` (score 0-5, not
#'   inverted), its Aqueduct category (`water_stress_label_<period>`), and the
#'   trajectory from baseline (`baseline_water_stress`,
#'   `baseline_water_stress_label`, `water_stress_change`,
#'   `water_stress_change_direction`, `water_stress_rank_change`,
#'   `water_stress_label_changed`, `water_stress_label_transition`)
#'
#' @export
calc_waterrights_exposure <- function(
    water_supplies,
    stress_data = wava_data_path("rapid",
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
    stress_data <- read_csv(stress_data, show_col_types = FALSE)
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
    filter(period == "baseline") %>%
    distinct(wsd_source_id, .keep_all = TRUE) %>%
    select(wsd_source_id,
           baseline_water_stress       = stress_score,
           baseline_water_stress_label = stress_label)
  
  # Worst-case RCP at the selected horizon — one row per water supply
  worst_future <- stress_data %>%
    filter(period == !!period, !is.na(stress_score)) %>%
    mutate(rcp = factor(as.character(rcp), levels = rcp_levels)) %>%
    group_by(wsd_source_id) %>%
    mutate(
      is_worst  = stress_score == max(stress_score),                       # flag the max-stress rcp(s)
      tied_rcps = paste(sort(unique(as.character(rcp[is_worst]))),         # always list low -> high
                        collapse = ", ")
    ) %>%
    filter(is_worst) %>%
    arrange(desc(rcp), .by_group = TRUE) %>%                               # highest forcing first
    slice(1) %>%                                                           # keep 1 row per water supply
    ungroup() %>%
    select(wsd_source_id, park_unit, park_name,
           worst_rcp = rcp, tied_rcps,
           raw_water_stress   = stress_score,
           water_stress_label = stress_label)
  
  if (nrow(worst_future) == 0) {
    stop("No rows found in stress_data for period == '", period, "'.")
  }
  
  # Trajectory from baseline to the selected end-state
  results <- worst_future %>%
    left_join(baseline_stress, by = "wsd_source_id") %>%
    mutate(
      worst_rcp = as.character(worst_rcp),
      
      # numeric change: future - baseline (positive = worsening)
      water_stress_change = round(raw_water_stress - baseline_water_stress, 3),
      
      water_stress_change_direction = factor(
        case_when(water_stress_change > 0 ~ "Worsened",
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
    rename_with(~ paste0(.x, "_", period),
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
    select(wsd_source_id, park_unit, park_name,
           worst_rcp, tied_rcps,
           all_of(c(raw_col, norm_col, rank_col, label_col)),
           baseline_water_stress, baseline_water_stress_label,
           water_stress_change, water_stress_change_direction,
           water_stress_rank_change, water_stress_label_changed,
           water_stress_label_transition) %>%
    arrange(wsd_source_id)
  
  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(output_dir,
                          paste0("waterrights_exposure_", period, "_", Sys.Date(), ".csv"))
    write_csv(final_results, filename)
    
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