# Shared fixtures for the rapid_app tests.
#
# testthat sources every helper-*.R file before running the tests, so
# make_fake_supplies() is available to all test files without an explicit
# source() call.

#' Build a synthetic water-supply table for calc_vulnerability_index()
#'
#' calc_vulnerability_index() operates on a plain data frame -- it never touches
#' geometry or coordinates -- so the fixture only needs `source_type` plus the
#' raw indicator / model-agreement columns the function reads. Values are
#' deterministic for a given `seed`.
#'
#' @param n            number of rows (water supplies)
#' @param source_type  single value recycled to every row, or a length-`n`
#'                     vector to set each row's type individually
#' @param na_cols      character vector of columns to force to all-NA (for the
#'                     "no data for this indicator" edge case)
#' @param seed         RNG seed, so a test can assert on exact values
make_fake_supplies <- function(n = 6,
                               source_type = "groundwater",
                               na_cols = character(0),
                               seed = 1) {
  set.seed(seed)

  exp_indicators <- c(
    "exp_runoff_change", "exp_precip_change",
    "exp_runoff_timing_change", "exp_precip_timing_change",
    "exp_drought_change", "exp_nearby_use_change", "exp_auc_change",
    "exp_inundation_slr", "exp_swi", "exp_storm_surge",
    "exp_fire_prob_change", "exp_temp_change"
  )
  sen_indicators <- c(
    "sen_visitation_trend", "sen_competition", "sen_auc_change",
    "sen_wildfire_hazard", "sen_flood_risk", "sen_inundation_current",
    "sen_temp_trend", "sen_runoff_trend", "sen_precip_trend",
    "sen_drought_trend", "sen_source_type"
  )
  # Model-agreement weights stored as percentages (0-100) in the source data.
  pct_weights <- c(
    "exp_runoff_model_agree", "exp_precip_model_agree",
    "exp_drought_model_agree", "exp_nearby_use_model_agree",
    "exp_fire_model_agree"
  )
  # Timing weights are stored as fractions (0-1).
  frac_weights <- c(
    "exp_runoff_timing_change_model_agree",
    "exp_precip_timing_change_model_agree"
  )

  df <- data.frame(
    wsd_source_id = sprintf("FAKE_%02d", seq_len(n)),
    source_type   = rep_len(source_type, n),
    stringsAsFactors = FALSE
  )
  for (col in c(exp_indicators, sen_indicators)) {
    df[[col]] <- round(stats::runif(n, -5, 5), 3)
  }
  for (col in pct_weights)  df[[col]] <- round(stats::runif(n, 0, 100), 1)
  for (col in frac_weights) df[[col]] <- round(stats::runif(n, 0, 1), 3)

  df[na_cols] <- NA_real_
  df
}
