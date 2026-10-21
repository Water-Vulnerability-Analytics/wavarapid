#' @title Fire exposure: projected change in fire probability
#'
#' @description
#' Assigns each water supply the projected mid-century change in fire
#' probability of its pyrome, from Gao et al. (2021). For each pyrome, the
#' climate models closest to the 10th, 50th and 90th percentiles of the
#' change are kept, along with the percentage of models that show an
#' increase. A supply outside every pyrome (e.g. outside the contiguous US)
#' gets the nearest pyrome. A supply that touches several pyromes gets one
#' row per pyrome.
#'
#' @param water_supplies sf object. The water supplies, as points with
#'   `wsd_source_id`, `park_unit` and `park_name` columns, e.g. from
#'   [wavacore::get_wsd()]. Required.
#' @param fire_prob_data Data frame. The change in fire probability by pyrome
#'   and climate model, with columns `PYROME`, `GCM` and `delta_fp` (percent
#'   change), as made in `rapid_data_processing.Rmd` and saved to
#'   `rapid/inputs/fire_probability_change.csv` in the shared data folder.
#'   Required.
#' @param pyrome_data sf object. The pyrome polygons, with a `PYROME` column,
#'   e.g. `rapid/inputs/pyromes.gpkg` in the shared data folder. Required.
#' @param save Logical. Whether to write the results to a CSV file.
#'   * `TRUE`: writes `fire_exposure_<date>.csv`, named with today's date, to
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
#'   * `TRUE`: shows a progress bar over the water supplies and messages with
#'     the run time and the path of the saved file.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return A tibble, one row per water supply and pyrome, with
#'   `wsd_source_id`, `park_unit`, `park_name`, `PYROME` and:
#'   * `GCM_p10`, `GCM_p50`, `GCM_p90`: the climate model closest to each
#'     percentile;
#'   * `raw_fire_prob_change_p10`, `_p50`, `_p90`: that model's percent change
#'     in fire probability;
#'   * `raw_fire_models_showing_increase`: the percentage of models with an
#'     increase;
#'   * `norm_*` and `rank_*` versions of the three changes from
#'     [add_standardized_indicators()]. A larger increase gets a higher score.
#' @references
#' Gao, P. et al. (2021). Robust projections of future fire probability for
#'   the conterminous United States. Science of the Total Environment, 789,
#'   147872.
#' @export
#' @seealso [calc_fire_sensitivity()] for current wildfire hazard;
#'   [calc_vulnerability_index()].
#' @examplesIf file.exists(wavacore::wava_data_path("rapid", "inputs", "pyromes.gpkg"))
#' ws <- wavacore::get_wsd(park = "ZION")
#' pyromes <- sf::st_read(wavacore::wava_data_path("rapid", "inputs", "pyromes.gpkg"),
#'                        quiet = TRUE)
#' fire_prob <- readr::read_csv(wavacore::wava_data_path("rapid", "inputs",
#'                                                       "fire_probability_change.csv"),
#'                              show_col_types = FALSE)
#' fire_exp <- calc_fire_exposure(ws, fire_prob_data = fire_prob, pyrome_data = pyromes)
#' fire_exp[, c("wsd_source_id", "PYROME", "raw_fire_prob_change_p50",
#'              "raw_fire_models_showing_increase")]

calc_fire_exposure <- function(water_supplies,
                               fire_prob_data,
                               pyrome_data,
                               save = FALSE,
                               output_dir = NULL,
                               verbose = TRUE) {
  
  # Input validation
  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir)
  )
  
  # Start timing
  start_time <- Sys.time()
  
  # convert both objects to terra for faster filtering
  pyrome_data <- terra::vect(pyrome_data)
  water_supplies <- terra::vect(water_supplies)
  
  # Calc percentiles and model agreement for all pyromes
  # Pull 10th, 50th and 90th percentile values for each pyrome
  pyromes_selected <- fire_prob_data %>%
    dplyr::group_by(PYROME) %>%
  dplyr::slice(c(
    which.min(abs(delta_fp - stats::quantile(delta_fp, 0.1))),
    which.min(abs(delta_fp - stats::quantile(delta_fp, 0.5))),
    which.min(abs(delta_fp - stats::quantile(delta_fp, 0.9)))
  )) %>%
  dplyr::mutate(percentile = c("p10", "p50", "p90")) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(
      names_from = percentile,
      values_from = c(GCM, delta_fp),
      names_sep = "_"
    )
  
  # calculate fraction of models showing increase
  model_fractions <- fire_prob_data %>%
    dplyr::group_by(PYROME) %>%
    dplyr::summarize(
      fire_models_showing_increase = 100 * sum(delta_fp > 0, na.rm = TRUE) / dplyr::n(),
      .groups = "drop"
    )
  
  # join percentiles and model agreements
  fire_prob_cleaned <- dplyr::left_join(pyromes_selected, model_fractions, by = "PYROME")
  
  # Create progress bar
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "Processing fire exposure [:bar] :percent (:current/:total) ETA: :eta",
      total = nrow(water_supplies),
      clear = FALSE,
      width = 80
    )
  }
  
  # Process each water supply
  results <- purrr::map_dfr(seq_len(nrow(water_supplies)), function(i) {
    if (verbose) pb$tick()
    
    # Get water supply pyrome
    ws_pyrome <- pyrome_data[terra::relate(pyrome_data, water_supplies[i,], "intersects"), ]
    
    # Fall back to nearest pyrome for water supplies outside pyrome coverage (e.g., non-CONUS)
    if (nrow(ws_pyrome) == 0) {
      nearest_idx <- terra::nearest(water_supplies[i, ], pyrome_data)$to_id
      ws_pyrome <- pyrome_data[nearest_idx, ]
    }
    
    # Extract fire probability data for this pyrome
    fp_data <- fire_prob_cleaned %>% 
      dplyr::filter(PYROME %in% ws_pyrome$PYROME) %>%
      dplyr::mutate(
        wsd_source_id = water_supplies[i, ]$wsd_source_id,
        park_unit = water_supplies[i, ]$park_unit,
        park_name = water_supplies[i, ]$park_name
      ) %>%
      dplyr::select(wsd_source_id, park_unit, park_name, dplyr::everything())
    
    return(fp_data)
  })
  
  # Calculate standardized indicators
  final_results <- results %>%
    # Raw values (renamed for consistency)
    dplyr::rename(
      raw_fire_prob_change_p10 = delta_fp_p10,
      raw_fire_prob_change_p50 = delta_fp_p50,
      raw_fire_prob_change_p90 = delta_fp_p90,
      raw_fire_models_showing_increase = fire_models_showing_increase
    ) %>%
    # Use helper function to add normalized and ranked columns
    add_standardized_indicators("raw_fire_prob_change")
  
  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(output_dir, paste0("fire_exposure_", Sys.Date(), ".csv"))
    readr::write_csv(final_results, filename)
    
    if (verbose) {
      message("Results saved to: ", filename)
    }
  }
  
  # Print timing
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message("\nFire exposure processing completed in ", 
            round(total_time, 2), " ", attr(total_time, "units"))
  }
  
  return(final_results)
}
