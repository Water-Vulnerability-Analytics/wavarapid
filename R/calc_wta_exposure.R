#' @title WTA exposure: projected change in watershed water stress
#'
#' @description
#' Calculates the change in the withdrawal-to-availability ratio (WTA) of each
#' water supply's watershed, 2040-2070 minus 2000-2020, from the Blumberg and
#' Warziniack WTA tables. WTA is withdrawals (public supply, irrigation,
#' thermoelectric) divided by modeled runoff. It has no units, so the change
#' is not divided by area. Each climate model and emission scenario pair is
#' one ensemble member; ssp126 is left out. For the 10th, 50th and 90th
#' percentiles, the member closest to each is reported (ties go to the first
#' model, then scenario, alphabetically).
#'
#' Each water supply gets the watershed it is in (the lowest HUC code when it
#' sits on a shared boundary), or the nearest one within
#' `max_watershed_distance_km`. The `norm_` and `rank_` columns are relative
#' to the water supplies passed in.
#'
#' @param water_supplies sf object. The water supplies, as points with
#'   `wsd_source_id`, `park_unit` and `park_name` columns, e.g. from
#'   [wavacore::get_wsd()]. Required.
#' @param data_dir String. The folder with `wta_ratios_2000_2020.csv` and
#'   `wta_ratios_2040_2070.csv`. Default:
#'   `Latest county-level data/water_demand_projections` in the shared data
#'   folder, see [wavacore::wava_data_path()].
#' @param watersheds String or sf object. The watershed polygons, with a
#'   `huc<huc_level>` column (e.g. `huc8`).
#'   * String: path to a WBD geodatabase; the `WBDHU<huc_level>` layer is
#'     read.
#'   * sf object: the polygons, already read (faster when calling the
#'     function more than once).
#'
#'   Default: `all/WBD_National_GDB.gdb` in the shared data folder.
#' @param huc_level Integer. The watershed size: `2`, `4`, `6` or `8` (HUC2
#'   to HUC8). Picks both the WTA rows and the watershed layer. Default:
#'   `8L`.
#' @param emission_scenarios Character vector or `NULL`. The emission
#'   scenarios to use.
#'   * A character vector, e.g. `"ssp245"`: keeps only those members.
#'   * `NULL`: keeps every scenario except ssp126.
#'
#'   Default: `NULL`.
#' @param climate_models Character vector or `NULL`. The climate models to
#'   use, in the same way as `emission_scenarios`. Default: `NULL` (all).
#' @param max_watershed_distance_km Numeric. How far, in km, a supply outside
#'   every watershed can be from the nearest one and still use it.
#'   * A positive distance: supplies within it use the nearest watershed
#'     (`watershed_match_method = "nearest"`); supplies farther away get `NA`.
#'   * `0`: only supplies inside a watershed get a value.
#'
#'   Default: `5`.
#' @param save Logical. Whether to write the results to a CSV file.
#'   * `TRUE`: writes `wta_exposure_huc<NN>_<date>.csv` (e.g.
#'     `wta_exposure_huc08_2026-10-21.csv`) to `output_dir`. The folder must
#'     already exist; `output_dir` can't be `NULL`.
#'   * `FALSE`: writes nothing.
#'
#'   Default: `FALSE`.
#' @param output_dir String or `NULL`. The folder for the CSV file when
#'   `save = TRUE`; ignored otherwise. Default: `NULL`.
#' @param verbose Logical. Whether to print messages.
#'   * `TRUE`: prints how many supplies got a finite change, and the path of
#'     the saved file.
#'   * `FALSE`: prints nothing.
#'
#'   Default: `TRUE`.
#'
#' @return A tibble, one row per water supply, with `wsd_source_id`,
#'   `park_unit`, `park_name`, the watershed (`huc`, `watershed_match_method`,
#'   `watershed_distance_km`), `baseline_wta` and:
#'   * `raw_wta_change_p10`, `_p50`, `_p90`: the change in WTA of the member
#'     closest to that percentile, with `raw_wta_percent_change_*`,
#'     `future_wta_*`, `climate_model_*` and `emission_scenario_*` for the
#'     same member;
#'   * `raw_wta_models_showing_increase_pct`: the percentage of members with
#'     an increase;
#'   * `n_models`, `n_future_records`, `n_valid_models` and `data_status`
#'     (`"ok"`, `"partial_ensemble"`, `"unmatched_watershed"`,
#'     `"missing_historical_watershed"`, `"missing_historical_ratio"`,
#'     `"missing_future_watershed"` or `"missing_future_ratio"`);
#'   * `norm_*` and `rank_*` versions of the change, percent change and
#'     model percentage columns from [add_standardized_indicators()];
#'   * `aggregation`, `baseline_period`, `future_period`, `wta_change_units`,
#'     `selected_emission_scenarios` and `selected_climate_models`.
#' @export
#' @seealso [calc_wta_sensitivity()] for the historical ratio;
#'   [calc_nearby_use_exposure()] for county withdrawals from the same
#'   projections.
#' @examplesIf interactive() && dir.exists(wavacore::wava_data_path("all", "WBD_National_GDB.gdb"))
#' # Reads the national HUC8 layer of the WBD geodatabase, so it takes a few
#' # minutes. Read it once and pass it in to reuse it.
#' huc8 <- sf::st_read(wavacore::wava_data_path("all", "WBD_National_GDB.gdb"),
#'                     layer = "WBDHU8", quiet = TRUE)
#' ws <- wavacore::get_wsd(park = c("ZION", "ROMO"))
#' wta_exp <- calc_wta_exposure(ws, watersheds = huc8)
#' wta_exp[, c("wsd_source_id", "huc", "raw_wta_change_p50", "data_status")]

calc_wta_exposure <- function(water_supplies,
                              data_dir = wavacore::wava_data_path("Latest county-level data", "water_demand_projections"),
                              watersheds = wavacore::wava_data_path("all", "WBD_National_GDB.gdb"),
                              huc_level = 8L,
                              emission_scenarios = NULL,
                              climate_models = NULL,
                              max_watershed_distance_km = 5,
                              save = FALSE,
                              output_dir = NULL,
                              verbose = TRUE) {

  aggregation <- sprintf("HUC%02d", huc_level)  # e.g. "HUC08"
  huc_format <- paste0("%0", huc_level, "d")    # pads leading zeros

  # 1. Read WTA ratios
  historical <- readr::read_csv(file.path(data_dir, "wta_ratios_2000_2020.csv"),
                         col_types = readr::cols(huc = "c", ratio = "d")) |>
    dplyr::filter(agg == aggregation) |>
    dplyr::transmute(huc = sprintf(huc_format, as.integer(huc)), baseline_wta = ratio,
              historical_record = TRUE)

  future <- readr::read_csv(file.path(data_dir, "wta_ratios_2040_2070.csv"),
                     col_types = readr::cols(huc = "c", ratio = "d")) |>
    dplyr::filter(agg == aggregation, emission_scenario != "ssp126") |>
    dplyr::transmute(huc = sprintf(huc_format, as.integer(huc)), climate_model,
              emission_scenario, future_wta = ratio, future_record = TRUE)
  if (!is.null(emission_scenarios)) future <- dplyr::filter(future, emission_scenario %in% emission_scenarios)
  if (!is.null(climate_models)) future <- dplyr::filter(future, climate_model %in% climate_models)
  members <- dplyr::distinct(future, climate_model, emission_scenario)

  # 2. Assign each supply to a watershed (lowest HUC on a shared boundary)
  if (!inherits(watersheds, "sf")) {
    watersheds <- sf::st_read(watersheds, layer = paste0("WBDHU", huc_level), quiet = TRUE)
  }
  watersheds <- watersheds |>
    dplyr::mutate(huc = sprintf(huc_format, as.integer(.data[[paste0("huc", huc_level)]]))) |>
    sf::st_make_valid() |>
    sf::st_transform(5070)
  points <- sf::st_transform(water_supplies, 5070)

  index <- sapply(sf::st_intersects(points, watersheds),
                  function(i) if (length(i)) i[order(watersheds$huc[i])][1] else NA_integer_)
  distance <- ifelse(is.na(index), NA_real_, 0)
  method <- ifelse(is.na(index), "unmatched", "intersects")

  outside <- which(is.na(index) & !sf::st_is_empty(points))
  if (length(outside) && max_watershed_distance_km > 0) {
    nearest <- sf::st_nearest_feature(points[outside, ], watersheds)
    km <- as.numeric(sf::st_distance(points[outside, ], watersheds[nearest, ], by_element = TRUE)) / 1000
    distance[outside] <- km
    index[outside[km <= max_watershed_distance_km]] <- nearest[km <= max_watershed_distance_km]
    method[outside[km <= max_watershed_distance_km]] <- "nearest"
  }

  locations <- sf::st_drop_geometry(water_supplies) |>
    dplyr::select(wsd_source_id, park_unit, park_name) |>
    dplyr::mutate(huc = watersheds$huc[index],
           watershed_match_method = method,
           watershed_distance_km = distance)

  # 3. Change per watershed and member (all members kept, to report coverage)
  changes <- tidyr::crossing(huc = unique(stats::na.omit(locations$huc)), members) |>
    dplyr::left_join(historical, by = "huc") |>
    dplyr::left_join(future, by = c("huc", "climate_model", "emission_scenario")) |>
    dplyr::mutate(wta_change = future_wta - baseline_wta,
           percent_change = dplyr::if_else(baseline_wta > 0, 100 * wta_change / baseline_wta, NA_real_)) |>
    dplyr::arrange(huc, climate_model, emission_scenario)

  # 4. Summarize the ensemble per watershed
  results <- changes |>
    dplyr::group_by(huc) |>
    dplyr::summarise(baseline_wta = dplyr::first(baseline_wta),
              historical_record = any(historical_record %in% TRUE),
              n_models = dplyr::n(),
              n_future_records = sum(future_record %in% TRUE),
              n_valid_models = sum(is.finite(wta_change)),
              raw_wta_models_showing_increase_pct = if (n_valid_models > 0) {
                100 * sum(wta_change > 0, na.rm = TRUE) / n_valid_models
              } else NA_real_) |>
    dplyr::mutate(data_status = dplyr::case_when(
      !historical_record ~ "missing_historical_watershed",
      !is.finite(baseline_wta) ~ "missing_historical_ratio",
      n_future_records == 0 ~ "missing_future_watershed",
      n_valid_models == 0 ~ "missing_future_ratio",
      n_valid_models < n_models ~ "partial_ensemble",
      TRUE ~ "ok")) |>
    dplyr::select(-historical_record)

  # Member nearest each percentile (a tiny tolerance ignores float roundoff in ties)
  for (p in c(10, 50, 90)) {
    nearest_member <- changes |>
      dplyr::filter(is.finite(wta_change)) |>
      dplyr::group_by(huc) |>
      dplyr::mutate(distance = abs(wta_change - stats::quantile(wta_change, p / 100))) |>
      dplyr::filter(distance <= min(distance) + 16 * .Machine$double.eps * max(1, abs(wta_change))) |>
      dplyr::slice(1) |>
      dplyr::ungroup() |>
      dplyr::transmute(huc, raw_wta_change = wta_change, raw_wta_percent_change = percent_change,
                future_wta, climate_model, emission_scenario) |>
      dplyr::rename_with(~ paste0(.x, "_p", p), -huc)
    results <- dplyr::left_join(results, nearest_member, by = "huc")
  }

  # 5. Join to supplies and standardize
  final <- locations |>
    dplyr::left_join(results, by = "huc") |>
    add_standardized_indicators("raw_wta_change") |>
    add_standardized_indicators("raw_wta_percent_change") |>
    add_standardized_indicators("raw_wta_models_showing_increase_pct") |>
    dplyr::mutate(aggregation = aggregation,
           baseline_period = "2000-2020",
           future_period = "2040-2070",
           wta_change_units = "absolute WTA ratio change",
           selected_emission_scenarios = paste(sort(unique(members$emission_scenario)), collapse = ","),
           selected_climate_models = paste(sort(unique(members$climate_model)), collapse = ","),
           data_status = dplyr::if_else(is.na(huc), "unmatched_watershed", data_status))

  if (save) {
    filename <- file.path(output_dir, paste0("wta_exposure_", tolower(aggregation), "_", Sys.Date(), ".csv"))
    readr::write_csv(final, filename)
    if (verbose) message("Results saved to: ", filename)
  }
  if (verbose) {
    message("WTA exposure: ", sum(is.finite(final$raw_wta_change_p50)), "/", nrow(final),
            " supplies with finite changes")
  }
  final
}
