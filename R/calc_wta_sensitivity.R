#' @title WTA sensitivity: historical watershed water stress
#'
#' @description
#' Gives each water supply the 2000-2020 withdrawal-to-availability ratio
#' (WTA) of its watershed, from the Blumberg and Warziniack WTA table. WTA is
#' withdrawals divided by modeled runoff; higher values mean more water
#' stress. It has no units, so it is not divided by area.
#'
#' Each water supply gets the watershed it is in (the lowest HUC code when it
#' sits on a shared boundary), or the nearest one within
#' `max_watershed_distance_km`. The `norm_` and `rank_` columns are relative
#' to the water supplies passed in.
#'
#' @param water_supplies sf object. The water supplies, as points with
#'   `wsd_source_id`, `park_unit` and `park_name` columns, e.g. from
#'   [wavacore::get_wsd()]. Required.
#' @param data_dir String. The folder with `wta_ratios_2000_2020.csv`.
#'   Default: `Latest county-level data/water_demand_projections` in the
#'   shared data folder, see [wavacore::wava_data_path()].
#' @param watersheds String or sf object. The watershed polygons, with a
#'   `huc<huc_level>` column (e.g. `huc8`).
#'   * String: path to a WBD geodatabase; the `WBDHU<huc_level>` layer is
#'     read.
#'   * sf object: the polygons, already read.
#'
#'   Default: `all/WBD_National_GDB.gdb` in the shared data folder.
#' @param huc_level Integer. The watershed size: `2`, `4`, `6` or `8` (HUC2
#'   to HUC8). Picks both the WTA rows and the watershed layer. Default:
#'   `8L`.
#' @param max_watershed_distance_km Numeric. How far, in km, a supply outside
#'   every watershed can be from the nearest one and still use it
#'   (`watershed_match_method = "nearest"`). Supplies farther away get `NA`
#'   and `data_status = "unmatched_watershed"`. Default: `5`.
#' @param save Logical. Whether to write the results to a CSV file.
#'   * `TRUE`: writes `wta_sensitivity_huc<NN>_<date>.csv` (e.g.
#'     `wta_sensitivity_huc08_2026-10-21.csv`) to `output_dir`. The folder
#'     must already exist; `output_dir` can't be `NULL`.
#'   * `FALSE`: writes nothing.
#'
#'   Default: `FALSE`.
#' @param output_dir String or `NULL`. The folder for the CSV file when
#'   `save = TRUE`; ignored otherwise. Default: `NULL`.
#' @param verbose Logical. Whether to print messages.
#'   * `TRUE`: prints how many supplies got a historical WTA, and the path of
#'     the saved file.
#'   * `FALSE`: prints nothing.
#'
#'   Default: `TRUE`.
#'
#' @return A tibble, one row per water supply, with `wsd_source_id`,
#'   `park_unit`, `park_name`, the watershed (`huc`, `watershed_match_method`,
#'   `watershed_distance_km`), `raw_wta` (the 2000-2020 ratio), `data_status`
#'   (`"ok"`, `"unmatched_watershed"` or `"missing_historical_ratio"`), the
#'   `norm_wta` and `rank_wta` columns from [add_standardized_indicators()],
#'   and `aggregation`, `historical_period` and `wta_units`.
#' @export
#' @seealso [calc_wta_exposure()] for the projected change;
#'   [calc_nearby_use_sensitivity()] for county withdrawals.
#' @examplesIf interactive() && dir.exists(wavacore::wava_data_path("all", "WBD_National_GDB.gdb"))
#' # Reads the national HUC8 layer of the WBD geodatabase, so it takes a few
#' # minutes. Read it once and pass it in to reuse it.
#' huc8 <- sf::st_read(wavacore::wava_data_path("all", "WBD_National_GDB.gdb"),
#'                     layer = "WBDHU8", quiet = TRUE)
#' ws <- wavacore::get_wsd(park = c("ZION", "ROMO"))
#' wta_sen <- calc_wta_sensitivity(ws, watersheds = huc8)
#' wta_sen[, c("wsd_source_id", "huc", "raw_wta", "data_status")]

calc_wta_sensitivity <- function(water_supplies,
                                 data_dir = wavacore::wava_data_path("Latest county-level data", "water_demand_projections"),
                                 watersheds = wavacore::wava_data_path("all", "WBD_National_GDB.gdb"),
                                 huc_level = 8L,
                                 max_watershed_distance_km = 5,
                                 save = FALSE,
                                 output_dir = NULL,
                                 verbose = TRUE) {

  aggregation <- sprintf("HUC%02d", huc_level)  # e.g. "HUC08"
  huc_format <- paste0("%0", huc_level, "d")    # pads leading zeros

  # 1. Read historical WTA ratios
  historical <- readr::read_csv(file.path(data_dir, "wta_ratios_2000_2020.csv"),
                         col_types = readr::cols(huc = "c", ratio = "d")) |>
    dplyr::filter(agg == aggregation) |>
    dplyr::transmute(huc = sprintf(huc_format, as.integer(huc)), raw_wta = ratio)

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
  if (length(outside)) {
    nearest <- sf::st_nearest_feature(points[outside, ], watersheds)
    km <- as.numeric(sf::st_distance(points[outside, ], watersheds[nearest, ], by_element = TRUE)) / 1000
    distance[outside] <- km
    index[outside[km <= max_watershed_distance_km]] <- nearest[km <= max_watershed_distance_km]
    method[outside[km <= max_watershed_distance_km]] <- "nearest"
  }

  # 3. Join ratios to supplies and standardize
  final <- sf::st_drop_geometry(water_supplies) |>
    dplyr::select(wsd_source_id, park_unit, park_name) |>
    dplyr::mutate(huc = watersheds$huc[index],
           watershed_match_method = method,
           watershed_distance_km = distance) |>
    dplyr::left_join(historical, by = "huc") |>
    dplyr::mutate(data_status = dplyr::case_when(is.na(huc) ~ "unmatched_watershed",
                                   is.na(raw_wta) ~ "missing_historical_ratio",
                                   TRUE ~ "ok")) |>
    add_standardized_indicators("raw_wta") |>
    dplyr::mutate(aggregation = aggregation,
           historical_period = "2000-2020",
           wta_units = "WTA ratio")

  if (save) {
    filename <- file.path(output_dir, paste0("wta_sensitivity_", tolower(aggregation), "_", Sys.Date(), ".csv"))
    readr::write_csv(final, filename)
    if (verbose) message("Results saved to: ", filename)
  }
  if (verbose) message("WTA sensitivity: ", sum(final$data_status == "ok"), "/", nrow(final), " supplies with historical WTA")
  final
}
