#' Calculate historical Withdrawal-to-Availability (WTA) sensitivity
#'
#' The 2000-2020 watershed WTA ratio for each water supply. WTA is a
#' dimensionless ratio, so it is not divided by area. Higher values mean more
#' historical water stress. norm_/rank_ columns are across supplies.
#'
#' @param water_supplies sf points with wsd_source_id, park_unit, park_name.
#' @param data_dir Folder with wta_ratios_2000_2020.csv.
#' @param watersheds sf polygons, or a WBD geodatabase path.
#' @param huc_level HUC width: 2, 4, 6, or 8 (default).
#' @param max_watershed_distance_km Nearest-watershed limit for points outside.
#' @param save,output_dir Save wta_sensitivity_hucNN_YYYY-MM-DD.csv to output_dir.
#' @param verbose Report coverage.
#' @return Tibble with one row per supply.
calc_wta_sensitivity <- function(water_supplies,
                                 data_dir = "data/Latest county-level data/water_demand_projections",
                                 watersheds = "data/all/WBD_National_GDB.gdb",
                                 huc_level = 8L,
                                 max_watershed_distance_km = 5,
                                 save = FALSE,
                                 output_dir = NULL,
                                 verbose = TRUE) {

  aggregation <- sprintf("HUC%02d", huc_level)  # e.g. "HUC08"
  huc_format <- paste0("%0", huc_level, "d")    # pads leading zeros

  # 1. Read historical WTA ratios
  historical <- read_csv(file.path(data_dir, "wta_ratios_2000_2020.csv"),
                         col_types = cols(huc = "c", ratio = "d")) |>
    filter(agg == aggregation) |>
    transmute(huc = sprintf(huc_format, as.integer(huc)), raw_wta = ratio)

  # 2. Assign each supply to a watershed (lowest HUC on a shared boundary)
  if (!inherits(watersheds, "sf")) {
    watersheds <- st_read(watersheds, layer = paste0("WBDHU", huc_level), quiet = TRUE)
  }
  watersheds <- watersheds |>
    mutate(huc = sprintf(huc_format, as.integer(.data[[paste0("huc", huc_level)]]))) |>
    st_make_valid() |>
    st_transform(5070)
  points <- st_transform(water_supplies, 5070)

  index <- sapply(st_intersects(points, watersheds),
                  function(i) if (length(i)) i[order(watersheds$huc[i])][1] else NA_integer_)
  distance <- ifelse(is.na(index), NA_real_, 0)
  method <- ifelse(is.na(index), "unmatched", "intersects")

  outside <- which(is.na(index) & !st_is_empty(points))
  if (length(outside)) {
    nearest <- st_nearest_feature(points[outside, ], watersheds)
    km <- as.numeric(st_distance(points[outside, ], watersheds[nearest, ], by_element = TRUE)) / 1000
    distance[outside] <- km
    index[outside[km <= max_watershed_distance_km]] <- nearest[km <= max_watershed_distance_km]
    method[outside[km <= max_watershed_distance_km]] <- "nearest"
  }

  # 3. Join ratios to supplies and standardize
  final <- st_drop_geometry(water_supplies) |>
    select(wsd_source_id, park_unit, park_name) |>
    mutate(huc = watersheds$huc[index],
           watershed_match_method = method,
           watershed_distance_km = distance) |>
    left_join(historical, by = "huc") |>
    mutate(data_status = case_when(is.na(huc) ~ "unmatched_watershed",
                                   is.na(raw_wta) ~ "missing_historical_ratio",
                                   TRUE ~ "ok")) |>
    add_standardized_indicators("raw_wta") |>
    mutate(aggregation = aggregation,
           historical_period = "2000-2020",
           wta_units = "WTA ratio")

  if (save) {
    filename <- file.path(output_dir, paste0("wta_sensitivity_", tolower(aggregation), "_", Sys.Date(), ".csv"))
    write_csv(final, filename)
    if (verbose) message("Results saved to: ", filename)
  }
  if (verbose) message("WTA sensitivity: ", sum(final$data_status == "ok"), "/", nrow(final), " supplies with historical WTA")
  final
}
