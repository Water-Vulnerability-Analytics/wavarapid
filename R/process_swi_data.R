#' @title Prepare the saltwater intrusion points (Adams et al. 2024)
#'
#' @description
#' Reads `output_data.txt` from the Adams et al. (2024) saltwater intrusion
#' dataset, keeps the points inside a box around the contiguous US
#' (longitude -125 to -66, latitude 24 to 50), and returns them as sf points.
#' [calc_swi_exposure()] uses the result.
#'
#' `Case_C` is the inland migration of the saltwater-freshwater interface
#' projected to 2100, in meters, from recharge change and sea level rise
#' combined. Positive values mean the interface moves inland.
#'
#' The point locations come from an older data pipeline (Sawyer/Zhou et al.)
#' and don't match any available HydroSHEDS watershed shapefile one to one.
#'
#' Citation: Adams, K.H., Reager, J.T., Buzzanga, B., David, C., Sawyer, A.,
#' and Hamlington, B. (2024). Climate-induced saltwater intrusion in 2100:
#' recharge-driven severity, sea level-driven prevalence. Geophysical Research
#' Letters.
#'
#' @param data_path String. Path to `output_data.txt` (space-separated, with
#'   `WatershedID`, `lat`, `lon` and `Case_C` columns). Stops with an error if
#'   the file doesn't exist. Default: `output_data.txt` in
#'   `all/coast_comp_data/adams_swi_2024/data` in the shared data folder, see
#'   [wavacore::wava_data_path()].
#' @param save Logical. Whether to write the points to a GeoPackage.
#'   * `TRUE`: writes `swi_points.gpkg` to `output_dir`, creating the folder
#'     if needed and replacing an existing file. Stops with an error if
#'     `output_dir` is `NULL`.
#'   * `FALSE`: writes nothing.
#'
#'   Default: `TRUE`.
#' @param output_dir String or `NULL`. The folder for `swi_points.gpkg`.
#'   * A folder path: where the file is written when `save = TRUE`. Ignored
#'     when `save = FALSE`.
#'   * `NULL`: only allowed with `save = FALSE`.
#'
#'   Default: `rapid/inputs/swi` in the shared data folder, which is where
#'   [calc_swi_exposure()] looks for it.
#'
#' @return An sf object of points (CRS 4326), one row per saltwater intrusion
#'   point, with columns `watershed_id`, `lat`, `lon` and `Case_C` (meters).
#' @export
#' @seealso [calc_swi_exposure()], which assigns these values to water
#'   supplies; [process_coastal_geometry()], which makes the distance-to-coast
#'   table that exposure function also needs.
#' @examplesIf file.exists(wavacore::wava_data_path("all", "coast_comp_data", "adams_swi_2024", "data", "output_data.txt"))
#' # Read the points without writing the GeoPackage
#' swi <- process_swi_data(save = FALSE)
#' nrow(swi)
#' summary(swi$Case_C)

process_swi_data <- function(
    data_path  = wavacore::wava_data_path("all", "coast_comp_data", "adams_swi_2024", "data", "output_data.txt"),
    save       = TRUE,
    output_dir = wavacore::wava_data_path("rapid", "inputs", "swi")) {

  stopifnot("output_data.txt not found \u2014 check data_path" = file.exists(data_path),
            "save is TRUE but output_dir not provided"    = !save || !is.null(output_dir))

  swi_pts <- utils::read.table(data_path, header = TRUE, sep = " ") %>%
    dplyr::as_tibble() %>%
    dplyr::rename(watershed_id = WatershedID) %>%
    # Case_C: meters of inland SWI migration by 2100 (combined recharge + SLR)
    dplyr::select(watershed_id, lat, lon, Case_C) %>%
    dplyr::filter(lon >= -125, lon <= -66, lat >= 24, lat <= 50) %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    pts_path <- file.path(output_dir, "swi_points.gpkg")
    sf::st_write(swi_pts, pts_path, delete_dsn = TRUE, quiet = TRUE)
  }

  return(swi_pts)

}
