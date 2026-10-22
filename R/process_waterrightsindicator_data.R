#' @title Prepare the Aqueduct 4.0 water stress table for the water rights indicators
#'
#' @description
#' Builds the long water stress table that [calc_waterrights_sensitivity()]
#' and [calc_waterrights_exposure()] read:
#' 1. Downloads the TIGER/Line state boundaries ([tigris::states()]) and
#'    keeps the contiguous US, so it needs internet.
#' 2. Reads the Aqueduct 4.0 `baseline_annual` and `future_annual` layers
#'    and keeps the sub-basins that touch the contiguous US.
#' 3. Attaches each water supply to its nearest sub-basin, separately for
#'    the baseline and future layers, because they use different HydroBASINS
#'    resolutions (about 3,137 and 1,011 sub-basins in the contiguous US).
#'    Both `pfaf_id_baseline` and `pfaf_id_future` are kept.
#' 4. Reshapes to one row per water supply, pathway and period: the baseline
#'    (1979-2019) plus RCP2.6, RCP7.0 and RCP8.5 at 2030, 2050 and 2080, so
#'    10 rows per supply.
#'
#' The water stress score is log-normalized from 0 to 5; higher means more
#' competition for water. "Arid and Low Water Use" scores 5, the same as
#' "Extremely High", and `stress_rank` gives both 5 so category changes read
#' correctly. `stress_raw` keeps Aqueduct's 9999 code for arid sub-basins;
#' recode it when plotting. State surface and groundwater doctrines are not
#' included; they are legal context kept separately.
#'
#' s2 geometry is turned off while the function runs (Aqueduct polygons are
#' large and often invalid) and turned back on when it ends.
#'
#' Data downloads (active as of August 2026):
#' - Sign-up form: <https://www.wri.org/data/aqueduct-global-maps-40-data>
#' - Direct link to the zip with the `.gdb`:
#'   <https://files.wri.org/aqueduct/aqueduct-4-0-water-risk-data.zip>
#' - Data dictionary:
#'   <https://github.com/wri/Aqueduct40/blob/master/data_dictionary_water-risk-atlas.md>
#'
#' Citation: Kuzma, S. et al. (2023) Aqueduct 4.0: Updated decision-relevant
#' global water risk indicators. World Resources Institute.
#'
#' @param water_supplies sf object. The water supplies, as points, with the
#'   WSD columns `wsd_source_id`, `park_unit`, `park_name`, `region`,
#'   `source_state`, `state`, `water_system_name`, `system_type`,
#'   `source_longitude` and `source_latitude`, e.g.
#'   `wavacore::get_wsd(rapid = TRUE)`. Supplies with no `source_state` use
#'   `state`; supplies outside the contiguous US are dropped, as are repeated
#'   `wsd_source_id`s. Required.
#' @param gdb_path String. Path to the Aqueduct 4.0 `.gdb`. Stops with an
#'   error if it doesn't exist. Default:
#'   `all/aqueduct_v4/Aqueduct40_waterrisk_download_Y2023M07D05/GDB/Aq40_Y2023D07M05.gdb`
#'   in the shared data folder, see [wavacore::wava_data_path()].
#' @param save Logical. Whether to write the table.
#'   * `TRUE`: writes `waterrights_stress_long.csv` to `output_dir`, creating
#'     the folder if needed and replacing an existing file. Stops with an
#'     error if `output_dir` is `NULL`.
#'   * `FALSE`: writes nothing, whatever `save_intermediate` is.
#'
#'   Default: `TRUE`.
#' @param output_dir String or `NULL`. The folder for the output files.
#'   * A folder path: where the files are written when `save = TRUE`.
#'     Ignored when `save = FALSE`.
#'   * `NULL`: only allowed with `save = FALSE`.
#'
#'   Default: `rapid/inputs/WaterRightsIndicator` in the shared data folder,
#'   which is where the water rights indicators look for the CSV.
#' @param save_intermediate Logical. Whether to also write the clipped
#'   Aqueduct layers (useful for maps). Only used when `save = TRUE`.
#'   * `TRUE`: also writes `CONUS_aqueduct_baseline_annual.gpkg` and
#'     `CONUS_aqueduct_future_annual.gpkg` to `output_dir`.
#'   * `FALSE`: writes only the CSV.
#'
#'   Default: `TRUE`.
#' @param verbose Logical. Whether to print progress.
#'   * `TRUE`: prints each step, the number of sub-basins and rows, the saved
#'     paths and the run time.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#' @param quiet_warnings Logical. Whether to hide the expected warnings and
#'   messages from the spatial steps.
#'   * `TRUE`: hides them (the GDAL `organizePolygons()` warning and the s2
#'     "assumes planar" notes). Errors still show.
#'   * `FALSE`: shows them.
#'
#'   Default: `TRUE`.
#'
#' @return A tibble, one row per water supply, pathway and period, with
#'   columns `wsd_source_id`, `park_unit`, `park_name`, `region`,
#'   `source_state`, `state_name`, `water_system_name`, `system_type`,
#'   `source_longitude`, `source_latitude`, `aq30_id`, `pfaf_id_baseline`,
#'   `pfaf_id_future`, `rcp` (`baseline`, `rcp26`, `rcp70`, `rcp85`),
#'   `period` (`baseline`, `2030`, `2050`, `2080`), `stress_raw`,
#'   `stress_score` (0-5), `stress_label` (ordered factor) and `stress_rank`
#'   (1-5).
#' @export
#' @seealso [calc_waterrights_sensitivity()] (uses the baseline rows) and
#'   [calc_waterrights_exposure()] (uses one future period).
#' @examplesIf interactive() && file.exists(wavacore::wava_data_path("all", "aqueduct_v4", "Aqueduct40_waterrisk_download_Y2023M07D05", "GDB", "Aq40_Y2023D07M05.gdb"))
#' # Downloads the state boundaries and reads the global Aqueduct layers, so
#' # it needs internet and takes several minutes
#' ws <- wavacore::get_wsd(park = c("ZION", "ROMO"))
#' stress <- process_waterrightsindicator_data(ws, save = FALSE)
#' table(stress$period, stress$rcp)

process_waterrightsindicator_data <- function(
    water_supplies,
    gdb_path          = wavacore::wava_data_path("all", "aqueduct_v4",
                             "Aqueduct40_waterrisk_download_Y2023M07D05", "GDB",
                             "Aq40_Y2023D07M05.gdb"),
    save              = TRUE,
    output_dir        = wavacore::wava_data_path("rapid", "inputs", "WaterRightsIndicator"),
    save_intermediate = TRUE,
    verbose           = TRUE,
    quiet_warnings    = TRUE) {

  stopifnot(
    "water_supplies must be an sf object"       = inherits(water_supplies, "sf"),
    "water_supplies must include wsd_source_id" = "wsd_source_id" %in% names(water_supplies),
    "Aqueduct .gdb not found - check gdb_path"  = file.exists(gdb_path),
    "save is TRUE but output_dir not provided"  = !save || !is.null(output_dir)
  )

  if (verbose) {
    message("Starting Aqueduct 4.0 water rights indicator processing...")
    start_time <- Sys.time()
  }

  # Silences only the expected noise (GDAL many-part polygon warning, s2 planar
  # notes); real errors still surface. Progress messages sit outside it.
  hush <- function(expr) {
    if (quiet_warnings) suppressWarnings(suppressMessages(expr)) else expr
  }

  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(TRUE), add = TRUE)

  # ---- 1. CONUS states + dissolved boundary -------------------------------
  if (verbose) message("Building CONUS boundary from TIGER/Line states...")

  non_conus <- c("Alaska", "Hawaii", "Puerto Rico", "American Samoa",
                 "Guam", "Northern Mariana Islands", "Virgin Islands",
                 "United States Virgin Islands",
                 "Commonwealth of the Northern Mariana Islands")

  conus_states <- hush(
    tigris::states(progress_bar = FALSE) %>%
      sf::st_transform(4326) %>%
      dplyr::select(STATEFP, GEOID, STUSPS, NAME) %>%
      dplyr::rename(state_name = NAME) %>%
      dplyr::filter(!state_name %in% non_conus)
  )

  conus_boundary <- hush(conus_states %>% sf::st_make_valid() %>% sf::st_union())

  # ---- 2. Water supplies, filtered to CONUS -------------------------------
  # Some water supplies have NA source_state; fall back to the "state" column.
  wsd_conus <- water_supplies %>%
    dplyr::mutate(
      source_state = dplyr::coalesce(source_state, state)
    ) %>%
    dplyr::filter(source_state %in% conus_states$STUSPS) %>%
    dplyr::distinct(wsd_source_id, .keep_all = TRUE) %>%
    dplyr::left_join(
      conus_states %>%
        sf::st_drop_geometry() %>%
        dplyr::select(source_state = STUSPS, state_name),
      by = "source_state"
    ) %>%
    dplyr::select(
      wsd_source_id, park_unit, park_name, region,
      source_state, state_name, water_system_name, system_type,
      source_longitude, source_latitude
    )
  
  if (verbose) {
    message("  ", nrow(wsd_conus), " CONUS water systems")
  }

  # ---- 3. Aqueduct annual layers, clipped to CONUS ------------------------
  stress_levels <- c("Low (<10%)", "Low - Medium (10-20%)",
                     "Medium - High (20-40%)", "High (40-80%)",
                     "Extremely High (>80%)", "Arid and Low Water Use")

  # same categories as spelled in the future_annual layer
  stress_levels_typo <- c("Low (<10%)", "Low-medium (10-20%)",
                          "Medium-high (20-40%)", "High (40-80%)",
                          "Extremely high (>80%)", "Arid and low water use")

  # severity rank: Arid ties with Extremely High (both score 5)
  stress_rank <- c("Low (<10%)" = 1, "Low - Medium (10-20%)" = 2,
                   "Medium - High (20-40%)" = 3, "High (40-80%)" = 4,
                   "Extremely High (>80%)" = 5, "Arid and Low Water Use" = 5)

  # Read + clip to CONUS (s2 off: Aqueduct polygons are large and often invalid)
  read_conus <- function(layer) {
    if (verbose) message("  reading layer: ", layer)

    out <- hush(
      sf::st_read(gdb_path, layer, quiet = TRUE) %>%
        sf::st_make_valid() %>%
        sf::st_transform(sf::st_crs(conus_boundary)) %>%
        sf::st_filter(conus_boundary, .predicate = sf::st_intersects)
    )

    if (verbose) message("    ", nrow(out), " sub-basins after clip")

    out
  }

  df_baseline <- read_conus("baseline_annual") %>%
    dplyr::select(aq30_id, pfaf_id, name_1, bws_raw, bws_score, bws_label) %>%
    dplyr::filter(pfaf_id != -9999, bws_score != -9999) %>%  # filter out "No Data" [-9999]
    dplyr::mutate(
      bws_label = factor(bws_label, levels = stress_levels, ordered = TRUE),
      dplyr::across(bws_raw, ~ round(.x, 3))
    ) %>%
    dplyr::rename_with(~ gsub("^bws_", "baseline_stress_", .x), dplyr::starts_with("bws_"))

  df_future <- read_conus("future_annual") %>%
    dplyr::select(pfaf_id, dplyr::matches(c("ws_x_r$", "ws_x_s$", "ws_x_l$"))) %>%
    dplyr::filter(dplyr::if_all(dplyr::matches("ws_x_r$"), ~ !is.na(.x))) %>%
    dplyr::mutate(
      dplyr::across(dplyr::matches("ws_x_l$"),
             ~ factor(.x, levels = stress_levels_typo,
                      labels = stress_levels, ordered = TRUE)),
      dplyr::across(dplyr::matches("_x_r$"), ~ round(.x, 3))
    ) %>%
    # bau/opt/pes -> rcp70/rcp26/rcp85; suffixes -> raw/score/label
    dplyr::rename_with(~ .x %>%
                  gsub("^bau", "rcp70_20", .) %>%
                  gsub("^opt", "rcp26_20", .) %>%
                  gsub("^pes", "rcp85_20", .) %>%
                  gsub("_ws_x_", "_stress_", .) %>%
                  gsub("_r$", "_raw", .) %>%
                  gsub("_s$", "_score", .) %>%
                  gsub("_l$", "_label", .))

  # ---- 4. Attach each water system to its sub-basin -----------------------
  # Separate nearest-feature joins: the two layers are at different HydroBASINS
  # resolutions, so a single join would mis-assign the finer baseline values.
  if (verbose) message("Joining water systems to Aqueduct sub-basins...")

  NPSWaterSystems <- hush(
    wsd_conus %>%
      sf::st_as_sf() %>%
      sf::st_make_valid() %>%
      sf::st_transform(sf::st_crs(df_baseline))
  )

  ws_baseline <- hush(sf::st_join(NPSWaterSystems, df_baseline,
                              join = sf::st_nearest_feature, left = TRUE))
  ws_future   <- hush(sf::st_join(NPSWaterSystems, df_future,
                              join = sf::st_nearest_feature, left = TRUE))

  future_attrs <- ws_future %>%
    sf::st_drop_geometry() %>%
    dplyr::select(wsd_source_id,
           pfaf_id_future = pfaf_id,
           dplyr::matches("^rcp\\d{2}_\\d{4}_stress_(raw|score|label)$"))

  combined <- ws_baseline %>%
    sf::st_drop_geometry() %>%
    dplyr::rename(pfaf_id_baseline = pfaf_id) %>%
    dplyr::left_join(future_attrs, by = "wsd_source_id")

  # ---- 5. Long form: one row per system x rcp x period --------------------
  stress_long <- combined %>%
    # give baseline a matching two-part prefix so the pivot regex catches it too
    dplyr::rename_with(~ sub("^baseline_stress_", "baseline_baseline_stress_", .x),
                .cols = dplyr::matches("^baseline_stress_")) %>%
    tidyr::pivot_longer(
      cols = dplyr::matches("_stress_(raw|score|label)$"),
      names_pattern = "^(rcp\\d{2}|baseline)_(\\d{4}|baseline)_stress_(raw|score|label)$",
      names_to = c("rcp", "period", ".value")
    ) %>%
    dplyr::rename(stress_raw   = raw,
           stress_score = score,
           stress_label = label) %>%
    dplyr::mutate(
      stress_label = factor(as.character(stress_label),
                            levels = stress_levels, ordered = TRUE),
      stress_rank  = unname(stress_rank[as.character(stress_label)]),
      rcp          = factor(rcp, levels = c("baseline", "rcp26", "rcp70", "rcp85")),
      period       = factor(period, levels = c("baseline", "2030", "2050", "2080"))
    ) %>%
    dplyr::select(wsd_source_id, park_unit, park_name, region, source_state, state_name,
           water_system_name, system_type, source_longitude, source_latitude,
           aq30_id, pfaf_id_baseline, pfaf_id_future,
           rcp, period, stress_raw, stress_score, stress_label, stress_rank) %>%
    dplyr::arrange(wsd_source_id, period, rcp)

  if (verbose) {
    message("  ", nrow(stress_long), " rows (",
            dplyr::n_distinct(stress_long$wsd_source_id), " water systems x ",
            dplyr::n_distinct(paste(stress_long$rcp, stress_long$period)),
            " scenarios)")
  }

  # ---- 6. Save ------------------------------------------------------------
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    stress_file <- file.path(output_dir, "waterrights_stress_long.csv")
    readr::write_csv(stress_long, stress_file)

    if (verbose) message("Saved: ", stress_file)

    if (save_intermediate) {
      baseline_file <- file.path(output_dir, "CONUS_aqueduct_baseline_annual.gpkg")
      future_file   <- file.path(output_dir, "CONUS_aqueduct_future_annual.gpkg")

      hush({
        sf::st_write(df_baseline, baseline_file, delete_dsn = TRUE, quiet = TRUE)
        sf::st_write(df_future,   future_file,   delete_dsn = TRUE, quiet = TRUE)
      })

      if (verbose) {
        message("Saved: ", baseline_file)
        message("Saved: ", future_file)
      }
    }
  }

  if (verbose) {
    elapsed <- Sys.time() - start_time
    message("\nCompleted in ", round(elapsed, 2), " ", attr(elapsed, "units"))
  }

  return(stress_long)
}
