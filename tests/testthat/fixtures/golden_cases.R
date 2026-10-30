# Cases for the golden comparison between the old code in
# nps_water_vulnerability (before the split) and wavarapid.
#
# Used by make_golden.R, which writes golden/<id>.rds from the old code, and by
# test-golden.R, which checks that wavarapid returns the identical object.
# Each case has these fields:
# `id`: the golden file name.
# `fun`: the function to call, by name.
# `args`: a list of arguments, or a function of the saved inputs (`inputs`,
#   below) that returns the list. The arguments must work with both the old
#   and the new function.
# `needs`: "data" = the shared WaVA data folder, "internet" = a live download.
# `reads`: for a "data" case, the files or folders under the data folder it
#   reads; the case skips without them, and their checksums go in
#   golden_env.csv.
# `downloads`: what an "internet" case fetches, named in its skip reason.
# `keep`: optional function applied to the result before it is saved or
#   compared, for results identical() can't compare as they are.
# `stop_on`: optional pattern. A warning that matches it means a download
#   failed inside a function that carries on anyway, so golden_run() turns it
#   into an error: a service that is down is an error, not a pass or a
#   mismatch.
#
# Not covered (see README.md): process_groundwater_depth() and
# process_inundation_data().

# Paths in the shared data folder, for arguments with no default. make_golden.R
# sets options(wava.data_dir) to the old code's `data` link, so both runs read
# the same files.
golden_data_path <- function(...) wavacore::wava_data_path(...)

# The folder a case that writes files writes into, inside tempdir(). Nothing
# else may be written (see golden_run()).
golden_out_dir <- function(id) file.path(tempdir(), "wavarapid_golden", id)

# Every parquet file a case wrote, by file name, for its `keep`
golden_read_parquet <- function(id) {
  files <- sort(list.files(golden_out_dir(id), pattern = "\\.parquet$",
                           full.names = TRUE))
  stats::setNames(lapply(files, arrow::read_parquet), basename(files))
}

# The water supplies the cases run on: 8 supplies in 4 parks, chosen
# (2026-10-27) so the rarer branches get real values.
golden_site_ids <- c(
  "ROMO_WS_002_01", # river/stream, snowmelt-fed
  "ROMO_WS_005_01", # shallow groundwater with a well depth (114 ft)
  "ZION_WS_003_01", # spring with no well depth (water table depth fallback)
  "ZION_WS_005_01", # deep groundwater, 925 ft well
  "PORE_WS_004_01", # river/stream on the Pacific coast, 2.8 km from shore
  "PORE_WS_005_01", # deep groundwater 14 m from the shore
  "CAHA_WS_001_01", # Atlantic barrier island, 64 m from the shore
  "CAHA_WS_002_01"  # shallow groundwater with a well depth, 0.9 km from shore
)

# Inputs, made once by make_golden.R with the OLD code and saved in
# golden/inputs.rds; test-golden.R reads them from there. So both runs get the
# identical objects, and the test doesn't call get_wsd() (wavacore's golden
# tests cover it). Each one is a function of the inputs made before it.
golden_inputs <- list(
  # Every rapid water supply. The water rights functions need all of them:
  # they stop unless water_supplies has exactly the supplies in the stress
  # table, which process_waterrightsindicator_data() made for the full set.
  ws_rapid = function(inputs) get_wsd(rapid = TRUE),
  ws = function(inputs) {
    inputs$ws_rapid[inputs$ws_rapid$wsd_source_id %in% golden_site_ids, ]
  },
  # As in rapid_workflow.Rmd, for the fire sensitivity
  ws_buffer = function(inputs) {
    sf::st_buffer(sf::st_transform(inputs$ws, crs = 32616), dist = 4000)
  },
  # As in rapid_workflow.Rmd: the 4 km inundation scenarios, with every water
  # supply joined on (those without a scenario get proj_file NA)
  ws_inundation = function(inputs) {
    scenarios <- sf::st_read(
      golden_data_path("rapid", "inputs", "coast_comp",
                       "buff_ws_inun_regions.gpkg"),
      quiet = TRUE
    )
    scenarios <- dplyr::filter(scenarios, !is.na(proj_file))
    dplyr::right_join(
      dplyr::select(scenarios, -park_unit),
      dplyr::select(sf::st_drop_geometry(inputs$ws), wsd_source_id, park_unit),
      by = "wsd_source_id"
    )
  },
  pyromes = function(inputs) {
    sf::st_read(golden_data_path("rapid", "inputs", "pyromes.gpkg"),
                quiet = TRUE)
  },
  fire_prob_data = function(inputs) {
    readr::read_csv(golden_data_path("rapid", "inputs",
                                     "fire_probability_change.csv"),
                    show_col_types = FALSE)
  },
  # make_fake_supplies() is in helper-fixtures.R. Every source type, since
  # rainwater and ocean have their own branches in the index.
  vuln_dat = function(inputs) {
    make_fake_supplies(8, source_type = c(
      "deep groundwater", "shallow groundwater", "unknown groundwater",
      "spring", "river/stream", "lake/reservoir", "rainwater", "ocean"
    ))
  },
  vuln_dat_na = function(inputs) {
    make_fake_supplies(6, na_cols = c("exp_swi", "sen_flood_risk"))
  }
)

# What the inputs read in the shared data folder, for golden_env.csv
golden_input_files <- c(
  file.path("Water_Supply_Systems", "NPS_Water_Systems_Database_Joined.xlsx"),
  file.path("rapid", "inputs", "coast_comp", "buff_ws_inun_regions.gpkg"),
  file.path("rapid", "inputs", "pyromes.gpkg"),
  file.path("rapid", "inputs", "fire_probability_change.csv")
)

# The cases. Built in local() so the shorthands below don't land in the old
# code's global environment, where make_golden.R sources this file.
golden_cases <- local({
  historical <- c(1980, 2010)
  future <- c(2035, 2065)

  # Shorthands for `reads`
  inputs_dir <- function(...) file.path("rapid", "inputs", ...)
  coast_csv <- inputs_dir("coast_comp", "coastal_site_geometry.csv")
  precip_dir <- inputs_dir("loca2_precipitation_water_supplies")
  pet_dir <- inputs_dir("loca2_pet_water_supplies")
  runoff_dir <- inputs_dir("loca2_runoff_water_supplies")
  swe_dir <- inputs_dir("loca2_swe_water_supplies")
  temp_dir <- inputs_dir("loca2_temperature_water_supplies")
  periods_dir <- inputs_dir("loca_cmip6_precip_periods_water_supplies")
  stress_csv <- inputs_dir("WaterRightsIndicator",
                           "waterrights_stress_long.csv")
  demand_dir <- file.path("Latest county-level data",
                          "water_demand_projections")
  wbd_gdb <- file.path("all", "WBD_National_GDB.gdb")
  storm_surge_gpkg <- inputs_dir("storm_surge",
                                 "storm_surge_RP10_muis_2023.gpkg")
  whp_dir <- file.path("all", "usfs_fire", "whp2023_GeoTIF")
  # The inundation tiles the two coastal parks use (PORE: CA_SFBay, CAHA:
  # NC_Middle2), not the whole 3.7 GB folder
  inundation_tiles <- inputs_dir("coast_comp", "Extent_Rasters", c(
    "CA/CA_SFBay_connectRaster_0_0.tif", "CA/CA_SFBay_connectRaster_2_0.tif",
    "NC/NC_Middle2_connectRaster_0_0.tif", "NC/NC_Middle2_connectRaster_2_5.tif"
  ))
  loca2_mwbm_nc <- file.path(
    "all", "loca2_mwbm", "snow_MWBM_LOCA2.CanESM5.ssp585.r1i1p1f1_2015-2100.nc"
  )
  # The training (Jay Alder) file of a LOCA2 variable, e.g. "jay_alder_pr"
  training_file <- function(dir, prefix) {
    file.path(dir, paste0(prefix,
                          "_MWBM_LOCA2_training_195001-201412.nc.parquet"))
  }
  visitation_download <- paste("NPS visitation for ROMO, ZION, PORE and CAHA",
                               "(irmaservices.nps.gov)")
  pr_metrics_nc <- file.path(
    "all", "loca_cmip6_precip_products", "historical",
    "NorESM2-MM.historical.r1i1p1f1.1950-2014.LOCA_16thdeg_pr_metrics.nc"
  )

  # A small table for add_standardized_indicators(): a column with an NA, an
  # all-NA column (its own branch) and a column the prefix doesn't match
  raw_table <- data.frame(
    wsd_source_id = c("A", "B", "C", "D", "E"),
    raw_a = c(1, NA, 3, 10, 5),
    raw_b = NA_real_,
    other = 1:5
  )

  list(
    # No data, no internet
    list(id = "minmax_norm",
         fun = "minmax_norm",
         args = list(x = c(3, NA, 7, 1, 7))),
    list(id = "minmax_norm_keep_na",
         fun = "minmax_norm",
         args = list(x = c(3, NA, 7, 1, 7), na.rm = FALSE)),
    list(id = "minmax_norm_constant",
         fun = "minmax_norm",
         args = list(x = c(2, 2, NA))),
    list(id = "add_standardized_indicators",
         fun = "add_standardized_indicators",
         args = list(data = raw_table, prefix = "raw_")),
    list(id = "add_standardized_indicators_invert",
         fun = "add_standardized_indicators",
         args = list(data = raw_table, prefix = "raw_", invert = TRUE)),
    list(id = "calc_vulnerability_index",
         fun = "calc_vulnerability_index",
         args = function(inputs) list(dat = inputs$vuln_dat)),
    list(id = "calc_vulnerability_index_na_cols",
         fun = "calc_vulnerability_index",
         args = function(inputs) list(dat = inputs$vuln_dat_na)),

    # Saved inputs only
    list(id = "calc_fire_exposure",
         fun = "calc_fire_exposure",
         args = function(inputs) list(water_supplies = inputs$ws,
                                      fire_prob_data = inputs$fire_prob_data,
                                      pyrome_data = inputs$pyromes)),

    # Internet: FEMA flood zones. calc_flood_risk() turns a failed request
    # into NaN for that supply with only a GDAL warning (issues.md Q39), so
    # `stop_on` makes that warning an error. These come early, so a FEMA
    # outage stops make_golden.R before the long cases.
    list(id = "calc_flood_risk_4km",
         fun = "calc_flood_risk",
         args = function(inputs) list(data = inputs$ws, dist = 4000),
         needs = "internet",
         downloads = "FEMA flood zones (hazards.fema.gov NFHL)",
         stop_on = "GDAL Error"),
    # The default dist = 0 buffers to empty polygons, so every supply gets NaN
    # (todo.md); kept to show the default behaves the same
    list(id = "calc_flood_risk_default",
         fun = "calc_flood_risk",
         args = function(inputs) list(data = inputs$ws),
         needs = "internet",
         downloads = "FEMA flood zones (hazards.fema.gov NFHL)",
         stop_on = "GDAL Error"),

    # Shared data folder: small files
    list(id = "calc_fire_sensitivity_cat",
         fun = "calc_fire_sensitivity",
         args = function(inputs) list(
           water_supplies = inputs$ws_buffer,
           whp_raster = golden_data_path(whp_dir, "whp2023_cls_conus.tif"),
           whp_type = "cat"
         ),
         needs = "data",
         reads = file.path(whp_dir, "whp2023_cls_conus.tif")),
    list(id = "calc_fire_sensitivity_raw",
         fun = "calc_fire_sensitivity",
         args = function(inputs) list(
           water_supplies = inputs$ws_buffer,
           whp_raster = golden_data_path(whp_dir, "whp2023_cnt_conus.tif"),
           whp_type = "raw"
         ),
         needs = "data",
         reads = file.path(whp_dir, "whp2023_cnt_conus.tif")),
    list(id = "calc_groundwater_depth_sensitivity",
         fun = "calc_groundwater_depth_sensitivity",
         args = function(inputs) list(
           water_supplies = inputs$ws,
           groundwater_data = golden_data_path(inputs_dir("groundwater_depth",
                                                          "wtd_all_sites.csv"))
         ),
         needs = "data",
         reads = inputs_dir("groundwater_depth", "wtd_all_sites.csv")),
    list(id = "calc_storm_surge_exposure",
         fun = "calc_storm_surge_exposure",
         args = function(inputs) list(
           water_supplies = inputs$ws,
           storm_surge_stations = golden_data_path(storm_surge_gpkg)
         ),
         needs = "data",
         reads = c(storm_surge_gpkg, coast_csv)),
    list(id = "calc_storm_surge_sensitivity",
         fun = "calc_storm_surge_sensitivity",
         args = function(inputs) list(
           water_supplies = inputs$ws,
           storm_surge_stations = golden_data_path(storm_surge_gpkg)
         ),
         needs = "data",
         reads = c(storm_surge_gpkg, coast_csv)),
    list(id = "calc_swi_exposure",
         fun = "calc_swi_exposure",
         args = function(inputs) list(water_supplies = inputs$ws),
         needs = "data",
         reads = c(inputs_dir("swi", "swi_points.gpkg"), coast_csv)),
    # All rapid water supplies (see ws_rapid above)
    list(id = "calc_waterrights_exposure_2030",
         fun = "calc_waterrights_exposure",
         args = function(inputs) list(water_supplies = inputs$ws_rapid,
                                      period = "2030"),
         needs = "data",
         reads = stress_csv),
    list(id = "calc_waterrights_exposure_2050",
         fun = "calc_waterrights_exposure",
         args = function(inputs) list(water_supplies = inputs$ws_rapid),
         needs = "data",
         reads = stress_csv),
    list(id = "calc_waterrights_exposure_2080",
         fun = "calc_waterrights_exposure",
         args = function(inputs) list(water_supplies = inputs$ws_rapid,
                                      period = "2080"),
         needs = "data",
         reads = stress_csv),
    list(id = "calc_waterrights_sensitivity",
         fun = "calc_waterrights_sensitivity",
         args = function(inputs) list(water_supplies = inputs$ws_rapid),
         needs = "data",
         reads = stress_csv),
    list(id = "calc_dry_period_exposure",
         fun = "calc_dry_period_exposure",
         args = function(inputs) list(water_supplies = inputs$ws,
                                      historical = historical, future = future),
         needs = "data",
         reads = periods_dir),
    list(id = "calc_dry_period_sensitivity",
         fun = "calc_dry_period_sensitivity",
         args = function(inputs) list(water_supplies = inputs$ws,
                                      historical = historical),
         needs = "data",
         reads = periods_dir),
    list(id = "calc_wet_period_exposure",
         fun = "calc_wet_period_exposure",
         args = function(inputs) list(water_supplies = inputs$ws,
                                      historical = historical, future = future),
         needs = "data",
         reads = periods_dir),
    list(id = "calc_wet_period_sensitivity",
         fun = "calc_wet_period_sensitivity",
         args = function(inputs) list(water_supplies = inputs$ws,
                                      historical = historical),
         needs = "data",
         reads = periods_dir),
    list(id = "calc_drought_exposure",
         fun = "calc_drought_exposure",
         args = function(inputs) list(water_supplies = inputs$ws,
                                      historical = historical, future = future),
         needs = "data",
         reads = inputs_dir(
           "loca2_drought_monthly_1980-2010_2035-2065_hindcast-TRUE.parquet"
         )),
    # The sensitivities read only the training (Jay Alder) file of each variable
    list(id = "calc_drought_sensitivity",
         fun = "calc_drought_sensitivity",
         args = function(inputs) list(water_supplies = inputs$ws,
                                      historical = historical),
         needs = "data",
         reads = c(training_file(precip_dir, "jay_alder_pr"),
                   training_file(pet_dir, "Jay_Alder_pet"))),
    list(id = "calc_precip_sensitivity",
         fun = "calc_precip_sensitivity",
         args = function(inputs) list(water_supplies = inputs$ws,
                                      historical = historical),
         needs = "data",
         reads = training_file(precip_dir, "jay_alder_pr")),
    list(id = "calc_runoff_sensitivity",
         fun = "calc_runoff_sensitivity",
         args = function(inputs) list(water_supplies = inputs$ws,
                                      historical = historical),
         needs = "data",
         reads = training_file(runoff_dir, "Jay_Alder_runoff")),
    list(id = "calc_swe_sensitivity",
         fun = "calc_swe_sensitivity",
         args = function(inputs) list(water_supplies = inputs$ws,
                                      historical = historical),
         needs = "data",
         reads = training_file(swe_dir, "Jay_Alder_snow")),
    list(id = "calc_temp_sensitivity",
         fun = "calc_temp_sensitivity",
         args = function(inputs) list(water_supplies = inputs$ws,
                                      historical = historical),
         needs = "data",
         reads = training_file(temp_dir, "jay_alder_tas")),
    list(id = "calc_inundation_exposure",
         fun = "calc_inundation_exposure",
         args = function(inputs) list(ws_scenarios = inputs$ws_inundation),
         needs = "data",
         reads = c(inundation_tiles, coast_csv)),
    list(id = "process_swi_data",
         fun = "process_swi_data",
         args = list(save = FALSE),
         needs = "data",
         reads = file.path("all", "coast_comp_data", "adams_swi_2024", "data",
                           "output_data.txt")),
    list(id = "process_storm_surge_data",
         fun = "process_storm_surge_data",
         args = function(inputs) list(
           data_dir = golden_data_path("all", "coast_comp_data",
                                       "storm_surge_muis_2023_RP10"),
           save = FALSE
         ),
         needs = "data",
         reads = file.path("all", "coast_comp_data",
                           "storm_surge_muis_2023_RP10")),
    # process_loca2_mwbm() returns nothing; its output is the parquet it writes
    list(id = "process_loca2_mwbm",
         fun = "process_loca2_mwbm",
         args = function(inputs) list(
           nc_path = golden_data_path(loca2_mwbm_nc),
           water_supplies = inputs$ws,
           var = "snow",
           output_dir = golden_out_dir("process_loca2_mwbm")
         ),
         needs = "data",
         reads = loca2_mwbm_nc,
         keep = function(x) {
           list(value = x, files = golden_read_parquet("process_loca2_mwbm"))
         }),
    list(id = "process_loca2_mwbm_nearest",
         fun = "process_loca2_mwbm",
         args = function(inputs) list(
           nc_path = golden_data_path(loca2_mwbm_nc),
           water_supplies = inputs$ws,
           var = "snow",
           output_dir = golden_out_dir("process_loca2_mwbm_nearest"),
           grab_nearest = TRUE
         ),
         needs = "data",
         reads = loca2_mwbm_nc,
         keep = function(x) {
           list(value = x,
                files = golden_read_parquet("process_loca2_mwbm_nearest"))
         }),
    # process_precip_period_data() returns the path of the parquet it writes,
    # which is in tempdir(); keep the file's name and contents instead
    list(id = "process_precip_period_data",
         fun = "process_precip_period_data",
         args = function(inputs) list(
           nc_path = golden_data_path(pr_metrics_nc),
           water_supplies = inputs$ws,
           output_dir = golden_out_dir("process_precip_period_data"),
           overwrite = TRUE
         ),
         needs = "data",
         reads = pr_metrics_nc,
         keep = function(x) {
           list(file = basename(x), data = arrow::read_parquet(x))
         }),

    # Shared data folder: the LOCA2 files of every model (several GB each)
    list(id = "calc_swe_exposure",
         fun = "calc_swe_exposure",
         args = function(inputs) list(water_supplies = inputs$ws,
                                      historical = historical, future = future),
         needs = "data",
         reads = swe_dir),
    list(id = "calc_precip_exposure",
         fun = "calc_precip_exposure",
         args = function(inputs) list(water_supplies = inputs$ws,
                                      historical = historical, future = future),
         needs = "data",
         reads = precip_dir),
    list(id = "calc_precip_timing_exposure",
         fun = "calc_precip_timing_exposure",
         args = function(inputs) list(water_supplies = inputs$ws,
                                      historical = historical, future = future),
         needs = "data",
         reads = precip_dir),
    list(id = "calc_runoff_exposure",
         fun = "calc_runoff_exposure",
         args = function(inputs) list(water_supplies = inputs$ws,
                                      historical = historical, future = future),
         needs = "data",
         reads = runoff_dir),
    list(id = "calc_runoff_timing_exposure",
         fun = "calc_runoff_timing_exposure",
         args = function(inputs) list(water_supplies = inputs$ws,
                                      historical = historical, future = future),
         needs = "data",
         reads = runoff_dir),
    list(id = "calc_temp_exposure",
         fun = "calc_temp_exposure",
         args = function(inputs) list(water_supplies = inputs$ws,
                                      historical = historical, future = future),
         needs = "data",
         reads = temp_dir),
    # One variable for each way the summary aggregates a year: snow, a mean
    # (temp), and a sum (runoff)
    list(id = "calc_loca2_change_summary",
         fun = "calc_loca2_change_summary",
         args = function(inputs) list(water_supplies = inputs$ws,
                                      variables = c("runoff", "snow", "temp")),
         needs = "data",
         reads = c(runoff_dir, swe_dir, temp_dir)),
    list(id = "process_drought_data",
         fun = "process_drought_data",
         args = function(inputs) list(water_supplies = inputs$ws,
                                      historical = historical, future = future,
                                      save = FALSE),
         needs = "data",
         reads = c(precip_dir, pet_dir)),

    # Shared data folder and internet
    list(id = "calc_nearby_use_exposure",
         fun = "calc_nearby_use_exposure",
         args = function(inputs) list(water_supplies = inputs$ws),
         needs = c("data", "internet"),
         reads = demand_dir,
         downloads = "the 2020 US county boundaries (tigris)"),
    list(id = "calc_nearby_use_sensitivity",
         fun = "calc_nearby_use_sensitivity",
         args = function(inputs) list(water_supplies = inputs$ws),
         needs = c("data", "internet"),
         reads = demand_dir,
         downloads = "the 2020 US county boundaries (tigris)"),
    list(id = "calc_wta_exposure",
         fun = "calc_wta_exposure",
         args = function(inputs) list(water_supplies = inputs$ws),
         needs = "data",
         reads = c(demand_dir, wbd_gdb)),
    list(id = "calc_wta_sensitivity",
         fun = "calc_wta_sensitivity",
         args = function(inputs) list(water_supplies = inputs$ws),
         needs = "data",
         reads = c(demand_dir, wbd_gdb)),
    list(id = "process_waterrightsindicator_data",
         fun = "process_waterrightsindicator_data",
         args = function(inputs) list(water_supplies = inputs$ws, save = FALSE),
         needs = c("data", "internet"),
         reads = file.path("all", "aqueduct_v4",
                           "Aqueduct40_waterrisk_download_Y2023M07D05", "GDB",
                           "Aq40_Y2023D07M05.gdb"),
         downloads = "the US state boundaries (tigris)"),
    # save = FALSE: with save = TRUE it also writes into the CUSP input folder
    list(id = "process_coastal_geometry",
         fun = "process_coastal_geometry",
         args = function(inputs) list(water_supplies = inputs$ws, save = FALSE),
         needs = c("data", "internet"),
         reads = file.path("all", "coast_comp_data", "NOAA_CUSP"),
         downloads = "point elevations from the USGS EPQS service (elevatr)"),

    # Internet only
    list(id = "calc_visitation_exposure",
         fun = "calc_visitation_exposure",
         args = function(inputs) list(water_supplies = inputs$ws),
         needs = "internet",
         downloads = visitation_download),
    list(id = "calc_visitation_sensitivity",
         fun = "calc_visitation_sensitivity",
         args = function(inputs) list(water_supplies = inputs$ws,
                                      start_year = 1980, end_year = 2025),
         needs = "internet",
         downloads = visitation_download)
  )
})

# Run one case: build its arguments, make sure it can only write into
# golden_out_dir(id), call it, and apply its `keep`. Used by make_golden.R
# (old code) and test-golden.R (wavarapid).
golden_run <- function(case, inputs) {
  args <- if (is.function(case$args)) case$args(inputs) else case$args
  # A process_* function writes unless save = FALSE (their defaults write into
  # the shared data folder); so does save = TRUE anywhere
  writes <- isTRUE(args$save) ||
    (startsWith(case$fun, "process_") && !isFALSE(args$save))
  if (writes) {
    out_dir <- golden_out_dir(case$id)
    if (!identical(args$output_dir, out_dir)) {
      stop("golden case ", case$id, " writes files, so its output_dir must ",
           "be golden_out_dir(\"", case$id, "\")", call. = FALSE)
    }
    unlink(out_dir, recursive = TRUE)
    dir.create(out_dir, recursive = TRUE)
  }
  value <- if (is.null(case$stop_on)) {
    do.call(case$fun, args)
  } else {
    withCallingHandlers(
      do.call(case$fun, args),
      warning = function(w) {
        if (grepl(case$stop_on, conditionMessage(w))) {
          stop("golden case ", case$id, ": a download failed (",
               conditionMessage(w), "), so there is nothing to compare. ",
               "The service is down or unreliable; try again later.",
               call. = FALSE)
        }
      }
    )
  }
  if (!is.null(case$keep)) value <- case$keep(value)
  value
}

# Files and folders in the shared data folder that the "data" cases read.
# make_golden.R records their checksums, and test-golden.R reports it when one
# has changed.
golden_data_files <- unique(unlist(lapply(golden_cases, function(case) {
  case$reads
})))

# md5 of a file, or one md5 over every file in a folder (names and contents)
golden_checksum <- function(path) {
  if (!file.exists(path)) return("missing")
  if (!dir.exists(path)) return(unname(tools::md5sum(path)))
  files <- sort(list.files(path, recursive = TRUE, all.files = TRUE))
  sums <- tools::md5sum(file.path(path, files))
  unname(tools::md5sum(bytes = charToRaw(paste(files, sums, collapse = "\n"))))
}

# Check the fields of every case, so a typo (e.g. `need = "data"`) stops here
# instead of silently dropping a skip.
local({
  fields <- c("id", "fun", "args", "needs", "reads", "downloads", "keep",
              "stop_on")
  for (case in golden_cases) {
    id <- if (is.character(case$id) && length(case$id) == 1) case$id else "?"
    bad <- function(...) stop("golden case ", id, ": ", ..., call. = FALSE)
    extra <- setdiff(names(case), fields)
    if (length(extra) > 0) bad("unknown field(s) ", toString(extra))
    if (id == "?") bad("`id` must be one string")
    if (!is.character(case$fun) || length(case$fun) != 1) {
      bad("`fun` must be one function name")
    }
    if (!is.list(case$args) && !is.function(case$args)) {
      bad("`args` must be a list or a function of the inputs")
    }
    if (!is.null(case$needs) && !all(case$needs %in% c("data", "internet"))) {
      bad("`needs` must be \"data\", \"internet\" or both")
    }
    if ("data" %in% case$needs &&
        !(is.character(case$reads) && length(case$reads) >= 1)) {
      bad("a \"data\" case needs `reads` (what it reads in the data folder)")
    }
    if (!"data" %in% case$needs && !is.null(case$reads)) {
      bad("`reads` is only for \"data\" cases")
    }
    if ("internet" %in% case$needs &&
        !(is.character(case$downloads) && length(case$downloads) == 1)) {
      bad("an \"internet\" case needs `downloads` (what it fetches)")
    }
    if (!is.null(case$keep) && !is.function(case$keep)) {
      bad("`keep` must be a function")
    }
    if (!is.null(case$stop_on) &&
        !(is.character(case$stop_on) && length(case$stop_on) == 1)) {
      bad("`stop_on` must be one pattern")
    }
  }
  ids <- vapply(golden_cases, function(case) case$id, character(1))
  if (anyDuplicated(ids)) {
    stop("golden case ids must be unique: ", toString(ids[duplicated(ids)]),
         call. = FALSE)
  }
})
