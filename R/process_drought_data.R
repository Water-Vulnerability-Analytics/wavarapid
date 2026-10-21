#' @title Build the monthly precipitation and PET table for drought (SPEI)
#'
#' @description
#' Reads the processed LOCA2 precipitation and PET parquet files (from
#' [process_loca2_mwbm()]), averages the grid cells to one value per water
#' supply, month and model run, joins precipitation and PET, and computes the
#' climatic water balance (precipitation minus PET). The result is the input
#' for [calc_drought_exposure()].
#'
#' Each historical and future window is padded with about two years of
#' earlier months, more than a 12-month SPEI needs, so
#' [calc_drought_exposure()] can trim to the accumulation scale it uses.
#'
#' Within each model, scenario and run, months from the primary historical
#' source come first; the matching SSP run fills historical months the
#' primary source doesn't have.
#'
#' @param water_supplies sf object. The water supplies, with a
#'   `wsd_source_id` column, e.g. from [wavacore::get_wsd()]. Only their IDs
#'   are used, to pick rows from the files. Required.
#' @param precip_data String. The folder of processed precipitation parquet
#'   files. Stops with an error if no files are left after the `hindcast`
#'   filter. Default: `rapid/inputs/loca2_precipitation_water_supplies` in
#'   the shared data folder, see [wavacore::wava_data_path()].
#' @param pet_data String. The folder of processed PET parquet files, used
#'   the same way. Default: `rapid/inputs/loca2_pet_water_supplies` in the
#'   shared data folder.
#' @param hindcast Logical. Where the historical months come from.
#'   * `TRUE`: each model's own historical run (training and Livneh files are
#'     left out). `source_type` is `"gcm_historical"`.
#'   * `FALSE`: the training run, copied to every model, scenario and run
#'     (only precipitation files whose names contain 2015, 2045, 2075 or
#'     "training", and PET files with 2015-2100 or "training", are read).
#'     `source_type` is `"training"`. Stops with an error if there are no
#'     training rows.
#'
#'   Default: `TRUE`.
#' @param historical Character or numeric vector of 2 years, e.g.
#'   `c("1980", "2010")`. The first and last year of the historical period;
#'   only the first 4 characters of each are used. Required.
#' @param future Character or numeric vector of 2 years, e.g.
#'   `c("2035", "2065")`. The first and last year of the future period.
#'   Required.
#' @param water_year Logical. How months are grouped into years.
#'   * `TRUE`: water years, October to September, labelled by the year they
#'     end in. The padding starts in March, two years before the first year.
#'   * `FALSE`: calendar years. The padding starts in June, two years before
#'     the first year.
#'
#'   Default: `TRUE`.
#' @param output_dir String. The folder for the parquet file when
#'   `save = TRUE`. Default: `rapid/inputs` in the shared data folder, which
#'   is where [calc_drought_exposure()] looks by default.
#' @param save Logical. Whether to write the table.
#'   * `TRUE`: writes
#'     `loca2_drought_monthly_<historical>_<future>_hindcast-<TRUE|FALSE>.parquet`
#'     (e.g. `loca2_drought_monthly_1980-2010_2035-2065_hindcast-TRUE.parquet`)
#'     to `output_dir`, creating the folder if needed.
#'   * `FALSE`: writes nothing.
#'
#'   Default: `TRUE`.
#' @param verbose Logical. Whether to print progress.
#'   * `TRUE`: prints each file as it is read, the saved path and the run
#'     time.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return A data.table, one row per water supply, model, scenario, run,
#'   period and month, with columns `wsd_source_id`, `park_unit`, `ym` (first
#'   day of the month), `year`, `month`, `analysis_year` (water or calendar
#'   year), `drought_period` (`"historical"` or `"future"`), `input_window`
#'   (`"reporting"` or `"padding"`), `reporting_period` (logical), `model`,
#'   `ssp`, `run`, `monthly_precipitation`, `monthly_pet`,
#'   `climatic_water_balance` (mm), `source_type`, `source_priority`,
#'   `source_model`, `source_ssp` and `source_run`.
#' @export
#' @seealso [calc_drought_exposure()], which calculates SPEI from this table;
#'   [process_loca2_mwbm()], which makes the parquet files it reads;
#'   [calc_drought_sensitivity()] for the historical drought indicator.
#' @examplesIf interactive() && dir.exists(wavacore::wava_data_path("rapid", "inputs", "loca2_pet_water_supplies"))
#' # Reads about 9 GB of processed precipitation and PET files, so it takes
#' # a long time
#' ws <- wavacore::get_wsd(park = "ZION")
#' drought_monthly <- process_drought_data(ws,
#'                                         historical = c("1980", "2010"),
#'                                         future = c("2035", "2065"),
#'                                         save = FALSE)
#' head(drought_monthly)

process_drought_data <- function(water_supplies,
                                 precip_data = wavacore::wava_data_path("rapid", "inputs", "loca2_precipitation_water_supplies/"),
                                 pet_data = wavacore::wava_data_path("rapid", "inputs", "loca2_pet_water_supplies/"),
                                 hindcast = TRUE,
                                 historical,
                                 future,
                                 water_year = TRUE,
                                 output_dir = wavacore::wava_data_path("rapid", "inputs/"),
                                 save = TRUE,
                                 verbose = TRUE) {
  
  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "water_supplies must include wsd_source_id" = "wsd_source_id" %in% names(water_supplies),
    "historical must be length 2" = length(historical) == 2,
    "future must be length 2" = length(future) == 2
  )
  
  if (verbose) {
    message("Starting drought data processing...")
    start_time <- Sys.time()
  }
  
  # Date boundaries
  hist_start_year   <- as.integer(substr(historical[1], 1, 4))
  hist_end_year     <- as.integer(substr(historical[2], 1, 4))
  future_start_year <- as.integer(substr(future[1], 1, 4))
  future_end_year   <- as.integer(substr(future[2], 1, 4))
  
  if (anyNA(c(hist_start_year, hist_end_year, future_start_year, future_end_year))) {
    stop("Could not parse years from historical/future inputs.")
  }
  
  # Padding is intentionally wider than SPEI-12 needs so calc_drought_exposure()
  # can subset more tightly based on user-selected spei_scale.
  if (water_year) {
    hist_input_start <- data.table::as.IDate(paste0(hist_start_year - 2L, "-03-01"))
    hist_input_end   <- data.table::as.IDate(paste0(hist_end_year, "-09-01"))
    fut_input_start  <- data.table::as.IDate(paste0(future_start_year - 2L, "-03-01"))
    fut_input_end    <- data.table::as.IDate(paste0(future_end_year, "-09-01"))
  } else {
    hist_input_start <- data.table::as.IDate(paste0(hist_start_year - 2L, "-06-01"))
    hist_input_end   <- data.table::as.IDate(paste0(hist_end_year, "-12-01"))
    fut_input_start  <- data.table::as.IDate(paste0(future_start_year - 2L, "-06-01"))
    fut_input_end    <- data.table::as.IDate(paste0(future_end_year, "-12-01"))
  }
  
  overall_start <- min(hist_input_start, fut_input_start)
  overall_end   <- max(hist_input_end, fut_input_end)
  
  ws_ids <- unique(water_supplies$wsd_source_id)
  
  precip_files <- list.files(precip_data, full.names = TRUE, pattern = "\\.parquet$")
  pet_files    <- list.files(pet_data, full.names = TRUE, pattern = "\\.parquet$")
  
  if (hindcast) {
    precip_files <- grep("training|livneh", precip_files, invert = TRUE, value = TRUE)
    pet_files    <- grep("training|livneh", pet_files, invert = TRUE, value = TRUE)
  } else {
    precip_files <- grep("2015|2045|2075|training", precip_files, value = TRUE)
    pet_files    <- grep("2015-2100|training", pet_files, value = TRUE)
  }
  
  if (!length(precip_files)) {
    stop("No processed precipitation files found.")
  }
  
  if (!length(pet_files)) {
    stop("No processed PET files found.")
  }
  
  if (verbose) {
    message("Processing ", length(precip_files), " precip + ", length(pet_files), " PET files")
  }
  
  # Parses model/ssp/run/time once per file. This avoids reading and parsing
  # the nc_file column for every row in massive parquet files.
  #
  # Expected examples:
  # pr.ACCESS-CM2.historical.r1i1p1f1.1950-2014.LOCA_16thdeg_v20240915.monthly.nc.parquet
  # pet_MWBM_LOCA2.ACCESS-CM2.historical.r1i1p1f1_1950-2014.nc.parquet
  # Jay_Alder_pr_MWBM_LOCA2_training_195001-201412.nc.parquet
  
  parse_file_meta <- function(files) {
    base <- basename(files)
    
    x <- sub("\\.parquet$", "", base)
    
    is_training <- grepl("training|livneh", x, ignore.case = TRUE)
    
    model <- rep(NA_character_, length(x))
    ssp   <- rep(NA_character_, length(x))
    run   <- rep(NA_character_, length(x))
    time  <- rep(NA_character_, length(x))
    
    model[is_training] <- "historical"
    ssp[is_training]   <- "training"
    run[is_training]   <- "training"
    time[is_training]  <- "1950-2100"
    
    idx <- which(!is_training)
    
    if (length(idx)) {
      y <- x[idx]
      
      y <- sub("^pr\\.", "", y)
      y <- sub("^pet_MWBM_LOCA2\\.", "", y)
      y <- sub("^pet\\.", "", y)
      
      parts <- data.table::tstrsplit(y, ".", fixed = TRUE)
      
      if (length(parts) < 3L) {
        stop("Could not parse model/ssp/run from filenames. Example bad name: ", y[1])
      }
      
      model[idx] <- parts[[1]]
      ssp[idx]   <- parts[[2]]
      
      run_token <- parts[[3]]
      has_under <- grepl("_", run_token, fixed = TRUE)
      
      run[idx] <- ifelse(
        has_under,
        sub("_.*$", "", run_token),
        run_token
      )
      
      time[idx] <- ifelse(
        has_under,
        sub("^.*_", "", run_token),
        if (length(parts) >= 4L) parts[[4]] else NA_character_
      )
      
      time[idx] <- sub("\\.nc$|\\.parquet$", "", time[idx])
    }
    
    meta <- data.table::data.table(
      file = files,
      model = model,
      ssp = ssp,
      run = run,
      time = time
    )
    
    if (anyNA(meta$model) || anyNA(meta$ssp) || anyNA(meta$run)) {
      print(meta[is.na(model) | is.na(ssp) | is.na(run)])
      stop("Some file metadata could not be parsed. Check filenames above.")
    }
    
    meta[]
  }
  
  # Process one climate variable at a time
  
  process_var <- function(files, value_col) {
    file_meta <- parse_file_meta(files)
    
    out <- vector("list", nrow(file_meta))
    
    for (i in seq_len(nrow(file_meta))) {
      f <- file_meta$file[i]
      
      if (verbose) {
        message(
          "  ", value_col, ": ",
          basename(f), " (", i, "/", nrow(file_meta), ")"
        )
      }
      
      # Read only the truly required columns. 
      dt <- arrow::open_dataset(f, format = "parquet") %>%
        dplyr::select(wsd_source_id,
                      date,
                      value,
                      park_unit) %>%
        dplyr::filter(wsd_source_id %in% ws_ids,
                      date >= as.character(overall_start),
                      date <= as.character(overall_end)) %>%
        dplyr::collect() %>%
        data.table::as.data.table() %>%
        .[, date := data.table::as.IDate(date)]
      
      if (!nrow(dt)) {
        out[[i]] <- NULL
        next
      }
      
      dt[, ym := date - (data.table::mday(date) - 1L)]
      
      dt <- dt[
        (ym >= hist_input_start & ym <= hist_input_end) |
          (ym >= fut_input_start & ym <= fut_input_end)
      ]
      
      if (!nrow(dt)) {
        out[[i]] <- NULL
        rm(dt)
        next
      }
      
      # Add file-level metadata as scalars, not row-parsed strings.
      dt[
        ,
        `:=`(
          model = file_meta$model[i],
          ssp   = file_meta$ssp[i],
          run   = file_meta$run[i]
        )
      ]
      
      dt[, year  := data.table::year(ym)]
      dt[, month := data.table::month(ym)]
      
      if (water_year) {
        dt[, analysis_year := year + data.table::fifelse(month >= 10L, 1L, 0L)]
      } else {
        dt[, analysis_year := year]
      }
      
      by_cols <- c("wsd_source_id",
                   "park_unit",
                   "ym",
                   "model",
                   "ssp",
                   "run",
                   "year",
                   "month",
                   "analysis_year")
      
      # Aggregate immediately to reduce memory pressure.
      ans <- dt[
        ,
        .(value = mean(value, na.rm = TRUE)),
        by = by_cols
      ]
      
      ans[!is.finite(value), value := NA_real_]
      data.table::setnames(ans, "value", value_col)
      
      out[[i]] <- ans
      
      rm(dt, ans)
    }
    
    keep <- vapply(
      out,
      function(x) !is.null(x) && data.table::is.data.table(x) && nrow(x) > 0L,
      logical(1)
    )
    
    out <- out[keep]
    
    if (!length(out)) {
      return(data.table::data.table())
    }
    
    out <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
    
    # Final collapse across files in case multiple source files overlap the same
    # site/month/model key.
    by_cols <- c("wsd_source_id",
                 "park_unit",
                 "ym",
                 "model",
                 "ssp",
                 "run",
                 "year",
                 "month",
                 "analysis_year")
    
    out <- out[
      ,
      .(value = mean(get(value_col), na.rm = TRUE)),
      by = by_cols
    ]
    
    out[!is.finite(value), value := NA_real_]
    data.table::setnames(out, "value", value_col)
    
    out[]
  }
  
  # Read and process precip + PET finals
  
  precip_monthly <- process_var(precip_files, "monthly_precipitation")
  pet_monthly    <- process_var(pet_files, "monthly_pet")
  
  if (!nrow(precip_monthly)) {
    stop("No precipitation rows remained after filtering.")
  }
  
  if (!nrow(pet_monthly)) {
    stop("No PET rows remained after filtering.")
  }
  
  # Join precip/PET
  
  join_keys <- c("wsd_source_id",
                 "park_unit",
                 "ym",
                 "model",
                 "ssp",
                 "run",
                 "year",
                 "month",
                 "analysis_year")
  
  data.table::setkeyv(precip_monthly, join_keys)
  data.table::setkeyv(pet_monthly, join_keys)
  
  climate <- merge(precip_monthly,
    pet_monthly,
    by = join_keys,
    all = TRUE)
  
  if (!nrow(climate)) {
    stop("No joinable precipitation/PET records were produced.")
  }
  
  climate[, climatic_water_balance := monthly_precipitation - monthly_pet]
  climate[!is.finite(climatic_water_balance), climatic_water_balance := NA_real_]
  
  # Build historical/future monthly rows
  scenario_combos <- unique(
    climate[
      !ssp %in% c("historical", "training") &
        !is.na(model) &
        !is.na(ssp) &
        !is.na(run),
      .(model, ssp, run)
    ]
  )
  
  if (!nrow(scenario_combos)) {
    print(climate[, .N, by = .(model, ssp, run)][order(model, ssp, run)])
    stop("No future model/SSP/run combinations found.")
  }
  
  scenario_combos[, combo_id := .I]
  
  if (hindcast) {
    # Model historical rows get attached to every matching future SSP combo for
    # the same model/run.
    hist_primary <- climate[
      ssp == "historical" &
        ym >= hist_input_start &
        ym <= hist_input_end
    ]
    
    hist_primary <- scenario_combos[
      hist_primary,
      on = .(model, run),
      allow.cartesian = TRUE,
      nomatch = 0L
    ]
    
    # data.table join keeps i.ssp from hist_primary if names collide.
    # We want output ssp to be the future scenario, but source_ssp to be historical.
    if ("i.ssp" %in% names(hist_primary)) {
      hist_primary[, source_ssp := i.ssp]
      hist_primary[, i.ssp := NULL]
    } else {
      hist_primary[, source_ssp := "historical"]
    }
    
    hist_primary[
      ,
      `:=`(
        source_model = model,
        source_run = run,
        drought_period = "historical",
        source_type = "gcm_historical",
        source_priority = 1L
      )
    ]
    
  } else {
    # Training rows get copied to every future scenario combo.
    hist_training <- climate[
      ssp == "training" &
        ym >= hist_input_start &
        ym <= hist_input_end
    ]
    
    if (!nrow(hist_training)) {
      stop("hindcast = FALSE, but no training rows were found in climate.")
    }
    
    # Cross join without bringing huge unnecessary columns.
    hist_training[, tmp_join_id := 1L]
    scenario_combos[, tmp_join_id := 1L]
    
    hist_primary <- scenario_combos[
      hist_training,
      on = "tmp_join_id",
      allow.cartesian = TRUE
    ]
    
    hist_primary[, tmp_join_id := NULL]
    scenario_combos[, tmp_join_id := NULL]
    hist_training[, tmp_join_id := NULL]
    
    # Normalize output model/ssp/run to scenario combo.
    if ("i.model" %in% names(hist_primary)) {
      hist_primary[, source_model := i.model]
      hist_primary[, i.model := NULL]
    } else {
      hist_primary[, source_model := "historical"]
    }
    
    if ("i.ssp" %in% names(hist_primary)) {
      hist_primary[, source_ssp := i.ssp]
      hist_primary[, i.ssp := NULL]
    } else {
      hist_primary[, source_ssp := "training"]
    }
    
    if ("i.run" %in% names(hist_primary)) {
      hist_primary[, source_run := i.run]
      hist_primary[, i.run := NULL]
    } else {
      hist_primary[, source_run := "training"]
    }
    
    hist_primary[
      ,
      `:=`(
        drought_period = "historical",
        source_type = "training",
        source_priority = 1L
      )
    ]
  }
  
  # Historical fill from matching SSP rows.
  hist_fill <- climate[
    !ssp %in% c("historical", "training") &
      ym >= hist_input_start &
      ym <= hist_input_end
  ]
  
  hist_fill <- scenario_combos[
    hist_fill,
    on = .(model, ssp, run),
    nomatch = 0L
  ]
  
  hist_fill[
    ,
    `:=`(
      source_model = model,
      source_ssp = ssp,
      source_run = run,
      drought_period = "historical",
      source_type = "ssp_historical_fill",
      source_priority = 2L
    )
  ]
  
  # Future rows from matching SSP rows.
  future_rows <- climate[
    !ssp %in% c("historical", "training") &
      ym >= fut_input_start &
      ym <= fut_input_end
  ]
  
  future_rows <- scenario_combos[
    future_rows,
    on = .(model, ssp, run),
    nomatch = 0L
  ]
  
  future_rows[
    ,
    `:=`(
      source_model = model,
      source_ssp = ssp,
      source_run = run,
      drought_period = "future",
      source_type = "ssp_future",
      source_priority = 1L
    )
  ]
  
  keep_cols <- c("wsd_source_id",
                 "park_unit",
                 "ym",
                 "year",
                 "month",
                 "analysis_year",
                 "drought_period",
                 "model",
                 "ssp",
                 "run",
                 "monthly_precipitation",
                 "monthly_pet",
                 "climatic_water_balance",
                 "source_type",
                 "source_priority",
                 "source_model",
                 "source_ssp",
                 "source_run")
  
  missing_keep <- setdiff(keep_cols, names(hist_primary))
  if (length(missing_keep)) {
    stop("hist_primary is missing expected columns: ", paste(missing_keep, collapse = ", "))
  }
  
  drought_monthly <- data.table::rbindlist(
    list(
      hist_primary[, ..keep_cols],
      hist_fill[, ..keep_cols],
      future_rows[, ..keep_cols]
    ),
    use.names = TRUE,
    fill = TRUE
  )
  
  if (!nrow(drought_monthly)) {
    stop("No monthly drought records were created.")
  }
  
  # Primary historical source wins over SSP fill when both exist.
  data.table::setorder(drought_monthly,
                       wsd_source_id,
                       park_unit,
                       model,
                       ssp,
                       run,
                       drought_period,
                       ym,
                       source_priority)
  
  drought_monthly <- drought_monthly[
    ,
    .SD[1L],
    by = .(wsd_source_id, park_unit, model, ssp, run, drought_period, ym)
  ]
  
  drought_monthly[
    ,
    reporting_period := (
      drought_period == "historical" &
        analysis_year >= hist_start_year &
        analysis_year <= hist_end_year
    ) |
      (
        drought_period == "future" &
          analysis_year >= future_start_year &
          analysis_year <= future_end_year
      )
  ]
  
  drought_monthly[
    ,
    input_window := data.table::fifelse(reporting_period, "reporting", "padding")
  ]
  
  preferred_order <- c("wsd_source_id",
                       "park_unit",
                       "ym",
                       "year",
                       "month",
                       "analysis_year",
                       "drought_period",
                       "input_window",
                       "reporting_period",
                       "model",
                       "ssp",
                       "run",
                       "monthly_precipitation",
                       "monthly_pet",
                       "climatic_water_balance",
                       "source_type",
                       "source_priority",
                       "source_model",
                       "source_ssp",
                       "source_run")
  
  data.table::setcolorder(drought_monthly,
                          c(preferred_order[preferred_order %in% names(drought_monthly)],
                            setdiff(names(drought_monthly), preferred_order)))
  
  data.table::setorder(drought_monthly,
                       wsd_source_id,
                       park_unit,
                       model,
                       ssp,
                       run,
                       drought_period,
                       ym)
  
  # Save results
  if (save) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    output_stub <- paste0(
      hist_start_year, "-", hist_end_year,
      "_",
      future_start_year, "-", future_end_year,
      "_hindcast-", hindcast
    )
    
    monthly_file <- file.path(
      output_dir,
      paste0("loca2_drought_monthly_", output_stub, ".parquet")
    )
    
    arrow::write_parquet(drought_monthly, monthly_file)
    
    if (verbose) {
      message("Saved joined monthly drought data to: ", monthly_file)
    }
  }
  
  if (verbose) {
    elapsed <- Sys.time() - start_time
    message("\nCompleted in ", round(elapsed, 2), " ", attr(elapsed, "units"))
  }
  
  drought_monthly
}
