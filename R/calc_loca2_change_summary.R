#' Summarize LOCA2 climate changes across the projection ensemble
#'
#' Calculates site-level changes for AET, PET, runoff, snow, storage,
#' precipitation, and temperature for one or more future periods relative to a
#' common historical period. Changes are calculated separately for each
#' GCMxSSP projection and then summarized using the ensemble minimum, 10th
#' percentile, median, 90th percentile, and maximum.
#'
#' Monthly grid-cell values follow the shared MWBM missing-data rule: every
#' finite cell in the extracted neighborhood is used and the site-month is
#' `NA` only when all cells are missing. `process_loca2_mwbm(...,
#' grab_nearest = TRUE)` should be used upstream when the nearest valid-cell
#' fallback is desired.
#'
#' AET, PET, runoff, and precipitation are summed to annual totals. Snow is
#' represented by annual peak SWE, while storage and temperature use annual
#' means. Ensemble changes are percent changes for all variables except
#' temperature, which uses absolute change in degrees C.
#'
#' @param water_supplies An `sf` object containing `wsd_source_id` and
#'   `park_unit`.
#' @param data_dirs Named character vector or list giving the processed parquet
#'   directory for each variable. Supported names are `aet`, `pet`, `runoff`,
#'   `snow`, `stor`, `precip`, and `temp`.
#' @param variables Character vector of variables to calculate. Defaults to all
#'   supported variables.
#' @param hindcast Logical. If `TRUE`, each GCM's historical run is paired with
#'   its future SSP run. If `FALSE`, the common LOCA2 training record is used as
#'   the historical baseline. Defaults to `TRUE`.
#' @param historical Two-element vector giving the historical start and end
#'   years. Defaults to `c(1980, 2010)`.
#' @param futures Named list of two-element future periods. Defaults to
#'   `list(mid_century = c(2035, 2065), late_century = c(2070, 2100))`.
#' @param water_year Logical. If `TRUE`, annual values use Oct-Sep water years;
#'   otherwise calendar years are used. Defaults to `TRUE`.
#' @param save Logical. Save the summary as a CSV file?
#' @param output_dir Output directory used when `save = TRUE`.
#' @param verbose Logical. Print progress messages?
#'
#' @return A tibble with one row per water supply, variable, and future period.
#'   `raw_historical_value` and `raw_future_value` are the ensemble means of
#'   the paired historical and future period climatologies. `change_min`,
#'   `change_p10`, `change_p50`, `change_p90`, and `change_max` summarize the
#'   GCM/SSP changes. `change_units` identifies whether the changes are
#'   percent or degrees C.
#'
#' @export
calc_loca2_change_summary <- function(water_supplies,
                                      data_dirs = c(aet = wavacore::wava_data_path("rapid", "inputs", "loca2_aet_water_supplies/"),
                                                    pet = wavacore::wava_data_path("rapid", "inputs", "loca2_pet_water_supplies/"),
                                                    runoff = wavacore::wava_data_path("rapid", "inputs", "loca2_runoff_water_supplies/"),
                                                    snow = wavacore::wava_data_path("rapid", "inputs", "loca2_swe_water_supplies/"),
                                                    stor = wavacore::wava_data_path("rapid", "inputs", "loca2_stor_water_supplies/"),
                                                    precip = wavacore::wava_data_path("rapid", "inputs", "loca2_precipitation_water_supplies/"),
                                                    temp = wavacore::wava_data_path("rapid", "inputs", "loca2_temperature_water_supplies/")),
                                      variables = names(data_dirs),
                                      hindcast = TRUE,
                                      historical = c(1980, 2010),
                                      futures = list(mid_century = c(2035, 2065),
                                                     late_century = c(2070, 2100)),
                                      water_year = TRUE,
                                      save = FALSE,
                                      output_dir = NULL,
                                      verbose = TRUE) {
  
  supported_variables <- c("aet", "pet", "runoff", "snow", "stor", "precip", "temp")
  
  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "water_supplies must include wsd_source_id" =
      "wsd_source_id" %in% names(water_supplies),
    "water_supplies must include park_unit" =
      "park_unit" %in% names(water_supplies),
    "historical must have two years" = length(historical) == 2L,
    "futures must be a non-empty list" = is.list(futures) && length(futures) > 0L,
    "save is TRUE but output_dir was not provided" =
      !save || !is.null(output_dir)
  )
  
  if (is.list(data_dirs)) {
    # Accept either a named list or named character vector for convenience.
    data_dirs <- unlist(data_dirs, use.names = TRUE)
  }
  
  if (is.null(names(data_dirs)) || any(names(data_dirs) == "")) {
    stop("data_dirs must be named with supported variable names.")
  }
  
  variables <- unique(tolower(variables))
  unknown_variables <- setdiff(variables, supported_variables)
  
  if (length(unknown_variables)) {
    stop("Unsupported variable(s): ", paste(unknown_variables, collapse = ", "),
         ". Supported variables are: ",
         paste(supported_variables, collapse = ", "), ".")
  }
  
  missing_dirs <- setdiff(variables, names(data_dirs))
  
  if (length(missing_dirs)) {
    stop("data_dirs is missing path(s) for: ",
         paste(missing_dirs, collapse = ", "), ".")
  }
  
  historical <- as.integer(historical)
  
  if (anyNA(historical) || historical[1] > historical[2]) {
    stop("historical must contain valid start and end years in ascending order.")
  }
  
  if (is.null(names(futures))) {
    # Period names become output labels, so create stable names when the caller
    # supplies an unnamed list.
    names(futures) <- paste0("future_", seq_along(futures))
  }
  
  if (any(names(futures) == "")) {
    blank_names <- which(names(futures) == "")
    names(futures)[blank_names] <- paste0("future_", blank_names)
  }
  
  futures <- lapply(futures, as.integer)
  invalid_futures <- vapply(futures, function(x) length(x) != 2L || anyNA(x) || x[1] > x[2],
                            logical(1))
  
  if (any(invalid_futures)) {
    stop(
      "Each future period must contain valid start and end years. Invalid: ",
      paste(names(futures)[invalid_futures], collapse = ", "), "."
    )
  }
  
  # Define the complete date window that must be read
  start_time <- Sys.time()
  site_ids <- unique(water_supplies$wsd_source_id)
  earliest_year <- min(c(historical[1], vapply(futures, `[`, integer(1), 1L)))
  latest_year <- max(c(historical[2], vapply(futures, `[`, integer(1), 2L)))
  
  if (water_year) {
    # Water year Y begins in October of Y - 1, so include the preceding fall.
    input_start <- as.Date(paste0(earliest_year - 1L, "-10-01"))
    input_end <- as.Date(paste0(latest_year, "-09-30"))
  } else {
    input_start <- as.Date(paste0(earliest_year, "-01-01"))
    input_end <- as.Date(paste0(latest_year, "-12-31"))
  }
  
  # Parse model metadata from processed-file names
  # The processed inputs use two filename families:
  #   var_MWBM_LOCA2.MODEL.SSP.RUN_YEARS.nc
  #   var.MODEL.SSP.RUN.YEARS....nc
  parse_file_metadata <- function(files) {
    parsed <- lapply(files, function(file) {
      stem <- basename(file)
      stem <- sub("\\.parquet$", "", stem, ignore.case = TRUE)
      stem <- sub("\\.nc$", "", stem, ignore.case = TRUE)
      is_training <- grepl("training|livneh|lusu", stem, ignore.case = TRUE)
      
      if (is_training) {
        # Training/observational products have no GCM/SSP/run identifiers. Give
        # them common keys so they can serve as one shared historical baseline.
        return(data.table::data.table(file = file,
                                      model = "training",
                                      ssp = "training",
                                      run = "training",
                                      is_training = TRUE))
      }
      
      normalized <- sub("^.*_MWBM_LOCA2\\.", "", stem)
      
      if (identical(normalized, stem)) {
        # The direct LOCA2 precipitation and temperature names begin with the
        # variable followed immediately by model, SSP, and run fields.
        normalized <- sub("^[^.]+\\.", "", stem)
      }
      
      parts <- strsplit(normalized, ".", fixed = TRUE)[[1]]
      
      if (length(parts) < 3L) {
        # Ignore unrelated parquet files whose names cannot supply the required
        # model, SSP, and run keys.
        return(NULL)
      }
      
      data.table::data.table(file = file, model = parts[1], ssp = parts[2],
                             run = sub("_.*$", "", parts[3]), is_training = FALSE)
    })
    
    data.table::rbindlist(parsed, use.names = TRUE, fill = TRUE)
  }
  
  # Read and summarize each requested climate variable:
  # Each list element will hold model-level changes for one variable across all
  # requested future periods.
  variable_changes <- vector("list", length(variables))
  
  for (variable_index in seq_along(variables)) {
    variable <- variables[variable_index]
    data_dir <- unname(data_dirs[variable])
    
    if (!dir.exists(data_dir)) {
      stop("Processed ", variable, " directory does not exist: ", data_dir)
    }
    
    files <- list.files(data_dir, full.names = TRUE,
                        pattern = "\\.parquet$", ignore.case = TRUE)
    
    if (!length(files)) {
      stop("No processed parquet files found for ", variable, " in ", data_dir)
    }
    
    file_metadata <- parse_file_metadata(files)
    
    if (!nrow(file_metadata)) {
      stop("No model metadata could be parsed for ", variable, ".")
    }
    
    if (hindcast) {
      # Retain modeled historical runs and modeled future scenarios. Training
      # products are intentionally excluded from a hindcast baseline.
      file_metadata <- file_metadata[
        !is_training & (ssp == "historical" | grepl("^ssp[0-9]+$", ssp))
      ]
    } else {
      # Prefer the consolidated Jay Alder training file when it exists so it is
      # not duplicated by older Livneh/LUSU files in the same directory.
      consolidated_training <- file_metadata[
        is_training & grepl("training", basename(file), ignore.case = TRUE)
      ]
      
      if (nrow(consolidated_training)) {
        file_metadata <- file_metadata[
          !is_training | grepl("training", basename(file), ignore.case = TRUE)
        ]
      }
      
      file_metadata <- file_metadata[
        is_training | grepl("^ssp[0-9]+$", ssp)
      ]
    }
    
    if (!nrow(file_metadata)) {
      stop("No ", variable, " files remain for hindcast = ", hindcast, ".")
    }
    
    if (verbose) {
      message(
        "Processing ", variable, " (", variable_index, "/",
        length(variables), "): ", nrow(file_metadata), " file(s)"
      )
    }
    
    # Store annual rather than monthly results from each file to limit memory
    # use during the full multi-variable, multi-model calculation.
    annual_files <- vector("list", nrow(file_metadata))
    
    for (file_index in seq_len(nrow(file_metadata))) {
      file <- file_metadata$file[file_index]
      
      dataset <- arrow::open_dataset(file, format = "parquet")
      required_columns <- c("wsd_source_id", "park_unit", "date", "value")
      absent_columns <- setdiff(required_columns, names(dataset))
      
      if (length(absent_columns)) {
        stop(
          basename(file), " is missing required column(s): ",
          paste(absent_columns, collapse = ", "), "."
        )
      }
      
      # Push site and date filters into Arrow before collecting into memory.
      raw <- dataset %>%
        dplyr::select(dplyr::all_of(required_columns)) %>%
        dplyr::filter(wsd_source_id %in% site_ids,
                      date >= as.character(input_start),
                      date <= as.character(input_end)) %>%
        dplyr::collect() %>%
        data.table::as.data.table()
      
      if (!nrow(raw)) {
        next
      }
      
      raw[, date := data.table::as.IDate(date)]
      raw[, value := as.numeric(value)]
      
      # Assign each monthly record to either an Oct-Sep water year or calendar
      # year, depending on the requested convention.
      raw[
        , analysis_year := if (water_year) {
          data.table::year(date) +
            data.table::fifelse(data.table::month(date) >= 10L, 1L, 0L)
        } else {
          data.table::year(date)
        }
      ]
      
      # Average all finite cells in the extracted neighborhood. Return NA only
      # when every neighborhood cell is non-finite.
      monthly <- raw[
        , .(site_value = if (any(is.finite(value))) {
          mean(value[is.finite(value)])
        } else {
          NA_real_
        }),
        by = .(wsd_source_id, park_unit, date, analysis_year)
      ]
      
      # File-level metadata are constant scalars, so attach them after spatial
      # aggregation instead of parsing them repeatedly for every raw row.
      monthly[
        , `:=`(
          model = file_metadata$model[file_index],
          ssp = file_metadata$ssp[file_index],
          run = file_metadata$run[file_index]
        )
      ]
      
      # Variable-specific temporal aggregation follows the existing RAPID
      # indicators: peak snow, mean storage/temperature, and summed fluxes.
      annual <- if (variable == "snow") {
        monthly[
          , .(annual_value = if (any(is.finite(site_value))) {
            max(site_value[is.finite(site_value)])
          } else {
            NA_real_
          }),
          by = .(wsd_source_id, park_unit, analysis_year, model, ssp, run)
        ]
      } else if (variable %in% c("stor", "temp")) {
        monthly[
          , .(annual_value = if (any(is.finite(site_value))) {
            mean(site_value[is.finite(site_value)])
          } else {
            NA_real_
          }),
          by = .(wsd_source_id, park_unit, analysis_year, model, ssp, run)
        ]
      } else {
        monthly[
          , .(annual_value = if (any(is.finite(site_value))) {
            sum(site_value[is.finite(site_value)])
          } else {
            NA_real_
          }),
          by = .(wsd_source_id, park_unit, analysis_year, model, ssp, run)
        ]
      }
      
      annual_files[[file_index]] <- annual
    }
    
    annual_data <- data.table::rbindlist(annual_files,
                                         use.names = TRUE,
                                         fill = TRUE)
    
    if (!nrow(annual_data)) {
      stop("No ", variable, " records remained after site/date filtering.")
    }
    
    # Direct temperature files contain separate tasmin and tasmax components;
    # averaging their annual means produces mean annual temperature. This also
    # deduplicates any overlapping, year-aligned input files.
    annual_data <- annual_data[
      , .(annual_value = if (any(is.finite(annual_value))) {
        mean(annual_value[is.finite(annual_value)])
      } else {
        NA_real_
      }),
      by = .(wsd_source_id, park_unit, analysis_year, model, ssp, run)
    ]
    
    # Each unique GCM x SSP x run is one ensemble projection.
    future_combinations <- unique(annual_data[grepl("^ssp[0-9]+$", ssp),
                                              .(model, ssp, run)])
    
    if (!nrow(future_combinations)) {
      stop("No future GCM/SSP/run combinations were found for ", variable, ".")
    }
    
    # Construct the historical baseline for every future projection
    if (hindcast) {
      # Calculate a separate historical climatology for each GCM/run. Joining
      # to future_combinations replicates that baseline across the model's SSPs.
      historical_climatology <- annual_data[
        ssp == "historical" &
          analysis_year >= historical[1] &
          analysis_year <= historical[2],
        .(historical_value = if (any(is.finite(annual_value))) {
          mean(annual_value[is.finite(annual_value)])
        } else {
          NA_real_
        }),
        by = .(wsd_source_id, park_unit, model, run)
      ]
      
      historical_climatology <- merge(historical_climatology,
                                      future_combinations,
                                      by = c("model", "run"),
                                      allow.cartesian = TRUE)
    } else {
      # Calculate one training-data climatology per site, then cross it with all
      # future projections so every projection uses the same observed baseline.
      historical_climatology <- annual_data[
        ssp == "training" &
          analysis_year >= historical[1] &
          analysis_year <= historical[2],
        .(historical_value = if (any(is.finite(annual_value))) {
          mean(annual_value[is.finite(annual_value)])
        } else {
          NA_real_
        }),
        by = .(wsd_source_id, park_unit)
      ]
      
      # A constant key provides an explicit Cartesian join.
      historical_climatology[, join_id := 1L]
      future_combinations[, join_id := 1L]
      historical_climatology <- merge(historical_climatology,
                                      future_combinations,
                                      by = "join_id",
                                      allow.cartesian = TRUE)
      historical_climatology[, join_id := NULL]
      future_combinations[, join_id := NULL]
    }
    
    if (!nrow(historical_climatology)) {
      stop(
        "No historical baseline could be paired with future projections for ",
        variable, "."
      )
    }
    
    #  Calculate projection-level changes for each future period
    period_changes <- vector("list", length(futures))
    
    for (period_index in seq_along(futures)) {
      future_years <- futures[[period_index]]
      
      # Reduce annual values to one future-period climatology for every
      # site/GCM/SSP/run combination.
      future_climatology <- annual_data[
        grepl("^ssp[0-9]+$", ssp) &
          analysis_year >= future_years[1] &
          analysis_year <= future_years[2],
        .(future_value = if (any(is.finite(annual_value))) {
          mean(annual_value[is.finite(annual_value)])
        } else {
          NA_real_
        }),
        by = .(wsd_source_id, park_unit, model, ssp, run)
      ]
      
      changes <- merge(
        historical_climatology,
        future_climatology,
        by = c("wsd_source_id", "park_unit", "model", "ssp", "run")
      )
      
      changes[, absolute_change := future_value - historical_value]
      changes[
        , change_value := if (variable == "temp") {
          # Temperature change is meaningful as an absolute difference in C.
          absolute_change
        } else {
          # Hydrologic variables follow existing exposure indicators and use
          # percent change. A zero baseline cannot produce a finite percentage.
          data.table::fifelse(
            is.finite(historical_value) & historical_value != 0,
            100 * absolute_change / historical_value,
            NA_real_
          )
        }
      ]
      changes[!is.finite(change_value), change_value := NA_real_]
      changes[
        , `:=`(
          variable = variable,
          period = names(futures)[period_index],
          historical_period = paste(historical, collapse = "-"),
          future_period = paste(future_years, collapse = "-"),
          change_units = if (variable == "temp") "degrees_C" else "percent"
        )
      ]
      
      period_changes[[period_index]] <- changes
    }
    
    variable_changes[[variable_index]] <- data.table::rbindlist(period_changes,
                                                                use.names = TRUE,
                                                                fill = TRUE)
  }
  
  # Combine variables only after every model-level calculation is complete.
  model_changes <- data.table::rbindlist(variable_changes,
                                         use.names = TRUE,
                                         fill = TRUE)
  
  # Summarize change across the GCM x SSP x run ensemble
  summary_results <- model_changes[
    , {
      # Report the mean raw period climatologies across the same paired
      # GCM/SSP/run combinations used to calculate change.
      valid_pair <- is.finite(historical_value) & is.finite(future_value)
      raw_historical_value <- if (any(valid_pair)) {
        mean(historical_value[valid_pair])
      } else {
        NA_real_
      }
      raw_future_value <- if (any(valid_pair)) {
        mean(future_value[valid_pair])
      } else {
        NA_real_
      }

      # Non-finite projections do not contribute to quantiles or model counts.
      valid_change <- change_value[is.finite(change_value)]
      
      if (!length(valid_change)) {
        list(
          raw_historical_value = raw_historical_value,
          raw_future_value = raw_future_value,
          change_min = NA_real_,
          change_p10 = NA_real_,
          change_p50 = NA_real_,
          change_p90 = NA_real_,
          change_max = NA_real_,
          n_projections = 0L
        )
      } else {
        # These are direct ensemble quantiles, not the nearest individual model
        # values. The endpoints are the ensemble minimum and maximum.
        quantiles <- as.numeric(
          stats::quantile(
            valid_change,
            probs = c(0, 0.10, 0.50, 0.90, 1),
            na.rm = TRUE,
            names = FALSE
          )
        )
        
        list(
          raw_historical_value = raw_historical_value,
          raw_future_value = raw_future_value,
          change_min = quantiles[1],
          change_p10 = quantiles[2],
          change_p50 = quantiles[3],
          change_p90 = quantiles[4],
          change_max = quantiles[5],
          n_projections = length(valid_change)
        )
      }
    },
    by = .(
      wsd_source_id,
      park_unit,
      variable,
      period,
      historical_period,
      future_period,
      change_units
    )
  ]
  
  # The model-level inner joins naturally omit sites without paired data. Build
  # the complete site x variable x period grid and join the summaries onto it so
  # those combinations remain explicit NA rows with n_projections = 0.
  sites <- unique(
    data.table::as.data.table(sf::st_drop_geometry(water_supplies))[
      , .(wsd_source_id, park_unit)
    ]
  )
  requested <- data.table::CJ(
    variable = variables,
    period = names(futures),
    unique = TRUE
  )
  requested[
    , `:=`(
      historical_period = paste(historical, collapse = "-"),
      future_period = vapply(
        period,
        function(x) paste(futures[[x]], collapse = "-"),
        character(1)
      ),
      change_units = data.table::fifelse(
        variable == "temp",
        "degrees_C",
        "percent"
      ),
      join_id = 1L
    )
  ]
  sites[, join_id := 1L]
  requested <- merge(sites, requested, by = "join_id", allow.cartesian = TRUE)
  requested[, join_id := NULL]
  
  join_columns <- c("wsd_source_id",
                    "park_unit",
                    "variable",
                    "period",
                    "historical_period",
                    "future_period",
                    "change_units")
  
  final_results <- merge(requested,
                         summary_results,
                         by = join_columns,
                         all.x = TRUE,
                         sort = FALSE)
  final_results[is.na(n_projections), n_projections := 0L]
  data.table::setcolorder(final_results, c("wsd_source_id",
                                           "park_unit",
                                           "variable",
                                           "period",
                                           "historical_period",
                                           "future_period",
                                           "change_units",
                                           "raw_historical_value",
                                           "raw_future_value",
                                           "change_min",
                                           "change_p10",
                                           "change_p50",
                                           "change_p90",
                                           "change_max",
                                           "n_projections"))
  data.table::setorder(final_results,
                       wsd_source_id,
                       park_unit,
                       variable,
                       period)
  
  final_results <- dplyr::as_tibble(final_results)
  
  if (save) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    output_file <- file.path(
      output_dir,
      paste0("loca2_climate_change_summary_", Sys.Date(), ".csv")
    )
    readr::write_csv(final_results, output_file)
    
    if (verbose) {
      message("Saved climate change summary to: ", output_file)
    }
  }
  
  if (verbose) {
    elapsed <- Sys.time() - start_time
    message(
      "LOCA2 change summary completed in ",
      round(elapsed, 2), " ", attr(elapsed, "units"), "."
    )
  }
  
  final_results
}
