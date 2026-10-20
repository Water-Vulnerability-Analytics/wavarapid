#' @keywords internal
"_PACKAGE"

## Functions from other packages that the code calls without a `pkg::` prefix.
## The playground got these from `library()` calls in setup.R.
#' @importFrom arrow read_parquet write_parquet
#' @importFrom data.table as.data.table data.table fifelse fread
#' @importFrom data.table setDT tstrsplit
#' @importFrom dplyr .data across all_of any_of arrange as_tibble bind_rows
#' @importFrom dplyr case_when coalesce count cross_join cume_dist desc
#' @importFrom dplyr distinct everything filter first group_by if_all if_else
#' @importFrom dplyr inner_join last_col left_join matches mutate n pull
#' @importFrom dplyr relocate rename rename_with select slice slice_max
#' @importFrom dplyr slice_min starts_with summarise summarize tibble
#' @importFrom dplyr transmute tribble ungroup
#' @importFrom lubridate month year ym
#' @importFrom ncdf4 nc_close nc_open ncatt_get ncvar_get
#' @importFrom purrr keep map map_chr map_dbl map_dfr map_int possibly pwalk
#' @importFrom purrr walk
#' @importFrom readr cols read_csv read_lines write_csv
#' @importFrom sf st_area st_as_sf st_bbox st_buffer st_crs st_distance
#' @importFrom sf st_drop_geometry st_filter st_intersection st_intersects
#' @importFrom sf st_is_empty st_join st_make_valid st_nearest_feature st_read
#' @importFrom sf st_transform st_union st_write
#' @importFrom stats coef lm median na.omit predict quantile sd weighted.mean
#' @importFrom stringr `%>%` str_extract str_remove str_replace str_split
#' @importFrom stringr str_starts
#' @importFrom terra crds crs rast relate vect
#' @importFrom tibble rowid_to_column
#' @importFrom tidyr crossing drop_na expand_grid pivot_longer pivot_wider
#' @importFrom utils data download.file head read.table write.csv
#' @importFrom wavacore getUnitVisitation wava_data_path
NULL

## Column names used inside dplyr and data.table calls.
utils::globalVariables(c(
  ":=", ".I", ".N", ".SD",
   ".", "..keep_cols", "1951-1980_best-fit",
   "1951-1980_high-bound-confidence-interval",
   "1951-1980_low-bound-confidence-interval", "1985-2014_best-fit",
   "1985-2014_high-bound-confidence-interval",
   "1985-2014_low-bound-confidence-interval", "2021-2050_best-fit",
   "2021-2050_high-bound-confidence-interval",
   "2021-2050_low-bound-confidence-interval", "abs_timing_change",
   "absolute_change", "agg", "alt_timing_p90", "analysis_year",
   "annual_total", "annual_value", "aq30_id", "baseline_water_stress",
   "baseline_water_stress_label", "baseline_wta", "buffer_km", "bws_label",
   "bws_raw", "bws_score", "Case_C", "cdd_change", "cdd_future",
   "cdd_historical", "cdd_historical_mean", "cdd_percent_change",
   "cdd_percent_slope", "cdd_slope", "cell_lat", "cell_lon", "change_units",
   "change_value", "climate_model", "climatic_water_balance", "coast_dist_m",
   "coast_elev_m", "combo_id", "county_area_km2", "cwd_change", "cwd_future",
   "cwd_historical", "cwd_historical_mean", "cwd_percent_change",
   "cwd_percent_slope", "cwd_slope", "data_status", "days_in_month",
   "delta_fp", "delta_fp_p10", "delta_fp_p50", "delta_fp_p90", "depth_ft",
   "depth_m", "depth_mm", "diff_from_alt_timing_p90", "diff_from_p10",
   "diff_from_p50", "diff_from_p90", "dist_weight", "distance_m",
   "drought_period", "dry_spell_max_length", "elev_weight", "elevation",
   "emission_scenario", "extraction_id", "factor_exp_drought",
   "factor_exp_slr", "factor_exp_wildfire", "factor_sen_drought",
   "factor_sen_flood", "factor_sen_slr", "factor_sen_wildfire", "fips",
   "FIPS", "fire_models_showing_increase", "flag_drought", "flag_fire",
   "flag_flood", "flag_slr", "FLD_ZONE", "flooded", "future_mean",
   "future_months_used", "future_period", "future_record", "future_value",
   "future_wta", "future_years_used", "GCM", "GEOID", "geom", "geometry",
   "has_data", "hist_months_used", "hist_years_used",
   "historic_withdrawal_mgd", "historical_period", "historical_record",
   "historical_value", "huc", "i.cell_lat", "i.cell_lon", "i.model",
   "i.park_unit", "i.point_id", "i.run", "i.ssp", "i.wsd_source_id", "ID",
   "in_pub_file", "info", "input_window", "inun_centroid_change",
   "inun_change", "inun_proj_closest_mm", "inun_proj_mm", "inun_region",
   "inun_state", "inun_surf_url", "irr", "is_gw_spring", "is_training",
   "is_worst", "join_id", "joint_weight", "label", "lat", "lon", "max_swe",
   "mean_swe", "mean_temperature", "missing_plant_records", "mk",
   "mk_p_value", "mk_tau", "model", "Month", "month_name", "monthly_balance",
   "monthly_percentage", "monthly_pet", "monthly_precip",
   "monthly_precipitation", "monthly_runoff", "n_complete_months",
   "n_earlier", "n_expected", "n_later", "n_months", "n_no_change",
   "n_plants", "n_projections", "n_public_months", "n_valid",
   "n_valid_models", "n_years", "NAME", "name_1", "nc_file",
   "norm_cell_change", "norm_coverage_pct_proj", "norm_coverage_pct_ref",
   "norm_pct_point_change", "p10", "p10_timing_change", "p50", "p90",
   "p90_timing_change", "park_name", "park_unit", "percent_change",
   "percentile", "period", "period_variant", "pfaf_id", "pfaf_id_baseline",
   "pfaf_id_future", "plant_id", "point_id", "population",
   "precip_alt_timing_change_direction", "precipitation_future",
   "precipitation_historical", "proj_file", "pub_irr_mgd", "pub_per_cap",
   "PYROME", "rank_cell_change", "rank_coverage_pct_proj",
   "rank_coverage_pct_ref", "rank_pct_point_change", "ratio",
   "raw_avg_nearby_use", "raw_avg_precip", "raw_avg_runoff", "raw_avg_swe",
   "raw_cell_change", "raw_coverage_pct_proj", "raw_coverage_pct_ref",
   "raw_file", "raw_flood_risk_area_m2", "raw_flood_risk_percent",
   "raw_flooded_cells_proj", "raw_flooded_cells_ref", "raw_pct_point_change",
   "raw_precip_alt_timing_change_magnitude",
   "raw_runoff_alt_timing_change_magnitude", "raw_slope_precip",
   "raw_slope_runoff", "raw_slope_swe", "raw_surge_change_max",
   "raw_surge_pct_change_max", "raw_surge_reference_m", "raw_swi_migration",
   "raw_total_cells_proj", "raw_total_cells_ref", "raw_water_stress",
   "raw_wildfire_hazard_mean", "rcp", "RecreationVisitors", "ref_file",
   "reference_file", "reference_mean", "region", "reporting_period",
   "rsl_total_IntHigh_percentiles=83_4", "run",
   "runoff_alt_timing_change_direction", "runoff_future", "runoff_historical",
   "score", "season", "selected_timing_change", "sen", "site_precip",
   "site_precipitation", "site_runoff", "site_temperature", "site_value",
   "slope_adjustment", "source_latitude", "source_longitude", "source_model",
   "source_priority", "source_run", "source_ssp", "source_state",
   "source_type", "spei_change", "spei_future",
   "spei_future_report_inf_months", "spei_future_report_neg_inf_months",
   "spei_future_report_pos_inf_months", "spei_hist_report_inf_months",
   "spei_hist_report_neg_inf_months", "spei_hist_report_pos_inf_months",
   "spei_historical", "spei_inf_months", "spei_neg_inf_months",
   "spei_pos_inf_months", "spei_raw", "spei_report_inf_months",
   "spei_report_months", "spei_report_neg_inf_months",
   "spei_report_pos_inf_months", "spei_total_fit_months", "spei_ym", "ssp",
   "st_agr<-", "state", "state_name", "STATEFP", "station_idx",
   "stress_label", "stress_raw", "stress_score", "STUSPS", "swe_future",
   "swe_historical", "system_type", "temperature_future",
   "temperature_historical", "tep_mgd", "tied_rcps", "time", "time_index",
   "timing_angle", "timing_angle_future", "timing_angle_historical",
   "timing_change", "timing_x", "timing_y", "tmp_join_id", "total",
   "total_mgd", "total_precip", "total_precipitation", "total_runoff",
   "total_visitation", "value", "water_stress_change",
   "water_stress_change_direction", "water_stress_label",
   "water_stress_label_changed", "water_stress_label_transition",
   "water_stress_rank_change", "water_system_name", "watershed_id",
   "WatershedID", "well_depth", "well_depth_m", "wet_spell_max_length",
   "withdrawal_mgd", "worst_rcp", "wsd_source_id", "wta_change",
   "wtd_point_m", "Year", "year_end", "year_single", "year_start"
))
