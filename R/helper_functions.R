# HELPER FUNCTIONS =========================================================

#' @title Add normalized and ranked versions of raw indicator columns
#'
#' @description
#' For every column whose name starts with `prefix`, adds two columns:
#' * `norm_*`: the values rescaled to 0-1 with [scales::rescale()] (all `NA`
#'   if the column is all `NA`).
#' * `rank_*`: the percentile rank, 0-100, from [dplyr::cume_dist()].
#'
#' The new names replace `raw_` with `norm_` and `rank_`, so
#' `raw_fire_prob_change_p50` gets `norm_fire_prob_change_p50` and
#' `rank_fire_prob_change_p50`. Every `calc_*()` indicator function calls this
#' on its output. Missing values stay missing. On a grouped data frame the
#' scaling and ranking are done within each group.
#'
#' @param data Data frame. The indicator table, with one or more `raw_*`
#'   columns. Required.
#' @param prefix String. The start of the names of the columns to transform,
#'   e.g. `"raw_fire_prob_change"` (matches `raw_fire_prob_change_p10`,
#'   `_p50` and `_p90`). The matching columns should start with `raw_`: the
#'   new names are made by replacing `raw_`, so a matching column without it
#'   is overwritten. Required.
#' @param invert Logical. The direction of the scale.
#'   * `TRUE`: the highest raw value gets the lowest score (`norm` is
#'     `1 - norm`, `rank` is `100 - rank`). Use it when a higher raw value
#'     means lower vulnerability, e.g. an increase in runoff.
#'   * `FALSE`: the highest raw value gets `norm = 1` and the highest rank.
#'
#'   Default: `FALSE`.
#'
#' @return `data` with a `norm_*` and a `rank_*` column added after the
#'   existing columns for each matching column.
#' @export
#' @seealso [minmax_norm()], the rescaling used by
#'   [calc_vulnerability_index()].
#' @examples
#' fire <- data.frame(
#'   wsd_source_id = c("A", "B", "C", "D"),
#'   raw_fire_prob_change_p50 = c(5, 10, 20, NA)
#' )
#'
#' # Higher fire probability change = higher score
#' add_standardized_indicators(fire, "raw_fire_prob_change")
#'
#' # Inverted: higher raw value = lower score
#' add_standardized_indicators(fire, "raw_fire_prob_change", invert = TRUE)

add_standardized_indicators <- function(data, prefix, invert = FALSE) {
  # Get columns that match the prefix
  cols_to_transform <- names(data)[stringr::str_starts(names(data), prefix)]

  # Create normalized columns 
  for (col in cols_to_transform) {
    norm_col <- stringr::str_replace(col, "raw_", "norm_")
    rank_col <- stringr::str_replace(col, "raw_", "rank_")

    data <- data %>%
      dplyr::mutate(
        !!norm_col := if (all(is.na(.data[[col]]))) NA_real_ else scales::rescale(.data[[col]], to = c(0, 1)),
        !!rank_col := dplyr::cume_dist(.data[[col]]) * 100
      )

    if (invert) {
      data <- data %>%
        dplyr::mutate(
          !!norm_col := 1 - .data[[norm_col]],
          !!rank_col := 100 - .data[[rank_col]]
        )
    }

  }

  return(data)
}
