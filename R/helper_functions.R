# HELPER FUNCTIONS =========================================================

#' Add standardized normalized and ranked columns to indicator data
#' 
#' @param data tibble with raw indicator columns
#' @param prefix character, prefix pattern to match columns (e.g., "raw_fire_prob")
#' @param invert logical, whether to invert the direction of scaling so that
#'   higher raw values correspond to lower normalized and ranked scores.
#' @return tibble with added norm_ and rank_ columns
#' 
add_standardized_indicators <- function(data, prefix, invert = FALSE) {
  # Get columns that match the prefix
  cols_to_transform <- names(data)[str_starts(names(data), prefix)]
  
  # Create normalized columns
  for (col in cols_to_transform) {
    norm_col <- str_replace(col, "raw_", "norm_")
    rank_col <- str_replace(col, "raw_", "rank_")
    
    data[[norm_col]] <- scales::rescale(data[[col]], to = c(0, 1))
    data[[rank_col]] <- cume_dist(data[[col]]) * 100
    
    if (invert) {
      data[[norm_col]] <- 1 - data[[norm_col]]
      data[[rank_col]] <- 100 - data[[rank_col]]
    }
    
  }
  
  return(data)
}