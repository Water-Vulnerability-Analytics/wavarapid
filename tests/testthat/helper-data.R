# Where the data folder path came from, for skip messages that need the
# shared WaVA data folder.
data_dir_source <- function() {
  if (!is.null(getOption("wava.data_dir"))) {
    "set by options(wava.data_dir)"
  } else if (nzchar(Sys.getenv("WAVA_DATA_DIR"))) {
    "set by WAVA_DATA_DIR"
  } else {
    "WAVA_DATA_DIR is not set, so it is \"data\" in the working directory"
  }
}
