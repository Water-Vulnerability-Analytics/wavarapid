# Write the golden files for test-golden.R from the OLD code in
# nps_water_vulnerability, as it was before the split.
#
# Run from the wavarapid repo root, in a fresh R session (not with wavarapid
# loaded):
#   Rscript tests/testthat/fixtures/make_golden.R ../nps_water_vulnerability
#
# The baseline is pinned to one commit, `golden_commit` below: upstream/main of
# nps_water_vulnerability when wavarapid was split out (2026-10-20). The golden
# files always come from that commit, whatever the playground clone is on, so
# later changes to the playground don't move the baseline. The script extracts
# that commit's setup.R, src/ and rapid_app/calc_vulnerability_index.R into a
# temporary folder (`git archive`; the clone itself is not touched), links the
# playground's `data` folder into it, and sources setup.R there as the old
# workflow did (about 70 packages plus every file in src/).
#
# It then sources rapid_app/calc_vulnerability_index.R. wavarapid's
# calc_vulnerability_index() and minmax_norm() came from that file (the app's
# version), not from the older src/rapid/calc_vulnerability_index.R that
# setup.R sources. Nothing else in the old src/ calls minmax_norm().
#
# Next it makes the inputs in golden_cases.R (the water supplies from the old
# get_wsd(), and so on) and saves them to golden/inputs.rds, which test-golden.R
# reads, so both runs use the identical inputs. Then it runs each case and
# writes one golden/<id>.rds per case (after the case's `keep`, if it has one),
# golden/golden_info.txt saying where they came from, and golden/golden_env.csv
# (R and package versions, checksums of the data files and folders) for
# test-golden.R to compare against when a case fails.
#
# Everything is written to a temporary folder first. golden/ is replaced only
# if every case ran, so a failed run leaves the old golden files as they were.
#
# This baseline is for checking the move out of the playground. Once
# development has moved off the playground, it is no longer the reference.
#
# Only regenerate the golden files on purpose (e.g. the source data changed and
# the team agreed); never to make a failing test pass.

golden_commit <- "b62f366"
index_file <- file.path("rapid_app", "calc_vulnerability_index.R")

args <- commandArgs(trailingOnly = TRUE)
old_repo <- normalizePath(
  if (length(args) > 0) args[1] else "../nps_water_vulnerability",
  mustWork = TRUE
)
fixtures <- normalizePath(file.path("tests", "testthat", "fixtures"),
                          mustWork = TRUE)
golden_dir <- file.path(fixtures, "golden")

if ("wavarapid" %in% loadedNamespaces()) {
  stop("wavarapid is loaded; run this in a fresh R session")
}
if (!dir.exists(file.path(old_repo, "data"))) {
  stop("no data folder or link in ", old_repo)
}

git <- function(...) {
  system2("git", c("-C", shQuote(old_repo), ...), stdout = TRUE)
}
golden_commit_full <- git("rev-parse", "--verify", "--quiet",
                          paste0(golden_commit, "^{commit}"))
if (length(golden_commit_full) == 0) {
  stop("commit ", golden_commit, " not found in ", old_repo,
       "; run git fetch upstream")
}
clone_head <- git("rev-parse", "--short", "HEAD")

# the old code exactly as it was at the pinned commit
old_code <- file.path(tempfile("golden_"), "nps_water_vulnerability")
dir.create(old_code, recursive = TRUE)
archive <- tempfile(fileext = ".tar")
git("archive", "--format=tar", "-o", shQuote(archive), golden_commit,
    "setup.R", "src", index_file)
utils::untar(archive, exdir = old_code)
file.symlink(normalizePath(file.path(old_repo, "data")),
             file.path(old_code, "data"))

# the new golden files go here until every case has run
out_dir <- tempfile("golden_out_")
dir.create(out_dir)

setwd(old_code)
source("setup.R")
source(index_file)

# Paths built with golden_data_path() resolve to the same `data` link
options(wava.data_dir = normalizePath("data"))

# The same run settings as test-golden.R, so the results don't depend on the
# session. Sort order "C", which testthat uses: list.files() sorts by the
# session's collation, which changes the row order of
# process_storm_surge_data(). Single-threaded arrow reads: with threads, arrow
# returns rows in a varying order, and calc_loca2_change_summary()'s sums
# change in the last bits from run to run (issues.md Q40).
Sys.setlocale("LC_COLLATE", "C")
options(arrow.use_threads = FALSE)

# make_fake_supplies(), for the index inputs
source(file.path(fixtures, "..", "helper-fixtures.R"))
source(file.path(fixtures, "golden_cases.R"))

info <- c(
  paste("Written:", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  paste("Old code: nps_water_vulnerability", golden_commit_full,
        "(pinned; setup.R, src/ and", index_file,
        "extracted from that commit)"),
  paste("calc_vulnerability_index() and minmax_norm() from", index_file,
        "(sourced after setup.R)"),
  paste("Playground clone was on:", clone_head, "(not used for the code)"),
  paste("Data folder:", normalizePath("data")),
  paste("Run settings: LC_COLLATE =", Sys.getlocale("LC_COLLATE"),
        "and options(arrow.use_threads = FALSE), as in test-golden.R"),
  paste("R:", R.version.string),
  "Package versions and data file checksums: golden_env.csv",
  ""
)

inputs <- list()
for (name in names(golden_inputs)) {
  message("Making input ", name)
  inputs[[name]] <- golden_inputs[[name]](inputs)
  info <- c(info, sprintf("input %s: %s, %d rows", name,
                          class(inputs[[name]])[1], NROW(inputs[[name]])))
}
saveRDS(inputs, file.path(out_dir, "inputs.rds"))
info <- c(info, "")

for (case in golden_cases) {
  message("Running ", case$id)
  warns <- character()
  value <- withCallingHandlers(
    golden_run(case, inputs),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    error = function(e) {
      message(case$id, " failed, so golden/ was not changed")
    }
  )
  saveRDS(value, file.path(out_dir, paste0(case$id, ".rds")))
  warn_text <- if (length(warns) > 0) {
    paste0(": ", paste(unique(warns), collapse = " | "))
  } else {
    ""
  }
  info <- c(info, sprintf("%s: %s, %d rows, %d warning(s)%s",
                          case$id, class(value)[1], NROW(value),
                          length(warns), warn_text))
}

writeLines(info, file.path(out_dir, "golden_info.txt"))

# Recorded after the cases ran, so packages they loaded are included
packages <- sort(loadedNamespaces())
data_files <- unique(c(golden_input_files, golden_data_files))
env <- rbind(
  data.frame(item = "R", name = "R", value = as.character(getRversion())),
  data.frame(item = "package", name = packages,
             value = vapply(packages, function(p) {
               as.character(utils::packageVersion(p))
             }, character(1))),
  data.frame(item = "data file", name = data_files,
             value = vapply(file.path("data", data_files),
                            golden_checksum, character(1), USE.NAMES = FALSE))
)
utils::write.csv(env, file.path(out_dir, "golden_env.csv"), row.names = FALSE)

# Every case ran: replace the old golden files with the new ones
unlink(list.files(golden_dir, full.names = TRUE))
dir.create(golden_dir, showWarnings = FALSE)
copied <- file.copy(list.files(out_dir, full.names = TRUE), golden_dir)
if (!all(copied)) stop("could not copy the golden files into ", golden_dir)
message("Golden files written to ", golden_dir)
