# Golden test fixtures

`test-golden.R` checks that wavarapid returns exactly what the old code in
`nps_water_vulnerability` returned before the split.

## The baseline: `nps_water_vulnerability` at `b62f366`

We chose one commit of the playground as the reference: `b62f366`, its
`upstream/main` when wavarapid was split out (2026-10-20). The golden files
always come from that commit, even after the playground changes.
`make_golden.R` extracts `setup.R`, `src/` and
`rapid_app/calc_vulnerability_index.R` from that commit with `git archive`,
whatever the local clone is on. The commit is set in one place,
`golden_commit` at the top of `make_golden.R`.

The vulnerability index has two old versions at that commit. wavarapid's
`calc_vulnerability_index()` and `minmax_norm()` came from the app's
version, `rapid_app/calc_vulnerability_index.R`, the same file as in
`nps_wsva_app`. `src/rapid/calc_vulnerability_index.R` is an older version
(it has an `id_cols` argument). `setup.R` sources it, so `make_golden.R`
sources the app's file afterwards. Nothing else in the old `src/` calls
`minmax_norm()`.

This is a temporary arrangement for checking the move out of the playground.
Once development has moved off the playground, `b62f366` stops being the
reference, and these tests either retire or become wavarapid's own
regression baseline.

The setup is copied from wavadeepdive's golden tests (which added `reads` and
`keep` to wavacore's), with three more pieces, described below:
- saved inputs;
- `args` that can be a function;
- `golden_run()`, which stops a case from writing outside `tempdir()`.

## Files

- `golden_cases.R`: the inputs, the calls to compare, and the helpers both
  scripts use.
  - **Inputs** (`golden_inputs`): the 8 water supplies the cases run on,
    with their buffered and inundation-scenario versions, the pyromes and
    fire probabilities, and the synthetic tables for the index (from
    `make_fake_supplies()` in `helper-fixtures.R`). `make_golden.R` makes
    them with the old code (the old `get_wsd(rapid = TRUE)`) and saves them
    to `golden/inputs.rds`; `test-golden.R` reads them from there. So both
    runs use the identical inputs, and the test doesn't call `get_wsd()`,
    which wavacore's golden tests cover. The sites are listed, with the
    reason for each, in `golden_site_ids`.
  - **Cases** (`golden_cases`): each has a function, its arguments, and
    whether it needs the shared data folder or internet. `args` is a list,
    or a function of the inputs that returns one. Two optional fields:
    - `reads`: for a data case, the files or folders it reads under the data
      folder. The case skips without them, and `make_golden.R` records their
      checksums. The inundation case lists the four tiles it reads rather
      than the whole 3.7 GB `Extent_Rasters` folder.
    - `keep`: a function applied to the result before it is saved or
      compared. Three cases use it: the two `process_loca2_mwbm()` cases,
      which return nothing and write a parquet file, and
      `process_precip_period_data()`, which returns the path of the parquet
      it writes in `tempdir()`. `keep` returns the file's contents instead.
    - `stop_on`: a pattern for warnings that mean a download failed inside
      a function that carries on anyway. `golden_run()` turns such a
      warning into an error. The two `calc_flood_risk()` cases use it:
      when a FEMA request fails (HTTP 500 or 504, or "Connection reset by
      peer", all seen on 2026-10-27), the function gives that supply `NaN`
      with only a GDAL warning (`issues.md` Q39). Without `stop_on`, an
      outage would show up as a mismatch, or be saved into the golden file.
      These cases run early, so an outage stops `make_golden.R` before the
      long cases.
  - **`golden_run()`** runs a case for both scripts. A case that writes files
    must write into `golden_out_dir(<id>)` inside `tempdir()`. Writing
    means `save = TRUE`, or a `process_*` function without `save = FALSE`,
    since their defaults write into the shared data folder. Otherwise it
    stops before calling the function.

  The file checks the fields of every case when it is sourced, so a typo
  stops with an error instead of silently dropping a skip.
- `make_golden.R`: runs those calls with the old code at `b62f366` and saves
  the results to `golden/<id>.rds`, with the inputs in `golden/inputs.rds`.
  It also writes:
  - `golden/golden_info.txt`: the date, the pinned commit, the commit the
    clone was on, the data folder, and row and warning counts;
  - `golden/golden_env.csv`: R and every package version, and a checksum of
    each data file or folder, including the database the inputs come from.

  It writes to a temporary folder first and replaces `golden/` only if every
  case ran, so a failed run (e.g. a service is down) leaves the old golden
  files as they were.

  Both `make_golden.R` and `test-golden.R` use the same two run settings, so
  results don't depend on the session:
  - **Sort order `"C"`** (`LC_COLLATE`), which testthat uses inside every
    test anyway. `list.files()` sorts by the session's collation, so under
    another one `process_storm_surge_data()` returns its rows in a different
    order. (Here `Rscript` runs with `C.UTF-8`, which R sorts with ICU, so
    `HadGEM3-GC31-HM_` comes before `HadGEM3-GC31-HM-SST`; under `"C"` it
    comes after.)
  - **Single-threaded arrow reads** (`options(arrow.use_threads = FALSE)`).
    With threads, arrow returns rows in a varying order, and
    `calc_loca2_change_summary()` changes in the last bits from run to run,
    in the old code and in wavarapid alike (`issues.md` Q40).

  Both were found on 2026-10-27, when the first golden run (made without these
  settings) failed on exactly those two cases.
- `golden/`: the saved results. Not in git for now; whether the outputs can go
  in this public repo is an open question for the team.

To create the golden files, from the wavarapid repo root in a fresh R session
(the clone needs commit `b62f366` and its `data` link; the script doesn't
change the clone):

```sh
Rscript tests/testthat/fixtures/make_golden.R ../nps_water_vulnerability
```

Most cases read the shared data folder. The LOCA2 cases read every model's
file for their variable, several GB each (precipitation 5.8 GB, temperature
9.1 GB, runoff 4.1 GB, PET 3.6 GB, SWE 2.4 GB), and the WTA cases read the
3.4 GB `WBD_National_GDB.gdb`. OneDrive downloads what isn't on the disk yet
the first time, about 27 GB. Use "Free up space" on those folders afterwards
to get the disk space back.

To run the comparison, point wavarapid at the same data folder the
playground's `data` link points to:

```sh
WAVA_DATA_DIR="/path/to/shared/data" Rscript -e 'devtools::test(filter = "golden")'
```

A few cases download live data, and those sources change: FEMA flood zones,
NPS visitation, county and state boundaries (tigris), and USGS point
elevations. Run the comparison soon after making the golden files. Later
failures may come from the source rather than the package (see below).

The comparison is `expect_identical()`. Every skip says which case it is and
what is missing:

- no golden files at all: one skip for the whole file, saying to run
  `make_golden.R` (this is what `R CMD check` shows, because `golden/` is in
  `.Rbuildignore`);
- no `inputs.rds`: one skip for the whole file, saying to run `make_golden.R`;
- one golden file missing: `<case>: golden file missing: <path>`;
- no data: `<case>: no <missing reads> in the data folder "<path>"`, plus
  whether the path came from the option, `WAVA_DATA_DIR`, or the `"data"`
  fallback;
- offline: `<case>: offline; this case downloads <what>` (the `downloads`
  field in `golden_cases.R`).

A service that is down while the machine is online is not skipped: the case
errors with the service's error, so an outage is never mistaken for a pass.
That is an error, not a mismatch. Two services can fail without an error:
- **FEMA:** `calc_flood_risk()` catches failed requests, and `stop_on`
  turns its warning back into an error.
- **NPS STATS:** `calc_visitation_*()` drop a park whose request fails,
  silently (`issues.md` Q19). If the server fails during a test run, the
  case shows a mismatch whose failure message names the download.

After making the golden files, check that the visitation results have values
for all four parks.

If a test fails, don't regenerate the golden file or loosen the comparison:
find out why the output changed first. A failure can come from wavarapid, or
from something else that changed since the golden run: an R or package
upgrade, a change to the shared data, or a live source. The failure message
lists the R, package and data changes it finds against `golden_env.csv`, and
names the download if the case has one. To find the cause:

1. Run `make_golden.R` from a scratch copy of wavarapid (so the real
   `golden/` stays as it is) to rebuild the golden files with today's
   packages, data and downloads.
2. Compare the fresh files with the saved ones (`identical()`).
3. If they are identical, the old code still gives the saved output, so the
   change is in wavarapid. If they differ, the change came from packages,
   data or a download, not from wavarapid.

## Not covered

- `process_groundwater_depth()`: its default raster is 36.6 GB and kept only
  online.
- `process_inundation_data()`: it overwrites intermediate files in the shared
  input folder (`coast_dir_path`), deletes `output_dir/Extent_Rasters`, and
  downloads several GB from NOAA.
- The `hindcast = FALSE` mode of the LOCA2 functions. The cases use
  `hindcast = TRUE`, as `rapid_workflow.Rmd` does. For precipitation,
  `hindcast = FALSE` needs Livneh files that aren't in the shared data
  folder.
- The `save = TRUE` branches, which write files.
