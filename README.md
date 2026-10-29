# wavarapid

Indicator functions for the rapid, nationwide Water Vulnerability Analytics
(WaVA) water supply vulnerability assessment of National Park Service (NPS)
units.

## What is in it

| Function family | What it does |
| --- | --- |
| `calc_*_exposure()` | Exposure indicators: change in a hazard or water supply condition at each water supply (runoff, precipitation, snow, drought, fire, storm surge, sea level rise, saltwater intrusion, nearby use, water rights, and others) |
| `calc_*_sensitivity()` | Sensitivity indicators: how sensitive a water supply is to those changes |
| `process_*()` | Prepare the input data that some indicators need (coastal geometry, storm surge, drought, groundwater depth, and others) |
| `add_standardized_indicators()` | Adds normalized and ranked columns to an indicator table |
| `calc_vulnerability_index()` | Combines the raw indicators into factors, components and the total vulnerability index. This is the same function the app uses. |

## Install

`wavarapid` depends on [`wavacore`](https://github.com/Water-Vulnerability-Analytics/wavacore).

```r
# install.packages("remotes")
remotes::install_github("Water-Vulnerability-Analytics/wavarapid")
```

## Data

The indicator functions read their inputs from the shared WaVA data folder,
which is not part of this repository. The default paths are built with
`wavacore::wava_data_path()`, which looks for the folder in this order:
`options(wava.data_dir = ...)`, the `WAVA_DATA_DIR` environment variable, then
a folder called `data` in the working directory. Every input and output path
is also an argument, so you can pass your own.

## AI Use Disclosure

AI-assisted coding tools (Claude Code, using Anthropic's Claude models) were used
during development of this repository, beginning in 2026, to help set up the
package, prepare it for release, write tests, and draft documentation. All
AI-assisted code and text was reviewed, tested, and validated by the authors,
who made all design and scientific decisions and remain responsible for the
content.

## License

MIT. See [LICENSE.md](LICENSE.md).
