## EventProjection 0.2.9.2

### Performance and memory

- Reworked prediction summaries in `predictEnrollment()` and `predictEvent()` to avoid large Cartesian expansions and reduce memory use for larger `nreps`.
- Added output controls for large subject-level simulation objects so the package can avoid returning bulky simulated subject/event tables when they are not needed.

### API and package surface

- Added size-control arguments for large simulation outputs, including `return_new_subjects`, `return_new_events`, `return_subject_data`, and `return_simulation_data`.
- Removed non-core simulation and testing utilities from the public package surface: `FH_test()`, `tte()`, `test_procedure()`, `test_plot()`, `b4pst()`, and `smed()`.

### Package maintenance

- Cleaned `NAMESPACE`, help files, and vignette references after removing the non-core utilities.
- Enhanced the packaged Shiny app with prior inputs for supported real-time models, separate observed-data and prior-informed fit views, and optional parameter summaries for event and dropout models.
- Updated package metadata for version `0.2.9.2`.

## EventProjection 0.2.9.1

## Major changes

- Forked from **eventPred** 0.2.9 and renamed to **EventProjection**.
- Merged cure-rate and delayed-treatment-effect functionality adapted from **EventPredInCure** into the main package workflow.
- Added support for cure-aware event projection with either estimated cure fractions or user-specified fixed cure rates.
- Consolidated fitting and prediction under a unified interface built around `fitEnrollment()`, `fitEvent()`, `fitDropout()`, `predictEnrollment()`, `predictEvent()`, and `getPrediction()`.

## Modeling updates

- Preserved standard non-cure, non-delay workflows so results remain as close as possible to **eventPred**.
- Added merged support for cure-aware event models inherited from **EventPredInCure**.
- Added support for `piecewise uniform` enrollment modeling in the merged workflow.
- Standardized internal handling of input column names by converting them to lower case.
- Added support for user-controlled random seeds through `seed.num` in the unified prediction workflow.

## Shiny app

- Updated the packaged Shiny app so it can work with the merged modeling workflow.
- Added support for importing study data from `.csv` and `.rda`/`.RData` files.
- Added `runShinyApp()` as the preferred public launcher for the packaged app.
- Retained `runShinyApp_eventPred()` for backward compatibility.
- Moved app-only dependencies to optional installation via `Suggests`, with runtime package checks before app launch.

## Package maintenance

- Updated package naming and documentation to reflect the **EventProjection** identity.
- Prepared the package for GitHub distribution and ongoing development outside CRAN.
- Maintained by **Haiming Zhou**.
