# River-size association of Palaeolithic findspots in Nagpur and Chandrapur districts, Central India

[![DOI](https://doi.org/10.5281/zenodo.23273489)

Code, tests, saved summary tables and figure-source tables for a spatial analysis of 197 Palaeolithic findspots
(59 source references) in the Wainganga-Wardha basin, Central India.

**Authors:** Sushant Begade, Prabash Sahu. [Department of Ancient Indian History Culture and Archaeology, RTMNU]

## Research question
Do Palaeolithic findspots lie closer to larger rivers than expected from survey effort and terrain? Distance is measured to the
nearest drainage line of WRIS order >= k, for k = 3 to 9 (a river-size threshold profile).

## Key design choices
- **Confirmatory test:** k = 6 on open-air sites (S1). Expected direction: intensity falls with distance. Decision rule: supported only if the
  40 km full-refit cluster-bootstrap 95% CI excludes 0 with beta < 0.
- **Model:** down-weighted Poisson point-process likelihood (presence weight 1e-6, quadrature weight 0.25 km2, 500 m grid), ridge prior
  (SD 2.5) on standardised slopes. Survey domain 16,116 km2 (64,464 cells) after removing ponds, reservoirs and quarries (+100 m).
- **Covariates:** elevation, log slope, plateau/hills, WorldCover crop fraction, built fraction, log distance to major road, log distance to settlement.
- **Effect scale:** log intensity ratio per IQR of log distance (beta).
- **Inference:** cluster bootstrap (40 km blocks primary, 1,000 full refits; 20 and 10 km sensitivity), spatial cross-validation,
  lineament-distance negative control, leave-one-reference-out, survey-bias and robustness checks. Exploratory analyses are labelled EXPLORATORY.
- **Sample sets:** S1 open-air (93), S2 all (197), S3 riverbed (88), S4 other (16).

## Headline outcome (see the paper for interpretation)
The pre-specified direction was not supported at k = 6 in S1: the association is positive (intensity rises with distance) and
depends on how sites are classified (exploratory analysis S13). Pooled sites (S2) show no robust association; riverbed sites (S3) are
at water by construction. Neither S1 nor S2 is a classification-independent test of river-size preference.
<!-- OPTIONAL, add numbers only when you are ready to disclose them before submission:
S1, k = 6: beta = 1.236 per IQR (95% CI 0.694 to 1.660; 40 km blocks, 1,000 full refits); intensity ratio near (q25) vs far (q75) 0.29 (0.19 to 0.50). -->

## Repository layout
    R/                      functions used by the scripts (modelling, covariates, diagnostics, figures, reproducibility checks)
    scripts/                numbered drivers S01-S19, run_all.R, freeze_manifest.R, check_repo_safe.R
    tests/testthat/         unit tests (124 pass, 2 skipped on the release machine)
   

## Pipeline
| Script | Purpose (abbreviated) |
|---|---|
| S01-S03 | clean site table, hydrography preparation, persistence index (wet/dry persistence was found unidentifiable and is not the primary analysis) |
| S04-S05b | external data download, regulation flags and cumulative regulation (dams) |
| S06-S08b | quadrature grid, terrain and geomorphology, detectability covariates (roads, places, land cover), OSM diagnostic |
| S09a-e | Step 0 simulations: power, ridge prior, block-size rule, profile recovery |
| S10 | site covariates (computed blind to outcomes) |
| S11 | primary model and confirmatory test |
| S12 | secondary diagnostics (sample sets, coordinate tiers, negative control, leave-one-out) |
| S13 | site-type composition (EXPLORATORY) |
| S14 | spatial cross-validation |
| S15-S16 | survey-bias and robustness checks |
| S17 | exploratory analyses (distance classes, periods) |
| S18 | figures and paper tables, every plotted value read from a saved table and audited |
| S19 | reproducibility checks (tests, deterministic reruns, freeze verification, coordinate scan) |

## Data (not included)
Exact site coordinates and all raw third-party layers are **not** redistributed. Obtain them from the original providers and cite them:
| Layer | Source | Notes |
|---|---|---|
| Drainage network and stream order | India WRIS | [FILL: layer name, version/date, terms] |
| Network extension outside the districts | HydroRIVERS | [FILL: version, licence] |
| Land cover | ESA WorldCover 2021, v200 | [FILL: licence statement from provider] |
| Roads and place nodes | OpenStreetMap via Overpass, snapshot 3 Oct 2026 | (c) OpenStreetMap contributors, ODbL; parsed with R/osm_xml.R |
| Dams | Global Dam Watch (GDW) | [FILL: version, citation] |
| Elevation/terrain | [FILL: DEM name, version, provider] | |
| Geomorphology, lineaments | [FILL: source] | |
[VARIANT A:] A 1 km coarsened site table is not part of this release. Contact the authors for access, subject to heritage-protection rules.
[VARIANT B instead:] `data/public/sites_public_1km.csv` holds site positions snapped to a 1 km grid; exact coordinates are available on request.

## Reproducing the results
Requires R 4.5.x and the raw inputs.

    1. Open the project (RStudio) and run renv::restore()
    2. Put your local data folder in .Renviron:  PALEO_DATA=<path>
    3. Quick check (about 1 minute):  source("scripts/run_all.R")
    4. Full rebuild in fresh sessions (hours; S11 alone takes about 26-29 minutes):
         RUN_MODE <- "full"; source("scripts/run_all.R")

Without the restricted data the quick check still runs tests, the freeze verification, the coordinate scan and S18 on the saved tables
(the Fig. 1 point layer needs the 1 km table).
Random seed is set in `R/config.R`. Each script writes only to `output/` and never overwrites raw data.

## Limitations (see paper)
Modern water bodies used as a proxy for Pleistocene water; legacy survey bias; mixed coordinate quality (tiers A-C); undated typology;
estimated dam capacities; no field validation.

## Data safety and ethics
Exact findspot coordinates are never published here. `scripts/check_repo_safe.R` blocks commits containing coordinate-like columns,
coordinate pairs, spatial file types or large files. [ADD YOUR ETHICS/HERITAGE STATEMENT: permissions, ASI/State Archaeology rules, reason for 1 km coarsening.]

## How to cite
Cite the archive (Zenodo version DOI for exact reproduction, concept DOI for the project) and the paper when published:

    Begade, S., Sahu, P. (2026). Code and derived results: river-size association of Palaeolithic findspots in the Wainganga-Wardha
    basin, Nagpur and Chandrapur districts, Central India (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.23273489

## Licence
Code: MIT. Derived tables and figures: CC BY 4.0. Third-party data keep their own licences.


## Contact
[Sushant Begade: sushantbegade@gmail.com]