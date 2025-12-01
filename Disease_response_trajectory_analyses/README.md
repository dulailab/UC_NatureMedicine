# Disease Response Trajectory Analyses

R scripts for trajectory visualization, survival analysis, and clustering examples. 

## Files
- `3TRTrajectoryPlot.R`: `TrajectoryVisualization()` function that builds PRO2 trajectory plots with LOESS smoothing, optional labels, and PDF export. Includes input validation and styling.
- `gemini time to event.R`: Kaplan–Meier survival curves for GEMINI (hospital/procedure outcomes) using `ggsurvfit`/`survminer`; saves PDFs with risk tables.
- `trajectory analysis example code.R`: Example of trajectory clustering with `gbmt` for Adalimumab/Vedolizumab/GEMINI PRO2 data; compares information criteria across cluster counts.

## Inputs and paths
- Example data paths in the scripts point to OneDrive locations (e.g., `C:\Users\...\gemini survival2.csv`, `pro2_wide.csv`, `longitudinal PO2 VDZ.csv`). Update to local paths before running.
- `TrajectoryVisualization()` expects a data frame with columns: `Subject`, `Cluster`, `Time`, `Value`, `Trajectory` (numeric `Time`/`Value`).

## Dependencies
- Common: `dplyr`, `ggplot2`, `Cairo`.
- Survival: `ggsurvfit`, `survminer`.
- Clustering example: `gbmt`.
Install in R (e.g., `install.packages(c("dplyr","ggplot2","Cairo","ggsurvfit","survminer","gbmt"))`).
