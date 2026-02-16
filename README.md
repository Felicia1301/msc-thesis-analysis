# Cardiac Responses to Feedback

Analysis scripts for the Master's thesis:

**"When Feedback Hits the Heart: Cardiac Responses to Feedback in Internalizing Disorders and Controls"**

This repository contains the R scripts used for preprocessing and statistical analysis of feedback-locked cardiac responses.

---

## Project Overview

The project investigates phasic heart rate responses to feedback in healthy controls and patients with internalizing psychopathology (HiTOP Distress and Fear subdomains).

Cardiac dynamics were analyzed using single-trial linear mixed-effects models across successive interbeat intervals (IBIs) surrounding feedback onset.

---

## Repository Structure

- `01_Demographics.R`  
  Descriptive statistics and group comparisons of demographic and clinical variables.

- `02_Preprocessing_Pipeline.R`  
  Marker cleaning, trial segmentation, IBI computation, feedback alignment, and construction of the final analysis table.

- `03_Analyses_LMM.R`  
  Linear mixed-effects models, assumption checks, estimated marginal means, and figure generation.

- `data/`  
  Not included due to data protection regulations.

- `outputs/`  
  Generated automatically when running the scripts.

---

## Reproducibility Notes

- Raw and processed data are not publicly available due to data protection regulations.
- Scripts are structured to run sequentially:
  
  1. `01_Demographics.R`
  2. `02_Preprocessing_Pipeline.R`
  3. `03_Analyses_LMM.R`

- Required R packages are loaded at the beginning of each script.

---

## Author

Felicia Sperber  
M.Sc. candidate in Clinical Psychology and Psychotherapy

