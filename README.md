# Pre-post intervention in stroke patients: TMS data analysis pipeline

From Alica. All rights reserved.

This pipeline runs on MELGS1 (for the correct Matlab version), and then on a local computer for the statistical analysis in R.

A copy of all of these files can be found under the main project directory on the server `/rri_disks/eugenia/meltzer_lab/TMS_stroke`.
All pertinent project information can be found on the meltzer_lab dropbox under `meltzer_lab_docs/other_docs/16-43_TMS_Stroke_motorintervention_Alica`

# TMS data preprocessing in Matlab

Matlab preprocessing scripts are written for each of the TMS protocols: SICI, LICI, RC, IHI, cSP, and iSP.
See `matlab_scripts_walkthrough` for a brief overview on how to run these scripts.

# Statistical analysis in R

Processing scripts for data cleaning and statistical analysis have been written in R for the SICI, LICI, and RC data. 
Scripts are thoroughly commented and should be easy to follow. See `statistical_analysis_walkthrough` for a brief overview of what the scripts do. 
