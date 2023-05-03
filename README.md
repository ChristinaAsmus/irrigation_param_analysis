# irrigation_param_analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7889385.svg)](https://doi.org/10.5281/zenodo.7889385)
[![test](https://github.com/christinaasmus/irrigation_param_analysis/actions/workflows/test.yaml/badge.svg)](https://github.com/christinaasmus/irrigation_param_analysis/actions/workflows/test.yaml)

.. image:: 
   :target: 
This repository contains the code for the analysis of the irrigation parameterization implemented into REMO2020-iMOVE. 

## Figure overview
Figure |  directory | plotting routine
--- | --- | --- |
Fig003 | experiment_setup | Fig03_Model_domain_irrigated_fraction_spatial.py
Fig004 | testing_schemes | Fig04_Test_schemes_irrigation_water_hourly_development.py
Fig005 | testing_schemes | Fig05_Test_schemes_irrigation_water_spatial.py
Fig006 | meteo_condition | Fig06_Meteo_condition_spatial.py
Fig007 | meteo_condition | Fig07_Meteo_condition_AMJJA_series.py
Fig008 | process_analysis | Fig08_Process_analysis_soil_surface_vars.py
Fig009 | process_analysis | Fig09_Process_analysis_soil_surface_fluxes_spatial.py
Fig010 | process_analysis | Fig10_Process_analysis_surface_energy_budget_diurnal_series.py
Fig011 | process_analysis | Fig11_Process_analysis_t2m_rh.py
Fig012 | process_analysis | Fig12_Process_analysis_precip.py
Fig013 | process_analysis | Fig13_Process_analysis_LAI_NPP_AMJ_spatial.py
Fig014 | process_analysis | Fig14_Process_analysis_LAI_NPP_series.py
Fig015 | delayed_effects | Fig15_Delayed_effects_t2m_spatial.py
Fig016 | delayed_effects | Fig16_Delayed_effects_heatwave_combined_series.py
FigA2 | appendix | A2_schemes_soil_moisture_runoff_drainage_series.py
FigA3 | appendix | A3_A4_schemes_effects_soil_moisture_surface_temp.py
FigA4 | appendix | A3_A4_schemes_effects_soil_moisture_surface_temp.py
FigB1 | appendix | B1_irrigation_Aug.py
FigC1 | appendix | C1_station_locations_spatial_scatter.py

## Installation

Create a new conda environment using, e.g.

```bash
conda env create -f environment.yaml
```

Afterwards, you can install the required plotting and `analysis_functions` using
```
pip install -e .
```

## Data access
To reproduce the analysis, please download the data from https://doi.org/10.5281/zenodo.7867329 using download_zenodo.py.
The observational data for the analysis in /compare_obs/Tab3_Eval_Scia_t2m.py and in /appendix/C1_station_locations_spatial_scatter.py can be downloaded from http://193.206.192.214/serverstazioni/stazioni400.php.
