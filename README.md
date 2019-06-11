
# Solar geoengineering may increase climate change droughts in drylands

## Article DOI


## Authors and affiliation

[Daniel R Schlaepfer](https://orcid.org/0000-0001-9973-2065),
[William K Lauenroth](https://orcid.org/0000-0002-3079-4484),
[John B Bradford](https://orcid.org/0000-0001-9257-6303), and
[Kyle A Palmquist](https://orcid.org/0000-0001-7665-4105)



## Code availability
The simulation model SOILWAT2 is available at
* https://github.com/DrylandEcology/SOILWAT2
* https://github.com/DrylandEcology/rSOILWAT2
* https://github.com/DrylandEcology/rSFSW2

and the code that we used to produce figures and tables that support the
findings of this study are available (__here__) at
* https://github.com/dschlaep/SolarGeoengineering_EcologicalDrought2019



## Short overview of code organization

The R script to analyze our simulation experiment and to produce figures and
tables consists of five files:
1. `Schlaepfer+_SolarGeoengineering_EcologicalDrought_0_Functions.R` which
   contains custom functions

2. `Schlaepfer+_SolarGeoengineering_EcologicalDrought_1_DescribeExperiment.R`
   which defines variables that describe our simulation experiment

3. `Schlaepfer+_SolarGeoengineering_EcologicalDrought_2_PrepareData.R` which
   loads and prepares simulation output

4. `Schlaepfer+_SolarGeoengineering_EcologicalDrought_3_Analysis.R` which
   carries out data analysis

5. `Schlaepfer+_SolarGeoengineering_EcologicalDrought_4_ProduceOutput.R` which
   produces figures and tables


## Data availability
All climate model and GeoMIP data used in this study are available through
the Earth System Grid Federation network (https://esgf.llnl.gov) and through
the Geoengineering Model Intercomparison Project
(http://climate.envsci.rutgers.edu/GeoMIP).

The simulation output produced by this study that support the findings of
this study are available in https://figshare.com with the identifier
https://doi.org/10.6084/m9.figshare.8258519.

Before the code can be run, data need to be downloaded and placed at the
correct paths:

* The folder `1_Simulation_Setup/` should contain the following five files
  in netCDF format which can be downloaded from
  https://doi.org/10.6084/m9.figshare.8258519:
  - `gridspec_fx_SOILWAT2_DAM_gn.nc`: specification of simulation gridcells
  - `sftlf_fx_SOILWAT2_DAM_gn.nc`: land mask of simulated area
  - `areacella_fx_SOILWAT2_DAM_gn.nc`: area (km2) of each gridcell
  - `climzone_fx_SOILWAT2_DAM_gn.nc`: assigning each gridcell to a climate zone
  - `region_fx_SOILWAT2_DAM_gn.nc`: assigning each gridcell to a geographic region

* The folder `2_Simulation_Output/` should contain the following files
  in netCDF format which can be downloaded from
  https://doi.org/10.6084/m9.figshare.8258519:
  - `All_LyrC_SOILWAT2-BNU-ESM_G3_gn_2040-2070.nc`
  - `All_LyrC_SOILWAT2-BNU-ESM_G4_gn_2040-2070.nc`
  - `All_LyrC_SOILWAT2-BNU-ESM_RCP45_gn_2040-2070.nc`
  - `All_LyrC_SOILWAT2-CanESM2_G4_gn_2040-2070.nc`
  - `All_LyrC_SOILWAT2-CanESM2_RCP45_gn_2040-2070.nc`
  - `All_LyrC_SOILWAT2-GISS-E2-R_G3_gn_2040-2070.nc`
  - `All_LyrC_SOILWAT2-GISS-E2-R_G4_gn_2040-2070.nc`
  - `All_LyrC_SOILWAT2-GISS-E2-R_RCP45_gn_2040-2070.nc`
  - `All_LyrC_SOILWAT2-HadGEM2-ES_G3_gn_2040-2070.nc`
  - `All_LyrC_SOILWAT2-HadGEM2-ES_G4_gn_2040-2070.nc`
  - `All_LyrC_SOILWAT2-HadGEM2-ES_RCP45_gn_2040-2070.nc`
  - `All_LyrC_SOILWAT2-historical_gn_1980-2010.nc`
  - `All_LyrC_SOILWAT2-IPSL-CM5A-LR_G3_gn_2040-2070.nc`
  - `All_LyrC_SOILWAT2-IPSL-CM5A-LR_RCP45_gn_2040-2070.nc`
  - `All_LyrC_SOILWAT2-MIROC-ESM_G4_gn_2040-2070.nc`
  - `All_LyrC_SOILWAT2-MIROC-ESM_RCP45_gn_2040-2070.nc`
  - `All_LyrC_SOILWAT2-MIROC-ESM-CHEM_G4_gn_2040-2070.nc`
  - `All_LyrC_SOILWAT2-MIROC-ESM-CHEM_RCP45_gn_2040-2070.nc`

  Each of these files contains output from the simulation model `SOILWAT2` that
  was forced by historical observations or by projected climate data
  from a GCM that was run under the climate change experiment RCP4.5
  or under a solar geoengineering experiment G3 or G4. The files contain values
  for three output variables:
  - "Shallow soil drought (days)"
  - "Deep soil drought (days)"
  - "Hot drought (days)"

  and two aggregated input variables:
  - "MAT (C)", mean annual temperature
  - "MAP (mm)", mean annual precipitation
