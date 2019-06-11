#!/usr/bin/env Rscript

################################################################################
#    Solar geoengineering may increase climate change droughts in drylands
#
#    Daniel R Schlaepfer1*, William K Lauenroth1, John B Bradford2,
#    Kyle A Palmquist3
#
#    * Corresponding author: daniel.schlaepfer@yale.edu
#
#    DRS: https://orcid.org/0000-0001-9973-2065
#    WKL: https://orcid.org/0000-0002-3079-4484
#    JBB: https://orcid.org/0000-0001-9257-6303
#    KAP: https://orcid.org/0000-0001-7665-4105
#
#    1 Yale School of Forestry and Environmental Studies, Yale University,
#      New Haven, CT 06511
#    2 Southwest Biological Science Center, U.S. Geological Survey,
#      Flagstaff, AZ 86011
#    3 Department of Biological Sciences, Marshall University,
#      Huntington, WV 25755
#
################################################################################
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################


if (!exists("SGEcologicalDrought_2_PrepareData")) {

#--- USER INPUTS
do_data_sharing <- TRUE # Produce data tables for sharing


#--- Load packages, custom functions, and description of simulation experiment
source("Schlaepfer+_SolarGeoengineering_EcologicalDrought_1_DescribeExperiment.R")

pkgs <- c("raster", "sp", "RNetCDF")
stopifnot(all(sapply(pkgs, requireNamespace)))


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#--- Data slices
kdatas <- c(val = "dats_MS1", ens = "dats_Ens")

#--- Extract data for MS1 analysis
mo <- seq_len(12)

vars_MS1 <- list(
  grid = c("site_id"),

  # Target response variables
  res = c(
    "DrySoilPeriods_SWPcrit3000kPa_NSadj_topLayers_AllLayersDry_Duration_LongestContinuous_days_mean",
    "DrySoilPeriods_SWPcrit3000kPa_NSadj_bottomLayers_AllLayersDry_Duration_LongestContinuous_days_mean",
    "SoilPeriods_WarmDry_allLayers_TcritPos35C_SWPcrit3000kPa_Count_days_mean"
  ),

  # Supplementary variables: climate
  clim = c("MAT_C_mean", "MAP_mm_mean")
)


vars_MS1[["res2"]] <- vars_MS1[["res"]][1:3]
vars_MS1[["MATMAP"]] <- vars_MS1[["clim"]][1:2]
vars_all <- unique(unlist(vars_MS1))
timeaggs_MS1_all <- sapply(strsplit(vars_all, split = "_"), function(x)
  x[length(x)])

vset_main <- "res2"
vset_clim <- "MATMAP"

names_MS1 <- list(
  grid = c("Site ID"),

  # Target response variables
  res = c(
    "Shallow soil drought",
    "Deep soil drought",
    "Hot drought"),

  # Supplementary variables: climate
  clim = c("MAT", "MAP")
)

units_MS1 <- list(
  grid = c("-"),

  # Target response variables
  res = c(
    "days",
    "days",
    "days"),

  # Supplementary variables: climate
  clim = c("C", "mm")
)

id <- match(vars_all, unlist(vars_MS1))
names_MS1_all <- unlist(names_MS1)[id]
units_MS1_all <- unlist(units_MS1)[id]

labs_MS1 <- list()
for (ns in names(names_MS1)) {
  labs_MS1[[ns]] <- paste0(names_MS1[[ns]],
    ifelse(nchar(units_MS1[[ns]]) > 0 & units_MS1[[ns]] != "-",
      paste0(" (", units_MS1[[ns]], ")"), ""))
}


labs_MS1[["res2"]] <- labs_MS1[["res"]][1:3]
labs_MS1[["MATMAP"]] <- labs_MS1[["clim"]][1:2]
labs_all <- unique(unlist(labs_MS1))


# 1, more is worse (drier/hotter); -1, more is better; NA, no direction
worse_direction <- list(
  grid = NA,
  res = c(1, 1, 1),
  clim = c(1, -1)
)

worse_direction[["res2"]] <- worse_direction[["res"]][1:3]
worse_direction[["MATMAP"]] <- worse_direction[["clim"]][1:2]
id <- match(vars_all, unlist(vars_MS1))
worse_direction_all <- unlist(worse_direction)[id]


#--- Load simulation data
print(paste(Sys.time(), "--", "Load simulation data"))

dats_MS1 <- list(
  vals = load_rSFSW2_data_for_analysis(
    meta = SFSW2_prj_meta,
    path = dir_prj_out,
    path_tmp = dir_res_data,
    whereClause = "Experimental_Label='Default'",
    variables = vars_all,
    sim_scenarios = sim_scens_id,
    sc_historical = sc_current,
    subset = id_ASA,
    write_to_netcdf = do_data_sharing,
    ftag_gatt = ftag_gatt, ftag_gatt2 = ftag_gatt2,
    var_names = names_MS1_all,
    var_units = units_MS1_all,
    timeaggs = timeaggs_MS1_all
  ))


} # endif not exists "SGEcologicalDrought_2_PrepareData"

SGEcologicalDrought_2_PrepareData <- TRUE


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
