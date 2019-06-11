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


if (!exists("SGEcologicalDrought_1_DescribeExperiment")) {


#--- Load custom functions
source("Schlaepfer+_SolarGeoengineering_EcologicalDrought_0_Functions.R")

pkgs <- c("raster", "sp")
stopifnot(all(sapply(pkgs, requireNamespace)))


#--- Define and create folder paths
dir_res <- if (requireNamespace("here")) here::here() else normalizePath(".")
dir_prj_in <- file.path(dir_res, "1_Simulation_Setup")
dir.create(dir_prj_in, recursive = TRUE, showWarnings = FALSE)

dir_prj_out <- file.path(dir_res, "2_Simulation_Output")
dir.create(dir_prj_out, recursive = TRUE, showWarnings = FALSE)

dir_res_data <- file.path(dir_res, "3_Prepared_Data")
dir.create(dir_res_data, recursive = TRUE, showWarnings = FALSE)

dir_res2 <- file.path(dir_res, "4_Analysis_Results")
dir.create(dir_res2, recursive = TRUE, showWarnings = FALSE)

dir_res_objMS1 <- dir_res2

ftag_MS1 <- "Obj_MS1_v3"

ftag_gatt <- list(
  title = "Solar geoengineering impacts on dryland water balance",
  mip_era = "CMIP5",
  activity_id = "GeoMIP",
  parent_experiment_id = "RCP4.5",
  source = paste("SOILWAT2 (6f4159755b7c557ceee0eadae1d588223181d262);",
    "rSOILWAT2 (2212524c5f11b2603d9f39f54f447e2f0be97c29);",
    "rSFSW2 (0e2624e79983d54021026c3a257567aba19be2d8)"),
  source_id = "SOILWAT2",
  source_type = "LAND",
  product = "model-output",
  further_info_url = "https://github.com/DrylandEcology/",
  contact = "daniel.schlaepfer@yale.edu",
  grid = "native geographic latitude-longitude grid (576x1152 latxlon)",
  grid_label = "gn",
  nominal_resolution = "50 km",
  realm = "land",
  version = "v20181206"
)

ftag_gatt2 <- list(
  experiment_id = "",
  parent_source_id = "",
  external_variables = "areacella"
)


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#--- LOAD DESCRIPTION OF SIMULATION EXPERIMENT (rSFSW2 object)
fmeta <- file.path(dir_prj_in, "SFSW2_project_descriptions.rds")
fname_gridspec <- file.path(dir_prj_in, "gridspec_fx_SOILWAT2_DAM_gn.nc")

if (file.exists(fmeta)) {
  library("rSFSW2")
  SFSW2_prj_meta <- readRDS(fmeta)

  if (!file.exists(fname_gridspec)) {
    r <- create_raster_from_variables(SFSW2_prj_meta,
      data = SFSW2_prj_meta[["sim_size"]][["runIDs_sites"]])
    create_netCDF_from_raster_with_variables(r,
      var_attributes = list(name = "gridspec",
        standard_name = "grid_specification",
        long_name = "Grid specification", units = "-",
        comment = "Value represents SiteID of dryland study simulation"),
      global_attributes = ftag_gatt,
      file = fname_gridspec)
  }

} else {
  warning("Simulation experiment description not available; script will only ",
    "work if intermediate objects are available in ",
    shQuote(basename(dir_prj_in)), " and ", shQuote(basename(dir_prj_out)))

  #---Create a rudimentary simulation experiment description
  SFSW2_prj_meta <- list()

  SFSW2_prj_meta[["fnames_out"]] <- list(
    dbOutput = ""
  )

  SFSW2_prj_meta[["sim_space"]] <- list()

  # Simulation grid
  SFSW2_prj_meta[["sim_space"]][["sim_raster"]] <- read_netCDF_to_raster(x =
      fname_gridspec, varname = "gridspec")

  # Simulated cells based on simulation grid
  temp <- data.frame(
    sp::coordinates(SFSW2_prj_meta[["sim_space"]][["sim_raster"]]),
    raster::as.data.frame(SFSW2_prj_meta[["sim_space"]][["sim_raster"]]))
  temp <- temp[order(temp[, 3], na.last = NA), 1:2]
  rownames(temp) <- NULL
  colnames(temp) <- c("X_WGS84", "Y_WGS84")
  SFSW2_prj_meta[["sim_space"]][["run_sites"]] <- sp::SpatialPoints(
    coords = temp,
    proj4string = raster::crs(SFSW2_prj_meta[["sim_space"]][["sim_raster"]]))

  SFSW2_prj_meta[["sim_size"]] <- list(
    runsN_sites = length(SFSW2_prj_meta[["sim_space"]][["run_sites"]])
  )

  SFSW2_prj_meta[["sim_scens"]] <- list(
    ambient = "Current",

    reqCSsPerM = list(
      c("RCP45", "RCP85", "G3", "G4"),
      c("RCP45", "RCP85", "G4"),
      c("RCP45", "RCP85", "G3", "G4"),
      c("RCP45", "RCP85", "G3", "G4"),
      c("RCP45", "RCP85", "G3"),
      c("RCP45", "RCP85", "G4"),
      c("RCP45", "RCP85", "G4")),

    reqMs = c("BNU-ESM", "CanESM2", "GISS-E2-R", "HadGEM2-ES", "IPSL-CM5A-LR",
      "MIROC-ESM", "MIROC-ESM-CHEM"),

    DeltaStr_yrs = "d60yrs",

    method_DS = "hybrid-delta-3mod",

    id = c("Current", "hybrid-delta-3mod.d60yrs.RCP45.BNU-ESM",
      "hybrid-delta-3mod.d60yrs.RCP45.CanESM2",
      "hybrid-delta-3mod.d60yrs.RCP45.GISS-E2-R",
      "hybrid-delta-3mod.d60yrs.RCP45.HadGEM2-ES",
      "hybrid-delta-3mod.d60yrs.RCP45.IPSL-CM5A-LR",
      "hybrid-delta-3mod.d60yrs.RCP45.MIROC-ESM",
      "hybrid-delta-3mod.d60yrs.RCP45.MIROC-ESM-CHEM",
      "hybrid-delta-3mod.d60yrs.RCP85.BNU-ESM",
      "hybrid-delta-3mod.d60yrs.RCP85.CanESM2",
      "hybrid-delta-3mod.d60yrs.RCP85.GISS-E2-R",
      "hybrid-delta-3mod.d60yrs.RCP85.HadGEM2-ES",
      "hybrid-delta-3mod.d60yrs.RCP85.IPSL-CM5A-LR",
      "hybrid-delta-3mod.d60yrs.RCP85.MIROC-ESM",
      "hybrid-delta-3mod.d60yrs.RCP85.MIROC-ESM-CHEM",
      "hybrid-delta-3mod.d60yrs.G3.BNU-ESM",
      "hybrid-delta-3mod.d60yrs.G3.GISS-E2-R",
      "hybrid-delta-3mod.d60yrs.G3.HadGEM2-ES",
      "hybrid-delta-3mod.d60yrs.G3.IPSL-CM5A-LR",
      "hybrid-delta-3mod.d60yrs.G4.BNU-ESM",
      "hybrid-delta-3mod.d60yrs.G4.CanESM2",
      "hybrid-delta-3mod.d60yrs.G4.GISS-E2-R",
      "hybrid-delta-3mod.d60yrs.G4.HadGEM2-ES",
      "hybrid-delta-3mod.d60yrs.G4.MIROC-ESM",
      "hybrid-delta-3mod.d60yrs.G4.MIROC-ESM-CHEM")
  )
}



#------------------------------------------------------------------------------#

#--- ASSIGNING GEOGRAPHIC REGIONS
fname_region <- file.path(dir_prj_in, "region_fx_SOILWAT2_DAM_gn.nc")

name_regions <- c("Global",
  "Australia", "Central Asia", "East Asia", "Mediterranean",
  "North America", "South America", "South Asia", "Southern Africa",
  "Subsaharan Africa", "Western Asia")

if (file.exists(fname_region)) {
  r <- read_netCDF_to_raster(x = fname_region, varname = "region")
  temp <- raster::extract(r, SFSW2_prj_meta[["sim_space"]][["run_sites"]])
  region_data <- factor(name_regions[-1][temp])

} else {
  SREX2012 <- ADA_SREX2012_regions(SFSW2_prj_meta)
  region_data <- SREX2012[["adjusted"]][["region_data"]]

  r <- create_raster_from_variables(SFSW2_prj_meta,
    data = as.integer(region_data))
  create_netCDF_from_raster_with_variables(r,
    var_attributes = list(name = "region",
      units = "-",
      comment = paste("Represents geographic region of cells included",
        "in dryland study")),
    global_attributes = ftag_gatt,
    file = fname_region)
}



#--- SIMULATION DEFINITION OF STUDY AREA
# Simulated cells = {cell ij:
#   AI(ij; mean(1961-1990)) in ]0.05, 0.5[ OR
#   AI(ij; mean(2071-2100)) in ]0.05, 0.5[
#   }
#
# based on data layers shared by Huang et al. 2015
#
# Huang, J., H. Yu, X. Guan, G. Wang, and R. Guo. 2015. Accelerated dryland
#   expansion under climate change. Nature Climate Change 6:166–171.
# Chen, M. Y., Xie, P. P., Janowiak, J. E. & Arkin, P. A. 2002. Global land
#   precipitation: A 50-yr monthly analysis based on gauge observations.
#   J. Hydrometeorol. 3, 249–266.
# Fan, Y. & Dool, H. V. D. 2008. A global monthly land surface air temperature
#   analysis for 1948–present. J. Geophys. Res. 113, D01103.

#--- ADJUSTED DEFINITION OF STUDY AREA
fname_landmask <- file.path(dir_prj_in, "sftlf_fx_SOILWAT2_DAM_gn.nc")

if (file.exists(fname_landmask)) {
  r <- read_netCDF_to_raster(x = fname_landmask, varname = "sftlf")
  temp <- raster::extract(r, SFSW2_prj_meta[["sim_space"]][["run_sites"]])
  id_ASA <- as.logical(temp)

} else {
  vars0 <- c("MAP_mm_mean", "UNAridityIndex_Normals_none_mean",
    "TeeriEtAl1976_NSadj_TempAirMin_7thMonth_C_mean")

  reqSCs <- c("Current", "RCP45")
  req_scens <- unlist(lapply(reqSCs, grep,
    x = SFSW2_prj_meta[["sim_scens"]][["id"]], value = TRUE))

  datSA <- extract_dbOut_to_df(SFSW2_prj_meta = SFSW2_prj_meta,
    variables = vars0, scenarios = req_scens,
    whereClause = "Experimental_Label='Default'",
    file = file.path(dir_res_data, "temp_dbOut", "StudyArea_v2.rds"),
    verbose = TRUE)

  # identify cells that are arid and/or semi-arid under current or
  # any GCM under RCP4.5 and that
  use_cond <- grep(paste(reqSCs, collapse = "|"), dimnames(datSA)[[3L]])
  id_ASA1 <- apply(datSA[, "UNAridityIndex_Normals_none_mean", use_cond], 1L,
    function(x) any(x > 0.05 & x <= 0.5))
  id_ASA2 <- apply(
    datSA[, "TeeriEtAl1976_NSadj_TempAirMin_7thMonth_C_mean", use_cond], 1L,
    function(x) all(x >= 0))
  id_ASA3 <- apply(datSA[, "MAP_mm_mean", use_cond], 1L,
    function(x) all(x < 10000))
  id_ASA4 <- !is.na(region_data) # 21 gridcells are not in any geographic region

  # subset to study area
  id_ASA <- id_ASA1 & id_ASA2 & id_ASA3 & id_ASA4

  r <- create_raster_from_variables(SFSW2_prj_meta,
    data = as.integer(id_ASA))
  create_netCDF_from_raster_with_variables(r,
    var_attributes = list(name = "sftlf",
      standard_name = "land_area_fraction",
      long_name = "Land Area Fraction", units = "%",
      comment = "Represents 0/1-mask of cells included in dryland study"),
    global_attributes = ftag_gatt,
    file = fname_landmask)
}



#--- AREAL EXTENT OF RASTER CELLS
fname_areacella <- file.path(dir_prj_in, "areacella_fx_SOILWAT2_DAM_gn.nc")

if (file.exists(fname_areacella)) {
  r <- read_netCDF_to_raster(x = fname_areacella, varname = "areacella")
  cells <- raster::extract(r, SFSW2_prj_meta[["sim_space"]][["run_sites"]])
  cells <- matrix(c(cells, cells / max(cells)), ncol = 2,
    dimnames = list(NULL, c("km2", "rel")))

} else {
  cells <- calculate_cell_area(
    sp_sites = SFSW2_prj_meta[["sim_space"]][["run_sites"]],
    sim_raster = SFSW2_prj_meta[["sim_space"]][["sim_raster"]])

  r <- create_raster_from_variables(SFSW2_prj_meta, data = cells[, "km2"])
  create_netCDF_from_raster_with_variables(r,
    var_attributes = list(name = "areacella",
      standard_name = "cell_area",
      long_name = "Atmosphere Grid-Cell Area", units = "km2"),
    global_attributes = ftag_gatt,
    file = fname_areacella)
}


if (is.null(ftag_gatt[["nominal_resolution"]]) ||
    nchar(ftag_gatt[["nominal_resolution"]]) == 0) {

  ftag_gatt[["nominal_resolution"]] <- calculate_nominal_resolution(
      grid = SFSW2_prj_meta[["sim_space"]][["sim_raster"]],
      sites = SFSW2_prj_meta[["sim_space"]][["run_sites"]],
      cell_areas_km2 = cells[, "km2"])
}



#--- CLIMATE ZONES: TEMPERATE, SUBTROPICAL, TROPICAL
fname_ClimZones <- file.path(dir_prj_in, "climzone_fx_SOILWAT2_DAM_gn.nc")

name_climzones <- c("Global",
  "boreal", "polar", "subtropical", "temperate", "tropical")

if (file.exists(fname_ClimZones)) {
  r <- read_netCDF_to_raster(x = fname_ClimZones, varname = "climzone")
  temp <- raster::extract(r, SFSW2_prj_meta[["sim_space"]][["run_sites"]])
  climzones_data <- factor(name_climzones[-1][temp])

} else {
  varsCZ <- c(paste0("TempAir_m", 1:12, "_C_mean"))

  datsTmon <- extract_dbOut_to_df(SFSW2_prj_meta = SFSW2_prj_meta,
    variables = varsCZ, scenarios = "Current",
    whereClause = "Experimental_Label='Default'",
    file = file.path(dir_res_data, "temp_dbOut", "MontlyTmean.rds"),
    verbose = TRUE)

  temp <- apply(datsTmon, c(1, 3), trewartha_climate)
  climzones_data <- factor(temp[, "Current"])

  r <- create_raster_from_variables(SFSW2_prj_meta,
    data = as.integer(climzones_data))
  create_netCDF_from_raster_with_variables(r,
    var_attributes = list(name = "climzone",
      units = "-",
      long_name = "Trewartha climate zones",
      comment = paste("Represents climate zones of cells under historical",
        "conditions included in dryland study")),
    global_attributes = ftag_gatt,
    file = fname_ClimZones)
}




#------------------------------------------------------------------------------#
#--- DEFINE SCENARIOS AND CLIMATE CONDITIONS OF SIMULATION EXPERIMENT

sc_current <- SFSW2_prj_meta[["sim_scens"]][["ambient"]]
reqCSs <- c("RCP45", "G3", "G4") #SFSW2_prj_meta[["sim_scens"]][["reqCSs"]]
reqACSs <- c(sc_current, reqCSs)
scen_SARs <- grep("G[[:digit:]]", reqCSs, value = TRUE)

deltas_byCS <- c(RCP45 = sc_current, RCP85 = sc_current,
  G3 = "RCP45", G4 = "RCP45")
deltas_byCS <- deltas_byCS[reqCSs]

req_scens <- unname(unlist(sapply(reqCSs, function(x)
  grep(x, SFSW2_prj_meta[["sim_scens"]][["id"]], value = TRUE))))

req_Downs <- unique(sapply(strsplit(req_scens, split = ".", fixed = TRUE),
  function(x) x[[1]]))
req_dTime <- c("d0yrs",
  unique(sapply(strsplit(req_scens, split = ".", fixed = TRUE),
    function(x) x[[2]])))

reqGCMs_perCS <- lapply(reqCSs, function(cs) {
  temp <- sapply(SFSW2_prj_meta[["sim_scens"]][["reqCSsPerM"]],
    function(x) any(x == cs))
  SFSW2_prj_meta[["sim_scens"]][["reqMs"]][temp]
})
names(reqGCMs_perCS) <- reqCSs
reqMs <- sort(unique(unlist(reqGCMs_perCS)))

temp <- c(sc_current, unlist(
  sapply(reqCSs, function(cs) sapply(reqGCMs_perCS[[cs]], function(x)
    grep(paste0(cs, ".", x, "$"), req_scens, value = TRUE)))))
temp <- SFSW2_prj_meta[["sim_scens"]][["id"]] %in% temp
sim_scens_id <- SFSW2_prj_meta[["sim_scens"]][["id"]][temp]



} # endif not exists "SGEcologicalDrought_1_DescribeExperiment"

SGEcologicalDrought_1_DescribeExperiment <- TRUE


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
