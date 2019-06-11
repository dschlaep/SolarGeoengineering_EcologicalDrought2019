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


#--- USER INPUTS

# Output to `Output_Background`
do_out_studyarea <- TRUE
do_out_mod_illustration <- TRUE

# Output to `Output_Background`: do_out_matmap
# Figure 1 | Cumulative density of change in mean annual climate inputs across global drylands.
# --> `Output_Background/Fig_dCDFs_Global_GCMs_by_Scenarios_MAPandMAT_refCurrent_EnsCW-mean_v2.pdf`
do_out_matmap <- TRUE

# Output to `Output_Ensemble`
do_out_ensembles <- TRUE

# Output to `Output_Response_Maps`: do_out_mod_maps
# Figure 2 | Drought responses under SSAI scenario G4 based on runs forced by BNU-ESM climate projections.
# --> `Output_Response_Maps/res2_G4/Fig_Map_dGiandMods_G4toRCP45vsCurrent_res2_GMods2c_tol0p500_BNU-ESM.png`
# Figure 3 | Drought responses under SSAI scenario G4 based on ensemble of runs.
# --> `Output_Response_Maps/res2_G4/Fig_Map_dGiandMods_G4toRCP45vsCurrent_res2_GMods2c_tol0p500_EnsCW_mean_mean.png`
do_out_mod_maps <- TRUE

# # Output to `Output_Response_Barplots`: do_out_mod_bp
# Figure 4 | Extent by dryland region of at least 50% modification by SSAI of baseline climate change.
# --> `Output_Response_Barplots/res2/Fig_Barplot_dGiandMods_G3G4toRCP45vsCurrent_res2_GMods2c_tol0p500_AllGCMs_v3_EnsCW_mean.pdf`
do_out_mod_bp <- TRUE

do_out_rel_bp <- TRUE

# Output to `Output_Response_CAFs`
do_out_CAFs <-  FALSE

# Output to `Output_VarPartitioning`
do_out_variancepartition <- TRUE

# Output to `Output_SpatialPatterns`
do_out_spatialpattern <- TRUE



#--- Load packages, custom functions, description of simulation experiment,
# prepared data, and analysis
source("Schlaepfer+_SolarGeoengineering_EcologicalDrought_3_Analysis.R")

stopifnot(requireNamespace("reshape2"))


#------------------------------------------------------------------------------#
#--------- PRODUCTS
#------------------------------------------------------------------------------#
print(paste(Sys.time(), "--", "Produce output"))

flag_main <- paste0("main=", vset_main)



#------------------------------------------------------------------------------#
#--- Ensemble GCMs
if (do_out_ensembles) {
  dir_ens <- file.path(dir_res_objMS1, "Output_Ensemble")
  dir.create(dir_ens, recursive = TRUE, showWarnings = FALSE)


  plot_GlobalEnsembleValues(data = vals_Ens2, data2 = valsCT_Ens2,
    ens_method = "relative", ens_wise = "EnsCW", ens_fun = "mean",
    ens_funCT = ens_funCT, ens_outs = ens_outs, ens_types = ens_types,
    labs_ens_types = labs_ens_types, scen_SARs = scen_SARs, reqMs = reqMs,
    path = dir_ens, ftag = flag_main)



  for (ew in ens_wises) {
    for (em in ens_methods) {
      write.csv(vals_Ens2[em, ew][[1]], file = file.path(dir_ens,
        paste0("Table_GlobalEnsembleValues_", em, "_", ew, "_",
          flag_main, ".csv")))

      write.csv(ranks_Ens2[em, ew][[1]], file = file.path(dir_ens,
        paste0("Table_GlobalEnsembleRanks_", em, "_", ew, "_",
          flag_main, ".csv")))

      if (em == "relative") {
        for (ef in ens_funCT) {
          write.csv(valsCT_Ens2[paste0(em, "_", ef), ew][[1]], file =
              file.path(dir_ens, paste0("Table_GlobalEnsembleValues_",
                em, "_", ew, "-", ef, "_", flag_main, ".csv")))
        }
      }
    }
  }

}


#------------------------------------------------------------------------------#
#--- Study area with regions and climate zones

if (do_out_studyarea) {
  dir_back <- file.path(dir_res_objMS1, "Output_Background")
  dir.create(dir_back, recursive = TRUE, showWarnings = FALSE)

  # Map of study area
  SREX2012 <- ADA_SREX2012_regions(SFSW2_prj_meta)

  plot_map_studyarea(
    meta = SFSW2_prj_meta, climzones_data = climzones_data,
    subset = id_ASA, add_SREX2012 = exists("SREX2012"),
    SREX2012 = if (exists("SREX2012")) SREX2012 else NULL,
    dtemp = dir_back)


  # Table of areal extent
  table_extent <- tabulate_by_region_x_climzone(x = cells, vars = "km2",
    subset = id_ASA, fmethod = "sum",
    region_data = region_data, climzones_data = climzones_data)

  table_extent2 <- reshape2::dcast(table_extent, Region ~ ClimateZone,
    value.var = "km2_sum")

  write.csv(table_extent2,
    file = file.path(dir_back, "Table_StudyArea_Extent.csv"), row.names = FALSE)


  #--- Extent of (arid and semiarid) study area (under current conditions)
  # in km2 and % of global terrestrial surface
  extent_km2 <- sum(cells[id_ASA, "km2"], na.rm = TRUE) # 37644769

  # Literature
  # Global terrestrial surface: http://en.wikipedia.org/wiki/Earth;
  # TODO: find better reference
  extent_terrestrial_km2 <- 148940000 # km2 from
  round(100 * extent_km2 / extent_terrestrial_km2, 2) # 25.28%
}



#------------------------------------------------------------------------------#
#--- Global MAP & MAT density plots

if (do_out_matmap) {
  dir_back <- file.path(dir_res_objMS1, "Output_Background")
  dir.create(dir_back, recursive = TRUE, showWarnings = FALSE)

  col_byCSs <- c(RCP45 = "orange", RCP85 = "firebrick4",
    G3 = "darkorchid4", G4 = "deeppink3")
  col_byCSs <- col_byCSs[reqCSs]

  ef <- "mean"
  vars <- vars_MS1[[vset_clim]]

  # Cumulative density plots
  temp <- rep(sc_current, length(reqCSs))
  names(temp) <- reqCSs
  used_deltas_byCS <- list(deltas_byCS, temp)

  for (k in seq_along(used_deltas_byCS)) {
    for (ew in ens_wises) {
      # Prepare data
      data_dCDFs <- dataEns_dCDFs <- list()

      for (sc in reqCSs) {
        tag <- paste0("delta_", used_deltas_byCS[[k]][sc], "_abs")

        iscens <- grep(sc, sim_scens_id)
        data_dCDFs[[sc]] <- dats_MS1[[tag]][, vars, iscens]

        iscen <- grep(paste(sc, ef, sep = "."),
          dimnames(dats_Ens[[ew]][[tag]])[[3]])
        dataEns_dCDFs[[sc]] <- dats_Ens[[ew]][[tag]][, vars, iscen]
      }

      plot_dCDFs_global_v2(data = data_dCDFs, dataEns = dataEns_dCDFs,
        vars = vars, var_labels = labs_MS1[[vset_clim]],
        cell_within_study_TF = id_ASA, reqCSs = reqCSs,
        deltas_byCS = used_deltas_byCS[[k]],
        reqMs = reqMs, col_byCSs = col_byCSs,
        path = dir_back,
        ftag = paste0("MAPandMAT_ref",
          paste(unique(used_deltas_byCS[[k]][reqCSs]), collapse = ""), "_",
          ew, "-", ef, "_v2"))
    }
  }

  # Table of ensemble median GCM
  if (all(c("EnsGW", "EnsRW") %in% names(dats_Ens))) {
    fname <- file.path(dir_back, paste0("Table_Climate_MATMAP_Ens-", ef,
      ".csv"))

    table_matmap <- NULL

    for (sc in reqCSs) {
      tag <- paste0("delta_", deltas_byCS[sc], "_abs")

      iscen <- grep(paste(sc, ef, sep = "."),
        dimnames(dats_Ens[["EnsGW"]][[tag]])[[3]])

      temp <- tabulate_by_region_x_climzone(
        x = dats_Ens[["EnsGW"]][[tag]][, vars, iscen],
        x_regional = dats_Ens[["EnsRW"]][[tag]][, vars, iscen],
        vars = vars, subset = id_ASA, fmethod = "quantile",
        area = cells[, "km2"],
        region_data = region_data, climzones_data = climzones_data)

      temp <- data.frame(Scenario = sc, Reference = deltas_byCS[sc], temp,
        stringsAsFactors = FALSE, row.names = NULL)

      table_matmap <- rbind(table_matmap, temp)
    }

    write.csv(table_matmap, file = fname, row.names = FALSE)
  }
}



#------------------------------------------------------------------------------#
#--- Maps of response variables and of modifications
#------------------------------------------------------------------------------#
if (do_out_mod_maps) {
  ef <- "mean"
  adhoc <- FALSE

  for (k in seq_along(scen_SARs)) {
    sc_name <- scen_SARs[k]
    ref1 <- deltas_byCS[sc_name]
    temp <- paste(paste0(sc_name, ".", reqGCMs_perCS[[sc_name]]),
      collapse = "|")
    gcms <- grep(temp, sim_scens_id, value = TRUE)
    ens <- paste("hybrid-delta-3mod.d60yrs", sc_name, ef, sep = ".")

    for (vset in c(vset_clim, vset_main)) { #names(vars_MS1)[-1]) {#
      dir_resp_out_tempM <- file.path(dir_res_objMS1, "Output_Response_Maps",
        paste0(vset, "_", sc_name))
      dir.create(dir_resp_out_tempM, recursive = TRUE, showWarnings = FALSE)

      for (ref0 in "Current") {# unique(deltas_byCS)) {
        dtemp1 <- paste0("delta_", ref0, "_abs")

        for (mwa in FALSE) { #c(FALSE, TRUE)) {#
          for (idt in seq_len(nrow(modtols))) {
            dtemp2 <- paste0("delta_", ref1, "_dir", modtols[idt, "tag"])
            ftemp <- paste0(sc_name, "to", ref1, "vs", ref0, "_", vset,
              if (mwa) "_GModsAll" else "_GMods2c", modtols[idt, "tag"])

            # Plots for each GCM
            plot_map_vars_dGi_and_modifier(
              x = dats_MS1[[dtemp1]][, vars_MS1[[vset]], gcms, drop = FALSE],
              xmod = dats_MS1[[dtemp2]][, vars_MS1[[vset]], gcms, drop = FALSE],
              meta = SFSW2_prj_meta, subset = id_ASA,
              GCMs = reqGCMs_perCS[[sc_name]], exp = sc_name,
              ref0 = ref0, ref1 = ref1,
              mod_with_all = mwa,
              mod_tol = modtols[idt, "tol"],
              expected_dy_lt_dx = worse_direction[[vset]],
              var_labels = labs_MS1[[vset]],
              path = dir_resp_out_tempM, ftag = ftemp)

            # Plots for ensemble across GCMs
            for (ew in ens_wises) {
              if (adhoc) {
                #TODO: use instead dats_Ens3 <- calculate_ensembles() once fixed
                temp <- dats_MS1[[dtemp2]][, vars_MS1[[vset]], gcms, drop = FALSE]
                xmod_ens <- calc_cellwise_ensemble(SFSW2_prj_meta,
                  data = temp, subset = id_ASA, funs = "majority",
                  na.rm = TRUE, verbose = FALSE)
                temp <- reshape2::melt(xmod_ens)
                temp[, "value"] <- as.character(temp[, "value"])
                xmod_ens <- reshape2::acast(temp, Var2 ~ Var3 ~ Var4 + Var1)
                dimnames(xmod_ens)[[3]] <- ens

              } else {
                xmod_ens <- dats_Ens[[ew]][[dtemp2]][, vars_MS1[[vset]], ens, drop = FALSE]
              }

              plot_map_vars_dGi_and_modifier(
                x = dats_Ens[[ew]][[dtemp1]][, vars_MS1[[vset]], ens, drop = FALSE],
                xmod = xmod_ens,
                meta = SFSW2_prj_meta, subset = id_ASA,
                GCMs = ef, exp = sc_name,
                ref0 = ref0, ref1 = ref1,
                mod_with_all = mwa,
                mod_tol = modtols[idt, "tol"],
                expected_dy_lt_dx = worse_direction[[vset]],
                var_labels = labs_MS1[[vset]],
                path = dir_resp_out_tempM,
                ftag = paste0(ftemp, "_", ew, "_", ef, if (adhoc) "-adhoc"))
            }
          }
        }
      }
    }
  }
}


if (do_out_mod_maps) {
  ef <- "mean"

  for (k in seq_along(1)) {
    sc_name <- "RCP45"
    ref1 <- deltas_byCS[sc_name]
    temp <- paste(paste0(sc_name, ".", reqGCMs_perCS[[sc_name]]),
      collapse = "|")
    gcms <- grep(temp, sim_scens_id, value = TRUE)
    ens <- paste("hybrid-delta-3mod.d60yrs", sc_name, ef, sep = ".")

    for (vset in c(vset_clim, vset_main)) { #names(vars_MS1)[-1]) {#
      dir_resp_out_tempM <- file.path(dir_res_objMS1, "Output_Response_Maps2",
        paste0(vset, "_", sc_name))
      dir.create(dir_resp_out_tempM, recursive = TRUE, showWarnings = FALSE)

      for (ref0 in "Current") {# unique(deltas_byCS)) {
        dtemp1 <- "vals" #paste0("delta_", ref0, "_abs")

        for (mwa in FALSE) { #c(FALSE, TRUE)) {#
          for (idt in 1) {
            dtemp2 <- paste0("delta_", ref1, "_abs")
            ftemp <- paste0(sc_name, "to", ref1, "vs", ref0, "_", vset)

            # Plots for each GCM: repeat historical
            x1 <- dats_MS1[[dtemp1]][, vars_MS1[[vset]], sc_current]
            rx1 <- array(NA, dim = c(dim(x1), length(gcms)),
              dimnames = c(dimnames(x1), list(gcms)))
            for (k in seq_along(gcms)) {
              rx1[, , k] <- x1
            }
            rx2 <- array(NA, dim = c(dim(x1), length(ens)),
              dimnames = c(dimnames(x1), ens))
            for (k in seq_along(ens)) {
              rx2[, , k] <- x1
            }

            plot_map_vars_dGi_and_modifier_v2(
              x1 = rx1,
              x2 = dats_MS1[[dtemp2]][, vars_MS1[[vset]], gcms, drop = FALSE],
              funs = c("vals", "delta"),
              meta = SFSW2_prj_meta, subset = id_ASA,
              GCMs = reqGCMs_perCS[[sc_name]], exp = sc_name,
              ref0 = ref0, ref1 = ref1,
              mod_with_all = mwa,
              mod_tol = modtols[idt, "tol"],
              expected_dy_lt_dx = worse_direction[[vset]],
              var_labels = labs_MS1[[vset]],
              path = dir_resp_out_tempM, ftag = ftemp)

            # Plots for ensemble across GCMs
            for (ew in ens_wises) {
              plot_map_vars_dGi_and_modifier_v2(
                x1 = rx2,
                x2 = dats_Ens[[ew]][[dtemp2]][, vars_MS1[[vset]], ens, drop = FALSE],
                funs = c("vals", "delta"),
                meta = SFSW2_prj_meta, subset = id_ASA,
                GCMs = ef, exp = sc_name,
                ref0 = ref0, ref1 = ref1,
                mod_with_all = mwa,
                mod_tol = modtols[idt, "tol"],
                expected_dy_lt_dx = worse_direction[[vset]],
                var_labels = labs_MS1[[vset]],
                path = dir_resp_out_tempM,
                ftag = paste0(ftemp, "_", ew, "_", ef))
            }
          }
        }
      }
    }
  }
}


#------------------------------------------------------------------------------#
#--- Barplots and tables of summed area of modifications by geoengineering
#    and of relative change from historical
#------------------------------------------------------------------------------#
if (do_out_mod_bp) {
  ef <- "mean"

  sc_names <- scen_SARs
  ref1u <- unique(deltas_byCS[sc_names])
  temp <- lapply(sc_names, function(x)
    paste(paste0(x, ".", reqGCMs_perCS[[x]]), collapse = "|"))
  gcms <- lapply(temp, function(x) grep(x, sim_scens_id, value = TRUE))
  ens <- paste("hybrid-delta-3mod.d60yrs", sc_names, ef, sep = ".")


  for (vset in c(vset_clim, vset_main)) { #names(vars_MS1)[-1]) {#
    for (ref0 in "Current") {# unique(deltas_byCS)) {
      for (mt in c("dir", "rdi")) {
        dir_resp_out_tempP <- file.path(dir_res_objMS1,
          "Output_ResponseMod_Barplots", paste0(vset, "_", mt))
        dir.create(dir_resp_out_tempP, recursive = TRUE, showWarnings = FALSE)
        dir_resp_out_tempT <- file.path(dir_res_objMS1,
          "Output_ResponseMod_Tabulation", paste0(vset, "_", mt))
        dir.create(dir_resp_out_tempT, recursive = TRUE, showWarnings = FALSE)

        if (mt == "dir") {
          is_ref <- paste(ref1u, collapse = "-")
          use_mod_with_all <- FALSE # c(FALSE, TRUE)
        } else {
          is_ref <- ref0
          use_mod_with_all <- FALSE
        }

        for (mwa in use_mod_with_all) {
          for (idt in seq_len(nrow(modtols))) {
            dtemp2 <- paste0("delta_", is_ref, "_", mt, modtols[idt, "tag"])
            ftemp <- paste0(mt, "_", paste(sc_names, collapse = ""), "to",
              is_ref, "vs", ref0, "_", vset,
              if (mwa) "_GModsAll" else "_GMods2c", modtols[idt, "tag"])

            # Tabulate for each GCM
            tabulate_vars_dGmodifiers(
              xmod = dats_MS1[[dtemp2]][, vars_MS1[[vset]], unlist(gcms), drop = FALSE],
              method = mt, vars = vars_MS1[[vset]], var_labs = labs_MS1[[vset]],
              subset = id_ASA, cell_areas = cells,
              region_data = region_data, name_regions = name_regions,
              GCMs = reqGCMs_perCS[sc_names], exp = sc_names,
              mod_with_all = mwa, mod_tol = modtols[idt, "tol"],
              expected_dy_lt_dx = if (mt == "dir") worse_direction[[vset]],
              path_table = dir_resp_out_tempT, ftag = ftemp)

            # Plots for each GCM
            barplot_vars_dGmodifiers(method = mt,
              vars = vars_MS1[[vset]], var_labs = labs_MS1[[vset]],
              name_regions = name_regions,
              scens = unlist(gcms), GCMs = reqGCMs_perCS[sc_names], exp = sc_names,
              mod_with_all = mwa, mod_tol = modtols[idt, "tol"],
              path_fig = dir_resp_out_tempP, path_table = dir_resp_out_tempT,
              ftag = ftemp)

            # Plots for ensemble across GCMs
            for (ew in ens_wises) {
              tabulate_vars_dGmodifiers(
                xmod = dats_Ens[[ew]][[dtemp2]][, vars_MS1[[vset]], ens, drop = FALSE],
                method = mt, vars = vars_MS1[[vset]], var_labs = labs_MS1[[vset]],
                subset = id_ASA, cell_areas = cells,
                region_data = region_data, name_regions = name_regions,
                GCMs = rep(list(ef), length(sc_names)), exp = sc_names,
                mod_with_all = mwa,
                mod_tol = modtols[idt, "tol"],
                expected_dy_lt_dx = if (mt == "dir") worse_direction[[vset]],
                path_table = dir_resp_out_tempT,
                ftag = paste0(ftemp, "_", ew, "_", ef))

              barplot_vars_dGmodifiers(method = mt,
                vars = vars_MS1[[vset]], var_labs = labs_MS1[[vset]],
                name_regions = name_regions,
                scens = ens, GCMs = rep(list(ef), length(sc_names)), exp = sc_names,
                mod_with_all = mwa, mod_tol = modtols[idt, "tol"],
                path_fig = dir_resp_out_tempP, path_table = dir_resp_out_tempT,
                ftag = paste0(ftemp, "_", ew, "_", ef))

              # Plots for all GCMs on top of ensemble across GCMs
              if (!mwa) {
                barplot_vars_dGmodifiers_v3(method = mt,
                  vars = vars_MS1[[vset]], var_labs = labs_MS1[[vset]],
                  name_regions = name_regions,
                  scens = unlist(gcms), scens_ens = ens,
                  GCMs = reqGCMs_perCS[sc_names],
                  GCMs_ens = rep(list(ef), length(sc_names)),
                  GCM_target = "BNU-ESM", ref = is_ref,
                  recalc_ens = TRUE,
                  mod_with_all = mwa, mod_tol = modtols[idt, "tol"],
                  path_fig = dir_resp_out_tempP, path_table = dir_resp_out_tempT,
                  ftag1 = ftemp, ftag2 = paste0(ew, "_", ef),
                  ftag3 = "AllGCMs_v4")
              }
            }
          }
        }
      }
    }
  }
}


#------------------------------------------------------------------------------#
#--- MAP RELATIVE AREA vs. DIRECTIONAL CHANGE (similar to CDFs)

if (do_out_CAFs) {
  # Global
  for (vset in c(vset_clim, vset_main)) { #names(vars_MS1)[-1]) {#
    for (sc in seq_along(scen_SARs)) {
      dir_resp_out_temp <- file.path(dir_res_objMS1, "Output_Response_CAFs",
        paste0(vset, "_", scen_SARs[sc]))
      dir.create(dir_resp_out_temp, recursive = TRUE, showWarnings = FALSE)

      plot_dCDF2s_regional(
        data_dvals = dats_MS1[["delta_RCP45_abs"]],
        data_dir = dats_MS1[["delta_RCP45_dir_tol0p500"]],
        vars = vars_MS1[[vset]], var_labs = labs_MS1[[vset]],
        SAR = scen_SARs[sc], sc_ref = "RCP45", sc_base = "Current",
        cat_dirs = define_direction_2deltas(),
        cell_within_study_TF = id_ASA, cell_areas = cells,
        region_data = region_data, name_regions = "Global",
        panel_indep = TRUE, path = dir_resp_out_temp,
        ftag = paste0("Global_", scen_SARs[sc], "_", vset))
    }
  }

  # By regions
  vset <- "res2"

  for (iv in seq_along(vars_MS1[[vset]])) {
    for (sc in seq_along(scen_SARs)) {
      dir_resp_out_temp <- file.path(dir_res_objMS1, "Output_Response_CAFs",
        paste0(vset, "_", scen_SARs[sc]))
      dir.create(dir_resp_out_temp, recursive = TRUE, showWarnings = FALSE)

      plot_dCDF2s_regional(
        data_dvals = dats_MS1[["delta_RCP45_abs"]],
        data_dir = dats_MS1[["delta_RCP45_dir_tol0p500"]],
        vars = vars_MS1[[vset]][iv], var_labs = labs_MS1[[vset]][iv],
        SAR = scen_SARs[sc], sc_ref = "RCP45", sc_base = "Current",
        cat_dirs = define_direction_2deltas(),
        cell_within_study_TF = id_ASA, cell_areas = cells,
        region_data = region_data,
        name_regions = grep("Global", name_regions, invert = TRUE,
          value = TRUE),
        panel_indep = TRUE, path = dir_resp_out_temp,
        ftag = paste0("Regional_", scen_SARs[sc], vars_MS1[[vset]][iv]))
    }
  }
}


#------------------------------------------------------------------------------#
#--- CALCULATE REGRESSION BETWEEN CHANGE IN RESPONSE AND CHANGE IN CLIMATE

if (do_out_variancepartition) {
  ef <- "mean"

  for (vset in c(vset_clim, vset_main)) { #names(vars_MS1)[-1]) {#
    dir_resp_out_temp <- file.path(dir_res_objMS1, "Output_VarPartitioning",
      vset)
    dir.create(dir_resp_out_temp, recursive = TRUE, showWarnings = FALSE)

    regress_predictor_contributions(
      data = dats_MS1,
      response = vars_MS1[[vset]], resp_labs = labs_MS1[[vset]],
      predictors = vars_MS1[[vset_clim]],
      ref_condition = "Current",
      cc_condition = "RCP45",
      SAR_conditions = scen_SARs,
      cell_within_study_TF = id_ASA, cell_areas = cells,
      path = dir_resp_out_temp, ftag = paste0(vset, "_data"))

    for (ew in ens_wises) {
      regress_predictor_contributions(
        data = dats_Ens[[ew]],
        response = vars_MS1[[vset]], resp_labs = labs_MS1[[vset]],
        predictors = vars_MS1[[vset_clim]],
        ref_condition = paste("Current", ef, sep = "."),
        cc_condition = paste("RCP45", ef, sep = "."),
        SAR_conditions = paste(scen_SARs, ef, sep = "."),
        cell_within_study_TF = id_ASA, cell_areas = cells,
        path = dir_resp_out_temp, ftag = paste0(vset, "_", ew, "-", ef))
    }
  }
}


#------------------------------------------------------------------------------#
#--- Tabulate SPATIAL STRUCTURE IN Gi modification categories


if (do_out_spatialpattern) {
  dir_resp_out_temp <- file.path(dir_res_objMS1, "Output_SpatialPatterns")
  dir.create(dir_resp_out_temp, recursive = TRUE, showWarnings = FALSE)

  ef <- "mean"
  ir_global <- "Global"

  ivars <- match(unlist(vars_MS1[c(vset_clim, vset_main)]), vars_all,
    nomatch = 0)

  Ns <- lengths(reqGCMs_perCS[scen_SARs])
  ids_disagree50 <- lapply(Ns, function(N) seq_len(N) / N <= 1 / 2)

  mwa <- TRUE
  for (idt in seq_len(nrow(modtols))) {
    # Spatial (dis)agreement
    ref1u <- unique(deltas_byCS[scen_SARs])
    dtemp2 <- paste0("delta_", ref1u, "_dir", modtols[idt, "tag"])
    sa_temp <- dats_ModAgree[[dtemp2]]["agreement", , vars_all[ivars], scen_SARs]

    x_sda <- sapply(scen_SARs, function(sc) {
        apply(sa_temp[, , sc], 2, function(x) {
          temp <- tapply(cells[, "km2"], INDEX = x, sum, na.rm = TRUE)
          temp <- temp / sum(temp)
          sum(temp[names(temp) %in% which(ids_disagree50[[sc]])])
        })
      })


    # Spatial structure: arrange data
    ss_temp <- NULL
    for (k in seq_along(scen_SARs)) {
      sc_name <- scen_SARs[k]

      for (kdata in kdatas) {
        if (kdata == kdatas["val"]) {
          used_ens_wises <- "data"
          scen_avail <- sim_scens_id
        } else if (kdata == kdatas["ens"]) {
          used_ens_wises <- ens_wises
          scen_avail <- grep(paste(ens_funCT, collapse = "|"),
            dimnames(dats_Ens[[1]][[1]])[[3]], value = TRUE)
        } else {
          stop("Unkown ", shQuote(kdata))
        }

        # Spatial structure
        for (ew in used_ens_wises) {
          temp <- res_SpatPat[ew, sc_name, modtols[idt, "tag"]][[1]]
          if (kdata == kdatas["ens"]) {
            id_ef <- grepl(paste0(".", ef, "$"), temp[, "Scenario"])
            temp <- temp[id_ef, ]
          }

          ss_temp <- rbind(ss_temp, temp,
            make.row.names = FALSE, stringsAsFactors = FALSE)
        }
      }
    }

    id_region <- ss_temp[, "Region"] == ir_global
    x_ss <- reshape2::dcast(ss_temp[id_region, ], Variable ~ Scenario,
      value.var = "Aggregation")

    # Save table to file
    ivars_sda <- match(vars_all[ivars], rownames(x_sda))
    ivars_ss <- match(vars_all[ivars], x_ss[, "Variable"])
    temp <- t(data.frame(
      Variable = labs_all[ivars],
      Disagreement50 = x_sda[ivars_sda, ],
      Spatial_Aggregation = x_ss[ivars_ss, -1],
      row.names = NULL))

    ftemp <- paste0(paste(scen_SARs, collapse = ""), "to",
      paste(ref1u, collapse = "-"), "_", ef, "_", ir_global,
      if (mwa) "_GModsAll" else "_GMods2c", modtols[idt, "tag"])

    write.csv(temp, file = file.path(dir_resp_out_temp,
      paste0("Table_SpatialPattern_", ftemp, ".csv")))


    # Plot: spatial structure
    if (FALSE) {
      par_prev <- par(mar = c(10, 10, 1, 1))
      image(as.matrix(x_ss[, -1]), col = blues9, axes = FALSE)
      axis(side = 1, at = seq(0, 1, length.out = 5),
        labels = x_ss[, 1], las = 2)
      axis(side = 2, at = seq(0, 1, length.out = 12),
        labels = colnames(x_ss[, -1]), las = 2)
      par(par_prev)
    }
  }
}


#------------------------------------------------------------------------------#

#--- Illustrate modification categories
if (do_out_mod_illustration) {
  dir_back <- file.path(dir_res_objMS1, "Output_Background")
  dir.create(dir_back, recursive = TRUE, showWarnings = FALSE)

  cols1 <- c("gray",
    sd = "lightblue", ri = "purple", li = "darkred",
    ld = "darkblue", rd = "mediumseagreen", si = "orange")
  cols2 <- c("gray", "purple", "aquamarine2")

  cl <- define_direction_2deltas(TRUE)
  cats2 <- c("More severe", "Less severe")

  n <- 100
  dx <- matrix(rnorm(n, sd = 2), ncol = 1, dimnames = list(NULL, "var"))
  dy <- matrix(dx + rnorm(n), ncol = 1, dimnames = list(NULL, "var"))

  n_panels <- c(2, nrow(modtols))

  png(units = "in", res = 150,
    height = 5 * n_panels[1], width = 7 * n_panels[1],
    file = file.path(dir_back, "Fig_IllustrateModifications.png"))
  par_prev <- par(mfcol = n_panels, mar = c(2.5, 2.5, 1, 1), mgp = c(1, 0, 0),
    tcl = 0.3, cex = 1)

  for (idt in seq_len(nrow(modtols))) {
    xy1 <- compare_direction_2deltas(dx, dy, cl = cl,
      tol = modtols[idt, "tol"])[, 1]
    fxy1 <- 1 + as.integer(xy1)

    xy2 <- factor(xy1, levels = cl[, "name"],
      labels = cats2[1 + as.integer(cl[, "expected"])])
    fxy2 <- 1 + as.integer(xy2)

    plot(dx[, 1], dy[, 1], asp = 1,
      xlab = "d(RCP45, historic)", ylab = "d(G3, historic)")
    abline(h = 0)
    abline(v = 0)
    abline(a = 0, b = 1)
    points(dx, dy, pch = 16, col = cols1[fxy1])
    mtext(side = 3, adj = 0, text =
        paste0("7 cat. | change > ", 100 * modtols[idt, "tol"], "%"))
    if (idt == 1) {
      legend("bottomright", legend = c("No/small modification", levels(xy1)),
        pch = 16, col = cols1, pt.cex = 2)
    }

    plot(dx[, 1], dy[, 1], asp = 1,
      xlab = "d(RCP45, historic)", ylab = "d(G3, historic)")
    abline(h = 0)
    abline(v = 0)
    abline(a = 0, b = 1)
    points(dx, dy, pch = 16, col = cols2[fxy2])
    mtext(side = 3, adj = 0, text =
        paste0("3 cat. | change > ", 100 * modtols[idt, "tol"], "%"))
    if (idt == 1) {
      legend("bottomright", legend = c("No/small modification", levels(xy2)),
        pch = 16, col = cols2, pt.cex = 2)
    }
  }

  par(par_prev)
  dev.off()
}


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
