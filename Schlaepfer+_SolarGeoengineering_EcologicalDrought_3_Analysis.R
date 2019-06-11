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


if (!exists("SGEcologicalDrought_3_Analysis")) {


#--- Load packages, custom functions, description of simulation experiment,
# and prepared data
source("Schlaepfer+_SolarGeoengineering_EcologicalDrought_2_PrepareData.R")

stopifnot(requireNamespace("reshape2"))


#------------------------------------------------------------------------------#
#--------- ANALYSIS (per GCM)
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#--- CALCULATE CHANGE IN RESPONSE VARIABLES (part 1)
print(paste(Sys.time(), "--", "Calculate absolute and relative change"))

# Absolute change
dats_MS1 <- c(dats_MS1,
  calc_response_change_from_reference(data = dats_MS1[["vals"]],
    ref_condition = c(sc_current, "RCP45"), method = "absolute"))

# Relative change
dats_MS1 <- c(dats_MS1,
  calc_response_change_from_reference(data = dats_MS1[["vals"]],
    ref_condition = c(sc_current, "RCP45"), method = "relative"))


# Make sure that we didn't introduce NAs
hasNAs <- sapply(dats_MS1, function(x) {
    temp <- grep(sc_current, dimnames(x)[[3]])
    anyNA(x[id_ASA, , -temp])
  })

if (any(hasNAs)) {
  stop("We have NAs in: ", paste(shQuote(names(hasNAs)[hasNAs]),
    collapse = ", "))
}



#------------------------------------------------------------------------------#
#--------- ENSEMBLES (per experimental family)
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#--- CALCULATE ENSEMBLE: (part 1)
print(paste(Sys.time(), "--", "Calculate ensemble values (part 1)"))

ens_funCT <- c("mean", "median")
ens_funs <- c("min", ens_funCT, "max") # other options: "span", "agreement"
ens_wises <- c("EnsCW", "EnsRW", "EnsGW")

rtemp <- region_data
rtemp[!id_ASA] <- NA

dats_Ens <- calculate_ensembles(SFSW2_prj_meta,
  data = dats_MS1,
  data_names = grep("dir|rdi", names(dats_MS1), invert = TRUE, value = TRUE),
  subset = id_ASA, cell_area_km2 = cells[, "km2"], id_region = rtemp,
  sc_historical = sc_current, req_Downs = req_Downs, req_dTime = req_dTime,
  ens_wises = ens_wises,
  ens_funs = ens_funs,
  path = file.path(dir_res_data, "temp_Ens"),
  ftag = ftag_MS1
)


#------------------------------------------------------------------------------#
#--- DETERMINE CONSERVATIVE/IDEAL GCM
print(paste(Sys.time(), "--", "Identify GCMs in ensemble"))

ens_outs <- data.frame(
  id = c(rep(1, length(vars_MS1[[vset_clim]])),
          rep(2, length(vars_MS1[[vset_main]]))),
  vars = c(vars_MS1[[vset_clim]], vars_MS1[[vset_main]]),
  labs = c(labs_MS1[[vset_clim]], labs_MS1[[vset_main]]),
  worse = c(worse_direction[[vset_clim]], worse_direction[[vset_main]]),
  stringsAsFactors = FALSE)

ens_methods <- c("absolute", "relative")
ens_types <- c("GCM_minGChange", "GCM_minGWorse", "GCM_maxGRevert")
labs_ens_types <- c("min{dG}", "min{dG:worse}", "max{dG:revert}")

reqCSs_ens1 <- reqMs_ens1 <- paste0(scen_SARs, ".", rep(ens_funCT, each = 2))

# Output containers
dats_Ens2 <- vector(mode = "list", length = length(ens_wises))
names(dats_Ens2) <- ens_wises
igcm_ens <- dats_Ens2

vals_Ens2 <- ranks_Ens2 <- matrix(list(),
  nrow = length(ens_methods), ncol = length(ens_wises),
  dimnames = list(ens_methods, ens_wises))

valsCT_Ens2 <- matrix(list(),
  nrow = length(ens_methods[2]) * length(ens_funCT), ncol = length(ens_wises),
  dimnames = list(paste0(ens_methods[2], "_", ens_funCT), ens_wises))


# Calculate ensembles
de2_template <- matrix(list(),
  nrow = length(ens_methods), ncol = length(ens_types),
  dimnames = list(ens_methods, ens_types))

for (ew in ens_wises) {
  dats_Ens2[[ew]] <- de2_template

  for (m in seq_along(ens_methods)) {
    print(paste(Sys.time(), dQuote(ew), "ensembles for", ens_methods[m]))

    mid <- substr(ens_methods[m], 1, 3)

    for (k in seq_along(ens_types)) {
      if (ens_types[k] == "GCM_minGChange") {
        # (i) GCM that minimizes how geoengineering modifies climate change:
        #     d(Gi, current) == d(RCP45, current), i.e., min of abs(d(Gi, RCP45))
        dtemp1 <- abs(dats_MS1[[paste0("delta_RCP45_", mid)]])
        dtemp2 <- abs(dats_Ens[[ew]][[paste0("delta_RCP45_", mid)]])

      } else if (ens_types[k] == "GCM_minGWorse") {
        # (ii) GCM that minimizes how geoengineering worsens climate change:
        #     min of [fworse * (d(Gi, RCP45)) > 0]
        dtemp1 <- dats_MS1[[paste0("delta_RCP45_", mid)]]
        dtemp2 <- dats_Ens[[ew]][[paste0("delta_RCP45_", mid)]]
        for (i in seq_along(ens_outs[, "vars"])) {
          var <- ens_outs[i, "vars"]
          dtemp1[, var, ] <- ens_outs[i, "worse"] * dtemp1[, var, ]
          dtemp2[, var, ] <- ens_outs[i, "worse"] * dtemp2[, var, ]
        }
        dtemp1[!is.na(dtemp1) & dtemp1 < 0] <- 0
        dtemp2[!is.na(dtemp2) & dtemp2 < 0] <- 0

      } else if (ens_types[k] == "GCM_maxGRevert") {
        # (iii) GCM that maximizes how geoengineering reverts climate change back
        #       to current: d(Gi, current) == 0, i.e., min of abs(d(Gi, current))
        dtemp1 <- abs(dats_MS1[[paste0("delta_Current_", mid)]])
        dtemp2 <- abs(dats_Ens[[ew]][[paste0("delta_Current_", mid)]])
      }

      # Remove infinite values found in relative data
      temp <- is.infinite(dtemp1)
      if (any(temp)) {
        dtemp1[temp] <- NA
      }
      temp <- is.infinite(dtemp2)
      if (any(temp)) {
        dtemp2[temp] <- NA
      }

      temp1 <- calc_extentwise_ensemble(SFSW2_prj_meta,
        data = dtemp1[, ens_outs[, "vars"], ], area = cells[, "km2"],
        fcentral = "mean", funs = NULL, probs = NULL, ties.method = "average",
        na.rm = TRUE)
      temp2 <- calc_extentwise_ensemble(SFSW2_prj_meta,
        data = dtemp2[, ens_outs[, "vars"], ], area = cells[, "km2"],
        fcentral = "mean", funs = NULL, probs = NULL, ties.method = "average",
        reqCSs = reqCSs_ens1, reqMs = reqMs_ens1, na.rm = TRUE)

      temp0a <- reshape2::melt(
        temp1[["ensemble_structure"]][, , , scen_SARs])
      temp0a[, "Var3"] <- as.character(temp0a[, "Var3"])
      temp0a[, "Var4"] <- as.character(temp0a[, "Var4"])

      temp0b <- NULL
      for (k2 in seq_along(reqCSs_ens1)) {
        ttemp <- reshape2::melt(
          temp2[["ensemble_structure"]][, , reqMs_ens1[k2], reqCSs_ens1[k2], drop = FALSE])
        ttemp[, "Var3"] <- as.character(ttemp[, "Var3"])
        ttemp[, "Var4"] <- substr(as.character(ttemp[, "Var4"]), 1, 2)
        temp0b <- rbind(temp0b, ttemp)
      }

      dats_Ens2[[ew]][ens_methods[m], ens_types[k]][[1]] <-
        reshape2::acast(rbind(temp0a, temp0b), Var1 ~ Var2 ~ Var3 ~ Var4)
    }
  }
}

# Put ensemble values into table
for (ew in ens_wises) {
  for (m in seq_along(ens_methods)) {
    for (k in seq_along(ens_types)) {
      temp <- reshape2::melt(
        dats_Ens2[[ew]][ens_methods[m], ens_types[k]][[1]]["fcentral", , , ])
      vals_Ens2[m, ew][[1]] <- rbind(vals_Ens2[m, ew][[1]],
        data.frame(ens_types = ens_types[k],
          reshape2::dcast(temp, Var3 + Var1 ~ Var2)))
    }
    vals_Ens2[m, ew][[1]][, "Var1"] <- ens_outs[, "labs"]
  }
}

# Aggregate relative values across variables (doesn't make sense for absolute)
em <- ens_methods[2]
for (ew in ens_wises) {
  temp <- vals_Ens2[em, ew][[1]]

  for (m in ens_funCT) {
    valsCT_Ens2[paste0(em, "_", m), ew][[1]] <-
      aggregate(temp[, -(1:3)], by = temp[, 1:2], m)
  }
}

# Aggregate ranks across variables
re2_template <- data.frame(
  out <- expand.grid(ens_types = ens_types, SAR = scen_SARs),
  matrix(NA, nrow = nrow(out), ncol = length(reqMs),
    dimnames = list(NULL, reqMs)))

for (ew in ens_wises) {
  for (m in seq_along(ens_methods)) {
    ranks_Ens2[m, ew][[1]] <- re2_template

    for (k in seq_along(ens_types)) {
      temp <- reshape2::melt(
        dats_Ens2[[ew]][ens_methods[m], ens_types[k]][[1]]["Rank", , reqMs, scen_SARs])
      temp <- reshape2::dcast(temp, Var3 + Var1 ~ Var2)
      temp <- aggregate(temp[, -c(1:2)], by = list(temp[, "Var3"]), median)

      ids <- ranks_Ens2[m, ew][[1]][, "ens_types"] == ens_types[k]
      ranks_Ens2[m, ew][[1]][ids, -c(1:2)] <- temp[, -1]
    }
  }
}

# Extract ensemble GCM names
for (ew in ens_wises) {
  ids <- unique(c(1, 2,
    which(!(colnames(valsCT_Ens2["relative_mean", ew][[1]]) %in% reqMs_ens1))))

  temp <- c(
    absolute = ranks_Ens2["absolute", ew],
    relative = ranks_Ens2["relative", ew],
    relative_mean = list(valsCT_Ens2["relative_mean", ew][[1]][, ids]),
    relative_median = list(valsCT_Ens2["relative_median", ew][[1]][, ids])
  )

  igcm_ens[[ew]] <- lapply(temp, function(x)
    data.frame(x[, c(1:2)], GCM = reqMs[apply(x[, -c(1:2)], 1, which.min)],
      stringsAsFactors = FALSE))
}


#------------------------------------------------------------------------------#
#--- CALCULATE CHANGE IN RESPONSE VARIABLES (part 2)
print(paste(Sys.time(), "--", "Calculate modification by geoengineering"))

modtols <- data.frame(tol = c(NA, 0.05, 0.2, 0.5))
modtols[, "tag"] <- paste0("_tol", sub(getOption("OutDec"), "p",
  formatC(ifelse(is.na(modtols[, "tol"]), 0, modtols[, "tol"]), format = "f",
    flag = "0", digits = 3), fixed = TRUE))

# Modification categories: direction of any/at least X% change by geoengineering
# (Gi) relative to climate change as is (RCP45) and relative to historical
dats_ModAgree <- list()
ref_conditions <- c(sc_current, "RCP45")
dnames_dir <- list()

for (ref_cond in ref_conditions) {
  method <- if (ref_cond == sc_current) "rdirection" else "direction"

  for (idt in seq_len(nrow(modtols))) {
    for (kdata in kdatas) {

      use_ens_wises <- if (kdata == kdatas["val"]) "data" else ens_wises

      for (ew in use_ens_wises) {
        fname_GtoScCond <- file.path(dir_res_data, paste0(ftag_MS1,
          "_DirectionChange_Gto", ref_cond, modtols[idt, "tag"], "_",
          ew, ".rds"))
        tag <- paste0("delta_", ref_cond, "_", substr(method, 1, 3),
          modtols[idt, "tag"])


        # Extract or calculate directional change
        if (file.exists(fname_GtoScCond)) {
          xdata <- readRDS(fname_GtoScCond)

        } else {
          use_dname <- if (method == "rdirection") "vals" else "delta_Current_abs"
          temp <- if (kdata == kdatas["val"]) {
              dats_MS1[[use_dname]]
            } else {
              dats_Ens[[ew]][[use_dname]]
            }

          xdata <- calc_response_change_from_reference(data = temp,
            ref_condition = ref_cond, method = method,
            tol = modtols[idt, "tol"], subset = id_ASA)

          saveRDS(xdata, file = fname_GtoScCond)
        }


        # Check that well formed: same number of dimensions; same number of
        #   sites and variables; same variables
        tt <- list(dats_MS1[["vals"]], xdata[[tag]])
        ds <- lapply(tt, function(x) list(d = dim(x), dn = dimnames(x)))
        rm(tt)
        stopifnot(
          identical(length(ds[[1]][["d"]]), length(ds[[2]][["d"]])),
          identical(ds[[1]][["d"]][1:2], ds[[2]][["d"]][1:2]),
          identical(ds[[1]][["dn"]][2], ds[[2]][["dn"]][2]))

        # Make sure that we didn't introduce NAs
        hasNAs <- sapply(xdata, function(x) {
          temp <- grep(sc_current, dimnames(x)[[3]])
          anyNA(x[id_ASA, , -temp])
        })

        if (any(hasNAs)) {
          stop("We have NAs in: ", paste(shQuote(names(hasNAs)[hasNAs]),
            collapse = ", "))
        }



        # Add directional change to data list and calculate agreement if needed
        if (kdata == kdatas["val"]) {
          dats_MS1 <- c(dats_MS1, xdata)

          # Calculate agreement among GCMs in modification categories
          if (method == "direction") {
            for (mwa in FALSE) { #c(FALSE, TRUE)) {#
              fname_GtoScCond_agree <- file.path(dir_res_data, paste0(ftag_MS1,
                "_DirectionChangeAgreement_Gto", ref_cond,
                if (mwa) "_GModsAll" else "_GMods2c", modtols[idt, "tag"],
                "_", ew, ".rds"))

              if (file.exists(fname_GtoScCond_agree)) {
                modagree <- readRDS(fname_GtoScCond_agree)

              } else {
                #--- Determine agreement in modification across GCMs
                dn <- dimnames(xdata[[1]])
                vars <- dn[[2]]
                scens <- dn[[3]]

                # Remap
                for (iv in seq_along(vars)) {
                  if (mwa) {
                    xdata[[1]][, iv, ] <- as.character(
                      factor_6directions(x[, iv, ]))
                  } else {
                    vdir <- as.logical(1 +
                        worse_direction_all[which(vars[iv] == vars_all)])

                    xdata[[1]][, iv, ] <- as.character(remap_directions_6to2(
                      as.vector(xdata[[1]][, iv, ]), expected_dy_lt_dx = vdir))
                  }
                }

                # Calculate agreement among GCMs in modification categories
                modagree <- calc_cellwise_ensemble(SFSW2_prj_meta,
                  data = xdata[[1]], funs = "agreement", na.rm = TRUE,
                  add_historic = FALSE, reqCSs = scen_SARs, verbose = TRUE)

                saveRDS(modagree, file = fname_GtoScCond_agree)
              }

              dats_ModAgree[[names(xdata)]] <- modagree
            }
          }

        } else {
          dats_Ens[[ew]] <- c(dats_Ens[[ew]], xdata)
        }

        rm(xdata)
      }
    }
  }

  temp <- paste0(ref_cond, "_", substr(method, 1, 3))
  dnames_dir[[ref_cond]] <- grep(temp, names(dats_MS1), value = TRUE)
}




#------------------------------------------------------------------------------#
#--- CALCULATE ENSEMBLE: (part 2)
print(paste(Sys.time(), "--", "Calculate ensemble values (part 2)"))

rtemp <- region_data
rtemp[!id_ASA] <- NA

#TODO: fix
dats_Ens3 <- calculate_ensembles(SFSW2_prj_meta,
  data = dats_MS1,
  data_names = grep("dir", names(dats_MS1), invert = FALSE, value = TRUE),
  subset = id_ASA, cell_area_km2 = cells[, "km2"], id_region = rtemp,
  sc_historical = sc_current, req_Downs = req_Downs, req_dTime = req_dTime,
  ens_wises = "EnsCW",
  ens_funs = "majority",
  path = file.path(dir_res_data, "temp_Ens"),
  ftag = ftag_MS1
)




#------------------------------------------------------------------------------#
#--- CALCULATE SPATIAL STRUCTURE IN Gi modification categories
print(paste(Sys.time(), "--", "Calculate spatial structure in modification"))


res_SpatPat <- array(list(),
  dim = c(1 + length(ens_wises), length(scen_SARs), nrow(modtols)),
  dimnames = list(c("data", ens_wises), scen_SARs, modtols[, "tag"]))


ivars <- match(unlist(vars_MS1[c(vset_clim, vset_main)]), vars_all,
  nomatch = 0)


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

  for (ew in used_ens_wises) {
    for (k in seq_along(scen_SARs)) {
      sc_name <- scen_SARs[k]
      ref1 <- deltas_byCS[sc_name]

      scens <- grep(paste(sc_name, collapse = "|"), scen_avail, value = TRUE)

      mwa <- FALSE
      for (idt in seq_len(nrow(modtols))) {
        dtemp2 <- paste0("delta_", ref1, "_dir", modtols[idt, "tag"])

        ftemp_SpatPat <- file.path(dir_res_data, paste0(ftag_MS1,
          "_SpatialPatterns-", kdata, "-", ew, "_", sc_name, "_",
          if (mwa) "GModsAll" else "GMods2c", "_", dtemp2, ".rds"))

        if (file.exists(ftemp_SpatPat)) {
          temp <- readRDS(ftemp_SpatPat)

        } else {
          # Calculate spatial patterns
          xtemp <- if (kdata == kdatas["val"]) {
              dats_MS1[[dtemp2]][, vars_all[ivars], scens, drop = FALSE]
            } else {
              dats_Ens[[ew]][[dtemp2]][, vars_all[ivars], scens, drop = FALSE]
            }

          temp <- range_spatialpattern(xtemp,
            meta = SFSW2_prj_meta,
            cell_within_study_TF = id_ASA,
            region_data = region_data, name_regions = name_regions,
            mod_with_all = mwa, expected_dy_lt_dx = worse_direction_all[ivars],
            path_temp = file.path(dir_res_data, "tmp_SpatialPatterns"),
            responses = c("Aggregation"), ftag = ftag_MS1, verbose = TRUE)

          saveRDS(temp, file = ftemp_SpatPat)
        }

        res_SpatPat[ew, sc_name, modtols[idt, "tag"]][[1]] <- temp
      }
    }
  }
}


} # endif not exists "SGEcologicalDrought_3_Analysis"

SGEcologicalDrought_3_Analysis <- TRUE


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
