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


#------------------------------------------------------------------------------#
#--------- CUSTOM FUNCTIONS
#------------------------------------------------------------------------------#

if (!exists("SGEcologicalDrought_0_Functions")) {

pkgs <- c("digest", "raster", "sp", "ncdf4", "geosphere", "Hmisc", "maptools",
  "maps", "reshape2", "viridisLite", "landscapemetrics")
stopifnot(all(sapply(pkgs, requireNamespace)))


#' Extract data from \var{dbOutput}
#'
#' This function is internally "memoised" per variable, i.e., it loads an
#' \code{array} from \code{file} if previously extracted from
#' \var{dbOutput}; otherwise, it extracts first data (and then saves to a file
#' disk).
#'
#' @param file A character string. The file path to the \var{rds}-file(s)
#'   containing the memoised data. If \code{NULL}, then data is not written/read
#'   from file, but extracted from \var{dbOutput}.
#'
#' @seealso \code{\link{get.SeveralOverallVariables_Scenario}},
#'   \code{\link{saveRDS}}
#' @return A three-dimensional numerical \code{array} with dimensions
#'   \code{runsN_sites}, \code{variables}, and \code{scenarios}.
#'
#' @export
extract_dbOut_to_df <- function(SFSW2_prj_meta, variables,
  MeanOrSD = "Mean", whereClause = NULL, runsN_sites = NULL, scenarios = NULL,
  fname_dbOut = NULL, file = NULL, verbose = FALSE) {

  if (is.null(runsN_sites)) {
    if (!is.null(SFSW2_prj_meta)) {
      runsN_sites <- SFSW2_prj_meta[["sim_size"]][["runsN_sites"]]
    } else {
      stop("Argument 'runsN_sites' is missing.")
    }
  }

  if (is.null(scenarios)) {
    if (!is.null(SFSW2_prj_meta)) {
      scenarios <- SFSW2_prj_meta[["sim_scens"]][["id"]]
    } else {
      stop("Argument 'scenarios' is missing.")
    }
  }
  scensN <- length(scenarios)

  if (is.null(fname_dbOut)) {
    if (!is.null(SFSW2_prj_meta)) {
      fname_dbOut <- SFSW2_prj_meta[["fnames_out"]][["dbOutput"]]
    } else {
      stop("Argument 'fname_dbOut' is missing.")
    }
  }

  if (!is.null(file)) {
    temp <- strsplit(fname_dbOut, split = .Platform$file.sep)[[1]]
    n <- length(temp)
    temp <- paste(temp[max(1, n - 2):n], collapse = .Platform$file.sep)
    temp <- list(fdbrSFSW2 = temp, MeanOrSD = MeanOrSD,
      whereClause = whereClause, runsN_sites = runsN_sites)
    call_id <- digest::digest(temp, algo = "sha1")

    dir_file <- dirname(file)
    dir.create(dir_file, recursive = TRUE, showWarnings = FALSE)
    base_file <- basename(file)
    ext_file <- ".rds"

    if (isTRUE(grepl(paste0(ext_file, "$"), base_file))) {
      temp <- substr(base_file, 1, nchar(base_file) - nchar(ext_file))

    } else {
      temp <- NULL
      dir_file <- file
    }

    tag_file <- paste("dbOut", temp, call_id, sep = "_")

    get_ftemp <- function(dir = dir_file, tag = tag_file, ext = ext_file,
      var, scen) {

      temp <- list(variable = as.character(var), scenario = as.character(scen))
      id <- digest::digest(temp, algo = "sha1")
      file.path(dir, paste0(tag, "-", id, ext))
    }

  } else {
    ftemp <- NA_character_
  }


  data <- array(NA, dim = c(runsN_sites, length(variables), scensN),
    dimnames = list(NULL, variables, scenarios))

  for (sc in seq_len(scensN)) {
    if (verbose) {
      print(paste(Sys.time(), scenarios[sc], whereClause))
    }

    # Load variables that have already been extracted for this scenario
    # and compile list of variables that still need to be extracted from dbOut
    vars_for_this_sc <- NULL

    for (iv in seq_along(variables)) {
      if (!is.null(file)) {
        ftemp <- get_ftemp(var = variables[iv], scen = scenarios[sc])
      }

      if (file.exists(ftemp)) {
        data[, variables[iv], scenarios[sc]] <- readRDS(ftemp)

      } else {
        vars_for_this_sc <- c(vars_for_this_sc, variables[iv])
      }
    }

    # Extract variables from dbOut and store in rds file
    if (length(vars_for_this_sc) > 0) {
      datv <- rSFSW2::get.SeveralOverallVariables_Scenario(
        fdbrSFSW2 = fname_dbOut, responseName = vars_for_this_sc,
        MeanOrSD = MeanOrSD, scenario = scenarios[sc],
        whereClause = whereClause)
      datv <- as.matrix(datv)

      if (identical(dim(datv), c(runsN_sites, length(vars_for_this_sc)))) {
        data[, vars_for_this_sc, scenarios[sc]] <- datv

        if (!is.null(file)) for (iv in seq_along(vars_for_this_sc)) {
          ftemp <- get_ftemp(var = vars_for_this_sc[iv], scen = scenarios[sc])
          dir.create(dirname(ftemp), recursive = TRUE, showWarnings = FALSE)
          saveRDS(data[, vars_for_this_sc[iv], scenarios[sc]], file = ftemp)
        }

      } else {
        stop("Dimension mismatch:",
          " requested = ", paste(dim(data[, , sc]), collapse = "/"),
          " extracted = ", paste(dim(datv), collapse = "/"),
          " for requested variables: ", paste(vars_for_this_sc, collapse = "/"),
          " and extracted variables: ", paste(dimnames(datv)[[2L]],
            collapse = "/"))
      }
    }
  }

  data
}



#' Calculate absolute or relative change in response variables
#'
#' Change in response variables between a climate reference
#' condition (or scenario family) and each climate condition is calculated as
#' \itemize{
#'   \item \code{x[condition k] - x[reference]}
#'          if \code{method} is \code{"absolute"}
#'   \item \code{(x[condition k] - x[reference]) / x[reference]}
#'          if \code{method} is \code{"relative"}
#'   \item the result of \code{\link{compare_direction_2deltas}}
#'          if \code{method} is \code{"direction"}
#'   \item the result of \code{(x[condition k] > tol, x[condition k] < tol)}
#'          if \code{method} is relative direction, i.e., \code{"rdirection"}
#' }
#'
#' If \code{reference} is the same for
#' each scenario, e.g., "current", then the same value is used for each
#' \code{condition k}; if \code{reference} is a scenario family, e.g., "RCP45",
#' then the appropriate \code{x} values per GCM are used.
#'
#' @param data A three-dimensional numerical \code{array} with dimensions
#'   \code{runsN_sites}, \code{variables}, and \code{scenarios}, e.g., the
#'   returned object from \code{\link{extract_dbOut_to_df}}.
#' @param ref_condition A vector of character strings. Each element is used
#'   as reference to calculate a separate set of change values.
#' @param method A character string. See description.
#' @param ... Additional arguments, e.g., \describe{
#'   \item{tol}{two-sided tolerance in relative units to identify
#'         "no change" for \code{method == "direction"} and
#'         \code{method == "rdirection"}}
#'   \item{subset}{a logical vector of length \code{runsN_sites} used by
#'         \code{method == "rdirection"} to identify "no change"}
#'   }
#'
#' @return A list of three-dimensional numerical \code{arrays} with dimensions
#'   \code{runsN_sites}, \code{variables}, and \code{scenarios}; each element
#'   is named as \var{delta_ref_met} where \var{ref} is one of
#'   \code{ref_condition} and \var{met} is the first three letters of
#'   \code{method} and the value of \code{tol} if provided.
#' @export
calc_response_change_from_reference <- function(data, ref_condition,
  method = c("absolute", "relative", "direction", "rdirection"), ...) {

  method <- match.arg(method)
  dots <- list(...)
  temp <- dimnames(data)
  scenarios <- temp[[3L]]

  rtol <- sqrt(.Machine$double.eps)

  tol <- name_tol <- subset <- NULL
  if (method %in% c("direction", "rdirection")) {
    variables <- temp[[2L]]
    cat_dirs <- if (method == "direction") {
        define_direction_2deltas()
      } else {
        define_direction_1rel()
      }

    if ("tol" %in% names(dots)) {
      has_tol <- isTRUE(is.finite(dots[["tol"]]))
      tol <- if (has_tol) dots[["tol"]]
      name_tol <- paste0("_tol", sub(getOption("OutDec"), "p",
        formatC(if (has_tol) tol else 0, format = "f", flag = "0", digits = 3),
        fixed = TRUE))
    }

    if (method == "rdirection") {
      if ("subset" %in% names(dots)) {
        subset <- dots[["subset"]]
        stopifnot(is.logical(subset), length(subset) == dim(data)[1],
          !anyNA(subset))
      } else {
        subset <- rep(TRUE, dim(data)[1])
      }
    }
  }

  res <- vector("list", length = length(ref_condition))
  names(res) <- paste0("delta_", ref_condition, "_", substr(method, 1, 3),
    name_tol)

  for (k in seq_along(ref_condition)) {
    #--- Do we have a reference that is the same for all GCMs or one that
    # varies?
    sc_ref <- grep(ref_condition[k], scenarios, value = TRUE)
    gcms_under_refs <- find_reqMs(sc_ref)
    ref_per_gcm <- length(sc_ref) > 1 || !identical(sc_ref, gcms_under_refs)

    x <- array(NA, dim = dim(data), dimnames = dimnames(data))

    for (sc in seq_along(scenarios)) {
      if (ref_per_gcm) {
        # Change from a GCM-specific reference per scenario
        temp <- which(find_reqMs(scenarios[sc]) == gcms_under_refs)
        if (length(temp) == 1L) {
          sc_gcms_under_refs <- which(sc_ref[temp] == scenarios)
        } else {
          # No GCM-specific scenario, e.g., "current"
          sc_gcms_under_refs <- NULL
        }

      } else {
        # Change from the same reference for each scenario
        sc_gcms_under_refs <- sc_ref
      }

      #--- Calculate absolute/relative/directional change
      if (length(sc_gcms_under_refs) == 1L) {
        x[, , scenarios[sc]] <-
          if (method == "absolute") {
            data[, , scenarios[sc]] - data[, , sc_gcms_under_refs]

          } else if (method == "relative") {
            temp <- (data[, , scenarios[sc]] - data[, , sc_gcms_under_refs]) /
              data[, , sc_gcms_under_refs]

            # Division by zero problems
            ids0 <- abs(data[, , sc_gcms_under_refs]) < rtol
            if (any(ids0, na.rm = TRUE)) {
              temp[ids0 & data[, , scenarios[sc]] > 0] <- Inf
              temp[ids0 & abs(data[, , scenarios[sc]]) < rtol] <- 0
              temp[ids0 & data[, , scenarios[sc]] < 0] <- -Inf
            }

            temp

          } else if (method == "rdirection") {
            as.matrix(compare_direction_1rel2ref(
              dx = data[, , sc_gcms_under_refs], dy = data[, , scenarios[sc]],
              vars = variables, cl = cat_dirs, tol = tol, subset = subset))

          } else if (method == "direction") {
            as.matrix(compare_direction_2deltas(
              dx = data[, , sc_gcms_under_refs], dy = data[, , scenarios[sc]],
              vars = variables, cl = cat_dirs, tol = tol))
          }

      }
    }

    res[[k]] <- x
  }

  res
}


#' Convert \code{\link{typeof}} to \code{\link[raster]{dataType}} types
#'
#' @references Relevant code adapted from \code{`raster:::dataType<-`}
get_raster_datatype <- function(data) {
  switch(EXPR = substr(toupper(typeof(data)), 1, 5),
    LOGIC = "LOG1S",
    BYTE = "INT1U",
    SMALL = "INT2S",
    INTEG = "INT4S",
    NUMER = , FLOAT = , SINGL = , REAL = "FLT4S",
    DOUBL = "FLT8S"
    )
}


#' Converts values into a Raster*
#'
#' @param SFSW2_prj_meta An environment
#' @param data A vector or two-dimensional object.
#' @param filename A character string. Passed to \code{\link[raster]{brick}}.
#'
#' @return A RasterLayer (if \code{data} is a vector) or RasterBrick (if
#'   \code{data} is two-dimensional).
create_raster_from_variables <- function(SFSW2_prj_meta, data, filename = "") {

  # prepare locations
  loc <- SFSW2_prj_meta[["sim_space"]][["run_sites"]]
  if (!raster::compareCRS(SFSW2_prj_meta[["sim_space"]][["crs_sites"]],
    SFSW2_prj_meta[["sim_space"]][["sim_crs"]])) {

    loc <- sp::spTransform(loc,
      CRS = SFSW2_prj_meta[["sim_space"]][["sim_crs"]])
  }

  # prepare data
  nl <- NCOL(data)
  cnames <- colnames(data)

  if (!is.numeric(data)) {
    if (nl > 1) {
      for (k in seq_len(nl)) {
        temp <- try(if (is.factor(data[, k])) {
            as.integer(data[, k])
          } else {
            as.double(data[, k])
          })
        stopifnot(!inherits(temp, "try-error"))
        data[, k] <- temp
      }

    } else {
      data <- try(if (is.factor(data)) {
          as.integer(data)
        } else {
          as.double(data)
        })
      stopifnot(!inherits(data, "try-error"))
    }
  }

  if (nl == 1) {
    data <- matrix(data, ncol = 1)
  }

  # create raster, init with NAs, and add data
  ids <- NULL
  rl <- list()
  if (nl > 1) {
    filenameks <- sapply(seq_len(nl), raster::rasterTmpFile)
  }

  for (k in seq_len(nl)) {
    rk <- raster::raster(SFSW2_prj_meta[["sim_space"]][["sim_raster"]])
    rk <- raster::init(rk, fun = function(x) rep(NA, x))
    if (k == 1) {
      ids <- raster::cellFromXY(rk, sp::coordinates(loc))
    }

    rk[ids] <- data[, k]

    if (nl > 1) {
      rk <- raster::writeRaster(rk, filename = filenameks[k])
    }

    rl <- c(rl, rk)
  }

  if (nl > 1) {
    names(rl) <- cnames
    # first convert list to stack before passing to brick because
    # raster v2.9.6 the list-method of brick ignores all ... arguments
    r <- raster::brick(raster::stack(rl), filename = filename)
    unlink(filenameks)

  } else {
    r <- rl[[1]]
  }


  # set datatype
  raster::dataType(r) <- get_raster_datatype(data)

  r
}


#' Convert raster where variables are organized in the third dimension to a
#' \var{netCDF} file
#'
#' @param x A Raster* object.
#' @param var_attributes A list of named character strings each of length
#'   equal to the number of variables, i.e., the third dimension of \code{x}.
#'   If not missing, then the variable attributes will be added to the
#'   \var{netCDF} file. If no variable names are provided, then the names of
#'   \code{x} are used; otherwise, these are added as \code{original_name}.
#' @param global_attributes A list of named character strings. If not missing,
#'   then the global attributes will be added to the \var{netCDF} file.
#' @param file A character string. The file path of the \var{netCDF} file to be
#'   created.
#' @param overwrite A logical value. If \code{TRUE}, file will be overwritten
#'   if it exists.
#' @return This function is used for the side-effect of creating a file.
create_netCDF_from_raster_with_variables <- function(x, time_bounds,
  var_attributes, global_attributes, file, force_v4 = TRUE,
  overwrite = FALSE) {

  stopifnot(requireNamespace("ncdf4"))
  if (force_v4) {
    # avoid "_FillValue" error in older versions of `raster` package
    stopifnot(utils::packageVersion("raster") >= "2.9.1")
  }

  nl <- raster::nlayers(x)
  var_names <- var_longnames <- names(x)
  var_units <- rep("", nl)

  ncdf4_datatype <- raster:::.getNetCDFDType(raster::dataType(x))

  NAflag <- raster::NAvalue(x)
  if (is.infinite(NAflag)) {
    # values: \code{\link[raster]{dataType}} and \code{`raster:::dataType<-`}
    NAflag <- switch(ncdf4_datatype,
      char = , byte = , short = -128L, integer = -2147483647L,
      float = -3.4e+38, double = -1.7e+308)
  }

  # Note: raster files are organized starting from NE corner
  xvals <- raster::xFromCol(x, seq_len(raster::ncol(x)))
  yvals <- raster::yFromRow(x, seq_len(raster::nrow(x)))
  grid_halfres <- raster::res(x) / 2

  has_time_central <- !missing(time_bounds)
  if (has_time_central) {
    stopifnot(length(time_bounds) == 2L)
    time_central <- mean(time_bounds)
  }

  if (!missing(var_attributes)) {
    if ("name" %in% names(var_attributes)) {
      var_names <- var_attributes[["name"]]
      var_attributes[["name"]] <- NULL
    }

    if ("long_name" %in% names(var_attributes)) {
      var_longnames <- var_attributes[["long_name"]]
      var_attributes[["long_name"]] <- NULL
    }

    if (is.null(var_longnames) && !is.null(var_names)) {
      var_longnames <- var_names
    }

    if ("units" %in% names(var_attributes)) {
      var_units <- var_attributes[["units"]]
      var_attributes[["units"]] <- NULL
    }
  }

  if (!missing(var_attributes)) {
    ns_att_vars <- names(var_attributes)
  }


  #--- setup of netCDF file
  var_chunksizes <- c(raster::ncol(x), raster::nrow(x))
  var_start <- c(1, 1)

  if (has_time_central) {
    var_chunksizes <- c(var_chunksizes, 1L)
    var_start <- c(var_start, 1)
  }

  # define dimensions
  bnddim <- ncdf4::ncdim_def(name = "bnds", units = "", vals = seq_len(2L),
    create_dimvar = FALSE)
  xdim <- ncdf4::ncdim_def(name = "lon", longname = "Longitude",
    units = "degrees_east", vals = xvals)
  ydim <- ncdf4::ncdim_def(name = "lat", longname = "Latitude",
    units = "degrees_north", vals = yvals)
  if (has_time_central) {
    tdim <- ncdf4::ncdim_def(name = "time", units = "Gregorian_year since 1900",
      vals = time_central)
  }

  # define variables
  var_dims <- list(xdim, ydim)
  if (has_time_central) {
    var_dims <- c(var_dims, list(tdim))
  }

  var_defs <- lapply(seq_len(nl), function(k)
    ncdf4::ncvar_def(name = var_names[k], units = var_units[k],
      dim = var_dims, chunksizes = var_chunksizes, missval = NAflag,
      longname = var_longnames[k], prec = ncdf4_datatype))

  crsdef <- ncdf4::ncvar_def(name = "crs", units = "", dim = list(),
    missval = NULL, prec = "integer")

  # define dimension bounds
  lonbnddef <- ncdf4::ncvar_def(name = "lon_bnds", units = "",
    dim = list(bnddim, xdim), missval = NULL,
    chunksizes = c(2L, var_chunksizes[1]))
  latbnddef <- ncdf4::ncvar_def(name = "lat_bnds", units = "",
    dim = list(bnddim, ydim), missval = NULL,
    chunksizes = c(2L, var_chunksizes[2]))

  if (has_time_central) {
    tbnddef <- ncdf4::ncvar_def(name = "time_bnds", units = "",
      dim = list(bnddim, tdim), missval = NULL, chunksizes = c(2L, 1L))
  }

  nc_dimvars <- list(lonbnddef, latbnddef)
  if (has_time_central) {
    nc_dimvars <- c(nc_dimvars, list(tbnddef))
  }

  #--- create netCDF file
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  nc <- ncdf4::nc_create(filename = file,
    vars = c(nc_dimvars, list(crsdef), var_defs), force_v4 = force_v4)
  on.exit(ncdf4::nc_close(nc))


  #--- write values of dimension bounds:
  try(ncdf4::ncvar_put(nc, varid = "lon_bnds",
    vals = rbind(xvals - grid_halfres[1], xvals + grid_halfres[1]),
    start = c(1, 1), count = c(2L, var_chunksizes[1])))

  try(ncdf4::ncvar_put(nc, varid = "lat_bnds",
    vals = rbind(yvals + grid_halfres[2], yvals - grid_halfres[2]),
    start = c(1, 1), count = c(2L, var_chunksizes[2])))

  if (has_time_central) {
    try(ncdf4::ncvar_put(nc, varid = "time_bnds",
      vals = time_bounds,
      start = c(1, 1), count = c(2L, 1L)))
  }


  #--- add attributes
  # add dimension attributes
  ncdf4::ncatt_put(nc, "lon", "axis", "X")
  ncdf4::ncatt_put(nc, "lon", "bounds", "lon_bnds")
  ncdf4::ncatt_put(nc, "lat", "axis", "Y")
  ncdf4::ncatt_put(nc, "lat", "bounds", "lat_bnds")

  if (has_time_central) {
    ncdf4::ncatt_put(nc, "time", "axis", "T")
    ncdf4::ncatt_put(nc, "time", "bounds", "time_bnds")
    ncdf4::ncatt_put(nc, "time", "calendar", "gregorian")
  }

  # add global attributes
  ncdf4::ncatt_put(nc, varid = 0, attname = "Conventions", attval = "CF-1.4")
  ncdf4::ncatt_put(nc, varid = 0, attname = "created_by",
    attval = paste0(R.version[["version.string"]], ", R packages ",
      if (requireNamespace("rSFSW2")) {
        paste0("rSFSW2 v", utils::packageVersion("rSFSW2"), " and ")
      },
      "ncdf4 v", utils::packageVersion("ncdf4"),
      ", and ", system2("nc-config", "--version", stdout = TRUE, stderr = TRUE))
    )
  ncdf4::ncatt_put(nc, varid = 0, attname = "creation_date",
    attval = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

  if (!missing(global_attributes)) {
    ns_att_glob <- names(global_attributes)

    for (natt in ns_att_glob) {
      ncdf4::ncatt_put(nc, varid = 0, attname = natt,
        attval = global_attributes[[natt]])
    }
  }

  if (!has_time_central) {
    ncdf4::ncatt_put(nc, varid = 0, attname = "time_label",
      attval = "None")
    ncdf4::ncatt_put(nc, varid = 0, attname = "time_title",
      attval = "No temporal dimensions ... fixed field")
  }

  # add coordinate system attributes
  prj <- raster::crs(x)
  if (!is.na(prj)) {
    ncdf4::ncatt_put(nc, varid = "crs", attname = "proj4",
      attval = as.character(prj))
  }


  #--- add variables values
  for (k in seq_len(nl)) {
    # write values to variable
    vals <- raster::getValues(if (nl > 1) raster::raster(x, layer = k) else x)

    try(ncdf4::ncvar_put(nc, varid = var_names[k],
      vals = matrix(vals, ncol = var_chunksizes[2]),
      start = var_start, count = var_chunksizes))

    # add variable attributes
    for (natt in ns_att_vars) {
      ncdf4::ncatt_put(nc, varid = var_names[k], attname = natt,
        attval = var_attributes[[natt]][k])
    }

    # Flush this step to the file so we dont lose it
    # if there is a crash or other problem
    ncdf4::nc_sync(nc)
  }

  invisible(TRUE)
}

#' @describeIn create_netCDF_from_raster_with_variables
create_netCDF_from_array_with_variables <- function(x, locations, grid,
  time_bounds, var_attributes, global_attributes, file, force_v4 = TRUE,
  overwrite = FALSE) {

  stopifnot(requireNamespace("ncdf4"))
  if (force_v4) {
    # avoid "_FillValue" error in older versions of `raster` package
    stopifnot(utils::packageVersion("raster") >= "2.9.1")
  }

  if (file.exists(file)) {
    if (overwrite) {
      unlink(file)
    } else {
      stop("File ", shQuote(basename(file)), " exists and 'overwrite' is FALSE")
    }
  }

  nl <- NCOL(x)
  if (nl == 1 && is.null(dim(x))) {
    x <- matrix(x, ncol = 1, dimnames = list(NULL, names(x)))
  }

  var_names <- var_longnames <- colnames(x)
  var_units <- rep("", nl)

  if (inherits(locations, "Spatial")) {
    locations <- sp::coordinates(locations)
  }

  ncdf4_datatype <- raster:::.getNetCDFDType(get_raster_datatype(x))

  # values: \code{\link[raster]{dataType}} and \code{`raster:::dataType<-`}
  NAflag <- switch(ncdf4_datatype,
    char = , byte = , short = -128L, integer = -2147483647L,
    float = -3.4e+38, double = -1.7e+308)

  # Note: raster files are organized starting from NE corner
  xvals <- raster::xFromCol(grid, seq_len(raster::ncol(grid)))
  yvals <- raster::yFromRow(grid, seq_len(raster::nrow(grid)))
  grid_halfres <- raster::res(grid) / 2

  has_time_central <- !missing(time_bounds)
  if (has_time_central) {
    stopifnot(length(time_bounds) == 2L)
    time_central <- mean(time_bounds)
  }

  if (!missing(var_attributes)) {
    if ("name" %in% names(var_attributes)) {
      var_names <- var_attributes[["name"]]
      var_attributes[["name"]] <- NULL
    }

    if ("long_name" %in% names(var_attributes)) {
      var_longnames <- var_attributes[["long_name"]]
      var_attributes[["long_name"]] <- NULL
    }

    if (is.null(var_longnames) && !is.null(var_names)) {
      var_longnames <- var_names
    }

    if ("units" %in% names(var_attributes)) {
      var_units <- var_attributes[["units"]]
      var_attributes[["units"]] <- NULL
    }
  }

  if (!missing(var_attributes)) {
    ns_att_vars <- names(var_attributes)
  }


  #--- setup of netCDF file
  var_chunksizes <- c(raster::ncol(grid), raster::nrow(grid))
  var_start <- c(1, 1)

  if (has_time_central) {
    var_chunksizes <- c(var_chunksizes, 1L)
    var_start <- c(var_start, 1)
  }

  # define dimensions
  bnddim <- ncdf4::ncdim_def(name = "bnds", units = "", vals = seq_len(2L),
    create_dimvar = FALSE)
  xdim <- ncdf4::ncdim_def(name = "lon", longname = "Longitude",
    units = "degrees_east", vals = xvals)
  ydim <- ncdf4::ncdim_def(name = "lat", longname = "Latitude",
    units = "degrees_north", vals = yvals)
  if (has_time_central) {
    tdim <- ncdf4::ncdim_def(name = "time", units = "Gregorian_year since 1900",
      vals = time_central)
  }

  # define variables
  var_dims <- list(xdim, ydim)
  if (has_time_central) {
    var_dims <- c(var_dims, list(tdim))
  }

  var_defs <- lapply(seq_len(nl), function(k)
    ncdf4::ncvar_def(name = var_names[k], units = var_units[k],
      dim = var_dims, chunksizes = var_chunksizes, missval = NAflag,
      longname = var_longnames[k], prec = ncdf4_datatype))

  crsdef <- ncdf4::ncvar_def(name = "crs", units = "", dim = list(),
    missval = NULL, prec = "integer")

  # define dimension bounds
  lonbnddef <- ncdf4::ncvar_def(name = "lon_bnds", units = "",
    dim = list(bnddim, xdim), missval = NULL,
    chunksizes = c(2L, var_chunksizes[1]))
  latbnddef <- ncdf4::ncvar_def(name = "lat_bnds", units = "",
    dim = list(bnddim, ydim), missval = NULL,
    chunksizes = c(2L, var_chunksizes[2]))

  if (has_time_central) {
    tbnddef <- ncdf4::ncvar_def(name = "time_bnds", units = "",
      dim = list(bnddim, tdim), missval = NULL, chunksizes = c(2L, 1L))
  }

  nc_dimvars <- list(lonbnddef, latbnddef)
  if (has_time_central) {
    nc_dimvars <- c(nc_dimvars, list(tbnddef))
  }

  #--- create netCDF file
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  nc <- ncdf4::nc_create(filename = file,
    vars = c(nc_dimvars, list(crsdef), var_defs), force_v4 = force_v4)
  on.exit(ncdf4::nc_close(nc))


  #--- write values of dimension bounds:
  try(ncdf4::ncvar_put(nc, varid = "lon_bnds",
    vals = rbind(xvals - grid_halfres[1], xvals + grid_halfres[1]),
    start = c(1, 1), count = c(2L, var_chunksizes[1])))

  try(ncdf4::ncvar_put(nc, varid = "lat_bnds",
    vals = rbind(yvals + grid_halfres[2], yvals - grid_halfres[2]),
    start = c(1, 1), count = c(2L, var_chunksizes[2])))

  if (has_time_central) {
    try(ncdf4::ncvar_put(nc, varid = "time_bnds",
      vals = time_bounds,
      start = c(1, 1), count = c(2L, 1L)))
  }


  #--- add attributes
  # add dimension attributes
  ncdf4::ncatt_put(nc, "lon", "axis", "X")
  ncdf4::ncatt_put(nc, "lon", "bounds", "lon_bnds")
  ncdf4::ncatt_put(nc, "lat", "axis", "Y")
  ncdf4::ncatt_put(nc, "lat", "bounds", "lat_bnds")

  if (has_time_central) {
    ncdf4::ncatt_put(nc, "time", "axis", "T")
    ncdf4::ncatt_put(nc, "time", "bounds", "time_bnds")
    ncdf4::ncatt_put(nc, "time", "calendar", "gregorian")
  }

  # add global attributes
  ncdf4::ncatt_put(nc, varid = 0, attname = "Conventions", attval = "CF-1.4")
  ncdf4::ncatt_put(nc, varid = 0, attname = "created_by",
    attval = paste0(R.version[["version.string"]], ", R packages ",
      if (requireNamespace("rSFSW2")) {
        paste0("rSFSW2 v", utils::packageVersion("rSFSW2"), " and ")
      },
      "ncdf4 v", utils::packageVersion("ncdf4"),
      ", and ", system2("nc-config", "--version", stdout = TRUE, stderr = TRUE))
  )
  ncdf4::ncatt_put(nc, varid = 0, attname = "creation_date",
    attval = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

  if (!missing(global_attributes)) {
    ns_att_glob <- names(global_attributes)

    for (natt in ns_att_glob) {
      ncdf4::ncatt_put(nc, varid = 0, attname = natt,
        attval = global_attributes[[natt]])
    }
  }

  if (!has_time_central) {
    ncdf4::ncatt_put(nc, varid = 0, attname = "time_label",
      attval = "None")
    ncdf4::ncatt_put(nc, varid = 0, attname = "time_title",
      attval = "No temporal dimensions ... fixed field")
  }

  # add coordinate system attributes
  prj <- raster::crs(grid)
  if (!is.na(prj)) {
    ncdf4::ncatt_put(nc, varid = "crs", attname = "proj4",
      attval = as.character(prj))
  }


  #--- prepare to put values to full grid format
  val_template <- rep(NA, raster::ncell(grid))
  val_ids <- raster::cellFromXY(grid, locations)

  #--- add variables values
  for (k in seq_len(nl)) {
    # put values into full grid format
    temp <- val_template
    temp[val_ids] <- x[, k]
    vals <- matrix(temp, ncol = var_chunksizes[2])

    # write values to variable
    try(ncdf4::ncvar_put(nc, varid = var_names[k], vals = vals,
      start = var_start, count = var_chunksizes))

    # add variable attributes
    for (natt in ns_att_vars) {
      ncdf4::ncatt_put(nc, varid = var_names[k], attname = natt,
        attval = var_attributes[[natt]][k])
    }

    # Flush this step to the file so we dont lose it
    # if there is a crash or other problem
    ncdf4::nc_sync(nc)
  }

  invisible(TRUE)
}



#' Enhance \code{\link[raster]{raster}} to handle added information written
#' to a \var{netCDF} file by function
#' \code{\link{create_netCDF_from_raster_with_variables}}
#'
#' @return A Raster* object.
read_netCDF_to_raster <- function(x, ...) {
  r <- raster::raster(x, ...)

  # Check whether projection was read correctly
  r_crs <- raster::crs(r)
  r_has_crs <- inherits(r_crs, "CRS") && !is.na(r_crs) &&
    rgdal::checkCRSArgs(raster::crs(r, asText = TRUE))[[1]]

  if (!r_has_crs) {
    nc <- RNetCDF::open.nc(x)
    nc_crs <- RNetCDF::att.get.nc(nc, variable = "crs", attribute = "proj4")
    RNetCDF::close.nc(nc)

    nc_crs <- raster::crs(nc_crs)
    if (inherits(nc_crs, "CRS") && !is.na(nc_crs) &&
        rgdal::checkCRSArgs(raster::crs(nc_crs, asText = TRUE))[[1]]) {
      raster::crs(r) <- nc_crs
    } else {
      warning("'read_netCDF_to_raster': could not locate a valid projection.")
    }
  }

  r
}

#' Load simulation output from \pkg{rSFSW2} into an array either from
#' a \var{sqlite3} database or from \var{netCDF} files
#'
load_rSFSW2_data_for_analysis <- function(meta, path, path_tmp,
  fname_dbOuts = NULL, whereClause = NULL, subprojects = NULL,
  variables, sim_scenarios, sc_historical, subset = NULL,
  write_to_netcdf = FALSE, ftag_gatt, ftag_gatt2, var_names, var_units,
  timeaggs) {

  # Check inputs
  has_subprojects <- !is.null(subprojects)
  if (!has_subprojects) {
    subprojects <- "Default"
  }

  if (is.null(subset)) {
    subset <- rep(TRUE, meta[["sim_size"]][["runsN_sites"]])
  }


  # Output container
  res <- if (has_subprojects) {
      array(NA,
        dim = c(meta[["sim_size"]][["runsN_sites"]],
          length(variables), length(sim_scenarios), length(subprojects)),
        dimnames = list(NULL, variables, sim_scenarios, subprojects))
    } else {
      array(NA,
        dim = c(meta[["sim_size"]][["runsN_sites"]],
          length(variables), length(sim_scenarios)),
        dimnames = list(NULL, variables, sim_scenarios))

    }

  # Check whether a previous call with write_to_netcdf = TRUE created netCDFs
  fname_datall <- list.files(path = path,
    pattern = "(^All_LyrC_SOILWAT2-)[[:print:]]+(nc$)", full.names = TRUE)


  if (length(fname_datall) == length(sim_scenarios) * length(subprojects)) {
    #--- Read from netCDF files and convert to array

    rdim <- dim(meta[["sim_space"]][["sim_raster"]])
    loc <- sp::coordinates(meta[["sim_space"]][["run_sites"]])
    id_map <- cbind(
      row = raster::colFromX(meta[["sim_space"]][["sim_raster"]], loc[, 1]),
      col = raster::rowFromY(meta[["sim_space"]][["sim_raster"]], loc[, 2])
    )

    for (sc in sim_scenarios) for (sprj in subprojects) {
      getM <- cur_to_hist(find_reqMs(sc))
      temp_pattern <- if (sc == sc_historical) {
          getM
        } else {
          paste0(getM, "_", find_reqCSs(sc))
        }
      if (has_subprojects) {
        temp_pattern <- paste0(temp_pattern, "_", sprj)
      }

      fname <- fname_datall[grep(temp_pattern, basename(fname_datall))]

      nc <- RNetCDF::open.nc(fname)
      dtemp <- RNetCDF::read.nc(nc)
      RNetCDF::close.nc(nc)

      stopifnot(
        length(dtemp[["lat"]]) == rdim[1],
        length(dtemp[["lon"]]) == rdim[2])

      for (iv in seq_along(variables)) {
        if (has_subprojects) {
          res[, variables[iv], sc, sprj] <- dtemp[[var_names[iv]]][id_map]
        } else {
          res[, variables[iv], sc] <- dtemp[[var_names[iv]]][id_map]
        }
      }

      # Check that data are correctly ordered
      if ("site_id" %in% variables) {
        stopifnot(identical(
          as.integer(if (has_subprojects) {
              res[, "site_id", sc, sprj]
            } else {
              res[, "site_id", sc]
            }),
          meta[["sim_size"]][["runIDs_sites"]]))
      }
    }

  } else {
    #--- Read from SQLite3 database and convert to array

    # Check inputs
    if (is.null(fname_dbOuts)) {
      fname_dbOuts <- meta[["fnames_out"]][["dbOutput"]]
      names(fname_dbOuts) <- subprojects
    }

    if (has_subprojects) {
      # Check that file names of dbOutput and subprojects have correct length
      if (length(fname_dbOuts) == 1 && length(subprojects) > 1) {
        fname_dbOuts <- rep_len(fname_dbOuts, length(subprojects))
      }

      stopifnot(length(fname_dbOuts) == length(subprojects))

      # Check that whereClause and subprojects do not interfere
      if (isTRUE(grepl("Experimental_Label", whereClause)) &&
          isTRUE(subprojects != "Default")) {
        stop("'whereClause' cannot set 'Experimental_Label' if there are ",
          "'subprojects'")
      }
    }


    # Extract data from simulation database `dbOutput`
    for (sprj in subprojects) {
      if (has_subprojects) {
        whereClause_temp <- paste0("Experimental_Label='", sprj, "'")
        if (!is.null(whereClause) && nchar(whereClause) > 0) {
          whereClause_temp <- paste(whereClause, "AND", whereClause_temp)
        }

       } else {
        whereClause_temp <- whereClause
      }

      fname_out <- if (has_subprojects) {
          file.path(path_tmp, "temp_dbOut", sprj)
        } else {
          file.path(path_tmp, "temp_dbOut")
        }

      temp <- extract_dbOut_to_df(meta,
        fname_dbOut = fname_dbOuts[sprj],
        variables = variables, scenarios = sim_scenarios,
        whereClause = whereClause_temp,
        file = fname_out, verbose = TRUE)

      # Subset to adjusted study area
      icol_siteid <- which("site_id" == dimnames(res)[[2]])
      temp[!subset, -icol_siteid, ] <- NA

      # Copy to output container
      if (has_subprojects) {
        res[, variables, , sprj] <- temp
      } else {
        res[, variables, ] <- temp
      }


      # Prepare simulation data for sharing
      if (write_to_netcdf) {
        for (sc in sim_scenarios) {
          ftag_gatt2_temp <- ftag_gatt2
          ftag_gatt2_temp[["experiment_id"]] <- cur_to_hist(find_reqMs(sc))
          if (sc != sc_historical) {
            ftag_gatt2_temp[["parent_source_id"]] <- find_reqCSs(sc)
          }

          time_bounds <- if (sc == sc_historical) {
            c(meta[["sim_time"]][["startyr"]],
              meta[["sim_time"]][["endyr"]])
          } else {
            fut_yrs <- meta[["sim_time"]][["future_yrs"]][find_reqDeltaYR(sc), ]
            c(fut_yrs[["DSfut_startyr"]], fut_yrs[["DSfut_endyr"]])
          }

          fname_nc <- file.path(path,
            paste0("All_", # all variables together
              "LyrC_", # Land-realm; year-climatology
              "SOILWAT2-",
              if (sc == sc_historical) {
                cur_to_hist(sc)
              } else {
                paste0(ftag_gatt2_temp[["experiment_id"]], "_",
                  ftag_gatt2_temp[["parent_source_id"]])
              },
              "_",
              if (has_subprojects) {
                paste0(sprj, "_")
              },
              "gn_", # grid-native
              paste(time_bounds, collapse = "-"),
              ".nc"))

          if (!file.exists(fname_nc)) {
            # Write raster to netCDF file
            create_netCDF_from_array_with_variables(
              x = if (has_subprojects) res[, , sc, sprj] else res[, , sc],
              locations = meta[["sim_space"]][["run_sites"]],
              grid = meta[["sim_space"]][["sim_raster"]],

              time_bounds = time_bounds,

              var_attributes = list(
                name = var_names,
                original_name = variables,
                units = var_units,
                cell_methods = paste("time:", timeaggs)
              ),

              global_attributes = c(ftag_gatt, ftag_gatt2_temp,
                if (has_subprojects) {
                  list(variant_info = paste("forcing:", sprj))
                },
                cell_measures = "data variables provided at cell centers"
              ),

              file = fname_nc
            )
          }
        }
      }
    }
  }

  res
}



#' Calculate cell-wise ensembles across `reqCSs`
#'
#' @param data A data.frame with three dimensions -- as generated by
#'   \code{\link{extract_dbOut_to_df}}.
#' @param subset A logical vector used to subset gridcells from \code{data}
#'   for calculation of ensemble.
#' @param funs A vector of character strings or of functions. The functions of
#'   the form \code{f(x, ...)} which will be used to calculate the ensemble
#'   values across \var{GCMs}.
#' @param ... Optional arguments passed to every \code{funs} function, e.g.,
#'   \code{na.rm = TRUE}.
#' @param variables A vector of character strings. A subset of names of the
#'   second dimension of \code{data} or \code{NULL} which indicates to use
#'   the full set of the second dimension of \code{data}.
#' @param reqCSs A vector of character strings. A subset of (climate) scenarios
#'   from \code{SFSW2_prj_meta[["sim_scens"]][["reqCSs"]]} for which ensembles
#'   are calculated.
#'
#' @export
calc_cellwise_ensemble <- function(SFSW2_prj_meta, data, subset = NULL,
  funs = c("mean", "min", "max"), ..., variables = NULL, add_historic = TRUE,
  reqCSs = NULL, reqMs = NULL, verbose = FALSE) {

  nfuns <- length(funs)

  dim_data <- dim(data)
  temp <- dimnames(data)
  has_sim_scens_id <- temp[[3]]
  vars_data <- temp[[2]]

  if (is.null(subset)) {
    subset <- rep(TRUE, dim_data[1])
  }

  # Check that data is well formed
  stopifnot(
    identical(length(dim_data), 3L),
    identical(dim_data[1], length(subset)),
    identical(dim_data[1], SFSW2_prj_meta[["sim_size"]][["runsN_sites"]])
    )

  # Locate scenarios
  if (is.null(reqCSs)) {
    reqCSs <- unique(find_reqCSs(has_sim_scens_id, SFSW2_prj_meta))
  }

  sc_current <- SFSW2_prj_meta[["sim_scens"]][["ambient"]]
  has_current <- sc_current %in% has_sim_scens_id
  reqCS2s <- if (has_current) c(sc_current, reqCSs) else reqCSs

  # Locate variables
  if (is.null(variables)) {
    variables <- vars_data
  } else {
    stopifnot(variables %in% vars_data)
  }

  # Prepare ensemble data object

  data_ens <- array(NA,
    dim = c(nfuns, SFSW2_prj_meta[["sim_size"]][["runsN_sites"]],
      length(variables), length(reqCS2s)),
    dimnames = list(funs, NULL, variables, reqCS2s))

  # Copy current values
  if (has_current && add_historic) for (k in seq_len(nfuns)) {
    data_ens[k, subset, , sc_current] <- data[subset, , sc_current]
  }

  # Calculate ensemble values
  for (sc in seq_along(reqCSs)) {
    if (verbose) {
      print(paste(Sys.time(), "'calc_cellwise_ensemble':", reqCSs[sc]))
    }

    isc <- grep(reqCSs[sc], has_sim_scens_id)

    if (length(isc) > 0) {
      data_ens[, subset, , reqCSs[sc]] <- apply(
        data[subset, , has_sim_scens_id[isc], drop = FALSE],
        MARGIN = 1:2,
        function(x, ...) sapply(funs, do.call, args = list(x = x, ...)))
    }
  }

  data_ens
}



#' Assign to cells the `inverse-ecdf` quantiles from extent-wide (area-weighted)
#' ranked, (equally-weighted) GCMs ensembles across `reqCSs`
#'
#' This ensemble approach, unlike cell-wise ensembles, accounts for physical
#' interdependence among grid-points within a GCM projection
#' (Madsen et al. 2017).
#'
#' @param data A data.frame with three dimensions -- as generated by
#'   \code{\link{extract_dbOut_to_df}}.
#' @param subset A logical vector used to subset gridcells from \code{data}
#'   for calculation of ensemble.
#' @param area A numeric vector. Its length corresponds to the first dimension
#'   of \code{data} and represents the cell areas (e.g., that may vary by
#'   latitude). If not \code{NULL} or not all equal, then a \code{fcentral}
#'   value "mean" uses \code{\link[stats]{weighted.mean}} instead of
#'   \code{\link[base]{mean}} and a value of "median" uses
#'   \code{\link[Hmisc]{wtd.quantile}} instead of \code{\link[stats]{median}}
#'   where \code{area} are used as non-random "reliability" weights.
#' @param fcentral A character string naming the function to calculate a central
#'   tendency across the spatial extent for each \var{GCMs};
#'   one of "mean" or "median". The function is used to calculate the
#'   region-wide (weighted) statistic on which the GCMs are ranked. See
#'   \code{area}.
#' @param funs A vector of character strings or of functions. The functions of
#'   the form \code{f(x, ...)} which will be used to calculate the ensemble
#'   values across \var{GCMs}. Note: \var{"min"}, \var{"median"}, \var{"max"}
#'   will be replaced by their quantile equivalents,
#'   i.e., \code{probs = 0, 0.5, or 1}  respectively.
#' @param probs A numeric vector of probabilities with values in \code{[0,1]}
#'   to determine which quantiles from among the \code{fcentral} values
#'   across \var{GCMs} are returned for each cell.
#' @param ... Optional arguments passed to \code{\link[stats]{quantile}}, e.g.,
#'   \code{na.rm = TRUE}.
#' @param ties.method A character string. Specifies how ties are treated. See
#'   \code{link[base]{rank}}. The method must produce unique ranks if
#'   \code{funs} or \code{probs} are not empty.
#' @param variables A vector of character strings. A subset of names of the
#'   second dimension of \code{data} or \code{NULL} which indicates to use the
#'   full set of the second dimension of \code{data}.
#' @param reqCSs A vector of character strings. A subset of (climate) scenarios
#'   from \code{SFSW2_prj_meta[["sim_scens"]][["reqCSs"]]} for which ensembles
#'   are calculated.
#' @param reqMs A vector of character strings. A subset of (climate) models
#'   from \code{SFSW2_prj_meta[["sim_scens"]][["reqMs"]]} across which ensembles
#'   are calculated.
#' @param verbose A logical value.
#'
#' @return A 4-dimensional, numeric array where the first dimension represents
#'   \code{probs}, the second the cells/sites (first dimension of \code{data}),
#'   the third the \code{variables}, and the fourth the scenarios (ambient plus
#'   \code{reqCSs}).
#'
#' @references Madsen, M. S., P. L. Langen, F. Boberg, and J. H. Christensen.
#'   2017. Inflated Uncertainty in Multimodel-Based Regional Climate
#'   Projections. Geophysical Research Letters 44:11606-11613.
#'
#' @export
calc_extentwise_ensemble <- function(SFSW2_prj_meta, data, subset = NULL,
  area = NULL, fcentral = c("median", "mean"), funs = c("min", "median", "max"),
  probs = c(0, 0.5, 1), ..., ties.method = "first",
  variables = NULL, add_historic = TRUE, reqCSs = NULL, reqMs = NULL,
  verbose = FALSE) {

  dim_data <- dim(data)
  dnames_data <- dimnames(data)
  if (is.null(area)) {
    area <- rep(1, SFSW2_prj_meta[["sim_size"]][["runsN_sites"]])
  }
  if (is.null(subset)) {
    subset <- rep(TRUE, dim_data[1])
  }


  # Check that data is well formed
  stopifnot(
    identical(length(dim_data), 3L),
    identical(length(area), SFSW2_prj_meta[["sim_size"]][["runsN_sites"]]),
    identical(dim_data[1], length(subset)),
    identical(dim_data[1L], SFSW2_prj_meta[["sim_size"]][["runsN_sites"]]))

  # Determine central tendency function
  fcentral <- fun_central(method = fcentral, area = area[subset], na.rm = TRUE)

  # Determine aggregation functions: quantiles and other functions
  names_aggs <- funs

  # convert equivalent funs to probs: 0% == min, 50% == median, 100% == max
  equivalent_p_f <- data.frame(
    probs = c(0, 0.5, 1),
    funs = c("min", "median", "max"),
    stringsAsFactors = FALSE)

  for (k in seq_len(nrow(equivalent_p_f))) {
    has_q <- abs(probs - equivalent_p_f[k, "probs"]) <= sqrt(.Machine$double.eps)
    has_f <- funs %in% equivalent_p_f[k, "funs"]
    if (any(has_f)) {
      funs <- funs[!has_f]
      if (!any(has_q)) {
        probs <- c(probs, equivalent_p_f[k, "probs"])
      }
    }
  }

  probs <- sort(unique(probs))
  nprobs <- length(probs)
  iout_probs <- match(equivalent_p_f[, "funs"], names_aggs, nomatch = 0)
  nfuns <- length(funs)
  naggs <- nprobs + nfuns

  # If funs include `span` and/or `agreement`, then we need 0, 0.5, and 1
  i_span <- match("span", names_aggs, nomatch = 0)
  has_span <- i_span > 0
  i_agreement <- match("agreement", names_aggs, nomatch = 0)
  has_agreement <- i_agreement > 0
  probs_temp <- if (has_span || has_agreement) c(0, 0.5, 1)
  nprobs_temp <- length(probs_temp)

  # Find "other" funs
  i_ofuns <- seq_along(names_aggs)[-c(iout_probs, i_span, i_agreement)]
  n_ofuns <- length(i_ofuns)


  # Check how to resolve ties in ranks
  if ((nprobs > 0 || nprobs_temp > 0) &&
      !(ties.method %in% c("first", "last"))) {
    warning("Unique ranks required: `ties.method` set to 'first'.")
    ties.method <- "first"
  }

  # Locate GCMs and scenarios
  sc_current <- SFSW2_prj_meta[["sim_scens"]][["ambient"]]
  has_sim_scens_id <- dnames_data[[3L]]
  has_current <- sc_current %in% has_sim_scens_id

  if (is.null(reqCSs)) {
    reqCSs <- unique(find_reqCSs(has_sim_scens_id, SFSW2_prj_meta))
  }

  reqCS2s <- if (has_current) c(sc_current, reqCSs) else reqCSs

  if (is.null(reqMs)) {
    reqMs <- unique(find_reqMs(has_sim_scens_id, SFSW2_prj_meta))
  }

  # Locate variables
  vars_data <- dnames_data[[2L]]
  if (is.null(variables)) {
    variables <- vars_data
  } else {
    stopifnot(variables %in% vars_data)
  }

  # Prepare ensemble data object
  data_ens <- array(NA,
    dim = c(naggs, SFSW2_prj_meta[["sim_size"]][["runsN_sites"]],
      length(variables), length(reqCS2s)),
    dimnames = list(names_aggs, NULL, variables, reqCS2s))

  data_ens2 <- array(NA,
    dim = c(2, length(variables), length(reqMs), length(reqCSs)),
    dimnames = list(c("fcentral", "Rank"), variables, reqMs, reqCSs))


  # Copy current values
  if (has_current && add_historic) for (k in seq_len(naggs)) {
    data_ens[k, subset, , sc_current] <- data[subset, , sc_current]
  }

  # Calculate ensemble values
  for (sc in seq_along(reqCSs)) {
    if (verbose) {
      print(paste(Sys.time(), "'calc_extentwise_ensemble':", reqCSs[sc]))
    }

    # Identify available GCMs in ensemble
    isc <- grep(reqCSs[sc], has_sim_scens_id)
    igcms <- sapply(reqMs, function(m) any(grepl(m, has_sim_scens_id[isc])))
    k <- length(isc)

    # Ensembles
    if (k > 0) {
      temp_data <- data[subset, , isc, drop = FALSE]

      # Calculate region-wide `fcentral` (area-weighted) values for each
      # variable and GCM
      temp_data_fin <- temp_data

      # Deal with infinite values that arose from relative deltas
      # (division by zero): ignore for the calculation of overall central value;
      # propagate otherwise
      ids_inf <- is.infinite(temp_data)
      if (any(ids_inf)) {
        temp_data_fin[ids_inf] <- NA
      }

      data_ens2["fcentral", , igcms, reqCSs[sc]] <-
        apply(temp_data_fin, 2:3, fcentral)

      # Rank available GMCs according to region-wide `fcentral`
      data_ens2["Rank", , igcms, reqCSs[sc]] <-
        if (k > 1) {
          t(apply(data_ens2["fcentral", , igcms, reqCSs[sc]], 1, rank,
            ties.method = ties.method))
        } else {
          rank(data_ens2["fcentral", , igcms, reqCSs[sc]],
            ties.method = ties.method)
        }


      # Identify (equally-weighted) GCMs by matching probs to the quantiles of
      # the 'Inverse of empirical distribution function' type
      if (nprobs > 0) {
        id_probs <- quantile(seq_len(k), probs = probs, type = 1L,
          names = FALSE)
      }
      if (nprobs_temp > 0) {
        id_probs_temp <- quantile(seq_len(k), probs = probs_temp, type = 1L,
          names = FALSE)
      }

      # Apply functions that use ranked values
      for (iv in seq_along(variables)) {
        if (nprobs > 0) {
          # Match probs to available ranks
          id_ranks <- match(id_probs,
            data_ens2["Rank", variables[iv], igcms, reqCSs[sc]])

          # Extract ensemble quantiles
          data_ens[iout_probs, subset, variables[iv], 1L + sc] <-
            t(temp_data[, variables[iv], id_ranks])
        }

        if (nprobs_temp > 0) {
          # Extract ensemble quantiles for `span` and/or `agreement`
          temp <- temp_data[, variables[iv], match(id_probs_temp,
            data_ens2["Rank", variables[iv], igcms, reqCSs[sc]])]

          # Calculate funs
          if (has_span) {
            data_ens[i_span, subset, variables[iv], 1L + sc] <-
              apply(temp, 1, span, ...)
          }

          if (has_agreement) {
            # count agreement with sign of cell values for median GCM
            temp2 <- cbind(median = temp[, 2], temp_data[, variables[iv], ])
            data_ens[i_agreement, subset, variables[iv], 1L + sc] <-
              apply(temp2, 1, function(x)
                agreement(x[-1], val = x[1], ...))
          }
        }
      }

      # Apply other functions across GCMs
      if (n_ofuns > 0) {
        data_ens[i_ofuns, subset, , 1L + sc] <-
          apply(temp_data, MARGIN = 1:2, function(x, ...)
            sapply(names_aggs[i_ofuns], do.call, args = list(x = x, ...)))
      }
    }
  }

  list(ensemble_values = data_ens, ensemble_structure = data_ens2)
}


#' Assign to cells the `inverse-ecdf` quantiles from region-wide ranked GCMs
#' ensembles across `reqCSs`
#'
#' For each \code{region}, the function \code{\link{calc_extentwise_ensemble}}
#' is called on the set of cells belonging to that region.
#'
#' @inheritParams calc_extentwise_ensemble
#' @param region A character or numeric vector of a length equal to the first
#'   dimension of \code{data}, i.e., cells.
#' @seealso \code{\link{calc_extentwise_ensemble}}
#'
#' @export
calc_regionwise_ensemble <- function(SFSW2_prj_meta, data, subset = NULL,
  region, area = NULL, fcentral = c("median", "mean"),
  funs = c("min", "median", "max"), probs = c(0, 0.5, 1), ...,
  ties.method = "first", variables = NULL, add_historic = TRUE, reqCSs = NULL,
  verbose = FALSE) {

  stopifnot(
    length(region) == dim(data)[1L],
    is.null(area) || length(area) == dim(data)[1L],
    is.null(subset) || length(subset) == dim(data)[1L])
  region_set <- na.exclude(unique(region))
  res <- NULL

  # Determine ensembles for each region seperately and combine cells back
  # together
  for (k in seq_along(region_set)) {
    if (verbose) {
      print(paste(Sys.time(), "'calc_regionwise_ensemble':", region_set[k]))
    }

    # Set all values not in region to NA
    id_region <- region %in% region_set[k]

    temp_data <- data
    temp_data[!id_region, , ] <- NA

    # Calculate region-wise ensemble
    x <- calc_extentwise_ensemble(SFSW2_prj_meta = SFSW2_prj_meta,
      data = temp_data, subset = subset, area = area,
      fcentral = fcentral, funs = funs, probs = probs, ...,
      ties.method = ties.method, variables = variables,
      add_historic = add_historic, reqCSs = reqCSs,
      verbose = verbose)[["ensemble_values"]]

    if (k == 1) {
      res <- array(NA, dim = dim(x), dimnames = dimnames(x))
    }

    # Store values for region
    res[, id_region, , ] <- x[, id_region, , ]
  }

  res
}



#' Ensemble from ranked GCM at the scale of interest (Madsen et al. 2017)
#' @references
#'   Madsen, M. S., P. L. Langen, F. Boberg, and J. H. Christensen. 2017.
#'   Inflated Uncertainty in Multimodel-Based Regional Climate Projections.
#'   Geophysical Research Letters 44:11606–11613.
#'
calculate_ensembles <- function(meta, data, data_names = names(data), subset,
  cell_area_km2, id_region,
  sc_historical, req_Downs, req_dTime,
  ens_wises = c("EnsCW", "EnsRW", "EnsGW"),
  ens_funs = c("min", "mean", "median", "max", "span", "agreement", "majority"),
  path, ftag) {

  ens_wises <- match.arg(ens_wises, several.ok = TRUE)
  ens_funs <- match.arg(ens_funs, several.ok = TRUE)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  if (missing(subset)) {
    subset <- rep(TRUE, dim(data[[1]])[1])
  }

  # Create unique file id number for memoization
  temp <- list(
    subset = as.vector(subset),
    cell_area_km2 = as.numeric(cell_area_km2),
    id_region = as.vector(id_region),
    sc_historical = as.character(sc_historical),
    req_Downs = as.character(req_Downs),
    req_dTime = as.character(req_dTime),
    ens_wises = as.character(ens_wises),
    ens_funs = as.character(ens_funs))
  cid <- digest::digest(temp, algo = "sha1")

  # Output container
  dats_Ens <- vector(mode = "list", length = length(ens_wises))
  names(dats_Ens) <- ens_wises

  for (ew in ens_wises) {
    dats_Ens[[ew]] <- vector("list", length = length(data_names))
    names(dats_Ens[[ew]]) <- data_names

    for (k in seq_along(data_names)) {
      print(paste(Sys.time(), dQuote(ew), "ensembles for", data_names[k]))

      temp <- list(data = data[[data_names[k]]],
        data_names = as.character(data_names[k]))
      fid <- digest::digest(temp, algo = "sha1")

      fname_EnsW <- file.path(path, paste0(ftag, "_", cid, "_", fid, ".rds"))

      if (file.exists(fname_EnsW)) {
        dats_Ens[[ew]][[data_names[k]]] <- readRDS(fname_EnsW)

      } else {
        temp <-
          if (ew == "EnsCW") {
            # (i) cell-wise ensembles: for each cell, foo across GCMs are
            # extracted from the GCM with the foo-rank based on values for each
            # cell independently
            calc_cellwise_ensemble(meta,
              data = data[[data_names[k]]], subset = subset, funs = ens_funs,
              na.rm = TRUE, verbose = TRUE)

          } else if (ew == "EnsRW") {
            # (ii) region-wise ensembles: for each cell, foo across GCMs are
            # extracted from the GCM with the foo-rank based on their regional
            # mean values
            calc_regionwise_ensemble(meta,
              data = data[[data_names[k]]], subset = subset, region = id_region,
              area = cell_area_km2, fcentral = "median", funs = ens_funs,
              na.rm = TRUE, verbose = TRUE)

          } else if (ew == "EnsGW") {
            # (iii) global ensembles: for each cell, foo across GCMs are
            # extracted from the GCM with the foo-rank based on their global
            # mean values
            calc_extentwise_ensemble(meta,
              data = data[[data_names[k]]], subset = subset,
              area = cell_area_km2, fcentral = "median", funs = ens_funs,
              na.rm = TRUE, verbose = TRUE)[["ensemble_values"]]
          }

        # Cast to same structure as `data`:
        # (i) fold 1st into 4th dimension
        dats_Ens[[ew]][[data_names[k]]] <- reshape2::acast(reshape2::melt(temp),
          Var2 ~ Var3 ~ Var4 + Var1)
        # (ii) fix naming scheme of scenario/3rd dimension
        sctemp <- dimnames(dats_Ens[[ew]][[data_names[k]]])[[3]]
        sctemp2 <- strsplit(sctemp, split = "_")
        itemp <- grep(sc_historical, sctemp)
        for (k2 in itemp) {
          sctemp[k2] <- paste(req_Downs, req_dTime[1], sctemp2[[k2]][1],
            sctemp2[[k2]][2], sep = ".")
        }
        for (k2 in seq_along(sctemp)[-itemp]) {
          sctemp[k2] <- paste(req_Downs, req_dTime[2], sctemp2[[k2]][1],
            sctemp2[[k2]][2], sep = ".")
        }
        dimnames(dats_Ens[[ew]][[data_names[k]]])[[3]] <- sctemp

        saveRDS(dats_Ens[[ew]][[data_names[k]]], file = fname_EnsW)
      }
    }


    # Make sure that we didn't introduce NAs
    hasNAs <- sapply(dats_Ens[[ew]], function(x) {
      temp <- grep(sc_historical, dimnames(x)[[3]])
      anyNA(x[subset, , -temp])
    })

    if (any(hasNAs)) {
      stop("We have NAs in: ", paste(shQuote(names(hasNAs)[hasNAs]),
        collapse = ", "))
    }
  }

  dats_Ens
}


#' Main climate classification according to Trewartha
#'
#' The main climate classifications based on mean monthly temperature, i.e.,
#' without group B--dry, include \describe{
#'   \item Group A (tropical): sum(months >= 18 C) == 12
#'   \item Group C (subtropical):
#'     sum(months >= 10 C) >= 8 & sum(months >= 18 C) < 12
#'   \item Group D (temperate and continental): sum(months >= 10 C) %in% 4:7
#'   \item Group E (boreal): sum(months >= 10 C) %in% 1:3
#'   \item Group F (polar): sum(months >= 10 C) == 0
#' }
#'
#' @param meanTmonthly_C A numeric vector of length 12. Mean monthly air
#'   temperature in degree Celsius.
#'
#' @seealso \code{\link[ClimClass](koeppen_geiger)}
#' @references Trewartha, G.T. and Lyle, H.H., 1980: An Introduction to Climate.
#'   MacGraw - Hill, 5th Ed. Appendix: Koeppen's Classification of Climates
trewartha_climate <- function(meanTmonthly_C) {
  stopifnot(length(meanTmonthly_C) == 12)

  temp18 <- sum(meanTmonthly_C >= 18)

  if (temp18 == 12) {
    "tropical"
  } else {
    temp10 <- sum(meanTmonthly_C >= 10)

    if (temp10 >= 8) {
      "subtropical"
    } else if (temp10 >= 4) {
      "temperate"
    } else if (temp10 >= 1) {
      "boreal"
    } else {
      "polar"
    }
  }
}




#' Create and plot a continuous color legend
#'
#' @param zlim A numeric vector of length two. The smallest and largest z-value
#'   of \code{grid} for which a complete color gradient and legend should be
#'   printed.
#' @param zextreme A numeric vector length two. The minimum and maximum z-value
#'   of \code{grid}. If \code{zextreme[1] < zlim[1] || zextreme[2] > zlim[2]},
#'   then the values in the ranges between each zextreme and zlim will be
#'   highlighted as extreme areas in a different color.
#' @param col_desc A list with the color names. A returned object from a call to
#'   the function \code{\link{make_colors}}.
#' @param grid A RasterLayer object for which the legend is created.
#' @param box A numeric vector of length four. The SW (lower left) and NE (upper
#'   right) corners of the color legend on the plotted \code{grid}.
#' @param whitebox A logical value. If \code{TRUE}, then the area of the color
#'   legend is printed white before the colors of the legend are added.
#' @param horiz A logical value. If \code{TRUE}, then the values of the color
#'   legend are printed below the legend box; if \code{FALSE}, then the values
#'   are printed on the right side of the box.
#' @param signed A numeric value. A scaling value multiplied with the z-values
#'   for determining legend values.
#' @param fun_inv_ens A function or \code{NULL}. If a function, then it will be
#'   applied to z-values for determining legend values.
#' @param at The points at which tick-marks are to be drawn. By default
#'   (when \code{NULL}) tickmark locations are computed,
#'   see \code{\link[graphics]{axis}}.
#' @param tick A logical value. Indicates whether or not tickmarks should be
#'   added to the legend color ramp.
#' @param labels The annotations to be written at the tick-marks,
#'   see \code{\link[graphics]{axis}}.
#' @param srt A numerical value. String rotation in degrees, see
#'   \code{\link[graphics]{text}}.
#' @param cex A numerical value. Numeric character expansion factor, see
#'   \code{\link[graphics]{text}}.
add_legend <- function(zlim, zextreme, col_desc, grid,
  box = c(-100, -97, -50, -10), whitebox = TRUE, horiz = FALSE, signed = 1,
  fun_inv_ens = NULL, at = NULL, tick = TRUE, labels = TRUE, srt = 90,
  cex = 1) {

  stopifnot(requireNamespace("raster"))
  if (is.null(fun_inv_ens)) fun_inv_ens <- function(x) x
  if (missing(zextreme) || is.null(zextreme)) zextreme <- zlim

  #--- Color legend
  # Create color legend
  zr <- raster::raster(xmn = box[1], xmx = box[2], ymn = box[3], ymx = box[4],
    crs = raster::projection(grid), resolution = raster::res(grid), vals = NULL)
  zr[] <- if (horiz) {
      rep(1:dim(zr)[2], times = dim(zr)[1])
    } else {
      rep(dim(zr)[1]:1, each = dim(zr)[2])
    }

  # Draw to plot
  if (whitebox) {
    raster::image(zr, col = "white", add = TRUE)
  }
  raster::image(zr, col = col_desc[["colors_label"]], add = TRUE)

  #--- Annotate
  # Calculate tickmarks on `zlim`-scale
  # Positions: generate default if `at` is NULL
  atz <- if (is.null(at)) pretty(zlim, n = 6L) else at

  # Beautify tickmarks
  temp2z <- atz <= zlim[1]
  if (any(temp2z)) {
    # adjust position(s) below the limit to one tickmark at the limit
    ids_temp2z <- which(temp2z)
    atz <- c(zlim[1], atz[!temp2z])
  }
  temp2z <- atz >= zlim[2]
  if (any(temp2z)) {
    # adjust position(s) above the limit to one tickmark at the limit
    ids_temp2z <- which(temp2z)
    atz <- c(atz[!temp2z], zlim[2])
  }

  datz <- diff(atz)
  temp <- 0.8 / (if (abs(srt) > 45) cex else 1)
  if (length(datz) > 0 && length(temp <- which(datz < temp * max(datz))) > 0) {
    # remove tickmarks if too close together, but not limits and not 0
    id_remove <- findInterval(x = 1 + temp, vec = seq_along(atz),
      all.inside = TRUE)
    temp <- which(0 == atz[id_remove])
    if (length(temp) > 0) {
      id_remove <- id_remove[-temp]
    }
    if (length(id_remove) > 0) {
      atz <- atz[-id_remove]
    }
  }
  if (length(datz) == 0) {
    atz <- zextreme
  }

  # Calculate tickmarks on `color`-scale
  ext <- raster::extent(zr)
  rtemp <- col_desc[["ncol_label_added"]] / length(col_desc[["colors_label"]])

  if (horiz) {
    xmin_orig <- ext@xmin
    if (col_desc[["added_below"]]) {
      xmin_orig <- xmin_orig + rtemp * (ext@xmax - ext@xmin)
    }

    xmax_orig <- ext@xmax
    if (col_desc[["added_above"]]) {
      xmax_orig <- xmax_orig - rtemp * (ext@xmax - ext@xmin)
    }

    temp <- (atz - min(zlim)) / diff(range(zlim))
    xs_orig <- xmin_orig + temp * (xmax_orig - xmin_orig)

  } else {
    ymin_orig <- ext@ymin
    if (col_desc[["added_below"]]) {
      ymin_orig <- ymin_orig + rtemp * (ext@ymax - ext@ymin)
    }

    ymax_orig <- ext@ymax
    if (col_desc[["added_above"]]) {
      ymax_orig <- ymax_orig - rtemp * (ext@ymax - ext@ymin)
    }

    temp <- (atz - min(zlim)) / diff(range(zlim))
    ys_orig <- ymin_orig + temp * (ymax_orig - ymin_orig)
  }

  # Draw ticks
  if (tick) {
    lwd_seg <- max(0.5, min(1, cex)) * par("lwd")
    if (horiz) {
      segments(x0 = xs_orig, x1 = xs_orig,
        y0 = ext@ymax - (ext@ymax - ext@ymin) / 3, y1 = ext@ymax, lwd = lwd_seg)

      if (col_desc[["added_below"]]) {
        segments(x0 = ext@xmin, x1 = ext@xmin, y0 = ext@ymin, y1 = ext@ymax,
          lwd = lwd_seg)
      }
      if (col_desc[["added_above"]]) {
        segments(x0 = ext@xmax, x1 = ext@xmax, y0 = ext@ymin, y1 = ext@ymax,
          lwd = lwd_seg)
      }

    } else {
      segments(x0 = ext@xmax - (ext@xmax - ext@xmin) / 3, x1 = ext@xmax,
        y0 = ys_orig, y1 = ys_orig, lwd = lwd_seg)

      if (col_desc[["added_below"]]) {
        segments(x0 = ext@xmin, x1 = ext@xmax, y0 = ext@ymin, y1 = ext@ymin,
          lwd = lwd_seg)
      }
      if (col_desc[["added_above"]]) {
        segments(x0 = ext@xmin, x1 = ext@xmax, y0 = ext@ymax, y1 = ext@ymax,
          lwd = lwd_seg)
      }
    }
  }

  # Decide whether to add tickmark labels
  if (is.logical(labels)) {
    add_automatic_labels <- labels
    add_arg_labels <- FALSE
  } else {
    add_automatic_labels <- FALSE
    labels <- grDevices::as.graphicsAnnot(labels)
    add_arg_labels <- length(labels) > 0 && length(labels) == length(atz)
  }
  add_labels <- add_automatic_labels || add_arg_labels

  # Tickmark labels
  if (add_labels) {
    if (add_automatic_labels) {
      # Create automatic tickmark labels
      ltxt <- prettyNum(signif(signed * fun_inv_ens(atz), 2))
      temp <- signif(signed * fun_inv_ens(zextreme), 2)
      ltext_extreme <- prettyNum(temp)
      if (any(temp < 0) && any(temp > 0)) {
        id <- temp > 0
        ltext_extreme[id] <- paste0("+", ltext_extreme[id])
      }
    } else {
      # generate argument-based tickmark labels
      ltxt <- labels
    }

    # Add tickmark labels to plot
    if (horiz) {
      if (abs(srt) > 45) {
        adj <- c(0.5, 1.3)
        ly <- ext@ymax
      } else {
        adj <- c(0.5, NA)
        ly <- ext@ymax + strheight(ltxt, units = "user", cex = cex * 1.05)
      }

      text(x = xs_orig, y = ly, labels = ltxt, srt = srt, adj = adj, cex = cex,
        xpd = TRUE)

      if (col_desc[["added_below"]] && !add_arg_labels) {
        text(x = ext@xmin, y = ext@ymin, labels = ltext_extreme[1], srt = 90,
          adj = c(0.1, -0.5), cex = cex, xpd = TRUE)
      }
      if (col_desc[["added_above"]] && !add_arg_labels) {
        text(x = ext@xmax, y = ext@ymin, labels = ltext_extreme[2], srt = 90,
          adj = c(0.1, 1.5), cex = cex, xpd = TRUE)
      }

    } else {
      if (abs(srt) > 45) {
        adj <- c(0.5, 1.3)
        lx <- ext@xmax
      } else {
        if (add_automatic_labels) {
          adj <- c(1, NA)
          temp1 <- if (ext@xmax > 0) -1 else +1
          temp2 <- if (cex < 0.5) 1.5 else 1.05
          lx <- ext@xmax + temp1 * max(strwidth(ltxt, units = "user",
            cex = cex * temp2))
        } else {
          adj <- c(-0.05, NA)
          lx <- ext@xmax
        }
      }

      text(x = lx, y = ys_orig, labels = ltxt, srt = srt, adj = adj, cex = cex,
        xpd = TRUE)

      if (col_desc[["added_below"]] && !add_arg_labels) {
        text(x = ext@xmax, y = ext@ymin, labels = ltext_extreme[1], srt = 0,
          adj = c(1, 1.3), cex = cex, xpd = TRUE)
      }
      if (col_desc[["added_above"]] && !add_arg_labels) {
        text(x = ext@xmax, y = ext@ymax, labels = ltext_extreme[2], srt = 0,
          adj = c(1, -0.5), cex = cex, xpd = TRUE)
      }
    }
  }

  invisible(TRUE)
}


#' Extend range by fact
extend_range <- function(x, fact = 1, na.rm = FALSE) {
  x + (fact - 1) * (x - mean(x, na.rm = na.rm))
}




#' Calculate areal extent of cells for each site
#'
#' @param sp_sites A \code{\link[sp:SpatialPoints-class]{sp::SpatialPoints}}
#'   object for which areal extent of cells will be computed.
#' @param sim_raster A \code{\link[raster:Raster-class]{raster::Raster}}
#'   object used for coordinate system and cells. The coordinate system must
#'   use angular coordinates (longitude/latitude).
#' @param tol A numerical value to inform the site to cell matching.
#'
#' @seealso \code{\link[geosphere::areaPolygon]{areaPolygon}}
#'
#' @return A \code{data.frame} with four columns: \var{\dQuote{Longitude}},
#'   \var{\dQuote{Latitude}}, \var{\dQuote{km2}}, and \var{\dQuote{rel}};
#'   and one row for each site of \code{sp_sites}.
#'   Cell area on ellipsoid, based on \code{sim_raster}, in square-kilometer
#'   \code{km2} and as fraction of maximal cell area on the equator \code{rel}.
#' @export
calculate_cell_area <- function(sp_sites, sim_raster,
  tol = sqrt(.Machine$double.eps)) {

  m2_to_km2 <- 1e-6

  temp <- sp::coordinates(sp_sites)
  colnames(temp) <- c("Longitude", "Latitude")

  cells <- data.frame(temp, km2 = NA, rel = NA)


  if (raster::isLonLat(sim_raster)) {
    # Use function `areaPolygon` which works on
    # angular coordinates (longitude/latitude) on an ellipsoid
    stopifnot(requireNamespace("geosphere"))

    unique_lats <- unique(cells[, "Latitude"])

    # Create empty raster except for cells at specified latitudes
    rtemp_init <- raster::init(sim_raster, fun = function(x) rep(NA, x))
    rtemp <- rtemp_init
    xy <- cbind(rep(0, length(unique_lats)), unique_lats)
    rtemp[raster::cellFromXY(rtemp, xy)] <- 1

    # Convert cells into polygons
    ptemp <- raster::rasterToPolygons(rtemp, dissolve = FALSE)

    # Calculate area of polyon for each cell
    for (k in seq_along(ptemp)) {
      icols <- abs(cells[, "Latitude"] - sp::coordinates(ptemp[k, ])[, 2]) < tol
      # Return value of `areaPolygon` is square-meter
      cells[icols, "km2"] <- m2_to_km2 * geosphere::areaPolygon(ptemp[k, ])
    }

    # Calculate area of maximal polyon for a cell on the equator
    rtemp <- rtemp_init
    rtemp[raster::cellFromXY(rtemp, matrix(c(0, 0), nrow = 1))] <- 1
    ptemp <- raster::rasterToPolygons(rtemp, dissolve = FALSE)
    cell_maxarea_km2 <- m2_to_km2 * geosphere::areaPolygon(ptemp)

  } else {
    # Use Euclidean area for projected/flat grids
    ar <- prod(raster::res(sim_raster))

    # Determine distance units: meters or kilometers?
    crs_txt <- raster::crs(sim_raster, asText = TRUE)
    temp <- lapply(strsplit(crs_txt, split = " ")[[1]], strsplit, split = "=")
    temp <- unlist(temp)
    crs_units <- temp[1 + which(temp == "+units")]

    ar_km2 <- ar * switch(crs_units, m = m2_to_km2, km2 = 1, NA)
    cells[, "km2"] <- ar_km2

    cell_maxarea_km2 <- ar_km2
  }

  cells[, "rel"] <- cells[, "km2"] / cell_maxarea_km2

  if (FALSE) {
    # Visualize cell area by latitude
    with(cells, plot(Latitude, km2))
    with(cells, plot(Latitude, rel))
    with(cells, image(MASS::kde2d(Latitude, km2)))

    # Shen, S. S. P. 2017. R Programming for Climate Data Analysis and
    # Visualization: Computing and plotting for NOAA data applications. 1st
    # revised edition. San Diego State University, San Diego, CA.
    # equation 2.2 ('area-weighted average'):
    aw_shen <- cos(cells[, "Latitude"] * pi / 180)

    ids <- match(unique_lats, cells[, "Latitude"])
    plot(cells[ids, "rel"], aw_shen[ids], pch = 46)
    abline(0, 1, col = "red", lty = 2)
    summary(lm(cells[, "rel"] ~ aw_shen, subset = ids))
    # --> Shen (2017) assumes a spherical earth which ever so slightly
    # underestimates cell areas at mid latitudes compared to cell areas on an
    # ellipsoid as used here
    diffs <- cells[, "rel"] - aw_shen
    plot(cell_maxarea_km2 * diffs[ids] ~ abs(cells[ids, "Latitude"]), pch = 46,
      xlab = "abs(Latitude)",
      ylab = "Cell area difference\n(ellipsoid - sphere; km2)")
    summary(lm(diffs ~ abs(cells[, "Latitude"]), subset = ids))
  }

  cells
}


#' Calculate "nominal resolution" of grid
#' @references CMIP6 Global Attributes, DRS, Filenames, Directory Structure, and CV’s
#'   10 September 2018 (v6.2.7)
#'   Appendix 2: Algorithms for Defining the "nominal_resolution" Attribute
#'   https://docs.google.com/document/d/1h0r8RZr_f3-8egBMMh7aqLwy3snpD6_MrDz1q8n5XUk/edit#bookmark=id.ibeh7ad2gpdi
calculate_nominal_resolution <- function(grid, sites, cell_areas_km2) {
  stopifnot(requireNamespace("geosphere"))
  # For a land surface model calculated on its own grid, include all land grid cells

  # For each grid cell, calculate the distance (in km) between each pair of
  # cell vertices and select the maximum distance ("dmax").
  # For latxlon grid cells, for example, dmax would be the diagonal distance.
  res <- raster::res(grid)

  if (raster::isLonLat(grid)) {
    xy <- sp::coordinates(sites)
    xy_lowerleft <- xy - res / 2
    xy_upperright <- xy + res / 2

    id_use <-
      xy_lowerleft[, 1] >= -180 & xy_lowerleft[, 1] <= 180 &
      xy_lowerleft[, 2] >= -90 & xy_lowerleft[, 2] <= 90 &
      xy_upperright[, 1] >= -180 & xy_upperright[, 1] <= 180 &
      xy_upperright[, 2] >= -90 & xy_upperright[, 2] <= 90

    dmax_km <- 1e-3 *
      geosphere::distGeo(xy_lowerleft[id_use, ], xy_upperright[id_use, ])

    # Calculate the mean over all cells of dmax, weighting each by the grid-cell's
    # area (A)
    mean_resolution_km <- weighted.mean(dmax_km, cell_areas_km2[id_use])

  } else {

    mean_resolution_km <- dmax_km <- 1e-3 * sqrt(sum(res ^ 2))
  }



  # Nominal resolution
  ifelse(mean_resolution_km < 0.72, "0.5 km",
    ifelse(mean_resolution_km < 1.6, "1 km",
      ifelse(mean_resolution_km < 3.6, "2.5 km",
        ifelse(mean_resolution_km < 7.2, "5 km",
          ifelse(mean_resolution_km < 16, "10 km",
            ifelse(mean_resolution_km < 36, "25 km",
              ifelse(mean_resolution_km < 72, "50 km",
                ifelse(mean_resolution_km < 160, "100 km",
                  ifelse(mean_resolution_km < 360, "250 km",
                    ifelse(mean_resolution_km < 720, "500 km",
                      ifelse(mean_resolution_km < 1600, "1000 km",
                        ifelse(mean_resolution_km < 3600, "2500 km",
                          ifelse(mean_resolution_km < 7200, "5000 km",
                            "10000 km")))))))))))))
}



#' Find the most common element of x
majority <- function(x, na.rm = FALSE) {
  ina <- is.na(x)
  if (sum(ina) >= length(x)) {
    res <- NA
  } else {
    useMode <- FALSE
    naVal <- NULL
    if (na.rm) {
      x <- na.exclude(x)
      useMode <- TRUE
    }
    if (!na.rm && sum(ina) > 0 && inherits(x, "numeric")) {
      naVal <- min(x, na.rm = TRUE) - 100
      x[ina] <- naVal
      useMode <- TRUE
    }
    if (useMode && requireNamespace("statip")) {
      res <- statip::mfv(x)
    } else {
      temp <- table(x, useNA = if (na.rm) "no" else "always")
      res <- names(temp)[which(c(temp) == max(c(temp)))]
    }

    if (length(res) > 1) res <- res[sample(x = length(res), size = 1)]
    if (useMode && !is.null(res) && isTRUE(res == naVal)) res <- NA
  }

  if (is.factor(x)) {
    factor(res, levels = levels(x), ordered = is.ordered(x))
  } else {
    as.vector(res, mode = typeof(x))
  }
}


setGeneric("span", function(x, na.rm = FALSE, ...) standardGeneric("span"))

setMethod("span", signature(x = "numeric"),
  function(x, na.rm = FALSE, ...) {
    res <- diff(range(x, na.rm = na.rm, ...))
    if (is.infinite(res)) {
      dots <- list(na.rm = na.rm, ...)
      if (isTRUE(dots[["na.rm"]]) || isTRUE(dots[["finite"]])) {
        res <- NA
      }
    }
    res
  })

setMethod("span", signature(x = "character"),
  function(x, na.rm = FALSE, ...) {
    if (na.rm) {
      x <- na.exclude(x)
    } else {
      if (anyNA(x)) return(NA_integer_)
    }

    length(unique(x, ...))
  })

setMethod("span", signature(x = "factor"),
  function(x, na.rm = FALSE, ...) {
    if (na.rm) {
      x <- factor(x, exclude = NA)
    } else {
      if (anyNA(x) && !anyNA(levels(x))) return(NA_integer_)
    }

    nlevels(x)
  })


setGeneric("agreement", function(x, val = NULL, ...) standardGeneric("agreement"))

setMethod("agreement", signature(x = "numeric"),
  function(x, val = NULL, ...) {
    if (is.null(val)) {
      val <- median(x, ...)

    } else {
      if (length(val) > 1) {
        val <- median(val, ...)
      }
    }

    if (is.na(val) || all(is.na(x))) {
      NA
    } else {
      sum(sign(x) == sign(val), ...)
    }
  })

setMethod("agreement", signature(x = "character"),
  function(x, val = NULL, ...) {
    if (is.null(val)) {
      val <- majority(x, ...)

    } else {
      if (length(val) > 1) {
        val <- majority(val, ...)
      }
    }

    if (is.na(val) || all(is.na(x))) {
      NA
    } else {
      sum(x == val, ...)
    }
  })

setMethod("agreement", signature(x = "factor"),
  selectMethod("agreement", signature = "character"))



quantiles_areaweighted <- function(x, area, probs, na.rm = FALSE) {
  stopifnot(requireNamespace("Hmisc"))
  # non-random reliability weights: normwt = TRUE
  Hmisc::wtd.quantile(x, weights = area, probs = probs, normwt = TRUE,
    na.rm = na.rm)
}

fun_central <- function(method = c("median", "mean"), area = NULL,
  na.rm = FALSE) {

  method <- match.arg(method)

  if (!is.null(area) && isTRUE(var(area) > 0)) {
    switch(method,
      median = function(x) quantiles_areaweighted(x, area = area, probs = 0.5,
        na.rm = na.rm),
      mean = function(x) stats::weighted.mean(x, w = area, na.rm = na.rm))

  } else {
    switch(method,
      median = function(x) median(x, na.rm = na.rm),
      mean = function(x) mean(x, na.rm = na.rm))
  }
}


define_direction_1rel <- function() {
  data.frame(
    name = c("Increase", "No change", "Decrease"),
    rx_vs_tol = c(">", NA, "<"),
    stringsAsFactors = FALSE
  )
}


define_direction_2deltas <- function(expected_dy_lt_dx = TRUE) {
  # Define categories
  cl <- data.frame(
      name = c("No change",
        "Smaller decrease", "Reverse to increase", "Larger increase",
        "Larger decrease", "Reverse to decrease", "Smaller increase"),
      dy_vs_dx = c(NA, ">", ">", ">", "<", "<", "<"),
      dx_vs_0 = c(NA, "<", "<", ">", "<", ">", ">"),
      dy_vs_0 = c(NA, "<", ">", ">", "<", "<", ">"),
      stringsAsFactors = FALSE
    )

  temp <- cl[, "dy_vs_dx"] == "<"
  cl[, "expected"] <- if (expected_dy_lt_dx) temp else !temp

  new_levels <- c("More severe", "Less severe")
  new_labels <- new_levels[1 + as.integer(cl[, "expected"])]
  new_labels[is.na(new_labels)] <- "No change"
  cl[, "name3"] <- new_labels

  cl
}

#' Convert data to factor according to
#' \code{define_direction_2deltas()[, "name"]}
factor_6directions <- function(data, cl = NULL) {
  if (is.null(cl)) {
    cl <- define_direction_2deltas()
  }

  if (is.factor(data)) {
    stopifnot(cl[, "name"] == levels(data))
    data <- data
  } else if (is.integer(data)) {
    stopifnot(na.exclude(data) %in% seq_along(cl[, "name"]))
    data <- factor(data, levels = seq_along(cl[, "name"]),
      labels = cl[, "name"])
  } else if (is.character(data)) {
    stopifnot(na.exclude(unique(data)) %in% cl[, "name"])
    data <- factor(data, levels = cl[, "name"])
  } else {
    stop("Data error")
  }

  data
}

#' Remap from six categories of directional modification to two categories
#'
#' @param data Integer, factor, or character vector coded according to
#'   \code{define_direction_2deltas()[, "name"]}
#' @return An object like \code{data} but coded according to
#'   \code{define_direction_2deltas()[, "expected"]}
remap_directions_6to2 <- function(data, expected_dy_lt_dx = TRUE, cl = NULL) {
  if (is.na(expected_dy_lt_dx)) {
    return(rep(NA, length(data)))
  }

  if (is.null(cl)) {
    cl <- define_direction_2deltas(expected_dy_lt_dx)
  }

  data <- factor_6directions(data, cl = cl)


  factor(data, levels = cl[, "name"], labels = cl[, "name3"])
}

compare_direction_2deltas <- function(dx, dy, vars = NULL, cl = NULL,
  tol = NULL) {

  if (length(tol) != 1 || isFALSE(is.finite(tol))) {
    tol <- sqrt(.Machine$double.eps)
  }
  stopifnot(tol < 1 && tol >= 0)

  # Identify variables to compare
  if (is.null(vars)) {
    vars <- intersect(colnames(dx), colnames(dy))
  }
  stopifnot(vars %in% colnames(dx), vars %in% colnames(dy))

  if (is.null(cl)) {
    cl <- define_direction_2deltas()
  }
  ncls <- seq_len(nrow(cl))

  cl[, "ftol"] <- ifelse(cl[, "dy_vs_dx"] == "<", -1, 1)

  # Prepare output container
  n <- nrow(dx)
  res <- data.frame(array(NA, dim = c(n, length(vars)),
    dimnames = list(NULL, vars)))

  # Identify categories
  for (iv in seq_along(vars)) {
    dx2 <- dx[, vars[iv]]
    dy2 <- dy[, vars[iv]]

    tmp <- rep(NA, n)
    id_nochange <- which(is.na(cl[, "dy_vs_dx"]))
    ids_hasdata <- is.finite(dx2) & is.finite(dy2)

    for (k in ncls[-id_nochange]) {
      tol1 <- 1 + cl[k, "ftol"] * sign(dx2) * tol
      ids <-
        match.fun(cl[k, "dy_vs_dx"])(dy2, tol1 * dx2) &
        match.fun(cl[k, "dx_vs_0"])(dx2, 0) &
        match.fun(cl[k, "dy_vs_0"])(dy2, 0)

      tmp[ids_hasdata & ids] <- k
    }

    tmp[ids_hasdata & is.na(tmp)] <- id_nochange

    res[, vars[iv]] <- factor(tmp, levels = ncls, labels = cl[, "name"])
  }

  res
}

#' @param dx Reference data
compare_direction_1rel2ref <- function(dx, dy, vars = NULL, cl = NULL,
  tol = NULL, subset = NULL) {

  n <- nrow(dx)

  if (length(tol) != 1 || isFALSE(is.finite(tol))) {
    tol <- sqrt(.Machine$double.eps)
  }
  stopifnot(tol < 1 && tol >= 0)

  if (is.null(subset)) {
    subset <- rep(TRUE, n)
  }

  # Identify variables to compare
  if (is.null(vars)) {
    vars <- intersect(colnames(dx), colnames(dy))
  }
  stopifnot(vars %in% colnames(dx), vars %in% colnames(dy))

  if (is.null(cl)) {
    cl <- define_direction_1rel()
  }
  ncls <- seq_len(nrow(cl))


  # Prepare output container
  res <- array(NA, dim = c(n, length(vars)), dimnames = list(NULL, vars))

  # Identify categories
  delta <- dy - dx
  rel <- sweep(delta, MARGIN = 1:2, STATS = dx, FUN = "/")

  ids_has_vals <- !is.na(dx) & !is.na(dy)

  tol2 <- sqrt(.Machine$double.eps)
  ids_ref_zero <- abs(dx) < tol2

  # no change: smaller absolute change than tol OR (in case rel is infinite:)
  #            both ref and res are zero
  ids_no_change <- abs(rel) <= tol | (ids_ref_zero & abs(dy) < tol2)
  res[ids_has_vals & ids_no_change] <- cl[is.na(cl[, "rx_vs_tol"]), "name"]

  # increase: larger change than tol OR (in case rel is infinite:)
  #           ref is zero and res is positive
  ids_inc <- rel > tol | (ids_ref_zero & dy > 0)
  res[ids_has_vals & !ids_no_change & ids_inc] <-
    cl[cl[, "rx_vs_tol"] %in% ">", "name"]

  # decrease: smaller change than -tol OR (in case rel is infinite:)
  #           ref is zero and res is negative
  ids_dec <- rel < -tol | (ids_ref_zero & dy < 0)
  res[ids_has_vals & !ids_no_change & ids_dec] <-
    cl[cl[, "rx_vs_tol"] %in% "<", "name"]

  stopifnot(!anyNA(res[subset, ]))

  res
}


find_scen_part <- function(n_elem_offset, sim_scens, meta_part = NULL) {
  temp <- strsplit(sim_scens, split = ".", fixed = TRUE)
  # old rSOILWAT2 wrapper had 3 nelems; new version (rSFSW2) with multiple
  # future time periods has 4 nelems
  nelems <- max(lengths(temp))
  temp <- na.exclude(sapply(temp, function(x) x[nelems - n_elem_offset]))

  if (is.null(meta_part)) {
    temp
  } else {
    meta_part[meta_part %in% temp]
  }
}

find_reqDS <- function(sim_scens, meta = NULL) {
  find_scen_part(n_elem_offset = 3L, sim_scens,
    meta_part = if (!is.null(meta)) meta[["sim_scens"]][["method_DS"]])
}

find_reqDeltaYR <- function(sim_scens, meta = NULL) {
  find_scen_part(n_elem_offset = 2L, sim_scens,
    meta_part = if (!is.null(meta))
      unique(meta[["sim_scens"]][["DeltaStr_yrs"]]))
}


find_reqCSs <- function(sim_scens, meta = NULL) {
  find_scen_part(n_elem_offset = 1L, sim_scens,
    meta_part = if (!is.null(meta)) meta[["sim_scens"]][["reqCSs"]])
}

find_reqMs <- function(sim_scens, meta = NULL) {
  find_scen_part(n_elem_offset = 0L, sim_scens,
    meta_part = if (!is.null(meta)) meta[["sim_scens"]][["reqMs"]])
}

#' @references Convert Latitude/Longitude to UTM [https://www.wavemetrics.com/code-snippet/convert-latitudelongitude-utm] (attributed to Chuck Gantz).
get_UTM_Zone <- function(longitude, latitude) {
  Long <- mean(longitude)
  Lat <- mean(latitude)

  # Make sure longitude is between -180.00 .. 179.9
  LongTemp <- Long - floor((Long + 180) / 360) * 360

  ZoneNumber <- floor((LongTemp + 180) / 6) + 1;

  if (Lat >= 56 && Lat < 64 && LongTemp >= 3 && LongTemp < 12) {
    ZoneNumber <- 32
  }

  # Special zones for Svalbard
  if (Lat >= 72 && Lat < 84) {
    if (LongTemp >= 0 && LongTemp < 9) {
      ZoneNumber <- 31
    } else if (LongTemp >= 9 && LongTemp < 21) {
      ZoneNumber <- 33
    } else if (LongTemp >= 21 && LongTemp < 33) {
      ZoneNumber <- 35
    } else if (LongTemp >= 33 && LongTemp < 42) {
      ZoneNumber <- 37
    }
  }

  ZoneNumber
}



#' Define (adjusted) analysis regions based on IPCC 2012 SREX
SREX2012_regions <- function(adjusted = FALSE) {
  slist_regions <- list()

  # IPCC (2012) Managing the Risks of Extreme Events and Disasters to Advance Climate
  # Change Adaptation. A Special Report of Working Groups I and II of the
  # Intergovernmental Panel on Climate Change  [C. B. Field, V. Baros, T. F. Stocker,
  # D. Qin, D. J. Dokken, K. L. Ebi, M. D. Mastrandrea, K .J. Mach, G.-K. Plattner,
  # S. K. Allen, M. Tignor and P. M. Midgley (eds.)], Cambridge University Press,
  # Cambridge, United Kingdom, and New York, NY, USA.

  # Figure 3-1 and Table 3.A-1

  #--- North America
  slist <- list()
  #   3 WNA = (28.566N, 105.000W) (28.566N, 130.000W) (60.000N, 130.000W) (60.000N, 105.000W),
  slist[["3_WNA"]] <- sp::Polygon(matrix(c(
    -105.000, 28.566,
    -130.000, 28.566,
    -130.000, 60.000,
    -105.000, 60.000,
    -105.000, 28.566
  ), ncol = 2, byrow = TRUE))
  #   4 CNA = (50.000N, 85.000W) (28.566N, 85.000W) (28.566N, 105.000W) (50.000N, 105.000W),
  slist[["4_CNA"]] <- sp::Polygon(matrix(c(
    -85.000, 50.000,
    -85.000, 28.566,
    -105.000, 28.566,
    -105.000, 50.000,
    -85.000, 50.000
  ), ncol = 2, byrow = TRUE))
  #   6 CAM = (11.439N, 68.800W) (1.239S, 79.729W) (28.566N, 118.323W) (28.566N, 90.315W)
  slist[["6_CAM"]] <- sp::Polygon(matrix(c(
    -68.800, 11.439,
    -79.729, -1.239,
    -118.323, 28.566,
    -90.315, 28.566,
    -68.800, 11.439
  ), ncol = 2, byrow = TRUE))

  if (adjusted) {
    # ADJUSTEMENT: add Caribbean islands to North America
    slist[["NA_Adj"]] <- sp::Polygon(matrix(c(
      -85.000 - 1, 25.000,
      -60.000, 25.000,
      -60.000, 11.439,
      -68.800, 11.439,
      -85.000 - 1, 25.000
    ), ncol = 2, byrow = TRUE))
  }

  slist_regions[["North America"]] <- sp::Polygons(slist, ID = "North America")


  #--- South America
  slist <- list()
  #   7 AMZ = (20.000S, 66.377W) (1.239S, 79.729W) (11.439N, 68.800W) (11.439N, 50.000W) (20.000S, 50.000W),
  slist[["7_AMZ"]] <- sp::Polygon(matrix(c(
    -66.377, -20.000,
    -79.729, -1.239,
    -68.800, 11.439,
    -50.000, 11.439,
    -50.000, -20.000,
    -66.377, -20.000
  ), ncol = 2, byrow = TRUE))
  #   8 NEB = (20.000S, 34.000W) (20.000S, 50.000W) (0.000N, 50.000W) (0.000N, 34.000W),
  slist[["8_NEB"]] <- sp::Polygon(matrix(c(
    -34.000, -20.000,
    -50.000, -20.000,
    -50.000, 0.000,
    -34.000, 0.000,
    -34.000, -20.000
  ), ncol = 2, byrow = TRUE))
  #   9 WSA = (1.239S, 79.729W) (20.000S, 66.377W) (50.000S, 72.141W) (56.704S, 67.348W) (56.704S, 82.022W) (0.530N, 82.022W),
  slist[["9_WSA"]] <- sp::Polygon(matrix(c(
    -79.729, -1.239,
    -66.377, -20.000,
    -72.141, -50.000,
    -67.348, -56.704,
    -82.022, -56.704,
    -82.022, 0.530,
    -79.729, -1.239
  ), ncol = 2, byrow = TRUE))
  #  10 SSA = (20.000S, 39.376W) (56.704S, 39.376W) (56.704S, 67.348W) (50.000S, 72.141W) (20.000S, 66.377W)
  slist[["10_SSA"]] <- sp::Polygon(matrix(c(
    -39.376, -20.000,
    -39.376, -56.704,
    -67.348, -56.704,
    -72.141, -50.000,
    -66.377, -20.000,
    -39.376, -20.000
  ), ncol = 2, byrow = TRUE))
  if (adjusted) {
    # ADJUSTEMENT: add Galapogos islands to South America
    slist[["SA_Adj"]] <- sp::Polygon(matrix(c(
      -82.022, 1.000,
      -92.000, 1.000,
      -92.000, -2.000,
      -82.022, -2.000,
      -82.022, 1.000
    ), ncol = 2, byrow = TRUE))
  }

  slist_regions[["South America"]] <- sp::Polygons(slist, ID = "South America")


  #--- Southern Africa
  slist <- list()
  #   SREX region 17 SAF = (35.000S, 10.000W) (11.365S, 10.000W) (11.365S, 51.990E) (35.000S, 51.990E)
  slist[["17_SAF"]] <- sp::Polygon(matrix(c(
    -10.000, -35.000,
    -10.000, -11.365,
    51.990, -11.365,
    51.990, -35.000,
    -10.000, -35.000
  ), ncol = 2, byrow = TRUE))

  slist_regions[["Southern Africa"]] <- sp::Polygons(slist, ID = "Southern Africa")


  #--- Sahara and subsaharan Africa
  slist <- list()
  #   14 SAH = (15.000N, 20.000W) (30.000N, 20.000W) (30.000N, 40.000E) (15.000N, 40.000E),
  slist[["14_SAH"]] <- sp::Polygon(matrix(c(
    # ADJUSTEMENT: add all of Marocco to Mediterranean instead of Subsaharan Africa
    -20.000, 15.000,
    -20.000, if (adjusted) 22.000 else 30.000,
    40.000, if (adjusted) 22.000 else 30.000,
    40.000, 15.000,
    -20.000, 15.000
  ), ncol = 2, byrow = TRUE))
  #   15 WAF = (11.365S, 20.000W) (15.000N, 20.000W) (15.000N, 25.000E) (11.365S, 25.000E),
  slist[["15_WAF"]] <- sp::Polygon(matrix(c(
    -20.000, -11.365,
    -20.000, 15.000,
    25.000, 15.000,
    25.000, -11.365,
    -20.000, -11.365
  ), ncol = 2, byrow = TRUE))
  #   16 EAF = (11.365S, 25.000E) (15.000N, 25.000E) (15.000N, 51.990E) (11.365S, 51.990E)
  slist[["16_EAF"]] <- sp::Polygon(matrix(c(
    # ADJUSTEMENT: add Yemen to Western Asia instead of Subsaharan Africa
    25.000, -11.365,
    25.000, 15.000,
    if (adjusted) c(40.000, 15.000),
    if (adjusted) c(40.000, 17.150),
    if (adjusted) c(44.500, 11.150),
    51.990, 15.000,
    51.990, -11.365,
    25.000, -11.365
  ), ncol = 2, byrow = TRUE))

  slist_regions[["Subsaharan Africa"]] <- sp::Polygons(slist, ID = "Subsaharan Africa")


  #--- Mediterranean Basin
  slist <- list()
  #   SREX region 13 MED = (30.000N, 10.000W) (45.000N, 10.000W) (45.000N, 40.000E) (30.000N, 40.000E)
  slist[["13_MED"]] <- sp::Polygon(matrix(c(
    # ADJUSTEMENT: add all of Marocco to Mediterranean instead of Subsaharan Africa
    -10.000, if (adjusted) 22.000 else 30.000,
    -10.000, 45.000,
    40.000, 45.000,
    40.000, if (adjusted) 22.000 else 30.000,
    -10.000, if (adjusted) 22.000 else 30.000
  ), ncol = 2, byrow = TRUE))
  if (adjusted) {
    # ADJUSTEMENT: add Azores, Madeira, and Canary Islands to Mediterranean Basin
    slist[["MB_Adj1"]] <- sp::Polygon(matrix(c(
      -32.000, 22.000,
      -32.000, 45.000,
      -10.000, 45.000,
      -10.000, 22.000,
      -32.000, 22.000
    ), ncol = 2, byrow = TRUE))
    # ADJUSTEMENT: add areas just adjacent to the North to the Mediterranean Basin
    slist[["MB_Adj2"]] <- sp::Polygon(matrix(c(
      -10.000, 45.000,
      -10.000, 50.000,
      40.000, 50.000,
      40.000, 45.000,
      -10.000, 45.000
    ), ncol = 2, byrow = TRUE))
  }

  slist_regions[["Mediterranean"]] <- sp::Polygons(slist, ID = "Mediterranean")


  #--- Western Asia
  slist <- list()
  #   SREX region 19 WAS = (15.000N, 40.000E) (50.000N, 40.000E) (50.000N, 60.000E) (15.000N, 60.000E)
  slist[["19_WAS"]] <- sp::Polygon(matrix(c(
    # ADJUSTEMENT: add Yemen to Western Asia instead of Subsaharan Africa
    40.000, if (adjusted) 17.150 else 15.000,
    40.000, 50.000,
    60.000, 50.000,
    60.000, 15.000,
    if (adjusted) c(51.990, 15.000),
    if (adjusted) c(44.500, 11.150),
    40.000, if (adjusted) 17.150 else 15.000
  ), ncol = 2, byrow = TRUE))
  if (adjusted) {
    # ADJUSTEMENT: add areas just adjacent to the North to Western Asia
    slist[["WA_Adj"]] <- sp::Polygon(matrix(c(
      40.000, 50.000,
      40.000, 55.000,
      60.000, 55.000,
      60.000, 50.000,
      40.000, 50.000
    ), ncol = 2, byrow = TRUE))
  }

  slist_regions[["Western Asia"]] <- sp::Polygons(slist, ID = "Western Asia")


  #--- Central Asia
  slist <- list()
  #   20 CAS = (30.000N, 60.000E) (50.000N, 60.000E) (50.000N, 75.000E) (30.000N, 75.000E)
  slist[["20_CAS"]] <- sp::Polygon(matrix(c(
    60.000, 30.000,
    60.000, 50.000,
    75.000, 50.000,
    75.000, 30.000,
    60.000, 30.000
  ), ncol = 2, byrow = TRUE))
  #   21 TIB = (30.000N, 75.000E) (50.000N, 75.000E) (50.000N, 100.000E) (30.000N, 100.000E)
  slist[["21_TIB"]] <- sp::Polygon(matrix(c(
    75.000, 30.000,
    75.000, 50.000,
    100.000, 50.000,
    100.000, 30.000,
    75.000, 30.000
  ), ncol = 2, byrow = TRUE))
  if (adjusted) {
    # ADJUSTEMENT: add areas just adjacent to the North to Central Asia
    slist[["CA_Adj"]] <- sp::Polygon(matrix(c(
      60.000, 50.000,
      60.000, 55.000,
      100.000, 55.000,
      100.000, 50.000,
      60.000, 50.000
    ), ncol = 2, byrow = TRUE))
  }

  slist_regions[["Central Asia"]] <- sp::Polygons(slist, ID = "Central Asia")


  #--- East Asia
  slist <- list()
  #   22 EAS = (20.000N, 100.000E) (50.000N, 100.000E) (50.000N, 145.000E) (20.000N, 145.000E)
  slist[["22_EAS"]] <- sp::Polygon(matrix(c(
    100.000, 20.000,
    100.000, 50.000,
    145.000, 50.000,
    145.000, 20.000,
    100.000, 20.000
  ), ncol = 2, byrow = TRUE))
  if (adjusted) {
    # ADJUSTEMENT: add areas just adjacent to the North to East Asia
    slist[["EA_Adj"]] <- sp::Polygon(matrix(c(
      100.000, 50.000,
      100.000, 55.000,
      145.000, 55.000,
      145.000, 50.000,
      100.000, 50.000
    ), ncol = 2, byrow = TRUE))
  }

  slist_regions[["East Asia"]] <- sp::Polygons(slist, ID = "East Asia")


  #--- South Asia
  slist <- list()
  #   23 SAS = (5.000N, 60.000E) (30.000N, 60.000E) (30.000N, 100.000E)
  #            (20.000N, 100.000E) (20.000N, 95.000E) (5.000N, 95.000E)
  slist[["23_SAS"]] <- sp::Polygon(matrix(c(
    60.000, 5.000,
    60.000, 30.000,
    100.000, 30.000,
    100.000, 20.000,
    95.000, 20.000,
    95.000, 5.000,
    60.000, 5.000
  ), ncol = 2, byrow = TRUE))

  slist_regions[["South Asia"]] <- sp::Polygons(slist, ID = "South Asia")


  #--- Australia = SREX regions 25 NAU, 26 SAU
  slist <- list()
  #   25 NAU = (30.000S, 110.000E) (10.000S, 110.000E) (10.000S, 155.000E) (30.000S, 155.000E)
  slist[["25_NAU"]] <- sp::Polygon(matrix(c(
    110.000, -30.000,
    110.000, -10.000,
    155.000, -10.000,
    155.000, -30.000,
    110.000, -30.000
  ), ncol = 2, byrow = TRUE))
  #   26 SAU = (50.000S, 110.000E) (30.000S, 110.000E) (30.000S, 180.000E) (50.000S, 180.000E)
  slist[["26_SAU"]] <- sp::Polygon(matrix(c(
    110.000, -50.000,
    110.000, -30.000,
    180.000, -30.000,
    180.000, -50.000,
    110.000, -50.000
  ), ncol = 2, byrow = TRUE))

  slist_regions[["Australia"]] <- sp::Polygons(slist, ID = "Australia")

  #--- Convert to SpatialPolygons
  wgs84 <- sp::CRS("+init=epsg:4326")
  spoly_regions <- sp::SpatialPolygons(slist_regions, proj4string = wgs84)

  region_names <- names(spoly_regions)
  spoly_regions2 <- NULL

  for (k in seq_along(spoly_regions)) {
    # Correct for 'orphaned holes'
    temp <- spoly_regions[k, ]
    slot(temp, "polygons") <- lapply(slot(temp, "polygons"), maptools::checkPolygonsHoles)

    # Merge individual pieces
    temp0 <- raster::aggregate(temp)
    temp0 <- sp::spChFIDs(temp0, region_names[k])

    # Put pieces back together into a SpatialPolygons object
    spoly_regions2 <- if (k > 1) {
        maptools::spRbind(spoly_regions2, temp0)
      } else {
        temp0
      }
  }

  spoly_regions2
}


ADA_SREX2012_regions <- function(meta) {
  sp_regions_orig <- SREX2012_regions(adjusted = FALSE)
  sp_regions_adj <- SREX2012_regions(adjusted = TRUE)

  temp <- sp::spTransform(meta[["sim_space"]][["run_sites"]],
    CRSobj = sp::proj4string(sp_regions_orig))
  data_orig <- sp::over(temp, sp_regions_orig)
  data_adj <- sp::over(temp, sp_regions_adj)

  list(
    original = list(
      region_data = factor(names(sp_regions_orig)[unname(data_orig)]),
      spoly = sp_regions_orig),
    adjusted = list(
      region_data = factor(names(sp_regions_adj)[unname(data_adj)]),
      spoly = sp_regions_adj)
  )
}


plot_densities <- function(meta, data, variables, var_labels, var_colors,
  xlim, region_categories, region_data, deltas = TRUE, dir_res_out, fname_tag) {

  reqCSs <- meta[["sim_scens"]][["reqCSs"]]
  reqGs <- if (deltas) grep("G", reqCSs, value = TRUE) else NULL
  reqCSs2 <- c(reqCSs, reqGs)

  sc_current <- meta[["sim_scens"]][["ambient"]]
  deltas_byScenario <- if (deltas) {
      c(rep(sc_current, length(reqCSs)), rep("RCP45", length(reqGs)))
    } else NULL

  n_panels <- c(length(reqCSs2), length(variables))
  yt <- seq(0, 1, length.out = 100L)
  var_colors <- rep_len(var_colors, length(variables))

  for (kr in seq_along(region_categories)) {
    id_region <- if (region_categories[kr] == "Global") {
        rep(TRUE, length(region_data))
      } else {
        region_data %in% region_categories[kr]
      }

    path_fig <- file.path(dir_res_out, "Figures",
      paste0("Densities_", if (deltas) "Changes_", "GCMs_by_Scenarios_", fname_tag, "_",
      region_categories[kr], ".png"))

    h.panel <- 40 / 25.4; h.edgeL <- 0.25; h.edgeU <- 0.25 # Heigth of panel and lower/upper edge
    w.panel <- 60 / 25.4; w.edgeL <- 0.25; w.edgeR <- 0.05 # Width of panel and left/right edge

    png(
      height = h.edgeL + h.panel * n_panels[1L] + h.edgeU,
      width = w.edgeL + w.panel * n_panels[2L] + w.edgeR,
      units = "in", res = 600L, file = path_fig)

    temp <- rep(0L, (2L + n_panels[1L]) * (2L + n_panels[2L]))
    lmat <- matrix(temp, nrow = 2L + n_panels[1L], ncol = 2L + n_panels[2L], byrow = TRUE)

    stemp2 <- seq_len(n_panels[2L])
    for (k in seq_len(n_panels[1L])) { # byrow = TRUE
      lmat[k + 1L, -c(1L, 2L + n_panels[2L])] <- (k - 1L) * n_panels[2L] + stemp2
    }

    layout(lmat,
      heights = c(h.edgeU, rep(h.panel, n_panels[1L]), h.edgeL),
      widths = c(w.edgeL, rep(w.panel, n_panels[2L]), w.edgeR))

    par_prev <- par(mar = c(1, 1, 1, 0.1), mgp = c(1, 0, 0), tcl = 0.3)

    for (sc in seq_along(reqCSs2)) {
      sc_name <- reqCSs2[sc]
      isc <- grep(sc_name, meta[["sim_scens"]][["id"]])

      if (deltas) {
        ref <- deltas_byScenario[sc]
        x <- data[["data"]][[paste0("wCO2_delta", ref)]]
        xno <- data[["data"]][[paste0("noCO2_delta", ref)]]
        xens <- data[["dataEns"]][[paste0("wCO2_delta", ref)]]
        xnoens <- data[["dataEns"]][[paste0("noCO2_delta", ref)]]
      } else {
        ref <- ""
        x <- data[["data"]][["wCO2"]]
        xno <- data[["data"]][["noCO2"]]
        xens <- data[["dataEns"]][["wCO2"]]
        xnoens <- data[["dataEns"]][["noCO2"]]
      }

      for (iv in seq_along(variables)) {
        vlab <- paste("delta", var_labels[iv])
        print(paste(Sys.time(), region_categories[kr], sc_name, ref, variables[iv],
          sep = "-"))

        if (length(isc) > 0) {
          xlim_used <- if (is.null(xlim)) {
              quantile(c(x[, variables[iv], ], xno[, variables[iv], ],
                xens["mean", , variables[iv], ], xnoens["mean", , variables[iv], ]),
                probs = c(0.025, 0.975), na.rm = TRUE)
            } else xlim

          plot(NA, type = "n", xlim = xlim_used, ylim = c(0, 1),
            xlab = "", ylab = if (iv == 1) "Fn(x)" else "", main = "", xpd = NA,
            col = var_colors[iv])
          temp <- if (deltas) paste(sc_name, "to", ref) else sc_name
          mtext(side = 3, text = temp, font = 2, cex = par("cex"),
            xpd = NA)

          if (sc == length(reqCSs2)) {
            temp <- gregexpr("\n", vlab)[[1L]]
            nlines <- 1L + sum(temp > 0)
            mtext(side = 1, line = nlines, text = vlab, cex = par("cex"), xpd = NA)
          }

          for (k in seq_along(isc)) {
            temp <- x[id_region, variables[iv], isc[k]]
            if (!all(is.na(temp))) {
              lines(ecdf(temp), do.points = FALSE, col = var_colors[iv])
            }

            temp <- xno[id_region, variables[iv], isc[k]]
            if (!all(is.na(temp))) {
              lines(ecdf(temp), do.points = FALSE, col = "gray", lty = 2)
            }
          }

          temp <- xens["mean", id_region, variables[iv], sc_name]
          if (!all(is.na(temp))) {
            dem <- ecdf(temp)
            lines(dem, do.points = FALSE, col = var_colors[iv], lwd = 2)
            if (FALSE) {
              # it seems that the lty argument is broken for the 'stepfun' methods of plot and lines
              lines(dem, do.points = FALSE, col = adjustcolor("black", alpha.f = 0.7),
                lty = 2, lwd = 2)
            } else {
              lines(x = quantile(dem, probs = yt), y = yt,
                col = adjustcolor("black", alpha.f = 0.7), lty = 2, lwd = 2)
            }
          }

          temp <- xnoens["mean", id_region, variables[iv], sc_name]
          if (!all(is.na(temp))) {
            demno <- ecdf(temp)
            lines(demno, do.points = FALSE, col = "gray", lwd = 2)
            if (FALSE) {
              # it seems that the lty argument is broken for the 'stepfun' methods of plot and lines
              lines(demno, do.points = FALSE, col = adjustcolor("black", alpha.f = 0.7),
                lty = 4, lwd = 2)
            } else {
              lines(x = quantile(demno, probs = yt), y = yt,
                col = adjustcolor("black", alpha.f = 0.7), lty = 4, lwd = 2)
            }
          }

          abline(v = 0, lty = 2)
          abline(h = 0.5, lty = 2)

        } else {
          plot.new()
        }
      }
    }

    par(par_prev)
    dev.off()
  }
}

#' Tabulate values by region and by climate zone
#'
#' @param is_x_per_area A logical value. If \code{TRUE}, then \code{x} values
#'   are interpreted to represent rates per unit area, i.e,. they
#'   are weighted by \code{area} (if \code{fmethod} is \code{quantile}
#'   or \code{mean}) or multiplied (if \code{fmethod} is \code{sum}). If
#'   \code{FALSE}, then \code{x} values are aggregated as is.
tabulate_by_region_x_climzone <- function(x, x_regional, vars, subset,
  fmethod = c("quantile", "mean", "sum"), is_x_per_area = TRUE, area,
  probs = c(0, 0.05, 0.5, 0.95, 1),
  region_data, name_regions = NULL, climzones_data, name_climzones = NULL) {

  # Organize data: global
  if (isTRUE(length(dim(x)) == 2)) {
    vars <- vars[vars %in% colnames(x)]
    x_glo <- x[, vars, drop = FALSE]

  } else {
    vars <- vars[1]
    x_glo <- matrix(x, nrow = length(x), ncol = 1,
      dimnames = list(NULL, vars))
  }

  # Organize data: regional if present; if not present, use `x`
  if (!missing(x_regional)) {
    if (isTRUE(length(dim(x_regional)) == 2)) {
      vars <- vars[vars %in% colnames(x_regional)]
      x_reg <- x_regional[, vars, drop = FALSE]

    } else {
      vars <- vars[1]
      x_reg <- matrix(x_regional, nrow = length(x_regional), ncol = 1,
        dimnames = list(NULL, vars))
    }

  } else {
    x_reg <- x_glo
  }


  fmethod <- match.arg(fmethod)
  if (fmethod %in% c("mean", "sum")) {
    probs <- -1
  }

  # Region and climate zone subsets
  if (is.null(name_regions)) {
    name_regions <- c("Global", sort(as.character(unique(region_data))))
  }
  if (is.null(name_climzones)) {
    name_climzones <- c("Global", sort(as.character(unique(climzones_data))))
  }

  #--- Prepare table
  subsets_rc <- expand.grid(Region = name_regions, ClimateZone = name_climzones,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  N <- nrow(subsets_rc)

  temp <- matrix(NA, nrow = N, ncol = length(vars) * length(probs))
  if (all(probs >= 0)) {
    colnames(temp) <- paste0(rep(vars, each = length(probs)), "_q",
      formatC(round(100 * probs), width = 3, flag = "0"))
  } else {
    colnames(temp) <- paste0(vars, "_", fmethod)
  }
  res <- cbind(subsets_rc, temp, stringsAsFactors = FALSE)

  icols_vars <- unlist(lapply(vars, function(x) grep(x, colnames(res))))


  #--- Prepare subsets
  id_all <- subset

  ids_region_subsets <- lapply(name_regions, function(region) {
    if (region == "Global") {
      id_all
    } else {
      id_all & (region_data %in% region)
    }
  })
  names(ids_region_subsets) <- name_regions

  ids_climzone_subsets <- lapply(name_climzones, function(cz) {
    if (cz == "Global") {
      id_all
    } else {
      id_all & (climzones_data %in% cz)
    }
  })
  names(ids_climzone_subsets) <- name_climzones


  #--- Summarize data into table
  for (k in seq_len(N)) {
    # Identify data slice
    id_use <- ids_region_subsets[[res[k, "Region"]]] &
      ids_climzone_subsets[[res[k, "ClimateZone"]]]

    # Identify global/regional dataset
    xdat <- if (res[k, "Region"] == "Global") x_glo else x_reg

    # Aggregate for each variable
    temp <- if (fmethod == "quantile") {
        apply(xdat[id_use, vars, drop = FALSE], 2, quantiles_areaweighted,
          area = area[id_use], probs = probs, na.rm = TRUE)
      } else {
        apply(xdat[id_use, vars, drop = FALSE], 2, match.fun(fmethod),
          na.rm = TRUE)
      }

    res[k, icols_vars] <- as.vector(temp)
  }

  res
}


add_panel_map <- function(x, meta, subset, ...) {
  # Prepare plotting arguments
  dots <- list(...)
  dots_names <- names(dots)

  id_dots_NAs <- match(c("useNA", "valNA"), dots_names, nomatch = 0)
  has_dots_NAs <- id_dots_NAs > 0

  if ("ylim" %in% dots_names && !("asp" %in% dots_names)) {
    # Adjust aspect ratio to correct for longitude/latitude
    dots[["asp"]] <- 1 / cos((mean(dots[["ylim"]]) * pi) / 180)
  }

  # Create raster from data
  if (!missing(subset)) {
    x[!subset] <- NA
    if (all(has_dots_NAs)) {
      x[subset & is.na(x)] <- dots[["valNA"]]
    }
  }
  rdata <- create_raster_from_variables(meta, data = x)

  # Clean up argument list
  if (any(has_dots_NAs)) {
    dots <- dots[-id_dots_NAs]
  }

  # Plot map
  do.call(selectMethod("image", class(rdata)), args = c(list(x = rdata), dots))

  if (requireNamespace("maps")) {
    maps::map(add = TRUE, resolution = 0, lwd = 0.5) # Country border
  }
  abline(h = 0, col = "darkgray", lty = 2, lwd = 1) # Equator

  invisible(rdata)
}


define_value_colorramp <- function(zlim, Nc = 255,
  cols = c("gray", "gold", "red")) {

  cramp <- colorRampPalette(cols)(Nc)

  list(
    panel = list(cols = cramp, zlim = zlim),
    legend = list(cols = cramp, zlim = zlim,
      desc = list(
        colors_label = cramp,
        ncol_label_added = 0,
        added_below = FALSE,
        added_above = FALSE))
  )
}

define_delta_colorramp <- function(zlim, zlim_trim = NULL, Nc = 255,
  use_rev_cols = FALSE, add_extremes_to_legend = FALSE) {

  col_zero <- "gray"
  cols_blues <- c("royalblue", "cornflowerblue", "lightblue")
  col_blue_extreme <- adjustcolor(cols_blues[1],
    green.f = 0.5, red.f = 0.5, blue.f = 0.5)
  cols_reds <- c("lightsalmon", "mediumorchid1", "purple", "purple4")
  col_red_extreme <- adjustcolor(cols_reds[length(cols_reds)],
    green.f = 0.5, red.f = 0.5, blue.f = 0.5)

  colramp <- if (use_rev_cols) {
      list(
        neg = colorRampPalette(rev(cols_reds)),
        neg_extreme = col_red_extreme,
        pos = colorRampPalette(rev(cols_blues)),
        pos_extreme = col_blue_extreme)
    } else {
      list(
        neg = colorRampPalette(cols_blues),
        neg_extreme = col_blue_extreme,
        pos = colorRampPalette(cols_reds),
        pos_extreme = col_red_extreme)
    }

  zlim_legend <- if (is.null(zlim_trim)) zlim else zlim_trim
  Nc_per_unit <- Nc / sum(abs(zlim_legend))

  if (is.finite(Nc_per_unit)) {
    ntemp <- abs(round(Nc_per_unit * zlim_legend))
    ctemp1 <- colramp[["neg"]](ntemp[1])
    ctemp2 <- colramp[["pos"]](ntemp[2])
    ctemp0 <- c(ctemp1, ctemp2)
    ctemp0[ntemp[1] + c(0, 1)] <- col_zero

    if (add_extremes_to_legend) {
      temp <- ceiling(0.05 * sum(ntemp))
      ncol_label_added <- if (temp > 0) max(5L, temp) else 0L
      cramp_legend <- c(
        rep(colramp[["neg_extreme"]], ncol_label_added),
        ctemp0,
        rep(colramp[["pos_extreme"]], ncol_label_added))

    } else {
      ncol_label_added <- 0
      cramp_legend <- ctemp0
    }

    ntemp2 <- abs(round(Nc_per_unit * (zlim - zlim_legend)))
    cramp <- c(
      rep(colramp[["neg_extreme"]], ntemp2[1]),
      ctemp0,
      rep(colramp[["pos_extreme"]], ntemp2[2]))

  } else {
    ncol_label_added <- 0
    cramp <- cramp_legend <- NA
  }

  list(
    panel = list(cols = cramp, zlim = zlim),
    legend = list(cols = cramp_legend, zlim = zlim_legend,
      desc = list(
        colors_label = cramp_legend,
        ncol_label_added = ncol_label_added,
        added_below = add_extremes_to_legend,
        added_above = add_extremes_to_legend))
  )
}


plot_map_studyarea <- function(meta, climzones_data, subset,
  add_SREX2012 = FALSE, SREX2012 = NULL, dtemp) {

  x <- factor(climzones_data)
  cols <- c("royalblue", "skyblue", "coral3", "orange", "darkolivegreen4")
  xlim <- c(-130, 160)
  ylim <- c(-60, 55)

  if (add_SREX2012) {
    spoly_labels <- unname(sapply(names(SREX2012[["adjusted"]][["spoly"]]),
      function(x) {
        x <- strsplit(x, split = " ", fixed = TRUE)[[1]]
        if (length(x) > 1) {
          paste(sapply(x, function(x) substr(x, 1, 2)), collapse = "")
        } else {
          substr(x, 1, 3)
        }
      }))
    spoly_label_pos <- t(sapply(seq_along(SREX2012[["adjusted"]][["spoly"]]),
      function(k) {
        temp <- slot(SREX2012[["adjusted"]][["spoly"]][k, ], "bbox")
        if (abs(temp["y", "max"] >= abs(temp["y", "min"]))) {
          # mainly north of equator
          temp <- c(temp["x", "min"], temp["y", "max"], -0.05, 1.1)
        } else {
          temp <- c(temp["x", "min"], temp["y", "min"], -0.05, -0.3)
        }
        temp <- c(min(xlim[2], max(xlim[1], temp[1])),
          min(ylim[2], max(ylim[1], temp[2])),
          temp[3:4])
      }))
    tag <- "SREX2012&"

  } else {
    tag <- ""
  }

  fexp <- 1.5
  w.panel <- fexp * 160 / 25.4
  h.panel <- sum(abs(ylim)) / sum(abs(xlim)) * w.panel

  pdf(height = h.panel, width = w.panel,
    file = file.path(dtemp, paste0("Fig_StudyArea_", tag, "ClimateZones.pdf")))
  par_prev <- par(mar = c(1.25, 1.25, 0.25, 0.25), mgp = c(1, 0, 0), tcl = 0.3,
    cex = 0.75 * fexp)

  add_panel_map(x = x, meta = meta, subset = subset,
    col = cols, xlim = xlim, ylim = ylim)

  if (add_SREX2012) {
    spoly <- SREX2012[["adjusted"]][["spoly"]]
    selectMethod("lines", class(spoly))(spoly, lwd = 1,
      col = adjustcolor("black", alpha.f = 0.5))

    for (k in seq_len(nrow(spoly_label_pos))) {
      text(x = spoly_label_pos[k, 1], y = spoly_label_pos[k, 2],
        adj = spoly_label_pos[k, 3:4],
        labels = spoly_labels[k], cex = 0.5, font = 2)
    }
  }

  legend("bottomleft", ncol = 1, cex = 0.85, inset = c(0.01, 0.025),
    legend = levels(x), fill = cols)

  par(par_prev)
  dev.off()
}



plot_map_vars_dGi_and_modifier <- function(x, xmod, meta, subset, GCMs, exp,
  ref0, ref1, var_labels, mod_with_all = TRUE, mod_tol = NA,
  expected_dy_lt_dx = TRUE, panel_indep = FALSE, path, ftag) {

  funs <- c("vals", "modifier")
  temp <- dimnames(x)
  vars <- temp[[2]]
  stopifnot(dim(x) == dim(xmod), all.equal(temp, dimnames(xmod)),
    sapply(GCMs, function(x) any(grepl(x, temp[[3]]))))

  n_vars <- length(vars)
  n_GCMs <- length(GCMs)
  n_funs <- length(funs)

  expected_dy_lt_dx <- if (is.logical(expected_dy_lt_dx)) {
      rep_len(expected_dy_lt_dx, n_vars)
    } else {
      expected_dy_lt_dx > 0
    }
  stopifnot(length(expected_dy_lt_dx) == n_vars)

  # Legend labels
  cl <- define_direction_2deltas(TRUE)
  lab_small_mod <- if (is.na(mod_tol)) {
      "No modification"
    } else {
      paste0("<", round(100 * mod_tol), "% modification")
    }
  legend_labs_mod <- if (mod_with_all) {
      temp <- cl[, "name"]
      temp[temp == "No change"] <- lab_small_mod
      temp
    } else {
      temp <- unique(cl[, "name3"])
      temp[temp == "No change"] <- lab_small_mod
      temp
    }

  # Re-order legend
  ibp_labs <- if (mod_with_all) {
      seq_along(legend_labs_mod)
    } else {
      match(rev(c("More severe", lab_small_mod, "Less severe")), legend_labs_mod)
    }

  # Colors
  col_funs <- list(
    median = NA,
    modifier = if (mod_with_all) {
        list(n = 7,
          cols = c("gray",
            sd = "lightblue", ri = "purple", li = "darkred",
            ld = "darkblue", rd = "mediumseagreen", si = "orange"))
      } else {
        list(n = 3,
          cols = c("gray", "purple", "aquamarine2"))
      }
  )


  xlim <- c(-130, 160)
  ylim <- c(-60, 55)

  n_panels <- c(n_vars, n_funs)

  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  # Figure size
  w.panel <- 160 / 25.4
  w.edgeL <- 0.1; w.edgeI <- 0; w.edgeR <- 0.05 # Width of panel and left/interior/right edge
  h.panel <- sum(abs(ylim)) / sum(abs(xlim)) * w.panel
  h.edgeL <- 0.1; h.edgeI <- if (panel_indep) 0.1 else 0; h.edgeU <- 0.1 # Heigth of panel and lower/interior/upper edge

  # Figure layout
  lmat <- matrix(0L,
    nrow = 2 + 2 * n_panels[1] - 1L, ncol = 2 + 2 * n_panels[2] - 1,
    byrow = TRUE)

  stemp2 <- seq_len(n_panels[2L])
  for (k in seq_len(n_panels[1L])) { # byrow = TRUE
    #lmat[k + 1L, -c(1L, 2L + n_panels[2L])] <- (k - 1L) * n_panels[2L] + stemp2
    lmat[(k - 1) * 2 + 2, seq_len(ncol(lmat)) %% 2 == 0] <- (k - 1L) * n_panels[2L] + stemp2
  }

  temp <- rep(c(h.panel, h.edgeI), n_panels[1L])
  layout_heights <- c(h.edgeU, temp[-length(temp)], h.edgeL)
  temp <- rep(c(w.panel, w.edgeI), n_panels[2L])
  layout_widths <- c(w.edgeL, temp[-length(temp)], w.edgeR)


  # Figure device
  fexp <- 1.5

  for (igcm in seq_len(n_GCMs)) {
    i <- 1
    fname <- file.path(path,
      paste0("Fig_Map_dGiandMods_", ftag, "_", GCMs[igcm], ".png"))

    png(units = "in", res = 150,
      height = fexp * (h.edgeL + h.panel * n_panels[1] + h.edgeI * (n_panels[1] - 1) + h.edgeU),
      width = fexp * (w.edgeL + w.panel * n_panels[2] + w.edgeI * (n_panels[2] - 1) + w.edgeR),
      file = fname)

    layout(lmat, heights = layout_heights, widths = layout_widths)

    mar <- if (panel_indep) c(0.5, 0.75, 0.75, 0.5) else rep(0.5, 4)
    par_prev <- par(mar = mar, mgp = c(1, 0, 0), tcl = 0.3,
      cex = fexp)

    for (iv in seq_len(n_vars)) {
      for (k in seq_len(n_funs)) {
        if (funs[k] == "modifier") {
          data <- xmod[, vars[iv], igcm]

          legend_title <- paste0(exp, " modifies d(", cur_to_hist(ref1), "):")
          legend_horiz <- FALSE
          legend_tick <- FALSE
          legend_adj <- c(0, NA)
          legend_pos <- NULL
          legend_cex <- 0.65
          if (mod_with_all) {
            legend_title_xy <- matrix(c(-135, 2), nrow = 1)
            legend_box <- c(-135, -129, -55, 0)
          } else {
            legend_title_xy <- matrix(c(-135, -8), nrow = 1)
            legend_box <- c(-135, -129, -40, -10)
          }

        } else {
          data <- x[, vars[iv], igcm]

          legend_title <- paste0("d(", exp, ", ", cur_to_hist(ref0), "): ",
            var_labels[iv])
          legend_horiz <- TRUE
          legend_box <- c(-35, 95, -55, -49)
          legend_tick <- TRUE
          legend_title_xy <- matrix(c(35, -47), nrow = 1)
          legend_cex <- 0.75
          legend_adj <- NULL
          legend_pos <- 3
        }

        # Specify legend
        ftemp <- col_funs[[funs[k]]]

        if (funs[k] %in% c("agreement", "modifier")) {
          ctemp <- define_value_colorramp(zlim = c(1, ftemp[["n"]]),
            Nc = ftemp[["n"]], cols = ftemp[["cols"]])
          ctemp[["legend"]][["zlim"]][2] <- ftemp[["n"]] + 1

          legend_ats <- seq_len(ftemp[["n"]]) + 0.5
          useNA <- TRUE

          if (funs[k] == "modifier") {
            data <- factor_6directions(data, cl = cl)
            if (!mod_with_all) {
              data <- remap_directions_6to2(data, expected_dy_lt_dx[iv])
            }
            data <- as.integer(data)

            # Re-order legend labels
            legend_labs <- legend_labs_mod[ibp_labs]
            ctemp[["legend"]][["desc"]][["colors_label"]] <-
              ctemp[["legend"]][["desc"]][["colors_label"]][ibp_labs]

          } else {
            legend_labs <- c("", seq_len(n_GCMs))
          }

        } else {
          legend_ats <- NULL
          legend_labs <- TRUE
          useNA <- FALSE
          qtemp <- quantile(data, probs = c(0, 0.01, 0.99, 1), na.rm = TRUE)
          #qtemp <- quantile(data, probs = c(0, 0, 1, 1), na.rm = TRUE)

          if (funs[k] %in% "span") {
            ctemp <- define_value_colorramp(zlim = qtemp[c(1, 4)],
              #zlim_trim = qtemp[c(2, 3)],
              Nc = ftemp[["n"]], cols = ftemp[["cols"]])

          } else {
            ctemp <- define_delta_colorramp(zlim = qtemp[c(1, 4)],
              zlim_trim = qtemp[c(2, 3)],
              add_extremes_to_legend = diff(qtemp[c(1, 4)]) >
                1.1 * diff(qtemp[c(2, 3)]))
          }
        }

        # Plot map to panel
        rdata <- add_panel_map(x = data, meta = meta, subset = subset,
          xlim = xlim, ylim = ylim, zlim = ctemp[["panel"]][["zlim"]],
          col = ctemp[["panel"]][["cols"]], useNA = useNA, valNA = 0,
          axes = FALSE, ann = FALSE)

        box()
        axis(side = 1, labels = panel_indep || k %in% c(3, 4))
        axis(side = 2, labels = panel_indep || k %in% c(1, 3))

        # Add legend to panel
        do_extremes <- ctemp[["legend"]][["desc"]][["added_below"]] ||
          ctemp[["legend"]][["desc"]][["added_above"]]

        add_legend(
          zlim = ctemp[["legend"]][["zlim"]], zextreme =
            if (do_extremes) ctemp[["panel"]][["zlim"]] else NULL,
          col_desc = ctemp[["legend"]][["desc"]], grid = rdata,
          box = legend_box, horiz = legend_horiz,
          at = legend_ats, tick = legend_tick, labels = legend_labs,
          srt = 0, cex = legend_cex)
        text(x = legend_title_xy, labels = legend_title, pos = legend_pos,
          offset = 0.75, adj = legend_adj, cex = legend_cex)


        # Add panel identifier
        mtext(side = 3, adj = 0, text = letters[i], font = 2, cex = fexp)
        i <- i + 1
      }
    }

    par(par_prev)
    dev.off()
  }
}


#' @param funs Is either of length two and consists of
#'   c("vals", "delta", "modifier")
plot_map_vars_dGi_and_modifier_v2 <- function(x1, x2,
  funs = c("delta", "modifier"), meta, subset, GCMs, exp,
  ref0, ref1, var_labels, mod_with_all = TRUE, mod_tol = NA,
  expected_dy_lt_dx = TRUE, panel_indep = FALSE, path, ftag) {

  temp <- dimnames(x1)
  vars <- temp[[2]]
  stopifnot(dim(x1) == dim(x2), all.equal(temp[-1], dimnames(x2)[-1]),
    sapply(GCMs, function(x1) any(grepl(x1, temp[[3]]))))

  x <- list(x1, x2)

  n_vars <- length(vars)
  n_GCMs <- length(GCMs)
  n_funs <- length(funs)

  expected_dy_lt_dx <- if (is.logical(expected_dy_lt_dx)) {
    rep_len(expected_dy_lt_dx, n_vars)
  } else {
    expected_dy_lt_dx > 0
  }
  stopifnot(length(expected_dy_lt_dx) == n_vars)

  # Legend labels
  cl <- define_direction_2deltas(TRUE)
  lab_small_mod <- if (is.na(mod_tol)) {
    "No modification"
  } else {
    paste0("<", round(100 * mod_tol), "% modification")
  }
  legend_labs_mod <- if (mod_with_all) {
    temp <- cl[, "name"]
    temp[temp == "No change"] <- lab_small_mod
    temp
  } else {
    temp <- unique(cl[, "name3"])
    temp[temp == "No change"] <- lab_small_mod
    temp
  }

  # Re-order legend
  ibp_labs <- if (mod_with_all) {
    seq_along(legend_labs_mod)
  } else {
    match(rev(c("More severe", lab_small_mod, "Less severe")), legend_labs_mod)
  }

  # Colors
  col_funs <- list(
    median = NA,
    modifier = if (mod_with_all) {
      list(n = 7,
        cols = c("gray",
          sd = "lightblue", ri = "purple", li = "darkred",
          ld = "darkblue", rd = "mediumseagreen", si = "orange"))
    } else {
      list(n = 3,
        cols = c("gray", "purple", "aquamarine2"))
    }
  )


  xlim <- c(-130, 160)
  ylim <- c(-60, 55)

  n_panels <- c(n_vars, n_funs)

  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  # Figure size
  w.panel <- 160 / 25.4
  w.edgeL <- 0.1; w.edgeI <- 0; w.edgeR <- 0.05 # Width of panel and left/interior/right edge
  h.panel <- sum(abs(ylim)) / sum(abs(xlim)) * w.panel
  h.edgeL <- 0.1; h.edgeI <- if (panel_indep) 0.1 else 0; h.edgeU <- 0.1 # Heigth of panel and lower/interior/upper edge

  # Figure layout
  lmat <- matrix(0L,
    nrow = 2 + 2 * n_panels[1] - 1L, ncol = 2 + 2 * n_panels[2] - 1,
    byrow = TRUE)

  stemp2 <- seq_len(n_panels[2L])
  for (k in seq_len(n_panels[1L])) { # byrow = TRUE
    #lmat[k + 1L, -c(1L, 2L + n_panels[2L])] <- (k - 1L) * n_panels[2L] + stemp2
    lmat[(k - 1) * 2 + 2, seq_len(ncol(lmat)) %% 2 == 0] <- (k - 1L) * n_panels[2L] + stemp2
  }

  temp <- rep(c(h.panel, h.edgeI), n_panels[1L])
  layout_heights <- c(h.edgeU, temp[-length(temp)], h.edgeL)
  temp <- rep(c(w.panel, w.edgeI), n_panels[2L])
  layout_widths <- c(w.edgeL, temp[-length(temp)], w.edgeR)


  # Figure device
  fexp <- 1.5

  for (igcm in seq_len(n_GCMs)) {
    i <- 1
    fname <- file.path(path,
      paste0("Fig_Map_dGiandMods_", ftag, "_", GCMs[igcm], ".png"))

    png(units = "in", res = 150,
      height = fexp * (h.edgeL + h.panel * n_panels[1] + h.edgeI * (n_panels[1] - 1) + h.edgeU),
      width = fexp * (w.edgeL + w.panel * n_panels[2] + w.edgeI * (n_panels[2] - 1) + w.edgeR),
      file = fname)

    layout(lmat, heights = layout_heights, widths = layout_widths)

    mar <- if (panel_indep) c(0.5, 0.75, 0.75, 0.5) else rep(0.5, 4)
    par_prev <- par(mar = mar, mgp = c(1, 0, 0), tcl = 0.3,
      cex = fexp)

    for (iv in seq_len(n_vars)) {
      for (k in seq_len(n_funs)) {
        data <- x[[k]][, vars[iv], igcm]

        if (funs[k] == "modifier") {
          legend_title <- paste0(exp, " modifies d(", cur_to_hist(ref1), "):")
          legend_horiz <- FALSE
          legend_tick <- FALSE
          legend_adj <- c(0, NA)
          legend_pos <- NULL
          legend_cex <- 0.65
          if (mod_with_all) {
            legend_title_xy <- matrix(c(-135, 2), nrow = 1)
            legend_box <- c(-135, -129, -55, 0)
          } else {
            legend_title_xy <- matrix(c(-135, -8), nrow = 1)
            legend_box <- c(-135, -129, -40, -10)
          }

        } else {
          legend_title <- if (funs[k] == "delta") {
              paste0("d(", exp, ", ", cur_to_hist(ref0), "): ", var_labels[iv])
            } else if (funs[k] == "vals") {
              paste0(cur_to_hist(ref0), ": ", var_labels[iv])
            }
          legend_horiz <- TRUE
          legend_box <- c(-35, 95, -55, -49)
          legend_tick <- TRUE
          legend_title_xy <- matrix(c(35, -47), nrow = 1)
          legend_cex <- 0.75
          legend_adj <- NULL
          legend_pos <- 3
        }

        # Specify legend
        ftemp <- col_funs[[funs[k]]]

        if (funs[k] %in% c("agreement", "modifier")) {
          ctemp <- define_value_colorramp(zlim = c(1, ftemp[["n"]]),
            Nc = ftemp[["n"]], cols = ftemp[["cols"]])
          ctemp[["legend"]][["zlim"]][2] <- ftemp[["n"]] + 1

          legend_ats <- seq_len(ftemp[["n"]]) + 0.5
          useNA <- TRUE

          if (funs[k] == "modifier") {
            data <- factor_6directions(data, cl = cl)
            if (!mod_with_all) {
              data <- remap_directions_6to2(data, expected_dy_lt_dx[iv])
            }
            data <- as.integer(data)

            # Re-order legend labels
            legend_labs <- legend_labs_mod[ibp_labs]
            ctemp[["legend"]][["desc"]][["colors_label"]] <-
              ctemp[["legend"]][["desc"]][["colors_label"]][ibp_labs]

          } else {
            legend_labs <- c("", seq_len(n_GCMs))
          }

        } else {
          legend_ats <- NULL
          legend_labs <- TRUE
          useNA <- FALSE
          qtemp <- quantile(data, probs = c(0, 0.01, 0.99, 1), na.rm = TRUE)
          #qtemp <- quantile(data, probs = c(0, 0, 1, 1), na.rm = TRUE)

          if (funs[k] %in% "span") {
            ctemp <- define_value_colorramp(zlim = qtemp[c(1, 4)],
              #zlim_trim = qtemp[c(2, 3)],
              Nc = ftemp[["n"]], cols = ftemp[["cols"]])

          } else {
            ctemp <- define_delta_colorramp(zlim = qtemp[c(1, 4)],
              zlim_trim = qtemp[c(2, 3)],
              add_extremes_to_legend = diff(qtemp[c(1, 4)]) >
                1.1 * diff(qtemp[c(2, 3)]))
          }
        }

        # Plot map to panel
        rdata <- add_panel_map(x = data, meta = meta, subset = subset,
          xlim = xlim, ylim = ylim, zlim = ctemp[["panel"]][["zlim"]],
          col = ctemp[["panel"]][["cols"]], useNA = useNA, valNA = 0,
          axes = FALSE, ann = FALSE)

        box()
        axis(side = 1, labels = panel_indep || k %in% c(3, 4))
        axis(side = 2, labels = panel_indep || k %in% c(1, 3))

        # Add legend to panel
        do_extremes <- ctemp[["legend"]][["desc"]][["added_below"]] ||
          ctemp[["legend"]][["desc"]][["added_above"]]

        add_legend(
          zlim = ctemp[["legend"]][["zlim"]], zextreme =
            if (do_extremes) ctemp[["panel"]][["zlim"]] else NULL,
          col_desc = ctemp[["legend"]][["desc"]], grid = rdata,
          box = legend_box, horiz = legend_horiz,
          at = legend_ats, tick = legend_tick, labels = legend_labs,
          srt = 0, cex = legend_cex)
        text(x = legend_title_xy, labels = legend_title, pos = legend_pos,
          offset = 0.75, adj = legend_adj, cex = legend_cex)


        # Add panel identifier
        mtext(side = 3, adj = 0, text = letters[i], font = 2, cex = fexp)
        i <- i + 1
      }
    }

    par(par_prev)
    dev.off()
  }
}



tabulate_vars_dGmodifiers <- function(xmod, method = c("dir", "rdi"),
  vars, var_labs, subset, cell_areas, region_data, name_regions, GCMs, exp,
  mod_with_all = TRUE, mod_tol = NA, expected_dy_lt_dx = TRUE,
  path_table, ftag) {

  method <- match.arg(method)

  if (is.null(name_regions)) {
    name_regions <- c("Global", sort(as.character(unique(region_data))))
  }

  cl <- switch(EXPR = method,
    dir = define_direction_2deltas(TRUE),
    rdi = define_direction_1rel(),
    NULL)

  lab_nochange <- "No change"
  legend_labs <- if (method == "dir") {
      if (mod_with_all) cl[, "name"] else unique(cl[, "name3"])
    } else {
      cl[, "name"]
    }

  temp <- dimnames(xmod)
  stopifnot(identical(vars, temp[[2]]))
  n_vars <- length(vars)

  scens <- temp[[3]]
  all_GCMs <- unique(unlist(GCMs))
  temp <- find_reqCSs(scens)
  reqCSs_per_M <- lapply(all_GCMs,
    function(x) temp[grep(paste0(x, "$"), scens)])
  names(reqCSs_per_M) <- all_GCMs

  dir.create(path_table, recursive = TRUE, showWarnings = FALSE)

  if (method != "rdi" && !is.null(expected_dy_lt_dx)) {
    expected_dy_lt_dx <- if (is.logical(expected_dy_lt_dx)) {
      rep_len(expected_dy_lt_dx, n_vars)
    } else {
      expected_dy_lt_dx > 0
    }
    stopifnot(length(expected_dy_lt_dx) == n_vars)
  }

  #--- Calculate region areas
  darea <- cell_areas[, "km2"]
  darea[!subset] <- NA
  region_areas <- aggregate(darea, by = list(region_data), sum, na.rm = TRUE)
  names(region_areas) <- c("Region", "Area_km2")
  region_areas <- rbind(region_areas,
    data.frame(Region = "Global", Area_km2 = sum(darea, na.rm = TRUE)))


  #--- Prepare table
  res <- NULL
  for (ie in seq_along(exp)) {
    temp <- expand.grid(GCM = GCMs[[ie]], Variable = vars,
      Region = name_regions, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    res <- rbind(res, cbind(Exp = exp[ie], temp))
  }
  N <- nrow(res)
  resX <- matrix(NA, nrow = N, ncol = 1 + length(legend_labs),
    dimnames = list(NULL, c("RegionArea_km2", legend_labs)))

  ids <- match(res[, "Region"], region_areas[, "Region"], nomatch = 0)
  resX[ids > 0, "RegionArea_km2"] <- region_areas[ids, "Area_km2"]

  #--- Calculate table
  for (k in seq_len(N)) {
    iexp <- res[k, "Exp"]
    igcm <- grep(paste(iexp, res[k, "GCM"], sep = "."), scens, value = TRUE)[1]
    ivar <- res[k, "Variable"]
    iused <- if (res[k, "Region"] == "Global") {
        subset
      } else {
        subset & region_data %in% res[k, "Region"]
      }

    data <- xmod[, ivar, igcm]

    if (method == "dir") {
      data <- factor_6directions(data, cl = cl)

      if (!mod_with_all) {
        # Remap
        data <- remap_directions_6to2(data,
          expected_dy_lt_dx[which(vars == res[k, "Variable"])])
      }
    } else if (method == "rdi") {
      data <- factor(data, levels = cl[, "name"])
    }

    xagg <- aggregate(darea[iused], by = list(data[iused]), sum, na.rm = TRUE,
      drop = FALSE)
    xagg[, 1] <- as.character(xagg[, 1])
    if (isTRUE(is.null(unlist(xagg[, 2])))) {
      # Set extent to 0 if no categories exist
      xagg[, 2] <- 0
    }
    if (anyNA(xagg[, 2])) {
      # Set extent to 0 if some categories don't exist
      xagg[is.na(xagg[, 2]), 2] <- 0
    }

    id <- match(xagg[, 1], levels(data), nomatch = 0)
    resX[k, levels(data)[id > 0]] <- xagg[id, 2]

    # Check that parts sum up to region extent
    temp <- sum(resX[k, legend_labs])
    ttemp <- all.equal(temp, resX[k, "RegionArea_km2"], check.attributes = FALSE)
    if (!isTRUE(ttemp)) {
      stop(paste0("`tabulate_vars_dGmodifiers` failed for: ", shQuote(ftag),
        " [row=", k, "; ", res[k, "Region"], "]: categories (",
        paste(shQuote(legend_labs), "=", round(resX[k, legend_labs]), "km2",
          collapse = ", "),
        ") do not sum up (", round(temp), " km2) to region area (",
        round(resX[k, "RegionArea_km2"]), " km2)"))
    }
  }

  fname <- file.path(path_table, paste0("Table_dGiandMods_", ftag, ".csv"))
  write.csv(cbind(res, resX), file = fname, row.names = FALSE)
}



barplot_vars_dGmodifiers <- function(method = c("dir", "rdi"), vars, var_labs,
  name_regions, scens, GCMs, exp, mod_with_all = TRUE, mod_tol = NA,
  path_fig, path_table, ftag) {

  method <- match.arg(method)

  cl <- switch(EXPR = method,
    dir = define_direction_2deltas(TRUE),
    rdi = define_direction_1rel(),
    NULL)

  legend_labs <- if (method == "dir") {
      if (mod_with_all) cl[, "name"] else unique(cl[, "name3"])
    } else {
      cl[, "name"]
    }

  n_vars <- length(vars)

  all_GCMs <- unique(unlist(GCMs))
  temp <- find_reqCSs(scens)
  reqCSs_per_M <- lapply(all_GCMs,
    function(x) temp[grep(paste0(x, "$"), scens)])
  names(reqCSs_per_M) <- all_GCMs

  dir.create(path_fig, recursive = TRUE, showWarnings = FALSE)


  #--- Load tables
  res_icols <- c("Exp", "GCM", "Variable", "Region")
  resX_icols <- c("RegionArea_km2", legend_labs)

  fname <- file.path(path_table, paste0("Table_dGiandMods_", ftag, ".csv"))
  temp <- read.csv(file = fname)

  res <- temp[, res_icols]
  resX <- temp[, make.names(resX_icols)]


  #--- Plot bars
  # Categories
  ibp_noc <- -1 # remove "no/small modification"
  ibp_labs0 <- seq_along(legend_labs[ibp_noc])

  ibp_labs <- if (mod_with_all) {
      tempn <- length(legend_labs) / 2
      temp <- seq_len(tempn)
      1 + c(temp, 0, tempn + temp)
    } else {
      if (method == "dir") {
        match(c("More severe", "No change", "Less severe"), legend_labs)
      } else {
        match(c("Increase", "No change", "Decrease"), legend_labs)
      }
    }

  mod_cols <- if (mod_with_all) {
      c("gray",
        sd = "lightblue", ri = "purple", li = "darkred",
        ld = "darkblue", rd = "mediumseagreen", si = "orange")
    } else {
      if (method == "dir") {
        c(noc = "gray", more = "purple", less = "aquamarine2")
      } else {
        c(inc = "purple", noc = "gray", dec = "aquamarine2")
      }
    }

  # Figure size
  fexp <- 1
  w.panel <- 3
  w.edgeL <- 0.6; w.edgeI <- 0; w.edgeR <- 0 # Width of panel and left/interior/right edge
  h.panel <- 2.5
  h.edgeL <- 1; h.edgeI <- 0.1; h.edgeU <- 0.1 # Heigth of panel and lower/interior/upper edge

  # Barplots for each GCM
  for (igcm in seq_along(all_GCMs)) {
    exp_per_m <- reqCSs_per_M[[all_GCMs[igcm]]]

    i <- 1
    fname <- file.path(path_fig,
      paste0("Fig_Barplot_dGiandMods_", ftag, "_", all_GCMs[igcm], ".pdf"))

    # Figure layout
    n_panels <- c(n_vars, length(exp_per_m))
    lmat <- matrix(0L,
      nrow = 2 + 2 * n_panels[1] - 1, ncol = 2 + 2 * n_panels[2] - 1,
      byrow = FALSE)

    stemp2 <- seq_len(n_panels[2L])
    for (k in seq_len(n_panels[1L])) { # byrow = FALSE
      lmat[(k - 1) * 2 + 2, seq_len(ncol(lmat)) %% 2 == 0] <- (stemp2 - 1L) * n_panels[1L] + k
    }

    # Figure size
    temp <- rep(c(h.panel, h.edgeI), n_panels[1L])
    layout_heights <- c(h.edgeU, temp[-length(temp)], h.edgeL)
    temp <- rep(c(w.panel, w.edgeI), n_panels[2L])
    layout_widths <- c(w.edgeL, temp[-length(temp)], w.edgeR)

    pdf(
      height = fexp * (h.edgeL + h.panel * n_panels[1] + h.edgeI * (n_panels[1] - 1) + h.edgeU),
      width = fexp * (w.edgeL + w.panel * n_panels[2] + w.edgeI * (n_panels[2] - 1) + w.edgeR),
      file = fname)

    layout(lmat, heights = layout_heights, widths = layout_widths)

    par_prev <- par(mar = c(0.5, 0.5, 0.5, 0.1), mgp = c(1, 0, 0), tcl = 0.3,
      cex = fexp)

    for (ie in seq_along(exp_per_m)) {
      for (iv in seq_len(n_vars)) {
        ks <- which(res[, "Exp"] %in% exp_per_m[ie] &
          res[, "GCM"] %in% all_GCMs[igcm] &
          res[, "Variable"] %in% vars[iv])
        stopifnot(length(ks) == length(name_regions))
        ig <- which(res[ks, "Region"] == "Global")

        dat_rarea <- resX[ks, "RegionArea_km2"]
        dat <- resX[ks, -1] / dat_rarea
        wreg <- dat_rarea / max(dat_rarea)

        # add bars for global average
        barplot(matrix(dat[ig, ibp_labs], ncol = 1),
          col = adjustcolor(mod_cols[ibp_labs], alpha.f = 0.4),
          xlim = c(0, 1.2), ylim = c(0, 1), axisnames = FALSE, axes = FALSE)

        # add bars for regions
        btemp <- barplot(t(dat[-ig, ibp_labs]), width = wreg[-ig],
          col = mod_cols[ibp_labs], xlim = c(0, 1.2), ylim = c(0, 1),
          axisnames = FALSE, axes = FALSE, add = TRUE)
        pos0 <- btemp[1] - wreg[-ig][1] / 2

        # add lines for global average
        gdat <- cumsum(dat[ig, ibp_labs, drop = TRUE])
        for (k in ibp_labs0) {
          segments(x0 = pos0, x1 = 1.2, y0 = gdat[k],
            col = adjustcolor(mod_cols[ibp_noc][k],
              offset = c(rep(-0.4, 3), 0)),
            lwd = 2, lend = 1)
        }

        # annotate and add axes
        axis(side = 2, labels = ie == 1)

        if (ie == 1) {
          mtext(side = 2, line = 1, xpd = NA,
            text = paste0( "Modification area (%) in\n", var_labs[iv]))
        }

        if (iv == n_vars) {
          text(x = btemp, y = 0, labels = as.character(res[ks, "Region"][-ig]),
            srt = 45, adj = c(1, 1), xpd = NA)
        }

        # Add panel identifier
        mtext(side = 3, adj = 0.05, text = letters[i], font = 2, cex = fexp)
        i <- i + 1
      }
    }

    par(par_prev)
    dev.off()
  }
}


barplot_vars_dGmodifiers_v3 <- function(method = c("dir", "rdi"),
  vars, var_labs, name_regions, scens, scens_ens, GCMs, GCMs_ens, GCM_target,
  ref, recalc_ens = FALSE, mod_with_all = TRUE, mod_tol = NA,
  path_fig, path_table, ftag1, ftag2, ftag3) {

  method <- match.arg(method)
  ref <- cur_to_hist(ref)

  cl <- switch(EXPR = method,
    dir = define_direction_2deltas(TRUE),
    rdi = define_direction_1rel(),
    NULL)


  if (method == "dir") {
    legend_labs <- if (mod_with_all) cl[, "name"] else unique(cl[, "name3"])
    icol_name_plus <- "More.severe"

  } else {
    legend_labs <- cl[, "name"]
    icol_name_plus <- "Increase"
  }

  icol_name_neg <- "No.change" # because of cumsum!


  n_vars <- length(vars)

  all_GCMs <- unique(unlist(GCMs))
  id_GCM_target <- if (missing(GCM_target)) 0 else which(all_GCMs == GCM_target)

  all_CSs <- find_reqCSs(scens)
  reqCSs_per_M <- lapply(all_GCMs,
    function(x) all_CSs[grep(paste0(x, "$"), scens)])
  names(reqCSs_per_M) <- all_GCMs
  all_CSs <- unique(all_CSs)

  all_GCMs_ens <- unique(unlist(GCMs_ens))
  all_CSs_ens <- find_reqCSs(scens_ens)
  reqCSs_per_M_ens <- lapply(all_GCMs_ens,
    function(x) all_CSs_ens[grep(paste0(x, "$"), scens_ens)])
  names(reqCSs_per_M_ens) <- all_GCMs_ens

  stopifnot(all_CSs_ens %in% all_CSs)

  dir.create(path_fig, recursive = TRUE, showWarnings = FALSE)


  #--- Load tables
  res_icols <- c("Exp", "GCM", "Variable", "Region")
  resX_icols <- c("RegionArea_km2", legend_labs)

  fname <- file.path(path_table, paste0("Table_dGiandMods_", ftag1, ".csv"))
  temp <- read.csv(file = fname)

  res <- temp[, res_icols]
  resX <- temp[, make.names(resX_icols)]

  fname <- file.path(path_table, paste0("Table_dGiandMods_", ftag1, "_", ftag2,
    ".csv"))
  temp <- read.csv(file = fname)

  res_ens <- temp[, res_icols]

  if (recalc_ens) {
    # Re-calculate ensemble for figure-panel consistency
    temp <- res
    temp[, "GCM"] <- "mean"
    temp <- aggregate(resX, by = as.list(temp), mean)

    id <- match(apply(res_ens, 1, paste, collapse = ""),
      table = apply(temp[, colnames(res_ens)], 1, paste, collapse = ""))
    resX_ens <- temp[id, colnames(resX)]

  } else {
    # Use pre-calculated ensemble for overall consistency
    resX_ens <- temp[, make.names(resX_icols)]
  }


  #--- Plot bars
  # Categories
  ibp_noc <- -1 # remove "no/small modification"
  ibp_labs0 <- seq_along(legend_labs[ibp_noc])


  ibp_labs <- if (mod_with_all) {
      tempn <- length(legend_labs) / 2
      temp <- seq_len(tempn)
      1 + c(temp, 0, tempn + temp)
      stop("Function not capable of processing 'mod_with_all'")
    } else {
      if (method == "dir") {
        match(c("More severe", "No change", "Less severe"), legend_labs)
      } else {
        match(c("Increase", "No change", "Decrease"), legend_labs)
      }
    }

  mod_cols <- if (mod_with_all) {
      c("gray",
        sd = "lightblue", ri = "purple", li = "darkred",
        ld = "darkblue", rd = "mediumseagreen", si = "orange")
    } else {
      if (method == "dir") {
        c(noc = "gray", more = "purple", less = "aquamarine2")
      } else {
        c(more = "purple", noc = "gray", less = "aquamarine2")
      }
    }

  tol_lab <- paste0(">", round(100 * if (is.na(mod_tol)) 0 else mod_tol), "%")


  mod_cols_dark1 <- adjustcolor(mod_cols, offset = c(rep(-0.2, 3), 0))
  mod_cols_dark2 <- adjustcolor(mod_cols, offset = c(rep(-0.5, 3), 0))
  names(mod_cols_dark1) <- names(mod_cols_dark2) <- names(mod_cols)

  # Barplots with all GCM
  fname <- file.path(path_fig,
    paste0("Fig_Barplot_dGiandMods_", ftag1, "_", ftag3, "_", ftag2, ".pdf"))

  # Figure size
  fexp <- 1
  w.panel <- 4
  w.edgeL <- 0.6; w.edgeI <- 0; w.edgeR <- 0 # Width of panel and left/interior/right edge
  h.panel <- 2.5
  h.edgeL <- 1; h.edgeI <- 0.1; h.edgeU <- 0.1 # Heigth of panel and lower/interior/upper edge

  # Figure layout
  n_panels <- c(n_vars, length(all_CSs))
  lmat <- matrix(0L,
    nrow = 2 + 2 * n_panels[1] - 1, ncol = 2 + 2 * n_panels[2] - 1,
    byrow = FALSE)

  stemp2 <- seq_len(n_panels[2L])
  for (k in seq_len(n_panels[1L])) { # byrow = FALSE
    lmat[(k - 1) * 2 + 2, seq_len(ncol(lmat)) %% 2 == 0] <- (stemp2 - 1L) * n_panels[1L] + k
  }

  # Figure size
  temp <- rep(c(h.panel, h.edgeI), n_panels[1L])
  layout_heights <- c(h.edgeU, temp[-length(temp)], h.edgeL)
  temp <- rep(c(w.panel, w.edgeI), n_panels[2L])
  layout_widths <- c(w.edgeL, temp[-length(temp)], w.edgeR)

  pdf(
    height = fexp * (h.edgeL + h.panel * n_panels[1] + h.edgeI * (n_panels[1] - 1) + h.edgeU),
    width = fexp * (w.edgeL + w.panel * n_panels[2] + w.edgeI * (n_panels[2] - 1) + w.edgeR),
    file = fname)

  layout(lmat, heights = layout_heights, widths = layout_widths)

  par_prev <- par(mar = c(0.5, 0.5, 0.5, 0.5), mgp = c(1, 0, 0), tcl = 0.3,
    cex = fexp)

  i <- 1 # panel counter
  pos0f <- 0
  pos1f <- 1.35
  wglob <- 0.1
  nids <- 4 # sub-bar elements

  # give more space between global and regions
  ig <- which(name_regions == "Global")
  bspace <- rep(0.2, length(name_regions))
  bspace[ig + 1] <- 0.5 # "amount of space ... left before each bar"

  # columns organized by scenarios:
  for (ie in seq_along(all_CSs)) {
    # rows organized by variables:
    for (iv in seq_len(n_vars)) {

      # panel: based on ensemble
      iens <- 1
      if (!(all_CSs[ie] %in% reqCSs_per_M_ens[[all_GCMs_ens[iens]]])) {
        next
      }

      ks_ens <- which(res_ens[, "Exp"] %in% all_CSs[ie] &
          res_ens[, "GCM"] %in% all_GCMs_ens[iens] &
          res_ens[, "Variable"] %in% vars[iv])
      stopifnot(length(ks_ens) == length(name_regions))

      id_regions <- match(name_regions, res_ens[ks_ens, "Region"])
      ks_ens2 <- ks_ens[id_regions]
      dat_rarea <- resX_ens[ks_ens2, "RegionArea_km2"]
      dat <- resX_ens[ks_ens2, 1 + ibp_labs] / dat_rarea

      wreg <- dat_rarea / max(dat_rarea)
      wreg[ig] <- wglob # fix global bar width

      # add ensemble bars for global average and for regions
      btemp <- barplot(t(dat), width = wreg, space = bspace,
        col = mod_cols[ibp_labs], border = NA,
        xlim = c(pos0f, pos1f), ylim = c(0, 1), axes = FALSE,
        names.arg = NULL, axisnames = FALSE)

      # positions of regional bars
      pos0r <- btemp - wreg / 2
      pos1r <- btemp + wreg / 2
      dposr <- (pos1r - pos0r) / nids

      # highlight bar for global average
      rect(pos0r[ig], 0, pos1r[ig], 1, lwd = 1.5)

      # panel: add among-GCM span as error bars
      ks <- which(res[, "Exp"] %in% all_CSs[ie] &
          res[, "GCM"] %in% all_GCMs &
          res[, "Variable"] %in% vars[iv])
      dat_rarea <- resX[ks, "RegionArea_km2"]


      # calculate error bars for regions
      temp <- apply(resX[ks, 1 + ibp_labs] / dat_rarea, 1, cumsum)
      temp <- aggregate(t(temp), by = list(res[ks, "Region"]), range,
        drop = FALSE)
      dat_ms <- temp[match(name_regions, temp[, "Group.1"]), ]

      # add error bars of "more severe" areas for regions
      dat_range_more <- dat_ms[, icol_name_plus]
      id <- 2
      rect(xleft = pos0r + dposr * (id - 0.75), ybottom = dat_range_more[, 1],
        xright = pos0r + dposr * (id - 0.25), ytop = dat_range_more[, 2],
        col = mod_cols_dark1["more"], border = NA)

      # add error bars of "less severe" areas for regions
      dat_range_less <- dat_ms[, icol_name_neg]
      id <- 3
      rect(xleft = pos0r + dposr * (id - 0.75), ybottom = dat_range_less[, 1],
        xright = pos0r + dposr * (id - 0.25), ytop = dat_range_less[, 2],
        col = mod_cols_dark1["less"], border = NA)

      # panel: add target GCM
      igcm <- id_GCM_target
      if (id_GCM_target > 0 && all_CSs[ie] %in% reqCSs_per_M[[all_GCMs[igcm]]]) {

        ks <- which(res[, "Exp"] %in% all_CSs[ie] &
            res[, "GCM"] %in% all_GCMs[igcm] &
            res[, "Variable"] %in% vars[iv])
        stopifnot(length(ks) == length(name_regions))
        id_regions <- match(name_regions, res[ks, "Region"])
        ks2 <- ks[id_regions]

        dat_rarea <- resX[ks2, "RegionArea_km2"]
        temp <- resX[ks2, 1 + ibp_labs] / dat_rarea
        dat_gcm <- apply(temp, 1, cumsum)

        # add points for regions
        id <- c(2, 3) #c(1, 4)
        for (k in ibp_labs0) {
          points(pos0r + dposr * (id[k] - 0.5), dat_gcm[k, ],
            pch = 8, font = 1, cex = 0.75, xpd = NA,
            col = mod_cols_dark2[ibp_noc][k])
        }
      }

      # annotate and add axes
      axis(side = 2, labels = ie == 1)

      if (ie == 1) {
        if (iv == ceiling(n_vars / 2)) {
          temp <- if (method == "dir") {
              paste0("SSAI modifies ", ref, "-climate by ", tol_lab, " in")
            } else {
              paste("SSAI-climate diverges from", ref, "by", tol_lab, "in")
            }
          mtext(side = 2, line = 2, xpd = NA,
            text = paste("Geographic extent (%): ", temp))
        }

        mtext(side = 2, line = 1, xpd = NA, text = var_labs[iv])
      }

      if (iv == n_vars) {
        text(x = btemp, y = 0, labels = c(name_regions[ig], name_regions[-ig]),
          srt = 45, adj = c(1, 1), xpd = NA)
      }

      # Add panel identifier
      mtext(side = 3, adj = 0.05, text = letters[i], font = 2, cex = fexp)
      i <- i + 1

    }
  }

  par(par_prev)
  dev.off()
}




plot_dCDFs_global_v2 <- function(data, dataEns, vars, var_labels,
  cell_within_study_TF, reqCSs, deltas_byCS, reqMs,
  col_byCSs, horiz = TRUE, path, ftag) {

  get_scen_colors <- function(k, col, reqMs. = reqMs) {
    adjustcolor(col, offset = c(rep((k / length(reqMs.) - 0.5) / 3, 3), 0))
  }

  n_panels <- if (horiz) c(1, length(vars)) else c(length(vars), 1)
  yt <- seq(0, 1, length.out = 500L)

  fname_fig <- file.path(path,
    paste0("Fig_dCDFs_Global_GCMs_by_Scenarios_", ftag, ".pdf"))

  h.panel <- 50 / 25.4; h.edgeL <- 0.20; h.edgeU <- 0.15 # Heigth of panel and lower/upper edge
  w.panel <- 60 / 25.4; w.edgeL <- 0.25; w.edgeR <- 0.05 # Width of panel and left/right edge

  lmat <- matrix(0, nrow = 2L + n_panels[1L], ncol = 2L + n_panels[2L],
    byrow = FALSE)

  stemp1 <- seq_len(n_panels[1L])
  for (k in seq_len(n_panels[2L])) { # byrow = FALSE
    lmat[-c(1L, 2L + n_panels[1L]), k + 1L] <- (k - 1L) * n_panels[1L] + stemp1
  }

  pdf(
    height = 1.25 * (h.edgeL + h.panel * n_panels[1L] + h.edgeU),
    width = 1.25 * (w.edgeL + w.panel * n_panels[2L] + w.edgeR),
    file = fname_fig)

  layout(lmat,
    heights = c(h.edgeU, rep(h.panel, n_panels[1L]), h.edgeL),
    widths = c(w.edgeL, rep(w.panel, n_panels[2L]), w.edgeR))

  par_prev <- par(mar = c(1, 0.5, 0.1, 0.1), mgp = c(1, 0, 0), tcl = 0.3,
    cex = 1)

  ltys_by_reqMs <- seq_along(reqMs)
  k0 <- 1

  for (iv in seq_along(vars)) {
    vlab <- paste("delta", var_labels[iv])

    xlim_used <- range(sapply(seq_along(reqCSs), function(sc) {
      quantile(data[[sc]][, vars[iv], ],
        probs = c(0.025, 0.975), na.rm = TRUE)
    }), na.rm = TRUE)

    for (sc in seq_along(reqCSs)) {
      sc_name <- reqCSs[sc]

      used_scens <- dimnames(data[[sc_name]])[[3]]
      ids_used_scens <- match(find_reqMs(used_scens), reqMs)
      nscen <- length(ids_used_scens)

      used_ltys <- ltys_by_reqMs[ids_used_scens]
      used_cols <- sapply(ids_used_scens, get_scen_colors, col = col_byCSs[sc])

      if (sc == 1) {
        plot(NA, type = "n", xlim = xlim_used, ylim = c(0, 1),
          ann = FALSE, axes = FALSE, frame.plot = TRUE, xpd = NA,
          col = col_byCSs[sc])
        axis(side = 1)
        do_axis2 <- !horiz || (horiz && iv == 1)
        axis(side = 2, labels = do_axis2)
        if (do_axis2) title(ylab = "Fn(x)", xpd = NA)

        temp <- gregexpr("\n", vlab)[[1L]]
        nlines <- 1L + sum(temp > 0)
        mtext(side = 1, line = nlines, text = vlab, cex = par("cex"), xpd = NA)

        mtext(side = 3, adj = 0.025, text = letters[k0], font = 2, xpd = NA)
        k0 <- k0 + 1
      }

      temp <- dataEns[[sc_name]][cell_within_study_TF, vars[iv]]
      if (!all(is.na(temp))) {
        dem <- ecdf(temp)
        lines(dem, do.points = FALSE, col = col_byCSs[sc], lty = 1, lwd = 2)
      }

      for (isc in seq_len(nscen)) {
        temp <- data[[sc_name]][cell_within_study_TF, vars[iv], isc]
        if (!all(is.na(temp))) {
          dem <- ecdf(temp)
          if (FALSE) {
            # it seems that the lty argument is broken for the
            # 'stepfun' methods of plot and lines
            lines(dem, do.points = FALSE,
              col = adjustcolor(used_cols[isc], alpha.f = 0.7),
              lty = used_ltys[isc], lwd = 1)
          } else {
            lines(x = quantile(dem, probs = yt), y = yt,
              col = adjustcolor(used_cols[isc], alpha.f = 0.7),
              lty = used_ltys[isc], lwd = 1)
          }
        }
      }

      abline(v = 0, lty = 2)
      abline(h = 0.5, lty = 2)
      abline(h = c(0, 1), lty = 2, col = "gray")

      if (sc == length(reqCSs) && iv == 2) {
        legend("bottomright",
          bty = "o", bg = "white", box.col = "white", inset = 0.025, cex = 0.55,
          legend = c(
            paste0("delta(", reqCSs, ", ", cur_to_hist(deltas_byCS), ")"),
            c(reqMs, paste("Ensemble"))),
          fill = temp <- c(col_byCSs, rep(NA, length(reqMs) + 1)),
          border = temp,
          col = c(rep(NA, length(reqCSs)),
            sapply(seq_along(reqMs), get_scen_colors, col = "gray25"), "black"),
          lty = c(rep(NA, length(reqCSs)), ltys_by_reqMs, 1),
          lwd = c(rep(NA, length(reqCSs)), rep(1, length(reqMs)), 2),
          merge = TRUE)
      }
    }
  }

  par(par_prev)
  dev.off()
}



#' Cummulative area by directional change and region plots
#' @param cat_dirs Directional change categories
plot_dCDF2s_regional <- function(data_dvals, data_dir, vars, var_labs,
  SAR, sc_ref, sc_base, cat_dirs = NULL,
  cell_within_study_TF, cell_areas, region_data, name_regions = NULL,
  panel_indep = FALSE, path, ftag) {

  xscale <- 1e6
  xdigits <- 2

  temp1 <- dimnames(data_dvals)
  temp2 <- dimnames(data_dir)
  vars_avail <- intersect(temp1[[2]], temp2[[2]])
  scens_avail <- intersect(temp1[[3]], temp2[[3]])

  stopifnot(vars %in% vars_avail)
  if (is.null(name_regions)) {
    name_regions <- c("Global", sort(as.character(unique(region_data))))
  }

  if (length(name_regions) > 1 && length(vars) > 1) {
    stop("Either one region (or global) and multiple variables OR ",
      "several regions and one variable")
  }

  if (is.null(cat_dirs)) {
    cat_dirs <- define_direction_2deltas()
  }

  col_cats <- c("lightblue", "mediumseagreen", "darkblue",
    "darkred", "purple", "orange")

  # Identify GCMs participating in scenario experiment
  scens <- grep(SAR, scens_avail, value = TRUE)
  gcms_under_scen <- find_reqMs(scens)

  # Determine used cells
  ids <- if (any(name_regions == "Global")) {
      cell_within_study_TF
    } else {
      cell_within_study_TF & match(region_data, name_regions, nomatch = 0) > 0
    }

  dr <- data_dir[, vars, scens, drop = FALSE]
  stopifnot(is.na(dr[ids, , ]) | dr[ids, , ] %in% cat_dirs[, "name"])

  #--- Prepare plot
  temp <- max(length(name_regions), length(vars))
  temp2 <- min(temp, 2)
  n_panels <- c(ceiling(temp / temp2), temp2)

  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  fname <- paste0("Fig_dCDFs_DirectionalChange_", ftag, ".pdf")
  fname <- file.path(path, fname)

  # Figure size
  w.panel <- 3
  w.edgeL <- 0.1; w.edgeI <- 0; w.edgeR <- 0.05 # Width of panel and left/interior/right edge
  h.panel <- 2
  h.edgeL <- 0.1; h.edgeI <- if (panel_indep) 0.1 else 0; h.edgeU <- 0.1 # Heigth of panel and lower/interior/upper edge

  # Figure layout
  lmat <- matrix(0L,
    nrow = 2 + 2 * n_panels[1] - 1L, ncol = 2 + 2 * n_panels[2] - 1,
    byrow = TRUE)

  stemp2 <- seq_len(n_panels[2L])
  for (k in seq_len(n_panels[1L])) { # byrow = TRUE
    #lmat[k + 1L, -c(1L, 2L + n_panels[2L])] <- (k - 1L) * n_panels[2L] + stemp2
    lmat[(k - 1) * 2 + 2, seq_len(ncol(lmat)) %% 2 == 0] <- (k - 1L) * n_panels[2L] + stemp2
  }

  temp <- rep(c(h.panel, h.edgeI), n_panels[1L])
  layout_heights <- c(h.edgeU, temp[-length(temp)], h.edgeL)
  temp <- rep(c(w.panel, w.edgeI), n_panels[2L])
  layout_widths <- c(w.edgeL, temp[-length(temp)], w.edgeR)


  # Figure device
  fexp <- 1.5; cexp <- 1
  pdf(
    height = fexp * (h.edgeL + h.panel * n_panels[1] + h.edgeI * (n_panels[1] - 1) + h.edgeU),
    width = fexp * (w.edgeL + w.panel * n_panels[2] + w.edgeI * (n_panels[2] - 1) + w.edgeR),
    file = fname)

  layout(lmat, heights = layout_heights, widths = layout_widths)

  mar <- if (panel_indep) c(1.5, 1.75, 0.75, 0.5) else rep(0.5, 4)
  par_prev <- par(mar = mar, mgp = c(1, 0, 0), tcl = 0.3,
    cex = cexp)

  #--- Plot data on panels
  k0 <- 1
  for (iv in seq_along(vars)) {
    dx <- data_dvals[, vars[iv], scens, drop = FALSE]
    data <- data.frame(area_km2 = cell_areas[, "km2"],
      dr = dr[, vars[iv], ], dx = dx[, 1, ], stringsAsFactors = FALSE)
    ylim <- range(dx[ids, , ], na.rm = TRUE)

    for (ir in seq_along(name_regions)) {
      id_region <- if (name_regions[ir] == "Global") {
          cell_within_study_TF
        } else {
          cell_within_study_TF & region_data %in% name_regions[ir]
        }

      area_region_km2 <- sum(data[id_region, "area_km2"])

      res <- expand.grid(GCM = gcms_under_scen, cat = cat_dirs[, "name"],
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      Ns <- nrow(res)
      vals <- matrix(vector(mode = "list", length = Ns), ncol = 1)

      for (k1 in seq_len(Ns)) {
        # Subset data to region, directional change category, and GCM
        id_cols <- grep(make.names(res[k1, "GCM"]), colnames(data))
        temp <- data[id_region, c(1, id_cols)]
        id_cat <- temp[, 2] %in% res[k1, "cat"]
        temp <- temp[id_cat, ]

        # Sort data: increasing/decreasing values of for pos/neg delta x
        is_neg <- all(temp[, 3] < 0, na.rm = TRUE)
        ords <- sort.list(temp[, 3], decreasing = is_neg)
        temp2 <- temp[ords, ]

        # Calculate cummulative area
        vals[[k1]] <- list(
          cum_area_km2 = cumsum(temp2[, "area_km2"]),
          x = temp2[, 3])
      }

      # Plot
      xlim <- c(1e3, max(sapply(vals, function(x) max(x[["cum_area_km2"]])),
        na.rm = TRUE))
      if (any(is.infinite(xlim))) {
        plot.new()
      } else {
        plot(NA, xlim = xlim / xscale, ylim = ylim, type = "n", log = "",
          xlab = paste0("Cumulative area (", name_regions[ir], " = ",
            round(area_region_km2 / xscale, xdigits), "; 10^6 km2)"), xpd = NA,
          ylab = paste0("d(", SAR, ", ", sc_ref, "): ", var_labs[iv]))
        abline(h = 0, col = "darkgray")

        for (k1 in seq_len(Ns)) {
          lines(vals[[k1]][["cum_area_km2"]] / xscale, vals[[k1]][["x"]],
            col = col_cats[which(cat_dirs[, "name"] %in% res[k1, "cat"])],
            lty = 1 + which(gcms_under_scen %in% res[k1, "GCM"]))
        }

        if (k0 == 1) {
          legend("topleft", inset = 0.025,
            title = paste0("d(", SAR, ", ", substr(sc_base, 1, 3), ") vs d(",
              sc_ref, ", ", substr(sc_base, 1, 3), ")"),
            legend = cat_dirs[, "name"], cex = 0.65 * cexp,
            col = col_cats, lwd = 2)
        }
      }

      mtext(side = 3, adj = -0.05, text = paste0("(", letters[k0], ")"),
        font = 2, cex = cexp)
      k0 <- k0 + 1
    }
  }

  par(par_prev)
  dev.off()
}




add_panel_vp <- function(x, xlab = "Variance (%)", plot = TRUE) {
  temp <- rev(x)
  vp <- t(cbind(temp, 1 - temp))

  mp <- barplot(vp,
    xlab = xlab, names.arg = rep("", length(temp)), axes = FALSE,
    plot = plot, beside = FALSE, xlim = c(0, 1), horiz = TRUE, las = 1)

  if (plot) {
    axis(1, labels = nchar(xlab) > 0, xpd = NA)
    title(xlab = xlab, xpd = NA)
  }

  space <- 0.2
  ylim <- c(space, ncol(vp) * (1 + space))

  invisible(list(mp = rev(mp), ylim = ylim))
}

# Plot variance partitioning
plot_vp_predictor_contributions <- function(data, response, resp_labs,
  GCMs, SARs, path, ftag) {

  # Figure size
  n_panels <- c(length(GCMs), length(SARs))
  w.panel <- 3
  w.edgeL <- 2.25; w.edgeI <- 0.1; w.edgeR <- 0.05 # Width of panel and left/interior/right edge
  h.panel <- 0.6 * length(response)
  h.edgeL <- 0.4; h.edgeI <- 0; h.edgeU <- 0.1 # Heigth of panel and lower/interior/upper edge

  # Figure layout
  lmat <- matrix(0L,
    nrow = 2 + 2 * n_panels[1] - 1L, ncol = 2 + 2 * n_panels[2] - 1,
    byrow = TRUE)

  stemp2 <- seq_len(n_panels[2L])
  for (k in seq_len(n_panels[1L])) { # byrow = TRUE
    lmat[(k - 1) * 2 + 2, seq_len(ncol(lmat)) %% 2 == 0] <- (k - 1L) * n_panels[2L] + stemp2
  }

  temp <- rep(c(h.panel, h.edgeI), n_panels[1L])
  layout_heights <- c(h.edgeU, temp[-length(temp)], h.edgeL)
  temp <- rep(c(w.panel, w.edgeI), n_panels[2L])
  layout_widths <- c(w.edgeL, temp[-length(temp)], w.edgeR)

  # Figure device
  fexp <- 1
  png(units = "in", res = 150,
    height = fexp * (h.edgeL + h.panel * n_panels[1] + h.edgeI * (n_panels[1] - 1) + h.edgeU),
    width = fexp * (w.edgeL + w.panel * n_panels[2] + w.edgeI * (n_panels[2] - 1) + w.edgeR),
    file = file.path(path, paste0("Fig_VariancePartition_", ftag, ".png")))

  layout(lmat, heights = layout_heights, widths = layout_widths)

  par_prev <- par(xaxs = "i", mar = rep(0.5, 4), mgp = c(1, 0, 0), tcl = 0.2,
    cex = fexp)

  i <- 1
  bp <- add_panel_vp(rep(0.5, length(response)), plot = FALSE)

  do_ann_GCMSAR <- all(n_panels > 1)

  for (sc in seq_along(GCMs)) {
    for (sar in SARs) {
      x <- data[sar, GCMs[sc], "R2", ]

      if (anyNA(x)) {
        plot(bp[["mp"]], xlim = c(0, 1), ylim = bp[["ylim"]],
          type = "n", axes = FALSE, ann = FALSE)

      } else {
        needs_xlab <- GCMs[sc] == GCMs[length(GCMs)]
        add_panel_vp(x, xlab = if (needs_xlab) "Variance (%)" else "")

        mtext(side = 3, text = letters[i], font = 2, adj = 0.025,
          cex = par("cex"))
        if (do_ann_GCMSAR) {
          mtext(side = 3, text = paste(sar, "|", GCMs[sc]), adj = 0.9)
        }
      }

      if (sar == SARs[[1]]) {
        text(x = -0.05, y = bp[["mp"]], labels = resp_labs, adj = 1, xpd = NA)
      }

      i <- i + 1
    }
  }

  par(par_prev)
  dev.off()
}



regress_predictor_contributions <- function(data, response, resp_labs,
  predictors, ref_condition, cc_condition, SAR_conditions, cell_within_study_TF,
  cell_areas, path, ftag) {

  temp <- dimnames(data[["vals"]])
  vars_avail <- temp[[2]]
  scens_avail <- temp[[3]]

  stopifnot(c(response, predictors) %in% vars_avail,
    any(grepl(ref_condition, scens_avail)),
    any(grepl(cc_condition, scens_avail)),
    sapply(SAR_conditions, function(gi) any(grepl(gi, scens_avail))))

  # Identify GCMs participating in both conditions
  sc_ref <- grep(ref_condition, scens_avail, value = TRUE)
  gcms_under_refs <- find_reqMs(sc_ref)
  sc_cc <- grep(cc_condition, scens_avail, value = TRUE)
  gcms_under_ccs <- find_reqMs(sc_cc)
  sc_sars <- sapply(SAR_conditions, function(gi)
    grep(gi, scens_avail, value = TRUE))
  gcms_under_sars <- lapply(sc_sars, find_reqMs)

  gcms_under_both <- lapply(gcms_under_sars, intersect, y = gcms_under_ccs)
  gcms_under_any <- names(sort(sapply(unique(unlist(gcms_under_both)),
    function(x) grep(x, scens_avail)[1])))

  scens_cc_used <- scens_sar_used <- list()
  for (sar in SAR_conditions) {
    temp <- match(gcms_under_both[[sar]], gcms_under_ccs, nomatch = 0)
    scens_cc_used[[sar]] <- sc_cc[temp]
    temp <- match(gcms_under_both[[sar]], gcms_under_sars[[sar]], nomatch = 0)
    scens_sar_used[[sar]] <- sc_sars[[sar]][temp]
  }

  # Formula
  temp <- paste(sapply(c("Ref", "ClimC", "GeoEng"), function(x)
    paste(paste(x, predictors, sep = "."), collapse = " * ")),
    collapse = " + ")
  vp_formula <- paste0('y ~ Region * (', temp, ')')

  # Output container
  m_term_labels <- labels(terms(as.formula(vp_formula)))
  np <- 2 + length(m_term_labels)
  res <- array(NA,
    dim = c(length(SAR_conditions), length(gcms_under_any), np, length(response)),
    dimnames = list(SAR_conditions, gcms_under_any,
      c("R2", paste0("SS_", c(m_term_labels, "Residuals"))), response))

  # Prepare predictor data
  x_ref <- reshape2::melt(data[["vals"]][, predictors, sc_ref, drop = FALSE])
  x_ref <- reshape2::dcast(x_ref, Var3 + Var1 ~ Var2, )
  colnames(x_ref) <- c("fGCM", "ID", predictors)

  dat_preds_template <- data.frame(
    Region = region_data,
    Ref = x_ref[, -c(1:2)])

  idata_ref <- paste0("delta_",
    strsplit(ref_condition, split = ".", fixed = TRUE)[[1]][1], "_abs")
  idata_cc <- paste0("delta_",
    strsplit(cc_condition, split = ".", fixed = TRUE)[[1]][1], "_abs")

  #--- Calculate
  for (sar in SAR_conditions) {
    for (sc in seq_along(gcms_under_both[[sar]])) {
      gcm <- gcms_under_both[[sar]][sc]

      x_cc <- reshape2::melt(
        data[[idata_ref]][, predictors, scens_cc_used[[sar]][sc],
        drop = FALSE])
      x_cc <- reshape2::dcast(x_cc, Var3 + Var1 ~ Var2, )
      colnames(x_cc) <- c("fGCM", "ID", predictors)

      x_gi <- reshape2::melt(
        data[[idata_ref]][, predictors, scens_sar_used[[sar]][sc],
        drop = FALSE])
      x_gi <- reshape2::dcast(x_gi, Var3 + Var1 ~ Var2, )
      colnames(x_gi) <- c("fGCM", "ID", predictors)

      dx <- cbind(dat_preds_template,
        ClimC = x_cc[, -c(1:2)],
        GeoEng = x_gi[, -c(1:2)])

      for (iv in seq_along(response)) {
        # Prepare response data
        y <- reshape2::melt(
          data[[idata_cc]][, response[iv], scens_sar_used[[sar]][sc],
          drop = FALSE])
        y <- reshape2::dcast(y, Var3 + Var1 ~ Var2)[, 3]

        # Variance partitioning / R2
        m <- lm(as.formula(vp_formula), weights = cell_areas[, "rel"],
          data = dx, subset = cell_within_study_TF)

        res[sar, gcm, "R2", iv] <- summary(m)[["r.squared"]]

        sm <- anova(m)
        res[sar, gcm, -1, iv] <- sm[, "Sum Sq"]

        # Check that R2 == sum(SS of residuals) / sum(SS total)
        stopifnot(all.equal(res[sar, gcm, 1, iv],
          1 - sum(res[sar, gcm, np, iv]) / sum(res[sar, gcm, -1, iv])))
      }
    }
  }


  #---Write vp to spreadsheet file
  write.csv(reshape2::melt(res), file = file.path(path,
    paste0("Table_VariancePartition_", ftag, "_All.csv")))


  #--- Plot variance partitioning for individual GCMs
  for (gcm in gcms_under_any) {
    plot_vp_predictor_contributions(data = res,
      response = response,resp_labs = resp_labs,
      GCMs = gcm,
      SARs = SAR_conditions[sapply(gcms_under_both, function(x) gcm %in% x)],
      path = path, ftag = paste0(ftag, "_", gcm))
  }

  #--- Plot variance partitioning for all G-GCM combinations
  plot_vp_predictor_contributions(data = res,
    response = response, resp_labs = resp_labs,
    GCMs = gcms_under_any, SARs = SAR_conditions,
    path = path, ftag = paste0(ftag, "_All"))


  invisible(res)
}





add_panel_bp <- function(x_rel, x_rel_mean, x_rel_median, xlim = NULL,
  xlab = "Change (%)", names, axes = TRUE, xann = TRUE, yann = TRUE,
  cols = NULL) {

  pt_cex <- 0.75 * par("cex")

  # Convert to percentages
  x_rel <- 100 * x_rel
  xlim <- if (!is.null(xlim)) 100 * xlim else range(x_rel)
  x_rel_mean <- 100 * x_rel_mean
  x_rel_median <- 100 * x_rel_median

  # Plot
  ats <- seq_len(ncol(x_rel))
  ys <- array(rep(ats, each = nrow(x_rel)), dim = dim(x_rel))

  plot(as.vector(unlist(x_rel)), y = as.vector(ys),
    xlim = xlim, ylim = range(extend_range(ats, 1.2)),
    log = "x", cex = pt_cex, axes = FALSE, ann = FALSE,
    frame.plot = TRUE)

  if (axes) {
    axis(side = 1, labels = xann)
  }

  if (xann) {
    title(xlab = xlab, xpd = NA)
  }

  if (yann) {
    opar <- par(mgp = c(1, 0.5, 0))
    on.exit(par(opar))
    axis(side = 2, at = ats, labels = names, las = 1, xpd = NA,
      cex = pt_cex)
  }

  # add variables
  if (is.null(cols)) {
    pch <- 1
    bg <- par("bg")
  } else {
    pch <- 21
    bg <- cols
  }
  cols <- "black"

  points(x = unlist(x_rel), y = rep(ats, each = nrow(x_rel)), pch = pch,
    col = cols, bg = bg, cex = pt_cex)

  # add means and medians
  points(x = x_rel_mean, y = ats, pch = 4, col = "red", lwd = 2,
    cex = pt_cex)
  points(x = x_rel_median, y = ats, pch = 5, col = "orange", lwd = 2,
    cex = pt_cex)
}


plot_GlobalEnsembleValues <- function(data, data2, ens_method = "relative",
  ens_wise = "EnsCW", ens_fun = "mean", ens_funCT, ens_outs, ens_types,
  labs_ens_types, scen_SARs, reqMs, path, ftag) {

  # Prepare data
  x1 <- data[ens_method, ens_wise][[1]]
  x2 <- data2[paste0(ens_method, "_", ens_funCT), ens_wise]
  var_labs <- unique(x1[, "Var1"])

  # Find reqMs (GCMs)
  ims1 <- match(reqMs, gsub(".", "-", colnames(x1), fixed = TRUE))
  ims2 <- lapply(x2, function(x)
    match(reqMs, gsub(".", "-", colnames(x), fixed = TRUE)))

  # Find ensembles based on `ens_fun`
  ies1 <- match(paste(scen_SARs, ens_fun, sep = "."), colnames(x1))
  ies2 <- lapply(x2, function(x)
    match(paste(scen_SARs, ens_fun, sep = "."), colnames(x)))

  #icols1 <- rev(seq_len(ncol(x1))[-c(1:3)])
  #icols2 <- lapply(x2, function(x) rev(seq_len(ncol(x))[-c(1:2)]))

  # Figure size
  n_panels <- c(length(scen_SARs), length(ens_types))
  w.panel <- 3
  w.edgeL <- 1.75; w.edgeI <- 0; w.edgeR <- 0.05 # Width of panel and left/interior/right edge
  h.panel <- 0.25 * length(reqMs)
  h.edgeL <- 0.3; h.edgeI <- 0; h.edgeU <- 0.1 # Heigth of panel and lower/interior/upper edge

  # Point/variable colors
  #ens_var_id <- factor(ens_var_id)
  #col_vars <- gray.colors(n = nlevels(ens_var_id),
  #  start = 0.8, end = 0.4)[ens_var_id]
  col_vars <- viridisLite::viridis(n = nrow(ens_outs))

  # Figure layout
  lmat <- matrix(0L,
    nrow = 2 + 2 * n_panels[1] - 1, ncol = 2 + 2 * n_panels[2] - 1,
    byrow = FALSE)

  stemp2 <- seq_len(n_panels[2L])
  for (k in seq_len(n_panels[1L])) { # byrow = FALSE
    lmat[(k - 1) * 2 + 2, seq_len(ncol(lmat)) %% 2 == 0] <- (stemp2 - 1L) * n_panels[1L] + k
    #lmat[(k - 1) * 2 + 2, seq_len(ncol(lmat)) %% 2 == 0] <- (k - 1L) * n_panels[2L] + stemp2
  }

  temp <- rep(c(h.panel, h.edgeI), n_panels[1L])
  layout_heights <- c(h.edgeU, temp[-length(temp)], h.edgeL)
  temp <- rep(c(w.panel, w.edgeI), n_panels[2L])
  layout_widths <- c(w.edgeL, temp[-length(temp)], w.edgeR)

  # Figure device
  fexp <- 1.5
  fname <- file.path(path,
    paste0("Fig_GlobalEnsembleValues_All_", ftag, ".png"))

  png(units = "in", res = 150,
    height = fexp * (h.edgeL + h.panel * n_panels[1] + h.edgeI * (n_panels[1] - 1) + h.edgeU),
    width = fexp * (w.edgeL + w.panel * n_panels[2] + w.edgeI * (n_panels[2] - 1) + w.edgeR),
    file = fname)

  layout(lmat, heights = layout_heights, widths = layout_widths)

  par_prev <- par(mar = rep(0.5, 4), mgp = c(1, 0, 0), tcl = 0.2, las = 1,
    cex = fexp)

  i <- 1

  for (k1 in seq_along(ens_types)) {
    ids1_k1 <- x1[, "ens_types"] == ens_types[k1]
    ids2_k1 <- lapply(x2, function(x) x[, "ens_types"] == ens_types[k1])

    xlim <- range(
      x1[ids1_k1, c(ims1, ies1)],
      unlist(lapply(seq_along(x2), function(j)
        x2[[j]][ids2_k1[[j]], c(ims2[[j]], ies2[[j]])])),
      na.rm = TRUE)

    for (k2 in seq_along(scen_SARs)) {
      ids1 <- which(ids1_k1 & x1[, "Var3"] == scen_SARs[k2])
      ids2 <- lapply(seq_along(x2), function(j) which(ids2_k1[[j]] &
        x2[[j]][, "Var3"] == scen_SARs[k2]))

      needs_xann <- k2 == length(scen_SARs)
      needs_yann <- k1 == 1

      add_panel_bp(
        x_rel = x1[ids1, rev(c(ims1, ies1[k2]))],
        x_rel_mean = x2[["relative_mean"]][ids2[[1]], rev(c(ims2[[1]], ies2[[1]][k2]))],
        x_rel_median = x2[["relative_median"]][ids2[[2]], rev(c(ims2[[2]], ies2[[2]][k2]))],
        xlim = xlim, xlab = "Remaining change (%)",
        names = rev(c(reqMs, paste("Ensemble", ens_fun))),
        axes = TRUE, xann = needs_xann, yann = needs_yann, cols = col_vars)

      if (k1 == 1 && k2 == 1) {
        legend("bottomright", inset = 0.05,
          legend = c(var_labs, ens_funCT), cex = 0.45,
          col = c(rep("black", length(var_labs)), "red", "orange"),
          pch = c(rep(21, length(var_labs)), 4, 5),
          pt.bg = c(col_vars, NA, NA))
      }

      mtext(side = 3, text = letters[i], font = 2, adj = 0.025,
        cex = par("cex"))
      mtext(side = 3, text = paste(scen_SARs[k2], "|", labs_ens_types[k1]),
        adj = 0.9)
      i <- i + 1
    }
  }

  par(par_prev)
  dev.off()
}


cur_to_hist <- function(x) {
  x <- gsub("current", "historical", x, ignore.case = TRUE)
  gsub("cur", "hist", x, ignore.case = TRUE)
}



#' Calculate range of variogram in kilometers
#'
#' @param x A RasterLayer in un-projected coordinates (i.e.,
#'   longitude/latitude)
variogram_range <- function(x, project_to_utm = TRUE,
  sub_samplepoints_N = NULL) {

  stopifnot(requireNamespace("sp"), requireNamespace("raster"),
    requireNamespace("automap"), requireNamespace("gstat"))

  # Create spatial points from raster
  points <- raster::rasterToPoints(x, spatial = TRUE)

  if (!is.null(sub_samplepoints_N)) {
    set.seed(2017)
    temp <- sample(seq_len(nrow(points)), sub_samplepoints_N,
      replace = FALSE)
    points <- points[temp, ]
  }

  if (project_to_utm) {
    # if region, then project to UTM coordinates with distance units of km
    # --> variogram will calculate Euclidean distances in map units
    xytemp <- sp::coordinates(points)
    utm_zone <- get_UTM_Zone(xytemp[, 1], xytemp[, 2])
    utm_NS <- if (mean(xytemp[, 2]) > 0) " +north" else " +south"
    points_prj <- sp::spTransform(points,
      CRSobj = sp::CRS(paste0("+proj=utm +zone=", utm_zone, utm_NS,
        " +datum=WGS84 +units=km +no_defs")))

  } else {
    # if global, then don't project, i.e. use long/lat with WGS84 datum
    # and distance units of meters
    # --> variogram will calculate great circle distances in kilometers(!)
    points_prj <- sp::spTransform(points,
      CRSobj = sp::CRS("+proj=longlat +datum=WGS84 +units=m +no_defs"))
  }

  names(points_prj) <- "target"

  # determine variogram; see gstat::variogram
  fittedVar <- automap::autofitVariogram(target ~ 1,
    input_data = points_prj,
    miscFitOptions = list(merge.small.bins = TRUE))

  # variogram range in km:
  # estimate is only valid if some variation in data (i.e., sill > 0)
  psill_sum <- sum(fittedVar[["var_model"]][, "psill"])

  if (abs(psill_sum) > sqrt(.Machine$double.eps)) {
    fittedVar[["var_model"]][2, "range"]
  } else NaN
}


# Estimates spatial patterns
range_spatialpattern <- function(x, meta, cell_within_study_TF, region_data,
  name_regions = NULL, mod_with_all = TRUE, expected_dy_lt_dx = TRUE,
  responses = NULL, path_temp = getwd(), ftag = "", verbose = FALSE) {

  dir.create(path_temp, recursive = TRUE, showWarnings = FALSE)

  stopifnot(requireNamespace("landscapemetrics"),
    packageVersion("landscapemetrics") >= "1.0.1")

  cl <- define_direction_2deltas(TRUE)

  if (is.null(name_regions)) {
    name_regions <- c("Global", sort(as.character(unique(region_data))))
  }

  temp <- dimnames(x)
  vars <- temp[[2]]
  n_vars <- length(vars)
  GCMs <- temp[[3]]
  n_GCMs <- length(GCMs)

  expected_dy_lt_dx <- if (is.logical(expected_dy_lt_dx)) {
      rep_len(expected_dy_lt_dx, n_vars)
    } else {
      expected_dy_lt_dx > 0
    }
  expected_dy_lt_dx <- as.vector(expected_dy_lt_dx)
  stopifnot(length(expected_dy_lt_dx) == n_vars)

  # Prepare data (part 1)
  x2 <- array(NA_integer_, dim = dim(x), dimnames = dimnames(x))
  for (iv in seq_along(vars)) {
    if (mod_with_all) {
      x2[, iv, ] <- as.integer(factor_6directions(x[, iv, ], cl = cl))
    } else {
      # Remap
      x2[, iv, ] <- as.integer(remap_directions_6to2(as.vector(x[, iv, ]),
        expected_dy_lt_dx[iv]))
    }
  }
  x <- x2
  rm(x2)

  # Possible response variables
  resp <- c("MoranI", "GearyC", "VariogramRange_km", "Aggregation",
    "Contagion", "Division", "LikeAdjacency")

  if (is.null(responses)) {
    responses <- resp
  }

  # Output container
  res <- expand.grid(Region = name_regions, Variable = vars, Scenario = GCMs,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  res <- cbind(res,
    data.frame(matrix(NA, nrow = nrow(res), ncol = length(resp),
      dimnames = list(NULL, resp)), stringsAsFactors = FALSE))

  # Calculate
  for (k in seq_len(nrow(res))) {
    # Skip if no expectation for variable
    iv <- which(res[k, "Variable"] == vars)
    if (is.na(expected_dy_lt_dx[iv])) next

    # Prepare data (part 2)
    data <- x[, res[k, "Variable"], res[k, "Scenario"]]

    id_region <- if (res[k, "Region"] == "Global") {
        cell_within_study_TF
      } else {
        cell_within_study_TF & region_data %in% res[k, "Region"]
      }
    data[!id_region] <- NA

    # Skip if not enough data
    n_values <- sum(!is.na(data))
    if (n_values < 1) next

    # Check if already done and stored on disk
    temp <- list(data = data, id_region = id_region,
      Variable = res[k, "Variable"], Scenario = res[k, "Scenario"],
      Region = res[k, "Region"], cl = cl, mod_with_all = mod_with_all,
      expected_dy_lt_dx = expected_dy_lt_dx[iv])
    run_id <- digest::digest(temp, algo = "sha1")

    get_fname <- function(tag) {
      file.path(path_temp, paste0(ftag, "_SpatialPatterns_Temp-", run_id,
        "_", tag, ".rds"))
    }

    raster::removeTmpFiles(h = 0.5)

    if (verbose) {
      print(paste(Sys.time(), "calculate spatial patterns for",
        shQuote(res[k, "Region"]), shQuote(res[k, "Variable"]),
        shQuote(res[k, "Scenario"])))
    }

    # Create raster from data
    rasterLayer <- create_raster_from_variables(meta, data)
    rasterLayer <- raster::trim(rasterLayer)

    # Moran's I and Geary's C (Queen's case with row-standardization)
    for (id in responses) {
      fname <- get_fname(id)

      if (file.exists(fname)) {
        temp <- readRDS(fname)
      } else {
        n_crit <- switch(id,
            MoranI = 0,
            GearyC = 0,
            VariogramRange_km = 3,
            Aggregation = 1,
            Contagion = 1,
            Division = 1,
            LikeAdjacency = 1
          )

        id_fun <- switch(id,
            MoranI = raster::Moran,
            GearyC = raster::Geary,
            VariogramRange_km = variogram_range,
            Aggregation = landscapemetrics::lsm_l_ai,
            Contagion = landscapemetrics::lsm_l_contag,
            Division = landscapemetrics::lsm_l_division,
            LikeAdjacency = landscapemetrics::lsm_l_pladj
          )

        ns <- try(getNamespaceName(environment(id_fun)), silent = TRUE)
        id_args <- if (isTRUE(ns == "landscapemetrics")) {
            list(landscape = rasterLayer)
          } else {
            list(x = rasterLayer)
          }

        if (id == "VariogramRange_km") {
          id_args <- c(id_args, project_to_utm = res[k, "Region"] != "Global")
        } else if (id == "Division") {
          id_args <- c(id_args, directions = 8)
        }

        temp <- if (n_values > n_crit) {
            xtemp <- try(do.call(what = id_fun, args = id_args),
              silent = !verbose)
            if (inherits(temp, "try-error")) {
              NaN
            } else if (isTRUE(ns == "landscapemetrics") && inherits(xtemp, "tbl")) {
              as.numeric(xtemp[1, "value"])
            } else {
              xtemp
            }
          } else {
            NaN
          }

        saveRDS(temp, file = get_fname(id))
      }

      res[k, id] <- temp
    }
  }

  res
}



} # endif not exists "SGEcologicalDrought_0_Functions"

SGEcologicalDrought_0_Functions <- TRUE
