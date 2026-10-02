#---------------------------------------------------------------
# Name:        nbmat (swash package)
# Purpose:     Neighborhood matrix and spatial statistics
# Author:      Thomas Wieland 
#              ORCID: 0000-0001-5168-9846
#              mail: geowieland@googlemail.com
# Version:     2.0.0
# Last update: 2026-09-26 14:02
# Copyright (c) 2022-2026 Thomas Wieland
#---------------------------------------------------------------


# Neighborhood matrix and spatial statistics functions:

# Class nbmatrix:
setClass(
  "nbmatrix",
  slots = list(
    nb = "ANY",
    nbmat = "data.frame",
    weights = "listw",
    nbcount = "data.frame",
    regions_col_name = "character",
    geodata = "sf",
    N = "numeric",
    N_neighbors = "numeric",
    polygon_data = "sf",
    nbmat_data = "data.frame",
    nbmat_data_aggregate = "data.frame",
    func = "character",
    moran = "list",
    getisord = "list",
    gstar = "list",
    data_col = "character",
    config = "list",
    timestamp = "list"
  ),
  prototype = list(
    polygon_data = st_sf(
      unit = character(0),
      geometry = st_sfc()
    ),
    nbmat_data = data.frame(),
    nbmat_data_aggregate = data.frame(),
    func = character(),
    moran = list(),
    getisord = list(),
    gstar = list(),
    data_col = character(),
    config = list(),
    timestamp = list()
  )
)


# Methods for class nbmatrix:

setMethod(
  "summary",
  "nbmatrix",
  function(object) {
    
    data_col <- object@data_col
    N <- object@N
    N_neighbors <- object@N_neighbors
    
    cat(nbmatrix_description, "\n\n")
    
    cat("Input data\n")
    cat(sprintf("  Unique units total     : %s\n", N[2]))
    cat(sprintf("  Indicator              : %s\n", data_col))
    
    cat("Neighborhood statistics\n")
    cat(sprintf("  Units with neighbors   : %s\n", N_neighbors[1]))
    cat(sprintf("  Total no. of neighbors : %s\n", N_neighbors[2]))
    cat(sprintf("  Av. no. of neighbors   : %.3f\n", N_neighbors[3]))
    
    moran <- object@moran
    if(length(moran) > 0) {
      cat("\n")
      if(moran$config$method == "Moran I test under randomisation") {
        cat(paste0(spatial_statistics_descriptions[["moran"]], " under randomization"), "\n")
      } else {
        cat(paste0(spatial_statistics_descriptions[["moran"]], " under normality"), "\n")
      }
      cat(sprintf("  Estimate               : %.3f\n", object@moran$I_estimate))
      cat(sprintf("  Expected value         : %.3f\n", object@moran$I_expected))
      cat(sprintf("  p value                : %.3f\n", object@moran$I_p))
    }
    
    getisord <- object@getisord
    if(length(getisord) > 0) {
      cat("\n")
      cat(paste0(spatial_statistics_descriptions[["getisord"]], " under randomization"), "\n")
      cat(sprintf("  Estimate               : %.3f\n", object@getisord$G_estimate))
      cat(sprintf("  Expected value         : %.3f\n", object@getisord$G_expected))
      cat(sprintf("  p value                : %.3f\n", object@getisord$G_p))
    }
    
    gstar <- object@gstar
    if(length(gstar) > 0) {
      
      results_col_names <- gstar$results_col_names
      Gi_estimates <- gstar$Gi_estimates
      Gi_estimates <- 
        .rename_cols(
          Gi_estimates,
          results_col_names
        )
      Gi_estimates <- Gi_estimates[, 1:(ncol(Gi_estimates) - 2)]
      
      cat("\n")
      cat(paste0(spatial_statistics_descriptions[["gstar"]]), "\n")
      cat(paste0("(Showing first and last 5 of ", nrow(Gi_estimates), " cases)"), "\n")
      
      print(head(Gi_estimates, 5))
      cat("...", "\n")
      print(tail(Gi_estimates, 5))
      
      cat("\n")
      cat("Use your_nbmatriobject@gstar$Gi_estimates to access the full table", "\n")
      
    }
    
    invisible(object)
    
  }
)

setMethod(
  "print", 
  "nbmatrix", 
  function(x) {
    
    cat(paste0(nbmatrix_description, " with ", x@N[2], " valid spatial units"), "\n")
    cat ("Use summary() for details")
    
    invisible(x)
    
  }
)

setMethod(
  "show", 
  "nbmatrix", 
  function(object) {
    
    cat(paste0(nbmatrix_description, " with ", object@N[2], " valid spatial units"), "\n")
    cat ("Use summary() for details")
    
    invisible(object)
    
  }
)

setGeneric(
  "nbstat",
  function(
    object,
    link_data,
    ID_col,
    data_col,
    func = "sum",
    verbose = FALSE
  ) {
    standardGeneric("nbstat")
  }
)

setMethod(
  "nbstat", 
  "nbmatrix", 
  function(
    object,
    link_data, 
    ID_col, 
    data_col, 
    func = "sum",
    verbose = FALSE
  ) {
    
    nbmat <- object@nbmat
    N <- object@N
    
    if(isTRUE(verbose)) {
      cat("Loading and checking link data ... ")
    }
    
    if("geometry" %in% colnames(link_data)) {
      link_data$geometry <- NULL
    }
    
    if(
      !ID_col %in% colnames(link_data)
      | !data_col %in% colnames(link_data)) {
      stop("Specified link data does not include unique ID column and/or data column")
    }
    
    if(!is.numeric(link_data[[data_col]])) {
      stop(paste0("Specified data column '", data_col, "' is not numeric"))
    }
    
    link_data_original_nrow <- nrow(link_data)
    
    link_data <- link_data[(!is.na(link_data[[ID_col]])) & ((!is.na(link_data[[data_col]]))),]
    
    geodata <- object@geodata
    regions_col_name <- object@regions_col_name
    
    polygon_data <-
      merge(
        geodata[c(regions_col_name, "geometry")],
        link_data[c(ID_col, data_col)],
        by.x = regions_col_name,
        by.y = ID_col
      )
    
    if(isTRUE(verbose)) {
      
      cat("OK", "\n")
      
      if(nrow(link_data) < link_data_original_nrow) {
        warning(paste0("From the link data, ", link_data_original_nrow-nrow(link_data), " rows were deleted because of NA values in unique ID or data col"))
      }
      
    }
    
    if(isTRUE(verbose)) {
      cat(paste0("Defined column '", data_col,"' as data column to be analyzed"), "\n")
    }
    
    nbmat_data <- 
      merge (
        nbmat, 
        link_data[c(ID_col, data_col)],
        by.x = "ID2_nb", 
        by.y = ID_col
      )
    
    nbmat_data_aggregate <- data.frame()
    
    if(!is.null(func)) {
      
      if(isTRUE(verbose)) {
        cat(paste0("Calculating neighborhood statistic for data column '", data_col,"' with function '", func, "' ... "))
      }
      
      nbmat_data_aggregate <- 
        aggregate(
          nbmat_data[[data_col]],
          by = list(nbmat_data$ID2),
          FUN = func,
          na.rm = TRUE
        )
      
      colnames(nbmat_data_aggregate) <- c("ID2", paste0(data_col, "_", func))
      
      if (isTRUE(verbose)) {
        cat("OK", "\n")
      }      
      
    } else {
      func <- character()
    }
    
    nbmatrix_object <-
      new(
        "nbmatrix",
        nb = object@nb,
        nbmat = object@nbmat,
        weights = object@weights,
        nbcount = object@nbcount,
        regions_col_name = object@regions_col_name,
        geodata = object@geodata,
        N = object@N,
        N_neighbors = object@N_neighbors,
        polygon_data = polygon_data,
        nbmat_data = nbmat_data,
        nbmat_data_aggregate = nbmat_data_aggregate,
        func = func,
        moran = object@moran,
        getisord = object@getisord,
        gstar = object@gstar,
        data_col = data_col,
        config = object@config,
        timestamp = object@timestamp
      )
    
    nbmatrix_object <- add_timestamp(
      nbmatrix_object,
      function_or_method = "nbstat",
      process = paste0("Defined column '", data_col,"' as data column to be analyzed")
    )
    
    if(length(func) > 0) {
      nbmatrix_object <- add_timestamp(
        nbmatrix_object,
        function_or_method = "nbstat",
        process = paste0("Calculated neighborhood statistic for data column '", data_col,"' with function '", func, "'")
      )
    }
    
    invisible(nbmatrix_object)
    
  }
)

setGeneric(
  "moran",
  function(
    object,
    randomization = TRUE,
    alternative = "greater",
    verbose = FALSE
  ) {
    standardGeneric("moran")
  }
)

setMethod(
  "moran", 
  "nbmatrix", 
  function(
    object,
    randomization = TRUE,
    alternative = "greater",
    verbose = FALSE
  ) { 
    
    data_col <- object@data_col
    
    if(length(data_col) == 0) {
      stop("No definition for data column to be analyzed. Run nbstat() first")
    }
    
    if(isTRUE(verbose)) {
      cat(paste0("Calculating ", spatial_statistics_descriptions$moran, " for data column '", data_col, "' ... "))
    }
    
    polygon_data <- object@polygon_data
    weights <- object@weights
    weights_zeropolicy <- object@config$weights_zeropolicy
    
    moran_results <- 
      moran.test(
        x = polygon_data[[data_col]], 
        listw = weights, 
        randomisation = randomization,
        alternative = alternative,
        zero.policy = weights_zeropolicy
      )
    
    moran <- 
      list(
        I_estimate = moran_results$estimate[[1]],
        I_expected = moran_results$estimate[[2]],
        I_variance = moran_results$estimate[[3]],
        I_p = moran_results$p.value,
        config = list(
          alternative = moran_results$alternative,
          method = moran_results$method,
          randomization = randomization
        )
      )
    
    if(isTRUE(verbose)) {
      cat("OK", "\n")
    }
    
    nbmatrix_object <-
      new(
        "nbmatrix",
        nb = object@nb,
        nbmat = object@nbmat,
        weights = object@weights,
        nbcount = object@nbcount,
        regions_col_name = object@regions_col_name,
        geodata = object@geodata,
        N = object@N,
        N_neighbors = object@N_neighbors,
        polygon_data = object@polygon_data,
        nbmat_data = object@nbmat_data,
        nbmat_data_aggregate = object@nbmat_data_aggregate,
        func = object@func,
        moran = moran,
        getisord = object@getisord,
        gstar = object@gstar,
        data_col = object@data_col,
        config = object@config,
        timestamp = object@timestamp
      )
    
    nbmatrix_object <- add_timestamp(
      nbmatrix_object,
      function_or_method = "moran",
      process = paste0("Calculated ", spatial_statistics_descriptions$moran, " for data column '", data_col, "'")
    )
    
    invisible(nbmatrix_object)
    
  }
)

setGeneric(
  "getisord",
  function(
    object,
    alternative = "greater",
    verbose = FALSE
  ) {
    standardGeneric("getisord")
  }
)

setMethod(
  "getisord",
  "nbmatrix",
  function(
    object,
    alternative = "greater",
    verbose = FALSE
  ) {
    
    data_col <- object@data_col
    
    if(length(data_col) == 0) {
      stop("No definition for data column to be analyzed. Run nbstat() first")
    }
    
    if(isTRUE(verbose)) {
      cat(paste0("Calculating ", spatial_statistics_descriptions$getisord, " for data column '", data_col, "' ... "))
    }
    
    polygon_data <- object@polygon_data
    weights <- object@weights
    weights_zeropolicy <- object@config$weights_zeropolicy
    
    getis_ord_results <- 
      globalG.test(
        x = polygon_data[[data_col]], 
        listw = weights, 
        alternative = alternative,
        zero.policy = weights_zeropolicy
      )
    
    getisord <- 
      list(
        G_estimate = getis_ord_results$estimate[[1]],
        G_expected = getis_ord_results$estimate[[2]],
        G_variance = getis_ord_results$estimate[[3]],
        G_p = getis_ord_results$p.value,
        config = list(
          alternative = getis_ord_results$alternative,
          method = getis_ord_results$method
        )
      )
    
    if(isTRUE(verbose)) {
      cat("OK", "\n")
    }
    
    nbmatrix_object <-
      new(
        "nbmatrix",
        nb = object@nb,
        nbmat = object@nbmat,
        weights = object@weights,
        nbcount = object@nbcount,
        regions_col_name = object@regions_col_name,
        geodata = object@geodata,
        N = object@N,
        N_neighbors = object@N_neighbors,
        polygon_data = object@polygon_data,
        nbmat_data = object@nbmat_data,
        nbmat_data_aggregate = object@nbmat_data_aggregate,
        func = object@func,
        moran = object@moran,
        getisord = getisord,
        gstar = object@gstar,
        data_col = object@data_col,
        config = object@config,
        timestamp = object@timestamp
      )
    
    nbmatrix_object <- add_timestamp(
      nbmatrix_object,
      function_or_method = "getisord",
      process = paste0("Calculated ", spatial_statistics_descriptions$getisord, " for data column '", data_col, "'")
    )
    
    invisible(nbmatrix_object)
    
  }
)

setGeneric(
  "gstar",
  function(
    object,
    alternative = "greater",
    nsim = 999,
    verbose = FALSE
  ) {
    standardGeneric("gstar")
  }
)

setMethod(
  "gstar",
  "nbmatrix",
  function(
    object,
    alternative = "greater",
    nsim = 999,
    verbose = FALSE
  ) {
    
    data_col <- object@data_col
    
    if(length(data_col) == 0) {
      stop("No definition for data column to be analyzed. Run nbstat() first")
    }
    
    if(isTRUE(verbose)) {
      cat(paste0("Calculating ", spatial_statistics_descriptions$gstar, " for data column '", data_col, "' ... "))
    }
    
    polygon_data <- object@polygon_data
    regions_col_name <- object@regions_col_name
    
    weights <- object@weights
    nb <- object@nb
    
    gstar_results <-
      local_gstar_perm(
        x = polygon_data[[data_col]],
        nb = nb,
        wt = weights,
        alternative = alternative,
        nsim = nsim
      )
    
    gstar_results[[regions_col_name]] <- polygon_data[[regions_col_name]]
    
    gstar_results <-
      gstar_results[c(regions_col_name, colnames(gstar_results)[1:length(colnames(gstar_results))-1])]
    
    results_cols_names <- gstar_df_colnames
    
    gstar <- 
      list(
        Gi_estimates = gstar_results,
        results_col_names = results_cols_names,
        config = list(
          alternative = alternative,
          nsim = nsim
        )
      )
    
    if(isTRUE(verbose)) {
      cat("OK", "\n")
    }
    
    nbmatrix_object <-
      new(
        "nbmatrix",
        nb = object@nb,
        nbmat = object@nbmat,
        weights = object@weights,
        nbcount = object@nbcount,
        regions_col_name = object@regions_col_name,
        geodata = object@geodata,
        N = object@N,
        N_neighbors = object@N_neighbors,
        polygon_data = object@polygon_data,
        nbmat_data = object@nbmat_data,
        nbmat_data_aggregate = object@nbmat_data_aggregate,
        func = object@func,
        moran = object@moran,
        getisord = object@getisord,
        gstar = gstar,
        data_col = object@data_col,
        config = object@config,
        timestamp = object@timestamp
      )
    
    nbmatrix_object <- add_timestamp(
      nbmatrix_object,
      function_or_method = "gstar",
      process = paste0("Calculated ", spatial_statistics_descriptions$gstar, " for data column '", data_col, "'")
    )
    
    invisible(nbmatrix_object)
    
  }
)


setMethod(
  "plot",
  "nbmatrix", 
  function(
    x,
    y = NULL,
    statistic = "nbstat",
    # "nbstat" or "gstar"
    attribute = NULL,
    # if statistic = "gstar"
    verbose = FALSE,
    ...
    # further params for sf::plot
  ) {
    
    data_col <- x@data_col
    
    if(length(data_col) == 0) {
      stop("No definition of data column to be analyzed. Run nbstat() first")
    }
    
    geodata <- x@geodata
    regions_col_name <- x@regions_col_name
    
    if(statistic == "nbstat") {
    # Neighborhood statistic (e.g., "sum", "mean")

      func <- x@func
      
      if(length(func) == 0) {
        stop("No definition of function to calculate the nbstat indicator. Run nbstat() first")
      }
      
      data_col_stat <- paste0(data_col, "_", func)
      
      nbmat_data_aggregate <- x@nbmat_data_aggregate
      
      geodata_stat <-
        merge(
          geodata,
          nbmat_data_aggregate,
          by.x = regions_col_name,
          by.y = "ID2"
        )
  
    } else if(statistic == "gstar") {
    # Results of gstar()
      
      if(is.null(attribute)) {
        stop(paste0("Parameter statistic was set to 'gstar'. Parameter attribute must be specified."))
      }
      
      gstar <- x@gstar
      if(length(gstar) == 0) {
        stop(paste0("nbmatrix object does not include calculation of ", spatial_statistics_descriptions[["gstar"]], ". Run gstar() first."))
      }
      
      i <- 0
      data_col_stat <- NULL
      for(i in 1:length(gstar_df_colnames)) {
        if(attribute == gstar_df_colnames[[i]]) {
          data_col_stat <- names(gstar_df_colnames)[i]
        }
      }
      
      if(is.null(data_col_stat)) {
        stop(paste0("The specified attribute '", attribute, "' is not included in the ", spatial_statistics_descriptions[["gstar"]], " data"))
      }
      
      Gi_estimates <- gstar$Gi_estimates
      
      if(!data_col_stat %in% colnames(Gi_estimates)) {
        stop(paste0("The gstar dataframe does not include the column '", data_col_stat, "'"))
      }
      
      regions_col_name <- colnames(Gi_estimates)[1]
      
      geodata_stat <-
        merge(
          geodata,
          Gi_estimates[c(regions_col_name, data_col_stat)],
          by.x = regions_col_name,
          by.y = regions_col_name
        )
      
    } else {
      stop(paste0("Specified parameter statistic = '", statistic, "' is unknown."))
    }
    
    if(isTRUE(verbose)) {
      cat(paste0("Attribute to be plot was set to '", attribute, "' in column '", data_col_stat, "'."), "\n")
    }
    
    plot(
      geodata_stat[data_col_stat],
      ...
    )
    
    invisible(geodata_stat)
    
  }
)


# Function nbmatrix to construct an instance of the nbmatrix class:
nbmatrix <- 
  function(
    polygon_sf, 
    ID_col,
    row.names = NULL,
    snap = NULL,
    queen = TRUE,
    weights_style = "W",
    weights_zeropolicy = NULL,
    verbose = FALSE
  ) {
    
    if (!all(st_geometry_type(polygon_sf) %in% c("POLYGON", "MULTIPOLYGON"))) {
      
      nrow_polygonsf <- nrow(polygon_sf)
      
      polygon_sf <- polygon_sf[st_geometry_type(polygon_sf) %in% c("POLYGON", "MULTIPOLYGON"),]
      
      if(nrow(polygon_sf) == 0) {
        stop("Stated sf object does not contain any polygons or multipolygons.")
      }
      
      if(nrow(polygon_sf) < nrow_polygonsf) {
        warning(paste0("From the stated sf object, ", nrow_polygonsf-nrow(polygon_sf), " non-polygon features were skipped."))
      }
      
    }
    
    geodata_clean <- 
      .clean_geodata(
        polygon_sf,
        unit_col = ID_col
      )
    
    polygon_sf <- geodata_clean$geodata
    N <- geodata_clean$N
    
    N_total <- N[[1]]
    N_with_uid <- N[[2]]
    N_valid <- N[[3]]
    
    if(isTRUE(verbose)) {
      cat(paste0("Constructing neighbors list from ", N_valid, " valid polygons ... "))
    }
    
    nb <- poly2nb(
      polygon_sf, 
      row.names = row.names,
      snap = snap,
      queen = queen
    )
    
    if(isTRUE(verbose)) {
      cat("OK", "\n")
      cat("Calculating spatial weights ... ")
    }
    
    weights <- 
      nb2listw(
        nb, 
        style = weights_style, 
        zero.policy = weights_zeropolicy
      )
    
    if(isTRUE(verbose)) {
      cat("OK", "\n")
      cat("Creating neighborhood matrix ... ")
    }
    
    polygon_sf$ID <- rownames(polygon_sf)
    
    poly_no <- nrow(polygon_sf)
    
    polys_nb <- data.frame(matrix(ncol = 3))
    colnames(polys_nb) <- c("ID", "ID_nb", "nb")
    
    i <- 0
    
    for (i in 1:poly_no) {
      
      neighbors <- nb[[i]]
      
      poly_id <- polygon_sf$ID[i]
      
      poly_nb <- data.frame(rep(poly_id, length(neighbors)), neighbors)
      colnames (poly_nb) <- c("ID", "ID_nb")
      poly_nb$nb <- 1
      
      polys_nb <- rbind(polys_nb, poly_nb)
      
    }
    
    polys_nb <- polys_nb[!is.na(polys_nb$ID),]
    
    polygon_sf_IDs <- cbind(polygon_sf$ID, polygon_sf[[ID_col]])
    colnames(polygon_sf_IDs) <- c("ID", "ID2")
    
    polys_nb <- 
      merge (
        polys_nb, 
        polygon_sf_IDs, 
        by.x = "ID", 
        by.y = "ID"
      )
    
    polys_nb <- 
      merge (
        polys_nb, 
        polygon_sf_IDs, 
        by.x = "ID_nb", 
        by.y = "ID"
      )
    
    polys_nb <- polys_nb[c(2, 4, 1, 5, 3)]  
    
    colnames(polys_nb) <- c("ID", "ID2", "ID_nb", "ID2_nb", "nb")
    
    if(isTRUE(verbose)) {
      cat("OK", "\n")
      cat("Calculating number of neighboring regions ... ")
    }
    
    neighbors_count <- aggregate(
      polys_nb$nb, 
      by=list(polys_nb$ID2), 
      FUN=sum
    )
    
    colnames(neighbors_count) <-
      c(
        "ID2",
        "nb_count"
      )
    
    if(isTRUE(verbose)) {
      cat("OK", "\n")
    }
    
    N <- 
      c(
        N_total,
        # No. of all features
        N_with_uid
        # No. of all features with unique ID
      )
    
    N_neighbors <-
      c(
        nrow(neighbors_count),
        # No. of features with neighbors
        N_neighbors <- nrow(polys_nb),
        # No. of all neighborhoods
        N_neighbors_av <- nrow(polys_nb)/N_with_uid
        # Average no. of neighbors
      )
    
    config <- 
      list(
        row.names = row.names,
        snap = snap,
        queen = queen,
        weights_style = weights_style,
        weights_zeropolicy = weights_zeropolicy
      )
    
    nbmatrix_object <-
      new(
        "nbmatrix",
        nb = nb,
        nbmat = polys_nb,
        weights = weights,
        nbcount = neighbors_count,
        regions_col_name = ID_col,
        geodata = polygon_sf,
        N = N,
        N_neighbors = N_neighbors,
        config = config
      )
    
    nbmatrix_object <- add_timestamp(
      nbmatrix_object,
      function_or_method = "nbmatrix",
      process = paste0("Construction of nbmatrix instance based on ", N_total, " features (", N_with_uid, " with unique ID)")
    )
    
    if(isTRUE(verbose)) {
      paste0("Constructed nbmatrix instance based on ", N_total, " features (", N_with_uid, " with unique ID)")
    }
    
    invisible(nbmatrix_object)
    
  }