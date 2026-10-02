#---------------------------------------------------------------
# Name:        infpan (swash package)
# Purpose:     Class infpan (Infections panel data)
#              Functions, classes and corresponding methods
# Author:      Thomas Wieland 
#              ORCID: 0000-0001-5168-9846
#              mail: geowieland@googlemail.com
# Version:     2.0.0
# Last update: 2026-10-02 19:03
# Copyright (c) 2021-2026 Thomas Wieland
#---------------------------------------------------------------



# Class infpan:
setClass(
  "infpan",
  slots = list(
    input_data = "data.frame",
    data_statistics = "numeric",
    index_col_names = "character",
    cases_col_name = "character",
    other_cols = "character",
    time_format = "character",
    time_unit = "character",
    geodata = "sf",
    timestamp = "list"
  ),
  prototype = list(
    other_cols = character(0),
    time_format = "%Y-%m-%d",
    time_unit = "days",
    geodata = st_sf(
      unit = character(0),
      geometry = st_sfc()
    ),
    timestamp = list()
  )
)

# Methods of class infpan:
setMethod(
  "summary",
  "infpan",
  function(object) {
    
    cat("Infections panel data\n\n")
    
    cat("Input data\n")
    cat(sprintf("  Units       : %s\n", object@data_statistics[1]))
    cat(sprintf("  Time points : %s\n", paste0(object@data_statistics[2], " ", object@time_unit)))
    cat(sprintf(
      "  Balanced    : %s\n",
      ifelse(object@data_statistics[4], "YES", "NO")
    )
    )
    
    cat("Main columns\n")
    cat(sprintf("  Units       : %s\n", object@index_col_names[1]))
    cat(sprintf("  Time points : %s\n", object@index_col_names[2]))
    cat(sprintf("  Cases       : %s\n", object@cases_col_name))
    
    if (length(object@other_cols) > 0) {
      
      cat("Other columns\n")
      
      for (name in names(object@other_cols)) {
        
        if (name %in% permitted_other_cols) {
          cat(sprintf("  %-11s : %s\n", name, object@other_cols[[name]]))
        }
        
      }
    }
    
    if(nrow(object@geodata) > 0) {
      cat("Geodata\n")
      cat(paste0("  EPSG        : ", st_crs(object@geodata)$input), "\n") 
    }
    
    invisible(object)
    
  }
)

setMethod(
  "print", 
  "infpan", 
  function(x) {
    
    cat(paste0("Infections panel data with ", x@data_statistics[1], " spatial units and ", x@data_statistics[2], " time points"), "\n")
    cat ("Use summary() for details")
    
    invisible(x)
    
  }
)

setMethod(
  "show", 
  "infpan", 
  function(object) {
    
    cat(paste0("Infections panel data with ", object@data_statistics[1], " spatial units and ", object@data_statistics[2], " time points"), "\n")
    cat ("Use summary() for details")
    
    invisible(object)
    
  }
)

setGeneric(
  "calculate_Rt", 
  function(
    object, 
    GP = 4,
    correction = FALSE,
    col_name = NULL,
    overwrite = FALSE,
    verbose = FALSE
  ) {
    standardGeneric("calculate_Rt")
  }
)

setMethod(
  "calculate_Rt",
  "infpan",
  function(
    object,
    GP = 4,
    correction = FALSE,
    col_name = NULL,
    overwrite = FALSE,
    verbose = FALSE
  ) {
    
    input_data <- object@input_data
    
    other_cols <- object@other_cols
    
    if ("R_t" %in% names(other_cols)) {
      
      col_R_t <- other_cols[["R_t"]]
      
      if (isFALSE(overwrite)) {
        
        warning(paste0("Effective reproduction number already included in column '", col_R_t, "' and overwrite is set to FALSE"), "\n")
        
        return(invisible(object))
        
      } else {
        
        if(isTRUE(verbose)) {
          message(paste0("Effective reproduction number included in column '", col_R_t, "' will be overwritten"), "\n")  
        }
        
        input_data[col_R_t] <- NULL
        
      }
      
    }
    
    col_region <- object@index_col_names[1]
    col_date <- object@index_col_names[2]
    col_cases <- object@cases_col_name
    
    N <- object@data_statistics[1]
    N_names <- levels(as.factor(input_data[[col_region]]))
    
    if (is.null(col_name)) {
      col_name <- "R_t"
    }
    
    data_Rt <- 
      data.frame(
        matrix(ncol = 2)
      )
    colnames(data_Rt) <- 
      c(
        paste0(col_region, "_", col_date),
        col_name
      )
    
    if (isTRUE(verbose)) {
      cat(paste0("Calculating effective reproduction number from cases in column '", col_cases, "' with GI=", GP, " for ", N, " regions in column '", col_name, "' ... "))
    }
    
    i <- 0
    
    for (i in 1:N) {
      
      input_data_i <- 
        input_data[input_data[[col_region]] == N_names[i],]
      
      input_data_i <- input_data_i[order(input_data_i[[col_region]], input_data_i[[col_date]]),]
      
      input_data_i_Rt <- 
        R_t(
          infections = input_data_i[[col_cases]],
          GP = GP,
          correction = correction
        )
      
      data_Rt_i <- data.frame(
        input_data_i[[paste0(col_region, "_", col_date)]],
        input_data_i_Rt[1]
      )
      colnames(data_Rt_i) <- colnames(data_Rt)
      
      data_Rt <- rbind(
        data_Rt,
        data_Rt_i
      )
      
    }
    
    data_Rt <- data_Rt[!is.na(data_Rt[[paste0(col_region, "_", col_date)]]),]
    
    input_data_Rt <-
      merge(
        input_data,
        data_Rt,
        by.x = paste0(col_region, "_", col_date),
        by.y = paste0(col_region, "_", col_date)
      )
    
    other_cols[["R_t"]] <- col_name
    
    infpan_object <- 
      new(
        "infpan", 
        input_data = data.frame(input_data_Rt),
        data_statistics = object@data_statistics,
        index_col_names = object@index_col_names,
        cases_col_name = object@cases_col_name,
        other_cols = other_cols,
        time_format = object@time_format,
        time_unit = object@time_unit,
        geodata = object@geodata,
        timestamp = object@timestamp
      )
    
    infpan_object <- add_timestamp(
      infpan_object,
      function_or_method = "calculate_Rt(infpan)",
      process = paste0("Calculated effective reproduction number from cases in column '", col_cases, "' with GI=", GP, " for ", N, " regions in column '", col_name, "'")
    )
    
    if (isTRUE(verbose)) {
      cat("OK", "\n")
    }
    
    invisible(infpan_object)
    
  }
)

setGeneric(
  "calculate_rollmean", 
  function(
    object,
    k = 7,
    align = "center",
    fill = NA,
    col_name = NULL,
    overwrite = FALSE,
    verbose = FALSE
  ) {
    standardGeneric("calculate_rollmean")
  }
)

setMethod(
  "calculate_rollmean",
  "infpan",
  function(
    object,
    k = 7,
    align = "center",
    fill = NA,
    col_name = NULL,
    overwrite = FALSE,
    verbose = FALSE
  ) {
    
    input_data <- object@input_data
    
    other_cols <- object@other_cols
    
    if (permitted_other_cols[5] %in% names(other_cols)) {
      
      col_rollmean <- other_cols[[permitted_other_cols[5]]]
      
      if (isFALSE(overwrite)) {
        
        warning(paste0("Rolling mean already included in column '", col_rollmean, "' and overwrite is set to FALSE"), "\n")
        
        return(invisible(object))
        
      } else {
        
        if(isTRUE(verbose)) {
          message(paste0("Rolling mean included in column '", col_rollmean, "' will be overwritten"), "\n")  
        }
        
        input_data[col_rollmean] <- NULL
        
      }
      
    }
    
    col_region <- object@index_col_names[1]
    col_date <- object@index_col_names[2]
    col_cases <- object@cases_col_name
    
    N <- object@data_statistics[1]
    N_names <- levels(as.factor(input_data[[col_region]]))
    
    if (is.null(col_name)) {
      col_name <- paste0(col_cases, "_rm")
    }
    
    data_rollmean <- 
      data.frame(
        matrix(ncol = 2)
      )
    colnames(data_rollmean) <- 
      c(
        paste0(col_region, "_", col_date),
        col_name
      )
    
    if (isTRUE(verbose)) {
      cat(paste0("Calculating rolling mean from cases in column '", col_cases,"' for ", N, " regions with k=", k, " and fill=", fill, " in column '", col_name, "' ... "))
    }
    
    i <- 0
    
    for (i in 1:N) {
      
      input_data_i <- 
        input_data[input_data[[col_region]] == N_names[i],]
      
      input_data_i <- input_data_i[order(input_data_i[[col_region]], input_data_i[[col_date]]),]
      
      input_data_i_rollmean <- 
        rollmean(
          input_data_i[[col_cases]],
          k = k, 
          fill = fill,
          align = align
        )
      
      data_rollmean_i <- data.frame(
        input_data_i[[paste0(col_region, "_", col_date)]],
        input_data_i_rollmean
      )
      colnames(data_rollmean_i) <- colnames(data_rollmean)
      
      data_rollmean <- rbind(
        data_rollmean,
        data_rollmean_i
      )
      
    }
    
    data_rollmean <- data_rollmean[!is.na(data_rollmean[[paste0(col_region, "_", col_date)]]),]
    
    input_data_rollmean <-
      merge(
        input_data,
        data_rollmean,
        by.x = paste0(col_region, "_", col_date),
        by.y = paste0(col_region, "_", col_date)
      )
    
    other_cols <- object@other_cols
    other_cols[[permitted_other_cols[5]]] <- col_name
    
    infpan_object <- 
      new(
        "infpan", 
        input_data = data.frame(input_data_rollmean),
        data_statistics = object@data_statistics,
        index_col_names = object@index_col_names,
        cases_col_name = object@cases_col_name,
        other_cols = other_cols,
        time_format = object@time_format,
        time_unit = object@time_unit,
        geodata = object@geodata,
        timestamp = object@timestamp
      )
    
    infpan_object <- add_timestamp(
      infpan_object,
      function_or_method = "calculate_rollmean(infpan)",
      process = paste0("Calculated rolling mean from cases in column '", col_cases,"' for ", N, " regions with k=", k, " and fill=", fill, " in column '", col_name, "'")
    )
    
    if (isTRUE(verbose)) {
      cat("OK", "\n")
    }
    
    invisible(infpan_object)
    
  }
)

setGeneric(
  "calculate_rollsum", 
  function(
    object,
    k = 7,
    align = "center",
    fill = NA,
    col_name = NULL,
    overwrite = FALSE,
    verbose = FALSE
  ) {
    standardGeneric("calculate_rollsum")
  }
)

setMethod(
  "calculate_rollsum",
  "infpan",
  function(
    object,
    k = 7,
    align = "center",
    fill = NA,
    col_name = NULL,
    overwrite = FALSE,
    verbose = FALSE
  ) {
    
    input_data <- object@input_data
    
    other_cols <- object@other_cols
    
    if (permitted_other_cols[6] %in% names(other_cols)) {
      
      col_rollsum <- other_cols[[permitted_other_cols[6]]]
      
      if (isFALSE(overwrite)) {
        
        warning(paste0("Rolling sum already included in column '", col_rollsum, "' and overwrite is set to FALSE"), "\n")
        
        return(invisible(object))
        
      } else {
        
        if(isTRUE(verbose)) {
          message(paste0("Rolling sum included in column '", col_rollsum, "' will be overwritten"), "\n")  
        }
        
        input_data[col_rollsum] <- NULL
        
      }
      
    }
    
    col_region <- object@index_col_names[1]
    col_date <- object@index_col_names[2]
    col_cases <- object@cases_col_name
    
    N <- object@data_statistics[1]
    N_names <- levels(as.factor(input_data[[col_region]]))
    
    if (is.null(col_name)) {
      col_name <- paste0(col_cases, "_rs")
    }
    
    data_rollsum <- 
      data.frame(
        matrix(ncol = 2)
      )
    colnames(data_rollsum) <- 
      c(
        paste0(col_region, "_", col_date),
        col_name
      )
    
    if (isTRUE(verbose)) {
      cat(paste0("Calculating rolling sum from cases in column '", col_cases,"' for ", N, " regions with k=", k, " and fill=", fill, " in column '", col_name, "' ... "))
    }
    
    i <- 0
    
    for (i in 1:N) {
      
      input_data_i <- 
        input_data[input_data[[col_region]] == N_names[i],]
      
      input_data_i <- input_data_i[order(input_data_i[[col_region]], input_data_i[[col_date]]),]
      
      input_data_i_rollsum <- 
        rollsum(
          input_data_i[[col_cases]],
          k = k, 
          fill = fill,
          align = align
        )
      
      data_rollsum_i <- data.frame(
        input_data_i[[paste0(col_region, "_", col_date)]],
        input_data_i_rollsum
      )
      colnames(data_rollsum_i) <- colnames(data_rollsum)
      
      data_rollsum <- rbind(
        data_rollsum,
        data_rollsum_i
      )
      
    }
    
    data_rollsum <- data_rollsum[!is.na(data_rollsum[[paste0(col_region, "_", col_date)]]),]
    
    input_data_rollsum <-
      merge(
        input_data,
        data_rollsum,
        by.x = paste0(col_region, "_", col_date),
        by.y = paste0(col_region, "_", col_date)
      )
    
    other_cols <- object@other_cols
    other_cols[[permitted_other_cols[6]]] <- col_name
    
    infpan_object <- 
      new(
        "infpan", 
        input_data = data.frame(input_data_rollsum),
        data_statistics = object@data_statistics,
        index_col_names = object@index_col_names,
        cases_col_name = object@cases_col_name,
        other_cols = other_cols,
        time_format = object@time_format,
        time_unit = object@time_unit,
        geodata = object@geodata,
        timestamp = object@timestamp
      )
    
    infpan_object <- add_timestamp(
      infpan_object,
      function_or_method = "calculate_rollsum(infpan)",
      process = paste0("Calculated rolling sum from cases in column '", col_cases,"' for ", N, " regions with k=", k, " and fill=", fill, " in column '", col_name, "'")
    )
    
    if (isTRUE(verbose)) {
      cat("OK", "\n")
    }
    
    invisible(infpan_object)
    
  }
)

setGeneric(
  "calculate_cum", 
  function(
    object,
    col_name = NULL,
    overwrite = FALSE,
    verbose = FALSE
  ) {
    standardGeneric("calculate_cum")
  }
)

setMethod(
  "calculate_cum",
  "infpan",
  function(
    object,
    col_name = NULL,
    overwrite = FALSE,
    verbose = FALSE
  ) {
    
    input_data <- object@input_data
    
    other_cols <- object@other_cols
    
    if ("Cum. cases" %in% names(other_cols)) {
      
      col_cum_cases <- other_cols[["Cum. cases"]]
      
      if (isFALSE(overwrite)) {
        
        warning(paste0("Cumulative infections already included in column '", col_cum_cases, "' and overwrite is set to FALSE"), "\n")
        
        return(invisible(object))
        
      } else {
        
        message(paste0("Cumulative infections included in column '", col_cum_cases, "' will be overwritten"), "\n")
        
        input_data[col_cum_cases] <- NULL
        
      }
      
    }
    
    col_region <- object@index_col_names[1]
    col_date <- object@index_col_names[2]
    col_cases <- object@cases_col_name
    
    N <- object@data_statistics[1]
    N_names <- levels(as.factor(input_data[[col_region]]))
    
    if (is.null(col_name)) {
      col_name <- paste0(col_cases, "_cum")
    }
    
    data_cum <- 
      data.frame(
        matrix(ncol = 2)
      )
    colnames(data_cum) <- 
      c(
        paste0(col_region, "_", col_date),
        col_name
      )
    
    if (isTRUE(verbose)) {
      cat(paste0("Calculating cumulative infections from cases in column '", col_cases, "' for ", N, " regions in column '", col_name, "' ... "))
    }
    
    i <- 0
    
    for (i in 1:N) {
      
      input_data_i <- 
        input_data[input_data[[col_region]] == N_names[i],]
      
      input_data_i <- input_data_i[order(input_data_i[[col_region]], input_data_i[[col_date]]),]
      
      input_data_i_cum <- 
        data.frame(
          cumsum(input_data_i[[col_cases]]), 
          input_data_i[[col_date]]
        )
      colnames(input_data_i_cum) <- c("y", "t")
      
      data_cum_i <- data.frame(
        input_data_i[[paste0(col_region, "_", col_date)]],
        input_data_i_cum$y
      )
      colnames(data_cum_i) <- colnames(data_cum)
      
      data_cum <- rbind(
        data_cum,
        data_cum_i
      )
      
    }
    
    data_cum <- data_cum[!is.na(data_cum[[paste0(col_region, "_", col_date)]]),]
    
    input_data_cum <-
      merge(
        input_data,
        data_cum,
        by.x = paste0(col_region, "_", col_date),
        by.y = paste0(col_region, "_", col_date)
      )
    
    other_cols <- object@other_cols
    other_cols[["Cum. cases"]] <- col_name
    
    infpan_object <- 
      new(
        "infpan", 
        input_data = data.frame(input_data_cum),
        data_statistics = object@data_statistics,
        index_col_names = object@index_col_names,
        cases_col_name = object@cases_col_name,
        other_cols = other_cols,
        time_format = object@time_format,
        time_unit = object@time_unit,
        geodata = object@geodata,
        timestamp = object@timestamp
      )
    
    infpan_object <- add_timestamp(
      infpan_object,
      function_or_method = "calculate_cum(infpan)",
      process = paste0("Calculated cumulative infections from cases in column '", col_cases, "' for ", N, " regions in column '", col_name, "'")
    )
    
    if (isTRUE(verbose)) {
      cat("OK", "\n")
    }
    
    invisible(infpan_object)
    
  }
)

setGeneric(
  "calculate_incidence", 
  function(
    object,
    use_column = NULL,
    col_name = NULL,
    pop_factor = 100000,
    overwrite = FALSE,
    verbose = FALSE
  ) {
    standardGeneric("calculate_incidence")
  }
)

setMethod(
  "calculate_incidence",
  "infpan",
  function(
    object,
    use_column = "Cases",
    col_name = NULL,
    pop_factor = 100000,
    overwrite = FALSE,
    verbose = FALSE
  ) {
    
    input_data <- object@input_data
    
    col_cases <- object@cases_col_name
    other_cols <- object@other_cols
    
    if (!permitted_other_cols[4] %in% names(other_cols)) {
      stop("No population column defined. Calculation of incidence is not possible.")
    }
    
    if (permitted_other_cols[3] %in% names(other_cols)) {
      
      col_incidence <- other_cols[[permitted_other_cols[3]]]
      
      if (isFALSE(overwrite)) {
        
        warning(paste0("Incidence already included in column '", col_incidence, "' and overwrite is set to FALSE"), "\n")
        
        return(invisible(object))
        
      } else {
        
        if(isTRUE(verbose)) {
          message(paste0("Incidence included in column '", col_incidence, "' will be overwritten"), "\n")  
        }
        
        input_data[col_incidence] <- NULL
        
      }
      
    }
    
    col_numerator <- col_cases
    
    if(use_column != "Cases") {
      
      permitted_other_numerator_cols <- permitted_other_cols[c(2,5:6)]
      
      if (use_column %in% permitted_other_numerator_cols) {
        col_numerator <- other_cols[[use_column]]
      } else {
        warning(paste0("Column identifier '", use_column, "' is unknown. Permitted identifier for other columns in incidence calculation are: ", paste(permitted_other_numerator_cols, collapse = ", "), "."), "\n")
      }
      
    }
    
    col_denominator <- other_cols[[permitted_other_cols[4]]]
    
    col_region <- object@index_col_names[1]
    col_date <- object@index_col_names[2]
    
    N <- object@data_statistics[1]
    N_names <- levels(as.factor(input_data[[col_region]]))
    
    if (is.null(col_name)) {
      col_name <- paste0(col_cases, "_inc")
    }
    
    data_incidence <- 
      data.frame(
        matrix(ncol = 2)
      )
    colnames(data_incidence) <- 
      c(
        paste0(col_region, "_", col_date),
        col_name
      )
    
    if (isTRUE(verbose)) {
      cat(paste0("Calculating incidence from ", use_column, " in column '", col_numerator, "' with population in column '", col_denominator, "' for ", N, " regions in column '", col_name, "' ... "))
    }
    
    i <- 0
    
    for (i in 1:N) {
      
      input_data_i <- 
        input_data[input_data[[col_region]] == N_names[i],]
      
      input_data_i <- input_data_i[order(input_data_i[[col_region]], input_data_i[[col_date]]),]
      
      input_data_i_incidence <- 
        input_data_i[[col_numerator]]/input_data_i[[col_denominator]]*pop_factor
      
      data_incidence_i <- data.frame(
        input_data_i[[paste0(col_region, "_", col_date)]],
        input_data_i_incidence
      )
      colnames(data_incidence_i) <- colnames(data_incidence)
      
      data_incidence <- rbind(
        data_incidence,
        data_incidence_i
      )
      
    }
    
    data_incidence <- data_incidence[!is.na(data_incidence[[paste0(col_region, "_", col_date)]]),]
    
    input_data_incidence <-
      merge(
        input_data,
        data_incidence,
        by.x = paste0(col_region, "_", col_date),
        by.y = paste0(col_region, "_", col_date)
      )
    
    other_cols <- object@other_cols
    other_cols[["Incidence"]] <- col_name
    
    infpan_object <- 
      new(
        "infpan", 
        input_data = data.frame(input_data_incidence),
        data_statistics = object@data_statistics,
        index_col_names = object@index_col_names,
        cases_col_name = object@cases_col_name,
        other_cols = other_cols,
        time_format = object@time_format,
        time_unit = object@time_unit,
        geodata = object@geodata,
        timestamp = object@timestamp
      )
    
    infpan_object <- add_timestamp(
      infpan_object,
      function_or_method = "calculate_incidence(infpan)",
      process = paste0("Calculated incidence from ", use_column, " in column '", col_numerator, "' with population in column '", col_denominator, "' for ", N, " regions in column '", col_name, "'")
    )
    
    if (isTRUE(verbose)) {
      cat("OK", "\n")
    }
    
    invisible(infpan_object)
    
  }
)


setGeneric(
  "swash", 
  function(
    object,
    verbose = FALSE
  ) {
    standardGeneric("swash")
  }
)

setMethod(
  "swash",
  "infpan",
  function(
    object,
    verbose = FALSE
  ) { 
    
    infpan_object <- object
    
    sbm_object <- swash_backwash(
      infpan = infpan_object,
      verbose = verbose
    )
    
    invisible(sbm_object)
    
  }
)

setGeneric(
  "growth", 
  function(
    object,
    S_iterations = 10, 
    S_start_est_method = "bisect", 
    seq_by = 10,
    nls = TRUE,
    add_constant = 1,
    overwrite = FALSE,
    verbose = FALSE
  ) {
    standardGeneric("growth")
  }
)

setMethod(
  "growth", 
  "infpan",
  function(
    object,
    S_iterations = 10, 
    S_start_est_method = "bisect", 
    seq_by = 10,
    nls = TRUE,
    add_constant = 1,
    overwrite = FALSE,
    verbose = FALSE
  ) {
    
    col_cases <- object@cases_col_name
    other_cols <- object@other_cols
    
    if ("Cum. cases" %in% names(other_cols)) {
      
      col_cum_cases <- other_cols[["Cum. cases"]]
      
      if (isFALSE(overwrite)) {
        
        message(paste0("Cumulative infections already included in column '", col_cum_cases, "'"), "\n")
        
      } else {
        
        object <- 
          calculate_cum(
            object,
            overwrite = overwrite,
            col_name = col_cum_cases,
            verbose = verbose
          )
        
      }
      
    } else {
      
      message("Cumulative infections not yet included", "\n")
      
      col_cum_cases <- paste0(col_cases, "_cum")
      
      object <- 
        calculate_cum(
          object,
          col_name = col_cum_cases,
          verbose = verbose
        )
      
    }
    
    input_data <- object@input_data
    
    N <- object@data_statistics[1]
    
    col_region <- object@index_col_names[1]
    col_date <- object@index_col_names[2]
    
    time_format <- object@time_format
    
    N_names <- levels(as.factor(input_data[[col_region]]))
    
    logistic_growth_models <- list()
    
    results <- data.frame(matrix(ncol = 20))
    
    for (i in 1:N) {
      
      input_data_i <- 
        input_data[input_data[[col_region]] == N_names[i],]
      
      input_data_i <-
        input_data_i[order(input_data_i[[col_date]]),]
      
      input_data_cum <- 
        data.frame(
          input_data_i[[col_cum_cases]], 
          input_data_i[[col_date]]
        )
      colnames(input_data_cum) <- c("y", "t")
      
      min_date <- min(input_data_cum$t)
      max_date <- max(input_data_cum$t)
      
      logistic_growth_i <-
        logistic_growth(
          y = input_data_cum$y, 
          t = input_data_cum$t, 
          S = max(input_data_cum$y)*1.01,
          S_start = NULL, 
          S_end = NULL, 
          S_iterations = S_iterations, 
          S_start_est_method = S_start_est_method, 
          seq_by = seq_by,
          nls = nls,
          add_constant = add_constant,
          verbose = verbose
        )
      
      results[i,1] <- N_names[i]
      results[i,2] <- min_date
      results[i,3] <- max_date
      results[i,4] <- logistic_growth_i@GrowthModel_OLS$S
      results[i,5] <- logistic_growth_i@GrowthModel_OLS$r
      results[i,6] <- logistic_growth_i@GrowthModel_OLS$y_0
      results[i,7] <- logistic_growth_i@GrowthModel_OLS$ip
      results[i,8] <- logistic_growth_i@GrowthModel_OLS$t_ip
      results[i,9] <- logistic_growth_i@GrowthModel_OLS$fit_metrics$SQR
      results[i,10] <- logistic_growth_i@GrowthModel_OLS$fit_metrics$R2
      results[i,11] <- logistic_growth_i@GrowthModel_OLS$fit_metrics$MAPE
      
      if (isTRUE(logistic_growth_i@config$nls_estimation)) {
        
        results[i,12] <- logistic_growth_i@GrowthModel_NLS$S
        results[i,13] <- logistic_growth_i@GrowthModel_NLS$r
        results[i,14] <- logistic_growth_i@GrowthModel_NLS$y_0
        results[i,15] <- logistic_growth_i@GrowthModel_NLS$ip
        results[i,16] <- logistic_growth_i@GrowthModel_NLS$t_ip
        results[i,17] <- logistic_growth_i@GrowthModel_NLS$fit_metrics$SQR
        results[i,18] <- logistic_growth_i@GrowthModel_NLS$fit_metrics$R2
        results[i,19] <- logistic_growth_i@GrowthModel_NLS$fit_metrics$MAPE
        
      } else {
        
        results[i,20] <- logistic_growth_i@config$nls_error_message 
        
      }
      
      logistic_growth_models <- append(logistic_growth_models, logistic_growth_i)
      
      names(logistic_growth_models)[i] <- N_names[i]
      
    }
    
    colnames(results) <-
      c(
        "region",
        "min_date",
        "max_date",
        "S_OLS",
        "r_OLS",
        "y0_OLS",
        "ip_OLS",
        "tip_OLS",
        "SSQ_OLS",
        "R2_OLS",
        "MAPE_OLS",
        "S_NLS",
        "r_NLS",
        "y0_NLS",
        "ip_NLS",
        "tip_NLS",
        "SSQ_NLS",
        "R2_NLS",
        "MAPE_NLS",
        "nls_error_message"
      )
    
    growthmodels_object <- 
      new(
        "growthmodels", 
        results = results,
        growth_models = logistic_growth_models,
        model_type = "Logistic",
        results_cols = c(
          "region", 
          "r_OLS", 
          "tip_OLS", 
          "R2_OLS", 
          "r_NLS", 
          "tip_NLS", 
          "R2_NLS"
        ),
        results_cols_names = c(
          col_region, 
          "OLS Growth rate",
          "OLS Inflection point",
          "OLS R-squared",
          "NLS Growth rate",
          "NLS Inflection point",
          "NLS R-squared"
        ),
        data_statistics = object@data_statistics,
        time_format = object@time_format
      )
    
    growthmodels_object <- add_timestamp(
      growthmodels_object,
      function_or_method = "growth",
      process = paste0("Estimation of ", nrow(results), " logistic growth models and creation of growthmodels object")
    )
    
    invisible(growthmodels_object)
    
  }
)

setGeneric(
  "growth_initial", 
  function(
    object,
    time_units = 10,
    GI = 4,
    nls = TRUE,
    nls_start = list(a = 1, b = 0.1),
    add_constant = 1,
    verbose = FALSE
  ) {
    standardGeneric("growth_initial")
  }
)

setMethod(
  "growth_initial",
  "infpan", 
  function(
    object,
    time_units = 10,
    GI = 4,
    nls = TRUE,
    nls_start = list(a = 1, b = 0.1),
    add_constant = 1,
    verbose = FALSE
  ) {
    
    N <- object@data_statistics[1]
    
    col_cases <- object@cases_col_name
    col_region <- object@index_col_names[1]
    col_date <- object@index_col_names[2]
    
    time_format <- object@time_format
    
    input_data <- object@input_data
    
    N_names <- levels(as.factor(input_data[[col_region]]))
    
    exponential_growth_models <- list()
    
    results <- data.frame(matrix(ncol = 16))
    
    for (i in 1:N) {
      
      input_data_i <- 
        input_data[input_data[[col_region]] == N_names[i],]
      
      if (time_units > nrow(input_data_i)) {
        stop(paste0("Dataset of region ", N_names[i], " contains ", nrow(input_data_i), " time units but ", time_units, " were stated."))  
      }
      
      input_data_i <-
        input_data_i[order(input_data_i[[col_date]]),]
      
      input_data_cases <- 
        data.frame(
          input_data_i[[col_cases]], 
          input_data_i[[col_date]]
        )
      
      input_data_cases <- input_data_cases[1:time_units,]
      colnames(input_data_cases) <- c("y", "t")
      
      if (nrow(input_data_cases) > 1) {
        
        min_date <- min(input_data_cases$t)
        max_date <- max(input_data_cases$t)
        
        exponential_growth_i <-
          exponential_growth(
            y = input_data_cases$y, 
            t = input_data_cases$t, 
            GI = GI,
            nls = nls,
            nls_start = nls_start,
            add_constant = add_constant,
            verbose = verbose
          )
        
        results[i,1] <- N_names[i]
        results[i,2] <- as.Date(min_date, format = time_format)
        results[i,3] <- as.Date(max_date, format = time_format)
        results[i,4] <- exponential_growth_i@GrowthModel_OLS$exp_gr
        results[i,5] <- exponential_growth_i@GrowthModel_OLS$R0
        results[i,6] <- exponential_growth_i@GrowthModel_OLS$doubling
        results[i,7] <- exponential_growth_i@GrowthModel_OLS$fit_metrics$SQR
        results[i,8] <- exponential_growth_i@GrowthModel_OLS$fit_metrics$R2
        results[i,9] <- exponential_growth_i@GrowthModel_OLS$fit_metrics$MAPE
        
        if (isTRUE(exponential_growth_i@config$nls_estimation)) {
          
          results[i,10] <- exponential_growth_i@GrowthModel_NLS$exp_gr
          results[i,11] <- exponential_growth_i@GrowthModel_NLS$R0
          results[i,12] <- exponential_growth_i@GrowthModel_NLS$doubling
          results[i,13] <- exponential_growth_i@GrowthModel_NLS$fit_metrics$SQR
          results[i,14] <- exponential_growth_i@GrowthModel_NLS$fit_metrics$R2
          results[i,15] <- exponential_growth_i@GrowthModel_NLS$fit_metrics$MAPE
          
        } else {
          
          results[i,16] <- exponential_growth_i@config$nls_error_message 
          
        }
        
        exponential_growth_models <- 
          append(
            exponential_growth_models, 
            exponential_growth_i
          )
        
        exponential_growth_models[[N_names[i]]] <- exponential_growth_i
        
      } else {
        
        warning(paste0("No cases for region ", N_names[i], ". No exponential growth model built."))
        
        results[i,1] <- N_names[i]
        results[i,2:15] <- NA
        
      }
      
    }
    
    colnames(results) <-
      c(
        "region",
        "min_date",
        "max_date",
        "r_OLS",
        "R_0_OLS",
        "doubling_rate_OLS",
        "SSQ_OLS",
        "R2_OLS",
        "MAPE_OLS",
        "r_NLS",
        "R_0_NLS",
        "doubling_rate_NLS",
        "SSQ_NLS",
        "R2_NLS",
        "MAPE_NLS"
      )
    
    growthmodels_object <- 
      new(
        "growthmodels", 
        results = results,
        growth_models = exponential_growth_models,
        model_type = "Exponential",
        results_cols = c(
          "region", 
          "r_OLS", 
          "R_0_OLS", 
          "doubling_rate_OLS",
          "R2_OLS",
          "r_NLS", 
          "R_0_NLS", 
          "doubling_rate_NLS",
          "R2_NLS"
        ),
        results_cols_names = c(
          col_region, 
          "OLS Growth rate",
          "OLS Basic reproduction number",
          "OLS Doubling rate",
          "OLS R-squared",
          "NLS Growth rate",
          "NLS Basic reproduction number",
          "NLS Doubling rate",
          "NLS R-squared"
        ),
        data_statistics = object@data_statistics,
        time_format = object@time_format
      )
    
    growthmodels_object <- add_timestamp(
      growthmodels_object,
      function_or_method = "growth_initial",
      process = "Analysis of exponential growth and creation of growthmodels object"
    )
    
    invisible(growthmodels_object)
    
  }
)


setGeneric(
  "growth_breaks", 
  function(
    object,
    ln = FALSE,
    add_constant = 1,
    alpha = 0.05,
    verbose = FALSE
  ) {
    standardGeneric("growth_breaks")
  }
)

setMethod(
  "growth_breaks",
  "infpan", 
  function(
    object,
    ln = FALSE,
    add_constant = 1,
    alpha = 0.05,
    verbose = FALSE
  ) { 
    
    N <- object@data_statistics[1]
    
    col_cases <- object@cases_col_name
    col_region <- object@index_col_names[1]
    col_date <- object@index_col_names[2]
    
    time_format <- object@time_format
    
    input_data <- object@input_data
    
    N_names <- levels(as.factor(input_data[[col_region]]))
    
    breaks_growth_models <- list()
    
    results <- data.frame(matrix(ncol = 21))
    colnames(results) <-
      c(
        "region",
        "min_date",
        "max_date",
        "first", 
        "last", 
        "intercept", 
        "intercept_CI25", 
        "intercept_CI975", 
        "intercept_p", 
        "slope", 
        "slope_CI25", 
        "slope_CI975", 
        "slope_p",
        "SQR",
        "SAR",
        "SQT",
        "R2",
        "MSE",
        "RMSE",
        "MAE",
        "MAPE"
      )
    
    for (i in 1:N) {
      
      input_data_i <- 
        input_data[input_data[[col_region]] == N_names[i],]
      
      input_data_i <-
        input_data_i[order(input_data_i[[col_date]]),]
      
      input_data_cases <- 
        data.frame(
          input_data_i[[col_cases]], 
          input_data_i[[col_date]]
        )
      
      colnames(input_data_cases) <- c("y", "t")
      
      if (nrow(input_data_cases) > 1) {
        
        min_date <- min(input_data_cases$t)
        max_date <- max(input_data_cases$t)
        
        breaks_growth_i <-
          breaks_growth(
            y = input_data_cases$y, 
            t = input_data_cases$t,
            ln = ln,
            add_constant = add_constant,
            alpha = alpha,
            verbose = verbose
          )
        
        results_i <- data.frame(breaks_growth_i@GrowthModel_OLS$segments_models_df)
        
        results_i$region <- N_names[i]
        results_i$min_date <- as.Date(min_date, format = time_format)
        results_i$max_date <- as.Date(max_date, format = time_format)
        
        results_i <- results_i[c(19:21, 1:18)]
        
        results <-
          rbind(
            results,
            results_i
          )
        
        breaks_growth_models <- 
          append(
            breaks_growth_models, 
            breaks_growth_i
          )
        
        breaks_growth_models[[N_names[i]]] <- breaks_growth_i
        
      } else {
        
        warning(paste0("No cases for region ", N_names[i], ". No breaking points model built."))
        
        results[i,1] <- N_names[i]
        results[i,2:21] <- NA
        
      }
      
    }
    
    colnames(results) <-
      c(
        "region",
        "min_date",
        "max_date",
        "first", 
        "last", 
        "intercept", 
        "intercept_CI25", 
        "intercept_CI975", 
        "intercept_p", 
        "slope", 
        "slope_CI25", 
        "slope_CI975", 
        "slope_p",
        "SQR",
        "SAR",
        "SQT",
        "R2",
        "MSE",
        "RMSE",
        "MAE",
        "MAPE"
      )
    
    results <- results[!is.na(results$region),]
    
    growthmodels_object <- 
      new(
        "growthmodels", 
        results = results,
        growth_models = breaks_growth_models,
        model_type = "Breaking Points",
        results_cols = c(
          "region", 
          "first", 
          "last", 
          "intercept", 
          "intercept_p", 
          "slope", 
          "slope_p",
          "R2"
        ),
        results_cols_names = c(
          col_region, 
          "First of segment",
          "Last of segment",
          "Intercept",
          "Intercept p value",
          "Slope",
          "Slope p value",
          "R-squared"
        ),
        data_statistics = object@data_statistics,
        time_format = object@time_format
      )
    
    growthmodels_object <- add_timestamp(
      growthmodels_object,
      function_or_method = "growth_breaks",
      process = "Analysis of breaking points and creation of growthmodels object"
    )
    
    invisible(growthmodels_object)
    
  }
)


setGeneric(
  "growth_hawkes", 
  function(
    object,
    optim_method = "L-BFGS-B",
    verbose = FALSE
  ) {
    standardGeneric("growth_hawkes")
  }
)

setMethod(
  "growth_hawkes",
  "infpan", 
  function(
    object,
    optim_method = "L-BFGS-B",
    verbose = FALSE
  ) {
    
    N <- object@data_statistics[1]
    
    col_cases <- object@cases_col_name
    col_region <- object@index_col_names[1]
    col_date <- object@index_col_names[2]
    
    time_format <- object@time_format
    
    input_data <- object@input_data
    
    N_names <- levels(as.factor(input_data[[col_region]]))
    
    hawkes_growth_models <- list()
    
    results <- data.frame(matrix(ncol = 7))
    
    for (i in 1:N) {
      
      input_data_i <- 
        input_data[input_data[[col_region]] == N_names[i],]
      
      input_data_i <-
        input_data_i[order(input_data_i[[col_date]]),]
      
      input_data_cases <- 
        data.frame(
          input_data_i[[col_cases]], 
          input_data_i[[col_date]]
        )
      
      colnames(input_data_cases) <- c("y", "t")
      
      if (nrow(input_data_cases) > 1) {
        
        min_date <- min(input_data_cases$t)
        max_date <- max(input_data_cases$t)
        
        hawkes_growth_i <-
          hawkes_growth(
            y = input_data_cases$y,
            verbose = verbose
          )
        
        results[i,1] <- N_names[i]
        results[i,2] <- as.Date(min_date, format = time_format)
        results[i,3] <- as.Date(max_date, format = time_format)
        results[i,4] <- hawkes_growth_i@mu
        results[i,5] <- hawkes_growth_i@alpha
        results[i,6] <- hawkes_growth_i@beta
        results[i,7] <- hawkes_growth_i@br
        
        hawkes_growth_models <- 
          append(
            hawkes_growth_models, 
            hawkes_growth_i
          )
        
        hawkes_growth_models[[N_names[i]]] <- hawkes_growth_i
        
      } else {
        
        warning(paste0("No cases for region ", N_names[i], ". No Hawkes growth model built."))
        
        results[i,1] <- N_names[i]
        results[i,2] <- NA
        results[i,3] <- NA
        results[i,4] <- NA
        results[i,5] <- NA
        results[i,6] <- NA
        results[i,7] <- NA
        
      }
      
    }
    
    colnames(results) <-
      c(
        "region",
        "min_date",
        "max_date",
        "mu",
        "alpha",
        "beta",
        "br"
      )
    
    growthmodels_object <- 
      new(
        "growthmodels", 
        results = results,
        growth_models = hawkes_growth_models,
        model_type = model_descriptions[["hawkes"]],
        results_cols = c(
          "region", 
          "mu",
          "alpha",
          "beta",
          "br"
        ),
        results_cols_names = c(
          col_region, 
          "Baseline",
          "Excitation",
          "Decay",
          "Breaking ratio"
        ),
        data_statistics = object@data_statistics,
        time_format = object@time_format
      )
    
    growthmodels_object <- add_timestamp(
      growthmodels_object,
      function_or_method = "growth_hawkes",
      process = "Analysis of Hawkes growth and creation of growthmodels object"
    )
    
    invisible(growthmodels_object)
    
  }
)

setMethod(
  "plot", 
  signature(x = "infpan"),
  function(
    x,
    y = NULL,
    col = "red",
    lty = "solid",
    scale = FALSE,
    normalize_by_col = NULL,
    normalize_factor = 1,
    plot_rollmean = FALSE,
    rollmean_col = "blue",
    rollmean_lty = "solid",
    rollmean_k = 7,
    rollmean_align = "center",
    rollmean_fill = NA,
    growth_col = "orange",
    growth_lty = "solid",
    growth_per_time_unit = 1
  ) {
    
    object <- x
    
    par_old <- par(no.readonly = TRUE)
    on.exit(par(par_old))
    dev.new()
    
    N <- object@data_statistics[1]
    
    plot_cols <- 4
    plot_rows <- ceiling(N/plot_cols)
    par (mfrow = c(plot_rows, plot_cols)) 
    
    input_data <- object@input_data
    
    col_cases <- object@cases_col_name
    col_region <- object@index_col_names[1]
    col_date <- object@index_col_names[2]
    
    N_names <- levels(as.factor(input_data[[col_region]]))
    
    if (!is.null(normalize_by_col)) {
      
      input_data[[paste0(col_cases, "_normalized")]] <-
        input_data[[col_cases]]/input_data[[normalize_by_col]]*normalize_factor
      
      y_max <- max(input_data[[paste0(col_cases, "_normalized")]])*1.05
      
    } else {
      
      y_max <- max(input_data[[col_cases]])*1.05
      
    }
    
    i <- 0
    
    for (i in 1:N) {
      
      par(mar = c(2, 2, 1, 1))
      
      input_data_i <- 
        input_data[input_data[[col_region]] == N_names[i],]
      input_data_i <-
        input_data_i[order(input_data_i[[col_date]]),]
      
      if (!is.null(normalize_by_col)) {
        
        if (scale == FALSE) {
          y_max <- max(input_data_i[[paste0(col_cases, "_normalized")]])*1.05
        }
        
        plot(
          input_data_i[[col_date]], 
          input_data_i[[paste0(col_cases, "_normalized")]],
          type = "l",
          col = col,
          lty = lty,
          main = N_names[i],
          ylim = c(0, y_max)
        )
        
        if (isTRUE(plot_rollmean)) {
          
          input_data_i[[paste0(col_cases, "_normalized_rm")]] <- 
            rollmean(
              input_data_i[[paste0(col_cases, "_normalized")]],
              k = rollmean_k, 
              fill = rollmean_fill,
              align = rollmean_align
            )
          
          lines (
            input_data_i[[col_date]],
            input_data_i[[paste0(col_cases, "_normalized_rm")]], 
            col = rollmean_col, 
            lty = rollmean_lty
          )
          
        }
        
      }
      
      else {
        
        if (scale == FALSE) {
          y_max <- max(input_data_i[[col_cases]])*1.05
        }
        
        plot(
          input_data_i[[col_date]], 
          input_data_i[[col_cases]],
          col = col, 
          main = N_names[i],
          type = "l",
          ylim = c(0, y_max)
        )
        
        if (isTRUE(plot_rollmean)) {
          
          input_data_i[[paste0(col_cases, "_rm")]] <- 
            rollmean(
              input_data_i[[col_cases]],
              k = rollmean_k, 
              fill = rollmean_fill,
              align = rollmean_align
            )
          
          lines (
            input_data_i[[col_date]],
            input_data_i[[paste0(col_cases, "_rm")]], 
            col = rollmean_col, 
            lty = rollmean_lty
          )
          
        }
        
      }
      
    }
    
    invisible(object)
  }
)


setGeneric(
  "add_geodata", 
  function(
    object,
    sf_geodata,
    unit_col,
    verbose = FALSE
  ) {
    standardGeneric("add_geodata")
  }
)

setMethod(
  "add_geodata",
  "infpan", 
  function(
    object,
    sf_geodata,
    unit_col,
    verbose = FALSE
  ) { 
    
    geodata_clean <- 
      .clean_geodata(
        sf_geodata,
        unit_col = unit_col
      )
    
    geodata <- geodata_clean[[1]]
    N <- geodata_clean[[2]]
    geodata_crs <- geodata_clean[[3]]
    
    if(isTRUE(verbose)) {
      cat(paste0("Added geodata (CRS: ", geodata_crs, ") with ", N[[1]], " features (", N[[2]], " with valid unique ID, ", N[[3]], " with valid geometry) to infpan object"), "\n")
    }
    
    infpan_object <- 
      new(
        "infpan", 
        input_data = object@input_data,
        data_statistics = object@data_statistics,
        index_col_names = object@index_col_names,
        cases_col_name = object@cases_col_name,
        other_cols = object@other_cols,
        time_format = object@time_format,
        time_unit = object@time_unit,
        geodata = geodata,
        timestamp = object@timestamp
      )
    
    infpan_object <- add_timestamp(
      infpan_object,
      function_or_method = "add_geodata(infpan)",
      process = paste0("Added geodata (CRS: ", geodata_crs, ") with ", N[[1]], " features (", N[[2]], " with valid unique ID, ", N[[3]], " with valid geometry) to infpan object")
    )
    
    invisible(infpan_object)
    
  }
)

setGeneric(
  "plot_map", 
  function(
    object,
    attribute = "Cases",
    timepoint = NULL,
    verbose = FALSE,
    ...
  ) {
    standardGeneric("plot_map")
  }
)

setMethod(
  "plot_map",
  "infpan", 
  function(
    object,
    attribute = "Cases",
    timepoint = NULL,
    verbose = FALSE,
    ...
  ) { 
    
    geodata <- object@geodata
    
    if(nrow(geodata) == 0) {
      stop("Infpan object does not include valid geodata. Use add_geodata() first")
    }
    
    geodata_region_col <- colnames(geodata)[1]
    
    cases_col_name <- object@cases_col_name
    other_cols <- object@other_cols
    
    input_data <- object@input_data
    
    col_date <- object@index_col_names[2]
    col_region <- object@index_col_names[1]
    
    if(is.null(timepoint)) {
      
      timepoint <- max(input_data[[col_date]])
      
      if(isTRUE(verbose)) {
        cat(paste0("No timepoint was specified. Timepoint for map plotting is set to last value = '", timepoint, "'"), "\n")
      }
      
    }
    
    input_data_timepoint <- input_data[input_data[[col_date]] == timepoint,]
    
    if(nrow(input_data_timepoint) == 0) {
      stop(paste0("Specified timepoint '", timepoint, "' is not included in infpan data."))
    }
    
    if(attribute == "Cases") {
      
      attribute_col <- cases_col_name
      
    } else {
      
      if (!attribute %in% permitted_cols) {
        stop(paste0("Specified attribute '", attribute, "' is unknown. Permitted columns are: ", paste(permitted_cols, collapse = ", "), "."))
      }
      
      attribute_col <- other_cols[[attribute]]
      
    }
    
    if(isTRUE(verbose)) {
      cat(paste0("Attribute to be plot was set to '", attribute, "' in column '", attribute_col, "'."), "\n")
    }
    
    geodata_with_attr <-
      merge(
        geodata,
        input_data_timepoint[c(col_region, attribute_col)],
        by.x = geodata_region_col,
        by.y = col_region
      )
    
    plot(
      geodata_with_attr[attribute_col],
      ...
    )
    
    invisible(object)
    
  }
)


setGeneric(
  "spatial_statistic", 
  function(
    object,
    attribute = "Cases",
    timepoint = NULL,
    statistic = "nbstat",
    func = NULL,
    row.names = NULL,
    snap = NULL,
    queen = TRUE,
    weights_style = "W",
    weights_zeropolicy = NULL,
    verbose = FALSE,
    ...
  ) {
    standardGeneric("spatial_statistic")
  }
)

setMethod(
  "spatial_statistic",
  "infpan", 
  function(
    object,
    attribute = "Cases",
    timepoint = NULL,
    statistic = "nbstat",
    func = NULL,
    row.names = NULL,
    snap = NULL,
    queen = TRUE,
    weights_style = "W",
    weights_zeropolicy = NULL,
    verbose = FALSE,
    ...
  ) {
    
    if(!statistic %in% names(spatial_statistics_descriptions)) {
      stop(paste0("Specified statistic '", statistic, "' is unknown. Permitted values are: ", paste(names(spatial_statistics_descriptions), collapse = ", "), "."))
    }
    
    geodata <- object@geodata
    
    if(nrow(geodata) == 0) {
      stop("Infpan object does not include valid geodata. Use add_geodata() first")
    }
    
    geodata_region_col <- colnames(geodata)[1]
    
    nbmatrix_object <- nbmatrix(
      polygon_sf = geodata, 
      ID_col = geodata_region_col,
      row.names = row.names,
      snap = snap,
      queen = queen,
      weights_style = weights_style,
      weights_zeropolicy = weights_zeropolicy,
      verbose = verbose
    )
    
    cases_col_name <- object@cases_col_name
    other_cols <- object@other_cols
    
    input_data <- object@input_data
    
    col_date <- object@index_col_names[2]
    col_region <- object@index_col_names[1]
    
    if(is.null(timepoint)) {
      
      timepoint <- max(input_data[[col_date]])
      
      if(isTRUE(verbose)) {
        cat(paste0("No timepoint was specified. Timepoint for map plotting is set to last value = '", timepoint, "'"), "\n")
      }
      
    }
    
    input_data_timepoint <- input_data[input_data[[col_date]] == timepoint,]
    
    if(nrow(input_data_timepoint) == 0) {
      stop(paste0("Specified timepoint '", timepoint, "' is not included in infpan data."))
    }
    
    if(attribute == "Cases") {
      
      attribute_col <- cases_col_name
      
    } else {
      
      if (!attribute %in% permitted_cols) {
        stop(paste0("Specified attribute '", attribute, "' is unknown. Permitted columns are: ", paste(permitted_cols, collapse = ", "), "."))
      }
      
      attribute_col <- other_cols[[attribute]]
      
    }
    
    if(isTRUE(verbose)) {
      cat(paste0("Attribute for calculating the spatial statistic was set to '", attribute, "' in column '", attribute_col, "'."), "\n")
    }

    nbmatrix_object <- nbstat(
      nbmatrix_object,
      link_data = input_data_timepoint, 
      ID_col = col_region, 
      data_col = attribute_col, 
      verbose = verbose
    )
        
    if(statistic == "moran") {
      nbmatrix_object <- moran(
        nbmatrix_object,
        verbose = verbose,
        ...
      )
    } else if(statistic == "getisord") {
      nbmatrix_object <- getisord(
        nbmatrix_object,
        verbose = TRUE,
        ...
      )
    } else if(statistic == "gstar") {
      nbmatrix_object <- gstar(
        nbmatrix_object,
        verbose = TRUE,
        ...
      )
    }
    
  invisible(nbmatrix_object)
    
  }
)


# Loading infections panel data (data.frame)
# and creating an object of class infpan:
load_infections_paneldata <-
  function(
    data,
    col_cases, 
    col_date, 
    col_region,
    other_cols = NULL,
    time_format = "%Y-%m-%d",
    time_unit = "days",
    verbose = FALSE
  ) {
    
    col_errors <- character(0)
    if (!col_cases %in% colnames(data))
    {
      col_errors <- c(col_errors, col_cases)
    }
    if (!col_date %in% colnames(data))
    {
      col_errors <- c(col_errors, col_date)
    }
    if (!col_region %in% colnames(data))
    {
      col_errors <- c(col_errors, col_region)
    }
    if (length(col_errors) > 0) {
      stop(paste0("Import failed. Columns not in data.frame: ", paste(col_errors, collapse = ", ")))
    }
    
    data[[col_date]] <- as.Date(data[[col_date]], format = time_format)
    
    data <- data[order(data[[col_region]], data[[col_date]]),]
    
    N <- nlevels(as.factor(data[[col_region]]))
    N_names <- levels(as.factor(data[[col_region]]))
    N_withoutcases <- 0
    TP <- nlevels(as.factor(data[[col_date]]))
    TP_t <- levels(as.factor(data[[col_date]]))
    
    data_check_balanced <- 
      is_balanced(
        data,
        col_cases = col_cases,
        col_region = col_region,
        col_date = col_date
      )
    data_balanced <- data_check_balanced$data_balanced
    
    data_statistics <- c(N, TP, N_withoutcases, data_balanced)
    
    if (isTRUE(verbose)) {
      
      cat(paste0("Infections panel data include ", N, " regions in column '", col_region, "' and ", TP, " time points in column '", col_date, "', with cases in column '", col_cases, "'. "))
      
      if (isTRUE(data_balanced)) {
        cat("The data is balanced.")
      } else {
        cat("The data is not balanced.")
      }
      cat("\n")
      
    }
    
    if (!is.null(other_cols)) {
      
      for (name in names(other_cols)) {
        
        if (name %in% permitted_other_cols) {
          other_cols <- c(other_cols, name)
        } else {
          warning(paste0("Column identifier '", name, "' is unknown. Permitted identifier for other columns are: ", paste(permitted_other_cols, collapse = ", "), "."), "\n")
        }
      }
      
    } else {
      other_cols <- character(0)
    }
    
    data[paste0(col_region, "_", col_date)] <-
      paste0(
        as.character(data[[col_region]]),
        "_",
        as.character(data[[col_date]])
      )
    
    infpan_object <- 
      new(
        "infpan", 
        input_data = data.frame(data),
        data_statistics = data_statistics,
        index_col_names = c(col_region, col_date),
        cases_col_name = col_cases,
        other_cols = other_cols,
        time_format = time_format,
        time_unit = time_unit
      )
    
    infpan_object <- add_timestamp(
      infpan_object,
      function_or_method = "load_infections_paneldata",
      process = "Import of infections panel data and creation of infpan object"
    )
    
    return(infpan_object)
    
  }
