#---------------------------------------------------------------
# Name:        helper (swash package)
# Purpose:     Helper functions for the swash package
# Author:      Thomas Wieland 
#              ORCID: 0000-0001-5168-9846
#              mail: geowieland@googlemail.com
# Version:     1.1.0
# Last update: 2026-10-13 21:29
# Copyright (c) 2025-2026 Thomas Wieland
#---------------------------------------------------------------


is_balanced <-
  function (
    data,
    col_cases, 
    col_date, 
    col_region,
    as_balanced = TRUE,
    fill_missing = 0,
    drop_duplicates = TRUE
  ) {
    
    if(!"data.frame" %in% class(data)) {
      stop("Specified object is not a data.frame")
    }
    
    required_cols <- c(col_cases, col_date, col_region)
    missing_cols <- setdiff(required_cols, colnames(data))
    if (length(missing_cols) > 0) {
      stop(
        paste("The following columns are missing:", paste(missing_cols, collapse = ", ")),
        call. = FALSE
      )
    }
    
    if(isTRUE(drop_duplicates)) {
      
      nrow_data <- nrow(data)
      
      data <- data[!duplicated(data[c(col_region, col_date)]), , drop = FALSE]
      
      if(nrow(data) < nrow_data) {
        warning(paste0("User-stated data frame included ", nrow_data-nrow(data), " rows with duplicated region/time combinations, which were excluded"))  
      }
      
    }
    
    N <- nlevels(as.factor(data[[col_region]]))
    TP <- nlevels(as.factor(data[[col_date]]))
    
    if (nrow(data) != (TP*N)) { 
      
      data_balanced <- FALSE
      
    } else {
      
      if (((length(unique(table(data[[col_date]]))) == 1) == FALSE) |
          (length(unique(table(data[[col_region]]))) == 1) == FALSE) {
        
        data_balanced <- FALSE
        
      } else {
        
        if (any (is.na(data[[col_cases]]))) {
          
          data_balanced <- FALSE
          
        } else {
          
          data_balanced <- TRUE
          
        }
      }
      
    }
    
    if ((data_balanced == FALSE) & (as_balanced == TRUE)) {
      
      data <-
        as_balanced(
          data,
          col_cases, 
          col_date, 
          col_region,
          fill_missing = fill_missing
        )
      
      data_balanced <- TRUE
    }
    
    results <- 
      list(
        data_balanced = data_balanced, 
        data = data
      )
    
    return (results)
    
  }


as_balanced <-
  function(
    data,
    col_cases, 
    col_date, 
    col_region,
    fill_missing = 0,
    drop_duplicates = TRUE
  ) {
    
    if(!"data.frame" %in% class(data)) {
      stop("Specified object is not a data.frame")
    }
    
    required_cols <- c(col_cases, col_date, col_region)
    missing_cols <- setdiff(required_cols, colnames(data))
    if (length(missing_cols) > 0) {
      stop(
        paste("The following columns are missing:", paste(missing_cols, collapse = ", ")),
        call. = FALSE
      )
    }
    
    if(isTRUE(drop_duplicates)) {
      
      nrow_data <- nrow(data)
      
      data <- data[!duplicated(data[c(col_region, col_date)]), , drop = FALSE]
      
      if(nrow(data) < nrow_data) {
        warning(paste0("User-stated data frame included ", nrow_data-nrow(data), " rows with duplicated region/time combinations, which were excluded"))  
      }
      
    }
    
    N <- nlevels(as.factor(data[[col_region]]))
    TP <- nlevels(as.factor(data[[col_date]]))
    N_names <- as.character(levels(as.factor(data[[col_region]])))
    TP_t <- as.character(levels(as.factor(data[[col_date]])))
    
    N_x_TPt <- merge (N_names, TP_t)
    colnames(N_x_TPt) <- c(paste0("__", col_region, "__"), paste0("__", col_date, "__"))
    N_x_TPt[[paste0("__", col_region, "_x_", col_date, "__")]] <-
      paste0(N_x_TPt[[paste0("__", col_region, "__")]], "_x_", N_x_TPt[[paste0("__", col_date, "__")]])
    
    data[[paste0("__", col_region, "_x_", col_date, "__")]] <-
      paste0(data[[col_region]], "_x_", data[[col_date]])
    
    data <-
      merge (
        data,
        N_x_TPt,
        by.x = paste0("__", col_region, "_x_", col_date, "__"),
        by.y = paste0("__", col_region, "_x_", col_date, "__"),
        all.y = TRUE
      )
    
    data[[col_region]] <- data[[paste0("__", col_region, "__")]]
    data[[col_date]] <- data[[paste0("__", col_date, "__")]]
    if (nrow(data[is.na(data[[col_cases]]),]) > 0) {
      data[is.na(data[[col_cases]]),][[col_cases]] <- fill_missing
    }
    
    data[[paste0("__", col_region, "_x_", col_date, "__")]] <- NULL
    data[[paste0("__", col_region, "__")]] <- NULL
    data[[paste0("__", col_date, "__")]] <- NULL
    
    return(data)
    
  }


setGeneric(
  "add_timestamp", 
  function(
    object, 
    function_or_method = "", 
    process = "",
    status = "OK"
  ) {
    standardGeneric("add_timestamp")
  }
)

setMethod(
  "add_timestamp",
  "ANY",
  function(
    object, 
    function_or_method = "", 
    process = "",
    status = "OK"
  ) {
    
    timestamp <- c(
      time = as.character(Sys.time()),
      package = paste0(package_name, " ", package_version),
      function_or_method = function_or_method,
      process = process,
      status = status
    )
    
    object@timestamp <- append(object@timestamp, list(timestamp))
    
    validObject(object)
    
    object
  }
)

setGeneric(
  "timestamps", 
  function(
    object
  ) {
    standardGeneric("timestamps")
  }
)

setMethod(
  "timestamps",
  "ANY",
  function(
    object
  ) {
    
    timestamp <- object@timestamp
    
    for (i in seq_along(timestamp)) {
      
      entry <- timestamp[[i]]
      
      cat(
        paste0(
          i, 
          " ", 
          entry[["time"]], 
          " | ", 
          entry[["package"]], 
          " | ", 
          entry[["function_or_method"]],
          " | ",
          entry[["process"]],
          " | ",
          entry[["status"]]
        ),
        "\n"
      )
      
    }
    
  }
)



.clean_geodata <- function(
    sf_geodata,
    unit_col,
    reduce_cols = TRUE
) {
  
  if(isTRUE(reduce_cols)) {
    sf_geodata <-
      sf_geodata[c(unit_col, "geometry")]  
  }
  
  nrow_geodata1 <- nrow(sf_geodata)
  
  sf_geodata <-
    sf_geodata[!is.na(sf_geodata[[unit_col]]),]
  
  nrow_geodata2 <- nrow(sf_geodata)
  
  if(nrow_geodata2 < nrow_geodata1) {
    warning(paste0("sf object includes ", nrow_geodata1-nrow_geodata2, " objects with NA values of '", unit_col, "', which are skipped"))
  }
  
  sf_geodata <-
    sf_geodata[!is.na(st_geometry(sf_geodata)) & st_is_valid(sf_geodata),]
  
  nrow_geodata3 <- nrow(sf_geodata)
  
  if(nrow_geodata3 < nrow_geodata2) {
    warning(paste0("sf object includes ", nrow_geodata2-nrow_geodata3, " objects with missing or invalid geometry, which are skipped"))
  }
  
  geodata_crs <- st_crs(sf_geodata)$input
  
  invisible(
    list(
      geodata = sf_geodata,
      N = c(
        nrow_geodata1, 
        nrow_geodata2,
        nrow_geodata3
      ),
      geodata_crs = geodata_crs
    )
  )
  
}


.rename_cols <-
  function(
    df,
    colnames_mapping
  ) {
    
    i <- 0
    
    for(i in 1:length(colnames_mapping)) {
      colnames(df)[names(df) == names(colnames_mapping)[i]] <- colnames_mapping[[i]]
    }
    
    return(df)
    
  }


.df_date_col <-
  function(
    df,
    col_date,
    time_format,
    verbose = FALSE
  ) {
    
    if(!"data.frame" %in% class(df)) {
      stop("User-specified data is not of class data.frame")
    }
    
    if(!col_date %in% colnames(df)) {
      stop(paste0("User-specified time column '", col_date, "' is not in data frame"))
    }
    
    if(isTRUE(verbose)) {
      cat(paste0("Checking time column '", col_date, "' from data frame with ", nrow(df), " rows ... "))
    }
    
    error_message <- NULL
    
    tryCatch(
      {
        
        df[[col_date]] <- 
          as.Date(
            df[[col_date]], 
            format = time_format
          )
        
      },
      
      error = function(cond) {
        
        error_message <- paste0("Time column '", col_date, "' does not seem to have a valid date format: ", conditionMessage(cond))
        stop(error_message)
        
      }
      
    )
    
    date_col_NA <- 
      nrow(df[is.na(df[[col_date]]),])
    
    if(date_col_NA > 0) {
      warning(paste0("There are ", date_col_NA, " NA values in time column '", col_date, "'"))
    }
    
    invisible(df)
   
  }