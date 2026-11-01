#---------------------------------------------------------------
# Name:        helper (swash package)
# Purpose:     Helper functions for the swash package
# Author:      Thomas Wieland 
#              ORCID: 0000-0001-5168-9846
#              mail: geowieland@googlemail.com
# Version:     1.2.0
# Last update: 2026-10-25 13:27
# Copyright (c) 2025-2026 Thomas Wieland
#---------------------------------------------------------------


is_balanced <- function(
  data,
  col_cases,
  col_date,
  col_region,
  as_balanced = TRUE,
  fill_missing = NA,
  drop_duplicates = TRUE,
  verbose = FALSE
) {

  if (!is.data.frame(data)) {
    stop("Specified object is not a data.frame", call. = FALSE)
  }

  if (isTRUE(verbose)) {
    cat("Checking for columns ... ")
  }

  required_cols <- c(col_cases, col_date, col_region)
  missing_cols <- setdiff(required_cols, colnames(data))

  if (length(missing_cols) > 0) {
    stop(
      paste(
        "The following columns are missing:",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (isTRUE(verbose)) {
    cat("OK\n")
  }

 
  if (isTRUE(verbose)) {
    cat("Checking for missing region/date keys ... ")
  }

  if (anyNA(data[[col_region]]) || anyNA(data[[col_date]])) {
    stop("Missing values are not allowed in the region or date columns.")
  }

  if (isTRUE(verbose)) {
    cat("OK\n")
  }

 
  if (isTRUE(verbose)) {
    cat("Checking for duplicates ... ")
  }

  nrow_original <- nrow(data)

  duplicate_rows <- duplicated(
    data[c(col_region, col_date)]
  )

  duplicates_count <- sum(duplicate_rows)
  duplicates <- duplicates_count > 0

  if (isTRUE(verbose)) {
    cat("OK\n")
  }

  if (duplicates) {

    if (isTRUE(drop_duplicates)) {

      data <- data[!duplicate_rows, , drop = FALSE]

      warning(
        paste0(
          "User-stated data frame included ",
          duplicates_count,
          " rows with duplicated region/time combinations, ",
          "which were excluded."
        ),
        call. = FALSE
      )

    } else {

      stop(
        paste0(
          "User-stated data frame includes ",
          duplicates_count,
          " rows with duplicated region/time combinations. ",
          "Use drop_duplicates = TRUE to clean up the dataset."
        ),
        call. = FALSE
      )

    }
  }


  check_balance <- function(x) {

    if (nrow(x) == 0) {
      return(FALSE)
    }

    if (
      anyNA(x[[col_region]]) ||
      anyNA(x[[col_date]])
    ) {
      return(FALSE)
    }

    regions <- unique(x[[col_region]])
    dates <- unique(x[[col_date]])

    panel_table <- table(
      factor(x[[col_region]], levels = regions),
      factor(x[[col_date]], levels = dates)
    )

    structurally_balanced <- all(panel_table == 1L)

    cases_complete <- !anyNA(x[[col_cases]])

    isTRUE(structurally_balanced && cases_complete)
  }

  
  if (isTRUE(verbose)) {
    cat("Checking whether data is balanced ... ")
  }

  data_balanced <- check_balance(data)

  if (isTRUE(verbose)) {
    cat("OK", "\n")
  }


  if (!data_balanced && isTRUE(as_balanced)) {

    balance_fun <- get(
      "as_balanced",
      mode = "function",
      inherits = TRUE
    )

    data <- balance_fun(
      data,
      col_cases,
      col_date,
      col_region,
      fill_missing = fill_missing,
      verbose = verbose
    )

    data_balanced <- check_balance(data)
    
  }

  results <- list(
    data_balanced = data_balanced,
    data = data
  )

  return(results)
}


as_balanced <-
  function(
    data,
    col_cases, 
    col_date, 
    col_region,
    fill_missing = NA,
    drop_duplicates = TRUE,
    verbose = FALSE
  ) {

    if (!is.data.frame(data)) {
      stop("Specified object is not a data.frame", call. = FALSE)
    }

    if (isTRUE(verbose)) {
      cat("Checking for columns ... ")
    }

    required_cols <- c(col_cases, col_date, col_region)
    missing_cols <- setdiff(required_cols, colnames(data))

    if (length(missing_cols) > 0) {
      stop(
        paste(
          "The following columns are missing:",
          paste(missing_cols, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    if (isTRUE(verbose)) {
      cat("OK", "\n")
    }

    if (anyNA(data[[col_region]]) || anyNA(data[[col_date]])) {
      stop("Missing values are not allowed in the region or date columns.")
    }

    if (isTRUE(verbose)) {
      cat("Checking for duplicates ... ")
    }

    nrow_data <- nrow(data)

    duplicate_rows <- duplicated(data[c(col_region, col_date)])
    duplicates_count <- sum(duplicate_rows)

    if (isTRUE(verbose)) {
      cat("OK", "\n")
    }

    if (duplicates_count > 0) {

      if (isTRUE(drop_duplicates)) {

        data <- data[!duplicate_rows, , drop = FALSE]

        warning(
          paste0("User-stated data frame included ", duplicates_count, " rows with duplicated region/time combinations, which were excluded."),
          call. = FALSE
        )

      } else {

        stop(
          paste0(
            "User-stated data frame includes ", duplicates_count, " rows with duplicated region/time combinations. Use drop_duplicates = TRUE to clean up the dataset.")
        )

      }
    }

    if (nrow(data) == 0) {
      stop("Cannot balance an empty data frame.")
    }

    if (isTRUE(verbose)) {
      cat(paste0("Balancing the dataset and filling missing values with ", paste(fill_missing, collapse = ", "), " ... "))
    }

    N_names <- unique(data[[col_region]])
    TP_t <- unique(data[[col_date]])

    N_x_TPt <- expand.grid(
      N_names,
      TP_t,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )

    colnames(N_x_TPt) <- c(col_region, col_date)

    data <- merge(
      data,
      N_x_TPt,
      by = c(col_region, col_date),
      all.y = TRUE,
      sort = FALSE
    )

    missing_cases <- is.na(data[[col_cases]])

    if (any(missing_cases)) {
      data[[col_cases]][missing_cases] <- fill_missing
    }

    if (isTRUE(verbose)) {
      cat("OK", "\n")
    }

    return(data)

  }


.is_cumulative <- function(y) {
  
  if(all(diff(y) >= 0)) {
    return (TRUE)
  } else {
    return (FALSE)
  }
  
}

.check_negative <- function(y) {
 
  if(any(y < 0)) {
    warning(paste0("Data contains ", length(y[y < 0]), " negative values which are set to NA"))
  }
   
  y[y < 0] <- NA
  
  invisible(y)
  
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
      package = paste0(package_name, " ", swash_package_version),
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


.clean_timeseries <- function(
    y,
    t,
    impute_missing_values = FALSE,
    accept_non_positive_values = TRUE,
    add_constant = 1,
    time_counter = TRUE,
    verbose = FALSE
) {


  if (length(y) != length(t)) {
    stop(paste0("Vectors 'y' and 't' differ in length: ", length(y), ", ", length(t)))
  }

  if (isTRUE(verbose)) {
    cat("Checking for correct types ... ")
  }

  type_errors <- character(0)

  if (!is.numeric(y)) {
    type_errors <- c(
      type_errors,
      "y vector must be of class 'numeric'."
    )
  }

  if (!is.numeric(t) && !inherits(t, "Date")) {
    type_errors <- c(
      type_errors,
      "Time vector must be of class 'numeric' or 'Date'."
    )
  }

  if (isTRUE(verbose)) {
    cat("OK", "\n")
  }

  if (length(type_errors) > 0) {
    stop(
      paste(
        "Type errors in input time series data:",
        paste(type_errors, collapse = " ")
      )
    )
  }

  if (length(y) == 0) {
    stop("Input time series is empty.", call. = FALSE)
  }


  if (isTRUE(verbose)) {
    cat("Checking time series for incomplete cases ... ")
  }

  yt_df <- data.frame(y = y, t = t)
  incomplete_cases <- sum(!complete.cases(yt_df))

  if (isTRUE(verbose)) {
    cat("OK", "\n")
  }

  if (incomplete_cases > 0) {

    message(
      paste0(
        "There are ", incomplete_cases,
        " incomplete rows in the specified time series."
      )
    )

    if (isTRUE(impute_missing_values)) {

      if (isTRUE(verbose)) {
        cat("Imputing missing y values ... ")
      }

      missing_times <- is.na(yt_df$t)
      incomplete_cases_t <- sum(missing_times)

      if (incomplete_cases_t > 0) {
        yt_df <- yt_df[!missing_times, , drop = FALSE]

        warning(
          paste0(
            incomplete_cases_t,
            " cases with missing time points were removed."
          ),
          call. = FALSE
        )
      }

      yt_df <- yt_df[
        order(yt_df$t),
        ,
        drop = FALSE
      ]

      if (anyDuplicated(yt_df$t)) {
        stop(
            "Duplicate time points were found. Resolve duplicate times before imputing missing values."
        )
      }

      if (sum(!is.na(yt_df$y)) < 2) {
        stop(
          paste(
            "At least two non-missing y values are required",
            "for interpolation."
          ),
          call. = FALSE
        )
      }

      yt_df$y <- na.approx(
        yt_df$y,
        x = as.numeric(yt_df$t),
        na.rm = FALSE
      )

      remaining_missing <- !complete.cases(yt_df)

      if (any(remaining_missing)) {

        n_removed <- sum(remaining_missing)

        yt_df <- yt_df[
          !remaining_missing,
          ,
          drop = FALSE
        ]

        warning(
          paste0(
            n_removed,
            " rows with non-interpolable values were removed."
          ),
          call. = FALSE
        )
      }

      y <- yt_df$y
      t <- yt_df$t

      if (isTRUE(verbose)) {
        cat("OK", "\n")
      }

    } else {

      if (isTRUE(verbose)) {
        cat("Removing incomplete rows ... ")
      }

      yt_df <- yt_df[
        complete.cases(yt_df),
        ,
        drop = FALSE
      ]

      y <- yt_df$y
      t <- yt_df$t

      if (isTRUE(verbose)) {
        cat("OK", "\n")
      }

      warning(
        paste0(
          incomplete_cases,
          " incomplete cases were removed."
        ),
        call. = FALSE
      )
    }
  }

  if (length(y) == 0) {
    stop(
      "No complete observations remain after cleaning.",
      call. = FALSE
    )
  }

  
  if (isTRUE(verbose)) {
    cat("Checking time points ... ")
  }

  time_order <- order(t)
  y <- y[time_order]
  t <- t[time_order]

  if (anyDuplicated(t)) {
    stop(
      paste(
        "Duplicate time points were found.",
        "Resolve duplicate times before further analysis."
      ),
      call. = FALSE
    )
  }

  if (isTRUE(verbose)) {
    cat("OK", "\n")
  }


  if (inherits(t, "Date") && isTRUE(time_counter)) {

    if (isTRUE(verbose)) {
      cat(
        "Time vector is of class 'Date'. ",
        "Calculating time counter ... "
      )
    }

    start_date <- min(t)

    t <- as.integer(t - start_date)

    if (isTRUE(verbose)) {
      cat("OK", "\n")
    }
  }

  
  if (!isTRUE(accept_non_positive_values) && any(y <= 0)) {

    if (isTRUE(verbose)) {
      cat("Correcting non-positive values ... ")
    }

    if (
      !is.null(add_constant) &&
      is.numeric(add_constant) &&
      length(add_constant) == 1 &&
      is.finite(add_constant) &&
      add_constant >= 0
    ) {

      y_new <- y + add_constant

      if (any(y_new <= 0)) {
        stop(
          paste(
            "add_constant is insufficient to make all",
            "y values strictly positive."
          ),
          call. = FALSE
        )
      }

      y <- y_new

      message(
        paste0(
          "Input time series data contains values <= 0. ",
          "All values are increased by a constant equal to ",
          add_constant, "."
        )
      )

    } else {

      keep <- y > 0

      y <- y[keep]
      t <- t[keep]

      warning(
        paste(
          "Input time series data contains values <= 0.",
          "These values and their corresponding time points",
          "are dropped."
        ),
        call. = FALSE
      )
    }

    if (isTRUE(verbose)) {
      cat("OK", "\n")
    }
  }

  
  if (length(y) != length(t)) {
    stop(
      "Internal error: cleaned vectors differ in length.",
      call. = FALSE
    )
  }

  if (length(y) == 0) {
    stop(
      "No observations remain after cleaning.",
      call. = FALSE
    )
  }

  if (anyNA(y) || anyNA(t)) {
    stop(
      "Missing values remain after cleaning.",
      call. = FALSE
    )
  }

  yt_vectors <- list(
    y = y,
    t = t
  )

  invisible(yt_vectors)
}



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
        
        if(isTRUE(verbose)) {
          cat("OK", "\n")
        }
        
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