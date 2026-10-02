#---------------------------------------------------------------
# Name:        swash (swash package)
# Purpose:     Swash-Backwash Model for the Single Epidemic Wave
#              Functions, classes and corresponding methods
# Author:      Thomas Wieland 
#              ORCID: 0000-0001-5168-9846
#              mail: geowieland@googlemail.com
# Version:     2.0.0
# Last update: 2026-09-08 19:13
# Copyright (c) 2022-2026 Thomas Wieland
#---------------------------------------------------------------



# Class sbm:
setClass(
  "sbm",
  slots = list(
    R_0A = "numeric",
    integrals = "numeric",
    velocity = "numeric",
    occ_regions = "data.frame",
    SIR_regions = "data.frame",
    cases_by_date = "data.frame",
    cases_by_region = "data.frame",
    input_data = "data.frame",
    data_statistics = "numeric",
    col_names = "character",
    timestamp = "list"
  ),
  prototype = list(
    timestamp = list()
  )  
)

# Methods for class sbm:
setMethod(
  "summary",
  "sbm",
  function(object) {
    
    cat(model_descriptions[["sbm"]], "\n\n")
    
    cat("Integrals\n")
    cat(sprintf("  Susceptible areas : %.3f\n", object@integrals[1]))
    cat(sprintf("  Infected areas    : %.3f\n", object@integrals[2]))
    cat(sprintf("  Recovered areas   : %.3f\n\n", object@integrals[3]))
    
    cat("Velocity\n")
    cat(sprintf("  Leading edge      : %.3f\n", object@velocity[1]))
    cat(sprintf("  Following edge    : %.3f\n\n", object@velocity[2]))
    
    cat(sprintf("Spatial reproduction number : %.3f\n\n", object@R_0A))
    
    cat("Input data\n")
    cat(sprintf("  Units       : %s\n", object@data_statistics[1]))
    cat(sprintf("  No-case     : %s\n", object@data_statistics[3]))
    cat(sprintf("  Time points : %s\n", object@data_statistics[2]))
    cat(sprintf(
      "  Balanced    : %s\n",
      ifelse(object@data_statistics[4], "YES", "NO")
    ))
  }
)

setMethod(
  "print", 
  "sbm", 
  function(x) {
    
    cat(paste0(model_descriptions[["sbm"]], " with ", x@data_statistics[1], " spatial units and ", x@data_statistics[2], " time points"), "\n")
    cat ("Use summary() for results", "\n")
    
    invisible(x)
    
  }
)

setMethod(
  "show", 
  "sbm", 
  function(object) {
    
    cat(paste0(model_descriptions[["sbm"]], " with ", object@data_statistics[1], " spatial units and ", object@data_statistics[2], " time points"), "\n")
    cat ("Use summary() for results", "\n")
    
    invisible((object))
    
  }
)

setMethod(
  "plot", 
  "sbm", 
  function(
    x, 
    y = NULL, 
    col_edges = "blue",
    xlab_edges = "Time", 
    ylab_edges = "Regions",
    main_edges = "Edges",
    col_SIR = c("blue", "red", "green"),
    lty_SIR = c("solid", "solid", "solid"),
    lwd_SIR = c(1,1,1),
    xlab_SIR = "Time", 
    ylab_SIR = "Regions",
    main_SIR = "SIR integrals",
    col_cases = "red",
    lty_cases = "solid",
    lwd_cases = 1,
    xlab_cases = "Time", 
    ylab_cases = "Infections",
    main_cases = "Daily infections",
    xlab_cum = "Cases",
    ylab_cum = "Regions",
    main_cum = "Cumulative infections per region",
    horiz_cum = TRUE,
    separate_plots = FALSE
  ) {
    
    par_old <- par(no.readonly = TRUE)
    on.exit(par(par_old))
    dev.new()
    
    if (separate_plots == FALSE) {
      par(mfrow = c(2,2))
    }
    
    barplot(
      x@occ_regions$LE_FE,
      col = col_edges,
      xlab = xlab_edges, 
      ylab = ylab_edges,
      main = main_edges
    )
    
    plot (
      x = as.Date(x@SIR_regions$date), 
      y = x@SIR_regions$susceptible, 
      col = col_SIR[1],
      xlab = xlab_SIR,
      ylab = ylab_SIR,
      "l",
      lty = lty_SIR[1],
      lwd = lwd_SIR[1],
      ylim = c(0, x@data_statistics[1]+1),
      main = main_SIR
    )
    lines (
      x = as.Date(x@SIR_regions$date), 
      y = x@SIR_regions$infected, 
      col = col_SIR[2],
      lty = lty_SIR,
      lwd = lwd_SIR[2]
    )
    lines (
      x = as.Date(x@SIR_regions$date), 
      y = x@SIR_regions$recovered, 
      col = col_SIR[3],
      lty = lty_SIR[3],
      lwd = lwd_SIR[3]
    )
    
    plot(
      x@cases_by_date$date, 
      x@cases_by_date$cases,
      type = "l",
      lty = lty_cases,
      lwd = lwd_cases,
      col = col_cases,
      xlab = xlab_cases, 
      ylab = ylab_cases,
      main = main_cases)
    
    x@cases_by_region <- x@cases_by_region[order(x@cases_by_region$cases_cumulative), ]
    barplot(
      height = x@cases_by_region$cases_cumulative,
      horiz = horiz_cum,
      names.arg = x@cases_by_region$region,
      xlab = xlab_cum,
      ylab = ylab_cum,
      main = main_cum,
      las = 1)
    
    invisible(x)
  }
)

setMethod(
  "confint", 
  "sbm", 
  function(
    object,
    iterations = 100,
    samples_ratio = 0.8,
    alpha = 0.05,
    replace = TRUE,
    verbose = FALSE
  ) {
    
    N = object@data_statistics[1]
    TP = object@data_statistics[2]
    input_data = object@input_data
    obs <- nrow(object@input_data)
    
    regions <- as.character(levels(as.factor(object@input_data[[object@col_names[3]]])))
    regions_to_sample <- round(N*samples_ratio)
    
    bootstrap_config <- list(
      iterations = iterations,
      samples_ratio = samples_ratio,
      regions_to_sample = regions_to_sample,
      alpha = alpha,
      replace = replace)
    
    if(isTRUE(verbose)) {
      cat(paste0("Calculating bootstrap confidence intervals with ", iterations, " iterations and alpha = ", alpha))
    }
    
    i <- 0
    swash_bootstrap <- matrix(ncol = 9, nrow = iterations)
    
    for (i in 1:iterations) {
      
      regions_sample <- 
        sample(
          x = regions, 
          size = regions_to_sample,
          replace = replace)
      
      data_sample <- 
        input_data[as.character(input_data[[object@col_names[3]]]) %in% regions_sample,]
      
      data_sample_size <- nrow(data_sample)
      
      data_sample_swash <- 
        swash_backwash(
          data = data_sample,
          col_cases = object@col_names[1], 
          col_date = object@col_names[2], 
          col_region = object@col_names[3]
        )
      
      swash_bootstrap[i, 1] <- i
      swash_bootstrap[i, 2] <- data_sample_swash@integrals[1]
      swash_bootstrap[i, 3] <- data_sample_swash@integrals[2]
      swash_bootstrap[i, 4] <- data_sample_swash@integrals[3]
      swash_bootstrap[i, 5] <- data_sample_swash@velocity[1]
      swash_bootstrap[i, 6] <- data_sample_swash@velocity[2]
      swash_bootstrap[i, 7] <- data_sample_swash@velocity[3]
      swash_bootstrap[i, 8] <- data_sample_swash@R_0A
      swash_bootstrap[i, 9] <- data_sample_size
      
    }
    
    colnames(swash_bootstrap) <-
      c(
        "iteration", 
        "S_A", 
        "I_A", 
        "R_A", 
        "t_LE", 
        "t_FE", 
        "t_FE-t_LE", 
        "R_0A", 
        "sample_size"
      )
    
    ci_lower <- alpha/2
    ci_upper <- 1-(alpha/2)
    
    S_A_ci <- quantile(swash_bootstrap[,2], probs = c(ci_lower, ci_upper))
    I_A_ci <- quantile(swash_bootstrap[,3], probs = c(ci_lower, ci_upper))
    R_A_ci <- quantile(swash_bootstrap[,4], probs = c(ci_lower, ci_upper))
    integrals_ci <- list(S_A_ci = S_A_ci, I_A_ci = I_A_ci, R_A_ci= R_A_ci)
    
    t_LE_ci <- quantile(swash_bootstrap[,5], probs = c(ci_lower, ci_upper))
    t_FE_ci <- quantile(swash_bootstrap[,6], probs = c(ci_lower, ci_upper))
    t_FE_t_LE_ci <- quantile(swash_bootstrap[,7], probs = c(ci_lower, ci_upper))
    velocity_ci <- list(t_LE_ci = t_LE_ci, t_FE_ci = t_FE_ci, t_FE_t_LE_ci= t_FE_t_LE_ci)
    
    R_0A_ci <- quantile(swash_bootstrap[,8], probs = c(ci_lower, ci_upper))
    
    if(isTRUE(verbose)) {
      print("OK", "\n")
    }
    
    sbm_ci_object <- new(
      "sbm_ci", 
      R_0A = object@R_0A, 
      integrals = object@integrals, 
      velocity = object@velocity,
      occ_regions = object@occ_regions, 
      cases_by_date = object@cases_by_date,
      input_data = object@input_data,
      data_statistics = object@data_statistics,
      col_names = object@col_names,
      integrals_ci = integrals_ci,
      velocity_ci = velocity_ci,
      R_0A_ci = R_0A_ci,
      iterations = data.frame(swash_bootstrap),
      ci = c(ci_lower, ci_upper),
      config = bootstrap_config
    )
    
    invisible(sbm_ci_object)
    
  }
)


# Class sbm_ci:
setClass(
  "sbm_ci",
  slots = list(
    R_0A = "numeric",
    integrals = "numeric",
    velocity = "numeric",
    occ_regions = "data.frame",
    cases_by_date = "data.frame",
    cases_by_region = "data.frame",
    input_data = "data.frame",
    data_statistics = "numeric",
    col_names = "character",
    integrals_ci = "list", 
    velocity_ci = "list", 
    R_0A_ci = "numeric", 
    iterations = "data.frame",
    ci = "numeric",
    config = "list"
  )
)

# Methods of class sbm_ci:
setMethod(
  "summary", 
  "sbm_ci", 
  function(object) {
    
    cat("Confidence Intervals for Swash-Backwash Model\n\n")
    
    cat("Integrals\n")
    cat(
      sprintf(
        "  Susceptible areas : [%.3f, %.3f]\n",
        object@integrals_ci$S_A_ci[1],
        object@integrals_ci$S_A_ci[2])
    )
    cat(
      sprintf(
        "  Infected areas    : [%.3f, %.3f]\n",
        object@integrals_ci$I_A_ci[1], 
        object@integrals_ci$I_A_ci[2]
      )
    )
    cat(
      sprintf(
        "  Recovered areas   : [%.3f, %.3f]\n\n",
        object@integrals_ci$R_A_ci[1], 
        object@integrals_ci$R_A_ci[2]
      )
    )
    
    cat("Velocity\n")
    cat(
      sprintf(
        "  Leading edge      : [%.3f, %.3f]\n", 
        object@velocity_ci$t_LE_ci[1], 
        object@velocity_ci$t_LE_ci[2]
      )
    )
    cat(
      sprintf(
        "  Following edge    : [%.3f, %.3f]\n\n", 
        object@velocity_ci$t_FE_ci[1], 
        object@velocity_ci$t_FE_ci[2]
      )
    )
    
    cat(
      sprintf(
        "Spatial reproduction number : [%.3f, %.3f]\n\n",
        object@R_0A_ci[1], 
        object@R_0A_ci[2]
      )
    )
    
    cat("Configuration for confidence intervals\n")
    cat(
      sprintf(
        "  CI alpha    : %s\n", 
        object@config$alpha
      )
    )
    cat(sprintf("  Sample      : %s %% (%s units)\n", 
                object@config$samples_ratio*100, object@config$regions_to_sample))
    cat(sprintf("  Iterations  : %s\n", object@config$iterations))
    cat(sprintf("  Bootstrap   : %s\n", ifelse(object@config$replace, "YES", "NO")))
    
    cat("\nInput data\n")
    cat(sprintf("  Units       : %s\n", object@data_statistics[1]))
    cat(sprintf("  No-case     : %s\n", object@data_statistics[3]))
    cat(sprintf("  Time points : %s\n", object@data_statistics[2]))
    cat(sprintf("  Balanced    : %s\n", ifelse(object@data_statistics[4], "YES", "NO")))
  }
)

setMethod(
  "print", 
  "sbm_ci", 
  function(x) {
    
    cat(paste0("Confidence Intervals for Swash-Backwash Model with ", x@data_statistics[1], " spatial units and ", x@data_statistics[2], " time points"), "\n")
    cat(paste0("Resampling of ", x@config$regions_to_sample, " spatial units (", x@config$samples_ratio*100, " %) with ", x@config$iterations, " iterations"), "\n")
    cat ("Use summary() for results", "\n")
    
  })


setMethod(
  "show", 
  "sbm_ci", 
  function(object) {
    
    cat(paste0("Confidence Intervals for Swash-Backwash Model with ", object@data_statistics[1], " spatial units and ", object@data_statistics[2], " time points"), "\n")
    cat(paste0("Resampling of ", object@config$regions_to_sample, " spatial units (", object@config$samples_ratio*100, " %) with ", object@config$iterations, " iterations"), "\n")
    cat ("Use summary() for results", "\n")
    
  })


setMethod(
  "plot", 
  "sbm_ci", 
  function(
    x, 
    y = NULL, 
    col_bars = "grey",
    col_ci = "red"
  ) {
    
    par_old <- par(no.readonly = TRUE)
    on.exit(par(par_old)) 
    dev.new()
    
    par (mfrow = c(2,3))
    
    alpha <- x@config$alpha
    
    hist_ci (
      x@iterations$S_A,
      alpha = alpha,
      col_bars = col_bars, 
      col_ci = col_ci, 
      main = "Susceptible areas integral",
      xlab = "S_A",
      ylab = "Frequency"
    )
    hist_ci (
      x@iterations$I_A,
      alpha = alpha,
      col_bars = col_bars, 
      col_ci = col_ci, 
      main = "Infected areas integral",
      xlab = "I_A",
      ylab = "Frequency"
    )
    hist_ci (
      x@iterations$R_A, 
      alpha = alpha,
      col_bars = col_bars, 
      col_ci = col_ci, 
      main = "Recovered areas integral",
      xlab = "R_A",
      ylab = "Frequency"
    )
    
    hist_ci (
      x@iterations$t_LE,
      alpha = alpha,
      col_bars = col_bars, 
      col_ci = col_ci, 
      main = "Leading edge",
      xlab = "t_LE",
      ylab = "Frequency"
    )
    hist_ci (
      x@iterations$t_FE, 
      alpha = alpha,
      col_bars = col_bars, 
      col_ci = col_ci, 
      main = "Following edge",
      xlab = "t_FE",
      ylab = "Frequency"
    )
    
    hist_ci (
      x@iterations$R_0A, 
      alpha = alpha,
      col_bars = col_bars, 
      col_ci = col_ci, 
      main = "Spatial reproduction number",
      xlab = "R_0A",
      ylab = "Frequency"
    )
    
    invisible(x)
  }
)


# Class countries:
setClass(
  "countries",
  slots = list(
    sbm_ci1 = "sbm_ci",
    sbm_ci2 = "sbm_ci",
    D = "numeric",
    D_ci = "numeric",
    config = "list",
    country_names = "character",
    indicator = "character"
  )
)

# Methods for class countries:
setMethod(
  "plot", 
  "countries", 
  function(
    x, 
    y = NULL, 
    col_bars = "grey",
    col_ci = "red"
  ) {
    
    par_old <- par(no.readonly = TRUE)
    on.exit(par(par_old))
    dev.new()
    
    par (mfrow = c(2,2))
    
    alpha <- x@config$alpha
    
    indicator <- x@indicator
    
    hist_ci (
      x@sbm_ci1@iterations[[indicator]],
      alpha = alpha,
      col_bars = col_bars, 
      col_ci = col_ci, 
      main = paste0("Indicator ", indicator, " for ", x@country_names[1]),
      xlab = indicator,
      ylab = "Frequency"
    )
    
    hist_ci (
      x@sbm_ci2@iterations[[indicator]],
      alpha = alpha,
      col_bars = col_bars, 
      col_ci = col_ci, 
      main = paste0("Indicator ", indicator, " for ", x@country_names[2]),
      xlab = indicator,
      ylab = "Frequency"
    )
    
    country1 <- data.frame(
      x@sbm_ci1@iterations[[indicator]], 
      x@country_names[1]
    )
    colnames(country1) <-
      c(indicator, "country")
    country2 <- data.frame(
      x@sbm_ci2@iterations[[indicator]], 
      x@country_names[2]
    )
    colnames(country2) <-
      c(indicator, "country")
    
    iterations_countries <-
      rbind(
        country1,
        country2
      )
    
    boxplot(
      iterations_countries[[indicator]] ~ iterations_countries$country,
      xlab = "Country",
      ylab = indicator
    )
    
    hist_ci (
      x@D,
      alpha = alpha,
      col_bars = col_bars, 
      col_ci = col_ci, 
      main = paste0("Difference in indicator ", indicator),
      xlab = paste0("D (", indicator, " )"),
      ylab = "Frequency"
    )
    
    invisible(x)
  }
)

setMethod(
  "summary", 
  "countries", 
  function(object) {
    
    D_ci <- object@D_ci
    indicator <- object@indicator
    
    D_mean <- mean(object@D, na.rm = TRUE)
    D_median <- median(object@D, na.rm = TRUE)
    
    cat("Two-country comparison for Swash-Backwash Model\n\n")
    
    cat(
      sprintf(
        "Difference in %s\n", 
        indicator
      )
    )
    cat(
      sprintf(
        "  Mean   : %.3f\n", 
        D_mean
      )
    )
    cat(
      sprintf(
        "  Median : %.3f\n", 
        D_median
      )
    )
    cat(
      sprintf(
        "  Confidence interval : [%.3f, %.3f]\n\n", 
        D_ci[1], 
        D_ci[2]
      )
    )
    
    cat("Configuration for confidence intervals\n")
    cat(
      sprintf(
        "  CI alpha    : %s\n", 
        object@config$alpha
      )
    )
    cat(
      sprintf(
        "  Sample      : %s %%\n", 
        object@config$samples_ratio*100
      )
    )
    cat(
      sprintf(
        "  Iterations  : %s\n", 
        object@config$iterations
      )
    )
    cat(
      sprintf(
        "  Bootstrap   : %s\n", 
        ifelse(object@config$replace, "YES", "NO")
      )
    )
    
  }
)

setMethod(
  "show", 
  "countries", 
  function(object) {
    
    cat(paste0("Two-country comparison with Swash-Backwash Model"), "\n")
    cat ("Use summary() for results", "\n")
    
  })

setMethod(
  "print", 
  "countries", 
  function(x) {
    
    cat(paste0("Two-country comparison with Swash-Backwash Model"), "\n")
    cat ("Use summary() for results", "\n")
    
  })


# Function for Swash-Backwash Model:
swash_backwash <- 
  function(
    infpan = NULL,
    data = NULL,
    col_cases = NULL, 
    col_date = NULL, 
    col_region = NULL,
    time_format = "%Y-%m-%d",
    verbose = FALSE
  ) {
    
    if (is.null(infpan) && is.null(data)) {
      stop("Either an infpan object or a data.frame must be stated") 
    }
    if (!is.null(data) && (is.null(col_cases) || is.null(col_date) || is.null(col_region))) {
      stop("When specifying a data.frame, three parameters must be stated: 'col_cases', 'col_date', 'col_region'")
    }
    if (!is.null(infpan) && !is.null(data)) {
      message("Only a data.frame OR an infpan object must be stated. Parameter 'data' is ignored")
      data <- NULL
    }
    
    if (!is.null(data)) {
      
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
      
      if(isTRUE(verbose)) {
        cat("Reading and diagnosing infections panel data ... ")
      }
      
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
      
      if (isTRUE(verbose)) {
        cat("OK", "\n")
      }
      
    }
    
    if (!is.null(infpan)) {
      
      if (isTRUE(verbose)) {
        cat("Loading infections panel data object ... ")
      }
      
      data <- infpan@input_data
      
      col_region <- infpan@index_col_names[1]
      col_date <- infpan@index_col_names[2]
      col_cases <- infpan@cases_col_name
      
      time_format <- infpan@time_format
      
      N <- infpan@data_statistics[1]
      N_names <- levels(as.factor(data[[col_region]]))
      N_withoutcases <- infpan@data_statistics[3]
      TP <- infpan@data_statistics[2]
      TP_t <- levels(as.factor(data[[col_date]]))
      data_balanced <- infpan@data_statistics[4]
      
      if (isTRUE(verbose)) {
        cat("OK", "\n")
      }
      
    }
    
    data[[col_date]] <- 
      as.Date(
        data[[col_date]], 
        format = time_format
      )
    
    data <- data[order(data[[col_region]], data[[col_date]]),]
    
    if (isTRUE(verbose)) {
      
      cat(paste0("Infections panel data include ", N, " regions in column '", col_region, "' and ", TP, " time points in column '", col_date, "', with cases in column '", col_cases, "'. "))
      
      if (isTRUE(data_balanced)) {
        cat("The data is balanced.")
      } else {
        cat("The data is not balanced.")
      }
      cat("\n")
      
    }
    
    if(isTRUE(verbose)) {
      cat(paste0("Calculating Swash-Backwash Model for ", N, " regions and ", TP, " time points ... "))
    }
    
    first_occ_regions <- data.frame(matrix(ncol = N+1, nrow = TP))
    first_occ_regions[,1] <- TP_t
    colnames(first_occ_regions)[1] <- "date"
    colnames(first_occ_regions)[2:(N+1)] <- N_names
    
    last_occ_regions <- data.frame(matrix(ncol = N+1, nrow = TP))
    last_occ_regions[,1] <- TP_t
    colnames(last_occ_regions)[1] <- "date"
    colnames(last_occ_regions)[2:(N+1)] <- N_names
    
    i <- 0
    
    for (i in 1:N) {
      
      data_n <- data[data[[col_region]] == N_names[i],]
      data_n$occurence <- 0
      
      if (sum(data_n[[col_cases]], na.rm = TRUE) > 0) {
        data_n[data_n[[col_cases]] > 0,]$occurence <- 1
        
      } else {
        
        N_withoutcases <- N_withoutcases+1
      }
      
      first_occ <- which(data_n$occurence == 1)[1]
      data_n$LE <- 0
      data_n$LE[first_occ] <- 1
      
      first_occ_regions[,i+1] <- data_n$LE
      
      colnames(first_occ_regions)[i+1] <- paste0("Region_", N_names[i])
      
      last_occ <- which(data_n$occurence == 1)[length(which(data_n$occurence == 1))]
      data_n$FE <- 0
      data_n$FE[last_occ] <- 1
      
      last_occ_regions[,i+1] <- data_n$FE
      
      colnames(last_occ_regions)[i+1] <- paste0("Region_", N_names[i])
      
    }
    
    first_occ_regions$no_regions_LE <- rowSums(first_occ_regions[, 2:(N+1)])
    first_occ_regions$t <- seq (1:TP)
    first_occ_regions$t_x_nt <- first_occ_regions$t*first_occ_regions$no_regions_LE
    
    t_LE <- sum(first_occ_regions$t_x_nt)/N
    
    last_occ_regions$no_regions_FE <- rowSums(last_occ_regions[, 2:(N+1)])
    last_occ_regions$t <- seq (1:TP)
    last_occ_regions$t_x_nt <- last_occ_regions$t*last_occ_regions$no_regions_FE
    
    t_FE <- sum(last_occ_regions$t_x_nt)/N
    
    S_A <- (t_LE-1)/TP
    I_A <- (t_FE/TP)-S_A
    R_A <- 1-(S_A+I_A)
    R_0A <- (1-S_A)/(1-R_A)
    
    integrals <- 
      c(
        S_A = S_A, 
        I_A = I_A, 
        R_A = R_A
      )
    
    velocity <- 
      c(
        t_LE = t_LE, 
        t_FE = t_FE, 
        diff = t_FE-t_LE
      )
    
    occ_regions <- data.frame(first_occ_regions$date, first_occ_regions[,N+2], last_occ_regions[,N+2])
    colnames(occ_regions) <- c("date", "LE", "FE")
    occ_regions$LE_FE <- occ_regions$LE-occ_regions$FE
    
    SIR_regions <- data.frame(matrix(ncol = 4, nrow = nrow(occ_regions)))
    SIR_regions[,1] <- occ_regions$date
    SIR_regions[,2] <- N-cumsum(occ_regions$LE)
    SIR_regions[,3] <- cumsum(occ_regions$LE_FE)
    SIR_regions[,4] <- cumsum(occ_regions$FE)
    colnames(SIR_regions) <-
      c(
        "date",
        "susceptible",
        "infected",
        "recovered"
      )
    
    cases_by_date <- aggregate(data[[col_cases]], by = list(data[[col_date]]), FUN = sum)
    colnames(cases_by_date) <- c("date", "cases")
    
    cases_by_region <- aggregate(data[[col_cases]], by = list(data[[col_region]]), FUN = sum)
    colnames(cases_by_region) <- c("region", "cases_cumulative")
    
    data_statistics <- 
      c(
        N, 
        TP, 
        N_withoutcases, 
        data_balanced
      )
    
    col_names = c(col_cases, col_date, col_region)
    
    if(isTRUE(verbose)) {
      cat("OK", "\n")
    }
    
    sbm_object <- 
      new(
        "sbm", 
        R_0A = R_0A, 
        integrals = integrals, 
        velocity = velocity,
        occ_regions = occ_regions,
        SIR_regions = SIR_regions,
        cases_by_date = cases_by_date,
        cases_by_region = cases_by_region,
        input_data = data.frame(data),
        data_statistics = data_statistics,
        col_names = col_names
      )
    
    sbm_object <- add_timestamp(
      sbm_object,
      function_or_method = "swash_backwash",
      process = paste0("Calculation of Swash-Backwash Model for ", N, " regions and ", TP, " time points")
    )
    
    invisible(sbm_object)
    
  }


compare_countries <-
  function(
    sbm1,
    sbm2,
    country_names = c("Country 1", "Country 2"),
    indicator = "R_0A",
    iterations = 20,
    samples_ratio = 0.8,
    alpha = 0.05,
    replace = TRUE
  ) {
    
    country_names <- as.character(country_names)
    
    sbm1_confint <- 
      confint(
        sbm1,
        iterations = iterations,
        samples_ratio = samples_ratio,
        alpha = alpha,
        replace = replace
      )
    
    sbm2_confint <- 
      confint(
        sbm2,
        iterations = iterations,
        samples_ratio = samples_ratio,
        alpha = alpha,
        replace = replace
      )
    
    D <- sbm1_confint@iterations[[indicator]]-sbm2_confint@iterations[[indicator]]
    
    D_ci <- quantile_ci(
      x = D, 
      alpha = alpha
    )
    
    bootstrap_config <- list(
      iterations = iterations,
      samples_ratio = samples_ratio,
      alpha = alpha,
      replace = replace
    )
    
    new("countries", 
        sbm_ci1 = sbm1_confint, 
        sbm_ci2 = sbm2_confint,  
        D = D,
        D_ci = D_ci,
        config = bootstrap_config,
        country_names = country_names,
        indicator = indicator
    )
    
  }