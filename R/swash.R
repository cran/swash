#-------------------------------------------------------------------------------
# Name:        swash
# Purpose:     Swash-Backwash Model for the Single Epidemic Wave
# Author:      Thomas Wieland (geowieland@googlemail.com)
# Version:     1.1.0
# Last update: 22.02.2025 12:56
# Copyright (c) 2025 Thomas Wieland
#-------------------------------------------------------------------------------


setClass("sbm",
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
           col_names = "character"
         ))


setClass("sbm_ci",
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
         ))


setClass("countries",
         slots = list(
           sbm_ci1 = "sbm_ci",
           sbm_ci2 = "sbm_ci",
           D = "numeric",
           D_ci = "numeric",
           config = "list",
           country_names = "character",
           indicator = "character"
         ))


swash <- 
  function(
    data,
    col_cases, 
    col_date, 
    col_region
    ) {
    
    par_old <- par(no.readonly = TRUE)
    on.exit(par(par_old))
    
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
        col_date = col_date)
    data_balanced <- data_check_balanced$data_balanced

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

    integrals <- c (S_A = S_A, I_A = I_A, R_A = R_A)
    velocity <- c (t_LE = t_LE, t_FE = t_FE, diff = t_FE-t_LE)

    occ_regions <- data.frame(first_occ_regions$date, first_occ_regions[,N+2], last_occ_regions[,N+2])
    colnames(occ_regions) <- c("date", "LE", "FE")
    occ_regions$LE_FE <- occ_regions$LE-occ_regions$FE

    SIR_regions <- data.frame(matrix(ncol = 4, nrow = nrow(occ_regions)))
    SIR_regions[,1] <- occ_regions$date
    SIR_regions[,2] <- N-cumsum(occ_regions$LE)
    SIR_regions[,3] <- cumsum(occ_regions$LE_FE)
    SIR_regions[,4] <- cumsum(occ_regions$FE)
    colnames(SIR_regions) <-
      c("date",
        "susceptible",
        "infected",
        "recovered")
    
    
    cases_by_date <- aggregate(data[[col_cases]], by = list(data[[col_date]]), FUN = sum)
    colnames(cases_by_date) <- c("date", "cases")

    cases_by_region <- aggregate(data[[col_cases]], by = list(data[[col_region]]), FUN = sum)
    colnames(cases_by_region) <- c("region", "cases_cumulative")

    data_statistics <- c(N, TP, N_withoutcases, data_balanced)

    col_names = c(col_cases, col_date, col_region)

    new("sbm", 
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

  }


setMethod(
  "summary", 
  "sbm", 
  function(object) {
    
    cat("Swash-Backwash Model", "\n")
    
    results_df <- data.frame(matrix(ncol = 1, nrow = 10))
    results_df[1,1] <- ""
    results_df[2,1] <- round (object@integrals[1], 3)
    results_df[3,1] <- round (object@integrals[2], 3)
    results_df[4,1] <- round (object@integrals[3], 3)
    results_df[5,1] <- ""
    
    results_df[6,1] <- ""
    results_df[7,1] <- round (object@velocity[1], 3)
    results_df[8,1] <- round (object@velocity[2], 3)
    results_df[9,1] <- ""
    
    results_df[10,1] <- round (object@R_0A, 3)
    
    rownames(results_df) <- c(
      "Integrals",
      "Susceptible areas", 
      "Infected areas", 
      "Recovered areas",
      "",
      "Velocity",
      "Leading edge",
      "Following edge",
      " ",
      "Spatial reproduction number"
    )
    
    colnames(results_df) <- " "

    print(results_df)
    cat("\n")
    
    cat ("Input data", "\n")
    cat (paste0("Units       ", object@data_statistics[1]), "\n")
    cat (paste0("No-case     ", object@data_statistics[3]), "\n")
    cat (paste0("Time points ", object@data_statistics[2], "\n"))
    if (object@data_statistics[4] == TRUE) {
      cat("Balanced    YES", "\n")
    } else {
      cat("Balanced    NOPE", "\n")
    }
    
  })


setMethod(
  "print", 
  "sbm", 
  function(x) {
    
    cat(paste0("Swash-Backwash Model with ", x@data_statistics[1], " spatial units and ", x@data_statistics[2], " time points"), "\n")
    cat ("Use summary() for results")

})


setMethod(
  "show", 
  "sbm", 
  function(object) {
    
    cat(paste0("Swash-Backwash Model with ", object@data_statistics[1], " spatial units and ", object@data_statistics[2], " time points"), "\n")
    cat ("Use summary() for results")
    
  })


setMethod(
  "plot", 
  "sbm", 
  function(x, y = NULL, 
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
    
    par(par_old)
  }
)


setGeneric("plot_regions", function(
    object,
    col = "red",
    scale = FALSE,
    normalize_by_col = NULL,
    normalize_factor = 1
    ) {
    standardGeneric("plot_regions")
  })

setMethod(
  "plot_regions", 
  "sbm", 
  function(object,
           col = "red",
           scale = FALSE,
           normalize_by_col = NULL,
           normalize_factor = 1
           ) {

    par_old <- par(no.readonly = TRUE)
    on.exit(par(par_old))

    N <- object@data_statistics[1]

    plot_cols <- 4
    plot_rows <- ceiling(N/plot_cols)
    par (mfrow = c(plot_rows, plot_cols)) 

    input_data <- object@input_data

    col_names <- object@col_names
    col_cases <- col_names[1]
    col_date <- col_names[2]
    col_region <- col_names[3]

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
          col = col, 
          main = N_names[i],
          type = "l",
          ylim = c(0, y_max)
        )
        
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
        
      }
      
    }
    
    par(par_old) 
  }
)


setMethod(
  "confint", 
  "sbm", 
  function(object,
           iterations = 100,
           samples_ratio = 0.8,
           alpha = 0.05,
           replace = TRUE) {

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
        swash (
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
      c("iteration", "S_A", "I_A", "R_A", "t_LE", "t_FE", "t_FE-t_LE", "R_0A", "sample_size")
    
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

    new("sbm_ci", 
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

  }
)


setMethod(
  "summary", 
  "sbm_ci", 
  function(object) {
  
    cis <- names(object@integrals_ci$S_A_ci)
    

    ci_df <- data.frame(matrix(ncol = 2, nrow = 11))
    ci_df[1,] <- cis
    ci_df[2,1] <- round(object@integrals_ci$S_A_ci[1], 3)
    ci_df[2,2] <- round(object@integrals_ci$S_A_ci[2], 3)
    ci_df[3,1] <- round(object@integrals_ci$I_A_ci[1], 3)
    ci_df[3,2] <- round(object@integrals_ci$I_A_ci[2], 3)
    ci_df[4,1] <- round(object@integrals_ci$R_A_ci[1], 3)
    ci_df[4,2] <- round(object@integrals_ci$R_A_ci[2], 3)
    ci_df[5,1:2] <- " "
    
    ci_df[6,] <- cis
    ci_df[7,1] <- round(object@velocity_ci$t_LE_ci[1], 3)
    ci_df[7,2] <- round(object@velocity_ci$t_LE_ci[2], 3)
    ci_df[8,1] <- round(object@velocity_ci$t_FE_ci[1], 3)
    ci_df[8,2] <- round(object@velocity_ci$t_FE_ci[2], 3)
    ci_df[9,1:2] <- " "

    ci_df[10,] <- cis
    ci_df[11,1] <- round(object@R_0A_ci[1], 3)
    ci_df[11,2] <- round(object@R_0A_ci[2], 3)

    rownames(ci_df) <- c(
      "Integrals",
      "Susceptible areas", 
      "Infected areas", 
      "Recovered areas",
      "",
      "Velocity",
      "Leading edge",
      "Following edge",
      " ",
      "   ",
      "Spatial reproduction number"
      )
    
    colnames(ci_df) <- c(" "," ")

    ci_df <- format(ci_df, justify = "right")

    cat("Confidence Intervals for Swash-Backwash Model", "\n")
    
    print(as.matrix(ci_df), quote = FALSE, colnames = FALSE)

    cat ("\n")
    cat ("Configuration for confidence intervals", "\n")
    cat (paste0("CI alpha    ", object@config$alpha), "\n")
    cat (paste0("Sample      ", object@config$samples_ratio*100, " % (", object@config$regions_to_sample, " units)"), "\n")
    cat (paste0("Iterations  ", object@config$iterations, "\n"))
    if (object@config$replace == TRUE) {
      cat ("Bootstrap   YES", "\n")
    } else {
      cat ("Bootstrap   NO", "\n")
    }

    cat ("\n")
    cat ("Input data", "\n")
    cat (paste0("Units       ", object@data_statistics[1]), "\n")
    cat (paste0("No-case     ", object@data_statistics[3]), "\n")
    cat (paste0("Time points ", object@data_statistics[2], "\n"))
    if (object@data_statistics[4] == TRUE) {
      cat("Balanced    YES", "\n")
    } else {
      cat("Balanced    NO", "\n")
    } 
  })


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
  function(x, y = NULL, 
           col_bars = "grey",
           col_ci = "red"
  ) {

    par_old <- par(no.readonly = TRUE)
    on.exit(par(par_old)) 
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
    
    par(par_old)
  }
)


compare_countries <-
  function(
    sbm1,
    sbm2,
    country_names = c("Country 1", "Country 2"),
    indicator = "R_0A",
    iterations = 100,
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

    D_ci <- quantile_ci(x = D, alpha = alpha)

    
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


setMethod(
  "plot", 
  "countries", 
  function(x, y = NULL, 
           col_bars = "grey",
           col_ci = "red"
  ) {

    par_old <- par(no.readonly = TRUE)
    on.exit(par(par_old)) 
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

    par(par_old)
  }
)


setMethod(
  "summary", 
  "countries", 
  function(object) {
    
    D_ci <- object@D_ci
    cis <- names(D_ci)
    indicator <- object@indicator
    
    D_mean <- mean(object@D, na.rm = TRUE)
    D_median <- median(object@D, na.rm = TRUE)
    
    
    ci_df <- data.frame(matrix(ncol = 4, nrow = 1))
    ci_df[1,1] <- round(D_mean, 3)
    ci_df[1,2] <- round(D_median, 3)
    ci_df[1,3:4] <- round(D_ci, 3)
    colnames(ci_df) <- c("Mean", "Median", cis)
    rownames(ci_df) <- paste0("Difference in ", indicator)
    
    cat("Two-country comparison for Swash-Backwash Model", "\n")
    cat("\n")
    
    print(ci_df)
    
    cat ("\n")
    cat ("Configuration for confidence intervals", "\n")
    cat (paste0("CI alpha    ", object@config$alpha), "\n")
    cat (paste0("Sample      ", object@config$samples_ratio*100, " % "), "\n")
    cat (paste0("Iterations  ", object@config$iterations), "\n")
    if (object@config$replace == TRUE) {
      cat ("Bootstrap   YES", "\n")
    } else {
      cat ("Bootstrap   NO", "\n")
    }
    
  })



setMethod(
  "show", 
  "countries", 
  function(object) {
    
    cat(paste0("Two-country comparison with Swash-Backwash Model"), "\n")
    cat ("Use summary() for results", "\n")
    
  })


quantile_ci <-
  function(
    x,
    alpha = 0.05
  ) {
    
    ci_lower <- alpha/2
    ci_upper <- 1-(alpha/2)

    ci <- quantile(x, probs = c(ci_lower, ci_upper))
    
    return(ci)
    
  }



hist_ci <-
  function(x,
           alpha = 0.05,
           col_bars = "grey",
           col_ci = "red",
           ...) {
    
    ci <- quantile_ci(
      x = x, 
      alpha = alpha
      )
    
    hist(x, col = col_bars, ...)
    abline(v = ci[1], col = col_ci)
    abline(v = ci[2], col = col_ci)
    abline(v = median(x), col = col_ci)
    
  }


is_balanced <-
  function (
    data,
    col_cases, 
    col_date, 
    col_region,
    as_balanced = TRUE,
    fill_missing = 0
  ) {
    
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
    
    results <- list (data_balanced = data_balanced, data = data)
    
    return (results)
    
  }


as_balanced <-
  function(
    data,
    col_cases, 
    col_date, 
    col_region,
    fill_missing = 0) {
    
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
        by.y = paste0("__", col_region, "_x_", col_date, "__")
      )

    data[[col_region]] <- data[[paste0("__", col_region, "__")]]
    data[[col_date]] <- data[[paste0("__", col_date(), "__")]]
    if (nrow(data[is.na(data[[col_cases]]),]) > 0) {
      data[is.na(data[[col_cases]]),][[col_cases]] <- fill_missing
    }
    
    data[[paste0("__", col_region, "_x_", col_date, "__")]] <- NULL
    data[[paste0("__", col_region, "__")]] <- NULL
    data[[paste0("__", col_date, "__")]] <- NULL

    return(data)

  }