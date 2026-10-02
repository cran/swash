#---------------------------------------------------------------
# Name:        growthmodels (swash package)
# Purpose:     Logistic and exponential growth model
#              and Hawkes process
#              Functions, classes and corresponding methods
# Author:      Thomas Wieland 
#              ORCID: 0000-0001-5168-9846
#              mail: geowieland@googlemail.com
# Version:     2.0.0
# Last update: 2026-09-08 19:11
# Copyright (c) 2020-2026 Thomas Wieland
#---------------------------------------------------------------


# Class growthmodels:
setClass(
  "growthmodels",
  slots = list(
    results = "data.frame",
    growth_models = "list",
    model_type = "character",
    results_cols = "character",
    results_cols_names = "character",
    data_statistics = "numeric",
    time_format = "character",
    timestamp = "list"
  ),
  prototype = list(
    timestamp = list()
  )  
)

# Methods for class growthmodels:
setMethod(
  "summary",
  "growthmodels",
  function(object) {
    
    model_type  <- object@model_type
    time_format <- object@time_format
    results_cols <- object@results_cols
    results_cols_names <- object@results_cols_names
    
    cat(paste0(model_type, " Growth Model"), "\n\n")
    
    cat("Input data\n")
    cat(sprintf("  Units       : %s\n", object@data_statistics[1]))
    cat(sprintf("  Time points : %s\n", object@data_statistics[2]))
    cat(sprintf(
      "  Balanced    : %s\n\n",
      ifelse(object@data_statistics[4], "YES", "NO")
    )
    )
    
    results_show <- object@results[results_cols]
    
    results_show[, 2:ncol(results_show)] <- round(results_show[, 2:ncol(results_show)], 3)
    colnames(results_show) <- results_cols_names
    
    cat("Results per region")
    
    if (nrow(results_show) > 10) {
      
      cat(paste0(" (Showing first and last 5 of ", object@data_statistics[1], " cases)"), "\n\n")
      
      print(head(results_show, 5), row.names = FALSE)
      cat("...", "\n")
      print(tail(results_show, 5), row.names = FALSE)
      
      cat("\n")
      cat("Use your_growth_models@results to access the full table", "\n")
      
    } else {
      
      cat("\n\n")
      print(results_show[results_cols], row.names = FALSE)
      
    }
    
    invisible(object)
    
  }
)

setMethod(
  "print", 
  "growthmodels", 
  function(x) {
    
    cat(paste0(x@model_type, " Growth Models for ", x@data_statistics[1], " spatial units and ", x@data_statistics[2], " time points"), "\n")
    cat ("Use summary() for results", "\n")
    
    invisible(x)
    
  }
)

setMethod(
  "show", 
  "growthmodels", 
  function(object) {
    
    cat(paste0(object@model_type, " Growth Models for ", object@data_statistics[1], " spatial units and ", object@data_statistics[2], " time points"), "\n")
    cat ("Use summary() for results", "\n")
    
    invisible(object)
    
  }
)


setClass(
  "expgrowth",
  slots = list(
    GrowthModel_OLS = "list",
    GrowthModel_NLS = "list",
    y = "numeric",
    t = "numeric",
    config = "list"
  )
)

setMethod(
  "summary", 
  "expgrowth", 
  function(object) {
    
    GrowthModel_OLS <- object@GrowthModel_OLS
    GrowthModel_NLS <- object@GrowthModel_NLS
    config <- object@config
    
    cat(model_descriptions[["expgrowth"]], "\n\n")
    
    cat("OLS model\n")
    cat(sprintf("  Growth rate              : %.3f\n", GrowthModel_OLS$exp_gr))
    cat(sprintf("  Baseline                 : %.3f\n", GrowthModel_OLS$y_0))
    cat(sprintf("  Basic reproduction number: %.3f\n", GrowthModel_OLS$R0))
    cat(sprintf("  Doubling rate            : %.3f\n", GrowthModel_OLS$doubling))
    cat(sprintf("  R-squared                : %.3f\n\n", GrowthModel_OLS$fit_metrics$R2))
    
    cat("NLS model\n")
    cat(sprintf("  Growth rate              : %.3f\n", GrowthModel_NLS$exp_gr))
    cat(sprintf("  Baseline                 : %.3f\n", GrowthModel_NLS$y_0))
    cat(sprintf("  Basic reproduction number: %.3f\n", GrowthModel_NLS$R0))
    cat(sprintf("  Doubling rate            : %.3f\n", GrowthModel_NLS$doubling))
    cat(sprintf("  R-squared                : %.3f\n\n", GrowthModel_NLS$fit_metrics$R2))
    
    cat("Input\n")
    
    cat(sprintf("  Time points (not NaN)    : %.0f\n", as.integer(object@config$TP)))
    
    
    invisible(object)
    
  }
)

setMethod(
  "plot", 
  "expgrowth", 
  function(
    x, 
    y = NULL,
    cp_col = "lightblue",
    cp_pch = 19,
    cl_col = "red",
    x_lab = "Time",
    y_lab = "Daily infections",
    plot_title = "Exponential growth model",
    text_size = 1,
    anno_size = 0.7,
    bgrid = TRUE,
    bgrid_col = "white",
    bgrid_type = 1,
    bgrid_size = 1,
    bg_col = "white",
    show_model_results = TRUE
  ) {
    
    par_old <- par(no.readonly = TRUE)
    on.exit(par(par_old))
    dev.new()
    
    t <- x@t
    y <- x@y
    GrowthModel_NLS <- x@GrowthModel_NLS
    GrowthModel_OLS <- x@GrowthModel_OLS
    
    if (is.null(GrowthModel_NLS)) {
      
      y_pred <- GrowthModel_OLS$y_pred
      exp_gr <- GrowthModel_OLS$exp_gr
      y_0 <- GrowthModel_OLS$y_0
      R0 <- GrowthModel_OLS$R0
      doubling <- GrowthModel_OLS$doubling
      sum_of_squares <- GrowthModel_OLS$fit_metrics$SQT
      R2 <- GrowthModel_OLS$fit_metrics$R2
      MAPE <- GrowthModel_OLS$fit_metrics$MAPE
      
      model = "OLS"
      
    } else {
      
      y_pred <- GrowthModel_NLS$y_pred
      exp_gr <- GrowthModel_NLS$exp_gr
      y_0 <- GrowthModel_NLS$y_0
      R0 <- GrowthModel_NLS$R0
      doubling <- GrowthModel_NLS$doubling
      sum_of_squares <- GrowthModel_NLS$fit_metrics$SQT
      R2 <- GrowthModel_NLS$fit_metrics$R2
      MAPE <- GrowthModel_NLS$fit_metrics$MAPE
      
      model = "NLS"
      
    }
    
    par(mar=c(5.1, 5, 4.1, 4.1)) 
    
    plot(
      t, 
      y, 
      col = cp_col, 
      "p", 
      pch = cp_pch, 
      xlab = x_lab, 
      ylab = y_lab, 
      main = plot_title, 
      cex.axis = text_size, 
      cex.lab = text_size, 
      cex.main = text_size
    )
    
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = bg_col)
    
    if (bgrid == TRUE) {
      grid(
        col = bgrid_col, 
        lty = bgrid_type, 
        lwd = bgrid_size
      )
    }
    
    points(
      t, 
      y, 
      col = cp_col, 
      pch = 19, 
      xlab = x_lab, 
      ylab = y_lab, 
      main = plot_title, 
      cex.axis = text_size, 
      cex.lab = text_size, 
      cex.main = text_size
    )
    
    lines(
      t, 
      y_pred, 
      col = cl_col
    )
    
    if(isTRUE(show_model_results)) {
      text (paste0("r = ", round(exp_gr, 5)), x = 0, y = max(y), pos = 4, cex = anno_size)
      text (paste0("y_0 = ", round(y_0, 3)), x = 0, y = (max(y)-(max(y)*0.04*anno_size)), pos = 4, cex = anno_size)
      text (paste0("R0 = ", round(R0, 3)), x = 0, y = (max(y)-(max(y)*0.12*anno_size)), pos = 4, cex = anno_size)
      text (paste0("doubling = ", round(doubling, 3)), x = 0, y = (max(y)-(max(y)*0.16*anno_size)), pos = 4, cex = anno_size)
      text (paste0("SSQ = ", round(sum_of_squares, 3)), x = 0, y = (max(y)-(max(y)*0.20*anno_size)), pos = 4, cex = anno_size)
      text (paste0("R^2 = ", round(R2, 3)), x = 0, y = (max(y)-(max(y)*0.24*anno_size)), pos = 4, cex = anno_size)
      text (paste0("MAPE = ", round(MAPE, 3)), x = 0, y = (max(y)-(max(y)*0.28*anno_size)), pos = 4, cex = anno_size)
    }
    
    invisible(x)
  }
)

setMethod(
  "print", 
  "expgrowth", 
  function(x) {
    
    cat(model_descriptions[["expgrowth"]], "\n")
    cat ("Use summary() for results", "\n")
    
  }
)

setMethod(
  "show", 
  "expgrowth", 
  function(object) {
    
    cat(model_descriptions[["expgrowth"]], "\n")
    cat ("Use summary() for results", "\n")
    
  }
)


setClass(
  "loggrowth",
  slots = list(
    LinModel = "list", 
    GrowthModel_OLS = "list", 
    GrowthModel_NLS = "list",
    y = "numeric",
    t = "numeric",
    config = "list"
  )
)

setMethod(
  "summary", 
  "loggrowth", 
  function(object) {
    
    GrowthModel_OLS <- object@GrowthModel_OLS
    GrowthModel_NLS <- object@GrowthModel_NLS
    config <- object@config
    
    cat(model_descriptions[["loggrowth"]], "\n\n")
    
    cat("OLS model\n")
    cat(sprintf("  Saturation               : %.3f\n", GrowthModel_OLS$S))
    cat(sprintf("  Growth rate              : %.3f\n", GrowthModel_OLS$r))
    cat(sprintf("  Baseline                 : %.3f\n", GrowthModel_OLS$y_0))
    cat(sprintf("  Inflection point         : %.3f\n", GrowthModel_OLS$ip))
    cat(sprintf("  Time of inflection point : %.3f\n", GrowthModel_OLS$t_ip))
    cat(sprintf("  R-squared                : %.3f\n\n", GrowthModel_OLS$fit_metrics$R2))
    
    if (length(GrowthModel_NLS) > 0) {
      
      cat("NLS model\n")
      cat(sprintf("  Saturation               : %.3f\n", GrowthModel_NLS$S))
      cat(sprintf("  Growth rate              : %.3f\n", GrowthModel_NLS$r))
      cat(sprintf("  Baseline                 : %.3f\n", GrowthModel_NLS$y_0))
      cat(sprintf("  Inflection point         : %.3f\n", GrowthModel_NLS$ip))
      cat(sprintf("  Time of inflection point : %.3f\n", GrowthModel_NLS$t_ip))
      cat(sprintf("  R-squared                : %.3f\n\n", GrowthModel_NLS$fit_metrics$R2))
      
    }
    
    cat("Input\n")
    
    cat(sprintf("  Time points (not NaN)    : %.0f\n", as.integer(object@config$TP)))
    
    if (config$nls == FALSE) {
      
      cat(sprintf("  NLS estimation         : Not desired by user\n"))
      
    } else {
      
      if (isTRUE(config$nls_estimation)) {
        
        if (!is.null(config$S)) {
          
          cat(paste0("  Saturation estimation    : ", config$S_start_est_method), "\n")
          cat(paste0("  Saturation value         : ", config$S), "\n")
          
        } else {
          
          cat(paste0("  Saturation estimation    : ", config$S_start_est_method), "\n")
          cat(paste0("  Saturation start value   : ", config$S_start), "\n")
          cat(paste0("  Saturation end value     : ", config$S_end), "\n")
          
        }
        
      } else {
        
        cat("  NLS estimation         : Failed\n")
        
      }
    }
    
    invisible(object)
    
  }
  
)

setMethod(
  "plot", 
  "loggrowth", 
  function(
    x, 
    y = NULL,
    cp_col = "lightblue",
    cp_pch = 19,
    cl_col = "red",
    plot_d = TRUE,
    dl_col = "blue",
    x_lab = "Time",
    y_lab = "Cumulative infections",
    y2_lab = "dC/dt",
    plot_title = "Logistic growth model",
    text_size = 1,
    anno_size = 0.7,
    bgrid = TRUE,
    bgrid_col = "white",
    bgrid_type = 1,
    bgrid_size = 1,
    bg_col = "white",
    show_model_results = TRUE
  ) {
    
    par_old <- par(no.readonly = TRUE)
    on.exit(par(par_old))
    dev.new()
    
    t <- x@t
    y <- x@y
    GrowthModel_NLS <- x@GrowthModel_NLS
    GrowthModel_OLS <- x@GrowthModel_OLS
    
    if (is.null(GrowthModel_NLS)) {
      
      y_pred <- GrowthModel_OLS$y_pred
      dy_dt <- GrowthModel_OLS$dy_dt
      r <- GrowthModel_OLS$r
      y_0 <- GrowthModel_OLS$y_0
      S <- GrowthModel_OLS$S
      ip <- GrowthModel_OLS$ip
      t_ip <- GrowthModel_OLS$t_ip
      sum_of_squares <- GrowthModel_OLS$fit_metrics$SQT
      R2 <- GrowthModel_OLS$fit_metrics$R2
      MAPE <- GrowthModel_OLS$fit_metrics$MAPE
      
      model = "OLS"
      
    } else {
      
      y_pred <- GrowthModel_NLS$y_pred
      dy_dt <- GrowthModel_NLS$dy_dt
      r <- GrowthModel_NLS$r
      y_0 <- GrowthModel_NLS$y_0
      S <- GrowthModel_NLS$S
      ip <- GrowthModel_NLS$ip
      t_ip <- GrowthModel_NLS$t_ip
      sum_of_squares <- GrowthModel_NLS$fit_metrics$SQT
      R2 <- GrowthModel_NLS$fit_metrics$R2
      MAPE <- GrowthModel_NLS$fit_metrics$MAPE
      
      model = "NLS"
      
    }
    
    par(mar=c(5.1, 5, 4.1, 4.1)) 
    
    plot(
      t, 
      y, 
      col = cp_col, 
      "p", 
      pch = cp_pch, 
      xlab = x_lab, 
      ylab = y_lab, 
      main = plot_title, 
      cex.axis = text_size, 
      cex.lab = text_size, 
      cex.main = text_size
    )
    
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = bg_col)
    
    if (bgrid == TRUE) {
      grid(
        col = bgrid_col, 
        lty = bgrid_type, 
        lwd = bgrid_size
      )
    }
    
    points(
      t, 
      y, 
      col = cp_col, 
      pch = 19, 
      xlab = x_lab, 
      ylab = y_lab, 
      main = plot_title, 
      cex.axis = text_size, 
      cex.lab = text_size, 
      cex.main = text_size
    )
    
    
    lines(
      t, 
      y_pred, 
      col = cl_col
    )
    
    if(isTRUE(show_model_results)) {
      text (paste0("r = ", round(r, 5)), x = 0, y = max(y), pos = 4, cex = anno_size)
      text (paste0("y_0 = ", round(y_0, 3)), x = 0, y = (max(y)-(max(y)*0.04*anno_size)), pos = 4, cex = anno_size)
      text (paste0("S = ", round(S, 3)), x = 0, y = (max(y)-(max(y)*0.08*anno_size)), pos = 4, cex = anno_size)
      text (paste0("ip = ", round(ip, 3)), x = 0, y = (max(y)-(max(y)*0.12*anno_size)), pos = 4, cex = anno_size)
      text (paste0("t_ip = ", round(t_ip, 3)), x = 0, y = (max(y)-(max(y)*0.16*anno_size)), pos = 4, cex = anno_size)
      text (paste0("SSQ = ", round(sum_of_squares, 3)), x = 0, y = (max(y)-(max(y)*0.20*anno_size)), pos = 4, cex = anno_size)
      text (paste0("R^2 = ", round(R2, 3)), x = 0, y = (max(y)-(max(y)*0.24*anno_size)), pos = 4, cex = anno_size)
      text (paste0("MAPE = ", round(MAPE, 3)), x = 0, y = (max(y)-(max(y)*0.28*anno_size)), pos = 4, cex = anno_size)
    }
    
    if (isTRUE(plot_d)) {
      
      par(new = TRUE)
      
      plot (
        t, 
        dy_dt, 
        xaxt = "n", 
        yaxt = "n",
        ylab = "", 
        xlab = "", 
        col = dl_col, 
        type = "l",
        cex.axis = text_size, 
        cex.lab = text_size
      )
      
      axis(side = 4, cex.axis = text_size)
      mtext(y2_lab, side = 4, line = 3, cex = text_size)
      
    }
    
    invisible(x)
    
  }
)

setMethod(
  "print", 
  "loggrowth", 
  function(x) {
    
    cat(model_descriptions[["loggrowth"]], "\n")
    cat ("Use summary() for results", "\n")
    
    invisible(x)
    
  }
)

setMethod(
  "show", 
  "loggrowth", 
  function(object) {
    
    cat(model_descriptions[["loggrowth"]], "\n")
    cat ("Use summary() for results", "\n")
    
    invisible(object)
    
  }
)


setClass(
  "hawkes",
  slots = list(
    y = "numeric",
    t = "numeric",
    mu = "numeric",
    alpha = "numeric",
    beta = "numeric",
    br = "numeric",
    y_pred = "numeric",
    fit_metrics = "list",
    config = "list"
  )
)

setMethod(
  "summary", 
  "hawkes", 
  function(object) {
    
    cat(model_descriptions[["hawkes"]], "\n\n")
    
    mu <- object@mu
    alpha <- object@alpha
    beta <- object@beta
    br <- object@br
    fit_metrics <- object@fit_metrics
    config <- object@config
    
    cat("Model optimization\n")
    cat(sprintf("  Mu                       : %.3f\n", mu))
    cat(sprintf("  Alpha                    : %.3f\n", alpha))
    cat(sprintf("  Beta                     : %.3f\n", beta))
    cat(sprintf("  Breaking ratio           : %.3f\n\n", br))
    
    cat(sprintf("  Sum of squares           : %.3f\n", fit_metrics$SQR))
    cat(sprintf("  R-squared                : %.3f\n", fit_metrics$R2))
    cat(sprintf("  MAPE                     : %.3f\n\n", fit_metrics$MAPE))
    
    cat(sprintf("  Time points (not NaN)    : %.0f\n", as.integer(config$TP)))
    cat(paste0("  Optimization method      : ", as.character(config$optim_method)), "\n")
    
  }
)

setMethod(
  "plot", 
  "hawkes", 
  function(
    x, 
    y = NULL,
    cp_col = "lightblue",
    cp_pch = 19,
    cl_col = "red",
    x_lab = "Time",
    y_lab = "Daily infections",
    plot_title = "Hawkes process model",
    text_size = 1,
    anno_size = 0.7,
    bgrid = TRUE,
    bgrid_col = "white",
    bgrid_type = 1,
    bgrid_size = 1,
    bg_col = "white",
    show_model_results = TRUE
  ) {
    
    par_old <- par(no.readonly = TRUE)
    on.exit(par(par_old))
    dev.new()
    
    t <- x@t
    y <- x@y
    y_pred <- x@y_pred
    
    mu <- x@mu
    alpha <- x@alpha
    beta <- x@beta
    br <- x@br
    
    fit_metrics <- x@fit_metrics
    sum_of_squares <- fit_metrics$SQR
    R2 <- fit_metrics$R2
    MAPE <- fit_metrics$MAPE
    
    config <- x@config
    
    par(mar=c(5.1, 5, 4.1, 4.1)) 
    
    plot(
      t, 
      y, 
      col = cp_col, 
      "p", 
      pch = cp_pch, 
      xlab = x_lab, 
      ylab = y_lab, 
      main = plot_title, 
      cex.axis = text_size, 
      cex.lab = text_size, 
      cex.main = text_size
    )
    
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = bg_col)
    
    if (bgrid == TRUE) {
      grid(
        col = bgrid_col, 
        lty = bgrid_type, 
        lwd = bgrid_size
      )
    }
    
    points(
      t, 
      y, 
      col = cp_col, 
      pch = 19, 
      xlab = x_lab, 
      ylab = y_lab, 
      main = plot_title, 
      cex.axis = text_size, 
      cex.lab = text_size, 
      cex.main = text_size
    )
    
    lines(
      t, 
      y_pred, 
      col = cl_col
    )
    
    if(isTRUE(show_model_results)) {
      text (paste0("Mu = ", round(mu, 5)), x = 0, y = max(y), pos = 4, cex = anno_size)
      text (paste0("Alpha = ", round(alpha, 3)), x = 0, y = (max(y)-(max(y)*0.04*anno_size)), pos = 4, cex = anno_size)
      text (paste0("Beta = ", round(beta, 3)), x = 0, y = (max(y)-(max(y)*0.08*anno_size)), pos = 4, cex = anno_size)
      text (paste0("br = ", round(br, 3)), x = 0, y = (max(y)-(max(y)*0.12*anno_size)), pos = 4, cex = anno_size)
      text (paste0("SSQ = ", round(sum_of_squares, 3)), x = 0, y = (max(y)-(max(y)*0.20*anno_size)), pos = 4, cex = anno_size)
      text (paste0("R^2 = ", round(R2, 3)), x = 0, y = (max(y)-(max(y)*0.24*anno_size)), pos = 4, cex = anno_size)
      text (paste0("MAPE = ", round(MAPE, 3)), x = 0, y = (max(y)-(max(y)*0.28*anno_size)), pos = 4, cex = anno_size)
    }
    
    invisible(x)
  }
)

setMethod(
  "print", 
  "hawkes", 
  function(x) {
    
    cat(model_descriptions[["hawkes"]], "\n")
    cat ("Use summary() for results", "\n")
    
  }
)

setMethod(
  "show", 
  "hawkes", 
  function(object) {
    
    cat(model_descriptions[["hawkes"]], "\n")
    cat ("Use summary() for results", "\n")
    
  }
)


setClass(
  "breaksgrowth",
  slots = list(
    GrowthModel_OLS = "list",
    y = "numeric",
    t = "numeric",
    config = "list"
  )
)

setMethod(
  "summary", 
  "breaksgrowth", 
  function(object) {
    
    GrowthModel_OLS <- object@GrowthModel_OLS
    config <- object@config
    
    cat("Time Series Model with Breakpoints for Infections\n\n")
    
    cat("OLS model\n")
    
    cat(sprintf("  No. of breakpoints       : %.0f\n", GrowthModel_OLS$bp_obj_breakpoints_no))
    cat(sprintf("  No. of segments          : %.0f\n", GrowthModel_OLS$bp_obj_segments_no))
    cat(paste0("  Breakpoints at positions : ", paste(GrowthModel_OLS$bp_obj_breakpoints, collapse = ", ")), "\n")
    
    cat("Input\n")
    
    cat(sprintf("  Time points (not NaN)    : %.0f\n", as.integer(object@config$TP)))
    if(isTRUE(object@config$ln)) {
      cat("  Natural logarithm        : YES", "\n")
    } else {
      cat("  Natural logarithm        : NO", "\n")
    }
    
    invisible(object)
    
  }
)

setMethod(
  "print", 
  "breaksgrowth", 
  function(x) {
    
    cat(paste0("Time Series Model with Breakpoints for Infections"), "\n")
    cat ("Use summary() for results", "\n")
    
  }
)

setMethod(
  "show", 
  "breaksgrowth", 
  function(object) {
    
    cat(paste0("Time Series Model with Breakpoints for Infections"), "\n")
    cat ("Use summary() for results", "\n")
    
  }
)

setMethod(
  "plot", 
  "breaksgrowth", 
  function(
    x, 
    y = NULL,
    line.col = "chocolate",
    ci.col = c(
      rgb(0, 0, 255, maxColorValue = 255, alpha = 0.5),
      rgb(0, 255, 0, maxColorValue = 255, alpha = 0.5),
      rgb(255, 0, 0, maxColorValue = 255, alpha = 0.5)
    ),
    legend.show = TRUE, 
    legend.pos = c(
      "topright", 
      "topright", 
      "topright", 
      "topright"
    ),
    xlab = "Time", 
    ylab = "Y", 
    ylim = NULL,
    plot.main = c(
      "Breakpoints", 
      "Segment model fits", 
      "BIC and Residual Sum of Squares", 
      "Model without breaks (one segment)"
    ),
    output.full = TRUE,
    separate_plots = FALSE,
    show_model_results = TRUE,
    ...
  ) {
    
    par_old <- par(no.readonly = TRUE)
    on.exit(par(par_old))
    dev.new()
    
    config <- x@config
    ln <- config$ln
    TP <- config$TP
    add_constant <- config$add_constant
    alpha <- config$alpha
    
    segments_models_df <- x@GrowthModel_OLS$segments_models_df
    segment_model_ci_points <- x@GrowthModel_OLS$segment_model_ci_points
    bp_obj_coef <- x@GrowthModel_OLS$bp_obj_coef
    bp_obj_breakpoints <- x@GrowthModel_OLS$bp_obj_breakpoints
    bp_obj_breakpoints_no <- x@GrowthModel_OLS$bp_obj_breakpoints_no
    bp_obj_segments_no <- x@GrowthModel_OLS$bp_obj_segments_no
    bp_obj_breakpoints_ci <- x@GrowthModel_OLS$bp_obj_breakpoints_ci
    bp_obj <- x@GrowthModel_OLS$bp_obj
    
    y <- as.numeric(x@y)
    t <- as.numeric(x@t)
    dataset <- data.frame(y, t)
    
    if (isTRUE(ln)) {
      bp_formula <- as.formula("log(y) ~ t")
    } else {
      bp_formula <- as.formula("y ~ t")
    }
    model_no_bp <- lm(
      bp_formula, 
      data = dataset
    )
    
    ci_lower <- round((alpha/2)*100, 1)
    ci_upper <- round((1-(alpha/2))*100, 1)
    CI <- paste0("CI (", ci_lower, ";", ci_upper, ")")
    
    if ((output.full == TRUE) & (separate_plots == FALSE)) {
      par (mfrow = c(2, 2))
    }
    
    
    col_ci <- ci.col[1]
    
    if (is.null(ylim)) {
      
      plot(
        model_no_bp$model[order(model_no_bp$model[,2]),][,2], 
        model_no_bp$model[order(model_no_bp$model[,2]),][,1], 
        lwd = 1.5, 
        type = "l", 
        col = line.col, 
        xlab = xlab, 
        ylab = ylab,
        main = plot.main[1]
      )
      
    } else {
      
      plot(
        model_no_bp$model[order(model_no_bp$model[,2]),][,2], 
        model_no_bp$model[order(model_no_bp$model[,2]),][,1], 
        lwd = 1.5, 
        type = "l", 
        col = line.col, 
        xlab = xlab, 
        ylab = ylab,
        main = plot.main[1], 
        ylim = ylim
      )
      
    }
    
    if (bp_obj_breakpoints_no > 0) {
      
      i <- 0
      
      for (i in 1:bp_obj_breakpoints_no) {
        
        rect(
          bp_obj_breakpoints_ci[[1]][i], 
          par("usr")[3], 
          bp_obj_breakpoints_ci[[1]][i+((2/3)*length(bp_obj_breakpoints_ci[[1]]))], 
          par("usr")[4], 
          col = col_ci, 
          lty = 0
        )
        
        abline (
          v = bp_obj_breakpoints_ci[[1]][i+((1/3)*length(bp_obj_breakpoints_ci[[1]]))], 
          col = "blue", 
          lty = 1, 
          lwd = 1.5
        )
        
      }
      
      i <- 0
      
      for (i in 1:bp_obj_segments_no) {
        
        text (
          x = mean(c(segments_models_df[i,1], segments_models_df[i,2])), 
          y = min(model_no_bp$model[,1]),
          (round(segments_models_df[i, 7], 3)), 
          cex = 1, 
          col = "black"
        )    
        
      }
      
      
      if (legend.show == TRUE) {
        
        legend(
          legend.pos[1], 
          legend = c("Breakpoint", CI), 
          col = c("blue", col_ci),
          lty = c(1, NA), 
          pch = c(NA, 15), 
          lwd = 1.5, 
          bty = "n", 
          bg = "white"
        )
        
        text(
          min(model_no_bp$model[,2])-max(model_no_bp$model[,2])*0.1, 
          min(model_no_bp$model[,1]), 
          "Slope", 
          xpd = NA
        )
        
      }
      
    }
    
    
    if (output.full == TRUE) {
      
      col_ci <- ci.col[2]
      
      if (is.null(ylim)) {
        
        plot(
          model_no_bp$model[,2], 
          model_no_bp$model[,1], 
          pch = 19, 
          col = line.col, 
          cex = 0.8,
          xlab = xlab, 
          ylab = ylab,
          main = plot.main[2]
        )
        
      }
      
      else {
        
        plot(
          model_no_bp$model[,2], 
          model_no_bp$model[,1], 
          pch = 19, 
          col = line.col, 
          cex = 0.8,
          xlab = xlab, 
          ylab = ylab, 
          main = plot.main[2], 
          ylim = ylim
        )
        
      }
      
      i <- 0
      
      for (i in 1:bp_obj_segments_no) {
        
        pol_coords <- 
          data.frame(
            cbind(
              model_no_bp$model[(segments_models_df[i,1]):(segments_models_df[i,2]), 2],
              rev(model_no_bp$model[(segments_models_df[i,1]):(segments_models_df[i,2]), 2]),
              segment_model_ci_points[[i]][,2]), 
            rev(segment_model_ci_points[[i]][,3]
            )
          )
        colnames(pol_coords) <- c("x1", "x2", "y1", "y2")
        
        polygon(
          c(pol_coords$x1, pol_coords$x2), 
          c(pol_coords$y1, pol_coords$y2),
          col = col_ci, 
          border = NA
        )
        
        lines (
          x = (model_no_bp$model[(segments_models_df[i,1]):(segments_models_df[i,2]), 2]),
          y = (segment_model_ci_points[[i]][,1]), col = "darkgreen", 
          lty = 2, 
          lwd = 1.5
        )
        
        text (
          x = mean(c(segments_models_df[i,1], segments_models_df[i,2])), 
          y = min(model_no_bp$model[,1]),
          (round(segments_models_df[i, 14], 3)), 
          cex = 1, 
          col = "black"
        )
        
      }
      
      if (legend.show == TRUE) {
        
        legend(
          legend.pos[2], 
          legend = c("Model fit", CI), 
          col = c("darkgreen", col_ci),
          lty = c(2, NA), 
          pch = c(NA, 15), 
          lwd = 1.5, 
          bty = "n", 
          bg = "white"
        )
        
        text(
          min(model_no_bp$model[,2])-max(model_no_bp$model[,2])*0.1, 
          min(model_no_bp$model[,1]), 
          expression(R^2), 
          xpd=NA
        )
        
      }
      
      plot(
        bp_obj, 
        lwd = 1.5, 
        main = plot.main[3]
      )
      
      if (is.null(ylim)) {
        
        plot(
          model_no_bp$model[,2], 
          model_no_bp$model[,1], 
          pch = 19, 
          col = line.col, 
          cex = 0.8,
          xlab = xlab, 
          ylab = ylab,
          main = plot.main[4]
        )
      }
      
      else {
        
        plot(
          model_no_bp$model[,2], 
          model_no_bp$model[,1], 
          pch = 19, 
          col = line.col, 
          cex = 0.8,
          xlab = xlab, 
          ylab = ylab,
          main = plot.main[4], 
          ylim = ylim
        )
        
      }
      
      model_no_bp_ci_points <- 
        predict(
          model_no_bp, 
          interval = "confidence", 
          level = 1-alpha
        )
      
      col_ci <- ci.col[3]
      
      pol_coords <- 
        data.frame(
          cbind(
            model_no_bp$model[,2],
            rev(model_no_bp$model[,2]),
            model_no_bp_ci_points[,2], 
            rev(model_no_bp_ci_points[,3])
          )
        )
      
      colnames(pol_coords) <- c("x1", "x2", "y1", "y2")
      
      polygon(
        c(pol_coords$x1, pol_coords$x2), 
        c(pol_coords$y1, pol_coords$y2),
        col = col_ci, 
        border = NA
      )
      
      lines(
        x = model_no_bp$model[,2], 
        y = model_no_bp$fitted.values, 
        lwd = 1.5, 
        lty = 2, 
        col = "red"
      )
      
      if (legend.show == TRUE) {
        
        legend(
          legend.pos[4], 
          legend = c("Model fit", CI),
          col = c("red", col_ci), 
          lty = c(2, NA), 
          pch = c(NA, 15), 
          lwd = 1.5, 
          bty = "n", 
          bg = "white"
        )
        
      }
      
    }
    
    invisible(x)
    
  }
)


exponential_growth <-
  function(
    y, 
    t,
    GI = 4,
    nls = TRUE,
    nls_start = list(a = 1, b = 0.1),
    add_constant = 1,
    verbose = FALSE
  ) {
    
    if (!is.numeric(y)) {
      stop ("Infections vector must be of class 'numeric'.")
    }
    
    if (!is.numeric(t) & !is.Date(t)) {
      stop ("Time vector must be of class 'numeric' or 'Date'.")
    }
    
    if (length(y) != length(t)) {
      stop (paste0("Vectors 'y' and 't' differ in length: ", length(y), ", ", length(t)))
    }
    
    yt_df <- data.frame(y, t)
    yt_df_length_clean <- yt_df[complete.cases(yt_df), ]
    
    if (nrow(yt_df_length_clean) < nrow(yt_df)) {
      
      y <- yt_df$y
      t <- yt_df$t
      
      nan_cases <- nrow(yt_df)-nrow(yt_df_length_clean)
      
      warning(paste0(nan_cases, " NaN values were removed"))
      
    }
    
    TP <- length(y)
    
    class_t <- "numeric"
    
    if (is.Date(t)) {
      
      class_t <- "Date"
      
      message("Time vector is of class 'Date'. Calculating time counter.")
      
      start_date <- min(t)
      
      time_counter <- as.integer(t-start_date)
      
      t <- time_counter
    }
    
    if (any(y <= 0)) {
      
      if (!is.null(add_constant) && is.numeric(add_constant)) {
        
        y_lin <- y+add_constant
        
        message(paste0("Input infections data contains values <= 0. All values are increased by a constant equal to ", add_constant, " for linear estimation."))
        
      } else {
        
        y_lin <- y[y > 0]
        
        warning(paste0("Input infections data contains values <= 0. These values are skipped for linear estimation."))
        
      }
      
    } else {
      
      y_lin <- y
      
    }
    
    if(isTRUE(verbose)) {
      cat("Performing OLS estimation ... ")
    }
    
    log_y_lin <- log(y_lin)
    
    linexpmodel <- lm (log_y_lin ~ t)
    
    y_0 <- linexpmodel$coefficients[1]
    
    exp_gr <- linexpmodel$coefficients[2]
    
    R0 <- exp (exp_gr*GI)
    
    doubling <- log(2)/exp_gr
    
    y_pred <- exp(predict(linexpmodel))
    
    fit_metrics <- metrics(
      observed = y_lin,
      expected = y_pred,
      plot = FALSE
    )
    
    fit_metrics <- fit_metrics[[1]]
    
    model_growth_ols_list <- 
      list (
        exp_gr = exp_gr,
        y_0 = y_0,
        R0 = R0, 
        doubling = doubling,
        y_pred = y_pred,
        model_data = linexpmodel,
        fit_metrics = fit_metrics
      )
    
    if(isTRUE(verbose)) {
      cat("OK", "\n")
    }
    
    model_growth_nls_list <- list()
    nls_estimation <- TRUE
    nls_error_message <- NULL
    
    if (isTRUE(nls)) {
      
      if(isTRUE(verbose)) {
        cat("Trying nonlinear estimation ... ")
      }
      
      tryCatch(
        {
          
          expmodel <- 
            nls(
              y~a*exp(b*t),
              start = nls_start
            )
          
          expmodel_summary <- summary(expmodel)
          
          y_0_NLS <- expmodel_summary$coefficients[1]
          
          exp_gr_NLS <- expmodel_summary$coefficients[2]
          
          R0_NLS <- exp (exp_gr_NLS*GI)
          
          doubling_NLS <- log(2)/exp_gr_NLS
          
          y_pred_NLS <- predict(expmodel)
          
          fit_metrics_NLS <- metrics(
            observed = y,
            expected = y_pred_NLS,
            plot = FALSE
          )
          
          fit_metrics_NLS <- fit_metrics_NLS[[1]]
          
          model_growth_nls_list <- 
            list (
              exp_gr = exp_gr_NLS,
              y_0 = y_0_NLS,
              R0 = R0_NLS, 
              doubling = doubling_NLS,
              y_pred = y_pred_NLS,
              model_data = expmodel,
              fit_metrics = fit_metrics_NLS
            )
          
        },
        
        error = function(cond) {
          
          nls_error_message <<- paste0("Nonlinear estimation failed: ", conditionMessage(cond))
          message(nls_error_message)
          
        }
        
        
      )
      
      if (length(model_growth_nls_list) == 0) {
        nls_estimation <- FALSE
      }
      
      if(isTRUE(verbose)) {
        cat("OK", "\n")
      }
      
    }
    
    config <-
      list (
        TP = TP,
        GI = GI,
        nls = nls,
        nls_start = nls_start,
        nls_estimation = nls_estimation,
        nls_error_message = nls_error_message
      )
    
    new(
      "expgrowth",
      GrowthModel_OLS = model_growth_ols_list, 
      GrowthModel_NLS = model_growth_nls_list,
      y = y,
      t = t,
      config = config
    )
    
  }


logistic_growth <- 
  function (
    y, 
    t, 
    S = NULL,
    S_start = NULL, 
    S_end = NULL, 
    S_iterations = 10, 
    S_start_est_method = "bisect", 
    seq_by = 10,
    nls = TRUE,
    add_constant = 1,
    verbose = FALSE
  ) 
  {
    
    loggrowth_saturation <- 
      function (
    y, 
    t, 
    S_start, 
    S_end, 
    S_start_est_method = "bisect", 
    S_iterations = 10, 
    seq_by = 10,
    verbose = FALSE
      ) {
        
        if(isTRUE(verbose)) {
          cat(paste0("Estimating saturation value via method: ", S_start_est_method, " ... "))
        }
        
        if (S_start_est_method == "trialanderror") {
          
          values_seq <- seq (S_start, S_end, by = seq_by)
          values_no <- length (values_seq)
          values_ssq <- matrix (ncol = 2, nrow = values_no)
          values_ssq[,1] <- values_seq
          
          i <- 0
          
          for (i in 1:values_no) {
            
            y_new <- log ((1/y)-(1/values_seq[i]))
            model_lin <- lm (y_new ~ t)
            model_lin_summary <- summary(model_lin)
            sum_of_squares_lin <- sum(model_lin_summary$residuals^2)
            
            values_ssq[i,2] <- sum_of_squares_lin 
            
          }
          
          values_ssq_order <- values_ssq[order(values_ssq[,2])]
          
          S_est <- values_ssq_order[1]
          
          plot(values_ssq[,1], values_ssq[,2], "l")
          
        }
        
        else {    
          
          i <- 0
          
          interval_m <- vector()
          
          for (i in 1:S_iterations)
          {
            
            interval_m[i] <- (S_start+S_end)/2
            
            y_new_start <- log ((1/y)-(1/S_start))
            
            model_lin_start <- lm (y_new_start ~ t)
            model_lin_start_summary <- summary(model_lin_start)
            sum_of_squares_lin_start <- sum(model_lin_start_summary$residuals^2)
            
            y_new_m <- log ((1/y)-(1/interval_m[i]))
            model_lin_m <- lm (y_new_m ~ t)
            model_lin_m_summary <- summary(model_lin_m)
            sum_of_squares_lin_m <- sum(model_lin_m_summary$residuals^2)
            
            y_new_end <- log ((1/y)-(1/S_end))
            model_lin_end <- lm (y_new_m ~ t)
            model_lin_end_summary <- summary(model_lin_m)
            sum_of_squares_lin_end <- sum(model_lin_end_summary$residuals^2)
            
            if (sum_of_squares_lin_start < sum_of_squares_lin_end)
              
            {
              
              S_start <- S_start
              S_end <- interval_m[i]
              
            }
            
            else
            {
              S_start <- interval_m[i]
              S_end <- S_end
            }
            
          }
          S_est <- interval_m[i]
          
        } 
        
        if(isTRUE(verbose)) {
          cat("OK", "\n")
        }
        
        return(S_est)
        
      }
    
    if (!is.numeric(y)) {
      stop("Cumulative infections vector must be of class 'numeric'.")
    }
    
    if (!is.numeric(t) & !is.Date(t)) {
      stop("Time vector must be of class 'numeric' or 'Date'.")
    }
    
    if (length(y) != length(t)) {
      stop (paste0("Vectors 'y' and 't' differ in length: ", length(y), ", ", length(t)))
    }
    
    yt_df <- data.frame(y, t)
    yt_df_length_clean <- yt_df[complete.cases(yt_df), ]
    
    if (nrow(yt_df_length_clean) < nrow(yt_df)) {
      
      y <- yt_df$y
      t <- yt_df$t
      
      nan_cases <- nrow(yt_df)-nrow(yt_df_length_clean)
      warning(paste0(nan_cases, " NaN values were removed"))
      
    }
    
    if (any(y <= 0)) {
      
      if (!is.null(add_constant) && is.numeric(add_constant)) {
        
        y_lin <- y+add_constant
        
        message(paste0("Input infections data contains values <= 0. All values are increased by a constant equal to ", add_constant, " for linear estimation."))
        
      } else {
        
        y_lin <- y[y > 0]
        
        warning(paste0("Input infections data contains values <= 0. These values are skipped for linear estimation."))
      }
      
    } else {
      
      y_lin <- y
      
    }
    
    TP <- length(y_lin)
    
    class_t <- "numeric"
    
    if (is.Date(t)) {
      
      class_t <- "Date"
      
      message("Time vector is of class 'Date'. Calculating time counter.")
      
      start_date <- min(t)
      
      time_counter <- as.integer(t-start_date)
      
      t <- time_counter
    }
    
    if (is.null(S) & (is.null(S_start) | is.null(S_end))) {
      stop ("Saturation value or start and end values are required for estimation")
    }
    
    if (!is.null(S)) {
      
      if (any(y <= 0) && !is.null(add_constant)) {
        
        S <- S+(add_constant*TP)
        
        message(paste0("Input infections data contains values <= 0 and constant is set to ", add_constant, ". Saturation is increased to ", S, " for linear estimation"))
        
      }
      
      if(S < max(y_lin)) {
        stop (paste0("Saturation value must be above or equal to the the maximum of y (", max(y_lin), ")"))
      }
      
    }
    
    if ((!is.null(S_start) & (!is.null(S_end)))) {
      
      if (any(y <= 0) && !is.null(add_constant)) {
        
        S_start <- S_start+(add_constant*TP)
        S_end <- S_end+(add_constant*TP)
        
        message(paste0("Input infections data contains values <= 0 and constant is set to ", add_constant, ". Saturation start and end values are increased to ", S_start, " and ", S_end, " for linear estimation"))
        
      }
      
      if ((S_start < max(y_lin)) || (S_end < max(y_lin))) {
        stop (paste0("Start and end values of saturation must be above or equal to the the maximum of y (", max(y_lin), ")"))
      }
      
      S <- 
        loggrowth_saturation(
          y = y_lin, 
          t = t, 
          S_start = S_start, 
          S_end = S_end, 
          S_iterations = S_iterations, 
          S_start_est_method = S_start_est_method, 
          seq_by = seq_by,
          verbose = verbose
        )
      
    }
    
    if(isTRUE(verbose)) {
      cat("Performing OLS estimation ... ")
    }
    
    y_new <- log ((1/y_lin)-(1/S))
    
    model_lin <- lm (y_new ~ t)
    
    model_lin_summary <- summary(model_lin)
    
    b <- model_lin_summary$coefficients[1]
    m <- model_lin_summary$coefficients[2]
    
    sum_of_squares_lin <- sum(model_lin_summary$residuals^2)
    
    model_lin_ols_list <- 
      list(
        b = b, 
        m = m, 
        sum_of_squares = sum_of_squares_lin
      )
    
    r <- -m/S
    
    y_0 <- S/(1+S*exp(m*t[1]+b))
    
    ip <- S/2
    c <- -log (y_0/(S-y_0))
    t_ip <- c/(r*S)
    
    y_pred <- S/(1+S*exp(m*t+b))
    
    dy_dt <- r*y_pred*(1-(y_pred/S))
    
    d2y_dt2 <- r*dy_dt*(1-((2*y_pred)/S))
    
    fit_metrics_ols <- metrics(
      observed = y_lin,
      expected = y_pred,
      plot = FALSE
    )
    
    fit_metrics_ols <- fit_metrics_ols[[1]]
    
    model_growth_ols_list <- 
      list (
        S = S, 
        r = r, 
        y_0 = y_0, 
        ip = ip, 
        t_ip = t_ip, 
        y_pred = y_pred,
        fit_metrics = fit_metrics_ols,
        dy_dt = dy_dt,
        d2y_dt2 = d2y_dt2
      )
    
    if(isTRUE(verbose)) {
      cat("OK", "\n")
    }
    
    model_growth_nls_list <- list()
    nls_estimation <- TRUE
    nls_error_message <- NULL
    
    if (nls == TRUE) {
      
      if(isTRUE(verbose)) {
        cat("Trying nonlinear estimation ... ")
      }
      
      tryCatch(
        {
          model_nls <- 
            nls(
              y ~ y_0 * S / (y_0 + (S - y_0) * exp(-r * S * t)), 
              start = list (y_0 = y_0, S = S, r = r), 
              control = list(maxiter = 500)
            )
          
          y_0_nls <-  model_nls$m$getPars()[1]
          
          S_nls <- model_nls$m$getPars()[2]
          
          r_nls <- model_nls$m$getPars()[3]
          
          ip_nls <- S_nls/2
          c_nls <- -log (y_0_nls/(S_nls-y_0_nls))
          t_ip_nls <- c_nls/(r_nls*S_nls)
          
          y_pred <- model_nls$m$predict()
          
          dy_dt <- r_nls*y_pred*(1-(y_pred/S_nls))
          
          d2y_dt2 <- r_nls*dy_dt*(1-((2*y_pred)/S_nls))
          
          fit_metrics_nls <- metrics(
            observed = y,
            expected = y_pred,
            plot = FALSE
          )
          
          fit_metrics_nls <- fit_metrics_nls[[1]]
          
          model_growth_nls_list <- 
            list (
              S = S_nls, 
              r = r_nls, 
              y_0 = y_0_nls, 
              ip = ip_nls, 
              t_ip = t_ip_nls, 
              y_pred = y_pred,
              fit_metrics = fit_metrics_nls,
              dy_dt = dy_dt,
              d2y_dt2 = d2y_dt2
            )
          
          if(isTRUE(verbose)) {
            cat("OK", "\n")
          }
          
        },
        
        error = function(cond) {
          
          nls_error_message <<- paste0("Nonlinear estimation failed: ", conditionMessage(cond))
          message(nls_error_message)
          
        }
        
      )
      
      if (length(model_growth_nls_list) == 0) {
        
        nls_estimation <- FALSE
        
      }
    }
    
    config <-
      list (
        S = S, 
        S_start = S_start,
        S_end = S_end,
        S_iterations = S_iterations,
        S_start_est_method = S_start_est_method,
        seq_by = seq_by,
        nls = nls,
        class_t = class_t,
        TP = TP,
        nls_estimation = nls_estimation,
        nls_error_message = nls_error_message
      )
    
    new(
      "loggrowth", 
      LinModel = model_lin_ols_list, 
      GrowthModel_OLS = model_growth_ols_list, 
      GrowthModel_NLS = model_growth_nls_list,
      y = y,
      t = t,
      config = config
    )
    
  }


hawkes_growth <- 
  function (
    y,
    optim_method = "L-BFGS-B",
    verbose = FALSE
  ) 
  {
    
    if (!is.numeric(y)) {
      stop ("Daily infections vector must be of class 'numeric'.")
    }
    
    TP <- length(y)
    
    loglik <- 
      function(par) {
        
        mu <- par[1]
        alpha <- par[2]
        beta <- par[3]
        
        lambda <- sapply(1:TP, function(t)
          mu + sum(alpha * exp(-beta*(t-1:(t-1)))*y[1:(t-1)])
        )
        
        -sum(dpois(y, lambda, log=TRUE))
      }
    
    fit <- 
      optim(
        c(
          mu=mean(y), 
          alpha=0.2, 
          beta=0.2
        ), 
        loglik, 
        method=optim_method,
        lower=c(
          1e-6,
          1e-6,
          1e-6
        )
      )
    
    mu    <- fit$par["mu"]
    alpha <- fit$par["alpha"]
    beta  <- fit$par["beta"]
    br <- alpha/beta
    
    y_pred <- numeric(TP)
    
    for(t in 1:TP){
      if(t == 1){
        y_pred[t] <- mu
      } else {
        y_pred[t] <- mu + sum(alpha * exp(-beta * (t-1:(t-1))) * y[1:(t-1)])
      }
    }
    
    fit_metrics <- metrics(
      observed = y,
      expected = y_pred,
      plot = FALSE
    )
    
    fit_metrics <- fit_metrics[[1]]
    
    config <-
      list(
        TP = TP,
        optim_method = optim_method
      )
    
    hawkes_object <-
      new(
        "hawkes",
        y = y,
        t = 1:TP,
        mu = mu,
        alpha = alpha,
        beta = beta,
        br = br,
        y_pred = y_pred,
        fit_metrics = fit_metrics,
        config = config
      )
  }


breaks_growth <-
  function(
    y, 
    t,
    ln = FALSE,
    add_constant = 1,
    alpha = 0.05,
    ...,
    verbose = FALSE
  ) {
    
    if (!is.numeric(y)) {
      stop("Infections vector must be of class 'numeric'.")
    }
    
    if (!is.numeric(t) & !is.Date(t)) {
      stop("Time vector must be of class 'numeric' or 'Date'.")
    }
    
    if (length(y) != length(t)) {
      stop (paste0("Vectors 'y' and 't' differ in length: ", length(y), ", ", length(t)))
    }
    
    yt_df <- data.frame(y, t)
    yt_df_length_clean <- yt_df[complete.cases(yt_df), ]
    
    if (nrow(yt_df_length_clean) < nrow(yt_df)) {
      
      y <- yt_df$y
      t <- yt_df$t
      
      nan_cases <- nrow(yt_df)-nrow(yt_df_length_clean)
      warning(paste0(nan_cases, " NaN values were removed"))
      
    }
    
    if (any(y <= 0)) {
      
      if (!is.null(add_constant) && is.numeric(add_constant)) {
        
        y_lin <- y+add_constant
        
        message(paste0("Input infections data contains values <= 0. All values are increased by a constant equal to ", add_constant, "."))
        
      } else {
        
        y_lin <- y[y > 0]
        
        warning(paste0("Input infections data contains values <= 0. These values are skipped."))
      }
      
    } else {
      
      y_lin <- y
      
    }
    
    TP <- length(y_lin)
    
    class_t <- "numeric"
    
    if (is.Date(t)) {
      
      class_t <- "Date"
      
      message("Time vector is of class 'Date'. Calculating time counter.")
      
      start_date <- min(t)
      
      time_counter <- as.integer(t-start_date)
      
      t <- time_counter
    }
    
    bp_formula <- as.formula("y ~ t")
    if (isTRUE(ln)) {
      bp_formula <- as.formula("log(y) ~ t")
    }
    
    if (isTRUE(verbose)) {
      if (isTRUE(ln)) { 
        cat(paste0("Calculating breakpoints for time series with ", TP, " time points, nat. log. of y values ... "))
      } else {
        cat(paste0("Calculating breakpoints for time series with ", TP, " time points ... "))
        
      }
    }
    
    bp_obj <- 
      breakpoints(
        formula = bp_formula, 
        data = yt_df_length_clean, 
        ...
      )
    
    bp_obj_coef <- coef(bp_obj)
    bp_obj_breakpoints <- bp_obj$breakpoints
    bp_obj_breakpoints_no <- length(bp_obj_breakpoints)
    bp_obj_segments_no <- bp_obj_breakpoints_no+1
    
    if (bp_obj_breakpoints_no > 0) {
      bp_obj_breakpoints_ci <- confint(
        bp_obj, 
        level = 1-alpha
      )  
    }
    
    if(isTRUE(verbose)) {
      cat("OK", "\n")
      cat(paste0(bp_obj_breakpoints_no, " breakpoints were found, which leads to ", bp_obj_segments_no, " segments"), "\n")
      cat(paste0("Calculating segment models ... "))
    }
    
    
    bp_obj_segments <- 
      matrix(
        ncol = 2, 
        nrow = bp_obj_segments_no
      )
    
    bp_obj_segments[1,] <- c(1, bp_obj_breakpoints[1])
    
    i <- 0
    
    for (i in 2:(bp_obj_segments_no-1)) {
      bp_obj_segments[i,] <- 
        c(
          (bp_obj_breakpoints[i-1]+1), 
          bp_obj_breakpoints[i]
        )
    }
    
    bp_obj_segments[bp_obj_segments_no,] <- 
      c(
        bp_obj_breakpoints[bp_obj_breakpoints_no], 
        length(y)
      )
    
    
    i <- 0
    
    slopes <- matrix(ncol = 3, nrow = bp_obj_segments_no)
    slopes_p <- vector()
    intercepts <- matrix(ncol = 3, nrow = bp_obj_segments_no)
    intercepts_p <- vector()
    SQR <- vector()
    SAR <- vector()
    SQT <- vector()
    R2 <- vector()
    MSE <- vector()
    RMSE <- vector()
    MAE <- vector()
    MAPE <- vector()
    segment_model_ci_points <- list()
    
    for (i in 1:bp_obj_segments_no) {
      
      yt_df_length_clean_segment <- yt_df_length_clean[bp_obj_segments[i,1]:bp_obj_segments[i,2],]
      
      segment_model <- 
        lm(
          bp_formula, 
          data = yt_df_length_clean_segment
        )
      
      slopes[i, 1] <- segment_model$coefficients[2]
      intercepts[i, 1] <- segment_model$coefficients[1]
      
      segment_model_ci_est <- 
        confint(
          segment_model,
          level = 1-alpha
        )
      slopes[i, 2] <- segment_model_ci_est[2]
      slopes[i, 3] <- segment_model_ci_est[4]
      intercepts[i, 2] <- segment_model_ci_est[1]
      intercepts[i, 3] <- segment_model_ci_est[3]
      
      segment_model_summary <- summary(segment_model)
      
      intercepts_p[i] <- segment_model_summary$coefficients[7]
      slopes_p[i] <- segment_model_summary$coefficients[8]
      
      segment_model_y <- yt_df_length_clean_segment$y
      segment_model_y_pred <- predict(segment_model)
      if (isTRUE(ln)) {
        segment_model_y_pred <- exp(segment_model_y_pred)
      }
      
      segment_model_fit_metrics <- metrics(
        segment_model_y,
        segment_model_y_pred,
        plot = FALSE
      )
      segment_model_fit_metrics <- segment_model_fit_metrics[[1]]
      SQR[i] <- segment_model_fit_metrics$SQR
      SAR[i] <- segment_model_fit_metrics$SAR
      SQT[i] <- segment_model_fit_metrics$SQT
      R2[i] <- segment_model_fit_metrics$R2
      MSE[i] <- segment_model_fit_metrics$MSE
      RMSE[i] <- segment_model_fit_metrics$RMSE
      MAE[i] <- segment_model_fit_metrics$MAE
      MAPE[i] <- segment_model_fit_metrics$MAPE
      segment_model_fit_metrics <- NULL
      
      segment_model_ci_points[[i]] <- 
        predict(
          segment_model, 
          interval = "confidence", 
          level = 1-alpha
        )
      
    }
    
    segments_models_df <- 
      cbind(
        bp_obj_segments, 
        intercepts, 
        intercepts_p, 
        slopes, 
        slopes_p,
        SQR,
        SAR,
        SQT,
        R2,
        MSE,
        RMSE,
        MAE,
        MAPE
      )
    
    colnames(segments_models_df) <- 
      c(
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
    
    model_breakpoints_list <- 
      list (
        segments_models_df = segments_models_df,
        segment_model_ci_points = segment_model_ci_points,
        bp_obj_coef = bp_obj_coef,
        bp_obj_breakpoints = bp_obj_breakpoints,
        bp_obj_breakpoints_no = bp_obj_breakpoints_no,
        bp_obj_segments_no = bp_obj_segments_no,
        bp_obj_breakpoints_ci = bp_obj_breakpoints_ci,
        bp_obj = bp_obj
      )
    
    config <-
      list (
        TP = TP,
        ln = ln,
        add_constant = add_constant,
        alpha = alpha
      )
    
    if(isTRUE(verbose)) {
      cat("OK", "\n") 
    }
    
    new(
      "breaksgrowth",
      GrowthModel_OLS = model_breakpoints_list, 
      y = y,
      t = t,
      config = config
    )
    
  }


R_t <- 
  function (
    infections, 
    GP = 4,
    correction = FALSE
  ) {
    
    if(!is.numeric(infections)) {
      stop("Vector 'infections' must be of type 'numeric'")
    }
    
    i <- 0
    
    infections_daily_A <- vector()
    infections_daily_B <- vector()
    
    R_t <- vector()
    R_t[1:(GP-1)] <- NA
    
    for (i in (GP*2):length(infections)) {
      
      infections_daily_A[i] <- sum (infections[(i-(GP-1)):i])
      infections_daily_B[i] <- sum (infections[(i-((GP*2)-1)):(i-GP)]) 
      
      if (correction == TRUE) {
        
        if (infections_daily_B[i] < 1) {
          
          infections_daily_B[i] <- 1
          
        }
      }
      
      R_t[i] <- as.numeric(infections_daily_A[i]/infections_daily_B[i])
      
    }
    
    results <- list (
      R_t = R_t,
      infections_data = cbind.data.frame(infections_daily_A, infections_daily_B, R_t)
    )
    
    return(results)
    
  }