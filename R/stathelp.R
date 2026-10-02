#---------------------------------------------------------------
# Name:        stathelp (swash package)
# Purpose:     Statistical helper functions for the swash package
# Author:      Thomas Wieland 
#              ORCID: 0000-0001-5168-9846
#              mail: geowieland@googlemail.com
# Version:     1.0.0
# Last update: 2026-09-12 10:13
# Copyright (c) 2022-2026 Thomas Wieland
#---------------------------------------------------------------



quantile_ci <-
  function(
    x,
    alpha = 0.05
  ) {
    
    ci_lower <- alpha/2
    ci_upper <- 1-(alpha/2)
    
    ci <- 
      quantile(
        x, 
        probs = c(ci_lower, ci_upper)
      )
    
    return(ci)
    
  }


hist_ci <-
  function(
    x,
    alpha = 0.05,
    col_bars = "grey",
    col_ci = "red",
    ...
  ) {
    
    ci <- quantile_ci(
      x = x, 
      alpha = alpha
    )
    
    hist(x, col = col_bars, ...)
    abline(v = ci[1], col = col_ci)
    abline(v = ci[2], col = col_ci)
    abline(v = median(x), col = col_ci)
    
  }



plot_coef_ci <- function(
    point_estimates,
    confint_lower,
    confint_upper,
    coef_names,
    p = NULL,
    estimate_colors = NULL,
    confint_colors = NULL,
    auto_color = FALSE,
    alpha = 0.05,
    set_estimate_colors = c("red", "grey", "green"),
    set_confint_colors = c("#ffcccb", "lightgray", "#CCFFCC"),
    skipvars = NULL,
    plot.xlab = "Independent variables",
    plot.main = "Point estimates with CI",
    axis.at = seq(-30, 40, by = 5),
    pch = 15,
    cex = 2,
    lwd = 5,
    y.cex = 0.8
) {
  
  if (
    length(coef_names) == length(point_estimates) &&
    length(point_estimates) == length(confint_lower) &&
    length(confint_lower) == length(confint_upper)
  ) {
    
    data_plot <- 
      data.frame(
        coef_names, 
        point_estimates,
        confint_lower,
        confint_upper
      )
    
  } else {
    
    stop("Vectors coef_names, point_estimates, confint_lower, and confint_upper differ in length.")
    
  }
  
  if (!is.null(p)) {
    
    if (length(p) == nrow(data_plot)) {
      
      data_plot$p <- p
      
    } else {
      
      stop("Vector p differs in length from coef_names, point_estimates, confint_lower, and confint_upper")
      
    }
    
  }
  
  if (!is.null(estimate_colors) & !is.null(confint_colors)) {
    
    data_plot$col_point <- estimate_colors
    data_plot$col_confint <- confint_colors
    
  } else {
    
    if (is.null(p)) {
      
      if (auto_color == FALSE) {
        
        stop("If estimate_colors and confint_colors are NULL and auto_color is FALSE, then p must be stated as numeric vector")
        
      } else {
        
        data_plot$col_point <- set_estimate_colors[2]
        data_plot$col_confint <- set_confint_colors[2]
        
        data_plot[
          (data_plot$point_estimates < 0)
          & (data_plot$confint_lower < 0)
          & (data_plot$confint_upper < 0)
          ,]$col_point <- set_estimate_colors[1]
        data_plot[
          (data_plot$point_estimates < 0)
          & (data_plot$confint_lower < 0)
          & (data_plot$confint_upper < 0)
          ,]$col_confint <- set_confint_colors[1]
        
        data_plot[
          (data_plot$point_estimates > 0)
          & (data_plot$confint_lower > 0)
          & (data_plot$confint_upper > 0)
          ,]$col_point <- set_estimate_colors[3]
        data_plot[
          (data_plot$point_estimates > 0)
          & (data_plot$confint_lower > 0)
          & (data_plot$confint_upper > 0)
          ,]$col_confint <- set_confint_colors[3]
        
      }
      
    }
    
  }
  
  if (!is.null(skipvars)) {
    data_plot <- data_plot[!data_plot$coef_names %in% skipvars,]
  }
  
  if (!is.null(p)) {
    
    data_plot$col_point <- set_estimate_colors[2]
    data_plot$col_confint <- set_confint_colors[2]
    
    if (nrow(data_plot[(data_plot$point_estimates < 0) & (data_plot$p < alpha),]) > 0) {
      data_plot[(data_plot$point_estimates < 0) & (data_plot$p < alpha),]$col_point <- set_estimate_colors[1]
    }
    if (nrow(data_plot[(data_plot$point_estimates > 0) & (data_plot$p < alpha),]) > 0) {
      data_plot[(data_plot$point_estimates > 0) & (data_plot$p < alpha),]$col_point <- set_estimate_colors[3]
    }
    
    if (nrow(data_plot[data_plot$col_point == set_estimate_colors[1],]) > 0) {
      data_plot[data_plot$col_point == set_estimate_colors[1],]$col_confint <- set_confint_colors[1]
    }
    if (nrow(data_plot[data_plot$col_point == set_estimate_colors[3],]) > 0) {
      data_plot[data_plot$col_point == set_estimate_colors[3],]$col_confint <- set_confint_colors[3]
    }
    
  } 
  
  data_plot$order <- 1:nrow(data_plot)
  data_plot <- data_plot[order(-data_plot$order),]
  
  par(mar=c(4,15,2,1))
  
  plot(
    x = data_plot$point_estimates, 
    y = (1:nrow(data_plot)), 
    xlim = c(
      -max(abs(data_plot$confint_lower)), 
      max(abs(data_plot$confint_upper))
    ), 
    yaxt = "n", 
    ylab = "", 
    xaxt = "n", 
    cex = 0.1, 
    xlab = plot.xlab, 
    main = plot.main
  )
  
  axis(
    1, 
    at = axis.at, 
    tck = 1, 
    lty = 2, 
    col = "gray"
  )
  
  par(las=1)
  
  axis (
    side = 2, 
    at = 1:nrow(data_plot), 
    labels = data_plot$coef_names, 
    cex.axis = y.cex, 
    tick = FALSE
  )
  
  abline(h = 0)
  abline(v = 0)
  
  i <- 0
  
  for (i in 1:nrow(data_plot)) {
    
    lines (
      x = c(
        data_plot$confint_lower[i], 
        data_plot$confint_upper[i]
      ), 
      y = c(i,i), 
      lwd = lwd,
      col = data_plot$col_confint[i]
    )
    
  }
  
  points (
    x = data_plot$point_estimates, 
    y = 1:nrow(data_plot), 
    pch = pch, 
    cex = cex, 
    col = data_plot$col_point
  )
  
}


metrics <- 
  function(
    observed,
    expected,
    plot = TRUE,
    plot.main = "Observed vs. expected",
    xlab = "Observed",
    ylab = "Expected",
    point.col = "blue",
    point.pch = 19,
    line.col = "red",
    plot_residuals.main = "Residuals",
    legend.cex = 0.7
  ) {
    
    if (length(observed) != length(expected)) {
      stop("Vectors 'observed' and 'expected' differ in length")
    }
    
    observed_expected <- data.frame(observed, expected)
    
    observed_expected$residuals <- observed_expected$expected-observed_expected$observed
    observed_expected$residuals_abs <- abs(observed_expected$expected-observed_expected$observed)
    observed_expected$residuals_sq <- (observed_expected$expected-observed_expected$observed)^2
    observed_expected$residuals_rel <- observed_expected$residuals/observed_expected$observed*100
    observed_expected$residuals_rel_abs <- abs(observed_expected$residuals_rel)
    
    SQR <- sum(observed_expected$residuals_sq)
    SAR <- sum(observed_expected$residuals_abs)
    SQT <- sum((observed_expected$observed - mean(observed_expected$observed))^2)
    R2 <- (1-(SQR/SQT))
    MSE <- mean((observed_expected$observed - observed_expected$expected)^2)
    RMSE <- sqrt(MSE)
    MAE <- sum(observed_expected$observed-observed_expected$expected)
    MAPE <- mean(observed_expected$residuals_rel_abs)
    
    if (plot == TRUE) {
      
      min = min(observed_expected$observed)
      max = max(observed_expected$observed)
      
      plot(
        observed_expected$observed, 
        observed_expected$expected, 
        xlab = xlab,
        ylab = ylab,
        pch = point.pch, 
        col = point.col,
        main = plot.main,
        xlim = c((min-min*0.1), (max*1.1)),
        ylim = c((min-min*0.1), (max*1.1))
      )
      
      legend(
        "topleft", 
        legend = c(
          bquote(R^2 == .(round(R2, 2))),
          paste0("MSE = ", round(MSE, 2)), 
          paste0("RMSE = ", round(RMSE, 2)), 
          paste0("MAE = ", round(MAE, 2)),
          paste0("MAPE = ", round(MAPE, 2))
        ), 
        cex = legend.cex
      )
      
      abline(coef = c(0,1), col = line.col)
      
      Y_residuals_stat <- 
        data.frame(
          table(
            cut(
              observed_expected$residuals_rel, 
              breaks = c(-Inf, seq(-90, 90, 10), Inf)
            )
          )
        )
      colnames(Y_residuals_stat) <- c("Range", "Freq")
      
      Y_residuals_stat$Freq_rel <- round(Y_residuals_stat$Freq/sum(Y_residuals_stat$Freq)*100, 2)
      
      barplot (
        Y_residuals_stat$Freq_rel, 
        names.arg = Y_residuals_stat$Range, 
        ylim = c(0, max(Y_residuals_stat$Freq_rel)*1.5),
        main = plot_residuals.main
      )
      
      legend(
        "topleft", 
        legend = c(
          paste0("+/- 30 %: ", round(sum(Y_residuals_stat[8:13,]$Freq_rel), 1), " %"),
          paste0("+/- 20 %: ", round(sum(Y_residuals_stat[9:12,]$Freq_rel), 1), " %"),
          paste0("+/- 10 %: ", round(sum(Y_residuals_stat[10:11,]$Freq_rel), 1), " %")
        ),
        cex = legend.cex
      )
      
    }
    
    fit_metrics <- list(
      SQR = SQR,
      SAR = SAR,
      SQT = SQT,
      R2 = R2,
      MSE = MSE,
      RMSE = RMSE,
      MAE = MAE,
      MAPE = MAPE
    )
    
    invisible(
      list(
        fit_metrics = fit_metrics, 
        observed_expected = observed_expected
      )
    )
    
  }


binary_metrics <- function(
    observed, 
    expected,
    no_information_rate = "negative"
) {
  
  if (length(observed) != length(expected)) {
    stop("Vectors 'observed' and 'expected' differ in length")
  }
  
  if (!all(observed %in% c(0,1)) || !all(expected %in% c(0,1))) {
    stop("Observed and/or expected values are not binary")
  }
  
  observed_expected <- data.frame(observed, expected)
  observed_expected$hit <- 0
  observed_expected[observed_expected$observed == observed_expected$expected,]$hit <- 1
  
  tab <- table(observed, expected)
  
  expected_required <- c("0", "1")
  expected_NA <- setdiff(expected_required, colnames(tab))
  
  if (length(expected_NA) > 0) {
    
    warning(paste0("No expected category for ", expected_NA), "\n")
    
    for (col in expected_NA) {
      
      tab <- cbind(tab, rep(0, nrow(tab)))
      
      colnames(tab)[ncol(tab)] <- col
      
    }
    
    tab <- tab[, expected_required]
  }
  
  sens <- tab[2,2] / sum(tab[2,])
  spec <- tab[1,1] / sum(tab[1,])
  acc  <- sum(diag(tab)) / sum(tab)
  
  if (no_information_rate == "positive") {
    nir <- length(observed[observed == 1])/length(observed)
  } else {
    nir <- length(observed[observed == 0])/length(observed)
  }
  
  fit_metrics <- 
    list(
      sens = sens,
      spec = spec,
      acc = acc,
      nir = nir
    )
  
  invisible(
    list(
      fit_metrics = fit_metrics,
      observed_expected = observed_expected
    )
  )
  
}


binary_metrics_glm <- function(
    logit_model, 
    threshold = 0.5
) {
  
  y_pred <- 
    ifelse(
      predict(
        logit_model, 
        type = "response"
      ) > threshold, 1, 0
    )
  
  logit_metrics <- binary_metrics(
    observed = logit_model$y,
    expected = y_pred
  )
  
  invisible(logit_metrics)
  
}
