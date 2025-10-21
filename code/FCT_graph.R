bootstrapITS <- function(fit, B, data, nmonths = 12, formula, group = NULL){
  
  BS <- list()
  length(BS) <- B
  
  for (b in 1:B) {
    # Sample 95% of the IDs
    sampled_ids <- sample(unique(data$ID), size = ceiling(length(unique(data$ID)) * 0.95))
    
    # Subset the data based on sampled IDs
    dataBS <- data[data$ID %in% sampled_ids, ]
    
    updatefit <- feols(
      as.formula(formula),
      data = dataBS,
      warn = FALSE,
      notes = FALSE
    )
    
    BS[[b]] <- resultsITS(data = dataBS, updatefit, nmonths, group = group)
    cat(b, "...")
  }
  BS <- bind_rows(BS)
  return(BS)
}


resultsITS <- function(fit, data, nmonths = 12, group = NULL){
  
  data <- data  |>
    mutate(
      pred = predict(fit, data),
      counterfactual = predict(fit, mutate(data, post = 0, post.2))
    )
  
  absolutechange <- (sum(data[data$time == 1,]$pred - data[data$time == 1,]$counterfactual)) /
    nrow(data[data$time == 1,])
  
  absolutechange12 <- (sum(data[data$time %in% 1:nmonths,]$pred -
                             data[data$time %in% 1:nmonths,]$counterfactual)) /
    nrow(data[data$time == 1,])
  
  relativechange = sum(data[data$time == 1,]$pred - data[data$time == 1,]$counterfactual) /
    sum(data[data$time == 1,]$counterfactual)
  
  relativechange12 <- sum(data[data$time %in% 1:nmonths,]$pred -
                            data[data$time %in% 1:nmonths,]$counterfactual) /
    sum(data[data$time %in% 1:nmonths,]$counterfactual)
  
  if(!is.null(group)){
    # Differences are computed for each group. Then the differences of the differences are taken and used for
    # testing whether they are statistically significant.
    dataG1 <- data[data[[group]] == 1, ]
    dataG0 <- data[data[[group]] == 0, ]
    
    absolutechangeG0 <- (sum(dataG0[dataG0$time == 1, ]$pred - dataG0[dataG0$time == 1,]$counterfactual)) /
      nrow(dataG0[dataG0$time == 1,])
    absolutechangeG1 <- (sum(dataG1[dataG1$time == 1, ]$pred - dataG1[dataG1$time == 1,]$counterfactual)) /
      nrow(dataG1[dataG1$time == 1,])
    
    absolutechange12_G0 <- (sum(dataG0[dataG0$time %in% 1:nmonths,]$pred -
                                  dataG0[dataG0$time %in% 1:nmonths,]$counterfactual)) /
      nrow(dataG0[dataG0$time == 1,])
    absolutechange12_G1 <- (sum(dataG1[dataG1$time %in% 1:nmonths,]$pred -
                                  dataG1[dataG1$time %in% 1:nmonths,]$counterfactual)) /
      nrow(dataG1[dataG1$time == 1,])
    
    relativechangeG0 = sum(dataG0[dataG0$time == 1,]$pred -
                             dataG0[dataG0$time == 1,]$counterfactual) /
      sum(dataG0[dataG0$time == 1,]$counterfactual)
    relativechangeG1 = sum(dataG1[dataG1$time == 1,]$pred -
                             dataG1[dataG1$time == 1,]$counterfactual) /
      sum(dataG1[dataG1$time == 1,]$counterfactual)
    
    relativechange12_G0 = sum(dataG0[dataG0$time %in% 1:nmonths,]$pred -
                                dataG0[dataG0$time %in% 1:nmonths,]$counterfactual) /
      sum(dataG0[dataG0$time %in% 1:nmonths,]$counterfactual)
    relativechange12_G1 = sum(dataG1[dataG1$time %in% 1:nmonths,]$pred -
                                dataG1[dataG1$time %in% 1:nmonths,]$counterfactual) /
      sum(dataG1[dataG1$time %in% 1:nmonths,]$counterfactual)
  }
  
  if (is.null(group)) {
    result_list <- list(
      "absolutechange" = absolutechange,
      "absolutechange12" = absolutechange12,
      "relativechange" = relativechange,
      "relativechange12" = relativechange12
    )
  } else {
    result_list <- list(
      "absolutechange_G0" = absolutechangeG0,
      "absolutechange_G1" = absolutechangeG1,
      "absolutechange12_G0" = absolutechange12_G0,
      "absolutechange12_G1" = absolutechange12_G1,
      "absolutechange_diff"   = absolutechangeG0 - absolutechangeG1,
      "absolutechange12_diff" = absolutechange12_G0 - absolutechange12_G1,
      
      "relativechange_G0" = relativechangeG0,
      "relativechange_G1" = relativechangeG1,
      "relativechange12_G0" = relativechange12_G0,
      "relativechange12_G1" = relativechange12_G1,
      "relativechange_group"   = relativechangeG0 - relativechangeG1,
      "relativechange12_group" = relativechange12_G0 - relativechange12_G1
    )
    
    
    
    
  }
  return(result_list)
}


bootstrapITS_ngroups <- function(fit, B, data, nmonths = 12, formula, group = NULL){
  
  BS <- list()
  length(BS) <- B
  
  for (b in 1:B) {
    # Sample 95% of the IDs
    sampled_ids <- sample(unique(data$ID), size = ceiling(length(unique(data$ID)) * 0.95))
    
    # Subset the data based on sampled IDs
    dataBS <- data[data$ID %in% sampled_ids, ]
    
    updatefit <- feols(
      as.formula(formula),
      data = dataBS,
      warn = FALSE,
      notes = FALSE
    )
    
    BS[[b]] <- resultsITS_ngroups(data = dataBS, updatefit, nmonths, group = group)
    cat(b, "...")
  }
  
  BS <- list(
    "effect" = bind_rows(lapply(BS, function(x) x[[1]])),
    "differences" = bind_rows(lapply(BS, function(x) x[[2]]))
  )
  return(BS)
}


# Bootstrap for more than two groups
resultsITS_ngroups <- function(fit, data, nmonths = 12, group = NULL) {
  
  data <- data |>
    mutate(
      pred = predict(fit, data),
      counterfactual = predict(fit, mutate(data, post = 0, post.2))
    )
  
  stopifnot(!is.null(group))
  
  # Check levels of the group variable
  group_levels <- unique(data[[group]])
  group_levels <- as.character(group_levels)
  
  # Initialize list to store results per group
  group_results <- lapply(group_levels, function(g) {
    
    # Filter on group
    group_data <- data[data[[group]] == g, ]
    
    # Compute effects
    data.frame(
      "group" = as.character(g),
      
      "absolutechange" = sum(group_data[group_data$time == 1,]$pred - group_data[group_data$time == 1,]$counterfactual) /
        nrow(group_data[group_data$time == 1,]),
      
      "absolutechange12" = sum(group_data[group_data$time %in% 1:nmonths,]$pred -
                                 group_data[group_data$time %in% 1:nmonths,]$counterfactual) /
        nrow(group_data[group_data$time == 1,]),
      
      "relativechange" = sum(group_data[group_data$time == 1,]$pred - group_data[group_data$time == 1,]$counterfactual) /
        sum(group_data[group_data$time == 1,]$counterfactual),
      
      "relativechange12" = sum(group_data[group_data$time %in% 1:nmonths,]$pred -
                                 group_data[group_data$time %in% 1:nmonths,]$counterfactual) /
        sum(group_data[group_data$time %in% 1:nmonths,]$counterfactual)
    )
  })
  
  # Compute pairwise differences between groups
  # combn takes m elements from the x object and applies a function to these two elements
  pairwise_differences <- combn(x = group_levels, m = 2, function(pair) {
    # Ensure consistent group order
    sorted_pair <- sort(pair)
    g1 <- sorted_pair[1]
    g2 <- sorted_pair[2]
    
    group1 <- group_results[[which(group_levels == g1)]]
    group2 <- group_results[[which(group_levels == g2)]]
    
    data.frame(
      "groups" = paste(g1, g2, sep = "_vs_"),
      "absolutechange_diff" = group1$absolutechange - group2$absolutechange,
      "absolutechange12_diff" = group1$absolutechange12 - group2$absolutechange12,
      "relativechange_diff" = group1$relativechange - group2$relativechange,
      "relativechange12_diff" = group1$relativechange12 - group2$relativechange12
    )
  }, simplify = FALSE)
  
  result_list <- list(
    per_group = bind_rows(group_results),
    pairwise_differences = bind_rows(pairwise_differences)
  )
  
  return(result_list)
}





# Plot ITS

plotITS <- function(model, y, x, add_vline = NULL, ylabel, counterfactual = TRUE) {
  
  d <- model$model
  colnames(d)[colnames(d) == y] = "y"
  colnames(d)[colnames(d) ==  "as.factor(MONAT_NR_STA)"] = "MONAT_NR_STA"
  colnames(d)[colnames(d) ==  "as.factor(JAHR_STA)"] = "JAHR_STA"
  
  cov  <- colnames(d)[!colnames(d) %in% c("y")]
  cov[cov == "as.factor(MONAT_NR_STA)"]  <- "MONAT_NR_STA"
  cov[cov == "as.factor(JAHR_STA)"]      <- "JAHR_STA"
  
  postIntervention <- as.numeric(d[d[[x]] == min(d[[x]]), x][1])
  
  d <- data.table(d)
  
  dfgraph <- d[, lapply(.SD, mean),
               by = c(cov),
               .SDcols = "y"
  ]
  dfgraph[x == 0, "post"] <- 1
  
  # Intervention time
  dfgraph <- dfgraph %>%
    mutate(group_intervention = case_when(dfgraph[[x]] <= -1 ~ 1,
                                          dfgraph[[x]] >= 0 & dfgraph[[x]] <= 14 ~ 2,
                                          .default = 3))
  
  dfgraph$predicted <- predict(model, dfgraph)
  dfgraph <- arrange(as_tibble(dfgraph),dfgraph[[x]])
  
  
  # Plot overall structure
  structurePlot <- ggplot(dfgraph, aes(y = y, x = .data[[x]])) +
    geom_line(aes(color = "Observed"), linewidth = 1) +
    labs(
      x = "Time to/from the interventions in months",
      y = ylabel
    ) +
    theme_classic() +
    scale_x_continuous(
      breaks = c(-50, -39, -27, -15, -3, 0, 9, 15, 21),
      labels = c(
        "Jan.\n2017", "Jan.\n2018", "Jan.\n2019",
        "Jan.\n2020", "Jan.\n2021", "Apr.\n2021",
        "Jan.\n2022", "Jun.\n2022", "Jan.\n2023"
      )
    ) +
    annotate("rect", xmin = -13, xmax = -9, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "#1C77D2") +
    geom_vline(xintercept = 0, color = "#E22920", lty = 20, linewidth = 1.5) +
    geom_vline(xintercept = 15, color = "#675482", lty = 20, linewidth = 1.5) +
    theme(legend.background = element_rect(fill = "#EEEEEE"),
          legend.position = c(0.1, 0.9),
          panel.grid.major.x = element_line(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12)
    )
  
  structurePlot <- structurePlot +
    geom_line(aes(x = time, y = as.matrix(dfgraph[, "y"]), color = "Observed", linetype = "Observed", size = "Observed")) +
    geom_line(aes(
      x = time,
      y = as.matrix(dfgraph[, "predicted"]),
      color = "Predicted",
      linetype = "Predicted",
      size = "Predicted"
    ))
  
  if (counterfactual) {
    # Counterfactual predicted values
    dfgraph$counterfactual <- dfgraph %>%
      mutate(post = 0, post.2 = 0) %>%
      predict(model, .)
    
    counterfactualM <- as.matrix(dfgraph[, "counterfactual"])
    counterfactualM[1:-postIntervention, ] <- NA
    
    # Plot
    structurePlot <- structurePlot +
      geom_line(aes(x = time, y = counterfactualM, color = "Counterfactual", linetype = "Counterfactual", size = "Counterfactual")) +
      scale_linetype_manual(
        name = element_blank(),
        values = c(
          "Observed" = 1,
          "Predicted" = 1,
          "Counterfactual" = 4
        ),
        labels = c("Counterfactual", "Observed", "Predicted")
      ) +
      scale_size_manual(
        name = element_blank(),
        values = c(
          "Observed" = 0.8,
          "Predicted" = 1,
          "Counterfactual" = 1.25
        ),
        labels = c("Counterfactual", "Observed", "Predicted")
      ) +
      scale_color_manual(
        name = element_blank(),
        values = c(
          "Observed" = "black",
          "Predicted" = "#E22920",
          "Counterfactual" = "#01BBA8"
        ),
        labels = c(
          "Counterfactual", "Observed", "Predicted"
        )
      ) +
      theme(
        legend.background = element_rect(fill = "white"),
        legend.position = c(0.9, 0.96)
      )
  } else {
    structurePlot <- structurePlot +
      scale_linetype_manual(
        name = element_blank(),
        values = c(
          "Observed" = 1,
          "Predicted" = 1
        ),
        labels = c("Observed", "Predicted")
      ) +
      scale_size_manual(
        name = element_blank(),
        values = c(
          "Observed" = 0.8,
          "Predicted" = 1
        ),
        labels = c("Observed", "Predicted")
      ) +
      scale_color_manual(
        name = element_blank(),
        values = c(
          "Observed" = "black",
          "Predicted" = "#E22920"
        ),
        labels = c("Observed", "Predicted")
      ) +
      theme(
        legend.background = element_rect(fill = "white"),
        legend.position = c(0.9, 0.96)
      )
  }
  
  if (!is.null(add_vline)) {
    for (i in add_vline) {
      structurePlot <-
        structurePlot +
        geom_vline(
          xintercept = i, color = SWICAcols[8], lty = 70,
          linewidth = 0.6
        )
    }
  }
  return(structurePlot)
}




# Plot overall structure
graphITSforSMR <- function(dfgraph, counterfactualM, out) {
  
  dfgraph <- if (out == "y_") {
    dfgraph <- dfgraph |>
      rename("y" = y_) |>
      mutate(y = y*100,
             pred = pred*100,
             counterfactualM = counterfactualM*100)
  } else if (out == "y1_") {
    dfgraph <- dfgraph |>
      rename("y" = y1_) |>
      mutate(y = y*100,
             counterfactualM = counterfactualM*100)
  } else if (out == "ycost_") {
    rename(dfgraph, "y" = ycost_)
  } else if (out == "y1cost_") {
    rename(dfgraph, "y" = y1cost_)
  }
  
  
  g <- ggplot(dfgraph, aes(y = y, x = time)) +
    geom_line(aes(color = "Observed", size = "Observed", linetype = "Observed",
                  alpha = ifelse(time <= 1, "before", "after"))
    ) +
    labs(
      x = "Time in months"
    ) +
    geom_vline(xintercept = 0, color = "#E22920", lty = 20, linewidth = 1) +
    geom_vline(xintercept = 14.5, color = "#E22920", lty = 20, linewidth = 1) +
    theme_classic() +
    geom_line(
      data = subset(dfgraph, time <= 0),
      aes(
        x = time,
        y = as.matrix(dfgraph[dfgraph$time <=0, "pred"]),
        col = "Predicted",
        size = "Predicted",
        linetype = "Predicted",
        alpha = "before"
      )
    ) +
    geom_line(
      data = subset(dfgraph, time >= 0),
      aes(
        x = time,
        y = as.matrix(dfgraph[dfgraph$time >= 0, "pred"]),
        col = "Predicted",
        size = "Predicted",
        linetype = "Predicted",
        alpha = "after"
      )
    ) +
    geom_line(
      aes(
        x = time,
        y = counterfactualM,
        col = "Counterfactual",
        size = "Counterfactual",
        linetype = "Counterfactual",
        alpha = ifelse(time <= 0, "before", "after")
      )
    ) +
    scale_linetype_manual(
      name = element_blank(),
      values = c(
        "Observed" = 1,
        "Predicted" = 1,
        "Counterfactual" = 4
      ),
      labels = c("Counterfactual", "Observed", "Predicted")
    ) +
    scale_size_manual(
      name = element_blank(),
      values = c("Observed" = 0.8,
                 "Predicted" = 1,
                 "Counterfactual" = 1.25),
      labels = c("Counterfactual", "Observed", "Predicted")
    ) +
    scale_color_manual(
      name = element_blank(),
      values = c(
        "Observed" = "grey40",
        "Predicted" = "#E22920",
        "Counterfactual" = "#01BBA8"
      ),
      labels = c(
        "Counterfactual", "Observed", "Predicted"
      )
    ) +
    scale_alpha_manual(
      name = element_blank(),
      values = c("before" = 0.5,
                 "after" = 1)
    ) +
    scale_x_continuous(
      breaks = c(-50, -39, -27, -15, -3, 0, 9, 15, 21),
      labels = c(
        "Jan. 2017", "Jan. 2018", "Jan. 2019",
        "Jan. 2020", "Jan. 2021", "  Apr. 2021",
        "Jan. 2022", "Jun. 2022", "Jan. 2023"
      ),
      limits = c(-39, 14)
    ) +
    guides(alpha = "none") +
    annotate("rect", xmin = -13, xmax = -9, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "#00857C") +
    geom_text(aes(x = -11, y = 1.8, label = "COVID19 \n Lockdown"), 
              color = "#00857C", size = 3) +
    theme(legend.background = element_rect(fill = "white"),
          legend.position = c(0.1, 0.95),
          panel.grid.major.x = element_line(),
          panel.grid.major.y = element_line())
  
  g <-  if (out == "y_") {
    g <- g +
      labs(
        y = "Average number of monthly VitaminD test \n per physician"
      )
  } else if (out == "y1_") {
    g <- g +
      labs(
        y = "Average number of monthly VitaminD test \n per physician per 100 consultations"
      )
  } else if (out == "ycost_") {
    g <- g +
      labs(
        y = "Average costs of monthly VitaminD test \n per physician"
      )
  } else if (out == "y1cost_") {
    g <- g +
      labs(
        y = "Average costs of monthly VitaminD test \n per physician per 100 consultations"
      )
  }
  
  g
}


# Graph for LIM
graphITSforLIM <- function(dfgraph, counterfactualM, out) {
  
  dfgraph <- if (out == "y_") {
    dfgraph <- dfgraph |>
      rename("y" = y_) |>
      mutate(y = y * 100,
             pred = pred * 100,
             counterfactualM = counterfactualM * 100)
  } else if (out == "y1_") {
    dfgraph <- dfgraph |>
      rename("y" = y1_) |>
      mutate(y = y * 100,
             pred = pred*100,
             counterfactualM = counterfactualM * 100)
  } else if (out == "ycost_") {
    rename(dfgraph, "y" = ycost_)
  } else if (out == "y1cost_") {
    rename(dfgraph, "y" = y1cost_)
  }
  
  g <- ggplot(dfgraph, aes(y = y, x = time)) +
    geom_line(aes(color = "Observed", size = "Observed", linetype = "Observed",
                  alpha = ifelse(time <= 0, "before", "after"))
    ) +
    labs(
      x = "Time to/from the interventions in months"
    ) +
    geom_vline(xintercept = 0, color = SWICAcols[10], lty = 20, linewidth = 1) +
    geom_vline(xintercept = -15, color = SWICAcols[10], lty = 20, linewidth = 1) +
    theme_classic() +
    geom_line(
      data = subset(dfgraph, time <= 0),
      aes(
        x = time,
        y = as.matrix(dfgraph[dfgraph$time <= 0, "pred"]),
        col = "Predicted",
        size = "Predicted",
        linetype = "Predicted",
        alpha = "before"
      )
    ) +
    geom_line(
      data = subset(dfgraph, time >= 0),
      aes(
        x = time,
        y = as.matrix(dfgraph[dfgraph$time >= 0, "pred"]),
        col = "Predicted",
        size = "Predicted",
        linetype = "Predicted",
        alpha = "after"
      )
    ) +
    geom_line(
      aes(
        x = time,
        y = counterfactualM,
        col = "Counterfactual",
        size = "Counterfactual",
        linetype = "Counterfactual",
        alpha = ifelse(time <= 0, "before", "after")
      )
    ) +
    scale_linetype_manual(
      name = element_blank(),
      values = c(
        "Observed" = 1,
        "Predicted" = 1,
        "Counterfactual" = 4
      ),
      labels = c("Counterfactual", "Observed", "Predicted")
    ) +
    scale_size_manual(
      name = element_blank(),
      values = c("Observed" = 0.8,
                 "Predicted" = 1,
                 "Counterfactual" = 1.25),
      labels = c("Counterfactual", "Observed", "Predicted")
    ) +
    scale_color_manual(
      name = element_blank(),
      values = c(
        "Observed" = SWICAcols[7],
        "Predicted" = SWICAcols[10],
        "Counterfactual" = SWICAcols[2]
      ),
      labels = c(
        "Counterfactual", "Observed", "Predicted"
      )
    ) +
    scale_alpha_manual(
      name = element_blank(),
      values = c("before" = 0.5,
                 "after" = 1)
    ) +
    scale_x_continuous(
      breaks = c(-50, -39, -27, -15, -3, 0, 9, 15, 21) - 15,
      labels = c(
        "Jan. 2017", "Jan. 2018", "Jan. 2019",
        "Jan. 2020", "Jan. 2021", "Apr. 2021",
        "Jan. 2022", "Jun. 2022", "Jan. 2023"
      ),
      limits = c(-50, 14)
    ) +
    guides(alpha = "none") +
    annotate("rect", xmin = -28, xmax = -24, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = SWICAcols[3]) +
    geom_text(aes(x = -26, y = 1.6, label = "COVID19"), color = SWICAcols[3], size = 3, alpha = 0.1) +
    theme(legend.background = element_rect(fill = "white"),
          legend.position = c(0.1, 0.95),
          panel.grid.major.x = element_line(),
          panel.grid.major.y = element_line())
  
  g <-  if (out == "y_") {
    g <- g +
      labs(
        y = "Average number of monthly VitaminD tests \n per physician"
      )
  } else if (out == "y1_") {
    g <- g +
      labs(
        y = "Average number of monthly VitaminD tests \n per physician per 100 consultations"
      )
  } else if (out == "ycost_") {
    g <- g +
      labs(
        y = "Average costs of monthly VitaminD tests \n per physician"
      )
  } else if (out == "y1cost_") {
    g <- g +
      labs(
        y = "Average costs of monthly VitaminD tests \n per physician per 100 consultations"
      )
  }
  
  g
}


graphITS_long_poster <- function(dfgraph, counterfactualM, out) {
  dfgraph <- if (out == "y_") {
    rename(dfgraph, "y" = y_)
  } else {
    rename(dfgraph, "y" = y1_)
  }
  
  # counterfactualM[40:54, ] <- NA
  
  ggplot(dfgraph, aes(y = y*100, x = time)) +
    geom_line(aes(color = "Observed",
                  size = "Observed",
                  linetype = "Observed")) +
    labs(
      x = "Time to/from the interventions in months",
      y = "Average number of monthly VitaminD perscriptions \n per physician per 100 consultations"
    ) +
    theme_classic() +
    geom_line(
      data = subset(dfgraph, time <= 0),
      aes(
        x = time,
        y = as.matrix(dfgraph[dfgraph$time <=0, "pred"])*100,
        col = "Predicted",
        size = "Predicted",
        linetype = "Predicted"
      )
    ) +
    geom_line(
      data = subset(dfgraph, time >= 0),
      aes(
        x = time,
        y = as.matrix(dfgraph[dfgraph$time >= 0, "pred"])*100,
        col = "Predicted",
        size = "Predicted",
        linetype = "Predicted"
      )
    ) +
    geom_line(
      aes(
        x = time,
        y = counterfactualM*100,
        col = "Counterfactual",
        size = "Counterfactual",
        linetype = "Counterfactual"
      ),
    ) +
    scale_linetype_manual(
      name = element_blank(),
      values = c(
        "Observed" = 1,
        "Predicted" = 1,
        "Counterfactual" = 4
      ),
      labels = c("Counterfactual", "Observed", "Predicted")
    ) +
    scale_size_manual(
      name = element_blank(),
      values = c("Observed" = x*0.8,
                 "Predicted" = x*1,
                 "Counterfactual" = x*1.25),
      labels = c("Counterfactual", "Observed", "Predicted")
    ) +
    scale_color_manual(
      name = element_blank(),
      values = c(
        "Observed" = SWICAcols[7],
        "Predicted" = SWICAcols[10],
        "Counterfactual" = SWICAcols[2]
      ),
      labels = c(
        "Counterfactual", "Observed", "Predicted"
      )
    ) +
    scale_x_continuous(
      breaks = c(-50, -39, -27, -15, -3, 0, 9, 15, 21)-15,
      labels = c(
        "Jan.\n2017", "Jan.\n2018", "Jan.\n2019",
        "Jan.\n2020", "Jan.\n2021", "Apr.\n2021",
        "Jan.\n2022", "Jun.\n2022", "Jan.\n2023"
      )
    ) +
    annotate("rect", xmin = -28, xmax = -24, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = SWICAcols[3]) +
    geom_text(aes(x = -26, y = 1.6, label = "Lockdown"), color = SWICAcols[3], size = x*3) +
    geom_vline(xintercept = -15, color = SWICAcols[10], lty = 20, linewidth = x*1.5) +
    geom_vline(xintercept = 0, color = "#675482", lty = 20, linewidth = x*1.5) +
    theme(legend.background = element_rect(fill = "#EEEEEE"),
          legend.position = c(0.1, 0.9),
          panel.grid.major.x = element_line(),
          panel.grid.major.y = element_line(),
          plot.background = element_rect(fill = "#EEEEEE"),
          # panel.background = element_rect(fill = "gray95"),
          # plot.background = element_rect(fill = "white"),
          plot.margin = margin(4, 10, 4, 10),  # Remove margins around the plot to avoid white space
          axis.text.x = element_text(size = x*12),
          axis.text.y = element_text(size = x*12),
          axis.title = element_text(size = x*12),
          legend.text = element_text(size = x*12)
    ) +
    coord_fixed(
      ratio = 4.5,
      xlim = c(-55, 11),
      ylim = c(1.5, 6.9),
      expand = FALSE)
}



graphITS_long_paper <- function(dfgraph, counterfactualM, out) {
  if (out == "y_") {
    dfgraph <- rename(dfgraph, "y" = y_)
  } else if (out == "y1_") {
    dfgraph <- rename(dfgraph, "y" = y1_)
  } else if (out == "SUPP_Anteil") {
    dfgraph <- rename(dfgraph, "y" = SUPP_Anteil)
  }
  
  ggplot(dfgraph, aes(y = y*100, x = time)) +
    geom_line(aes(color = "Observed",
                  size = "Observed",
                  linetype = "Observed")) +
    labs(
      x = "",
      y = "Average number of monthly VitaminD perscriptions \n per physician per 100 consultations"
    ) +
    theme_classic() +
    geom_line(
      data = subset(dfgraph, time <= 0),
      aes(
        x = time,
        y = as.matrix(dfgraph[dfgraph$time <=0, "pred"])*100,
        col = "Predicted",
        size = "Predicted",
        linetype = "Predicted"
      )
    ) +
    geom_line(
      data = subset(dfgraph, time >= 0),
      aes(
        x = time,
        y = as.matrix(dfgraph[dfgraph$time >= 0, "pred"])*100,
        col = "Predicted",
        size = "Predicted",
        linetype = "Predicted"
      )
    ) +
    geom_line(
      aes(
        x = time,
        y = counterfactualM*100,
        col = "Counterfactual",
        size = "Counterfactual",
        linetype = "Counterfactual"
      ),
    ) +
    scale_linetype_manual(
      name = element_blank(),
      values = c(
        "Observed" = 1,
        "Predicted" = 1,
        "Counterfactual" = 4
      ),
      labels = c("Counterfactual", "Observed", "Predicted")
    ) +
    scale_size_manual(
      name = element_blank(),
      values = c("Observed" = x*0.8,
                 "Predicted" = x*1,
                 "Counterfactual" = x*1.25),
      labels = c("Counterfactual", "Observed", "Predicted")
    ) +
    scale_color_manual(
      name = element_blank(),
      values = c(
        "Observed" = "grey40",
        "Predicted" = "#E22920",
        "Counterfactual" = "#01BBA8"
      ),
      labels = c(
        "Counterfactual", "Observed", "Predicted"
      )
    ) +
    scale_x_continuous(
      breaks = c(-50, -39, -27, -15, -3, 0, 9, 15, 21)-15,
      labels = c(
        "Jan.\n2017", "Jan.\n2018", "Jan.\n2019",
        "Jan.\n2020", "Jan.\n2021", "Apr.\n2021",
        "Jan.\n2022", "Jun.\n2022", "Jan.\n2023"
      )
    ) +
    annotate("rect", xmin = -28, xmax = -24, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "#00857C") +
    geom_text(aes(x = -26, y = 1.8, label = "COVID19 \nLockdown"), color = "#00857C", size = x*3) +
    geom_vline(xintercept = -15, color = "#E22920", lty = 20, linewidth = x*1.5) +
    geom_vline(xintercept = 0, color = "#675482", lty = 20, linewidth = x*1.5) +
    theme(legend.background = element_rect(fill = "#EEEEEE"),
          legend.position = c(0.1, 0.9),
          panel.grid.major.x = element_line(),
          panel.grid.major.y = element_line(),
          # plot.background = element_rect(fill = "#EEEEEE"),
          # panel.background = element_rect(fill = "gray95"),
          # plot.background = element_rect(fill = "white"),
          # plot.margin = margin(4, 10, 4, 10),  # Remove margins around the plot to avoid white space
          axis.text.x = element_text(size = x*12),
          axis.text.y = element_text(size = x*12),
          axis.title = element_text(size = x*12),
          legend.text = element_text(size = x*12)
    ) 
}

