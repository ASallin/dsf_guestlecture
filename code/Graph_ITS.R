# Import csv
library(readr)
library(dplyr)
library(ggplot2)
library(qs)

dataGraph <- qread("data/dataGraph.qs") 

dfgraph <- dataGraph$dfGraph_long
counterfactualM <- dataGraph$counterfactual_long
out <- dataGraph$out
x <- 1

# Create graph
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
