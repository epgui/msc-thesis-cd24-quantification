library("ggplot2")
library("mvoutlier")
library("HH")
library("plyr")
library("dplyr") # Always load after plyr to avoid problems!
library("ggpubr")
library("ggsci")
library("reshape2")
library("nlme")
library("lme4")
library("knitr")

setwd("/Users/gp/Documents/Maitrise/Results/N3L2/786 n1/")
data <- read.csv("Statistics.csv") # Relative proportions data frame

order.of.cell.lines <- c('786-O', 'VHL')

# Constants for bead quantification procedure
bead.volume.used       <-       5        # in uL
bead.concentration     <- 1000000 / 1000 # beads / uL
total.sample.volume    <-    1000        # uL
used.sample.volume     <-     100        # uL
normalized.to.x.beads  <-    5000


data$Temps   <- as.factor(data$Temps)
data$A23187  <- as.factor(data$A23187)
data$Lignée   = factor(data$Lignée, levels = order.of.cell.lines)

data$Norm.P.events     = (data$P.Events   / data$P.Beads) * normalized.to.x.beads
data$Norm.P.positive   = (data$P.Positive / data$P.Beads) * normalized.to.x.beads
data$Norm.P.high       = (data$P.High / data$P.Beads) * normalized.to.x.beads

data$P.RelativeConcentration = 1 / data$P.Dilution

pd <- position_dodge(0.1) # move the error bars .05 to the left and right


f <- function(param) { return(var(c((4251 - param) * 5, (3901 - param) * 10, (3807 - param) * 20, (3511 - param) * 50))) }
o <- optimize(f, interval=c(1000,5000))
o$minimum

# Nombre de microparticules?

ggplot(data, aes(x = data$P.RelativeConcentration, y = (Norm.P.high), fill = Lignée)) +
  labs(
    title    = "Étalonnage des microparticules",
    subtitle = "PMA 30 nM (n = 1)",
    x        = "Concentration",
    y        = "Microparticules") +
  theme(
    panel.background   = element_rect(fill = "#F4F4F4"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_line(colour = "#EEEEEE"),
    axis.ticks         = element_blank(),
    axis.title.x       = element_text(vjust = -2),
    axis.title.y       = element_text(angle = 90, vjust = 2)) +
  scale_y_continuous(
    limits       = c(0, 10000),
    expand       = c(0, 0),
    labels       = scales::scientific,
    breaks       = seq(0, 10000, 1000),
    minor_breaks = seq(0, 10000, 200)) +
  geom_bar(
    position = position_dodge(),
    stat     = "identity") +
  facet_grid(
    . ~ Lignée,
    labeller = as_labeller(
      c(`786-O` = "786-O (ATCC CRL-1932)",
        `VHL`   = "786-O/VHL")))
