# Call simpop functions
library(phangorn)
library(ape)
library(ggplot2)
library(ggtree)
library(rsimpop)
library(rtreefit)
library(gridExtra)

rm(list = ls())

source("~/rotation_fall2021/ultra_tree_lambda_fn.R")
#source("~/rotation_fall2021/init_growth_rate_with_and_without_zeroes.R")
source("~/rotation_fall2021/init_growth_rate_vec.R")

SEED=37774323
initSimPop(SEED,bForce = TRUE)

################# Post Conception ####################################
#birthRate_init <- 0.05 # per day, immediately following conception (vary this)
#df <- init_growth_rate_with_and_without_zeroes(birthRate_init, nIters = 10)

birthRate_initVec <- c(0.005, 0.01, 0.05, 0.1, 0.5, 1, 10) # per day (units irrelevant, birth only)
df <- init_growth_rate_vec(birthRate_initVec, nIters = 1000, nYears = 20)

ggplot(df) + geom_point(aes(x = method, y = log(lambda, base = 10), color = factor(iteration))) + 
  geom_hline(yintercept = c(log(birthRate_initVec, base = 10)), linetype = "dashed") +
  #geom_line(aes(x = method, y = lambda, group = factor(iteration)), linetype = "dashed", size = 0.1)+
  theme_bw()

plot_list <- list()
for (i in c(1:length(birthRate_initVec))) {
  p <- list(ggplot(df[df$actual_lambda == unique(df$actual_lambda)[i],], aes( x = method, y = lambda)) + 
    geom_violin() + geom_jitter(size = .1, width = .3, height = 0) + 
    theme_bw() + 
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1),
          plot.title = element_text(size = 10)) +
    geom_hline(yintercept = unique(df$actual_lambda)[i], linetype = "dashed") + 
    ylim(c(0, 2*unique(df$actual_lambda)[i])) + 
    ggtitle(paste0("Actual lambda: ", unique(df$actual_lambda)[i])))
  plot_list <- append(plot_list, p)
}

grid.arrange(grobs = plot_list, ncol = length(birthRate_initVec))
write.csv(df, "~/rotation_fall2021/simulated_distributions_06052022.csv")
