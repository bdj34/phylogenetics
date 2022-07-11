# Plot lambdas from jason with lambdas from phylofit without ACF
#NOTE: Need 9478 from jason
library(ggplot2)

rm(list = ls())

df_jason <- read.csv("~/rotation_fall2021/lambdas_from_jason_05292022.csv")
df_phylofit <- read.csv("~/rotation_fall2021/williams_clade_estimates_threshold_n_10_phylofit_no_ACF07062022.csv")

df_phylofit$driver_and_pid_and_age <- paste0(df_phylofit$pid, " ", df_phylofit$driver, " ", round(df_phylofit$age_at_sample, 2))

which(df_jason$X %in% df_phylofit$driver_and_pid_and_age)

df_combined <- df_jason
df_combined$phylofit_growthRate <- 0
df_combined$phylofit_lb_growthRate <- 0
df_combined$phylofit_ub_growthRate <- 0

for (i in 1:length(df_jason$X)) {
  
  df_phylofit$driver_and_pid_and_age %in% df_jason$X
  
  # Get the row from phylofit df and copy to combined
  row_phylofit <- which(df_phylofit$driver_and_pid_and_age == df_jason$X[i])
  
  df_combined$phylofit_growthRate[i] <- df_phylofit$phylofit_growthRate[row_phylofit]
  df_combined$phylofit_lb_growthRate[i] <- df_phylofit$phylofit_lb_growthRate[row_phylofit]
  df_combined$phylofit_ub_growthRate[i] <- df_phylofit$phylofit_ub_growthRate[row_phylofit]
}

df_combined$X <- factor(df_combined$X,  levels = unique(df_combined$X[order(df_combined$phylofit_growthRate)]))

pdf("~/rotation_fall2021/growthRates_phylofit_vs_ours_07062022.pdf")
ggplot(df_combined, aes(x = X)) + geom_point(aes(y = phylofit_growthRate), color = "blue", size = 2.5) + 
  geom_errorbar(aes(ymin=phylofit_lb_growthRate, ymax=phylofit_ub_growthRate), width=.3, color = "blue") + 
  geom_point(aes(y = our_lambda), color = "red", size = 2, alpha = .8) +
  geom_errorbar(aes(ymin = our_lambda_lower, ymax=our_lambda_upper), width = .2, color = "red")+
  #geom_point(aes(y = lambda_moment), color = "red", fill = "black", shape = 24, size = 2, alpha = .8) + 
  #geom_errorbar(aes(ymin = lambda_moment_lower, ymax=lambda_moment_upper), width = .2, color = "red")+
  #geom_point(aes(y = lambda_mle), color = "red", fill = "black", shape = 23, size = 2, alpha = .8) +
  #geom_errorbar(aes(ymin=our_growthRate_VCF_mut_Tree_minusSD, ymax=our_growthRate_VCF_mut_Tree_plusSD), width=.2, color = "red")+
  geom_text(aes(label = paste0("n=", N_tips), y = our_lambda_upper + 0.1), size = 4)+
  ylab("Net growth rates (per year)") + #ylim(0, 2) 
  theme_bw() +
  ggtitle("Method comparisons") +
  theme(axis.text.x = element_text(face = "bold", angle = 50, vjust = 1, hjust=1, size = 10),
        panel.grid.minor = element_blank(), axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(face = "bold", size = 10))
dev.off()

pdf("~/rotation_fall2021/growthRates_ABC_with_ACF_vs_ours_07062022.pdf")
ggplot(df_combined, aes(x = X)) + geom_point(aes(y = lambda_ABC), color = "blue", size = 2.5) + 
  geom_errorbar(aes(ymin=lower_lambda_ABC, ymax=upper_lambda_ABC), width=.3, color = "blue") + 
  geom_point(aes(y = our_lambda), color = "red", size = 2, alpha = .8) +
  geom_errorbar(aes(ymin = our_lambda_lower, ymax=our_lambda_upper), width = .2, color = "red")+
  #geom_point(aes(y = lambda_moment), color = "red", fill = "black", shape = 24, size = 2, alpha = .8) + 
  #geom_errorbar(aes(ymin = lambda_moment_lower, ymax=lambda_moment_upper), width = .2, color = "red")+
  #geom_point(aes(y = lambda_mle), color = "red", fill = "black", shape = 23, size = 2, alpha = .8) +
  #geom_errorbar(aes(ymin=our_growthRate_VCF_mut_Tree_minusSD, ymax=our_growthRate_VCF_mut_Tree_plusSD), width=.2, color = "red")+
  geom_text(aes(label = paste0("n=", N_tips), y = our_lambda_upper + 0.1), size = 4)+
  ylab("Net growth rates (per year)") + #ylim(0, 2) 
  theme_bw() +
  ggtitle("Method comparisons") +
  theme(axis.text.x = element_text(face = "bold", angle = 50, vjust = 1, hjust=1, size = 10),
        panel.grid.minor = element_blank(), axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(face = "bold", size = 10))
dev.off()  


df_cut <- df_combined[df_combined$N_tips >= 20,]

pdf("~/rotation_fall2021/growthRates_threshold_20_phylofit_vs_ours_07062022.pdf")
ggplot(df_cut, aes(x = X)) + geom_point(aes(y = phylofit_growthRate), color = "blue", size = 2.5) + 
  geom_errorbar(aes(ymin=phylofit_lb_growthRate, ymax=phylofit_ub_growthRate), width=.3, color = "blue") + 
  geom_point(aes(y = our_lambda), color = "red", size = 2, alpha = .8) +
  geom_errorbar(aes(ymin = our_lambda_lower, ymax=our_lambda_upper), width = .2, color = "red")+
  geom_text(aes(label = paste0("n=", N_tips), y = our_lambda_upper + 0.1), size = 4)+
  ylab("Net growth rates (per year)") + #ylim(0, 2) 
  theme_bw() +
  ggtitle("Method comparisons") +
  theme(axis.text.x = element_text(face = "bold", angle = 50, vjust = 1, hjust=1, size = 10),
        panel.grid.minor = element_blank(), axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(face = "bold", size = 10))
dev.off()

pdf("~/rotation_fall2021/growthRates_threshold_20_ABC_with_ACF_vs_ours_07062022.pdf")
ggplot(df_cut, aes(x = X)) + geom_point(aes(y = lambda_ABC), color = "blue", size = 2.5) + 
  geom_errorbar(aes(ymin=lower_lambda_ABC, ymax=upper_lambda_ABC), width=.3, color = "blue") + 
  geom_point(aes(y = our_lambda), color = "red", size = 2, alpha = .8) +
  geom_errorbar(aes(ymin = our_lambda_lower, ymax=our_lambda_upper), width = .2, color = "red")+
  geom_text(aes(label = paste0("n=", N_tips), y = our_lambda_upper + 0.1), size = 4)+
  ylab("Net growth rates (per year)") + #ylim(0, 2) 
  theme_bw() +
  ggtitle("Method comparisons") +
  theme(axis.text.x = element_text(face = "bold", angle = 50, vjust = 1, hjust=1, size = 10),
        panel.grid.minor = element_blank(), axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(face = "bold", size = 10))
dev.off() 
