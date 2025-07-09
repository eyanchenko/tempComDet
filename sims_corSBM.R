########## FUNCTIONS ##############
library(ggplot2)
library(dplyr)
library(RColorBrewer)

source("functions.R")

############## COLORS
myColors <- brewer.pal(5,"Set1")


############## SIMULATIONS
set.seed(12345)
####### Setting 1. Increasing delta (rho=0.25)
n.sims = 100

n = 1000
b = 0.01
delta.seq = seq(0, 0.010, 0.0005)

K = 2
C <- c(rep(1,0.80*n), rep(2, 0.20*n))

rho = 0.25
TT = 10

cnt = 1
df <- tibble(iter = rep(rep(1:n.sims, each=5), length(delta.seq)), 
             delta  = rep(delta.seq, each=n.sims*5),
             method = rep(c("int", "VS", "kappa=0.25", "kappa=0.50", "kappa=0.75"), 
                          n.sims*length(delta.seq)), 
             e=0)

for(delta in delta.seq){
  B = matrix(c(b+delta/2, b-delta/2, b-delta/2, b+delta/2), nrow=2, ncol=2)
  
  for(iter in 1:n.sims){
    # Generate networks
    A <- generate_SBM(n, B, C, rho, TT)
    
    pval <- spectral.adj.pval(A)
    
    df$e[cnt + 0] <- mean(p_to_e(pval, "int"))
    df$e[cnt + 1] <- mean(p_to_e(pval, "VS"))
    df$e[cnt + 2] <- mean(p_to_e(pval, "kappa", kappa=0.25))
    df$e[cnt + 3] <- mean(p_to_e(pval, "kappa", kappa=0.50))
    df$e[cnt + 4] <- mean(p_to_e(pval, "kappa", kappa=0.75))
    
    cnt = cnt + 5
  }
  save(df, file="~/Documents/Research/CommDetTemp/Flip_PtoE/Results/sims_corSBM_vary_delta_rho0.25.RData")
  print(delta)
}


load("~/Documents/Research/CommDetTemp/Flip_PtoE/Results/sims_corSBM_vary_delta_rho0.25.RData")
df$e[is.na(df$e)] <- Inf
df_plot = df %>% group_by(method, delta) %>% summarize(e=median(e))
df_plot$method <- factor(df_plot$method, levels = c("kappa=0.25", "kappa=0.50", "kappa=0.75",
                                                    "int", "VS"))

p1 <- ggplot(df_plot, aes(x=delta, y=e, color=method)) +
  geom_line() +
  scale_colour_manual(
    values = myColors,
    labels = c("int" = "Int", "VS" = "VS", "kappa=0.25" = expression(kappa==0.25),
               "kappa=0.50" = expression(kappa==0.50), "kappa=0.75" = expression(kappa==0.75))
  ) +
  guides(
    color = guide_legend(title = "Calibrator"),
    linetype = guide_legend(title = "Calibrator"),
  ) +
  xlab(expression(delta)) +
  ylab("e-value") +
  scale_y_continuous(trans = 'log10') +
  theme_bw()+
  theme(legend.position="bottom")
p1

ggsave("~/Documents/Research/CommDetTemp/Flip_PtoE/Results/sims_corSBM_vary_delta_rho0.25.pdf", 
       units="in", width = 5, height = 4)



# Variability plots showing each observation
ggplot(df, aes(x=delta, y=e, color=method))+
  geom_point()+
  scale_colour_manual(
    values = myColors,
    labels = c("int" = "Int", "VS" = "VS", "kappa=0.25" = expression(kappa==0.25),
               "kappa=0.50" = expression(kappa==0.50), "kappa=0.75" = expression(kappa==0.75))
  ) +
  guides(
    color = guide_legend(title = "Calibrator"),
    linetype = guide_legend(title = "Calibrator"),
  ) +
  xlab(expression(delta)) +
  ylab("e-value") +
  scale_y_continuous(limits = c(0.1, 100000), trans = 'log10') +  
  theme_bw()+
  theme(legend.position="bottom")

ggsave("~/Documents/Research/CommDetTemp/Flip_PtoE/Results/sims_corSBM_vary_delta_rho0.25_var.pdf", 
       units="in", width = 5, height = 4)




####### Setting 2. Increasing delta (rho=0.75)
n.sims = 100

n = 1000
b = 0.01
delta.seq = seq(0, 0.010, 0.0005)

K = 2
C <- c(rep(1,0.80*n), rep(2, 0.20*n))

rho = 0.75
TT = 10

cnt = 1
df <- tibble(iter = rep(rep(1:n.sims, each=5), length(delta.seq)), 
             delta  = rep(delta.seq, each=n.sims*5),
             method = rep(c("int", "VS", "kappa=0.25", "kappa=0.50", "kappa=0.75"), 
                          n.sims*length(delta.seq)), 
             e=0)

for(delta in delta.seq){
  B = matrix(c(b+delta/2, b-delta/2, b-delta/2, b+delta/2), nrow=2, ncol=2)
  
  for(iter in 1:n.sims){
    # Generate networks
    A <- generate_SBM(n, B, C, rho, TT)
    
    pval <- spectral.adj.pval(A)
    
    df$e[cnt + 0] <- mean(p_to_e(pval, "int"))
    df$e[cnt + 1] <- mean(p_to_e(pval, "VS"))
    df$e[cnt + 2] <- mean(p_to_e(pval, "kappa", kappa=0.25))
    df$e[cnt + 3] <- mean(p_to_e(pval, "kappa", kappa=0.50))
    df$e[cnt + 4] <- mean(p_to_e(pval, "kappa", kappa=0.75))
    
    cnt = cnt + 5
  }
  save(df, file="~/Documents/Research/CommDetTemp/Flip_PtoE/Results/sims_corSBM_vary_delta_rho0.75.RData")
  print(delta)
}


load("~/Documents/Research/CommDetTemp/Flip_PtoE/Results/sims_corSBM_vary_delta_rho0.75.RData")
df$e[is.na(df$e)] <- Inf

df_plot = df %>% group_by(method, delta) %>% summarize(e=median(e))
df_plot$method <- factor(df_plot$method, levels = c("kappa=0.25", "kappa=0.50", "kappa=0.75",
                                                    "int", "VS"))

p2 <- ggplot(df_plot, aes(x=delta, y=e, color=method)) +
  geom_line() +
  scale_colour_manual(
    values = myColors,
    labels = c("int" = "Int", "VS" = "VS", "kappa=0.25" = expression(kappa==0.25),
               "kappa=0.50" = expression(kappa==0.50), "kappa=0.75" = expression(kappa==0.75))
  ) +
  guides(
    color = guide_legend(title = "Calibrator"),
    linetype = guide_legend(title = "Calibrator")
  ) +
  xlab(expression(delta)) +
  ylab("e-value") +
  scale_y_continuous(trans = 'log10') +
  theme_bw()+
  theme(legend.position="bottom")
p2

ggsave("~/Documents/Research/CommDetTemp/Flip_PtoE/Results/sims_corSBM_vary_delta_rho0.75.pdf", 
       units="in", width = 5, height = 4)








####### Setting 3. Increasing n (rho=0.25)
n.sims = 100

n.seq = seq(100, 1000, 50)
b = 0.01
delta = 0.009
B = matrix(c(b+delta/2, b-delta/2, b-delta/2, b+delta/2), nrow=2, ncol=2)

K = 2
rho = 0.25
TT = 10

cnt = 1
df <- tibble(iter = rep(rep(1:n.sims, each=5), length(n.seq)), 
             n  = rep(n.seq, each=n.sims*5),
             method = rep(c("int", "VS", "kappa=0.25", "kappa=0.50", "kappa=0.75"), 
                          n.sims*length(n.seq)), 
             e=0)

for(n in n.seq){
  C <- c(rep(1,0.80*n), rep(2, 0.20*n))
  
  for(iter in 1:n.sims){
    # Generate networks
    A <- generate_SBM(n, B, C, rho, TT)
    
    pval <- spectral.adj.pval(A)
    
    df$e[cnt + 0] <- mean(p_to_e(pval, "int"))
    df$e[cnt + 1] <- mean(p_to_e(pval, "VS"))
    df$e[cnt + 2] <- mean(p_to_e(pval, "kappa", kappa=0.25))
    df$e[cnt + 3] <- mean(p_to_e(pval, "kappa", kappa=0.50))
    df$e[cnt + 4] <- mean(p_to_e(pval, "kappa", kappa=0.75))
    
    cnt = cnt + 5
  }
  save(df, file="~/Documents/Research/CommDetTemp/Flip_PtoE/Results/sims_corSBM_vary_n_rho0.25.RData")
  print(n)
}


load("~/Documents/Research/CommDetTemp/Flip_PtoE/Results/sims_corSBM_vary_n_rho0.25.RData")
df$e[is.na(df$e)] <- Inf

df_plot = df %>% group_by(method, n) %>% summarize(e=median(e))
df_plot$method <- factor(df_plot$method, levels = c("kappa=0.25", "kappa=0.50", "kappa=0.75",
                                                    "int", "VS"))

p3 <- ggplot(df_plot, aes(x=n, y=e, color=method)) +
  geom_line() +
  scale_colour_manual(
    values = myColors,
    labels = c("int" = "Int", "VS" = "VS", "kappa=0.25" = expression(kappa==0.25),
               "kappa=0.50" = expression(kappa==0.50), "kappa=0.75" = expression(kappa==0.75))
  ) +
  guides(
    color = guide_legend(title = "Calibrator"),
    linetype = guide_legend(title = "Calibrator")
  ) +
  xlab("n") +
  ylab("e-value") +
  scale_y_continuous(trans = 'log10') +
  theme_bw()+
  theme(legend.position="bottom")
p3

ggsave("~/Documents/Research/CommDetTemp/Flip_PtoE/Results/sims_corSBM_vary_n_rho0.25.pdf", 
       units="in", width = 5, height = 4)


# Variability plots showing each observation
ggplot(df, aes(x=n, y=e, color=method))+
  geom_point()+
  scale_colour_manual(
    values = myColors,
    labels = c("int" = "Int", "VS" = "VS", "kappa=0.25" = expression(kappa==0.25),
               "kappa=0.50" = expression(kappa==0.50), "kappa=0.75" = expression(kappa==0.75))
  ) +
  guides(
    color = guide_legend(title = "Calibrator"),
    linetype = guide_legend(title = "Calibrator"),
  ) +
  xlab("n") +
  ylab("e-value") +
  scale_y_continuous(limits = c(0.1, 100000), trans = 'log10') +
  theme_bw()+
  theme(legend.position="bottom")

ggsave("~/Documents/Research/CommDetTemp/Flip_PtoE/Results/sims_corSBM_vary_n_rho0.25_var.pdf", 
       units="in", width = 5, height = 4)



####### Setting 4. Increasing n (rho=0.25)
n.sims = 100

n.seq = seq(100, 1000, 50)
b = 0.01
delta = 0.009
B = matrix(c(b+delta/2, b-delta/2, b-delta/2, b+delta/2), nrow=2, ncol=2)

K = 2

rho = 0.75
TT = 10

cnt = 1
df <- tibble(iter = rep(rep(1:n.sims, each=5), length(n.seq)), 
             n  = rep(n.seq, each=n.sims*5),
             method = rep(c("int", "VS", "kappa=0.25", "kappa=0.50", "kappa=0.75"), 
                          n.sims*length(n.seq)), 
             e=0)

for(n in n.seq){
  C <- c(rep(1,0.80*n), rep(2, 0.20*n))
  
  for(iter in 1:n.sims){
    # Generate networks
    A <- generate_SBM(n, B, C, rho, TT)
    
    pval <- spectral.adj.pval(A)
    
    df$e[cnt + 0] <- mean(p_to_e(pval, "int"))
    df$e[cnt + 1] <- mean(p_to_e(pval, "VS"))
    df$e[cnt + 2] <- mean(p_to_e(pval, "kappa", kappa=0.25))
    df$e[cnt + 3] <- mean(p_to_e(pval, "kappa", kappa=0.50))
    df$e[cnt + 4] <- mean(p_to_e(pval, "kappa", kappa=0.75))
    
    cnt = cnt + 5
  }
  save(df, file="~/Documents/Research/CommDetTemp/Flip_PtoE/Results/sims_corSBM_vary_n_rho0.75.RData")
  print(n)
}


load("~/Documents/Research/CommDetTemp/Flip_PtoE/Results/sims_corSBM_vary_n_rho0.75.RData")
df$e[is.na(df$e)] <- Inf

df_plot = df %>% group_by(method, n) %>% summarize(e=median(e))
df_plot$method <- factor(df_plot$method, levels = c("kappa=0.25", "kappa=0.50", "kappa=0.75",
                                                    "int", "VS"))

p4 <- ggplot(df_plot, aes(x=n, y=e, color=method)) +
  geom_line() +
  scale_colour_manual(
    values = myColors,
    labels = c("int" = "Int", "VS" = "VS", "kappa=0.25" = expression(kappa==0.25),
               "kappa=0.50" = expression(kappa==0.50), "kappa=0.75" = expression(kappa==0.75))
  ) +
  guides(
    color = guide_legend(title = "Calibrator"),
    linetype = guide_legend(title = "Calibrator")
  ) +
  xlab("n") +
  ylab("e-value") +
  scale_y_continuous(trans = 'log10') +
  theme_bw()+
  theme(legend.position="bottom")
p4

ggsave("~/Documents/Research/CommDetTemp/Flip_PtoE/Results/sims_corSBM_vary_n_rho0.75.pdf", 
       units="in", width = 5, height = 4)


