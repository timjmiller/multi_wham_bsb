library(wham)
library(ggplot2)
library(dplyr)
library(patchwork)
source(here::here("R","plot_functions.R"))

theme_set(theme_bw())
theme_update(strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = rel(2)), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)))

fit <- readRDS(file.path("results","fits_no_M_re_rev.RDS"))[[2]] #M_1
proj_1 <- readRDS(here::here("results",paste0("m1_proj_",1,".RDS")))
proj_1_fjp <- readRDS(here::here("results",paste0("m1_proj_",1,"_sdrep_fjp.RDS")))
proj <- readRDS(here::here("results","m1_proj_3_R_opt1.RDS"))
pred_R_N <- make_pred_R_cv(proj_1, proj_1_fjp, type = 1)
pred_R_S <- make_pred_R_cv(proj_1, proj_1_fjp, type = 1, region = 2)


df_use <- cbind.data.frame(Year = proj$input$years_Ecov, name = "Ecov_RE", val = TMB:::as.list.sdreport(proj$sdrep, report = TRUE, what = "Est")$Ecov_x[,1], 
	se = TMB:::as.list.sdreport(proj$sdrep, report = TRUE, what = "Std")$Ecov_x[,1], Region = "North")
df_use <- rbind(df_use, cbind.data.frame(Year = proj$input$years_Ecov, name = "Ecov_RE", val = TMB:::as.list.sdreport(proj$sdrep, report = TRUE, what = "Est")$Ecov_x[,2], 
	se = TMB:::as.list.sdreport(proj$sdrep, report = TRUE, what = "Std")$Ecov_x[,2], Region = "South"))
df_use <- rbind(df_use, cbind.data.frame(Year = proj$input$years_Ecov, name = "Ecov_obs", val = proj$input$data$Ecov_obs[,1], se = exp(proj$parList$Ecov_obs_logsigma[,1]), Region = "North"))
df_use <- rbind(df_use, cbind.data.frame(Year = proj$input$years_Ecov, name = "Ecov_obs", val = proj$input$data$Ecov_obs[,2], se = exp(proj$parList$Ecov_obs_logsigma[,2]), Region = "South"))

df_use <- rbind(df_use, cbind.data.frame(Year = proj$years_full, name = "SSB", val = TMB:::as.list.sdreport(proj$sdrep, report = TRUE, what = "Est")$log_SSB[,1], 
	se = TMB:::as.list.sdreport(proj$sdrep, report = TRUE, what = "Std")$log_SSB[,1], Region = "North"))
df_use <- rbind(df_use, cbind.data.frame(Year = proj$years_full, name = "SSB", val = TMB:::as.list.sdreport(proj$sdrep, report = TRUE, what = "Est")$log_SSB[,2], 
	se = TMB:::as.list.sdreport(proj$sdrep, report = TRUE, what = "Std")$log_SSB[,2], Region = "South"))
df_use <- rbind(df_use, cbind.data.frame(Year = proj$years_full, name = "SSB", val = TMB:::as.list.sdreport(proj$sdrep, report = TRUE, what = "Est")$log_SSB_all, 
	se = TMB:::as.list.sdreport(proj$sdrep, report = TRUE, what = "Std")$log_SSB_all, Region = "Total"))

df_use <- rbind(df_use, cbind.data.frame(Year = proj$years_full, name = "Fbar", val = TMB:::as.list.sdreport(proj$sdrep, report = TRUE, what = "Est")$log_Fbar[,5], 
	se = TMB:::as.list.sdreport(proj$sdrep, report = TRUE, what = "Std")$log_Fbar[,5], Region = "North"))
df_use <- rbind(df_use, cbind.data.frame(Year = proj$years_full, name = "Fbar", val = TMB:::as.list.sdreport(proj$sdrep, report = TRUE, what = "Est")$log_Fbar[,6], 
	se = TMB:::as.list.sdreport(proj$sdrep, report = TRUE, what = "Std")$log_Fbar[,6], Region = "South"))
df_use <- rbind(df_use, cbind.data.frame(Year = proj$years_full, name = "Fbar", val = TMB:::as.list.sdreport(proj$sdrep, report = TRUE, what = "Est")$log_Fbar[,7], 
	se = TMB:::as.list.sdreport(proj$sdrep, report = TRUE, what = "Std")$log_Fbar[,7], Region = "Total"))

df_use <- rbind(df_use, cbind.data.frame(Year = proj$years_full, name = "R_RE", val = TMB:::as.list.sdreport(proj$sdrep, report = TRUE, what = "Est")$log_NAA_rep[1,1,,1], 
	se = TMB:::as.list.sdreport(proj$sdrep, report = TRUE, what = "Std")$log_NAA_rep[1,1,,1], Region = "North"))
df_use <- rbind(df_use, cbind.data.frame(Year = proj$years_full, name = "R_RE", val = TMB:::as.list.sdreport(proj$sdrep, report = TRUE, what = "Est")$log_NAA_rep[2,2,,1], 
	se = TMB:::as.list.sdreport(proj$sdrep, report = TRUE, what = "Std")$log_NAA_rep[2,2,,1], Region = "South"))

df_use <- rbind(df_use, cbind.data.frame(Year = proj$years_full, name = "R_exp", val = log(pred_R_N[,1]), se = pred_R_N[,2], Region = "North"))
df_use <- rbind(df_use, cbind.data.frame(Year = proj$years_full, name = "R_exp", val = log(pred_R_S[,1]), se = pred_R_S[,2], Region = "South"))
df_use$lo <- df_use$val + qnorm(0.025)*df_use$se
df_use$hi <- df_use$val + qnorm(0.975)*df_use$se


this_df <- subset(df_use, Year < 2022 & name %in% c("Ecov_obs", "Ecov_RE"))
this_df$name <- c("Random Effect", "Observed")[match(this_df$name, c("Ecov_RE","Ecov_obs"))]
this_df$name <- factor(this_df$name)

p_E_full <- ggplot(this_df, aes(x = Year, y = val, colour = name, fill = name)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + 
    labs(color = "Covariate Type", fill = "Covariate Type") + xlab("Year") +
    facet_grid(~Region) + ylab("Bottom Temperature\nAnomaly") + geom_vline(xintercept = 1989, linetype = 2, alpha = 0.4, linewidth = 1.5) +
    geom_point(data = subset(this_df, name == "Observed"), aes(x = Year, y = val, colour = name), size = 2) + 
    geom_pointrange(data = subset(this_df, name == "Observed"), aes(ymin = lo, ymax = hi, colour = name)) + 
    geom_point(data = subset(this_df, name == "Random Effect"), aes(x = Year, y = val, colour = name), alpha = 0.6, shape = 22, size = 4) + 
    geom_line(data = subset(this_df, name == "Random Effect"), aes(x = Year, y = val, colour = name), linetype = 2) + 
    geom_ribbon(data =  subset(this_df, name == "Random Effect"), aes(ymin=lo, ymax=hi, fill = name), alpha=0.4, linetype = 0) 
p_E_full

###########################################################################################
# Fig. S4
cairo_pdf(file.path("paper", "BTA_full_fig.pdf"), width = 20, height = 10)
p_E_full
dev.off()
###########################################################################################

this_df <- subset(df_use, Year > 1988 & Year < 2022 & name %in% c("Ecov_obs", "Ecov_RE"))
this_df$name <- c("Random Effect", "Observed")[match(this_df$name, c("Ecov_RE","Ecov_obs"))]
this_df$name <- factor(this_df$name)

p_E <- ggplot(this_df, aes(x = Year, y = val, colour = name, fill = name)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + 
    labs(color = "Covariate Type", fill = "Covariate Type") + xlab("Year") +
    facet_grid(~Region) + ylab("Bottom Temperature\nAnomaly") + 
    geom_point(data = subset(this_df, name == "Observed"), aes(x = Year, y = val, colour = name), size = 2) + 
    geom_pointrange(data = subset(this_df, name == "Observed"), aes(ymin = lo, ymax = hi, colour = name)) + 
    geom_point(data = subset(this_df, name == "Random Effect"), aes(x = Year, y = val, colour = name), alpha = 0.6, shape = 22, size = 4) + 
    geom_line(data = subset(this_df, name == "Random Effect"), aes(x = Year, y = val, colour = name), linetype = 2, linewidth = 1) + 
    geom_ribbon(data =  subset(this_df, name == "Random Effect"), aes(ymin=lo, ymax=hi, fill = name), alpha=0.4, linetype = 0) 
p_E

this_df <- subset(df_use, Year < 2022 & Year > 1988 & name %in% c("R_RE", "R_exp"))
this_df$name <- c("Random Effect", "Expected")[match(this_df$name, c("R_RE","R_exp"))]
this_df$name <- factor(this_df$name)
this_df$val <- exp(this_df$val)/1e3
this_df$hi <- exp(this_df$hi)/1e3
this_df$lo <- exp(this_df$lo)/1e3

cnms <- c(North = "", South = "")
plt_R <- ggplot(this_df, aes(x = Year, y = val, colour = name)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") +  
    labs(colour = "Recruitment Type", fill = "Recruitment Type") +
    geom_point(size = 2) + geom_line(linetype = 2)  + 
    ylab(expression(paste("Recruitment (",10^6,")"))) + 
    # ylab(expression(paste("Recruitment (",phantom(0)%*%10^6,")"))) + 
    xlab("Year") +
    facet_grid(~ Region, labeller = labeller(Region = cnms)) + coord_cartesian(ylim = c(0,1.25e2)) +
    geom_ribbon(ggplot2::aes(ymin=lo, ymax=hi, fill=name), alpha=0.3, linetype = 0)
plt_R

this_df <- subset(df_use, Year < 2022 & Year > 1988 & name %in% c("SSB"))
this_df$val <- exp(this_df$val)/1e3
this_df$hi <- exp(this_df$hi)/1e3
this_df$lo <- exp(this_df$lo)/1e3

cnms <- c(North = "", South = "", Total = "Total")
plt_SSB <- ggplot(this_df, aes(x = Year, y = val)) + 
    # scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") +  
    # labs(color = "Recruitment Type", fill = "Recruitment Type") +
    geom_point(size = 2) + geom_line(linetype = 2)  + ylab("SSB (kmt)") + xlab("Year") +
    facet_grid(~ Region, labeller = labeller(Region = cnms)) + #coord_cartesian(ylim = c(0,6e4)) +#ylim(0,1e5) + #xlim(1989,2021)
    geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.3, linetype = 0)
plt_SSB


this_df <- subset(df_use, Year < 2022 & Year > 1988 & name %in% c("Fbar"))
this_df$val <- exp(this_df$val)
this_df$hi <- exp(this_df$hi)
this_df$lo <- exp(this_df$lo)

cnms <- c(North = "", South = "", Total = "")
plt_F <- ggplot(this_df, aes(x = Year, y = val)) + 
    # scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") +  
    # labs(color = "Recruitment Type", fill = "Recruitment Type") +
    geom_point(size = 2) + geom_line(linetype = 2)  + ylab(expression(bar(italic(F))~~"(Ages 6-7)")) + xlab("Year") +
    facet_grid(~ Region, labeller = labeller(Region = cnms)) + #coord_cartesian(ylim = c(0,6e4)) +#ylim(0,1e5) + #xlim(1989,2021)
    geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.3, linetype = 0)
plt_F


###########################################################################################
# Fig. 4
cairo_pdf(file.path("paper", "E_R_SSB_F_fig.pdf"), width = 20, height = 16)
design <- c(area(1,1,1,778), area(2,1,2,778), area(3,1,3,1000), area(4,1,4,1000))
(p_E +  xlab("")  + theme(axis.title.x = element_blank(), strip.text=element_text(margin=margin(b = 5)), axis.text.x=element_blank(), plot.margin = margin(b = 1, t = 0))) + 
  (plt_R + xlab("")  + theme(strip.clip = "off", strip.text=element_text(margin=margin(b = 5, t=-20)), axis.title.x = element_blank(), axis.text.x=element_blank(), plot.margin = margin(t = 0, b = 1))) + 
  (plt_SSB + xlab("")  + theme(strip.clip = "off", strip.text=element_text(margin=margin(b = 5, t=-20)), axis.title.x = element_blank(), legend.position="none",axis.text.x=element_blank(), plot.margin = margin(t = 0, b = 1))) + 
  (plt_F + theme(legend.position="none", strip.clip = "off", strip.text = element_text(margin=margin(t=-15)), plot.margin = margin(b = 1, t = 0))) + 
  plot_layout(design = design)
dev.off()
###########################################################################################
