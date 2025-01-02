library(wham)
library(ggplot2)
library(dplyr)
library(patchwork)
source(here::here("R","plot_functions.R"))
theme_set(theme_bw())
theme_update(strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = rel(2)), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)))

fit1 <- readRDS(file.path("results","fits_no_M_re_rev.RDS"))[[2]] #M_1

proj_R1 <- lapply(1:3, \(x) readRDS(here::here("results",paste0("m1_proj_",x,"_R_opt1.RDS"))))
proj_R3 <- lapply(1:3, \(x) readRDS(here::here("results",paste0("m1_proj_",x,"_R_opt3.RDS"))))
proj <- proj_R1[[3]]
se_R1 <- lapply(proj_R1, \(x) TMB:::as.list.sdreport(x$sdrep, what = "Std", report = TRUE))
se_R3 <- lapply(proj_R3, \(x) TMB:::as.list.sdreport(x$sdrep, what = "Std", report = TRUE))

# fjp <- list() 
# for(i in 1:3) fjp[[i]] <- readRDS(here::here("results",paste0("m1_proj_",i,"_sdrep_fjp.RDS")))

#have to do se(predR) by hand
source(here::here("R","plot_functions.R"))
temp <- make_pred_R_cv(proj_R3[[1]], type = 1)
temp <- rbind(temp, make_pred_R_cv(proj_R3[[1]], type = 1, region = 2))
temp <- rbind(temp,make_pred_R_cv(proj_R3[[2]], type = 2))
temp <- rbind(temp,make_pred_R_cv(proj_R3[[2]], type = 2, region = 2))
temp <- rbind(temp,make_pred_R_cv(proj_R3[[3]], type = 3))
temp <- rbind(temp,make_pred_R_cv(proj_R3[[3]], type = 3, region = 2))

source(here::here("R","plot_functions.R"))
brp_results <- lapply(c(proj_R1,proj_R3), get.brp.status.results)


ny_full <- length(proj_R3[[1]]$years_full)
df_annual <- cbind.data.frame(Year = proj_R3[[1]]$years_full, 
  val = exp(c(sapply(c(proj_R1,proj_R3), \(x) c(x$rep$log_SSB_FXSPR)))),
  cv = c(sapply(c(se_R1,se_R3), \(x) c(x$log_SSB_FXSPR))),
  # val = exp(c(sapply(paste0("proj_",1:3,"_R", rep(c(1,3),each = 3)), \(x) c(get(x)$rep$log_SSB_FXSPR)))),
  name = "SSB(F40)",
  Region = rep(c("North","South","Total"), each = ny_full),
  Proj_type = rep(c("AR1", "Average", "Linear"), each = 3*ny_full), 
  XSPR_R_type = rep(c("Random Effect", "Expected"), each = 3*3*ny_full))
df_annual <- rbind(df_annual, cbind.data.frame(Year = proj_R3[[1]]$years_full, 
  val = exp(c(sapply(c(proj_R1,proj_R3), \(x) c(x$rep$log_Fbar_XSPR[,5:7])))),
  cv = c(sapply(c(se_R1,se_R3), \(x) c(x$log_Fbar_XSPR[,5:7]))),
  name = "F40",
  Region = rep(c("North","South","Total"), each = ny_full),
  Proj_type = rep(c("AR1", "Average", "Linear"), each = 3*ny_full), 
  XSPR_R_type = rep(c("Random Effect", "Expected"), each = 3*3*ny_full)))
df_annual <- rbind(df_annual, cbind.data.frame(Year = proj_R3[[1]]$years_full, 
  val = c(sapply(c(proj_R1,proj_R3), \(x) c(x$rep$SSB, apply(x$rep$SSB,1,sum)))),
  cv = c(sapply(c(se_R1,se_R3), \(x) c(x$log_SSB, x$log_SSB_all))),
  name = "SSB",
  Region = rep(c("North","South","Total"), each = ny_full),
  Proj_type = rep(c("AR1", "Average", "Linear"), each = 3*ny_full), 
  XSPR_R_type = rep(c("Random Effect", "Expected"), each = 3*3*ny_full)))
df_annual <- rbind(df_annual, cbind.data.frame(Year = proj_R3[[1]]$years_full, 
  val = c(sapply(c(proj_R1,proj_R3), \(x) c(x$rep$Fbar[,5:7]))),
  cv = c(sapply(c(se_R1,se_R3), \(x) c(x$log_Fbar[,5:7]))),
  name = "Fbar",
  Region = rep(c("North","South","Total"), each = ny_full),
  Proj_type = rep(c("AR1", "Average", "Linear"), each = 3*ny_full), 
  XSPR_R_type = rep(c("Random Effect", "Expected"), each = 3*3*ny_full)))

df_annual <- rbind(df_annual, cbind.data.frame(Year = proj_R3[[1]]$years_full, 
  val = c(sapply(c(proj_R1,proj_R3), \(x) c(x$rep$NAA[1,1,,1],x$rep$NAA[2,2,,1]))),
  cv = c(sapply(c(se_R1,se_R3), \(x) c(x$log_NAA_rep[1,1,,1],x$log_NAA_rep[2,2,,1]))),
  name = "Random Effect",
  Region = rep(c("North","South"), each = ny_full),
  Proj_type = rep(c("AR1", "Average", "Linear"), each = 2*ny_full), 
  XSPR_R_type = rep(c("Random Effect", "Expected"), each = 2*3*ny_full)))

# dim(subset(x, Proj_type == "Last" & Region == "North" & Year < 2022))

# mean(apply(proj_R1[[1]]$rep$FAA[1:2,35,6:7],2,sum))
# proj_R1[[1]]$rep$Fbar[35,5]
# se_R1[[1]]$log_F_tot


df_annual <- rbind(df_annual, cbind.data.frame(Year = proj_R3[[1]]$years_full, 
  val = temp[,1],
  cv = temp[,2],
  name = "Expected",
  Region = rep(c("North","South"), each = ny_full),
  Proj_type = rep(c("AR1", "Average", "Linear"), each = 2*ny_full), 
  XSPR_R_type = rep(c("Random Effect", "Expected"), each = 2*3*ny_full)))


df_annual <- rbind(df_annual, cbind.data.frame(Year = c(sapply(brp_results, \(x) c(x$df.status$Year))),
  val = exp(c(sapply(brp_results, \(x) x$df.status$val))),
  cv = c(sapply(brp_results, \(x) x$df.status$se)),
  name = c(sapply(brp_results, \(x) c("SSB_status", "F_status")[match(x$df.status$type, c("SSB","Fbar"))])),
  Region = c(sapply(brp_results, \(x) x$df.status$region)),
  Proj_type = rep(c("AR1", "Average", "Linear"), each = NROW(brp_results[[1]]$df.status)), 
  XSPR_R_type = rep(c("Random Effect", "Expected"), each = 3*NROW(brp_results[[1]]$df.status))))

df_annual$lo <- df_annual$val*exp(-qnorm(0.975)*df_annual$cv)
df_annual$hi <- df_annual$val*exp(qnorm(0.975)*df_annual$cv)

df_bta <- cbind.data.frame(Year = proj_R1[[3]]$input$years_Ecov,
  val = c(sapply(c(proj_R1,proj_R3), \(x) c(x$rep$Ecov_x))),
  cv = c(sapply(c(se_R1,se_R3), \(x) c(x$Ecov_x))),
  name = "BTA",
  Region = rep(c("North","South"), each = length(proj_R1[[3]]$input$years_Ecov)),
  Proj_type = rep(c("AR1", "Average", "Linear"), each = 2*length(proj_R1[[3]]$input$years_Ecov)), 
  XSPR_R_type = rep(c("Random Effect", "Expected"), each = 2*3*length(proj_R1[[3]]$input$years_Ecov)))

ind_x <- matrix(which(names(fit1$sdrep$value) == "Ecov_x"), length(fit1$input$years_Ecov), 2)
ind_x <- ind_x[which(fit1$input$years_Ecov %in% tail(fit1$years,5)),]

ind <- which(df_bta$Proj_type == "Average" & df_bta$Year>2021 & df_bta$Region == "North")
df_bta$val[ind] <- mean(fit1$sdrep$value[ind_x[,1]])
df_bta$cv[ind] <- sqrt(sum(fit1$sdrep$cov[ind_x[,1],ind_x[,1]]))/5
ind <- which(df_bta$Proj_type == "Average" & df_bta$Year>2021 & df_bta$Region == "South")
df_bta$val[ind] <- mean(fit1$sdrep$value[ind_x[,2]])
df_bta$cv[ind] <- sqrt(sum(fit1$sdrep$cov[ind_x[,2],ind_x[,2]]))/5

ind <- which(df_bta$Proj_type == "Linear" & df_bta$Year>2021 & df_bta$Region == "North")
df_bta$val[ind] <- proj_R1[[3]]$rep$Ecov_out_R[1,which(!proj_R1[[3]]$input$years_full %in% proj_R1[[3]]$input$years),1]
df_bta$cv[ind] <- 0
ind <- which(df_bta$Proj_type == "Linear" & df_bta$Year>2021 & df_bta$Region == "South")
df_bta$val[ind] <- proj_R1[[3]]$rep$Ecov_out_R[1,which(!proj_R1[[3]]$input$years_full %in% proj_R1[[3]]$input$years),2]
df_bta$cv[ind] <- 0

df_bta$lo <- df_bta$val -qnorm(0.975)*df_bta$cv
df_bta$hi <- df_bta$val + qnorm(0.975)*df_bta$cv

# df_annual <- rbind(df_annual, temp)

# df_proj <- subset(df_annual, Year %in% 2019:2031)

this_df <- subset(df_bta, XSPR_R_type == "Random Effect" &  Year %in% 2019:2031)
this_df <- subset(this_df, (Year >2021 & Proj_type != "AR1") | Proj_type == "AR1")
this_df <- subset(this_df, Region == "North")
# proj_df <- subset(this_df, Year >=2021)
# mod_df <- subset(this_df, Proj_type == "AR1")

BTA_proj <- ggplot(this_df, aes(x = Year, y = val, colour = Proj_type, fill = Proj_type)) + xlab("Year") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + 
    facet_grid(~Region) + 
    labs(color = "Projection Type", fill = "Projection Type") + 
    ylab("Bottom Temperature\nAnomaly") + geom_vline(xintercept = 2021, linetype = 2, alpha = 0.4, linewidth = 1.5) +
    geom_line() + 
    geom_point(shape = 21, size = 2) + 
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha=0.4, linetype = 0)
BTA_proj

this_df <- subset(df_annual, Year %in% 2019:2031 & name == "Random Effect" & XSPR_R_type == "Expected")
this_df$val <- this_df$val/1000
this_df$lo <- this_df$lo/1000
this_df$hi <- this_df$hi/1000
this_df$Proj_type <- factor(this_df$Proj_type)
model_df <- subset(this_df, Year <=2021)
# proj_df <- subset(this_df, Year >=2021)
proj_df <- subset(this_df, Year >2021)
north_ar1_df <- subset(this_df, Region == "North" & Proj_type == "AR1")
south_df <- subset(this_df, Region == "South" & Proj_type == "AR1")
north_proj_df <- subset(proj_df, Region == "North" & Proj_type != "AR1")

cnms <- c(North = "", South = "South")
R_RE_proj <- ggplot(this_df, aes(x = Year, y = val)) + xlab("Year") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + 
    labs(color = "Projection\nType", fill = "Projection\nType") + coord_cartesian(ylim = c(0,80)) +
    facet_grid(~Region, labeller = labeller(Region = cnms)) + 
    ylab(expression(Recruitment~~Random~~Effect~~(10^6))) + geom_vline(xintercept = 2021, linetype = 2, alpha = 0.4, linewidth = 1.5) +
    geom_line(data=south_df, mapping = aes(x = Year, y = val), colour = "black") + 
    geom_point(data=south_df, mapping = aes(x = Year, y = val), colour = "black", fill = "black", shape = 21, size = 2) + 
    geom_ribbon(data = south_df, aes(x = Year, ymin = lo, ymax = hi), fill = alpha("black", 0.4), linetype = 0) +
    geom_line(data=north_ar1_df, mapping = aes(x = Year, y = val, colour = Proj_type)) + 
    geom_point(data=north_ar1_df, mapping = aes(x = Year, y = val, colour = Proj_type, fill = Proj_type), shape = 21, size = 2) + 
    geom_ribbon(data = north_ar1_df, aes(ymin = lo, ymax = hi, fill = Proj_type), alpha=0.4, linetype = 0) +
    geom_line(data=north_proj_df, mapping = aes(x = Year, y = val, colour = Proj_type)) + 
    geom_point(data=north_proj_df, mapping = aes(x = Year, y = val, colour = Proj_type, fill = Proj_type), shape = 21, size = 2) + 
    geom_ribbon(data = north_proj_df, aes(ymin = lo, ymax = hi, fill = Proj_type), alpha=0.4, linetype = 0)
R_RE_proj

this_df <- subset(df_annual, Year %in% 2019:2031 & name == "Expected" & XSPR_R_type == "Expected")
this_df$val <- this_df$val/1000
this_df$lo <- this_df$lo/1000
this_df$hi <- this_df$hi/1000
# north_df <- subset(this_df, Region == "North")
# south_df <- subset(this_df, Region == "South" & Proj_type == "AR1")
model_df <- subset(this_df, Year <=2021)
# proj_df <- subset(this_df, Year >=2021)
proj_df <- subset(this_df, Year >2021)
north_ar1_df <- subset(this_df, Region == "North" & Proj_type == "AR1")
south_df <- subset(this_df, Region == "South" & Proj_type == "AR1")
north_proj_df <- subset(proj_df, Region == "North" & Proj_type != "AR1")
cnms <- c(North = "", South = "")
R_Exp_proj <- ggplot(this_df, aes(x = Year, y = val)) + xlab("Year") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + 
    labs(color = "Projection\nType", fill = "Projection\nType") + coord_cartesian(ylim = c(0,80)) + 
    facet_grid(~Region, labeller = labeller(Region = cnms)) + ylab(expression(Expected~~Recruitment~~(10^6))) + geom_vline(xintercept = 2021, linetype = 2, alpha = 0.4, linewidth = 1.5) +
    geom_line(data=south_df, mapping = aes(x = Year, y = val), colour = "black") + 
    geom_line(data=south_df, mapping = aes(x = Year, y = val), colour = "black") + 
    geom_point(data=south_df, mapping = aes(x = Year, y = val), colour = "black", fill = "black", shape = 21, size = 2) + 
    geom_ribbon(data = south_df, aes(x = Year, ymin = lo, ymax = hi), fill = alpha("black", 0.4), linetype = 0) +
    geom_line(data=north_ar1_df, mapping = aes(x = Year, y = val, colour = Proj_type)) + 
    geom_point(data=north_ar1_df, mapping = aes(x = Year, y = val, colour = Proj_type, fill = Proj_type), shape = 21, size = 2) + 
    geom_ribbon(data = north_ar1_df, aes(ymin = lo, ymax = hi, fill = Proj_type), alpha=0.4, linetype = 0) +
    geom_line(data=north_proj_df, mapping = aes(x = Year, y = val, colour = Proj_type)) + 
    geom_point(data=north_proj_df, mapping = aes(x = Year, y = val, colour = Proj_type, fill = Proj_type), shape = 21, size = 2) + 
    geom_ribbon(data = north_proj_df, aes(ymin = lo, ymax = hi, fill = Proj_type), alpha=0.4, linetype = 0)
R_Exp_proj

cairo_pdf(file.path("paper", "proj_ecov_Recruit_results.pdf"), width = 20, height = 20)
design <- c(area(1,1,1,572), area(2,1,2,1000), area(3,1,3,1000))
(BTA_proj +  xlab("")  + theme(axis.title.x = element_blank(), strip.text=element_text(margin=margin(b = 5)), axis.text.x=element_blank(), plot.margin = margin(b = 1, t = 0))) + 
  (R_RE_proj + xlab("")  + theme(strip.clip = "off", strip.text=element_text(margin=margin(b = 5, t=-20)), axis.title.x = element_blank(), legend.position="none", axis.text.x=element_blank(), plot.margin = margin(t = 0, b = 1))) + 
  (R_Exp_proj + xlab("Year")  + theme(legend.position="none", strip.clip = "off", strip.text = element_text(margin=margin(t=-15)), plot.margin = margin(b = 1, t = 0))) + 
  plot_layout(design = design)
dev.off()



this_df <- subset(df_annual, Year %in% 2019:2031 & name == "Fbar" & XSPR_R_type == "Expected" & Proj_type == "AR1")
# model_df <- subset(this_df, Year <=2021)
# proj_line_df <- subset(this_df, Year >=2021)
# proj_point_df <- subset(this_df, Year >2021)
Fbar_proj <- ggplot(this_df, aes(x = Year, y = val)) + xlab("Year") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + 
    labs(color = "Projection\nType", fill = "Projection\nType") + 
    facet_grid(~Region) + ylab(expression(bar(italic(F))~~"(Ages 6-7)")) + 
    geom_vline(xintercept = 2021, linetype = 2, alpha = 0.4, linewidth = 1.5) +
    geom_point(this_df, mapping = aes(x = Year, y = val), shape = 21, size = 2, fill = "black", colour = "black") + 
    geom_line(this_df, mapping = aes(x = Year, y = val), colour = "black")  + 
    geom_ribbon(this_df, mapping = aes(x = Year, y = val, ymin = lo, ymax = hi), colour = "black", alpha=0.4, linetype = 0)  
    # geom_point(proj_point_df, mapping = aes(x = Year, y = val, colour = Proj_type, fill = Proj_type), shape = 21, size = 2) + 
    # geom_line(proj_line_df, mapping = aes(x = Year, y = val, colour = Proj_type)) + geom_ribbon(proj_line_df, mapping = aes(x = Year, y = val, ymin = lo, ymax = hi, fill = Proj_type), alpha=0.4, linetype = 0)
Fbar_proj

cnms <- c(North = "", South = "", Total = "")
this_df <- subset(df_annual, Year %in% 2019:2031 & name == "SSB" & XSPR_R_type == "Expected")
this_df$val <- this_df$val/1000
this_df$lo <- this_df$lo/1000
this_df$hi <- this_df$hi/1000
model_df <- subset(this_df, Proj_type == "AR1")
# south_df <- subset(this_df, Region == "South" & Proj_type == "AR1")
# proj_df <- subset(this_df, Year >=2021)
proj_df <- subset(this_df, Year > 2021 & Proj_type != "AR1")
# south_df <- subset(this_df, Region == "South" & Proj_type == "AR1")
# not_south_ar1_df <- subset(this_df, Region != "South" & Proj_type == "AR1")
# not_south_proj_df <- subset(proj_df, Region != "South" & Proj_type != "AR1")
# proj_line_df <- subset(this_df, Year >2021)
# proj_point_df <- subset(this_df, Year >2021)
SSB_proj <- ggplot(this_df, aes(x = Year, y = val)) + xlab("Year") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + 
    labs(color = "Projection Type", fill = "Projection Type") + 
    facet_grid(~Region, labeller = labeller(Region = cnms)) + ylab("SSB (kmt)") + 
    geom_vline(xintercept = 2021, linetype = 2, alpha = 0.4, linewidth = 1.5) +
    geom_point(model_df, mapping = aes(x = Year, y = val, color = Proj_type, fill = Proj_type), shape = 21, size = 2) + 
    geom_line(model_df, mapping = aes(x = Year, y = val, color = Proj_type))  + 
    geom_ribbon(model_df, mapping = aes(x = Year, y = val, ymin = lo, ymax = hi, fill = Proj_type), alpha=0.4, linetype = 0) + 
    geom_point(proj_df, mapping = aes(x = Year, y = val, colour = Proj_type, fill = Proj_type), shape = 21, size = 2) + 
    geom_line(proj_df, mapping = aes(x = Year, y = val, colour = Proj_type)) + 
    geom_ribbon(proj_df, mapping = aes(x = Year, y = val, ymin = lo, ymax = hi, fill = Proj_type), alpha=0.4, linetype = 0)
SSB_proj

cairo_pdf(file.path("paper", "proj_F_SSB.pdf"), width = 20, height = 20)
design <- c(area(1,1,1,1), area(2,1,2,1))
(Fbar_proj +  xlab("")  + theme(axis.title.x = element_blank(), strip.text=element_text(margin=margin(b = 5)), axis.text.x=element_blank(), plot.margin = margin(b = 1, t = 0))) + 
  (SSB_proj + xlab("Year")  + theme(legend.position="none", strip.clip = "off", strip.text = element_text(margin=margin(t=-15)), plot.margin = margin(b = 1, t = 0))) + 
  plot_layout(design = design)
dev.off()

this_df <- subset(df_annual, Year %in% 2019:2031 & name == "SSB(F40)")
this_df$val <- this_df$val/1000
this_df$lo <- this_df$lo/1000
this_df$hi <- this_df$hi/1000
model_df <- subset(this_df, Proj_type == "AR1")
proj_df <- subset(this_df, Year > 2021 & Proj_type != "AR1")
SSB40_proj <- ggplot(this_df, aes(x = Year, y = val)) + xlab("Year") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + 
    labs(color = "Projection Type", fill = "Projection Type") + 
    facet_grid(Region~XSPR_R_type) + ylab(expression(SSB[40*"%"]~~(kmt))) + 
    geom_vline(xintercept = 2021, linetype = 2, alpha = 0.4, linewidth = 1.5) +
    geom_point(model_df, mapping = aes(x = Year, y = val, colour = Proj_type, fill = Proj_type), shape = 21, size = 2) + 
    geom_line(model_df, mapping = aes(x = Year, y = val, colour = Proj_type))  + 
    geom_ribbon(model_df, mapping = aes(x = Year, y = val, ymin = lo, ymax = hi, fill = Proj_type), alpha=0.4, linetype = 0) + 
    geom_point(proj_df, mapping = aes(x = Year, y = val, colour = Proj_type, fill = Proj_type), shape = 21, size = 2) + 
    geom_line(proj_df, mapping = aes(x = Year, y = val, colour = Proj_type)) + 
    geom_ribbon(proj_df, mapping = aes(x = Year, y = val, ymin = lo, ymax = hi, fill = Proj_type), alpha=0.4, linetype = 0)
SSB40_proj

cairo_pdf(file.path("paper", "proj_SSB40_results.pdf"), width = 20, height = 20)
SSB40_proj
dev.off()

this_df <- subset(df_annual, Year %in% 2019:2031 & name == "F40")
# model_df <- subset(this_df, Year <=2021)
model_df <- subset(this_df, Proj_type == "AR1")
proj_df <- subset(this_df, Year > 2021 & Proj_type != "AR1")
# proj_line_df <- subset(this_df, Year >=2021)
# proj_point_df <- subset(this_df, Year >2021)
F40_proj <- ggplot(this_df, aes(x = Year, y = val)) + xlab("Year") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + 
    labs(color = "Projection Type", fill = "Projection Type") + 
    facet_grid(Region~XSPR_R_type, scales="free_y") + ylab(expression(bar(italic(F))[40*"%"]~~"(Ages 6-7)")) + 
    geom_vline(xintercept = 2021, linetype = 2, alpha = 0.4, linewidth = 1.5) +
    geom_point(model_df, mapping = aes(x = Year, y = val, colour = Proj_type, fill = Proj_type), shape = 21, size = 2) + 
    geom_line(model_df, mapping = aes(x = Year, y = val, colour = Proj_type))  + 
    geom_ribbon(model_df, mapping = aes(x = Year, y = val, ymin = lo, ymax = hi, fill = Proj_type), alpha=0.4, linetype = 0) + 
    geom_point(proj_df, mapping = aes(x = Year, y = val, colour = Proj_type, fill = Proj_type), shape = 21, size = 2) + 
    geom_line(proj_df, mapping = aes(x = Year, y = val, colour = Proj_type)) + 
    geom_ribbon(proj_df, mapping = aes(x = Year, y = val, ymin = lo, ymax = hi, fill = Proj_type), alpha=0.4, linetype = 0)
F40_proj

cairo_pdf(file.path("paper", "proj_F40_results.pdf"), width = 20, height = 20)
F40_proj
dev.off()

cairo_pdf(file.path("paper", "proj_brps.pdf"), width = 24, height = 16)
design <- c(area(1,1,1,1), area(1,2,1,2))
(F40_proj +  xlab("Year")  + theme(legend.position="none", axis.text.x=element_text(size = rel(0.7)))) + 
  (SSB40_proj + xlab("Year") + theme(axis.text.x=element_text(size = rel(0.7)))) + 
  plot_layout(design = design)
dev.off()

this_df <- subset(df_annual, Year %in% 2019:2031 & name == "SSB_status")
model_df <- subset(this_df, Proj_type == "AR1")
proj_df <- subset(this_df, Year > 2021 & Proj_type != "AR1")
# model_df <- subset(this_df, Year <=2021)
# proj_line_df <- subset(this_df, Year >=2021)
# proj_point_df <- subset(this_df, Year >2021)
SSB_status_proj <- ggplot(this_df, aes(x = Year, y = val)) + xlab("Year") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + 
    labs(color = "Projection Type", fill = "Projection Type") + 
    facet_grid(Region~XSPR_R_type, scales="free_y") + ylab(bquote(SSB*"/"*SSB[40*"%"])) + 
    geom_vline(xintercept = 2021, linetype = 2, alpha = 0.4, linewidth = 1.5) + geom_hline(yintercept = 1, linetype = 1, colour = alpha("gray", 0.4), linewidth = 1.5) +
    geom_point(model_df, mapping = aes(x = Year, y = val, colour = Proj_type, fill = Proj_type), shape = 21, size = 2) + 
    geom_line(model_df, mapping = aes(x = Year, y = val, colour = Proj_type))  + 
    geom_ribbon(model_df, mapping = aes(x = Year, y = val, ymin = lo, ymax = hi, fill = Proj_type), alpha=0.4, linetype = 0) + 
    geom_point(proj_df, mapping = aes(x = Year, y = val, colour = Proj_type, fill = Proj_type), shape = 21, size = 2) + 
    geom_line(proj_df, mapping = aes(x = Year, y = val, colour = Proj_type)) + 
    geom_ribbon(proj_df, mapping = aes(x = Year, y = val, ymin = lo, ymax = hi, fill = Proj_type), alpha=0.4, linetype = 0)
SSB_status_proj

cairo_pdf(file.path("paper", "proj_SSB_status_results.pdf"), width = 20, height = 20)
SSB_status_proj
dev.off()



this_df <- subset(df_annual, Year %in% 2019:2031 & name == "F_status")
model_df <- subset(this_df, Proj_type == "AR1")
proj_df <- subset(this_df, Year > 2021 & Proj_type != "AR1")
# model_df <- subset(this_df, Year <=2021)
# proj_line_df <- subset(this_df, Year >=2021)
# proj_point_df <- subset(this_df, Year >2021)
F_status_proj <- ggplot(this_df, aes(x = Year, y = val)) + xlab("Year") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + 
    labs(color = "Projection Type", fill = "Projection Type") + 
    facet_grid(Region~XSPR_R_type, scales="free_y") + ylab(bquote(bar(italic(F))*"/"*bar(italic(F))[40*"%"]~~"(Ages 6-7)")) + 
    geom_vline(xintercept = 2021, linetype = 2, alpha = 0.4, linewidth = 1.5) +
    geom_point(model_df, mapping = aes(x = Year, y = val, colour = Proj_type, fill = Proj_type), shape = 21, size = 2) + 
    geom_line(model_df, mapping = aes(x = Year, y = val, colour = Proj_type))  + 
    geom_ribbon(model_df, mapping = aes(x = Year, y = val, ymin = lo, ymax = hi, fill = Proj_type), alpha=0.4, linetype = 0) + 
    geom_point(proj_df, mapping = aes(x = Year, y = val, colour = Proj_type, fill = Proj_type), shape = 21, size = 2) + 
    geom_line(proj_df, mapping = aes(x = Year, y = val, colour = Proj_type)) + 
    geom_ribbon(proj_df, mapping = aes(x = Year, y = val, ymin = lo, ymax = hi, fill = Proj_type), alpha=0.4, linetype = 0)
F_status_proj

cairo_pdf(file.path("paper", "proj_F_status_results.pdf"), width = 20, height = 20)
F_status_proj
dev.off()

this_df <- subset(df_annual, name %in% c("Expected", "Random Effect") & XSPR_R_type == "Expected")
model_df <- subset(this_df, (Region == "South" | (Region == "North" & Year <=2021)) & Proj_type == "AR1")
proj_line_df <- subset(this_df, Region == "North" & Year >=2021)
proj_point_df <- subset(this_df, Region == "North" & Year >2021)
plt_R_cv <- ggplot(this_df, aes(x = Year, y = cv)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", alpha = 0.5) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", alpha = 0.5) + 
    labs(color = "Projection Type", fill = "Projection Type")  + geom_vline(xintercept = 2021, linetype = 2, alpha = 0.4, linewidth = 1.5) + 
    geom_point(model_df, mapping = aes(x = Year, y = cv), shape = 21, size = 2, fill = "black", colour = "black") + 
    geom_line(model_df, mapping = aes(x = Year, y = cv), colour = "black")  + 
    geom_point(proj_point_df, mapping = aes(x = Year, y = cv, colour = Proj_type, fill = Proj_type), shape = 21, size = 2) + 
    geom_line(proj_line_df, mapping = aes(x = Year, y = cv, colour = Proj_type))  + 
    ylab("CV(Recruitment)") + xlab("Year") + facet_grid(name ~ Region)
plt_R_cv

cnms <- c(North = "", South = "", Total = "Total")
this_df <- subset(df_annual, name == "Fbar" & XSPR_R_type == "Expected")
model_df <- subset(this_df, Year <=2021 )
proj_line_df <- subset(this_df, Year >=2021)
proj_point_df <- subset(this_df, Year >2021)
plt_F_cv <- ggplot(this_df, aes(x = Year, y = cv)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", alpha = 0.5) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", alpha = 0.5) + 
    labs(color = "Projection Type", fill = "Projection Type") + geom_vline(xintercept = 2021, linetype = 2, alpha = 0.4, linewidth = 1.5) + 
    geom_point(model_df, mapping = aes(x = Year, y = cv), shape=21, size = 2, colour = "black", fill = "black") + 
    geom_line(model_df, mapping = aes(x = Year, y = cv), colour = "black")  + 
    geom_point(proj_point_df, mapping = aes(x = Year, y = cv, colour = Proj_type, fill = Proj_type), shape=21, size = 2) + 
    geom_line(proj_line_df, mapping = aes(x = Year, y = cv, colour = Proj_type))  + 
    ylab(expression(CV*"("*bar(italic(F))*")")) + xlab("Year") + facet_grid(~ Region, labeller = labeller(Region = cnms))
plt_F_cv

cnms <- c(North = "", South = "", Total = "")
this_df <- subset(df_annual, name == "SSB" & XSPR_R_type == "Expected")
model_df <- subset(this_df, Year <=2021)
proj_line_df <- subset(this_df, Year >=2021)
proj_point_df <- subset(this_df, Year >2021)
plt_SSB_cv <- ggplot(this_df, aes(x = Year, y = cv)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", alpha = 0.5) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", alpha = 0.5) + 
    labs(color = "Projection Type", fill = "Projection Type") + geom_vline(xintercept = 2021, linetype = 2, alpha = 0.4, linewidth = 1.5) + 
    geom_point(model_df, mapping = aes(x = Year, y = cv), shape=21, size = 2, colour = "black", fill = "black") + 
    geom_line(model_df, mapping = aes(x = Year, y = cv), colour = "black")  + 
    geom_point(proj_point_df, mapping = aes(x = Year, y = cv, colour = Proj_type, fill = Proj_type), shape=21, size = 2) + 
    geom_line(proj_line_df, mapping = aes(x = Year, y = cv, colour = Proj_type))  + 
    ylab("CV(SSB)") + xlab("Year") + facet_grid(~ Region, labeller = labeller(Region = cnms))
plt_SSB_cv


cairo_pdf(file.path("paper", "R_SSB_F_cv_results.pdf"), width = 20, height = 15)
design <- c(area(1,1,2,777), area(3,1,3,1000), area(4,1,4,1000))
(plt_R_cv +  xlab("")  + theme(axis.title.x = element_blank(), strip.text=element_text(margin=margin(b = 5)), axis.text.x=element_blank(), plot.margin = margin(b = 1, t = 0))) + 
  (plt_F_cv + xlab("")  + theme(strip.clip = "off", strip.text=element_text(margin=margin(b = 5, t=-20)), axis.title.x = element_blank(), legend.position="none",axis.text.x=element_blank(), plot.margin = margin(t = 0, b = 1))) + 
  (plt_SSB_cv + theme(legend.position="none", strip.clip = "off", strip.text = element_text(margin=margin(t=-15)), plot.margin = margin(b = 1, t = 0))) + 
  plot_layout(design = design)
dev.off()


this_df <- subset(df_annual, name %in% c("Expected", "Random Effect") & XSPR_R_type == "Expected")
model_df <- subset(this_df, (Region == "South" | (Region == "North" & Year <=2021)) & Proj_type == "AR1")
proj_line_df <- subset(this_df, Region == "North" & Year >=2021)
proj_point_df <- subset(this_df, Region == "North" & Year >2021)
plt_R_cv <- ggplot(this_df, aes(x = Year, y = cv)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", alpha = 0.5) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", alpha = 0.5) + 
    labs(color = "Projection Type", fill = "Projection Type")  + geom_vline(xintercept = 2021, linetype = 2, alpha = 0.4, linewidth = 1.5) + 
    geom_point(model_df, mapping = aes(x = Year, y = cv), shape = 21, size = 2, fill = "black", colour = "black") + 
    geom_line(model_df, mapping = aes(x = Year, y = cv), colour = "black")  + 
    geom_point(proj_point_df, mapping = aes(x = Year, y = cv, colour = Proj_type, fill = Proj_type), shape = 21, size = 2) + 
    geom_line(proj_line_df, mapping = aes(x = Year, y = cv, colour = Proj_type))  + 
    ylab("CV(Recruitment)") + xlab("Year") + facet_grid(name ~ Region)
plt_R_cv

this_df <- subset(df_annual, name == "F40")
model_df <- subset(this_df, Year <=2021)
proj_line_df <- subset(this_df, Year >=2021)
proj_point_df <- subset(this_df, Year >2021)
F40_cv_proj <- ggplot(this_df, aes(x = Year, y = cv, colour = Proj_type, fill = Proj_type)) + xlab("Year") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", alpha = 0.5) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", alpha = 0.5) + 
    labs(color = "Projection Type", fill = "Projection Type") + 
    facet_grid(Region~XSPR_R_type, scales="free_y") + ylab(bquote(CV*"("*bar(italic(F))[40*"%"]*")")) + 
    geom_vline(xintercept = 2021, linetype = 2, alpha = 0.4, linewidth = 1.5) +
    geom_point(model_df, mapping = aes(x = Year, y = cv), shape = 21, size = 2, fill = "black", colour = "black") + 
    geom_line(model_df, mapping = aes(x = Year, y = cv), colour = "black")  + 
    geom_point(proj_point_df, mapping = aes(x = Year, y = cv, colour = Proj_type, fill = Proj_type), shape = 21, size = 2) + 
    geom_line(proj_line_df, mapping = aes(x = Year, y = cv, colour = Proj_type))
F40_cv_proj

cairo_pdf(file.path("paper", "proj_F40_CV.pdf"), width = 20, height = 20)
F40_cv_proj
dev.off()

this_df <- subset(df_annual, name == "SSB(F40)")
model_df <- subset(this_df, Year <=2021)
proj_line_df <- subset(this_df, Year >=2021)
proj_point_df <- subset(this_df, Year >2021)
SSB40_cv_proj <- ggplot(this_df, aes(x = Year, y = cv, colour = Proj_type, fill = Proj_type)) + xlab("Year") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", alpha = 0.5) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", alpha = 0.5) + 
    labs(color = "Projection Type", fill = "Projection Type") + 
    facet_grid(Region~XSPR_R_type, scales="free_y") + ylab(bquote(CV*"("*SSB[40*"%"]*")")) + 
    geom_vline(xintercept = 2021, linetype = 2, alpha = 0.4, linewidth = 1.5) +
    geom_point(model_df, mapping = aes(x = Year, y = cv), shape = 21, size = 2, fill = "black", colour = "black") + 
    geom_line(model_df, mapping = aes(x = Year, y = cv), colour = "black")  + 
    geom_point(proj_point_df, mapping = aes(x = Year, y = cv, colour = Proj_type, fill = Proj_type), shape = 21, size = 2) + 
    geom_line(proj_line_df, mapping = aes(x = Year, y = cv, colour = Proj_type))
SSB40_cv_proj

cairo_pdf(file.path("paper", "proj_SSB40_CV.pdf"), width = 20, height = 20)
SSB40_cv_proj
dev.off()

# this_df <- subset(df_annual, Year %in% 2019:2031 & name == "F_status")
this_df <- subset(df_annual, name == "F_status")
model_df <- subset(this_df, Year <=2021)
proj_line_df <- subset(this_df, Year >=2021)
proj_point_df <- subset(this_df, Year >2021)
F_status_cv_proj <- ggplot(this_df, aes(x = Year, y = cv, colour = Proj_type, fill = Proj_type)) + xlab("Year") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", alpha = 0.5) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", alpha = 0.5) + 
    labs(color = "Projection Type", fill = "Projection Type") + 
    facet_grid(Region~XSPR_R_type, scales="free_y") + ylab(bquote(CV*"("*bar(italic(F))*"/"*bar(italic(F))[40*"%"]*")")) + 
    geom_vline(xintercept = 2021, linetype = 2, alpha = 0.4, linewidth = 1.5) +
    geom_point(model_df, mapping = aes(x = Year, y = cv), shape = 21, size = 2, fill = "black", colour = "black") + 
    geom_line(model_df, mapping = aes(x = Year, y = cv), colour = "black")  + 
    geom_point(proj_point_df, mapping = aes(x = Year, y = cv, colour = Proj_type, fill = Proj_type), shape = 21, size = 2) + 
    geom_line(proj_line_df, mapping = aes(x = Year, y = cv, colour = Proj_type))
F_status_cv_proj

cairo_pdf(file.path("paper", "proj_F_status_CV.pdf"), width = 20, height = 20)
F_status_cv_proj
dev.off()

this_df <- subset(df_annual, name == "SSB_status")
model_df <- subset(this_df, Year <=2021)
proj_line_df <- subset(this_df, Year >=2021)
proj_point_df <- subset(this_df, Year >2021)
SSB_status_cv_proj <- ggplot(this_df, aes(x = Year, y = cv, colour = Proj_type, fill = Proj_type)) + xlab("Year") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", alpha = 0.5) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", alpha = 0.5) + 
    labs(color = "Projection Type", fill = "Projection Type") + 
    facet_grid(Region~XSPR_R_type, scales="free_y") + ylab(bquote(CV*"("*SSB*"/"*SSB[40*"%"]*")")) + 
    geom_vline(xintercept = 2021, linetype = 2, alpha = 0.4, linewidth = 1.5) +
    geom_point(model_df, mapping = aes(x = Year, y = cv), shape = 21, size = 2, fill = "black", colour = "black") + 
    geom_line(model_df, mapping = aes(x = Year, y = cv), colour = "black")  + 
    geom_point(proj_point_df, mapping = aes(x = Year, y = cv, colour = Proj_type, fill = Proj_type), shape = 21, size = 2) + 
    geom_line(proj_line_df, mapping = aes(x = Year, y = cv, colour = Proj_type))
SSB_status_cv_proj

cairo_pdf(file.path("paper", "proj_SSB_status_CV.pdf"), width = 20, height = 20)
SSB_status_cv_proj
dev.off()

#######################################################################################
#old

proj <- list() 
for(i in 1:3) proj[[i]] <- readRDS(here::here("results",paste0("m1_proj_",i,".RDS")))
fjp <- list() 
for(i in 1:3) fjp[[i]] <- readRDS(here::here("results",paste0("m1_proj_",i,"_sdrep_fjp.RDS")))
for(i in 1:3) proj[[i]]$sdrep <- fjp[[i]]

temp <- cbind.data.frame(est3 = TMB:::as.list.sdreport(proj[[3]]$sdrep, report = TRUE, what = "Est")$Ecov_x[,1])
temp <- cbind(temp, se3 = TMB:::as.list.sdreport(proj[[3]]$sdrep, report = TRUE, what = "Std")$Ecov_x[,1])
temp <- cbind(est1 = NA, se1 = NA, est2 = NA, se2 = NA, temp)
temp$est1[which(proj[[3]]$input$years_Ecov %in% proj[[1]]$input$years_Ecov)] <- TMB:::as.list.sdreport(proj[[1]]$sdrep, report = TRUE, what = "Est")$Ecov_x[,1]
temp$se1[which(proj[[3]]$input$years_Ecov %in% proj[[1]]$input$years_Ecov)] <- TMB:::as.list.sdreport(proj[[1]]$sdrep, report = TRUE, what = "Std")$Ecov_x[,1]
temp$est2[which(proj[[3]]$input$years_Ecov %in% proj[[2]]$input$years_Ecov)] <- TMB:::as.list.sdreport(proj[[2]]$sdrep, report = TRUE, what = "Est")$Ecov_x[,1]
temp$se2[which(proj[[3]]$input$years_Ecov %in% proj[[2]]$input$years_Ecov)] <- TMB:::as.list.sdreport(proj[[2]]$sdrep, report = TRUE, what = "Std")$Ecov_x[,1]

temp <- temp[which(proj[[3]]$input$years_Ecov %in% proj[[2]]$input$years_full), ]
temp <- cbind(year = proj[[2]]$input$years_full, temp)
rownames(temp) <- proj[[2]]$input$years_full

fit1 <- readRDS(file.path("results","fits_no_M_re_rev.RDS"))[[2]]

ind <- matrix(which(names(fit1$sdrep$value) == "Ecov_x"), length(fit1$input$years_Ecov), 2)[,1]
ind <- ind[which(fit1$input$years_Ecov %in% tail(fit1$years,5))]

proj_est2 <- mean(fit1$sdrep$value[ind])
proj_est2
proj[[2]]$rep$Ecov_out_R[1,,1] #same

proj_se2 <- sqrt(sum(fit1$sdrep$cov[ind,ind]))/5

temp$est2[which(!proj[[2]]$input$years_full %in% fit1$input$years)] <- proj_est2
temp$se2[which(!proj[[2]]$input$years_full %in% fit1$input$years)] <- proj_se2

proj[[3]]$rep$Ecov_out_R[1,,1] 
proj_est3 <- proj[[3]]$rep$Ecov_out_R[1,which(!proj[[2]]$input$years_full %in% fit1$input$years),1]
temp$est3[which(!proj[[2]]$input$years_full %in% fit1$input$years)] <- proj_est3
temp$se3[which(!proj[[2]]$input$years_full %in% fit1$input$years)] <- 0 #input value, not estimated internally

temp$R1 <- proj[[1]]$rep$NAA[1,1,,1]/1000
temp$predR1 <- proj[[1]]$rep$pred_NAA[1,1,,1]/1000
temp$R1.cv <- TMB:::as.list.sdreport(proj[[1]]$sdrep, report = TRUE, what = "Std")$log_NAA_rep[1,1,,1]
temp$R2 <- proj[[2]]$rep$NAA[1,1,,1]/1000
temp$predR2 <- proj[[2]]$rep$pred_NAA[1,1,,1]/1000
temp$R2.cv <- TMB:::as.list.sdreport(proj[[2]]$sdrep, report = TRUE, what = "Std")$log_NAA_rep[1,1,,1]
temp$R3 <- proj[[3]]$rep$NAA[1,1,,1]/1000
temp$predR3 <- proj[[3]]$rep$pred_NAA[1,1,,1]/1000
temp$R3.cv <- TMB:::as.list.sdreport(proj[[3]]$sdrep, report = TRUE, what = "Std")$log_NAA_rep[1,1,,1]

#have to do se(predR) by hand
x <- list(make_pred_R_cv(proj[[1]], fjp[[1]], type = 1))
x[[2]] <- make_pred_R_cv(proj[[2]], fjp[[2]], type = 2)
x[[3]] <- make_pred_R_cv(proj[[3]], fjp[[3]], type = 3)
temp$predR1 
x[[1]][,"pred_R"]
temp$predR2 
x[[2]][,"pred_R"]
temp$predR3 
x[[3]][,"pred_R"]
temp$predR1.cv <- x[[1]][,"cv"]
temp$predR2.cv <- x[[2]][,"cv"]
temp$predR3.cv <- x[[3]][,"cv"]
# temp <- temp[which(proj[[2]]$input$years_full<2032),]

#lag is 0 so years match

temp$lo1 <- temp$est1 + qnorm(0.025)*temp$se1
temp$lo2 <- temp$est2 + qnorm(0.025)*temp$se2
temp$lo3 <- temp$est3 + qnorm(0.025)*temp$se3
temp$hi1 <- temp$est1 + qnorm(0.975)*temp$se1
temp$hi2 <- temp$est2 + qnorm(0.975)*temp$se2
temp$hi3 <- temp$est3 + qnorm(0.975)*temp$se3

temp$loR1 <- exp(log(temp$R1) + qnorm(0.025)*temp$R1.cv)
temp$loR2 <- exp(log(temp$R2) + qnorm(0.025)*temp$R2.cv)
temp$loR3 <- exp(log(temp$R3) + qnorm(0.025)*temp$R3.cv)
temp$hiR1 <- exp(log(temp$R1) + qnorm(0.975)*temp$R1.cv)
temp$hiR2 <- exp(log(temp$R2) + qnorm(0.975)*temp$R2.cv)
temp$hiR3 <- exp(log(temp$R3) + qnorm(0.975)*temp$R3.cv)

temp$lopredR1 <- exp(log(temp$predR1) + qnorm(0.025)*temp$predR1.cv)
temp$lopredR2 <- exp(log(temp$predR2) + qnorm(0.025)*temp$predR2.cv)
temp$lopredR3 <- exp(log(temp$predR3) + qnorm(0.025)*temp$predR3.cv)
temp$hipredR1 <- exp(log(temp$predR1) + qnorm(0.975)*temp$predR1.cv)
temp$hipredR2 <- exp(log(temp$predR2) + qnorm(0.975)*temp$predR2.cv)
temp$hipredR3 <- exp(log(temp$predR3) + qnorm(0.975)*temp$predR3.cv)

cairo_pdf(file.path("paper", "R_proj_results.pdf"), width = 16, height = 12)
df <- subset(temp, year %in% 2019:2031)
cols <- viridis::viridis_pal(begin = 0.2, end = 0.8, option = "H")(3)
poly_cols <- viridis::viridis_pal(begin = 0.2, end = 0.8, option = "H", alpha = 0.3)(3)
par(mfcol = c(2,2), oma = c(5,1,1,1), mar = c(1,6,1,1))
ylim <- range(subset(temp, year %in% 2019:2031)[c("lo1","lo2","lo3", "hi1", "hi2", "hi3")])
plot(temp$year, temp$est1, type = "n", xlim = c(2018,2031), lwd = 2, ylab = "", ylim = ylim, xaxt = "n", cex.axis = 2, cex.lab = 2, las = 1)
mtext(side = 2, outer = FALSE, "Bottom Temperature Anomaly", cex = 2, line = 4)
axis(1, labels = FALSE)
grid(lty =2, col = gray(0.7))
lines(temp$year, temp$est1, lwd = 2, col = cols[1])
proj_ind <- which(df$year>2020)
lines(df$year[proj_ind], df$est2[proj_ind], col = cols[2], lwd = 2)
lines(df$year[proj_ind], df$est3[proj_ind], col = cols[3], lwd = 2)
polygon(c(temp$year,rev(temp$year)), c(temp$lo1,rev(temp$hi1)), col = poly_cols[1], border = "transparent")
proj_ind <- which(df$year>2021)
polygon(c(df$year[proj_ind],rev(df$year[proj_ind])), c(df$lo2[proj_ind],rev(df$hi2[proj_ind])), col = poly_cols[2], border = "transparent")
polygon(c(df$year[proj_ind],rev(df$year[proj_ind])), c(df$lo3[proj_ind],rev(df$hi3[proj_ind])), col = poly_cols[3], border = "transparent")
abline(v=2021, lty = 3, lwd = 2)

ylim <- range(subset(temp, year %in% 2019:2031)[c("loR1","loR2","loR3", "hiR1", "hiR2", "hiR3")])
plot(temp$year, temp$R1, type = "n", xlim = c(2018,2031), lwd = 2, xlab = "", ylab = "", ylim = ylim, cex.axis = 2, cex.lab = 2, las = 1)
mtext(side = 2, outer = FALSE, expression(Recruitment~~(10^6)), cex = 2, line = 4)
grid(lty =2, col = gray(0.7))
lines(temp$year, temp$R1, lwd = 2, col = cols[1])
lines(temp$year, temp$predR1, lty = 2, lwd = 2, col = cols[1])
proj_ind <- which(df$year>2020)
lines(df$year[proj_ind], df$R2[proj_ind], col = cols[2], lwd = 2)
lines(df$year[proj_ind], df$predR2[proj_ind], col = cols[2], lty = 2, lwd = 2)
lines(df$year[proj_ind], df$R3[proj_ind], col = cols[3], lwd = 2)
lines(df$year[proj_ind], df$predR3[proj_ind], col = cols[3], lty = 2, lwd = 2)
polygon(c(temp$year,rev(temp$year)), c(temp$loR1,rev(temp$hiR1)), col = poly_cols[1], border = "transparent")
proj_ind <- which(df$year>2021)
polygon(c(df$year[proj_ind],rev(df$year[proj_ind])), c(df$loR2[proj_ind],rev(df$hiR2[proj_ind])), col = poly_cols[2], border = "transparent")
polygon(c(df$year[proj_ind],rev(df$year[proj_ind])), c(df$loR3[proj_ind],rev(df$hiR3[proj_ind])), col = poly_cols[3], border = "transparent")
abline(v=2021, lty = 3, lwd = 2)

ylim <- c(0,max(subset(temp, year %in% 2019:2031)[c("R1.cv","R2.cv","R3.cv")]))
plot(temp$year, temp$R1.cv, type = "n", xlim = c(2018,2031), lwd = 2, xlab = "", ylab = "", ylim = ylim, xaxt = "n", cex.axis = 2, cex.lab = 2, las = 1)
mtext(side = 2, outer = FALSE, "Random Effect CV", cex = 2, line = 4)
axis(1, labels = FALSE)
grid(lty =2, col = gray(0.7))
lines(temp$year, temp$R1.cv, lwd = 2, col = cols[1])
proj_ind <- which(df$year>2020)
lines(df$year[proj_ind], df$R2.cv[proj_ind], col = cols[2], lwd = 2)
lines(df$year[proj_ind], df$R3.cv[proj_ind], col = cols[3], lwd = 2)
abline(v=2021, lty = 3, lwd = 2)

ylim <- c(0,max(subset(temp, year %in% 2019:2031)[c("predR1.cv","predR2.cv","predR3.cv")]))
plot(temp$year, temp$predR1.cv, type = "n", xlim = c(2018,2031), lwd = 2, xlab = "", ylab = "", ylim = ylim, cex.axis = 2, cex.lab = 2, las = 1)
mtext(side = 2, outer = FALSE, "Expected CV", cex = 2, line = 4)
grid(lty =2, col = gray(0.7))
lines(temp$year, temp$predR1.cv, lwd = 2, col = cols[1])
proj_ind <- which(df$year>2020)
lines(df$year[proj_ind], df$predR2.cv[proj_ind], col = cols[2], lwd = 2)
lines(df$year[proj_ind], df$predR3.cv[proj_ind], col = cols[3], lwd = 2)
abline(v=2021, lty = 3, lwd = 2)

mtext(side = 1, line = 2, "Year", cex = 2, outer = TRUE)
dev.off()

brp_results <- lapply(c(proj_R1,proj_R3), get.brp.status.results)
# x <- brp_results[[1]]
# c("SSB_status", "F_status")[match(x$df.status$type, c("SSB","Fbar"))]
# x$df.status

df_annual <- rbind(df_annual, cbind.data.frame(Year = c(sapply(brp_results, \(x) c(x$df.status$Year))),
  val = exp(c(sapply(brp_results, \(x) x$df.status$val))),
  cv = c(sapply(brp_results, \(x) x$df.status$se)),
  name = c(sapply(brp_results, \(x) c("SSB_status", "F_status")[match(x$df.status$type, c("SSB","Fbar"))])),
  Region = c(sapply(brp_results, \(x) x$df.status$region)),
  Proj_type = rep(c("AR1", "Last", "Linear"), each = NROW(brp_results[[1]]$df.status)), 
  XSPR_R_type = rep(c("Random Effect", "Expected"), each = 3*NROW(brp_results[[1]]$df.status))))

df_annual$lo <- df_annual$val*exp(-qnorm(0.975)*df_annual$cv)
df_annual$hi <- df_annual$val*exp(qnorm(0.975)*df_annual$cv)

est <-TMB:::as.list.sdreport(proj[[1]]$sdrep, report = TRUE, what = "Est")
se <-TMB:::as.list.sdreport(proj[[1]]$sdrep, report = TRUE, what = "Std.")
proj_df <- cbind.data.frame(year = proj[[1]]$years_full, 
  val = exp(c(est$log_SSB, est$log_SSB_all)), cv = c(se$log_SSB, se$log_SSB_all),
  region = rep(c("North", "South", "Total"), each = length(proj[[1]]$years_full)),
  proj_type = "Continue")
est <-TMB:::as.list.sdreport(proj[[2]]$sdrep, report = TRUE, what = "Est")
se <-TMB:::as.list.sdreport(proj[[2]]$sdrep, report = TRUE, what = "Std.")
proj_df <- rbind(proj_df, cbind.data.frame(year = proj[[2]]$years_full, 
  val = exp(c(est$log_SSB, est$log_SSB_all)), cv = c(se$log_SSB, se$log_SSB_all),
  region = rep(c("North", "South", "Total"), each = length(proj[[2]]$years_full)),
  proj_type = "Average"))
est <-TMB:::as.list.sdreport(proj[[3]]$sdrep, report = TRUE, what = "Est")
se <-TMB:::as.list.sdreport(proj[[3]]$sdrep, report = TRUE, what = "Std.")
proj_df <- rbind(proj_df, cbind.data.frame(year = proj[[3]]$years_full, 
  val = exp(c(est$log_SSB, est$log_SSB_all)), cv = c(se$log_SSB, se$log_SSB_all),
  region = rep(c("North", "South", "Total"), each = length(proj[[3]]$years_full)),
  proj_type = "Trend"))

proj_df$lo <- proj_df$val*exp(qnorm(0.025)*proj_df$cv)
proj_df$hi <- proj_df$val*exp(qnorm(0.025)*proj_df$cv)

this_df <- subset(proj_df, year %in% 2019:2031)
ggplot(this_df, aes(x = year, y = val)) + geom_line(aes(colour = proj_type)) + facet_grid(~region)


temp$SSB_N <- proj[[1]]$rep$NAA[1,1,,1]/1000
temp$predR1 <- proj[[1]]$rep$pred_NAA[1,1,,1]/1000
temp$R1.cv <- TMB:::as.list.sdreport(proj[[1]]$sdrep, report = TRUE, what = "Std")$log_NAA_rep[1,1,,1]

temp$SSB <- 
#plot how R and pred R come together in projection period
#plot different projections of pred R
#plot time series of pred R and Ecov with color gradient f(Ecov)
