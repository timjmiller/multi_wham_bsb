library("wham", lib.loc = "c:/work/wham/old_packages/53e236b")
library(ggplot2)
library(dplyr)
library(patchwork)
source(here::here("R","plot_functions.R"))
theme_set(theme_bw())
theme_update(strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = rel(2)), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)))

proj_R1 <- lapply(1:3, \(x) readRDS(here::here("results",paste0("m1_proj_",x,"_R_opt1.RDS"))))
proj_R3 <- lapply(1:3, \(x) readRDS(here::here("results",paste0("m1_proj_",x,"_R_opt3.RDS"))))
se_R1 <- lapply(proj_R1, \(x) TMB:::as.list.sdreport(x$sdrep, what = "Std", report = TRUE))
se_R3 <- lapply(proj_R3, \(x) TMB:::as.list.sdreport(x$sdrep, what = "Std", report = TRUE))

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
  Proj_type = rep(c("AR1", "Last", "Linear"), each = 3*ny_full), 
  XSPR_R_type = rep(c("Random Effect", "Expected"), each = 3*3*ny_full))
df_annual <- rbind(df_annual, cbind.data.frame(Year = proj_R3[[1]]$years_full, 
  val = exp(c(sapply(c(proj_R1,proj_R3), \(x) c(x$rep$log_Fbar_XSPR[,5:7])))),
  cv = c(sapply(c(se_R1,se_R3), \(x) c(x$log_Fbar_XSPR[,5:7]))),
  name = "F40",
  Region = rep(c("North","South","Total"), each = ny_full),
  Proj_type = rep(c("AR1", "Last", "Linear"), each = 3*ny_full), 
  XSPR_R_type = rep(c("Random Effect", "Expected"), each = 3*3*ny_full)))
df_annual <- rbind(df_annual, cbind.data.frame(Year = proj_R3[[1]]$years_full, 
  val = c(sapply(c(proj_R1,proj_R3), \(x) c(x$rep$NAA[1,1,,1],x$rep$NAA[2,2,,1]))),
  cv = c(sapply(c(se_R1,se_R3), \(x) c(x$log_NAA_rep[1,1,,1],x$log_NAA_rep[2,2,,1]))),
  name = "Random Effect",
  Region = rep(c("North","South"), each = ny_full),
  Proj_type = rep(c("AR1", "Last", "Linear"), each = 2*ny_full), 
  XSPR_R_type = rep(c("Random Effect", "Expected"), each = 2*3*ny_full)))

df_annual <- rbind(df_annual, cbind.data.frame(Year = proj_R3[[1]]$years_full, 
  val = temp[,1],
  cv = temp[,2],
  name = "Expected",
  Region = rep(c("North","South"), each = ny_full),
  Proj_type = rep(c("AR1", "Last", "Linear"), each = 2*ny_full), 
  XSPR_R_type = rep(c("Random Effect", "Expected"), each = 2*3*ny_full)))

df_annual <- rbind(df_annual, cbind.data.frame(Year = c(sapply(brp_results, \(x) c(x$df.status$Year))),
  val = exp(c(sapply(brp_results, \(x) x$df.status$val))),
  cv = c(sapply(brp_results, \(x) x$df.status$se)),
  name = c(sapply(brp_results, \(x) c("SSB_status", "F_status")[match(x$df.status$type, c("SSB","Fbar"))])),
  Region = c(sapply(brp_results, \(x) x$df.status$region)),
  Proj_type = rep(c("AR1", "Last", "Linear"), each = NROW(brp_results[[1]]$df.status)), 
  XSPR_R_type = rep(c("Random Effect", "Expected"), each = 3*NROW(brp_results[[1]]$df.status))))

df_annual <- rbind(df_annual, cbind.data.frame(Year = proj_R3[[1]]$years_full, 
  val = c(sapply(c(proj_R1,proj_R3), \(x) c(x$rep$Fbar[,5:7]))),
  cv = c(sapply(c(se_R1,se_R3), \(x) c(x$log_Fbar[,5:7]))),
  name = "Fbar",
  Region = rep(c("North","South","Total"), each = ny_full),
  Proj_type = rep(c("AR1", "Last", "Linear"), each = 3*ny_full), 
  XSPR_R_type = rep(c("Random Effect", "Expected"), each = 3*3*ny_full)))

df_annual <- rbind(df_annual, cbind.data.frame(Year = proj_R3[[1]]$years_full, 
  val = c(sapply(c(proj_R1,proj_R3), \(x) c(x$rep$SSB, apply(x$rep$SSB,1,sum)))),
  cv = c(sapply(c(se_R1,se_R3), \(x) c(x$log_SSB,x$log_SSB_all))),
  name = "SSB",
  Region = rep(c("North","South","Total"), each = ny_full),
  Proj_type = rep(c("AR1", "Last", "Linear"), each = 3*ny_full), 
  XSPR_R_type = rep(c("Random Effect", "Expected"), each = 3*3*ny_full)))


df_annual$lo <- df_annual$val*exp(-qnorm(0.975)*df_annual$cv)
df_annual$hi <- df_annual$val*exp(qnorm(0.975)*df_annual$cv)

this_df <- subset(df_annual, Proj_type == "Last" & Year < 2022)
# cnms <- c(North = "", South = "", Total = "Total")
SSB40_df <- subset(this_df, name == "SSB(F40)")
SSB40_df$val <- SSB40_df$val/1000
SSB40_df$hi <- SSB40_df$hi/1000
SSB40_df$lo <- SSB40_df$lo/1000
plt_SSB40 <- ggplot(SSB40_df, aes(x = Year, y = val, colour = XSPR_R_type)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") +  
    labs(color = "Recruitment Type", fill = "Recruitment Type") +
    geom_point(size = 2) + geom_line()  + ylab(expression(SSB["40%"]~"(kmt)")) + xlab("Year") +
    facet_grid(~ Region) + coord_cartesian(ylim = c(0,6e1)) +#ylim(0,1e5) + #xlim(1989,2021)
    geom_ribbon(ggplot2::aes(ymin=lo, ymax=hi, fill=XSPR_R_type), alpha=0.3, linetype = 0)
plt_SSB40

cnms <- c(North = "", South = "", Total = "")
plt_SSB_status <- ggplot(subset(this_df, name == "SSB_status"), aes(x = Year, y = val, colour = XSPR_R_type)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") +  
    labs(color = "Recruitment Type", fill = "Recruitment Type") +
    geom_point(size = 2) + geom_line()  + ylab(bquote(SSB*"/"*SSB[40*"%"])) + xlab("Year") +
    facet_grid(~ Region, labeller = labeller(Region = cnms)) + #coord_cartesian(ylim = c(0,6e4)) +#ylim(0,1e5) + #xlim(1989,2021)
    geom_ribbon(ggplot2::aes(ymin=lo, ymax=hi, fill=XSPR_R_type), alpha=0.3, linetype = 0)
plt_SSB_status

plt_F40 <- ggplot(subset(this_df, name == "F40"), aes(x = Year, y = val, colour = XSPR_R_type)) + 
    # scale_colour_viridis_d(begin = 0.2, end = 0.8) + scale_fill_viridis_d(begin = 0.2, end = 0.8) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + 
    labs(color = "Recruitment Type", fill = "Recruitment Type") +
    geom_point(size = 2) + geom_line()  + ylab(expression(bar(italic(F))[ "40%"]~~"(Ages 6-7)")) + xlab("Year") +
    facet_grid(~ Region, labeller = labeller(Region = cnms)) + 
    geom_ribbon(ggplot2::aes(ymin=lo, ymax=hi, fill=XSPR_R_type), alpha=0.3, linetype = 0)
plt_F40

plt_F_status <- ggplot(subset(this_df, name == "F_status"), aes(x = Year, y = val, colour = XSPR_R_type)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") +  
    labs(color = "Recruitment Type", fill = "Recruitment Type") +
    geom_point(size = 2) + geom_line()  + ylab(bquote(bar(italic(F))*"/"*bar(italic(F))[40*"%"]~~"(Ages 6-7)")) + xlab("Year") +
    facet_grid(~ Region, labeller = labeller(Region = cnms)) + 
    geom_ribbon(ggplot2::aes(ymin=lo, ymax=hi, fill=XSPR_R_type), alpha=0.3, linetype = 0)
plt_F_status


brp_results_static <- lapply(c(proj_R1[1],proj_R3[1]), get.brp.status.results, static = TRUE)
kobe_results_static <- lapply(c(proj_R1[1],proj_R3[1]), kobe.plot, static = TRUE, status.years = length(proj_R1[[1]]$years))
source(here::here("R","plot_functions.R"))
poly.colors <- c(kobe_results_static[[1]]$poly.colors, ellipse = "transparent", center = "transparent")
df <- rbind(subset(brp_results_static[[1]]$df.ellipse, Year == 2021),
	subset(brp_results_static[[2]]$df.ellipse, Year == 2021))
df$XSPR_R_type <- rep(c("Random Effect", "Expected"), each = NROW(kobe_results_static[[1]]$df.ellipse))
df$x <- exp(df$x)
df$y <- exp(df$y)
borders <- kobe_results_static[[1]]$borders
df <- rbind(df, cbind.data.frame(Year = 2021, region = rep(c("North","South","Total"), each = NROW(borders)),
	ptype = borders$region, x = borders$x, y = borders$y, XSPR_R_type = NA))
df$ptype <- factor(df$ptype)
df$region <- factor(df$region)
df.ellipse <- subset(df,ptype == "ellipse")
df.center <- subset(df,ptype == "center")
df.colors <- subset(df,!ptype %in% c("ellipse","center"))
df.prob <- cbind.data.frame(Year = 2021, 
	probs = paste0("P = ", round(c(sapply(kobe_results_static, \(z) sapply(z$p.status, \(y) unlist(y[[1]])))),2)),
	region = rep(c("North","South","Total"), 4), x = rep(c(0,0,max(df.ellipse$x),max(df.ellipse$x)), each = 3),
	y = c(rep(c(0.1,max(df.ellipse$y),0.1, max(df.ellipse$y)), each = 3), rep(c(0,max(df.ellipse$y)-0.1,0, max(df.ellipse$y)-0.1), each = 3)))
df.prob$XSPR_R_type <- rep(c("Random Effect", "Expected"), each = 4*3)
plt_kobe <- ggplot(df, aes(x = x, y = y)) + facet_grid(~ region, labeller = labeller(region = cnms)) + 
    geom_polygon(df.colors, mapping = aes(x = x, y = y, fill = ptype)) + coord_cartesian(ylim = c(0, max(df.ellipse$y)), xlim = c(0,max(df.ellipse$x))) + 
    geom_polygon(subset(df.ellipse, XSPR_R_type == "Random Effect"), mapping = aes(x = x, y = y), fill = NA, colour = alpha("black", 0.4), linewidth = 1.5) + 
    geom_polygon(subset(df.ellipse, XSPR_R_type == "Expected"), mapping = aes(x = x, y = y), fill = NA, colour = alpha("gray48", 0.4), linewidth = 1.5) + 
    ylab(bquote(bar(italic(F))*"/"*bar(italic(F))[40*"%"]~~"(Ages 6-7)")) + xlab(bquote(SSB*"/"*SSB[40*"%"])) +
    geom_text(subset(df.center, XSPR_R_type == "Random Effect"), mapping  = aes(x = x , y = y, label = Year)) + 
    geom_text(subset(df.center, XSPR_R_type == "Expected"), color = "gray48", mapping  = aes(x = x , y = y, label = Year)) + 
    geom_text(subset(df.prob, XSPR_R_type == "Random Effect"), mapping  = aes(x = x , y = y, label = probs), hjust = "inward") + 
    geom_text(subset(df.prob, XSPR_R_type == "Expected"), mapping  = aes(x = x , y = y, label = probs), color = "gray48", hjust = "inward") + 
    scale_fill_manual(values = poly.colors) + theme(legend.position="none")
plt_kobe

theme_update(strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = rel(2)), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)))

#Fig 6
cairo_pdf(file.path("paper", "brp_status_results.pdf"), width = 15, height = 20)
design <- c(area(1,1,1,1), area(2,1,2,1), area(3,1,3,1), area(4,1,4,1), area(5,1,5,1))
	(plt_SSB40 +  xlab("")  + theme(legend.position="none", axis.title.x = element_blank(), strip.text=element_text(margin=margin(b = 5)), axis.text.x=element_blank(), plot.margin = margin(b = 1, t = 0))) + 
  (plt_SSB_status + xlab("")  + theme(legend.position="none", strip.clip = "off", strip.text=element_text(margin=margin(b = 5, t=-20)), plot.margin = margin(b = 1, t = 0), axis.title.x = element_blank(), axis.text.x=element_blank())) + 
  (plt_F40 + xlab("") + theme(legend.position="none", strip.clip = "off", strip.text = element_text(margin=margin(b = 5, t=-15)), plot.margin = margin(b = 1, t = 0), axis.title.x = element_blank(), axis.text.x=element_blank())) + 
  (plt_F_status + theme(legend.position="none", strip.clip = "off", strip.text = element_text(margin=margin(b = 5, t=-15)), plot.margin = margin(b = 1, t = 0))) + 
  (plt_kobe + theme(legend.position="none", strip.clip = "off", strip.text = element_text(margin=margin(b = 5, t=-15)), plot.margin = margin(b = 1, t = 0))) + 
  plot_layout(design = design)
dev.off()

#######################################
#CVs

plt_R_cv <- ggplot(subset(this_df, XSPR_R_type == "Random Effect" & name %in% c("Random Effect", "Expected")), aes(x = Year, y = cv, colour = name)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + 
    labs(color = "Recruitment Type") + geom_point(size = 2) + geom_line()  + ylab("CV(Recruitment)") + xlab("Year") + facet_grid(~ Region)
plt_R_cv

cnms <- c(North = "", South = "", Total = "Total")
plt_SSB40_cv <- ggplot(subset(this_df, name == "SSB(F40)"), aes(x = Year, y = cv, colour = XSPR_R_type)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + 
    labs(color = "Recruitment Type") + geom_point(size = 2) + geom_line()  + ylab(expression(CV(SSB(italic(F)[ "40%"])))) + xlab("Year") +
    facet_grid(~ Region, labeller = labeller(Region = cnms)) 
plt_SSB40_cv

cnms <- c(North = "", South = "", Total = "")
plt_F40_cv <- ggplot(subset(this_df, name == "F40"), aes(x = Year, y = cv, colour = XSPR_R_type)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + 
    labs(color = "Recruitment Type") + geom_point(size = 2) + geom_line()  + ylab(expression(CV(bar(italic(F))[ "40%"]))) + xlab("Year") +
    facet_grid(~ Region, labeller = labeller(Region = cnms))
plt_F40_cv


cairo_pdf(file.path("paper", "brp_cv_results.pdf"), width = 20, height = 12)
design <- c(area(1,1,1,778), area(2,1,2,1000), area(3,1,3,1000))
(plt_R_cv +  xlab("")  + theme(axis.title.x = element_blank(), strip.text=element_text(margin=margin(b = 5)), axis.text.x=element_blank(), plot.margin = margin(b = 1, t = 0))) + 
  (plt_SSB40_cv + xlab("")  + theme(strip.clip = "off", strip.text=element_text(margin=margin(b = 5, t=-20)), axis.title.x = element_blank(), legend.position="none",axis.text.x=element_blank(), plot.margin = margin(t = 0, b = 1))) + 
  (plt_F40_cv + theme(legend.position="none", strip.clip = "off", strip.text = element_text(margin=margin(t=-15)), plot.margin = margin(b = 1, t = 0))) + 
  plot_layout(design = design)
dev.off()

