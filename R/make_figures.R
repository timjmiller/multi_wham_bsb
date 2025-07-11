# load the appropriate R package libraries
library(rnaturalearth)
library(maps)
library(sp)
library(sf)
library("wham", lib.loc = "c:/work/wham/old_packages/53e236b")
# library(wham)
library(ggplot2)
library(dplyr)
library(patchwork)
# load some helper functions
source(here::here("R","plot_functions.R"))

#load fits
fits <- readRDS(file.path("results","fits_no_M_re_rev_1.RDS"))
fits_M_re <- readRDS(file.path("results","fits_M_re_best.RDS"))

##############################################################################################


#Fig 1.
#hudson
#hudson.n_UTM <- data.frame(x=c(357212.46,19314.01,658471.31,935303.78),y=c(4261346,4566676,5099985,4635884))
spRef_UTM_19 <- "+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
spRef_DD <- "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"
hudson.n_UTM <- data.frame(x=c(357212.46,19314.01),y=c(4261346,4566676))
hudson.n_UTM.rl <- rbind(hudson.n_UTM, hudson.n_UTM[1,])
hudson.n_UTM.l <- Line(hudson.n_UTM.rl)
hudson.n_UTM.l <- Lines(list(hudson.n_UTM.l), ID = 'n')
hudson.n_UTM.spL <- SpatialLines(list(hudson.n_UTM.l))
proj4string(hudson.n_UTM.spL) <- CRS(spRef_UTM_19)
hudson.n <- spTransform(hudson.n_UTM.spL, CRS(spRef_DD))
hudson.n.sf <- st_as_sf(hudson.n,coords = c("x", "y"))
world <- st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))
states <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))

# contours <- sf::st_as_sf(sf::read_sf(("c:/work/shapefiles/Contours"), quiet=TRUE))
# contour_400m <- subset(contours, elev_m == -400)
# saveRDS(contour_400m, here::here("data","contour_400m.RDS"))
contour_400m <- readRDS(here::here("data","contour_400m.RDS"))
lat.range <- c(34.5,42)
lon.range <- c(-78,-70)
lat.range <- c(35,45)
lon.range <- c(-76,-66)

state_prov <- rnaturalearth::ne_states(c("united states of america", "canada"))

# ggmap <- ggplot(states) +
ggmap <- ggplot(state_prov) +
  geom_sf(alpha = 0.5) +
  # geom_sf(data = states, alpha = 0.5) +
  geom_sf(data = contour_400m, linewidth = 0.1) +
  geom_sf(data = hudson.n.sf, linewidth = 1, colour = "red") +
  xlim(lon.range) +
  ylim(lat.range) +
  theme_bw() +
  annotate("text", x = -70, y = 37, label = "Atlantic Ocean") +
  annotate("text", x = -68, y = 43, label = "Gulf\nof\nMaine") +
  theme(legend.title=element_blank(), axis.title = element_blank())

cairo_pdf(file.path("paper", "map.pdf"), width = 7, height = 7)
ggmap
dev.off()

##############################################################################################

#Fig. 2 from power point 
##############################################################################################

#Fig. 3
#prior/posterior for movement parameter
dlogitnorm <- function(x, p, sigma){

  y <- log(x/(1-x))
  mu <- log(p/(1-p))
  fx <- dnorm(y, mu, sigma, log = F) /(x*(1-x)) 
  return(fx)
}

fit1 <- readRDS(here::here("results","fit_1.RDS")) #M_1

cairo_pdf(file.path("paper", "move_prior_post.pdf"), width = 7, height = 11)
par(mfrow = c(2,1), mar = c(4,4,1,1), oma = c(1,1,1,1))
ps <- seq(0,0.1,0.001)
#prior for both north-south and south-north
plot(ps, dlogitnorm(ps,0.02214863, 0.2), type = 'l', lwd = 2, xlab = expression(italic("\u03BC")["N" %->% phantom(0)* "S"]), ylab = expression(f(italic("\u03BC")["N " %->% phantom(0)* "S"])), ylim = c(0,300))
abline(v=0.02214863, lwd =2)
#posterior for north to south
lines(ps, dlogitnorm(ps, p = fit1$rep$mu[1,8,8,1,1,2], sigma = TMB:::as.list.sdreport(fit1$sdrep, "Std")$mu_prior_re[1,8,1,1]), lwd = 2, col = "red")
abline(v=fit1$rep$mu[1,8,8,1,1,2], col = "red", lwd = 2)

ps <- seq(0,0.6,0.001)
plot(ps, dlogitnorm(ps,0.3130358, 0.2), type = 'l', lwd = 2, xlab = expression(italic("\u03BC")["S" %->% phantom(0)* "N"]), ylab = expression(f(italic("\u03BC")["S " %->% phantom(0)* "N"])))
abline(v=0.3130358, lwd =2)
#posterior for north to south
lines(ps, dlogitnorm(ps, p = fit1$rep$mu[1,8,1,1,2,1], sigma = TMB:::as.list.sdreport(fit1$sdrep, "Std")$mu_prior_re[1,8,2,1]), lwd = 2, col = "red")
abline(v=fit1$rep$mu[1,8,1,1,2,1], col = "red", lwd = 2)
dev.off()

sapply(c(fits,fits_M_re), \(x) x$rep$mu[1,8,8,1,1,2])
sapply(c(fits,fits_M_re), \(x) x$rep$mu[1,8,1,1,2,1])


################################################################################

#Fig. 5

#plot selectivity at time and age for blocks with RE.
source(here::here("R","plot_functions.R"))
cairo_pdf(file.path("paper", "selectivity_re_plot.pdf"), width = 16, height = 12)
plot.selectivity(fit1, blocks = c(1:2,5:6), block.names = c("North Commercial Fleet", "North Recreational Fleet", "North Recreational Catch/Effort Index", "North Spring Index"))
dev.off()


################################################################################



#BT for north and south on same plot


cairo_pdf(file.path("paper", paste0("Ecov.pdf")), height = 10, width = 10)
plot.ecov(fit1)
legend("topleft", legend = c("North", "South"), col = mypalette(2), lty =  1, pch = 19)
dev.off()

library(latex2exp) #https://github.com/stefano-meschiari/latex2exp
################################################################################
#Fig. S7
cairo_pdf(file.path("paper", paste0("Ecov_M1_rel_M0.pdf")), height = 10, width = 10)
par(mar = c(4,7,1,1), oma = c(0,0,0,0))
plot(fits[[1]]$input$years_Ecov, (fits[[2]]$rep$Ecov_x[,1] - fits[[1]]$rep$Ecov_x[,1])/abs(fits[[1]]$rep$Ecov_x[,1]), xlab = "Year", ylab = TeX("$\\frac{\\widehat{X}\\left(M_1\\right) - \\widehat{X}\\left(M_0\\right)}{|\\widehat{X}\\left(M_0\\right)|}$" ))
grid(col = gray(0.7), lwd = 2, lty = 2)
dev.off()
################################################################################


###############################################
#SSB (and model-weighted SSB?)

summarize_res_fn <- function(fits, rep_name = "log_SSB", index = 1, age = 1){
  out <- list()
  if(rep_name == "log_NAA_rep"){
    out$Est <- sapply(fits, \(x) {
      return(TMB:::as.list.sdreport(x$sdrep, "Est", report = TRUE)[[rep_name]][index,index,,age])
    })
    out$SE <- sapply(fits, \(x) {
      return(TMB:::as.list.sdreport(x$sdrep, "Std", report = TRUE)[[rep_name]][index,index,,age])
    })
  } else {
    out$Est <- sapply(fits, \(x) {
      return(TMB:::as.list.sdreport(x$sdrep, "Est", report = TRUE)[[rep_name]][,index])
    })
    out$SE <- sapply(fits, \(x) {
      return(TMB:::as.list.sdreport(x$sdrep, "Std", report = TRUE)[[rep_name]][,index])
    })
  }
  out$lo <- out$Est + qnorm(0.025) * out$SE
  out$hi <- out$Est + qnorm(0.975) * out$SE
  # print(out$median_res)
  return(out)
}

x <- summarize_res_fn(c(fits, fits_M_re))

x <- summarize_res_fn(c(fits, fits_M_re), rep_name = "log_NAA_rep")

make_df <- function(fits, rep_name = "log_SSB", label = "SSB", age = 1){
  df <- cbind.data.frame(year = numeric(), model = character(), stock = character(), SSB = numeric(), lo = numeric(),hi = numeric())
  indices <- 1:2
  if(rep_name == "log_Fbar") indices <- indices + 4 #4 fleets, 2 regions, 1 total
  print(indices)
  summary_N <- summarize_res_fn(fits, rep_name = rep_name, index = indices[1])
  print(NCOL(summary_N$Est))
  summary_S <- summarize_res_fn(fits, rep_name = rep_name, index = indices[2])
  mods <- paste0("$\\textit{M}_{",1:NCOL(summary_N$Est) - 1,"}$")
  print(mods)
  for(i in 1:NCOL(summary_N$Est)){
    df <- rbind(df, cbind.data.frame(year = fits[[1]]$years, model = mods[i], stock = "North", type = label, est = exp(summary_N$Est[,i]), lo = exp(summary_N$lo[,i]),hi = exp(summary_N$hi[,i])))
    df <- rbind(df, cbind.data.frame(year = fits[[1]]$years, model = mods[i], stock = "South", type = label, est = exp(summary_S$Est[,i]), lo = exp(summary_S$lo[,i]),hi = exp(summary_S$hi[,i])))
  }
  df$model <- factor(df$model, levels = mods)
  return(df)
}

df <- make_df(c(fits,fits_M_re), label = "SSB (mt)")
df <- rbind(df, make_df(c(fits,fits_M_re), rep_name = "log_Fbar", label = "Avg. F (Ages 6-7)"))
temp_df <- make_df(c(fits,fits_M_re), rep_name = "log_NAA_rep", label = "Recruitment (Millions)")
temp_df[,c("est","lo","hi")] <-temp_df[,c("est","lo","hi")]/1000
df <- rbind(df,temp_df)
df_M1 <- subset(df, model == levels(df$model)[2])
plt <- ggplot(subset(df, model == levels(df$model)[2]), aes(x = year, y = est)) + scale_colour_viridis_d() + 
    geom_point(size = 1) + geom_line() +
    facet_grid(type ~ stock, scales = "free_y", switch = "y") + theme_bw() +
    xlab("Year") + ylab(NULL) + 
    theme(strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = rel(2)), 
      axis.title.x = element_text(size = rel(2)), axis.text = element_text(size = rel(2))) +#, strip.text.y = element_blank()) +
    geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.3, linetype = 0)
plt

cairo_pdf(file.path("paper", "M_1_SSB_F_R.pdf"), width = 12, height = 16)
plt
dev.off()

df_rel <- cbind(df, rel_M1_est = NA)
for(i in unique(df$model)) for(j in unique(df$type)) for(k in unique(df$stock)) {
  ind <- which(df_rel$model ==i & df_rel$type == j & df_rel$stock == k)
  ind1 <- which(df_rel$model ==levels(df_rel$model)[2] & df_rel$type == j & df_rel$stock == k)
  df_rel$rel_M1_est[ind] <- df_rel$est[ind]/df_rel$est[ind1]
}
df_rel$type <- factor(df_rel$type)

library(latex2exp) #https://github.com/stefano-meschiari/latex2exp
TeX("SSB/SSB($M_1$)")
type_labels <- c("frac(bar(italic(F)),bar(italic(F))(italic(M)[1]))", "frac(Recruitment,Recruitment(italic(M)[1]))", "frac(SSB,SSB(italic(M)[1]))")

df_rel$type_labels <- type_labels[match(df_rel$type, levels(df_rel$type))]

plt <- ggplot(df_rel, aes(x = year, y = rel_M1_est, colour = model)) + #scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + 
    geom_point(size = 1) + geom_line() +
    facet_grid(type_labels ~ stock, scales = "free_y", switch = "y", labeller = label_parsed) + theme_bw() +
    xlab("Year") + ylab(NULL) + labs(colour = "Model") +
    scale_color_viridis_d(begin = 0.2, end = 0.8, option = "turbo", labels=lapply(levels(df$model), TeX)) + 
    theme(strip.background = element_blank(), strip.placement = "outside", legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), strip.text = element_text(size = rel(2)), 
      axis.title.x = element_text(size = rel(2)))
plt

###############################################
#Fig. S6
cairo_pdf(file.path("paper", "SSB_F_R_rel_M1.pdf"), width = 12, height = 16)
plt
dev.off()

###############################################


#reference points
#R1 = annual NAA, R3 = annual pred_NAA

proj_R1 <- lapply(1:3, \(x) readRDS(here::here("results",paste0("m1_proj_",x,"_R_opt1.RDS"))))
proj_R3 <- lapply(1:3, \(x) readRDS(here::here("results",paste0("m1_proj_",x,"_R_opt3.RDS"))))

se_R1 <- lapply(proj_R1, \(x) TMB:::as.list.sdreport(x$sdrep, what = "Std", report = TRUE))
se_R3 <- lapply(proj_R3, \(x) TMB:::as.list.sdreport(x$sdrep, what = "Std", report = TRUE))
ny_full <- length(proj_R3[[1]]$years_full)
df_brps <- cbind.data.frame(Year = proj_R3[[1]]$years_full, 
  val = exp(c(sapply(c(proj_R1,proj_R3), \(x) c(x$rep$log_SSB_FXSPR)))),
  cv = c(sapply(c(se_R1,se_R3), \(x) c(x$log_SSB_FXSPR))),
  # val = exp(c(sapply(paste0("proj_",1:3,"_R", rep(c(1,3),each = 3)), \(x) c(get(x)$rep$log_SSB_FXSPR)))),
  name = "SSB(F40)",
  Region = rep(c("North","South","Total"), each = ny_full),
  Proj_type = rep(c("AR1", "Last", "Linear"), each = 3*ny_full), 
  XSPR_R_type = rep(c("Random Effect", "Expected"), each = 3*3*ny_full))
df_brps <- rbind(df_brps, cbind.data.frame(Year = proj_R3[[1]]$years_full, 
  val = exp(c(sapply(c(proj_R1,proj_R3), \(x) c(x$rep$log_Fbar_XSPR[,5:7])))),
  cv = c(sapply(c(se_R1,se_R3), \(x) c(x$log_Fbar_XSPR[,5:7]))),
  name = "F40",
  Region = rep(c("North","South","Total"), each = ny_full),
  Proj_type = rep(c("AR1", "Last", "Linear"), each = 3*ny_full), 
  XSPR_R_type = rep(c("Random Effect", "Expected"), each = 3*3*ny_full)))
df_brps <- rbind(df_brps, cbind.data.frame(Year = proj_R3[[1]]$years_full, 
  val = c(sapply(c(proj_R1,proj_R3), \(x) c(x$rep$NAA[1,1,,1],x$rep$NAA[2,2,,1]))),
  cv = c(sapply(c(se_R1,se_R3), \(x) c(x$log_NAA_rep[1,1,,1],x$log_NAA_rep[2,2,,1]))),
  name = "Random Effect",
  Region = rep(c("North","South"), each = ny_full),
  Proj_type = rep(c("AR1", "Last", "Linear"), each = 2*ny_full), 
  XSPR_R_type = rep(c("Random Effect", "Expected"), each = 2*3*ny_full)))

dim(subset(x, Proj_type == "Last" & Region == "North" & Year < 2022))

# fjp <- list() 
# for(i in 1:3) fjp[[i]] <- readRDS(here::here("results",paste0("m1_proj_",i,"_sdrep_fjp.RDS")))


# #have to do se(predR) by hand
# source(here::here("R","plot_functions.R"))
# temp <- make_pred_R_cv(proj_R3[[1]], fjp[[1]], type = 1)
# temp <- rbind(temp, make_pred_R_cv(proj_R3[[1]], fjp[[1]], type = 1, region = 2))
# temp <- rbind(temp,make_pred_R_cv(proj_R3[[2]], fjp[[2]], type = 2))
# temp <- rbind(temp,make_pred_R_cv(proj_R3[[2]], fjp[[2]], type = 2, region = 2))
# temp <- rbind(temp,make_pred_R_cv(proj_R3[[3]], fjp[[3]], type = 3))
# temp <- rbind(temp,make_pred_R_cv(proj_R3[[3]], fjp[[3]], type = 3, region = 2))

temp <- make_pred_R_cv(proj_R3[[1]], type = 1)
temp <- rbind(temp, make_pred_R_cv(proj_R3[[1]], type = 1, region = 2))
temp <- rbind(temp,make_pred_R_cv(proj_R3[[2]], type = 2))
temp <- rbind(temp,make_pred_R_cv(proj_R3[[2]], type = 2, region = 2))
temp <- rbind(temp,make_pred_R_cv(proj_R3[[3]], type = 3))
temp <- rbind(temp,make_pred_R_cv(proj_R3[[3]], type = 3, region = 2))

df_brps <- rbind(df_brps, cbind.data.frame(Year = proj_R3[[1]]$years_full, 
  val = temp[,1],
  cv = temp[,2],
  name = "Expected",
  Region = rep(c("North","South"), each = ny_full),
  Proj_type = rep(c("AR1", "Last", "Linear"), each = 2*ny_full), 
  XSPR_R_type = rep(c("Random Effect", "Expected"), each = 2*3*ny_full)))

df_brps$lo <- df_brps$val*exp(-qnorm(0.975)*df_brps$cv)
df_brps$hi <- df_brps$val*exp(qnorm(0.975)*df_brps$cv)

saveRDS(df_brps, here::here("results", "df_brps_rev.RDS"))


this_df <- subset(df_brps, Proj_type == "Last" & Year < 2022)
library(ggplot2)  
plt_R <- ggplot(subset(this_df, XSPR_R_type == "Random Effect" & name %in% c("Random Effect", "Expected")), aes(x = Year, y = val, colour = name)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") +  
    labs(color = "Recruitment Type", fill = "Recruitment Type") +
    geom_point(size = 2) + geom_line()  + ylab("Recruitment (1000s)") + xlab("Year") +
    facet_grid(~ Region) + coord_cartesian(ylim = c(0,1.25e5)) +#ylim(0,1e5) + #xlim(1989,2021)
    geom_ribbon(ggplot2::aes(ymin=lo, ymax=hi, fill=name), alpha=0.3, linetype = 0)
plt_R

cnms <- c(North = "", South = "", Total = "Total")
plt_SSB40 <- ggplot(subset(this_df, name == "SSB(F40)"), aes(x = Year, y = val, colour = XSPR_R_type)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") +  
    labs(color = "Recruitment Type", fill = "Recruitment Type") +
    geom_point(size = 2) + geom_line()  + ylab(expression(SSB(italic(F)[ "40%"])~"(mt)")) + xlab("Year") +
    facet_grid(~ Region, labeller = labeller(Region = cnms)) + coord_cartesian(ylim = c(0,6e4)) +#ylim(0,1e5) + #xlim(1989,2021)
    geom_ribbon(ggplot2::aes(ymin=lo, ymax=hi, fill=XSPR_R_type), alpha=0.3, linetype = 0)
plt_SSB40

cnms <- c(North = "", South = "", Total = "")
plt_F40 <- ggplot(subset(this_df, name == "F40"), aes(x = Year, y = val, colour = XSPR_R_type)) + 
    # scale_colour_viridis_d(begin = 0.2, end = 0.8) + scale_fill_viridis_d(begin = 0.2, end = 0.8) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + 
    labs(color = "Recruitment Type", fill = "Recruitment Type") +
    geom_point(size = 2) + geom_line()  + ylab(expression(bar(italic(F))[ "40%"])) + xlab("Year") +
    facet_grid(~ Region, labeller = labeller(Region = cnms)) + 
    geom_ribbon(ggplot2::aes(ymin=lo, ymax=hi, fill=XSPR_R_type), alpha=0.3, linetype = 0)
plt_F40

plt_R + plt_SSB40 + plot_layout(ncol = 2)


cairo_pdf(file.path("paper", "brp_results.pdf"), width = 20, height = 12)
design <- c(area(1,1,1,781), area(2,1,2,1000), area(3,1,3,1000))
(plt_R +  xlab("")  + theme(axis.title.x = element_blank(), strip.text=element_text(margin=margin(b = 5)), axis.text.x=element_blank(), plot.margin = margin(b = 1, t = 0))) + 
  (plt_SSB40 + xlab("")  + theme(strip.clip = "off", strip.text=element_text(margin=margin(b = 5, t=-20)), axis.title.x = element_blank(), legend.position="none",axis.text.x=element_blank(), plot.margin = margin(t = 0, b = 1))) + 
  (plt_F40 + theme(legend.position="none", strip.clip = "off", strip.text = element_text(margin=margin(t=-15)), plot.margin = margin(b = 1, t = 0))) + 
  plot_layout(design = design)
dev.off()

# source(here::here("R","plot_functions.R"))
# x <- kobe.plot(proj, status.years = length(proj$years), static = FALSE)
# df <- x$df.ellipse

# source(here::here("R","plot_functions.R"))
# x <- get.status.results(proj)
# y <- get.status.results(proj, static = TRUE)

# ggplot(x$borders, aes(x = x, y = y, fill = region)) + 
#     geom_polygon() + coord_cartesian(ylim = c(0, max(x$borders$y)), xlim = c(0,max(x$borders$x))) + 
#     scale_fill_manual(values = x$poly.colors) + theme(legend.position="none")
# df <- x$df.ellipse
# df$ptype <- factor(df$ptype)
# df$region <- factor(df$region)
#   gg.kobe <- ggplot(df, aes(x = x, y = y)) + facet_grid(~region)+ 
#     geom_polygon(subset(x$df.ellipse, ptype == "ellipse"), aes(x = x, y = y)) + 


###############################################
# Ecov effects on M and R
fit_3 <- readRDS(file.path("results","fits_M_re_best.RDS"))[[5]]
Ecov <- fit_3$rep$Ecov_out_M[1,1,1,,][,1]
M_out <- fit_3$rep$MAA[1,1,,1]
ord <- order(Ecov)

Ecov_plt <- seq(min(Ecov),max(Ecov),0.01)
M <- exp(fit_3$parList$Mpars[1,1,1] + fit_3$parList$Ecov_beta_M[1,1,1,1,1] * Ecov_plt)
pal <- viridis::viridis_pal(option = "H")(length(ord))

df <- cbind.data.frame(M=M, Temperature  = Ecov_plt, group = 1)
df <- rbind(df, cbind.data.frame(M = M_out, Temperature = Ecov, group = 2))
plt <- ggplot(df, aes(Temperature, M, colour = Temperature)) +
  geom_line(data = subset(df, group == 1), linewidth = 1.2) +
  geom_point(data = subset(df, group == 2)) + 
  scale_colour_gradientn(colours = pal)

fit0 <- readRDS(file.path("results","fits_no_M_re_rev_1.RDS"))[[1]]
fit1 <- readRDS(file.path("results","fits_no_M_re_rev_1.RDS"))[[2]]
R_out <- fit1$rep$NAA[1,1,,1]
Ecov <- fit1$rep$Ecov_out_R[1,,1]
ord <- order(Ecov)
Ecov_plt <- seq(min(Ecov),max(Ecov),0.01)
Expected_R <- exp(fit1$parList$mean_rec_pars[1,1] + fit1$parList$Ecov_beta_R[1,1,1] * Ecov_plt)
pal <- viridis::viridis_pal(begin = 0.2, end = 0.8, option = "H")(length(ord))
years <- 1:length(ord)
cbind(exp(fit1$parList$mean_rec_pars[1,1] + fit1$parList$Ecov_beta_R[1,1,1] * Ecov), R_out)

plot(log(fit0$rep$NAA[1,1,-1,1]), log(fit1$rep$NAA[1,1,-1,1]), ylab = "log(Recruitment) with BT Effect", xlab = "log(Recruitment) without BT Effect")
abline(0,1)
plot(fit1$years[-1], fit1$rep$pred_NAA[1,1,-1,1], ylab = "Expected Recruitment (1000s)", xlab = "Year", type = 'l')
lines(fit1$years[-1], fit0$rep$pred_NAA[1,1,-1,1], col = "gray")

df_brps <- readRDS(here::here("results", "df_brps_rev.RDS"))
R_df <- subset(df_brps, Proj_type == "Last" & Year < 2022 & XSPR_R_type == "Random Effect" & name %in% c("Random Effect", "Expected"))
R_df[,"BT"] <- c(fits[[2]]$rep$Ecov_out_R[1,,1],fits[[2]]$rep$Ecov_out_R[2,,2])
R_df_N <- subset(R_df, Region =="North")
R_df_N$name <- factor(R_df_N$name)

df <- cbind.data.frame(Year = fit1$years, Recruitment = fit1$rep$pred_NAA[1,1,,1], BT = Ecov, Type = "Expected")
df <- rbind(df, cbind.data.frame(Year = fit1$years, Recruitment = fit1$rep$NAA[1,1,,1], BT = Ecov, Type = "Random Effect"))
df[,"order_Ecov"] <- ord
df$order_Ecov <- factor(df$order_Ecov)
df$Type <- factor(df$Type)
df$lo <- c(subset(R_df_N, name == "Expected")$lo, subset(R_df_N, name == "Random Effect")$lo)
df$hi <- c(subset(R_df_N, name == "Expected")$hi, subset(R_df_N, name == "Random Effect")$hi)
plt <- ggplot(df, aes(Year, Recruitment, colour = BT, shape = Type)) + theme_bw() + 
  geom_line(data = subset(df,  Type == "Expected"), colour = gray(0.7, alpha = 0.4), linewidth = 2) +
  geom_point(data = subset(df, Type == "Expected"), size = 2) +
  geom_point(data = subset(df, Type == "Random Effect"), size = 2) + 
  ylab("Recruitment (1000s)") + labs(shape="Estimate Type", colour="Bottom Temperature\nAnomaly") + 
  theme(axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.title = element_text(size = rel(2), hjust=0.5), 
    legend.box.just = "center",legend.text = element_text(size = rel(1.5))) +
  scale_colour_gradientn(colours = pal)# + scale_fill_gradientn(colours = pal)
plt

###############################################
#Fig. S5
cairo_pdf(file.path("paper", "best_R_Ecov.pdf"), width = 14, height = 12)
plt
dev.off()
###############################################

# plt <- ggplot(df, aes(Year, Recruitment, colour = order_Ecov, shape = Type, ymin=lo, ymax=hi)) + theme_bw() + 
#   # geom_line(data = subset(df,  Type == "Expected"), colour = gray(0.7, alpha = 0.4), size = 2) +
#   geom_line(aes(Year, Recruitment, colour = Type), size = 1.5, inherit.aes = FALSE) +
#   geom_point(aes(fill=order_Ecov), size=4, shape=21, stroke=0) +
#   scale_fill_gradientn(colours = pal) +
#   ylab("Recruitment (1000s)") + labs(shape="Estimate Type", fill = "Estimate Type", colour="Bottom Temperature\nAnomaly") + 
#   theme(axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.title = element_text(size = rel(2), hjust=0.5), 
#     legend.box.just = "center",legend.text = element_text(size = rel(1.5))) +
#   geom_ribbon(aes(Year, Recruitment, ymin=lo, ymax=hi), alpha=0.3, linetype = 0, inherit.aes = FALSE)
# plt

#######################

############################################################################################
#Fig. S3
#self-test results
#source(here::here("R","plot_functions.R"))

#Fig. S3
summarize_res_fn <- function(sims, rep_name = "SSB", index = 1, fit){
  out <- list()
  if(rep_name == "SSB" & is.null(index)){ #total SSB
    out$resid <- sapply(sims, \(x) {
      if(!is.null(x[["SSB"]])) return(apply(x[["SSB"]],1,sum)/apply(fit$rep[["SSB"]],1,sum) - 1)
      else return(rep(NA, NROW(fit$rep[["SSB"]])))
    })
  }
  if(rep_name == "SSB" & !is.null(index)){ #stock SSB
    out$resid <- sapply(sims, \(x) {
      if(!is.null(x[[rep_name]])) return(x[[rep_name]][,index]/fit$rep[[rep_name]][,index] - 1)
      else return(rep(NA, NROW(fit$rep[[rep_name]])))
    })
  }
  if(rep_name == "F"){ #total SSB
    out$resid <- sapply(sims, \(x) {
      if(!is.null(x[["F"]])) return(x[["F"]]/exp(fit$rep[["log_F_tot"]]) - 1)
      else return(rep(NA, length(fit$rep[["log_F_tot"]])))
    })
  }
  out$n <- apply(out$resid,1, \(x) sum(!is.na(x))) #by year 
  out$median_res <- apply(out$resid, 1, median, na.rm = TRUE)
  # print(out$median_res)
  out$bnds <- cbind(qbinom(0.025, out$n, 0.5)/out$n, qbinom(0.975, out$n, 0.5)/out$n)
  out$quantiles <- sapply(1:length(out$n), \(y) quantile(out$resid[y,], probs = c(out$bnds[y,1],0.5, out$bnds[y,2]), na.rm = TRUE))
  return(out)
}

make_df <- function(results_files, types, fit){
  df <- cbind.data.frame(year = numeric(), type = character(), stock = character(), median_y = numeric(), lo_y = numeric(),hi_y = numeric(), 
    median_median = numeric(), median_all = numeric())
  for(i in 1:length(results_files)){
    res <- readRDS(results_files[i])
    # res <- readRDS(here::here("results","self_test", "fit_1_fix_reml", "self_test_results.RDS"))
    print(i)
    summary_N <- summarize_res_fn(res[[1]], rep_name = "SSB", index = 1, fit)
    summary_S <- summarize_res_fn(res[[1]], rep_name = "SSB", index = 2, fit)
    summary_all <- summarize_res_fn(res[[1]], rep_name = "SSB", index = NULL, fit)
    summary_F <- summarize_res_fn(res[[1]], rep_name = "F", index = NULL, fit)
    df_i <- cbind.data.frame(year = fit$years, type = types[i], stock = "SSB North", median_y = summary_N$median_res, lo_y = summary_N$quantiles[1,],hi_y = summary_N$quantiles[3,], 
      median_median = median(summary_N$median_res, na.rm = TRUE), median_all = median(summary_N$resid, na.rm = TRUE))
    df_i <- rbind(df_i, cbind.data.frame(year = fit$years, type = types[i], stock = "SSB South", median_y = summary_S$median_res, lo_y = summary_S$quantiles[1,],hi_y = summary_S$quantiles[3,], 
      median_median = median(summary_S$median_res, na.rm = TRUE), median_all = median(summary_S$resid, na.rm = TRUE)))
    df_i <- rbind(df_i, cbind.data.frame(year = fit$years, type = types[i], stock = "SSB Total", median_y = summary_all$median_res, lo_y = summary_all$quantiles[1,],hi_y = summary_all$quantiles[3,], 
      median_median = median(summary_all$median_res, na.rm = TRUE), median_all = median(summary_all$resid, na.rm = TRUE)))
    df_i <- rbind(df_i, cbind.data.frame(year = fit$years, type = types[i], stock = "F Total", median_y = summary_F$median_res, lo_y = summary_F$quantiles[1,],hi_y = summary_F$quantiles[3,], 
      median_median = median(summary_F$median_res, na.rm = TRUE), median_all = median(summary_F$resid, na.rm = TRUE)))
    df <- rbind(df, df_i)
  }
  return(df)
}

fit_1 <- readRDS(here::here("results","fit_1.RDS"))
res_files <- here::here("results","self_test", c("fit_1", "fit_1_fix"), "self_test_results.RDS")
self_test_df <- make_df(res_files, types = c("Index CV Est", "Index CV Fix"), fit_1)
# self_test_df$type <- factor(self_test_df$type, levels = c("ML", "REML", "Low Error"))
library(ggplot2)
plt <- ggplot(self_test_df, aes(x = year, y = median_y)) + scale_colour_viridis_d() + 
    geom_point(size = 2) + geom_line() +
    facet_grid(stock  ~ type) + theme_bw() +
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_errorbar(aes(ymin = lo_y, ymax = hi_y), width = .01, linewidth = 1) +
    geom_line(aes(x = year, y = median_median), alpha = 0.5, linewidth = 2, linetype = "dashed") +
    geom_line(aes(x = year, y = median_all),alpha = 0.5, linewidth = 2, linetype = "dashed", color = "red") +
    ylab(bquote(Relative~Error~(SSB))) + xlab("Year")

cairo_pdf(file.path("paper", "self_test_results.pdf"), width = 12, height = 16)
plt
dev.off()

#For median CIs
  # https://www-users.york.ac.uk/~mb55/intro/cicent.htm
  # https://stats.stackexchange.com/questions/122001/confidence-intervals-for-median
  # http://www.jstor.com/stable/2957563

########################################
#projections


########################################

#Fig. S8
cairo_pdf(file.path("paper", "selectivity_south_plot.pdf"), width = 16, height = 12)
plot.selectivity(fits[[2]], blocks = c(3:4,7:8), block.names = c("South Commercial Fleet", "South Recreational Fleet", "South Recreational Catch/Effort Index", "South Spring Index"))
dev.off()

source(here::here("R","plot_functions.R"))
theme_set(theme_bw())
theme_update(strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = rel(2)), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)))

cairo_pdf(file.path("paper", "catch_plot.pdf"), width = 16, height = 16)
plot.catch(fits[[2]])
dev.off()

###################################################################
#Ecov, R, SSB, F same figure

summarize_res_fn <- function(fits, rep_name = "log_SSB", index = 1, age = 1){
  out <- list()
  if(rep_name == "log_NAA_rep"){
    out$Est <- sapply(fits, \(x) {
      return(TMB:::as.list.sdreport(x$sdrep, "Est", report = TRUE)[[rep_name]][index,index,,age])
    })
    out$SE <- sapply(fits, \(x) {
      return(TMB:::as.list.sdreport(x$sdrep, "Std", report = TRUE)[[rep_name]][index,index,,age])
    })
  } else {
    out$Est <- sapply(fits, \(x) {
      return(TMB:::as.list.sdreport(x$sdrep, "Est", report = TRUE)[[rep_name]][,index])
    })
    out$SE <- sapply(fits, \(x) {
      return(TMB:::as.list.sdreport(x$sdrep, "Std", report = TRUE)[[rep_name]][,index])
    })
  }
  out$lo <- out$Est + qnorm(0.025) * out$SE
  out$hi <- out$Est + qnorm(0.975) * out$SE
  return(out)
}

make_df <- function(fits, rep_name = "log_SSB", label = "SSB", age = 1){
  df <- cbind.data.frame(year = numeric(), model = character(), stock = character(), SSB = numeric(), lo = numeric(),hi = numeric())
  summary_N <- summarize_res_fn(fits, rep_name = rep_name, index = 1)
  print(NCOL(summary_N$Est))
  summary_S <- summarize_res_fn(fits, rep_name = rep_name, index = 2)
  mods <- paste0("$\\textit{M}_{",1:NCOL(summary_N$Est) - 1,"}$")
  print(mods)
  for(i in 1:NCOL(summary_N$Est)){
    df <- rbind(df, cbind.data.frame(year = fits[[1]]$years, model = mods[i], stock = "North", type = label, est = exp(summary_N$Est[,i]), lo = exp(summary_N$lo[,i]),hi = exp(summary_N$hi[,i])))
    df <- rbind(df, cbind.data.frame(year = fits[[1]]$years, model = mods[i], stock = "South", type = label, est = exp(summary_S$Est[,i]), lo = exp(summary_S$lo[,i]),hi = exp(summary_S$hi[,i])))
  }
  df$model <- factor(df$model, levels = mods)
  return(df)
}

df <- make_df(c(fits,fits_M_re), label = "SSB (mt)")
df <- rbind(df, make_df(c(fits,fits_M_re), rep_name = "log_Fbar", label = "Avg. F (Ages 6-7)"))
temp_df <- make_df(c(fits,fits_M_re), rep_name = "log_NAA_rep", label = "Recruitment (Millions)")
temp_df[,c("est","lo","hi")] <-temp_df[,c("est","lo","hi")]/1000
df <- rbind(df,temp_df)

ecov_res <- plot.ecov(fits[[2]])

df_use <- subset(df, model == levels(df$model)[2])

df_use <- rbind(df_use, cbind.data.frame(year = ecov_res$Year, model = levels(df$model)[2], stock = ecov_res$region, type = ecov_res$type, 
  est = ecov_res$est, lo = ecov_res$lo, hi = ecov_res$hi))

df_ecov_all <- subset(df_use, type %in% c("obs","pred"))
df_ecov_all$Type <- factor(c("Observed", "Posterior")[match(df_ecov_all$type, c("obs","pred"))])
plt_cols <- scales::viridis_pal(option = "turbo")(2)
names(plt_cols) <- levels(df_ecov_all$Type)
plt_cols <- scale_colour_manual(name = "Type",values = scales::viridis_pal(option = "turbo")(2))
df_ecov_all_pred <- subset(df_ecov_all, type %in% c("pred"))
df_ecov_all_obs <- subset(df_ecov_all, type %in% c("obs"))
df_ecov <- subset(df_ecov_all, stock == "North")
df_ecov_pred <- subset(df_ecov, type %in% c("pred"))
df_ecov_obs <- subset(df_ecov, type %in% c("obs"))
p_ecov_N <- ggplot(df_ecov, aes(x = year, y = est, colour = Type, fill = Type)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + 
    labs(color = "Estimate Type", fill = "Estimate Type") + xlab("Year") +
    facet_grid(~stock) + ylab("Bottom Temperature Anomaly") + geom_vline(xintercept = min(fits[[1]]$years), linetype = 2, alpha = 0.4, linewidth = 1.5) +
    geom_point(data = df_ecov_obs, aes(x = year, y = est, colour = Type), size = 1) + geom_pointrange(data = df_ecov_obs, aes(ymin = lo, ymax = hi, colour = Type)) + 
    geom_line(data = df_ecov_pred, aes(x = year, y = est, colour = Type)) + geom_ribbon(data = df_ecov_pred, aes(ymin=lo, ymax=hi, fill = Type), alpha=0.4, linetype = 0) 
p_ecov_N


df_ecov <- subset(df_ecov_all, stock == "South")
df_ecov_pred <- subset(df_ecov, type %in% c("pred"))
df_ecov_obs <- subset(df_ecov, type %in% c("obs"))
p_ecov_S <- ggplot(df_ecov, aes(x = year, y = est, colour = Type, fill = Type)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + 
    labs(color = "Estimate Type", fill = "Estimate Type") + xlab("Year") +
    facet_grid(~stock) + ylab("Bottom Temperature Anomaly") + geom_vline(xintercept = min(fits[[1]]$years), linetype = 2, alpha = 0.4, linewidth = 1.5) +
    geom_point(data = df_ecov_obs, aes(x = year, y = est, colour = Type), size = 1) + geom_pointrange(data = df_ecov_obs, aes(ymin = lo, ymax = hi, colour = Type)) + 
    geom_line(data = df_ecov_pred, aes(x = year, y = est, colour = Type)) + geom_ribbon(data = df_ecov_pred, aes(ymin=lo, ymax=hi, fill = Type), alpha=0.4, linetype = 0) 
p_ecov_S

p_ecov <- ggplot(df_ecov_all, aes(x = year, y = est, colour = Type, fill = Type)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") + 
    labs(color = "Estimate Type", fill = "Estimate Type") + xlab("Year") +
    facet_grid(~stock) + ylab("Bottom Temperature Anomaly") + geom_vline(xintercept = min(fits[[1]]$years), linetype = 2, alpha = 0.4, linewidth = 1.5) +
    geom_point(data = df_ecov_all_obs, aes(x = year, y = est, colour = Type), size = 1) + geom_pointrange(data = df_ecov_all_obs, aes(ymin = lo, ymax = hi, colour = Type)) + 
    geom_line(data = df_ecov_all_pred, aes(x = year, y = est, colour = Type)) + geom_ribbon(data = df_ecov_all_pred, aes(ymin=lo, ymax=hi, fill = Type), alpha=0.4, linetype = 0) 
p_ecov

df_brps <- readRDS(here::here("results", "df_brps_rev.RDS"))

R_df <- subset(df_brps, Proj_type == "Last" & Year < 2022 & XSPR_R_type == "Random Effect" & name %in% c("Random Effect", "Expected"))


Ecov_N <- fits[[2]]$rep$Ecov_out_R[1,,1]
R_df[,"BT"] <- c(fits[[2]]$rep$Ecov_out_R[1,,1],fits[[2]]$rep$Ecov_out_R[2,,2])
R_df_N <- subset(R_df, Region =="North")
ord <- order(Ecov_N)
years <- 1:length(ord)

pal <- viridis::viridis_pal(option = "H")(length(ord))

R_df_N$name <- factor(R_df_N$name)
# plt <- ggplot(R_df_N, aes(Year, val)) + 
#   geom_line(data=subset(R_df_N, name == "Expected"), colour = gray(0.7), linewidth = 1.5) +  
#   scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") +
#   scale_colour_viridis_d(begin = 0.2, end = 0.8) + 
#   geom_ribbon(data=subset(R_df_N, name == "Expected"), aes(ymin=lo, ymax=hi, fill = name), fill= gray(0.7), alpha=0.3, linetype = 0) +
#   geom_point(aes(fill = BT), shape = 21, size = 4) +  
#   ylab("Recruitment (1000s)") + labs(shape="Estimate Type", colour="Bottom Temperature\nAnomaly", fill = "Estimate Type") #+
# plt

(p_ecov_N + theme(legend.position="none") + xlab("")) + (p_ecov_S  + ylab("") + xlab(""))


#Fig. S1
theme_set(theme_bw())
theme_update(strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = rel(2)), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)))

res_dir <- here::here("results","jitter", "fit_0_rev")
fit_0 <- readRDS(here::here("results", "fit_0_rev.RDS"))
fit_0_jit_res <- readRDS(here::here(res_dir, "fit_0_rev_jitter_results.RDS"))
jitters <- fit_0_jit_res[[1]]


y <- cbind(nll = sapply(jitters, \(x) x$obj), mgrad = sapply(jitters, \(x) max(abs(x$grad))))
ordy <- order(y[,1])
y[ordy,]
cbind(round(fit_0$opt$par,2), round(fit_0$opt$par - jitters[[ordy[2]]]$par,2),round(fit_0$opt$par - jitters[[ordy[3]]]$par,2),round(fit_0$opt$par - jitters[[ordy[4]]]$par,2),round(fit_0$opt$par - jitters[[ordy[5]]]$par,2))
# ordered fits 3 and 4 look identical

# pcol <- c("black","gray","red")[sapply(jitters, \(x) findInterval(max(abs(x$grad)), c(1e-20, 1e-10, 1)))]
pcol <- c("black","gray","red")[sapply(jitters, \(x) findInterval(y[,2]), c(1e-2, 1))+1)]
badgrads <- cbind(iter = which(y[,2]>1), y[y[,2]>1,])
char_badgrad <- format(badgrads[,3])
p1 <- as.numeric(substr(char_badgrad,1,3))
p2 <- as.numeric(substr(char_badgrad,nchar(char_badgrad)-1,nchar(char_badgrad)))
labs <- lapply(1:NROW(badgrads), \(i) bquote(.(p1[i])%*%10^.(p2[i])))

cairo_pdf(file.path("paper", "fit_0_jitter_plt.pdf"), width = 7, height = 7)
plot(sapply(jitters, \(x) x$obj), ylab = "Negative Log-Likelihood", xlab = "Jitter #", pch = 19, cex = 1.2, col = pcol, ylim =c(-980,-940))
#text(x = badgrads[,1], y = badgrads[,2], labels = bquote(.(p1)%*%10^.(p2)), adj = c(1,0))
lapply(1:NROW(badgrads), \(i) text(x = badgrads[i,1], y = badgrads[i,2]+0.5, labels = labs[[i]], adj = c(1,0)))
abline(h = fit_0$opt$obj, cex = 1)
dev.off()

##############################################################################################

#Fig. S2
res_dir <- here::here("results","jitter", "fit_1")
fit_1_jit_res <- readRDS(here::here(res_dir, "fit_1_jitter_results.RDS"))
fit_1 <- readRDS(here::here("results", "fit_1.RDS")))
jitters <- fit_1_jit_res[[1]]
sapply(jitters, \(x) max(abs(x$grad)))[1:10]
sort(sapply(jitters, \(x) x$obj))
y <- order(sapply(jitters, \(x) x$obj))
sapply(jitters, \(x) max(abs(x$grad)))[y]
badpar <- sapply(y[1:3], \(x) jitters[[x]]$par) - fit_1$opt$par
round(badpar,2)
pcol <- c("black","gray","red")[sapply(jitters, \(x) findInterval(max(abs(x$grad)), c(1e-6,1e-2, 10)))]
cairo_pdf(file.path("paper", "fit_1_jitter_plt.pdf"), width = 7, height = 7)
plot(sapply(jitters, \(x) x$obj), ylab = "Negative Log-Likelihood", xlab = "Jitter #", pch = 19, cex = 1.2, col = pcol)
abline(h = fit_1$opt$obj, cex = 1)
dev.off()

###############################################
