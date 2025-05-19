library("wham", lib.loc = "c:/work/wham/old_packages/53e236b")
# library("wham", lib.loc = "c:/work/wham/old_packages/24c8156")
library(here)
source(here::here("R","make_wham_inputs.R"))

#don't use prior
move$use_prior <- array(0, dim = c(2,length(seasons),2,1))

fit1 <- readRDS(here::here("results", paste0("fits_no_M_re_rev_1.RDS")))[[2]]


#don't use prior
input_1 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_1, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp)
input_1$fleet_names <- paste0(rep(c("North_", "South_"),each = 2), temp$fleet_names)
input_1$index_names <- paste0(rep(c("North_", "South_"),c(2,2)), temp$index_names)
input_1$par <- fit1$parList
input_1$par$mu_prior_re[] <- 0
input_1$map$trans_mu <- factor(rep(NA, length(input_1$par$trans_mu)))
# input_1$par <- fit_0$parList
# fit_1_1 <- fit_wham(input_1, do.fit = FALSE, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)
fit_1_1 <- fit_wham(input_1, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)
saveRDS(fit_1_1, here::here("results", "fit_1_fix_mu.RDS"))


#estimate median M (constant across stocks, regions)
source(here::here("R","make_wham_inputs.R"))
M <- list(mean_model = "estimate-M")
input_2 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_1, catch_info = catch_info, M = M,
	index_info = index_info, age_comp = age_comp)
input_2$fleet_names <- paste0(rep(c("North_", "South_"),each = 2), temp$fleet_names)
input_2$index_names <- paste0(rep(c("North_", "South_"),c(2,2)), temp$index_names)
input_2$par <- fits[[2]]$parList

#input_2$par <- fit_1_1$parList

# fit_1_2 <- fit_wham(input_2, do.fit = FALSE, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)
fit_1_2 <- fit_wham(input_2, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)
#estimated M is much larger: 0.72
exp(fit_1_2$opt$par["Mpars"])
fit_1_2$rep$mu[1,8,1,1,,] #not much different when estimating M
fits[[2]]$rep$mu[1,8,1,1,,]
saveRDS(fit_1_2, here::here("results", "fit_1_estimate_M.RDS"))

#no movement
source(here::here("R","make_wham_inputs.R"))
#move$use_prior <- array(0, dim = c(2,length(seasons),2,1))
basic_info$NAA_where <- NULL
input_3 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, ecov = ecov_1, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp)
input_3$fleet_names <- paste0(rep(c("North_", "South_"),each = 2), temp$fleet_names)
input_3$index_names <- paste0(rep(c("North_", "South_"),c(2,2)), temp$index_names)
input_3$par <- fit1$parList
input_3$par$mu_prior_re[] <- 0
input_3$par$trans_mu[] <- 0
# input_1$par <- fit_0$parList
# fit_1_3 <- fit_wham(input_3, do.fit = FALSE, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)
fit_1_3 <- fit_wham(input_3, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)
saveRDS(fit_1_3, here::here("results", "fit_1_no_move.RDS"))

ny <- length(fit1$years)
df <- cbind.data.frame(m = "italic(M)[1]", year = rep(fit1$years,4), type = rep(c("SSB", "Fbar"), each = 2*ny), val =c(fit1$rep$SSB/1000,fit1$rep$Fbar[,5:6]), region = rep(rep(c("North","South"), each = ny),2))
mods <- c("fit_1_1", "fit_1_2", "fit_1_3")
modnames <- c("Fixed movement", "Estimate M", "No movement")
for( i in 1:3) df <- rbind(df, 
	cbind.data.frame(m = modnames[i], year = rep(get(mods[i])$years,4), type = rep(c("SSB", "Fbar"), each = 2*ny), val = c(get(mods[i])$rep$SSB/1000, get(mods[i])$rep$Fbar[,5:6]), region = rep(rep(c("North","South"), each = ny),2)))
df$m <- factor(df$m, levels = c( "italic(M)[1]", modnames))

library(ggplot2)
library(patchwork)

theme_set(theme_bw())
theme_update(strip.text = element_text(size = rel(2)), strip.placement = "outside", strip.background = element_rect(), #fill = "transparent"), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)))
line.labs <- list(bquote(italic(M)[1]), "Fixed movement", "Estimate M", "No movement")

this_df <- subset(df, type == "SSB")
SSB_plt <- ggplot(this_df, aes(x = year, y = val, colour = m)) + scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", labels = line.labs) +
    geom_line(linewidth = 1.5) +
    facet_grid(~ region, switch = "y") +
    xlab("Year") + ylab("SSB (kmt)") + labs(colour = "Model") 
SSB_plt

cnms <- c(North = "", South = "")
this_df <- subset(df, type == "Fbar")
Fbar_plt <- ggplot(this_df, aes(x = year, y = val, colour = m)) + scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", labels = line.labs) +
    geom_line(linewidth = 1.5) +
    facet_grid(~ region, switch = "y", labeller = labeller(region = cnms)) +
    xlab("Year") + ylab(expression(bar(italic(F))~~"(Ages 6-7)")) + labs(colour = "Model") 
Fbar_plt

design <- c(area(1,1,1,1), area(2,1,2,1))
(SSB_plt +  xlab("")  + theme(axis.title.x = element_blank(), axis.text.x=element_blank())) + 
  (Fbar_plt  + theme(strip.text.x=element_blank(), legend.position = "none")) + 
  plot_layout(design = design)

cairo_pdf(here::here("paper", "SSB_F_sensitivity_plots.pdf"), width = 30*2/3, height = 20*2/3)
design <- c(area(1,1,1,1), area(2,1,2,1))
(SSB_plt +  xlab("")  + theme(axis.title.x = element_blank(), axis.text.x=element_blank())) + 
  (Fbar_plt  + theme(strip.text.x=element_blank(), legend.position = "none")) + 
  plot_layout(design = design)
dev.off()