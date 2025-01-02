library(wham)
library(ggplot2)
library(dplyr)
library(patchwork)
source(here::here("R","plot_functions.R"))

theme_set(theme_bw())
theme_update(strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = rel(2)), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)))

fits <- readRDS(file.path("results","fits_no_M_re_rev.RDS"))
fits_M_re <- readRDS(file.path("results","fits_M_re_better.RDS"))



res_dir <- here::here("results","jitter", "fit_0")
fit_0 <- readRDS(here::here("results", "fit_0_rev.RDS"))
fit_0_jit_res <- readRDS(here::here(res_dir, "fit_0_jitter_results.RDS"))
# jitter_files <- dir(here::here("results","jitter","fit_0"), full.names = TRUE, pattern = glob2rx("jitter_sim*.RDS"))
# jitters <- lapply(fit_0_jit_res, \(x) readRDS(x))
jitters <- fit_0_jit_res[[1]]
sort(sapply(jitters, \(x) x$obj))
y <- order(sapply(jitters, \(x) x$obj))
# sapply(jitters, \(x) max(abs(x$grad)))[y]
cbind(round(jitters[[y[3]]]$par,2), round(jitters[[y[1]]]$par - jitters[[y[3]]]$par,2))

temp <- jitter_wham(fit_RDS = here::here("results","fit_0_rev.RDS"), n_jitter = 1, initial_vals = rbind(fit_0_jit_res[[2]][y[1],]))
jit_mod <- readRDS(here::here("results", "fit_0_rev.RDS"))
jit_mod$env$last.par.best <- temp[[1]][[1]]$last.par.best
jit_mod$fn(temp[[1]][[1]]$par) #works!

jit_mod <- readRDS(here::here("results", "fit_0_rev.RDS"))
jit_mod$par <- fit_0_jit_res[[2]][y[1],]
temp1 <- fit_tmb(jit_mod, n.newton = 3, do.sdrep = FALSE)

jit_mod <- readRDS(here::here("results", "fit_0_rev.RDS"))
jit_mod$par <- fit_0_jit_res[[2]][y[1],]
temp11 <- fit_tmb(jit_mod, n.newton = 0, do.sdrep = FALSE)
jit_mod <- readRDS(here::here("results", "fit_0_rev.RDS"))
jit_mod$fn(temp11$opt$par)

temp2 <- nlminb(fit_0_jit_res[[2]][y[1],], fit_0$fn, fit_0$gr)
fit_0$par <- fit_0_jit_res[[2]][y[1],]

pcol <- c("black","gray","red")[sapply(jitters, \(x) findInterval(max(abs(x$grad)), c(1e-20, 1e-10, 1)))]
cairo_pdf(file.path("paper", "fit_0_jitter_plt.pdf"), width = 7, height = 7)
plot(sapply(jitters, \(x) x$obj), ylab = "Negative Log-Likelihood", xlab = "Jitter #", pch = 19, cex = 1.2, col = pcol)
abline(h = fit_0$opt$obj, cex = 1)
dev.off()


res_dir <- here::here("results","jitter", "fit_1")
fit_1_jit_res <- readRDS(here::here(res_dir, "fit_1_jitter_results.RDS"))

jitters <- fit_1_jit_res[[1]]
sapply(jitters, \(x) max(abs(x$grad)))[1:10]
sort(sapply(jitters, \(x) x$obj))
y <- order(sapply(jitters, \(x) x$obj))
sapply(jitters, \(x) max(abs(x$grad)))[y]
badpar <- sapply(y[1:3], \(x) jitters[[x]]$par) - fits[[2]]$opt$par
round(badpar,2)
pcol <- c("black","gray","red")[sapply(jitters, \(x) findInterval(max(abs(x$grad)), c(1e-6,1e-2, 1e4)))]
cairo_pdf(file.path("paper", "fit_1_jitter_plt.pdf"), width = 7, height = 7)
plot(sapply(jitters, \(x) x$obj), ylab = "Negative Log-Likelihood", xlab = "Jitter #", pch = 19, cex = 1.2, col = pcol)
abline(h = fits[[2]]$opt$obj, cex = 1)
dev.off()
###############################################
#BT for north and south on same plot


cairo_pdf(file.path("paper", paste0("Ecov.pdf")), height = 10, width = 10)
plot.ecov(fits[[2]])
legend("topleft", legend = c("North", "South"), col = mypalette(2), lty =  1, pch = 19)
dev.off()

library(latex2exp) #https://github.com/stefano-meschiari/latex2exp

cairo_pdf(file.path("paper", paste0("Ecov_M1_rel_M0.pdf")), height = 10, width = 10)
par(mar = c(4,7,1,1), oma = c(0,0,0,0))
plot(fits[[1]]$input$years_Ecov, (fits[[2]]$rep$Ecov_x[,1] - fits[[1]]$rep$Ecov_x[,1])/abs(fits[[1]]$rep$Ecov_x[,1]), xlab = "Year", ylab = TeX("$\\frac{\\widehat{X}\\left(M_1\\right) - \\widehat{X}\\left(M_0\\right)}{|\\widehat{X}\\left(M_0\\right)|}$" ))
grid(col = gray(0.7), lwd = 2, lty = 2)
dev.off()

cbind(fits[[1]]$input$years_Ecov, fits[[1]]$rep$Ecov_x)
plot(fits[[2]]$rep$Ecov_out_R[1,,1], log(fits[[2]]$rep$NAA[1,1,,1]))
abline(a = fits[[2]]$parList$mean_rec_pars[1,1], b = fits[[2]]$parList$Ecov_beta_R[1,1,1])
# plot(fits[[2]]$rep$Ecov_out_R[2,,2], log(fits[[2]]$rep$NAA[2,2,,1]))
# abline(a = fits[[2]]$parList$mean_rec_pars[2,1], b = fits[[2]]$parList$Ecov_beta_R[1,1,1])

# lm(log(fits[[2]]$rep$NAA[2,2,,1]) ~ poly(fits[[2]]$rep$Ecov_out_R[2,,2], order = 2))
# lm(log(fits[[2]]$rep$NAA[2,2,,1]) ~ fits[[2]]$rep$Ecov_out_R[2,,2])
# lm(log(fits[[2]]$rep$NAA[1,1,,1]) ~ fits[[2]]$rep$Ecov_out_R[1,,1])
c(fits[[2]]$parList$Ecov_beta_R[1,1,1], sapply(1:7, \(x) fits[[2]]$peels[[x]]$parList$Ecov_beta_R[1,1,1]))


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
# strip_labels <- list("$\\frac{\\overline{\\textit{F}}$}{\\overline{\\textit{F}}(M_1)}$", "$\\frac{Recruitment}{Recruitment(M_1)}$", "$\\frac{SSB}{SSB(M_1)}$")
# names(strip_labels) <- levels(df$typ)
# latex2exp::TeX(strip_labels[[1]])
# strip_labeller <- function(variable,value){
#   return(latex2exp::TeX(strip_labels[value]))
# }

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
cairo_pdf(file.path("paper", "SSB_F_R_rel_M1.pdf"), width = 12, height = 16)
plt
dev.off()


###############################################



###############################################
# Ecov effects on M and R
fit_3 <- readRDS(file.path("results","fits_M_re_better.RDS"))[[5]]
Ecov <- fit_3$rep$Ecov_out_M[1,1,1,,][,1]
M_out <- fit_3$rep$MAA[1,1,,1]
ord <- order(Ecov)

Ecov_plt <- seq(min(Ecov),max(Ecov),0.01)
M <- exp(fit_3$parList$Mpars[1,1,1] + fit_3$parList$Ecov_beta_M[1,1,1,1,1] * Ecov_plt)
# pal <- mypalette(length(ord))
pal <- viridis::viridis_pal(option = "H")(length(ord))
# plot(Ecov_plt, M, type = 'l', lwd = 2, ylim = c(0, max(M_out)))
# points(Ecov[ord], M_out[ord], col = pal, pch = 19)

df <- cbind.data.frame(M=M, Temperature  = Ecov_plt, group = 1)
df <- rbind(df, cbind.data.frame(M = M_out, Temperature = Ecov, group = 2))
plt <- ggplot(df, aes(Temperature, M, colour = Temperature)) +
  geom_line(data = subset(df, group == 1), size = 1.2) +
  geom_point(data = subset(df, group == 2)) + 
  scale_colour_gradientn(colours = pal)

#SSB <- seq(1,max(fits[[2]]$rep$SSB[,1]),100)
#R <- fits[[2]]$rep$pred_NAA[1,1,,1]
fit0 <- readRDS(file.path("results","fits_no_M_re_rev.RDS"))[[1]]
fit1 <- readRDS(file.path("results","fits_no_M_re_rev.RDS"))[[2]]
R_out <- fit1$rep$NAA[1,1,,1]
Ecov <- fit1$rep$Ecov_out_R[1,,1]
ord <- order(Ecov)
Ecov_plt <- seq(min(Ecov),max(Ecov),0.01)
Expected_R <- exp(fit1$parList$mean_rec_pars[1,1] + fit1$parList$Ecov_beta_R[1,1,1] * Ecov_plt)
pal <- viridis::viridis_pal(begin = 0.2, end = 0.8, option = "H")(length(ord))
# pal <- mypalette(length(ord))
years <- 1:length(ord)
#R <- sapply(1:NROW(ab), function(x) ab[x,1]*SSB/(1 + ab[x,2]*SSB))
# Rmax <- ab[,1]*max(SSB)/(1 + ab[,2]*max(SSB))
# matplot(SSB, R, type = 'l', col = pal, lty = 1, lwd = 2)
cbind(exp(fit1$parList$mean_rec_pars[1,1] + fit1$parList$Ecov_beta_R[1,1,1] * Ecov), R_out)

plot(log(fit0$rep$NAA[1,1,-1,1]), log(fit1$rep$NAA[1,1,-1,1]), ylab = "log(Recruitment) with BT Effect", xlab = "log(Recruitment) without BT Effect")
abline(0,1)
plot(fit1$years[-1], fit1$rep$pred_NAA[1,1,-1,1], ylab = "Expected Recruitment (1000s)", xlab = "Year", type = 'l')
lines(fit1$years[-1], fit0$rep$pred_NAA[1,1,-1,1], col = "gray")

df_brps <- readRDS(here::here("results", "df_brps.RDS"))
R_df <- subset(df_brps, Proj_type == "Last" & Year < 2022 & XSPR_R_type == "Random Effect" & name %in% c("Random Effect", "Expected"))
R_df[,"BT"] <- c(fits[[2]]$rep$Ecov_out_R[1,,1],fits[[2]]$rep$Ecov_out_R[2,,2])
R_df_N <- subset(R_df, Region =="North")
# Ecov_N <- fits[[2]]$rep$Ecov_out_R[1,,1]
# ord <- order(Ecov_N)
# years <- 1:length(ord)
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
cairo_pdf(file.path("paper", "best_R_Ecov.pdf"), width = 14, height = 12)
plt
dev.off()

# df <- cbind.data.frame(R=Expected_R, BT = Ecov_plt, group = 1)
# df <- rbind(df, cbind.data.frame(R = R_out, BT = Ecov, group = 2))
# plt <- ggplot(df, aes(BT, R, colour = BT)) +
#   geom_line(data = subset(df, group == 1), size = 1.2) +
#   geom_point(data = subset(df, group == 2)) + 
#   scale_colour_gradientn(colours = pal)

# plt <- ggplot(df, aes(Year, Recruitment, colour = BT, shape = Type, ymin=lo, ymax=hi, fill = Type)) + theme_bw() + 
plt <- ggplot(df, aes(Year, Recruitment, colour = order_Ecov, shape = Type, ymin=lo, ymax=hi)) + theme_bw() + 
  # geom_line(data = subset(df,  Type == "Expected"), colour = gray(0.7, alpha = 0.4), size = 2) +
  geom_line(aes(Year, Recruitment, colour = Type), size = 1.5, inherit.aes = FALSE) +
  geom_point(aes(fill=order_Ecov), size=4, shape=21, stroke=0) +
  scale_fill_gradientn(colours = pal) +
  # geom_point(size = 2) + scale_colour_gradientn(colours = pal) +
  # geom_point(data = subset(df, Type == "Expected"), size = 2) +
  # geom_point(data = subset(df, Type == "Random Effect"), size = 2) + 
  ylab("Recruitment (1000s)") + labs(shape="Estimate Type", fill = "Estimate Type", colour="Bottom Temperature\nAnomaly") + 
  theme(axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.title = element_text(size = rel(2), hjust=0.5), 
    legend.box.just = "center",legend.text = element_text(size = rel(1.5))) +
  geom_ribbon(aes(Year, Recruitment, ymin=lo, ymax=hi), alpha=0.3, linetype = 0, inherit.aes = FALSE)
  # geom_ribbon(data = subset(df, Type == "Expected"), aes(Year, Recruitment, ymin=lo, ymax=hi), fill = gray(0.7), alpha=0.3, linetype = 0, inherit.aes = FALSE)
  # geom_ribbon(fill = gray(0.7), alpha=0.3, linetype = 0) #+
  # geom_ribbon(alpha=0.3, linetype = 0) #+
  # scale_fill_manual(values = c("Expected"="gray", "Random Effect"="Transparent")) +
  # scale_colour_gradientn(colours = pal)# + scale_fill_gradientn(colours = pal)
plt
# df <- cbind.data.frame(year = rep(years,each = length(SSB)), SSB = SSB, R = c(R), GSI = rep(Ecov, each = length(SSB)), group = 1)
# df <- rbind.data.frame(df, cbind.data.frame(year = years, SSB = fit1$rep$SSB, R = fit1$rep$NAA[1,1,,1], GSI = Ecov, group = 2))
# plt <- ggplot(df, aes(SSB, R, colour = GSI, group = factor(year))) +
#   geom_line(data = subset(df, group == 1), size = 1.2) +
#   geom_point(data = subset(df, group == 2)) + 
#   scale_colour_gradientn(colours = pal)
# ggsave(file.path("paper", "best_R_ecov.png"), plt)

#######################
#prior/posterior for movement parameter
dlogitnorm <- function(x, p, sigma){

  y <- log(x/(1-x))
  mu <- log(p/(1-p))
  fx <- dnorm(y, mu, sigma, log = F) /(x*(1-x)) 
  return(fx)
}

#png(file.path("paper", "move_prior_post.png"), width = 7*144, height = 7*144, res = 144, pointsize = 12, type = "cairo-png")
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

#######################

#######################
#self-test results
make_df <- function(results_files, types, fit){
  df <- cbind.data.frame(year = numeric(), type = character(), stock = character(), median_y = numeric(), lo_y = numeric(),hi_y = numeric(), 
    median_median = numeric(), median_all = numeric())
  for(i in 1:length(results_files)){
    res <- readRDS(results_files[i])
    # res <- readRDS(here::here("results","self_test", "fit_1_fix_reml", "self_test_results.RDS"))
    summary_N <- summarize_res_fn(res[[1]], rep_name = "SSB", index = 1, fit)
    summary_S <- summarize_res_fn(res[[1]], rep_name = "SSB", index = 2, fit)
    df_i <- cbind.data.frame(year = fit$years, type = types[i], stock = "North", median_y = summary_N$median_res, lo_y = summary_N$quantiles[1,],hi_y = summary_N$quantiles[3,], 
      median_median = median(summary_N$median_res, na.rm = TRUE), median_all = median(summary_N$resid, na.rm = TRUE))
    df_i <- rbind(df_i, cbind.data.frame(year = fit$years, type = types[i], stock = "South", median_y = summary_S$median_res, lo_y = summary_S$quantiles[1,],hi_y = summary_S$quantiles[3,], 
      median_median = median(summary_S$median_res, na.rm = TRUE), median_all = median(summary_S$resid, na.rm = TRUE)))
    df <- rbind(df, df_i)
  }
  return(df)
}
fit_1 <- readRDS(file.path("results","fits_no_M_re_rev.RDS"))[[2]]
res_files <- here::here("results","self_test", c("fit_1_fix", "fit_1_fix_reml", "fit_1_low_error"), "self_test_results.RDS")
self_test_df <- make_df(res_files, types = c("ML", "REML", "Low Error"), fit = fit_1)
self_test_df$type <- factor(self_test_df$type, levels = c("ML", "REML", "Low Error"))
library(ggplot2)
plt <- ggplot(self_test_df, aes(x = year, y = median_y)) + scale_colour_viridis_d() + 
    geom_point(size = 2) + geom_line() +
    facet_grid(stock  ~ type) + theme_bw() +
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_errorbar(aes(ymin = lo_y, ymax = hi_y), width = .01, linewidth = 1) +
    geom_line(aes(x = year, y = median_median),linewidth = 2, linetype = "dashed") +
    geom_line(aes(x = year, y = median_all),linewidth = 2, linetype = "dashed", color = "red") +
    ylab(bquote(Relative~Error~(SSB))) + xlab("Year")

cairo_pdf(file.path("paper", "self_test_results.pdf"), width = 20, height = 12)
plt
dev.off()

ggsave(here::here("paper", "self_test_results.png"), plt, width = 20, height = 12, units = "in")

#For median CIs
  # https://www-users.york.ac.uk/~mb55/intro/cicent.htm
  # https://stats.stackexchange.com/questions/122001/confidence-intervals-for-median
  # http://www.jstor.com/stable/2957563

########################################
#projections


########################################

#plot selectivity at time and age for blocks with RE.
fit <- readRDS(file.path("results","fits_no_M_re_rev.RDS"))[[2]]
source(here::here("R","plot_functions.R"))
cairo_pdf(file.path("paper", "selectivity_re_plot.pdf"), width = 16, height = 12)
plot.selectivity(fits[[2]], blocks = c(1:2,5:6), block.names = c("North Commercial Fleet", "North Recreational Fleet", "North Recreational Catch/Effort Index", "North Spring Index"))
dev.off()


source(here::here("R","plot_functions.R"))
theme_set(theme_bw())
theme_update(strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = rel(2)), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)))

cairo_pdf(file.path("paper", "catch_plot.pdf"), width = 16, height = 16)
plot.catch(fits[[2]])
dev.off()


#reference points
#R1 = annual NAA, R3 = annual pred_NAA

proj_R1 <- lapply(1:3, \(x) readRDS(here::here("results",paste0("m1_proj_",x,"_R_opt1.RDS"))))
proj_R3 <- lapply(1:3, \(x) readRDS(here::here("results",paste0("m1_proj_",x,"_R_opt3.RDS"))))
# proj_1_R1 <- readRDS(here::here("results","m1_proj_1_R_opt1.RDS"))
# proj_2_R1 <- readRDS(here::here("results","m1_proj_2_R_opt1.RDS"))
# proj_3_R1 <- readRDS(here::here("results","m1_proj_3_R_opt1.RDS"))
# proj_1_R3 <- readRDS(here::here("results","m1_proj_1_R_opt3.RDS"))
# proj_2_R3 <- readRDS(here::here("results","m1_proj_2_R_opt3.RDS"))
# proj_3_R3 <- readRDS(here::here("results","m1_proj_3_R_opt3.RDS"))

se_R1 <- lapply(proj_R1, \(x) TMB:::as.list.sdreport(x$sdrep, what = "Std", report = TRUE))
se_R3 <- lapply(proj_R3, \(x) TMB:::as.list.sdreport(x$sdrep, what = "Std", report = TRUE))
# se_1_R1 <- TMB:::as.list.sdreport(proj_1_R1$sdrep, what = "Std", report = TRUE)
# se_2_R1 <- TMB:::as.list.sdreport(proj_2_R1$sdrep, what = "Std", report = TRUE)
# se_3_R1 <- TMB:::as.list.sdreport(proj_3_R1$sdrep, what = "Std", report = TRUE)
# se_1_R3 <- TMB:::as.list.sdreport(readRDS(here::here("results","m1_proj_1_sdrep_fjp.RDS")), what = "Std", report = TRUE)
# se_2_R3 <- TMB:::as.list.sdreport(readRDS(here::here("results","m1_proj_2_sdrep_fjp.RDS")), what = "Std", report = TRUE)
# se_3_R3 <- TMB:::as.list.sdreport(readRDS(here::here("results","m1_proj_3_sdrep_fjp.RDS")), what = "Std", report = TRUE)
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

fjp <- list() 
for(i in 1:3) fjp[[i]] <- readRDS(here::here("results",paste0("m1_proj_",i,"_sdrep_fjp.RDS")))


#have to do se(predR) by hand
source(here::here("R","plot_functions.R"))
temp <- make_pred_R_cv(proj_R3[[1]], fjp[[1]], type = 1)
temp <- rbind(temp, make_pred_R_cv(proj_R3[[1]], fjp[[1]], type = 1, region = 2))
temp <- rbind(temp,make_pred_R_cv(proj_R3[[2]], fjp[[2]], type = 2))
temp <- rbind(temp,make_pred_R_cv(proj_R3[[2]], fjp[[2]], type = 2, region = 2))
temp <- rbind(temp,make_pred_R_cv(proj_R3[[3]], fjp[[3]], type = 3))
temp <- rbind(temp,make_pred_R_cv(proj_R3[[3]], fjp[[3]], type = 3, region = 2))

df_brps <- rbind(df_brps, cbind.data.frame(Year = proj_R3[[1]]$years_full, 
  val = temp[,1],
  cv = temp[,2],
  name = "Expected",
  Region = rep(c("North","South"), each = ny_full),
  Proj_type = rep(c("AR1", "Last", "Linear"), each = 2*ny_full), 
  XSPR_R_type = rep(c("Random Effect", "Expected"), each = 2*3*ny_full)))

df_brps$lo <- df_brps$val*exp(-qnorm(0.975)*df_brps$cv)
df_brps$hi <- df_brps$val*exp(qnorm(0.975)*df_brps$cv)

saveRDS(df_brps, here::here("results", "df_brps.RDS"))

# df <- cbind.data.frame(Year = proj_1_R3$years_full, SSB = exp(c(sapply(paste0("proj_",1:3,"_R", rep(c(1,3),each = 3)), \(x) get(x)$rep$log_SSB_FXSPR[,1]))),
#   Proj_type = rep(c("AR1", "Last", "Linear"), each = length(proj_1_R3$years_full)), R_type = rep(c("Estimated", "Expected"), each = 3*length(proj_1_R3$years_full)))

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

source(here::here("R","plot_functions.R"))
x <- kobe.plot(proj, status.years = length(proj$years), static = FALSE)
df <- x$df.ellipse

source(here::here("R","plot_functions.R"))
x <- get.status.results(proj)
y <- get.status.results(proj, static = TRUE)

ggplot(x$borders, aes(x = x, y = y, fill = region)) + 
    geom_polygon() + coord_cartesian(ylim = c(0, max(x$borders$y)), xlim = c(0,max(x$borders$x))) + 
    scale_fill_manual(values = x$poly.colors) + theme(legend.position="none")
df <- x$df.ellipse
df$ptype <- factor(df$ptype)
df$region <- factor(df$region)
  gg.kobe <- ggplot(df, aes(x = x, y = y)) + facet_grid(~region)+ 
    geom_polygon(subset(x$df.ellipse, ptype == "ellipse"), aes(x = x, y = y)) + 

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
  # print(out$median_res)
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

df_brps <- readRDS(here::here("results", "df_brps.RDS"))

R_df <- subset(df_brps, Proj_type == "Last" & Year < 2022 & XSPR_R_type == "Random Effect" & name %in% c("Random Effect", "Expected"))

# R_out <- fits[[2]]$rep$NAA[1,1,,1]

Ecov_N <- fits[[2]]$rep$Ecov_out_R[1,,1]
R_df[,"BT"] <- c(fits[[2]]$rep$Ecov_out_R[1,,1],fits[[2]]$rep$Ecov_out_R[2,,2])
R_df_N <- subset(R_df, Region =="North")
ord <- order(Ecov_N)
years <- 1:length(ord)

# Ecov_plt <- seq(min(Ecov_N),max(Ecov_N),0.01)
# Expected_R <- exp(fits[[2]]$parList$mean_rec_pars[1,1] + fits[[2]]$parList$Ecov_beta_R[1,1,1] * Ecov_plt)
pal <- viridis::viridis_pal(option = "H")(length(ord))
# pal <- mypalette(length(ord))
#R <- sapply(1:NROW(ab), function(x) ab[x,1]*SSB/(1 + ab[x,2]*SSB))
# Rmax <- ab[,1]*max(SSB)/(1 + ab[,2]*max(SSB))
# matplot(SSB, R, type = 'l', col = pal, lty = 1, lwd = 2)
# cbind(exp(fit1$parList$mean_rec_pars[1,1] + fit1$parList$Ecov_beta_R[1,1,1] * Ecov), R_out)

# df_R_Ecov <- cbind.data.frame(Year = fits[[2]]$years, Recruitment = fits[[2]]$rep$pred_NAA[1,1,,1], BT = Ecov, Type = "Expected")
# df_R_Ecov <- rbind(df_R_Ecov, cbind.data.frame(Year = fits[[2]]$years, Recruitment = fits[[2]]$rep$NAA[1,1,,1], BT = Ecov, Type = "Random Effect"))
# df_R_Ecov$Type <- factor(df_R_Ecov$Type)
R_df_N$name <- factor(R_df_N$name)
plt <- ggplot(R_df_N, aes(Year, val)) + 
  geom_line(data=subset(R_df_N, name == "Expected"), colour = gray(0.7), linewidth = 1.5) +  
  #scale_fill_viridis_d(begin = 0.2, end = 0.8) +
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo") +
  scale_colour_viridis_d(begin = 0.2, end = 0.8) + 
  # geom_ribbon(data = subset(R_df_N, name == "Expected"), aes(ymin=lo, ymax=hi), fill = gray(0.7, alpha = 0.3), linetype = 0) +
  geom_ribbon(data=subset(R_df_N, name == "Expected"), aes(ymin=lo, ymax=hi, fill = name), fill= gray(0.7), alpha=0.3, linetype = 0) +
  geom_point(aes(fill = BT, shape = 21), size = 4) +  
  # geom_point(data = subset(R_df_N, name == "Expected"), size = 2) +
  # geom_point(data = subset(R_df_N, name == "Random Effect"), size = 2) + 
  ylab("Recruitment (1000s)") + labs(shape="Estimate Type", colour="Bottom Temperature\nAnomaly", fill = "Estimate Type") #+
  # scale_colour_gradientn(colours = BT) + scale_fill_gradientn(colours = BT)
plt

(p_ecov_N + theme(legend.position="none") + xlab("")) + (p_ecov_S  + ylab("") + xlab(""))
