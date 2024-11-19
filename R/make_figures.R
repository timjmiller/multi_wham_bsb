library(wham)
library(ggplot2)
library(dplyr)
source(here::here("R","plot_functions.R"))

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
plot.ecov(fits[[1]])
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

df <- cbind(df, rel_M1_est = NA)
for(i in unique(df$model)) for(j in unique(df$type)) for(k in unique(df$stock)) {
  df$rel_M1_est[which(df$model ==i & df$type == j & df$stock == k)] <- df$est[which(df$model ==i & df$type == j & df$stock == k)]/df$est[which(df$model ==levels(df$model)[2] & df$type == j & df$stock == k)]
}


library(latex2exp) #https://github.com/stefano-meschiari/latex2exp
plt <- ggplot(subset(df), aes(x = year, y = rel_M1_est, colour = model)) + scale_colour_viridis_d() + 
    geom_point(size = 1) + geom_line() +
    facet_grid(type ~ stock, scales = "free_y", switch = "y") + theme_bw() +
    xlab("Year") + ylab(NULL) + labs(colour = "Model") +
    scale_color_discrete(labels=lapply(levels(df$model), TeX)) + 
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
pal <- viridis::viridis_pal(option = "H")(length(ord))
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

df <- cbind.data.frame(Year = fit1$years, Recruitment = fit1$rep$pred_NAA[1,1,,1], BT = Ecov, Type = "Expected")
df <- rbind(df, cbind.data.frame(Year = fit1$years, Recruitment = fit1$rep$NAA[1,1,,1], BT = Ecov, Type = "Random Effect"))
df$Type <- factor(df$Type)
plt <- ggplot(df, aes(Year, Recruitment, colour = BT, shape = Type)) + theme_bw() + 
  geom_line(data = subset(df,  Type == "Expected"), colour = gray(0.7, alpha = 0.4), size = 2) +
  geom_point(data = subset(df, Type == "Expected"), size = 2) +
  geom_point(data = subset(df, Type == "Random Effect"), size = 2) + 
  ylab("Recruitment (1000s)") + labs(shape="Estimate Type", colour="Bottom Temperature\nAnomaly") + 
  theme(axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.title = element_text(size = rel(2), hjust=0.5), 
    legend.box.just = "center",legend.text = element_text(size = rel(1.5))) +
  scale_colour_gradientn(colours = pal) + scale_fill_gradientn(colours = pal)
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

proj <- list() 
for(i in 1:3) proj[[i]] <- readRDS(here::here("results",paste0("m1_proj_",i,".RDS")))
fjp <- list() 
for(i in 1:3) fjp[[i]] <- readRDS(here::here("results",paste0("m1_proj_",i,"_sdrep_fjp.RDS")))

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

temp$R1 <- proj[[1]]$rep$NAA[1,1,,1]
temp$predR1 <- proj[[1]]$rep$pred_NAA[1,1,,1]
temp$R1.cv <- TMB:::as.list.sdreport(proj[[1]]$sdrep, report = TRUE, what = "Std")$log_NAA_rep[1,1,,1]
temp$R2 <- proj[[2]]$rep$NAA[1,1,,1]
temp$predR2 <- proj[[2]]$rep$pred_NAA[1,1,,1]
temp$R2.cv <- TMB:::as.list.sdreport(proj[[2]]$sdrep, report = TRUE, what = "Std")$log_NAA_rep[1,1,,1]
temp$R3 <- proj[[3]]$rep$NAA[1,1,,1]
temp$predR3 <- proj[[3]]$rep$pred_NAA[1,1,,1]
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
temp <- temp[which(proj[[2]]$input$years_full<2032),]

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
cols <- viridis::viridis_pal(option = "H")(3)
poly_cols <- viridis::viridis_pal(option = "H", alpha = 0.3)(3)
par(mfcol = c(2,2), oma = c(5,1,1,1), mar = c(1,5,1,1))
ylim <- range(subset(temp, year %in% 2019:2031)[c("lo1","lo2","lo3", "hi1", "hi2", "hi3")])
plot(temp$year, temp$est1, type = "n", xlim = c(2018,2031), lwd = 2, ylab = "Bottom Temperature Anomaly", ylim = ylim, xaxt = "n", cex.axis = 2, cex.lab = 2)
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
plot(temp$year, temp$R1, type = "n", xlim = c(2018,2031), lwd = 2, xlab = "", ylab = "Recruitment (1000s)", ylim = ylim, cex.axis = 2, cex.lab = 2)
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
plot(temp$year, temp$R1.cv, type = "n", xlim = c(2018,2031), lwd = 2, xlab = "", ylab = "CV(Recruitment)", ylim = ylim, xaxt = "n", cex.axis = 2, cex.lab = 2)
axis(1, labels = FALSE)
grid(lty =2, col = gray(0.7))
lines(temp$year, temp$R1.cv, lwd = 2, col = cols[1])
proj_ind <- which(df$year>2020)
lines(df$year[proj_ind], df$R2.cv[proj_ind], col = cols[2], lwd = 2)
lines(df$year[proj_ind], df$R3.cv[proj_ind], col = cols[3], lwd = 2)
abline(v=2021, lty = 3, lwd = 2)

ylim <- c(0,max(subset(temp, year %in% 2019:2031)[c("predR1.cv","predR2.cv","predR3.cv")]))
plot(temp$year, temp$predR1.cv, type = "n", xlim = c(2018,2031), lwd = 2, xlab = "", ylab = "CV(Exp. Recruitment)", ylim = ylim, cex.axis = 2, cex.lab = 2)
grid(lty =2, col = gray(0.7))
lines(temp$year, temp$predR1.cv, lwd = 2, col = cols[1])
proj_ind <- which(df$year>2020)
lines(df$year[proj_ind], df$predR2.cv[proj_ind], col = cols[2], lwd = 2)
lines(df$year[proj_ind], df$predR3.cv[proj_ind], col = cols[3], lwd = 2)
abline(v=2021, lty = 3, lwd = 2)

mtext(side = 1, line = 2, "Year", cex = 2, outer = TRUE)
dev.off()

#plot how R and pred R come together in projection period
#plot different projections of pred R
#plot time series of pred R and Ecov with color gradient f(Ecov)

########################################

#plot selectivity at time and age for blocks with RE.
fit <- readRDS(file.path("results","fits_no_M_re_rev.RDS"))[[2]]
source(here::here("R","plot_functions.R"))
cairo_pdf(file.path("paper", "selectivity_re_plot.pdf"), width = 16, height = 12)
plot.selectivity(fit, blocks = c(1:2,5:6), block.names = c("North Commercial Fleet", "North Recreational Fleet", "North Recreational Catch/Effort Index", "North Spring Index"))
dev.off()
#reference points

