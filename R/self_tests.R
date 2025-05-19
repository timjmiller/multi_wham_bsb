
# library(wham)
library("wham", lib.loc = "c:/work/wham/old_packages/53e236b")
#seeds <- 1:10
fit_1 <- readRDS(here::here("results","fit_1.RDS"))
# seeds <- sample(0:1e9, 100, replace = TRUE)
res_dir <- here::here("results","self_test", "fit_1_rev")
dir.create(res_dir)
# self_test_1 <- self_test(fit_RDS = here::here("results","fit_1.RDS"), n = 50, seeds = seeds, res_dir = res_dir)
self_test_1 <- self_test(fit_RDS = here::here("results","fit_1.RDS"), n = 50, res_dir = res_dir)
saveRDS(self_test_1, file = file.path(res_dir, "self_test_results.RDS"))

self_test_1 <- readRDS(file = file.path(res_dir, "self_test_results.RDS"))

test_files <- dir(here::here("results","self_test","fit_1"), full.names = TRUE, pattern = glob2rx("cond_sim*.RDS"))
sims <- lapply(test_files, \(x) readRDS(x))

# sims <- self_test_1[[1]]

sd_scale_pars <- t(sapply(sims, \(x) x$par[names(x$par) =="log_index_sig_scale"]))

# use <- which(sapply(sims, \(x) !(is.na(max(abs(x$grad))) | max(abs(x$grad))>0.01)))
use <- 1:length(sims)

par(mfrow= c(2,2))
resid_SSB_N <- sapply(sims[use], \(x) x$SSB[,1]/fit_1$rep$SSB[,1] - 1)
med_res_SSB_N <- apply(resid_SSB_N, 1, median, na.rm = TRUE)
plot(fit_1$years, med_res_SSB_N)
abline(h = 0)
abline(h = median(med_res_SSB_N, na.rm = TRUE), lty =2)
abline(h = median(resid_SSB_N, na.rm = TRUE), lty =2, col = "red")

resid_SSB_S <- sapply(sims[use], \(x) x$SSB[,2]/fit_1$rep$SSB[,2] - 1)
med_res_SSB_S <- apply(resid_SSB_S, 1, median, na.rm = TRUE)
plot(fit_1$years, med_res_SSB_S)
abline(h = 0)
abline(h = median(med_res_SSB_S, na.rm = TRUE), lty =2)
abline(h = median(resid_SSB_S, na.rm = TRUE), lty =2, col = "red")

sapply(sims, \(x) max(abs(x$grad)))

resid_par <- sapply(sims, \(x) x$par - fit_1$opt$par)
med_res_par <- apply(resid_par, 1, median)
plot(med_res_par, xlab = "par")
abline(h = 0)

#seeds <- sample(0:1e9, 100, replace = TRUE)
seeds <- self_test_1$seeds
res_dir <- here::here("results","self_test", "fit_1_rev_fix")
dir.create(res_dir)
fit_1 <- readRDS(here::here("results","fit_1.RDS"))
source(here::here("R","plot_functions.R"))
map_change <- list(log_index_sig_scale = factor(rep(NA, length(fit_1$parList$log_index_sig_scale))))
self_test_fix <- self_test(fit_RDS = here::here("results","fit_1.RDS"), n = 50, seeds = seeds, res_dir = res_dir, map_change = map_change)
saveRDS(self_test_fix, file = file.path(res_dir, "self_test_results.RDS"))
# self_test_fix <- readRDS(file = file.path(res_dir, "self_test_results.RDS"))
sapply(self_test_fix[[1]], \(x) x$opt$conv)
sapply(self_test_fix[[1]], \(x) max(abs(x$grad)))
sapply(self_test_1[[1]], \(x) max(abs(x$grad)))
# self_test_fix <- readRDS(here::here("results","self_test", "fit_1_fix", "self_test_results.RDS"))

resid_SSB_N <- sapply(self_test_fix[[1]], \(x) x$SSB[,1]/fit_1$rep$SSB[,1] - 1)
med_res_SSB_N <- apply(resid_SSB_N, 1, median, na.rm = TRUE)
plot(fit_1$years, med_res_SSB_N)
abline(h = 0)
abline(h = median(med_res_SSB_N, na.rm = TRUE), lty =2)
abline(h = median(resid_SSB_N, na.rm = TRUE), lty =2, col = "red")

resid_SSB_S <- sapply(self_test_fix[[1]], \(x) x$SSB[,2]/fit_1$rep$SSB[,2] - 1)
med_res_SSB_S <- apply(resid_SSB_S, 1, median, na.rm = TRUE)
plot(fit_1$years, med_res_SSB_S)
abline(h = 0)
abline(h = median(med_res_SSB_S, na.rm = TRUE), lty =2)
abline(h = median(resid_SSB_S, na.rm = TRUE), lty =2, col = "red")


#median across years for each sim
med_sim_resid_SSB_S <- sapply(self_test_fix[[1]], \(x) median(x$SSB[,2]/fit_1$rep$SSB[,2] - 1, na.rm = TRUE))

omit <- which(names(fit_1$env$last.par.best) == "log_index_sig_scale")

resid_last_par <- sapply(self_test_fix[[1]], \(x) x$last.par.best -fit_1$env$last.par.best[-omit]) #DO NOT DO relative bias for these because parameters can be negative.
temp <- cor(cbind(med_sim_resid_SSB_S, t(resid_last_par))[-28,])
which(abs(temp[-1,1])>0.7) #5 and 401 (logit_q and log_NAA)
temp[rownames(temp) == "logit_q",1] 
temp[-1,1][c(5,401)]

temp <- cor(cbind(t(resid_SSB_S), t(resid_last_par))[-28,])
which(sapply(1:33, \(x) any(abs(temp[-c(1:33, which(rownames(temp) == "log_NAA")),x])>0.7)))
sapply(1:33, \(x) which(abs(temp[-c(1:33, which(rownames(temp) == "log_NAA")),x])>0.7))

which(abs(temp[-1,1])>0.7) #5 and 401 (logit_q and log_NAA)
temp[rownames(temp) == "logit_q",1] 
temp[-1,1][c(5,401)]

temp <- sapply(self_test_fix[[1]], \(x) x$parList$mu_prior_re[1,1,,] - fit_1$parList$mu_prior_re[1,1,,])
apply(temp,1,range, na.rm = TRUE)
fit_1$env$last.par.best[ind]

temp <- sapply(self_test_fix[[1]], \(x) x$rep$mu[1,8,1,1,2,1] - fit_1$rep$mu[1,8,1,1,2,1])
self_test_fix[[1]][[1]]$sim_input$data$mu[1,8,1,1,,]
self_test_fix[[1]][[1]]$rep$mu[1,8,1,1,,]
fit_1$rep$mu[1,8,1,1,,]

# ini <- self_test_fix[[1]][[38]]$sim_input #bad negative median across all years
# ini$random <- c(ini$random, c("mean_rec_pars","log_N1", "logit_q","F_pars","logit_selpars","Ecov_beta_R"))
# reml_fit <- fit_wham(ini, do.sdrep = FALSE, do.retro = FALSE, do.osa = FALSE, do.brps = FALSE)

# #reml provides lower bias for this sim
# median(reml_fit$rep$SSB[,2]/fit_1$rep$SSB[,2]-1) 
# median(self_test_fix[[1]][[38]]$rep$SSB[,2]/fit_1$rep$SSB[,2]-1)
# median(reml_fit$rep$SSB[,1]/fit_1$rep$SSB[,1]-1) 
# median(self_test_fix[[1]][[38]]$rep$SSB[,1]/fit_1$rep$SSB[,1]-1)
 
# reml_fit
# unique(names(self_test_fix[[1]][[1]]$par))

# res_dir <- here::here("results","self_test", "fit_1_fix_reml")
# dir.create(res_dir)
# fit_1 <- readRDS(here::here("results","fit_1.RDS"))
# source(here::here("R","plot_functions.R"))
# map_change <- list(log_index_sig_scale = factor(rep(NA, length(fit_1$parList$log_index_sig_scale))))
# self_test_reml <- self_test(fit_RDS = here::here("results","fit_1.RDS"), n = 50, seeds = seeds, res_dir = res_dir, map_change = map_change, reml = TRUE)
# saveRDS(self_test_reml, file = file.path(res_dir, "self_test_results.RDS"))
# #self_test_reml <- readRDS(here::here("results","self_test", "fit_1_fix_reml", "self_test_results.RDS"))
# test_files <- dir(here::here("results","self_test", "fit_1_fix_reml"), full.names = TRUE, pattern = glob2rx("cond_sim*.RDS"))
# sims <- lapply(test_files, \(x) readRDS(x))
# seeds_done <- sapply(test_files, \(x) readRDS(x)$seed)
# which_seeds_done <- match(seeds_done, seeds)

# which_seeds_left <- which(!seeds %in% seeds_done)
# x <- self_test_i(which_seeds_left[1], fit_RDS = here::here("results","fit_1.RDS"), seeds = seeds, map_change = map_change, res_dir = res_dir, reml = TRUE, MakeADFun.silent = FALSE)
# x <- self_test_i(which_seeds_left[2], fit_RDS = here::here("results","fit_1.RDS"), seeds = seeds, map_change = map_change, res_dir = res_dir, reml = TRUE, MakeADFun.silent = FALSE)
# self_test_left <- self_test(fit_RDS = here::here("results","fit_1.RDS"), which_seeds = which_seeds_left, seeds = seeds, res_dir = res_dir, map_change = map_change, reml = TRUE)
# seeds_left <- sapply(self_test_left[[1]], \(x) x$seed)
# self_test_reml <- list(self_test_results = list(), seeds = self_test_left$seeeds)
# self_test_reml[[1]][match(seeds_left, seeds)] <- self_test_left[[1]]
# sims <- lapply(test_files, \(x) readRDS(x))
# sim_seeds <- sapply(sims, \(x) x$seed)
# sims_done_before <- sims[!sim_seeds %in% seeds_left]
# seeds_done_before <- sapply(sims_done_before, \(x) x$seed)
# self_test_reml[[1]][match(seeds_done_before,seeds)] <- sims_done_before
# sapply(self_test_reml[[1]], \(x) x$seed) - seeds
# sapply(self_test_reml[[1]], \(x) length(x$opt$par))
# saveRDS(self_test_reml, file = file.path(res_dir, "self_test_results.RDS"))

# source(here::here("R","plot_functions.R"))
# plot_SSB(self_test_reml[[1]],fit = fit_1)
# x <- self_test_reml[[1]]

# temp <- SSB_res(sims,fit_1)
# #seeds <- sapply(self_test_reml[[1]], \(x) x$seed)
# #seeds_2 <- sapply(self_test_fix[[1]], \(x) x$seed)


# ##########################################
# #low age comp observation error self test
# res_dir <- here::here("results","self_test", "fit_1_low_error")
# dir.create(res_dir)
# fit_1_le <- readRDS(here::here("results","fit_1.RDS"))
# fit_1_le$parList$catch_paa_pars[2:4,1] <- log(0.1)
# # fit_1_le$parList$catch_paa_pars[1,1] <- log(5)
# fit_1_le$parList$catch_paa_pars[2:4,2] <- 0
# fit_1_le$parList$index_paa_pars[c(1,3:4),1] <- log(0.1)
# # fit_1_le$parList$index_paa_pars[2,1] <- log(5)
# fit_1_le$parList$index_paa_pars[2,1] <- log(10)
# fit_1_le$parList$index_paa_pars[c(1,3:4),2] <- 0
# # fit_1_le$map$catch_paa_pars <- factor(rep(NA,length(fit_1_le$parList$catch_paa_pars)))
# # fit_1_le$map$index_paa_pars <- factor(rep(NA,length(fit_1_le$parList$index_paa_pars)))
# saveRDS(fit_1_le, file = file.path(res_dir, "fit_1_low_error.RDS"))
# self_test_le <- self_test(fit_RDS = file.path(res_dir, "fit_1_low_error.RDS"), n = 50, seeds = seeds, res_dir = res_dir, map_change = map_change)
# saveRDS(self_test_le, file = file.path(res_dir, "self_test_results.RDS"))
# # test_files <- dir(res_dir, full.names = TRUE, pattern = glob2rx("cond_sim*.RDS"))
# # sims <- lapply(test_files, \(x) readRDS(x))
# # self_test_le <- readRDS(here::here("results","self_test", "fit_1_low_error", "fit_1_low_error.RDS"))
# sims <- self_test_le[[1]]
# sapply(sims, \(x) max(abs(x$grad)))
# apply(sapply(sims, \(x) x$parList$catch_paa_pars[,1]),1,range)

# source(here::here("R","plot_functions.R"))
# x <- SSB_res(sims, fit_1_le)

##########################################

# self_test_reml <- readRDS(here::here("results","self_test", "fit_1_fix_reml", "self_test_results.RDS"))
# seeds <- sapply(self_test_reml[[1]], \(x) x$seed)
# self_test_ml <- readRDS(here::here("results","self_test", "fit_1_fix", "self_test_results.RDS"))
# par(mfrow = c(1,2))

# reml_resid_SSB_S <- sapply(self_test_reml[[1]], \(x) x$SSB[,2]/fit_1$rep$SSB[,2] - 1)
# ml_resid_SSB_S <- sapply(self_test_ml[[1]], \(x) x$SSB[,2]/fit_1$rep$SSB[,2] - 1)
# med_reml_res_SSB_S <- apply(reml_resid_SSB_S, 1, median, na.rm = TRUE)
# med_ml_res_SSB_S <- apply(ml_resid_SSB_S, 1, median, na.rm = TRUE)
# plot(fit_1$years, med_reml_res_SSB_S, ylim = range(c(med_reml_res_SSB_S, med_ml_res_SSB_S), na.rm =T), ylab = "Rel Diff.")
# mtext(side = 3, "REML")
# abline(h = 0)
# #median of annual medians (n_years)
# abline(h = median(med_reml_res_SSB_S, na.rm = TRUE), lty =2)
# #median of all resids (n_sims x n_years)
# abline(h = median(reml_resid_SSB_S, na.rm = TRUE), lty =2, col = "red")

# plot(fit_1$years, med_ml_res_SSB_S, ylim = range(c(med_reml_res_SSB_S, med_ml_res_SSB_S), na.rm =T), ylab = "")
# mtext(side = 3, "ML")
# abline(h = 0)
# #median of annual medians (n_years)
# abline(h = median(med_ml_res_SSB_S, na.rm = TRUE), lty =2)
# #median of all resids (n_sims x n_years)
# abline(h = median(ml_resid_SSB_S, na.rm = TRUE), lty =2, col = "red")

# par(mfrow = c(1,2))
# reml_resid_SSB_N <- sapply(self_test_reml[[1]], \(x) x$SSB[,1]/fit_1$rep$SSB[,1] - 1)
# ml_resid_SSB_N <- sapply(self_test_ml[[1]], \(x) x$SSB[,1]/fit_1$rep$SSB[,1] - 1)
# med_reml_res_SSB_N <- apply(reml_resid_SSB_N, 1, median, na.rm = TRUE)
# med_ml_res_SSB_N <- apply(ml_resid_SSB_N, 1, median, na.rm = TRUE)
# plot(fit_1$years, med_reml_res_SSB_N, ylim = range(c(med_reml_res_SSB_N, med_ml_res_SSB_N), na.rm =T), ylab = "Rel Diff.")
# mtext(side = 3, "REML")
# abline(h = 0)
# #median of annual medians (n_years)
# abline(h = median(med_reml_res_SSB_N, na.rm = TRUE), lty =2)
# #median of all resids (n_sims x n_years)
# abline(h = median(reml_resid_SSB_N, na.rm = TRUE), lty =2, col = "red")

# plot(fit_1$years, med_ml_res_SSB_N, ylim = range(c(med_reml_res_SSB_N, med_ml_res_SSB_N), na.rm =T), ylab = "")
# mtext(side = 3, "ML")
# abline(h = 0)
# #median of annual medians (n_years)
# abline(h = median(med_ml_res_SSB_N, na.rm = TRUE), lty =2)
# #median of all resids (n_sims x n_years)
# abline(h = median(ml_resid_SSB_N, na.rm = TRUE), lty =2, col = "red")

# fit_1$rep$mu[1,8,1,1,,]
# fit_1$rep$mu[1,8,11,1,,]
# ml_resid_mu_SN <- sapply(self_test_ml[[1]], \(x) x$rep$mu[1,8,1,1,2,1] - fit_1$rep$mu[1,8,1,1,2,1])
# ml_resid_mu_NS <- sapply(self_test_ml[[1]], \(x) x$rep$mu[1,8,11,1,1,2] - fit_1$rep$mu[1,8,11,1,1,2])
# med_ml_res_mu_SN <- apply(ml_resid_mu_SN, 1, median, na.rm = TRUE)

