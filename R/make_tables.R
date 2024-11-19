
mod_names <- paste0("$M_{", 0:13,"}$")
col_names <- c("Model", "Temperature effect (north)", "Temperature effect (south)", "$M$ at age 1 random effects")

mod_table <- matrix(nrow = length(mod_names), ncol = length(col_names), dimnames = list(mod_names,col_names))

mod_table[,1] <- mod_names
mod_table[,2:4] <- "--"
mod_table[,4] <- rep(c("none", "time-varying"), each = 7)
mod_table[c(2,4),2] <- "Recruitment"
mod_table[c(3,4),3] <- "Recruitment"
mod_table[c(5,7),2] <- "$M$ at age 1"
mod_table[c(6,7),3] <- "$M$ at age 1"
mod_table[8:14,2:3] <- mod_table[1:7,2:3]


fits <- readRDS(file.path("results","fits_no_M_re_rev.RDS"))
fits_M_re <- readRDS(file.path("results","fits_M_re_better.RDS"))

# AIC

aic1 <- sapply(fits, function(x) {
	2*c(x$opt$obj + length(x$opt$par), sapply(x$peels, function(y) y$opt$obj + length(y$opt$par)))
})
aic2 <- sapply(fits_M_re, function(x) {
	2*c(x$opt$obj + length(x$opt$par), sapply(x$peels, function(y) y$opt$obj + length(y$opt$par)))
})

diff_aic <- apply(cbind(aic1,aic2),1, function(x) x - min(x))
aic_wts <- apply(diff_aic,2, function(x) exp(-x/2)/sum(exp(-x/2)))
rownames(diff_aic) <- rownames(aic_wts) <- mod_table[,1]
colnames(diff_aic) <- colnames(aic_wts) <- paste0("Peel ", 0:7)
saveRDS(diff_aic, file.path("results","diff_aic.RDS"))
saveRDS(aic_wts, file.path("results","aic_wts.RDS"))

#evidence for models with/without temporal variation in age 1 M.
apply(aic_wts, 2, \(x) c(sum(x[1:7]), sum(x[8:14])))
apply(diff_aic, 2, \(x) exp(-x[1:7]/2)/sapply(1:7, \(y) sum(exp(-x[c(y,y+7)]/2))))


nlminb_conv <- sapply(c(fits,fits_M_re), function(x) {
  x$opt$conv
})

nasd_conv <- sapply(c(fits,fits_M_re), function(x) {
  x$na_sdrep
})

maxgr_conv <- sapply(c(fits,fits_M_re), function(x) {
  max(abs(x$final_gradient))
})

rho_ssb <- rbind(t(sapply(fits, function(x) mohns_rho(x)$SSB)),
 t(sapply(fits_M_re, function(x) mohns_rho(x)$SSB)))
rho_Fbar <- rbind(t(sapply(fits, function(x) mohns_rho(x)$Fbar)),
 t(sapply(fits_M_re, function(x) mohns_rho(x)$Fbar)))
# rho_R <- rbind(t(sapply(fits, function(x) c(mohns_rho(x)$naa[1,1,1],mohns_rho(x)$naa[2,2,1]))),
#  t(sapply(fits_M_re, function(x) c(mohns_rho(x)$naa[1,1,1],mohns_rho(x)$naa[2,2,1]))))

rho_table <- cbind.data.frame(round(rho_ssb,3), round(rho_Fbar,3))#, round(rho_R,3))
rho_table <- cbind(mod_names, rho_table)
saveRDS(rho_table, file.path("results", "rho_table.RDS"))

########################################################
#beta and var parameter table
par_peels <- cbind(
  R_beta_N_1 = c(fits[[2]]$parList$Ecov_beta_R[1,1,1], sapply(1:7, \(x) fits[[2]]$peels[[x]]$parList$Ecov_beta_R[1,1,1])),
  R_beta_S_2 = c(fits[[3]]$parList$Ecov_beta_R[2,2,1], sapply(1:7, \(x) fits[[3]]$peels[[x]]$parList$Ecov_beta_R[2,2,1])),
  R_beta_N_3 = c(fits[[4]]$parList$Ecov_beta_R[1,1,1], sapply(1:7, \(x) fits[[4]]$peels[[x]]$parList$Ecov_beta_R[1,1,1])),
  R_beta_S_3 = c(fits[[4]]$parList$Ecov_beta_R[2,2,1], sapply(1:7, \(x) fits[[4]]$peels[[x]]$parList$Ecov_beta_R[2,2,1])),
  # M_beta_N_6 = c(fits[[7]]$parList$Ecov_beta_M[1,1,1,1,1], sapply(1:7, \(x) fits[[7]]$peels[[x]]$parList$Ecov_beta_M[1,1,1,1,1])),
  # M_beta_N_13 = c(fits_M_re[[7]]$parList$Ecov_beta_M[1,1,1,1,1], sapply(1:7, \(x) fits_M_re[[7]]$peels[[x]]$parList$Ecov_beta_M[1,1,1,1,1])),
  R_sig_N_0 = c(fits[[1]]$parList$log_NAA_sigma[1,1,1], sapply(1:7, \(x) fits[[1]]$peels[[x]]$parList$log_NAA_sigma[1,1,1])),
  R_sig_N_1 = c(fits[[2]]$parList$log_NAA_sigma[1,1,1], sapply(1:7, \(x) fits[[2]]$peels[[x]]$parList$log_NAA_sigma[1,1,1])),
  R_rho_N_0 = c(fits[[1]]$parList$trans_NAA_rho[1,1,3], sapply(1:7, \(x) fits[[1]]$peels[[x]]$parList$trans_NAA_rho[1,1,3])),
  R_rho_N_1 = c(fits[[2]]$parList$trans_NAA_rho[1,1,3], sapply(1:7, \(x) fits[[2]]$peels[[x]]$parList$trans_NAA_rho[1,1,3]))
)
par_peels[,5:6] <- exp(par_peels[,5:6])
par_peels[,7:8] <- -1 + 2/(1+exp(-par_peels[,7:8]))
par_peels <- cbind(par_peels,par_peels[,5:6]/sqrt(1 - par_peels[,7:8]^2))
colnames(par_peels)[9:10] <- c("R_marg_sig_N_0","R_marg_sig_N_1")
par_peels <- cbind(par_peels,par_peels[,6]/par_peels[,5], par_peels[,10]/par_peels[,9])
colnames(par_peels)[11:12] <- c("R_sig_N_ratio","R_marg_sig_N_ratio")

saveRDS(par_peels, file.path("results", "beta_sigma_par_peel_table.RDS"))

########################################################
#static BRP table

sapply(fits, \(x) x$rep$log_FXSPR_static)
sapply(fits, \(x) x$rep$log_SSB_FXSPR_static)
sapply(fits, \(x) x$rep$log_SPR0[,3])
exp(sapply(fits, \(x) x$rep$log_SSB_FXSPR[,3]))


fits[[1]]$rep$log_SPR0[,3])
