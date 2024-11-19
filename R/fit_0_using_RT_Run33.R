Run34 <- readRDS(file.path("c:/work/BSB.2023.RT.Modeling", "2023.RT.Runs","Run34", "fit.RDS"))
Run33_no_effect <- readRDS(file.path("c:/work/BSB.2023.RT.Modeling", "2023.RT.Runs","Run33", "fit_no_effect.RDS"))

library(wham)

library(here)

source(here("R","make_wham_inputs.R"))

input_0 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_0, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp)

nms <- sort(unique(names(Run33_no_effect$opt$par)))
temp <- input_0
temp$par[nms] <- Run33_no_effect$parList[nms]
temp$par$logit_selpars <- input_0$par$logit_selpars
temp$par$sel_repars <- input_0$par$sel_repars
temp$par$logit_selpars <- Run33_no_effect$parList$logit_selpars[-(5:8),]
temp$par$sel_repars <- Run33_no_effect$parList$sel_repars[-(5:8),]
temp$par$sel_repars[,2:3] <- temp$par$sel_repars[,2:3]*2 #previous version had yet to change transformation for selectivity correlation parameters
x <- temp$par$Ecov_process_pars
x[2,] <- x[2,] + 0.5*log(1 - (-1 + 2/(1+exp(-x[3,])))^2) # previous version estimated marginal sig, instead of conditional
temp$par$Ecov_process_pars <- x 
# temp$par$Ecov_beta_R[] <- 0
input_0_Run33 <- temp
fit_0_Run33 <- fit_wham(input_0_Run33, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)
saveRDS(fit_0_Run33, here::here("results", "fit_0_Run33.RDS"))
