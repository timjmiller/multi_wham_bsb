#like run 30, but fit models with and without temperature effects on recruitment

#multi-wham (12/12/2023)
#library("wham", lib.loc = "c:/work/wham/old_packages/4a44134")
#pkgbuild::compile_dll("c:/work/wham/wham", debug = FALSE)
#pkgload::load_all("c:/work/wham/wham")

Run34 <- readRDS(file.path("c:/work/BSB.2023.RT.Modeling", "2023.RT.Runs","Run34", "fit.RDS"))

#devel branch is now multi-wham
library(wham)

library(here)

source(here("R","make_wham_inputs.R")

# pkgload::load_all("c:/work/wham/wham")
input_1 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_1, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp)
input_1$fleet_names <- paste0(rep(c("North_", "South_"),each = 2), temp$fleet_names)
input_1$index_names <- paste0(rep(c("North_", "South_"),c(2,2)), temp$index_names)

nms <- sort(unique(names(Run34$opt$par)))
temp <- input_1
temp$par[nms] <- Run34$parList[nms]
temp$par$logit_selpars <- input_1$par$logit_selpars
temp$par$sel_repars <- input_1$par$sel_repars
temp$par$logit_selpars <- Run34$parList$logit_selpars[-(5:8),]
temp$par$sel_repars <- Run34$parList$sel_repars[-(5:8),]
temp$par$sel_repars[,2:3] <- temp$par$sel_repars[,2:3]*2 #previous version had yet to change transformation for selectivity correlation parameters
x <- temp$par$Ecov_process_pars
x[2,] <- x[2,] + 0.5*log(1 - (-1 + 2/(1+exp(-x[3,])))^2) # previous version estimated marginal sig, instead of conditional
temp$par$Ecov_process_pars <- x 


input_1 <- temp
# fit <- readRDS(here("results", "Run34_update_fit.RDS"))

fit_1 <- fit_wham(input_1, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)
fit_1 <- do_reference_points(fit_1)
fit_1 <- do_sdreport(fit_1)
fit_1 <- do_retro_peels(fit_1)

input_0 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_0, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp)
# input_0$par <- fit_1$parList
# input_0$par$Ecov_beta_R[] <- 0
fit_0 <- fit_wham(input_0, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)
#not quite as good as using better starting values
fit_0$fn(fit_1[[1]]$opt$par)
fit_0 <- do_retro_peels(fit_0)

input_2 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_2, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp)
input_2$par <- fit_0$parList
fit_2 <- fit_wham(input_2, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)

input_3 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_3, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp)
input_3$par <- fit_0$parList
fit_3 <- fit_wham(input_3, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)

input_4 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_4, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp)
input_4$par <- fit_0$parList
fit_4 <- fit_wham(input_4, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)

input_5 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_5, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp)
input_5$par <- fit_0$parList
fit_5 <- fit_wham(input_5, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)

input_6 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_6, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp)
input_6$par <- fit_0$parList
fit_6 <- fit_wham(input_6, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)

for(i in 0:6) saveRDS(get(paste0("fit_", i)), here("results", paste0("fit_", i, ".RDS")))

sapply(0:6, function(x) {
	mod <- get(paste0("fit_",x))
	2*(mod$opt$obj + length(mod$opt$par))
})

fits <- lapply(0:6, function(x) get(paste0("fit_",x)))
sapply(fits, function(x) 2*(x$opt$obj + length(x$opt$par)))

fits <- lapply(fits, do_retro_peels)
fits <- lapply(fits, do_reference_points)
fits <- lapply(fits, do_sdreport)
fits <- lapply(fits, make_osa_residuals)
saveRDS(fits, here("results", paste0("fits_no_M_re.RDS")))
parLists_fits <- sapply(fits, function(x) x$parList)
saveRDS(parLists_fits, here::here("results", paste0("parLists_no_M_re.RDS")))

input_M_re_0 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp, M = M_list)
nms <- names(input_M_re_0$par)
nms <- nms[which(nms != "M_repars")]
input_M_re_0$par[nms] <- fit_0$parList[nms]
# nofit <- fit_wham(input_M_re_0, do.fit = FALSE, do.brps = FALSE)
# x <- input_M_re_0$par$M_repars
# x[] <- input_M_re_0$map$M_repars
fit_M_re_0 <- fit_wham(input_M_re_0, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)

input_M_re_1 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_1, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp, M = M_list)
input_M_re_1$par[nms] <- fit_0$parList[nms]
fit_M_re_1 <- fit_wham(input_M_re_1, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)

input_M_re_2 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_2, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp, M = M_list)
input_M_re_2$par[nms] <- fit_0$parList[nms]
fit_M_re_2 <- fit_wham(input_M_re_2, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)

input_M_re_3 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_3, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp, M = M_list)
input_M_re_3$par[nms] <- fit_0$parList[nms]
fit_M_re_3 <- fit_wham(input_M_re_3, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)

input_M_re_4 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_4, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp, M = M_list)
input_M_re_4$par[nms] <- fit_0$parList[nms]
fit_M_re_4 <- fit_wham(input_M_re_4, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)

input_M_re_5 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_5, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp, M = M_list)
input_M_re_5$par[nms] <- fit_0$parList[nms]
fit_M_re_5 <- fit_wham(input_M_re_5, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)

input_M_re_6 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_6, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp, M = M_list)
input_M_re_6$par[nms] <- fit_0$parList[nms]
fit_M_re_6 <- fit_wham(input_M_re_6, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)

sapply(0:6, function(x) {
	mod <- get(paste0("fit_M_re_",x))
	2*(mod$opt$obj + length(mod$opt$par))
})

for(i in 0:6) saveRDS(get(paste0("fit_M_re_", i)), here("results", paste0("fit_M_re_", i, ".RDS")))
sapply(Run34$input$map, function(x) sum(!is.na(x)))
sapply(nofit$input$map, function(x) sum(!is.na(x)))



# fits <- lapply(0:6, function(i) readRDS(here("results", paste0("fit_",i,".RDS"))))

# pkgload::load_all("c:/work/wham/wham")
# inpt <- fits[[1]]$input
# inpt$par <- fits[[1]]$parList
# inpt$data$use_indices[tail(1:33,1),] <- 0
# inpt$data$use_agg_catch[tail(1:33,1),] <- 0
# inpt$data$use_catch_paa[tail(1:33,1),] <- 0
# inpt$data$use_index_paa[tail(1:33,1),] <- 0
# ind <- which(inpt$options$ecov$year >= tail(inpt$years,1))
# inpt$data$Ecov_use_obs[ind,] <- 0
# temp <- matrix(NA, 33, NCOL(inpt$par$F_pars))
# temp[1:32,] <- 1:(32*NCOL(inpt$par$F_pars))
# inpt$map$F_pars <- factor(temp)
# temp <- fit_wham(inpt, do.osa=FALSE, do.retro = FALSE, do.sdrep = FALSE, do.brps = FALSE)
# h <- stats::optimHess(temp$opt$par, temp$fn, temp$gr)

# inpt$par$F_pars[tail(1:33,7),] <- 0
# temp <- fit_wham(inpt, do.osa=FALSE, do.retro = FALSE, do.sdrep = FALSE, do.brps = FALSE)
# inpt$map <- inpt$map[names(inpt$map) != "F_pars"]
# inpt$random <- c(inpt$random, "F_pars")
# temp <- fit_wham(inpt, do.osa=FALSE, do.retro = FALSE, do.sdrep = FALSE, do.brps = FALSE)

# fits_M_re <- lapply(0:6, function(i) readRDS(here("results", paste0("fit_M_re_",i,".RDS"))))

nms <- names(fit_1$input$par)
nms <- nms[which(nms != "M_repars")]
fits_M_re <- lapply(0:6, function(i){
	input <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = get(paste0("ecov_",i)), catch_info = catch_info,
	index_info = index_info, age_comp = age_comp, M = M_list_2)
	input$par[nms] <- fits[[1]]$parList[nms]
	return(fit_wham(input, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE))
})

sapply(fits, function(x) mohns_rho(x)$SSB)
sapply(fits_M_re, function(x) mohns_rho(x)$SSB)
fits_M_re <- lapply(fits_M_re, do_retro_peels)
fits_M_re <- lapply(fits_M_re, do_reference_points)
fits_M_re <- lapply(fits_M_re, do_sdreport)
fits_M_re <- lapply(fits_M_re, make_osa_residuals)
saveRDS(fits_M_re, here("results", paste0("fits_M_re_better.RDS")))

sapply(fits_M_re, function(x) mohns_rho(x)$SSB)

# saveRDS(fits_M_re, here("results", paste0("fits_M_re_better.RDS")))
# fits_M_re <- lapply(fits_M_re, function(x) {
# 	y <- try(do_retro_peels(x))
# 	if(is.character(y)) return(x)
# 	else return(y)
# 	})
# sapply(fits_M_re, function(x) "peels" %in% names(x))

# temp <- fits_M_re[[1]]$input
# temp$data$q_upper[] <- 0.1
# q <- fits_M_re[[1]]$rep$q[1,]
# temp$par <- fits_M_re[[1]]$parList
# temp$par$logit_q <- wham:::gen.logit(q, 0, 0.1)
# temp <- fit_wham(temp, do.brps = TRUE, do.sdrep = TRUE, do.osa = FALSE, do.retro = FALSE)
# temp <- do_retro_peels(temp)
# peels <- list()
# for(i in 1:7){
# 	peels[[i]] <- try(fit_peel(i,input = temp$input, n.newton = 0))
# }
# x <- temp
# x$peels <- peels
# mohns_rho(x)$SSB


# temp2 <- fits_M_re[[1]]$input
# temp2$data$q_upper[] <- 0.1
# q <- fits_M_re[[3]]$rep$q[1,]
# temp2$par <- fits_M_re[[3]]$parList
# temp2$par$logit_q <- wham:::gen.logit(q, 0, 0.1)
# temp2 <- fit_wham(temp2, do.brps = TRUE, do.sdrep = TRUE, do.osa = FALSE, do.retro = FALSE)

gen.invlogit <- function(x,lo,hi,s){
	lo + (hi-lo)/(1+exp(-s*x))
}

gen.invlogit(-16.08 + 0.3, 0, 1000, 1)

sapply(fits_M_re, function(x) "peels" %in% names(x))


sapply(fits, function(x) {
	2*(x$opt$obj + length(x$opt$par))
})

aic1 <- sapply(fits, function(x) {
	2*c(x$opt$obj + length(x$opt$par), sapply(x$peels, function(y) y$opt$obj + length(y$opt$par)))
})
aic2 <- sapply(fits_M_re, function(x) {
	2*c(x$opt$obj + length(x$opt$par), sapply(x$peels, function(y) y$opt$obj + length(y$opt$par)))
})

aic <- t(apply(cbind(aic1,aic2),1, function(x) x - min(x)))
aic_wts <- t(round(apply(aic,1, function(x) exp(-x/2)/sum(exp(-x/2))),2))

sapply(fits_M_re, function(x) {
	2*(x$opt$obj + length(x$opt$par))
})

sapply(sort(grep("nll", names(nofit$rep), value = TRUE)), function(x) sum(nofit$rep[[x]]))
sapply(sort(grep("nll", names(Run34$rep), value = TRUE)), function(x) sum(Run34$rep[[x]]))

nms <- sort(unique(names(Run34$opt$par)))
x <- nofit$env$parList()
lapply(nms, function(i) x[[i]]- Run34$parList[[i]])

