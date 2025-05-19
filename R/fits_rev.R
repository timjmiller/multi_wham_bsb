#multi-wham (12/12/2023)
#library("wham", lib.loc = "c:/work/wham/old_packages/24c8156")

#pkgbuild::compile_dll("c:/work/wham/wham", debug = FALSE)
#pkgload::load_all("c:/work/wham/wham")

#install.packages("systemfonts", repos = "https://cloud.r-project.org", type = "binary")
#pak::pkg_install("timjmiller/wham@devel", upgrade = FALSE, ask = FALSE)

#devel branch is now multi-wham
#library(wham)
library("wham", lib.loc = "c:/work/wham/old_packages/53e236b")
# library("wham", lib.loc = "c:/work/wham/old_packages/fb8b089")

library(here)

source(here("R","make_wham_inputs.R"))

input_0 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_0, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp)
fit_0 <- fit_wham(input_0, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)
fit_0 <- do_sdreport(fit_0)
saveRDS(fit_0, here::here("results", "fit_0_rev.RDS"))

####################
#this will give the same initial values as used by jitter_wham
# cov <- fit_0$sdrep$cov.fixed
# chol.L <- t(chol(cov))
# set.seed(8675309)
# initial_vals <- t(sapply(1:50, function(x) fit_0$opt$par + chol.L %*% cbind(rnorm(n= NCOL(cov)))))
####################

#Do jitter for null model
res_dir <- here::here("results","jitter", "fit_0_rev")
dir.create(res_dir)
fit_0_jit_res <- jitter_wham(fit_RDS = here::here("results","fit_0_rev.RDS"), n_jitter = 50, res_dir = res_dir)#, do_parallel = FALSE)
saveRDS(fit_0_jit_res, here::here(res_dir, "fit_0_rev_jitter_results.RDS"))
#saveRDS(fit_0_jit_res, here::here(res_dir, "fit_0_jitter_results.RDS"))

jitters <- fit_0_jit_res[[1]]
y <- cbind(nll = sapply(jitters, \(x) x$obj), mgrad = sapply(jitters, \(x) max(abs(x$grad))))
ordy <- order(y[,1])
y[ordy,]
cbind(round(fit_0$opt$par,2), round(fit_0$opt$par - jitters[[ordy[2]]]$par,2),round(fit_0$opt$par - jitters[[ordy[3]]]$par,2),round(fit_0$opt$par - jitters[[ordy[4]]]$par,2),round(fit_0$opt$par - jitters[[ordy[5]]]$par,2))
# ordered fits 3 and 4 look identical

input_0_better <- fit_0$input
input_0_better$par <- fit_0$env$parList(jitters[[ordy[3]]]$par)
temp <- fit_wham(input_0_better, do.fit = FALSE, do.brps = FALSE)
fit_0_better <- fit_wham(input_0_better, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)
fit_0_better <- do_sdreport(fit_0_better)
saveRDS(fit_0_better, here::here("results", "fit_0_best.RDS"))

fit_0 <- fit_0_better
input_1 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_1, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp)
input_1$fleet_names <- paste0(rep(c("North_", "South_"),each = 2), temp$fleet_names)
input_1$index_names <- paste0(rep(c("North_", "South_"),c(2,2)), temp$index_names)
input_1$par <- fit_0$parList
fit_1 <- fit_wham(input_1, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)
fit_1 <- do_sdreport(fit_1)
saveRDS(fit_1, here::here("results", "fit_1.RDS"))
x <- readRDS(here::here("results", "fit_1.RDS"))
res_dir <- here::here("results","jitter", "fit_1_rev")
dir.create(res_dir)
fit_1_jit_res <- jitter_wham(fit_RDS = here::here("results","fit_1.RDS"), n_jitter = 50, res_dir = res_dir)
saveRDS(fit_1_jit_res, here::here(res_dir, "fit_1_jitter_results.RDS"))

for(i in 2:6){
	assign(paste0("input_",i), prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = get(paste0("ecov_", i)), catch_info = catch_info,
		index_info = index_info, age_comp = age_comp))
	eval(parse(text = paste0("input_",i, "$par <- fit_0$parList")))
	assign(paste0("fit_",i), fit_wham(get(paste0("input_", i)), do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE))
}

for(i in 2:6) saveRDS(get(paste0("fit_", i)), here("results", paste0("fit_", i, ".RDS")))

fits <- list(fit_0, fit_1)
fits[3:7] <- lapply(2:6, function(x) get(paste0("fit_",x)))
sapply(fits, function(x) 2*(x$opt$obj + length(x$opt$par)))

# x <- readRDS(here("results", paste0("fits_no_M_re_rev.RDS")))

fits <- lapply(fits, do_retro_peels)
fits <- lapply(fits, do_reference_points)
fits <- lapply(fits, do_sdreport)
fits <- lapply(fits, make_osa_residuals)
saveRDS(fits, here("results", paste0("fits_no_M_re_rev_1.RDS")))
saveRDS(fits, here("results", paste0("fits_no_M_re_rev.RDS")))
parLists_fits <- lapply(fits, function(x) x$parList)
saveRDS(parLists_fits, here::here("results", paste0("parLists_no_M_re_rev.RDS")))
saveRDS(fits[[2]], here::here("results", "fit_1.RDS"))

fit_0 <- readRDS(here::here("results","fit_0_best.RDS"))
input_M_re_0 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp, M = M_list)
nms <- names(input_M_re_0$par)
nms <- nms[which(nms != "M_repars")]
input_M_re_0$par[nms] <- fit_0$parList[nms]
fit_M_re_0 <- fit_wham(input_M_re_0, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)
saveRDS(get(paste0("fit_M_re_", 0)), here("results", paste0("fit_M_re_", 0, ".RDS")))

for(i in 1:6){
	assign(paste0("input_M_re_",i), prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = get(paste0("ecov_",i)), catch_info = catch_info,
	index_info = index_info, age_comp = age_comp, M = M_list))
	eval(parse(text = paste0("input_M_re_",i, "$par[nms] <- fit_0$parList[nms]")))
	assign(paste0("fit_M_re_",i), fit_wham(get(paste0("input_M_re_", i)), do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE))
	saveRDS(get(paste0("fit_M_re_", i)), here("results", paste0("fit_M_re_", i, ".RDS")))
}
for(i in 0:6) saveRDS(get(paste0("fit_M_re_", i)), here("results", paste0("fit_M_re_", i, ".RDS")))


fits_M_re <- lapply(0:6, function(i){
	input <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = get(paste0("ecov_",i)), catch_info = catch_info,
	index_info = index_info, age_comp = age_comp, M = M_list_2)
	input$par[nms] <- fit_0$parList[nms]
	return(fit_wham(input, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE))
})
#saveRDS(fits_M_re, here("results", paste0("fits_M_re_better.RDS")))
saveRDS(fits_M_re, here("results", paste0("fits_M_re_best.RDS")))

sapply(fits_M_re, \(y) y$opt$obj)
#x <- readRDS(here("results", paste0("fits_M_re_better.RDS")))
#sapply(x, \(y) y$opt$obj)

fits_M_re <- lapply(fits_M_re, do_retro_peels)
saveRDS(fits_M_re, here("results", paste0("fits_M_re_best.RDS")))
fits_M_re <- lapply(fits_M_re, do_reference_points)
saveRDS(fits_M_re, here("results", paste0("fits_M_re_best.RDS")))
fits_M_re <- lapply(fits_M_re, do_sdreport)
saveRDS(fits_M_re, here("results", paste0("fits_M_re_best.RDS")))
fits_M_re <- lapply(fits_M_re, make_osa_residuals)
saveRDS(fits_M_re, here("results", paste0("fits_M_re_best.RDS")))
#saveRDS(fits_M_re, here("results", paste0("fits_M_re_better.RDS")))


