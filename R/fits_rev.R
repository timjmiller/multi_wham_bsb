#multi-wham (12/12/2023)
#library("wham", lib.loc = "c:/work/wham/old_packages/4a44134")

#pkgbuild::compile_dll("c:/work/wham/wham", debug = FALSE)
#pkgload::load_all("c:/work/wham/wham")

#install.packages("systemfonts", repos = "https://cloud.r-project.org", type = "binary")
#pak::pkg_install("timjmiller/wham@devel", upgrade = FALSE, ask = FALSE)

#devel branch is now multi-wham
library(wham)

library(here)

source(here("R","make_wham_inputs.R"))

input_0 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_0, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp)
# input_0$par <- fit_1$parList
# input_0$par$Ecov_beta_R[] <- 0

fit_0 <- fit_wham(input_0, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)
saveRDS(fit_0, here::here("results", "fit_0_rev.RDS"))
# fit_0 <- readRDS(here::here("results", "fit_0_rev.RDS"))
res_dir <- here::here("results","jitter", "fit_0")
dir.create(res_dir)
fit_0_jit_res <- jitter_wham(fit_RDS = here::here("results","fit_0_rev.RDS"), n_jitter = 50, res_dir = res_dir)
saveRDS(fit_0_jit_res, here::here(res_dir, "fit_0_jitter_results.RDS"))

####################
#this will give the same initial values as used by jitter_wham
# cov <- fit_0$sdrep$cov.fixed
# chol.L <- t(chol(cov))
# set.seed(8675309)
# initial_vals <- t(sapply(1:50, function(x) fit_0$opt$par + chol.L %*% cbind(rnorm(n= NCOL(cov)))))
####################

jitters <- fit_0_jit_res[[1]]
sort(sapply(jitters, \(x) x$obj))
y <- order(sapply(jitters, \(x) x$obj))

input_0_better <- fit_0$input
input_0_better$par <- fit_0$env$parList(jitters[[13]]$par)
fit_0_better <- fit_wham(input_0_better, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)
saveRDS(fit_0_better, here::here("results", "fit_0_best.RDS"))

####################
#try starting at arbitrary 0 values for age comp dispersion parameters
# input_0_alt <- input_0
# input_0_alt$par$catch_paa_pars[] <- 0
# input_0_alt$par$index_paa_pars[] <- 0
# fit_0_alt <- fit_wham(input_0_alt, do.brps = FALSE, do.sdrep = TRUE, do.osa = FALSE, do.retro = FALSE, n.newton = 0)
# saveRDS(fit_0_alt, here::here("results", "fit_0_paa_pars_init0.RDS"))
####################

fit_0 <- fit_0_better

input_1 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov_1, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp)
input_1$fleet_names <- paste0(rep(c("North_", "South_"),each = 2), temp$fleet_names)
input_1$index_names <- paste0(rep(c("North_", "South_"),c(2,2)), temp$index_names)
input_1$par <- fit_0$parList
fit_1 <- fit_wham(input_1, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)
saveRDS(fit_1, here::here("results", "fit_1.RDS"))


for(i in 2:6){
	assign(paste0("input_",i), prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = get(paste0("ecov_", i)), catch_info = catch_info,
		index_info = index_info, age_comp = age_comp))
	eval(parse(text = paste0("input_",i, "$par <- fit_0$parList")))
	assign(paste0("fit_",i), fit_wham(get(paste0("input_", i)), do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE))
}

for(i in 2:6) saveRDS(get(paste0("fit_", i)), here("results", paste0("fit_", i, ".RDS")))

fits <- lapply(0:6, function(x) get(paste0("fit_",x)))
sapply(fits, function(x) 2*(x$opt$obj + length(x$opt$par)))

fits <- lapply(fits, do_retro_peels)
fits <- lapply(fits, do_reference_points)
fits <- lapply(fits, do_sdreport)
fits <- lapply(fits, make_osa_residuals)
saveRDS(fits, here("results", paste0("fits_no_M_re_rev.RDS")))
parLists_fits <- lapply(fits, function(x) x$parList)
saveRDS(parLists_fits, here::here("results", paste0("parLists_no_M_re_rev.RDS")))

saveRDS(fits[[2]], here::here("results", "fit_1.RDS"))

res_dir <- here::here("results","jitter", "fit_1")
dir.create(res_dir)
fit_1_jit_res <- jitter_wham(fit_RDS = here::here("results","fit_1.RDS"), n_jitter = 50, res_dir = res_dir)
saveRDS(fit_1_jit_res, here::here(res_dir, "fit_1_jitter_results.RDS"))

jitters <- fit_1_jit_res[[1]]
sort(sapply(jitters, \(x) x$obj))
y <- order(sapply(jitters, \(x) x$obj))
sapply(jitters, \(x) max(abs(x$grad)))[y]

input_M_re_0 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = ecov, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp, M = M_list)
nms <- names(input_M_re_0$par)
nms <- nms[which(nms != "M_repars")]
input_M_re_0$par[nms] <- fit_0$parList[nms]
fit_M_re_0 <- fit_wham(input_M_re_0, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)

for(i in 1:6){
	assign(paste0("input_M_re_",i), prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = get(paste0("ecov_",i)), catch_info = catch_info,
	index_info = index_info, age_comp = age_comp, M = M_list))
	eval(parse(text = paste0("input_M_re_",i, "$par[nms] <- fit_0$parList[nms]")))
	assign(paste0("fit_M_re_",i), fit_wham(get(paste0("input_M_re_", i)), do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE))
}
for(i in 0:6) saveRDS(get(paste0("fit_M_re_", i)), here("results", paste0("fit_M_re_", i, ".RDS")))


fits_M_re <- lapply(0:6, function(i){
	input <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re, basic_info = basic_info, move = move, ecov = get(paste0("ecov_",i)), catch_info = catch_info,
	index_info = index_info, age_comp = age_comp, M = M_list_2)
	input$par[nms] <- fits[[1]]$parList[nms]
	return(fit_wham(input, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE))
})

fits_M_re <- lapply(fits_M_re, do_retro_peels)
fits_M_re <- lapply(fits_M_re, do_reference_points)
fits_M_re <- lapply(fits_M_re, do_sdreport)
fits_M_re <- lapply(fits_M_re, make_osa_residuals)
saveRDS(fits_M_re, here("results", paste0("fits_M_re_better.RDS")))


