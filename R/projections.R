library(wham)
fit1 <-readRDS(file.path("results","fits_no_M_re_rev.RDS"))[[2]]

tdat <- cbind.data.frame(time = fit1$input$years_Ecov, bt_N = fit1$rep$Ecov_x[,1],bt_S = fit1$rep$Ecov_x[,2])
lmtdat <- lm(bt_N ~time, data = tdat)
newdat <- cbind.data.frame(time = c(tdat$time,tail(tdat$time,1)+1:20), bt_N = c(tdat$bt_N, rep(NA,20)), bt_S = c(tdat$bt_S, rep(NA,20)))
predlm <- predict(lmtdat, newdata = newdat, se.fit = TRUE)
newdat$pred_N <- predlm$fit
newdat$pred_lo_N <- newdat$pred_N + qnorm(0.025)*predlm$se.fit
newdat$pred_hi_N <- newdat$pred_N+ qnorm(0.975)*predlm$se.fit

plot(newdat$time, newdat$bt_N)
lines(newdat$time, newdat$pred_N)
polygon(c(newdat$time,rev(newdat$time)), c(newdat$pred_lo_N, rev(newdat$pred_hi_N)), border = "transparent", col = gray(0.7,alpha=0.3))

lmtdat <- lm(bt_S ~time, data = tdat)
predlm <- predict(lmtdat, newdata = newdat, se.fit = TRUE)
newdat$pred_S <- predlm$fit
newdat$pred_lo_S <- newdat$pred_S + qnorm(0.025)*predlm$se.fit
newdat$pred_hi_S <- newdat$pred_S + qnorm(0.975)*predlm$se.fit

plot(newdat$time, newdat$bt_S)
lines(newdat$time, newdat$pred_S)
polygon(c(newdat$time,rev(newdat$time)), c(newdat$pred_lo_S, rev(newdat$pred_hi_S)), border = "transparent", col = gray(0.7,alpha=0.3))

proj.opts.list <- list()
# proj.opts.list[[1]] <- list(n.yrs = 11, proj_F_opt = c(5,rep(3,10)), proj_Fcatch = rep(10000,11)) #continue Ecov process by default
# proj.opts.list[[2]] <- list(n.yrs = 11, proj_F_opt = c(5,rep(3,10)), proj_Fcatch = rep(10000,11),  cont.ecov=FALSE, use.last.ecov=FALSE, avg.ecov.yrs=tail(fit1$years,5)) #use average of last 5 years
# proj.opts.list[[3]] <- list(n.yrs = 11, proj_F_opt = c(5,rep(3,10)), proj_Fcatch = rep(10000,11), cont.ecov=FALSE, use.last.ecov=FALSE, 
# 	proj.ecov = cbind(newdat$pred_N,newdat$pred_S)[65:74,]) #use prediction from linear trend
proj.opts.list[[1]] <- list(n.yrs = 10, proj_F_opt = rep(3,10), proj_Fcatch = rep(10000,10)) #continue Ecov process by default
proj.opts.list[[2]] <- list(n.yrs = 10, proj_F_opt = rep(3,10), proj_Fcatch = rep(10000,10),  cont.ecov=FALSE, use.last.ecov=FALSE, avg.ecov.yrs=tail(fit1$years,5)) #use average of last 5 years
proj.opts.list[[3]] <- list(n.yrs = 10, proj_F_opt = rep(3,10), proj_Fcatch = rep(10000,10), cont.ecov=FALSE, use.last.ecov=FALSE, 
	proj.ecov = cbind(newdat$pred_N,newdat$pred_S)[65:74,]) #use prediction from linear trend
saveRDS(proj.opts.list, here::here("results", "m1_proj_opts_list.RDS"))
# temp <- project_wham(fit1, proj.opts = proj.opts.list[[1]], check.version = F, save.sdrep = FALSE)
proj.list <- list()
proj.list[[1]] <- project_wham(fit1, proj.opts = proj.opts.list[[1]], check.version = F, do.sdrep = FALSE)
saveRDS(proj.list[[1]], here::here("results","m1_proj_1.RDS"))
proj.list[[2]] <- project_wham(fit1, proj.opts = proj.opts.list[[2]], check.version = F, do.sdrep = FALSE)
saveRDS(proj.list[[2]], here::here("results","m1_proj_2.RDS"))
proj.list[[3]] <- project_wham(fit1, proj.opts = proj.opts.list[[3]], check.version = F, do.sdrep = FALSE)
saveRDS(proj.list[[3]], here::here("results","m1_proj_3.RDS"))

# temp <- readRDS(here::here("results","m1_proj_1.RDS"))

# plot(proj.list[[3]]$years_full, proj.list[[3]]$rep$NAA[1,1,,1], type = "l")
# lines(proj.list[[3]]$years_full, proj.list[[3]]$rep$pred_NAA[1,1,,1], col = "red")
# lines(proj.list[[1]]$years_full, proj.list[[1]]$rep$pred_NAA[1,1,,1], col = "blue")
# proj.list[[3]]$rep$Ecov_out_R[1,,1]
# fit1$rep$Ecov_out_R[1,,1]

full_JP1 <- TMB::sdreport(proj.list[[1]], getJointPrecision = TRUE)
saveRDS(full_JP1, here::here("results", "m1_proj_1_sdrep_fjp.RDS"))
full_JP2 <- TMB::sdreport(proj.list[[2]], getJointPrecision = TRUE)
saveRDS(full_JP2, here::here("results", "m1_proj_2_sdrep_fjp.RDS"))
full_JP3 <- TMB::sdreport(proj.list[[3]], getJointPrecision = TRUE)
saveRDS(full_JP3, here::here("results", "m1_proj_3_sdrep_fjp.RDS"))

# temp <- TMB:::as.list.sdreport(fit1$sdrep, report = TRUE, what = "Est")$Ecov_x

# ind <- which(names(fit1$sdrep$value) == "Ecov_x")
# temp - matrix(fit1$sdrep$value[ind], length(fit1$input$years_Ecov), 2)
# ind <- matrix(ind,length(fit1$input$years_Ecov), 2)[,1]
# Sig <- fit1$sdrep$cov[ind,ind]
# obj <- function(pars, x, Sig){
# 	mu = pars[1] + pars[2]*fit1$input$years_Ecov
# 	res <- x - mu
# 	return(-mvtnorm::dmvnorm(res,sigma = Sig, log = TRUE))
# }
# opt <- nlminb(c(0,0), obj, x = temp[,1], Sig = Sig)
# opt$hess <- optimHess(opt$par, fn = obj, x = temp[,1], Sig = Sig)

####################################
#change usage of recruitment for brps
library(wham)
fit1 <-readRDS(file.path("results","fits_no_M_re_rev.RDS"))[[2]]

input <- fit1$input
input$data$XSPR_R_opt <- 3 #use annual expected recruitment
input$par <- fit1$parList
temp <- fit_wham(input, do.fit = FALSE)
temp$opt <- fit1$opt #fake it to do sdreport on projection
saveRDS(temp, here::here("results", "m1_XSPR_opt_3.RDS"))
# temp <- readRDS(here::here("results", "m1_XSPR_opt_3.RDS"))
x <- TMB::sdreport(temp)
saveRDS(x,here::here("results", "m1_XSPR_opt_3_sdrep.RDS"))
proj.opts.list <- readRDS(here::here("results", "m1_proj_opts_list.RDS"))
proj_input_1 <- prepare_projection(temp, proj.opts = proj.opts.list[[1]])

# proj_1 <- TMB::MakeADFun(proj_input_1$data, proj_input_1$par, random = proj_input_1$random, map = proj_input_1$map)
# x <- proj_1$report()
# proj_1$fn()

# proj.opts.list[[3]]$n.yrs <- 10
# proj.opts.list[[3]]$proj_Fcatch <- rep(1000,10)
# proj.opts.list[[3]]$proj_F_opt <- c(5,rep(3,9))
# proj_input_3 <- prepare_projection(temp, proj.opts = proj.opts.list[[3]])

# pkgload::load_all("~/wham_research/wham")
# proj_3 <- TMB::MakeADFun(proj_input_3$data, proj_input_3$par, random = proj_input_3$random, map = proj_input_3$map)
# proj_3$fn()
# x <- proj_1$report()

proj.list <- list()
proj.list[[1]] <- project_wham(temp, proj.opts = proj.opts.list[[1]], check.version = FALSE, do.sdrep = FALSE)

saveRDS(proj.list[[1]], here::here("results","m1_proj_1_R_opt3.RDS"))
proj.list[[2]] <- project_wham(temp, proj.opts = proj.opts.list[[2]], check.version = F, do.sdrep = FALSE)
saveRDS(proj.list[[2]], here::here("results","m1_proj_2_R_opt3.RDS"))
proj.list[[3]] <- project_wham(temp, proj.opts = proj.opts.list[[3]], check.version = F, do.sdrep = FALSE)
saveRDS(proj.list[[3]], here::here("results","m1_proj_3_R_opt3.RDS"))

temp <- readRDS(here::here("results", "m1_XSPR_opt_3.RDS"))
input <- temp$input
input$data$XSPR_R_opt <- 1 #use annual random effects
temp2 <- fit_wham(input, do.fit = FALSE)
temp2$opt <- temp$opt #fake it to do sdreport on projection
saveRDS(temp2, here::here("results", "m1_XSPR_opt_1.RDS"))
# temp2 <- readRDS(here::here("results", "m1_XSPR_opt_3.RDS"))
proj.opts.list <- readRDS(here::here("results", "m1_proj_opts_list.RDS"))
proj.list <- list()
proj.list[[1]] <- project_wham(temp2, proj.opts = proj.opts.list[[1]], check.version = FALSE)
saveRDS(proj.list[[1]], here::here("results","m1_proj_1_R_opt1.RDS"))
proj.list[[2]] <- project_wham(temp2, proj.opts = proj.opts.list[[2]], check.version = F)
saveRDS(proj.list[[2]], here::here("results","m1_proj_2_R_opt1.RDS"))
proj.list[[3]] <- project_wham(temp2, proj.opts = proj.opts.list[[3]], check.version = F)
saveRDS(proj.list[[3]], here::here("results","m1_proj_3_R_opt1.RDS"))
