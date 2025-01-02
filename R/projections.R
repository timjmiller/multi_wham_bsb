# library(wham)
devtools::load_all("c:/work/wham/wham")
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

#Use terminal F in projections
proj.opts.list <- list()
proj.opts.list[[1]] <- list(n.yrs = 10, proj_F_opt = rep(1,10), proj_Fcatch = rep(10000,10)) #continue Ecov process by default
proj.opts.list[[2]] <- list(n.yrs = 10, proj_F_opt = rep(1,10), proj_Fcatch = rep(10000,10),  cont.ecov=FALSE, use.last.ecov=FALSE, avg.ecov.yrs=tail(fit1$years,5)) #use average of last 5 years
proj.opts.list[[3]] <- list(n.yrs = 10, proj_F_opt = rep(1,10), proj_Fcatch = rep(10000,10), cont.ecov=FALSE, use.last.ecov=FALSE, 
	proj.ecov = cbind(newdat$pred_N,newdat$pred_S)[65:74,]) #use prediction from linear trend
saveRDS(proj.opts.list, here::here("results", "m1_proj_opts_list.RDS"))
# temp <- project_wham(fit1, proj.opts = proj.opts.list[[1]], check.version = F, save.sdrep = FALSE)
proj.list <- list()
proj.list[[1]] <- project_wham(fit1, proj.opts = proj.opts.list[[1]], check.version = F, do.sdrep = TRUE, TMB.jointPrecision = TRUE)
saveRDS(proj.list[[1]], here::here("results","m1_proj_1_R_opt2.RDS"))
proj.list[[2]] <- project_wham(fit1, proj.opts = proj.opts.list[[2]], check.version = F, do.sdrep = TRUE, TMB.jointPrecision = TRUE)
saveRDS(proj.list[[2]], here::here("results","m1_proj_2_R_opt2.RDS"))
proj.list[[3]] <- project_wham(fit1, proj.opts = proj.opts.list[[3]], check.version = F, do.sdrep = TRUE, TMB.jointPrecision = TRUE)
saveRDS(proj.list[[3]], here::here("results","m1_proj_3_R_opt2.RDS"))

#change R_XSPR_opt
for(i in 1:3){
	print(i)
	temp <- readRDS(here::here("results",paste0("m1_proj_",i,"_R_opt2.RDS")))
	temp$input$data$XSPR_R_opt <- temp$env$data$XSPR_R_opt <- 3 #use annual expected recruitment
	temp <- do_sdreport(temp, TMB.jointPrecision = TRUE)#, save.sdrep = TRUE, TMB.bias.correct = FALSE, TMB.jointPrecision = FALSE))
	saveRDS(temp, here::here("results", paste0("m1_proj_",i,"_R_opt3.RDS")))
	temp$input$data$XSPR_R_opt <- temp$env$data$XSPR_R_opt <- 1 #use annual random effects recruitment 
	temp <- do_sdreport(temp, TMB.jointPrecision = TRUE)#, save.sdrep = TRUE, TMB.bias.correct = FALSE, TMB.jointPrecision = FALSE))
	saveRDS(temp, here::here("results", paste0("m1_proj_",i,"_R_opt1.RDS")))
	remove(temp)
}


####################################
#change usage of recruitment for brps
# library(wham)
# fit1 <-readRDS(file.path("results","fits_no_M_re_rev.RDS"))[[2]]

fit1$input$data$XSPR_R_opt <- fit1$env$data$XSPR_R_opt <- 3 #use annual expected recruitment
fit1 <- do_sdreport(fit1, TMB.jointPrecision = TRUE)#, save.sdrep = TRUE, TMB.bias.correct = FALSE, TMB.jointPrecision = FALSE))
saveRDS(fit1, here::here("results", "m1_XSPR_R_opt_3.RDS"))

fit1$input$data$XSPR_R_opt <- fit1$env$data$XSPR_R_opt <- 1 #use annual random effects recruitment
fit1 <- do_sdreport(fit1, TMB.jointPrecision = TRUE)#, save.sdrep = TRUE, TMB.bias.correct = FALSE, TMB.jointPrecision = FALSE))
saveRDS(fit1, here::here("results", "m1_XSPR_R_opt_1.RDS"))

proj <- project_wham(temp, proj.opts = proj.opts.list[[1]], check.version = FALSE, do.sdrep = TRUE)
saveRDS(proj, here::here("results","m1_proj_1_R_opt3.RDS"))
proj <- project_wham(temp, proj.opts = proj.opts.list[[2]], check.version = F, do.sdrep = TRUE)
saveRDS(proj, here::here("results","m1_proj_2_R_opt3.RDS"))
proj <- project_wham(temp, proj.opts = proj.opts.list[[3]], check.version = F, do.sdrep = TRUE)
saveRDS(proj, here::here("results","m1_proj_3_R_opt3.RDS"))


temp <- readRDS(here::here("results", "m1_XSPR_opt_3.RDS"))
input <- temp$input
input$data$XSPR_R_opt <- 1 #use annual random effects
temp2 <- fit_wham(input, do.fit = FALSE)
temp2$opt <- temp$opt #fake it to do sdreport on projection
saveRDS(temp2, here::here("results", "m1_XSPR_opt_1.RDS"))
# temp2 <- readRDS(here::here("results", "m1_XSPR_opt_1.RDS"))
x <- TMB::sdreport(temp2)
saveRDS(x,here::here("results", "m1_XSPR_opt_1_sdrep.RDS"))
# temp2 <- readRDS(here::here("results", "m1_XSPR_opt_1.RDS"))
proj.opts.list <- readRDS(here::here("results", "m1_proj_opts_list.RDS"))
proj.list <- list()
proj <- project_wham(temp2, proj.opts = proj.opts.list[[1]], check.version = FALSE, do.sdrep = TRUE)
saveRDS(proj, here::here("results","m1_proj_1_R_opt1.RDS"))
proj <- project_wham(temp2, proj.opts = proj.opts.list[[2]], check.version = F, do.sdrep = TRUE)
saveRDS(proj, here::here("results","m1_proj_2_R_opt1.RDS"))
proj <- project_wham(temp2, proj.opts = proj.opts.list[[3]], check.version = F, do.sdrep = TRUE)
saveRDS(proj, here::here("results","m1_proj_3_R_opt1.RDS"))


for(i in 1:3){
	temp <- readRDS(here::here("results",paste0("m1_proj_",i,"_R_opt1.RDS")))
	temp$env$data$proj_F_opt[] <- temp$input$data$proj_F_opt[] <- 1
	hess <- solve(temp$sdrep$cov.fixed)
	par <- temp$opt$par
	temp$rep <- temp$report(temp$env$last.par.best)
	temp$sdrep <- sdreport(temp, par.fixed = par, hessian.fixed = hess, getJointPrecision = TRUE)
	saveRDS(temp, here::here("results",paste0("m1_proj_",i,"_last_F_R_opt1.RDS")))
}

for(i in 1:3){
	temp <- readRDS(here::here("results",paste0("m1_proj_",i,"_R_opt3.RDS")))
	temp$env$data$proj_F_opt[] <- temp$input$data$proj_F_opt[] <- 1
	hess <- solve(temp$sdrep$cov.fixed)
	par <- temp$opt$par
	temp$rep <- temp$report(temp$env$last.par.best)
	temp$sdrep <- sdreport(temp, par.fixed = par, hessian.fixed = hess, getJointPrecision = TRUE)
	saveRDS(temp, here::here("results",paste0("m1_proj_",i,"_last_F_R_opt3.RDS")))
}

x <- readRDS(here::here("results",paste0("m1_proj_",3,"_last_F_R_opt1.RDS")))
names(x$sdrep)
y <- readRDS(here::here("results",paste0("m1_proj_",3,"_last_F_R_opt3.RDS")))
names(y$sdrep)
