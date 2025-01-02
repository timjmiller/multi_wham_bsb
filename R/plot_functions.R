mypalette = function(n){
  palette.fn <- colorRampPalette(c("dodgerblue","green","red"), space = "Lab")
  palette.fn(n)
}

plot.ecov <- function(mod, plot.pad = FALSE, do.tex=FALSE, do.png=FALSE, fontfam="", res=72, od){
  origpar <- par(no.readonly = TRUE)
  dat = mod$env$data
  ecov.pred = mod$rep$Ecov_x
  ecov.pred.low <- ecov.pred.high <- ecov.pred.se <- matrix(NA, NROW(ecov.pred), NCOL(ecov.pred))
  ecov.obs = dat$Ecov_obs[1:dat$n_years_Ecov,,drop=F]
  years <- seq(from=mod$input$years_Ecov[1], by=1, length.out=NROW(ecov.obs))
  years_full <- seq(from=mod$input$years_Ecov[1], by=1, length.out=NROW(ecov.pred))#dat$n_years_Ecov+dat$n_years_proj_Ecov)

  if(class(mod$sdrep)[1] == "sdreport"){
    temp <- list(TMB:::as.list.sdreport(mod$sdrep, what = "Est", report = T),
      TMB:::as.list.sdreport(mod$sdrep, what = "Std", report = T))
    ecov.pred.se[] <- temp[[2]]$Ecov_x #TMB:::as.list.sdreport(mod$sdrep, what = "Std.", report=TRUE)$Ecov_x
    ecov.pred.low[] <- ecov.pred - 1.96 * ecov.pred.se
    ecov.pred.high[] <- ecov.pred + 1.96 * ecov.pred.se
  }

  ecov.obs.sig = mod$rep$Ecov_obs_sigma # Ecov_obs_sigma is filled with fixed, or estimated values (fe or re) for each covariate depending on the respective options
  ecov.use = dat$Ecov_use_obs[1:dat$n_years_Ecov,,drop=F]
  ecov.obs.sig = ecov.obs.sig[1:dat$n_years_Ecov,,drop=F]
  ecov.obs.sig[ecov.use == 0] <- NA

  # default: don't plot the padded entries that weren't used in ecov likelihood
  if(!plot.pad) ecov.obs[ecov.use == 0] <- NA

  ecov.res = (ecov.obs - ecov.pred[1:dat$n_years_Ecov,]) / ecov.obs.sig # standard residual (obs - pred)

  ecovs <- 1:dat$n_Ecov
  plot.colors = mypalette(dat$n_Ecov)

  ecov.low <- ecov.obs - 1.96 * ecov.obs.sig
  ecov.high <- ecov.obs + 1.96 * ecov.obs.sig
  y.min <- ifelse(min(ecov.low,na.rm=T) < 0, 1.1*min(ecov.low,na.rm=T), 0.9*min(ecov.low,na.rm=T))
  y.max <- ifelse(max(ecov.high,na.rm=T) < 0, 0.9*max(ecov.high,na.rm=T), 1.1*max(ecov.high,na.rm=T))
  if(max(ecov.pred,na.rm=T) > y.max) y.max <- max(ecov.pred,na.rm=T)
  if(min(ecov.pred,na.rm=T) < y.min) y.min <- min(ecov.pred,na.rm=T)
  plot(years_full, ecov.pred[,1], type='n', xlab="Year", ylab="Bottom Temperature Anomaly",
       ylim=c(y.min, y.max))
  for (i in ecovs)
  {
    polygon(c(years_full,rev(years_full)), c(ecov.pred.low[,i], rev(ecov.pred.high[,i])), col=adjustcolor(plot.colors[i], alpha.f=0.4), border = "transparent")
    arrows(years, ecov.low[,i], years, ecov.high[,i], length=0, col = plot.colors[i])
    points(years, ecov.obs[,i], pch=19, col = plot.colors[i])
    lines(years_full, ecov.pred[,i], col=plot.colors[i], lwd=3)
    if(length(years_full) > length(years)) abline(v=tail(years,1), lty=2)
  }
  return(cbind.data.frame(Year = years_full, est = c(ecov.obs, ecov.pred), lo = c(ecov.low, ecov.pred.low), hi = c(ecov.high, ecov.pred.high), 
    region = rep(c("North", "South"), each = length(years_full)),
    type = rep(c("obs","pred"), each = 2*length(years_full))))
}

self_test_i <- function(i, fit_RDS = NULL, seeds = NULL, conditional = TRUE, map_change = NULL, res_dir = NULL, wham_location = NULL, test_dir = NULL, save_inputs = FALSE, 
  reml = FALSE, MakeADFun.silent = TRUE) {
  
  if(is.null(test_dir)) library(wham, lib.loc = wham_location)
  else pkgload::load_all(test_dir)
  fit <- readRDS(fit_RDS)
  sim_input <- fit$input
  sim_input$par <- fit$parList
  sim_input$data$do_SPR_BRPs[] <- 0
  sim_input$random <- NULL # so fit_wham doesn't try to do inner optimization
  if(conditional){
    temp <- grep(glob2rx("do_simulate_*_re"), names(sim_input$data), value = TRUE)
    # temp <- c("do_simulate_Ecov_re", "do_simulate_L_re", "do_simulate_M_re", "do_simulate_mu_prior_re", "do_simulate_mu_re", "do_simulate_N_re",
    #           "do_simulate_q_prior_re", "do_simulate_q_re", "do_simulate_sel_re")
    sim_input$data[temp] <- rep(list(0), length(temp))
    #sim_input$data[temp] <- lapply(temp, function(x) sim_input$data[[x]][] <- 0)
  }
  if(!is.null(map_change)) sim_input$map[names(map_change)] <- map_change
  sim_mod <- try(fit_wham(sim_input, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = MakeADFun.silent))
  if(!(is.character(sim_mod))){
    set.seed(seeds[i])
    sim_input$data <- sim_mod$simulate(complete=TRUE)
    sim_input$random <- fit$input$random #set random correctly for estimation
    if(reml){
      reml_random <- c("mean_rec_pars","logit_q","F_pars","trans_mu","log_N1","logit_selpars","Mpars", "Ecov_beta_R", "Ecov_beta_M", "Ecov_beta_q", "Ecov_beta_mu")
      #first col of L_repars is mu
      #first col of Ecov_process_pars is mu/Ecov1
      sim_input$random <- unique(c(sim_input$random, reml_random))
    }
    x <- try(fit_wham(sim_input, do.sdrep = FALSE, do.retro = FALSE, do.osa = FALSE, do.brps = FALSE, MakeADFun.silent = MakeADFun.silent))
    out <- list(obj = NA, 
      par = rep(NA,length(sim_mod$par)), 
      grad = rep(NA, length(sim_mod$par)), 
      SSB = matrix(NA,NROW(sim_mod$rep$SSB),NCOL(sim_mod$rep$SSB)), 
      F = rep(NA,length(sim_mod$rep$log_F_tot)), 
      NAA = array(NA, dim = dim(sim_mod$rep$NAA)))
    out$rep <- NULL
    out$parList <- NULL
    out$opt <- NULL
    out$last.par.best <- rep(NA, length(fit$env$last.par.best))
    out$sim_input <- sim_input
    out$seed <- seeds[i]
    if(!(is.character(x))) if(!is.null(x$opt)) if(!is.character(x$opt)){
      out$opt <- x$opt
      out$rep <- x$rep
      out$parList <- x$parList
      out$last.par.best <- x$env$last.par.best
      out$obj <- x$opt$obj
      out$par <- x$opt$par
      out$grad <- x$final_gradient
      out$SSB <- x$rep$SSB
      out$F <- exp(x$rep$log_F_tot)
      out$NAA <- x$rep$NAA
    }
    if(!is.null(res_dir)){
      sim_type <- "sim"
      if(conditional) sim_type <- "cond_sim"
      if(save_inputs) saveRDS(sim_input, file.path(res_dir, paste0(sim_type, "_input_", i, ".RDS")))
      saveRDS(out, file.path(res_dir, paste0(sim_type, "_", i, ".RDS")))
    }
    return(out)
  } else return(NULL)
}
self_test <- function(fit_RDS = NULL, n = 10, seeds = NULL, which_seeds = NULL, conditional = TRUE, map_change = NULL, do_parallel = TRUE, n_cores  = NULL, 
  res_dir = NULL, wham_location = NULL, test_dir = NULL, save_inputs = FALSE, reml = FALSE){
  
  if(is.null(fit_RDS)) stop("Provide fit_RDS, an RDS file name for a fitted WHAM model.")
  if(!is.null(res_dir)) {
    cat("res_dir is provided, so self test files will be saved to ", res_dir, ". \n")
  }
  #if(is.null(wham_location)) wham_location <- system.file(package="wham")
  if(is.null(which_seeds)) which_seeds <- 1:n
  is_snowfall <- nchar(system.file(package="snowfall"))>0
  is_parallel <- nchar(system.file(package="parallel"))>0
  
  # mod <- readRDS(fit_RDS)
  # sim_type <- "sim"
  # if(conditional) sim_type <- "cond_sim"

  if(is.null(seeds)) {
    set.seed(8675309)
    seeds <- sample(0:1e9, n)
  }
  if(!all(which_seeds %in% 1:length(seeds))) stop("some of which_seeds are outside 1 to n.")
  
  if(do_parallel){
    if(is_snowfall & is_parallel){
      if(is.null(n_cores)) n_cores <- parallel::detectCores()/2
      if(!is.null(res_dir)) snowfall::sfInit(parallel=TRUE, cpus=n_cores, slaveOutfile=file.path(res_dir,"self_test_log.txt"))
      snowfall::sfExportAll()
      sim_res <- try(snowfall::sfLapply(which_seeds, function(i) self_test_i(i = i, fit_RDS = fit_RDS, seeds = seeds, conditional = conditional, map_change = map_change, 
        res_dir = res_dir, wham_location = wham_location, test_dir = test_dir, save_inputs = save_inputs, reml = reml)))
      snowfall::sfStop()
    } else stop("To do self test fits in parallel, install the snowfall and parallel packages. Otherwise, set do_parallel = FALSE.")
  } else{
    if(!(is_snowfall & is_parallel)) cat("If snowfall and parallel packages are installed, self test can be fit in parallel. \n")
    sim_res <- list()
    for(i in 1:length(which_seeds)){
      # fit <- readRDS(fit_RDS)
      # sim_input <- fit$input
      # sim_input$par <- fit$parList
      # sim_input$data$do_SPR_BRPs[] <- 0
      # sim_input$random <- NULL # so fit_wham doesn't try to do inner optimization
      # if(conditional){
      #     temp <- grep(glob2rx("do_simulate_*_re"), names(sim_input$data), value = TRUE)
      #   # temp <- c("do_simulate_Ecov_re", "do_simulate_L_re", "do_simulate_M_re", "do_simulate_mu_prior_re", "do_simulate_mu_re", "do_simulate_N_re",
      #   #           "do_simulate_q_prior_re", "do_simulate_q_re", "do_simulate_sel_re")
      #   sim_input$data[temp] <- rep(list(0), length(temp))
      #   # sim_input$data[temp] <- lapply(temp, function(x) sim_input$data[[x]][] <- 0)
      # }
      # if(!is.null(map_change)) sim_input$map[names(map_change)] <- map_change
      # sim_mod <- fit_wham(sim_input, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)
      # set.seed(seeds[i])
      # sim_input$data <- sim_mod$simulate(complete=TRUE)
      # sim_input$random <- fit$input$random #set random correctly for estimation
      # if(reml){
      #   reml_random <- c("mean_rec_pars","logit_q","F_pars","trans_mu","log_N1","logit_selpars","Mpars", "Ecov_beta_R", "Ecov_beta_M", "Ecov_beta_q", "Ecov_beta_mu")
      #   #first col of L_repars is mu
      #   #first col of Ecov_process_pars is mu/Ecov1
      #   sim_input$random <- unique(c(sim_input$random, reml_random))
      # }
      # x <- try(fit_wham(sim_input, do.sdrep = FALSE, do.retro = FALSE, do.osa = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE))
      # out <- list(obj = NA, 
      #   par = rep(NA,length(sim_mod$par)), 
      #   grad = rep(NA, length(sim_mod$par)), 
      #   SSB = matrix(NA,NROW(sim_mod$rep$SSB),NCOL(sim_mod$rep$SSB)), 
      #   F = rep(NA,length(sim_mod$rep$log_F_tot)), 
      #   NAA = array(NA, dim = dim(sim_mod$rep$NAA)))
      # out$rep <- NULL
      # out$parList <- NULL
      # out$opt <- NULL
      # out$last.par.best <- rep(NA, length(fit$env$last.par.best))
      # out$sim_input <- sim_input
      # out$seed <- seeds[i]
      # if(!(is.character(x) | is.null(x$opt))){
      #   out$opt <- x$opt
      #   out$rep <- x$rep
      #   out$parList <- x$parList
      #   out$last.par.best <- x$env$last.par.best
      #   out$obj <- x$opt$obj
      #   out$par <- x$opt$par
      #   out$grad <- x$final_gradient
      #   out$SSB <- x$rep$SSB
      #   out$F <- exp(x$rep$log_F_tot)
      #   out$NAA <- x$rep$NAA
      # }
      # out$seed <- seeds[i]
      # sim_res[[i]] <- out
      # if(!is.null(res_dir)){
      #   if(save_inputs) saveRDS(sim_input, file.path(res_dir, paste0(sim_type, "_input_", i, ".RDS")))
      #   saveRDS(out, file.path(res_dir, paste0(sim_type, "_", i, ".RDS")))
      # }
      sim_res[[i]] <- try(self_test_i(i = i, fit_RDS = fit_RDS, seeds = seeds, conditional = conditional, map_change = map_change, 
        res_dir = res_dir, wham_location = wham_location, test_dir = test_dir, save_inputs = save_inputs, reml = reml))
    }
  }
  return(list(self_test_results = sim_res, seeds = seeds))
}

summarize_res_fn <- function(sims, rep_name = "SSB", index = 1, fit){
  out <- list()
  out$resid <- sapply(sims, \(x) {
    if(!is.null(x$rep[[rep_name]])) return(x$rep[[rep_name]][,index]/fit$rep[[rep_name]][,index] - 1)
    else return(rep(NA, NROW(fit$rep[[rep_name]])))
  })
  out$n <- apply(out$resid,1, \(x) sum(!is.na(x))) #by year 
  out$median_res <- apply(out$resid, 1, median, na.rm = TRUE)
  # print(out$median_res)
  out$bnds <- cbind(qbinom(0.025, out$n, 0.5)/out$n, qbinom(0.975, out$n, 0.5)/out$n)
  out$quantiles <- sapply(1:length(out$n), \(y) quantile(out$resid[y,], probs = c(out$bnds[y,1],0.5, out$bnds[y,2]), na.rm = TRUE))
  return(out)
}
# x <- summarize_res_fn(sims, fit = fit_1_le)

plot_SSB <- function(sims, rep_name = "SSB", indices = 1:2, index_names = c("North", "South"), fit, do_plot = TRUE){
  par(mfrow = c(1,length(indices)))
  sum_res <- lapply(indices, \(x) summarize_res_fn(sims=sims, rep_name = rep_name, index = x, fit=fit))
  print(length(sum_res))
  ylim = range(sapply(indices, \(x) range(sum_res[[x]]$quantiles, na.rm = TRUE)), na.rm = TRUE)
  for(i in 1:length(indices)){
    plot(fit$years, sum_res[[i]]$quantiles[2,], ylim = ylim, ylab = paste0("Rel Diff. ", rep_name), xlab = "Year")
    polygon(c(fit$years,rev(fit$years)), c(sum_res[[i]]$quantiles[1,], rev(sum_res[[i]]$quantiles[3,])), col=adjustcolor("black", alpha.f=0.4), border = "transparent")
    mtext(side = 3, index_names[i])
    abline(h = 0)
    #median of annual medians (n_years)
    abline(h = median(sum_res[[i]]$median_res, na.rm = TRUE), lty =2)
    #median of all resids (n_sims x n_years)
    abline(h = median(sum_res[[i]]$resid, na.rm = TRUE), lty =2, col = "red")
  }
}

whichList <- function(mod, par = mod$env$last.par.best) {
  nms <- names(mod$env$parameters)
  which_by_name <- function(x) {
    which_par <- which(names(par) == x)
    y <- attr(mod$env$parameters[[x]],"shape")
    if(is.null(y)) { #all estimated
      y <- mod$env$parameters[[x]] 
      y[] <- NA
      y[] <- which_par
    } else { #mapping used, not all estimated and/or some commonly estimated
      y[] <- NA
      f <- attr(mod$env$parameters[[x]], "map")
      ind_unique <- unique(f[which(f>=0)]) #location of each unique estimated par
      if(length(which_par)) for(l in 1:length(ind_unique)) {
        y[which(f == ind_unique[l])] <- which_par[ind_unique[l]+1] #+1 because starts at 0
      }
    }
    y
  }
  out <- lapply(nms,which_by_name)
  names(out) <- nms
  return(out)
}

make_pred_R_cv <- function(proj, fjp=NULL, type = 1, region = 1){

  if(is.null(fjp)) fjp <- proj$sdrep
  if(is.null(fjp$jointPrecision)) stop("no jointPrecision matrix")
  x <- TMB:::as.list.sdreport(fjp, report = FALSE, what = "Est")
  cov_full <- solve(fjp$jointPrecision)
  index <- whichList(proj, proj$env$last.par.best)
  if(type == 1){ #continue ecov process
    #theta = mean_rec_par, ecov_mu, Ecov_beta_R, Ecov_re
    theta <- c(x$mean_rec_pars[region,1], x$Ecov_process_pars[1,region], x$Ecov_beta_R[region,region,1], x$Ecov_re[,region])
    par_index <- c(index$mean_rec_pars[region,1], index$Ecov_process_pars[1,region], index$Ecov_beta_R[region,region,1], index$Ecov_re[,region])
    print(paste("type: ", type))
    print(par_index)
    #log_pred_R = mean_rec_pars + Ecov_beta_R * ecov_mu + Ecov_beta_R * Ecov_re
    log_pred_R <- theta[1] + theta[3] * (theta[2] + x$Ecov_re[,region])
    jac <- cbind(1,theta[3], theta[2]+x$Ecov_re[,region], diag(theta[3], length(x$Ecov_re[,region])))
    # if(region == 2) ind[3] <- NA #beta not estimated
    cov <- cov_full[par_index,par_index]
    if(region == 2) cov[is.na(cov)] <- 0
    sd <- sqrt(diag(jac %*% cov %*% t(jac)))
    log_pred_R <- log_pred_R[which(proj$input$years_Ecov %in% proj$years_full)]
    sd <- sd[which(proj$input$years_Ecov %in% proj$years_full)]
  }
  if(type %in% 2:3){ # use average ecov in recent years
    #theta = mean_rec_par, ecov_mu, Ecov_beta_R, Ecov_re
    ind_m_yrs <- which(proj$input$years_Ecov %in% proj$years) #which Ecov years are model years
    # ind_m_yrs <- which(proj$input$years_Ecov %in% proj$years[proj$input$data$avg_years_Ecov+1]) #which Ecov years are averaged
    print(ind_m_yrs)
    theta <- c(x$mean_rec_pars[region,1], x$Ecov_process_pars[1,region], x$Ecov_beta_R[region,region,1], x$Ecov_re[ind_m_yrs,region])
    par_index <- c(index$mean_rec_pars[region,1], index$Ecov_process_pars[1,region], index$Ecov_beta_R[region,region,1], index$Ecov_re[ind_m_yrs,region])
    #log_pred_R = mean_rec_pars + Ecov_beta_R * ecov_mu + Ecov_beta_R * Ecov_re
    log_pred_R <- theta[1] + theta[3] * (theta[2] + x$Ecov_re[ind_m_yrs,region])
    jac <- cbind(1,theta[3], theta[2]+x$Ecov_re[ind_m_yrs,region], diag(theta[3], length(ind_m_yrs)))
    print(dim(jac))
    # if(region == 2) ind[3] <- NA #beta not estimated
    print(paste("type: ", type))
    print(par_index)
    cov <- cov_full[par_index,par_index]
    if(region == 2) cov[is.na(cov)] <- 0
    sd <- sqrt(diag(jac %*% cov %*% t(jac)))
    ind_p_yrs <- which(proj$input$years_full > max(proj$years)) #which projection years of Ecov_out_R to use
    if(type == 2){
      ind_avg_yrs <- which(proj$input$years_Ecov %in% tail(proj$years,5)) #which model years were used for average
      print(ind_p_yrs)
      theta <- c(x$mean_rec_pars[region,1], x$Ecov_process_pars[1,region], x$Ecov_beta_R[region,region,1], x$Ecov_re[ind_avg_yrs,region])
      par_index <- c(index$mean_rec_pars[region,1], index$Ecov_process_pars[1,region], index$Ecov_beta_R[region,region,1], index$Ecov_re[ind_avg_yrs,region])
      # if(region == 2) ind[3] <- NA #beta not estimated
      log_pred_R <- c(log_pred_R, theta[1] + theta[3] * (theta[2] + rep(mean(x$Ecov_re[ind_avg_yrs,region]),length(ind_p_yrs))))
      jac <- c(1,theta[3], theta[2]+mean(x$Ecov_re[ind_p_yrs,region]), rep(theta[3]/5, length(ind_avg_yrs)))
      jac <- matrix(jac, nrow = length(ind_p_yrs), ncol = length(jac), byrow = TRUE)
    }
    if(type == 3){
      ind_p_yrs <- which(proj$input$years_full > max(proj$years)) #which projection years of Ecov_out_R to use
      print(ind_p_yrs)
      theta <- c(x$mean_rec_pars[region,1], x$Ecov_beta_R[region,region,1])
      par_index <- c(index$mean_rec_pars[region,1], index$Ecov_beta_R[region,region,1])
      # if(region == 2) ind[2] <- NA #beta not estimated
      log_pred_R <- c(log_pred_R, theta[1] + theta[2] * proj$rep$Ecov_out_R[region,ind_p_yrs,region])
      jac <- cbind(1,proj$rep$Ecov_out_R[region,ind_p_yrs,region])
    }
    print(paste("type: ", type))
    print(par_index)
    cov <- cov_full[par_index,par_index]
    if(region == 2) cov[is.na(cov)] <- 0
    sd <- c(sd, sqrt(diag(jac %*% cov %*% t(jac))))
  }
  print(sd)
  return(cbind(pred_R = exp(log_pred_R), cv = sd))
}

plot.catch <- function(mod, fleet_names = mod$input$fleet_names)
{
  years <- mod$years
  dat = mod$env$data
  years_full = mod$years_full
  pred_catch = mod$rep$pred_catch
  sigma = dat$agg_catch_sigma %*% diag(exp(mod$parList$log_catch_sig_scale), nrow = length(mod$parList$log_catch_sig_scale)) # dims: [ny,nf] x [nf]
  catch = dat$agg_catch
  fleets <- 1:dat$n_fleets
  cols <- viridis::viridis_pal(option = "H")(length(fleets))
  poly_cols <- viridis::viridis_pal(option = "H", alpha = 0.3)(length(fleets))
  df <- cbind.data.frame(year = years, catch = c(dat$agg_catch), fleet = rep(fleet_names, each = NROW(dat$agg_catch)), obs_cv = c(dat$agg_catch_sigma),
    pred_catch = c(mod$rep$pred_catch), region = rep(c("North", "South"), each = NROW(dat$agg_catch)*2))
  df$lo <- df$catch*exp(qnorm(0.025)*df$obs_cv)
  df$hi <- df$catch*exp(qnorm(0.975)*df$obs_cv)

  ggplot(df, aes(x = year, y = catch)) + 
  facet_grid(fleet~region) + ylab("Catch (mt)") + xlab("Year") +
  geom_point() + 
  geom_pointrange(aes(ymin = lo, ymax = hi)) + 
  geom_line(aes(x = year, y = pred_catch))
}

plot.selectivity <- function(mod, blocks = NULL, block.names,fontfam="", od){
  dat = mod$env$data
  rep = mod$rep
  years = mod$years
  n_years = length(years)
  n_ages = dat$n_ages
  ages <- 1:n_ages
  ages.lab = 1:n_ages
  if(!is.null(mod$ages.lab)) ages.lab = mod$ages.lab

  # selAA for all blocks using facet_wrap
    n_selblocks <- length(blocks)
    sel_mod <- c("age-specific","logistic","double-logistic","decreasing-logistic")[dat$selblock_models[blocks]]
    sel_re <- c("no","IID","AR1(age)","AR1(year)","2DAR1")[dat$selblock_models_re[blocks]]
    df.selAA <- data.frame(matrix(NA, nrow=0, ncol=n_ages+2))
    colnames(df.selAA) <- c(paste0("Age_",1:n_ages),"Year","Block")
    block.names <- paste0(block.names, ": ", sel_mod,"\n(",sel_re," random effects)")
    for(i in blocks) {#if(include.selblock[i]){
      tmp = as.data.frame(rep$selAA[[i]])
      tmp$Year <- years
      colnames(tmp) <- c(paste0("Age_",1:n_ages),"Year")
      tmp$Block = block.names[which(blocks == i)]
      df.selAA <- rbind(df.selAA, tmp)
    }
    print(blocks)
    print(dim(df.selAA))
    df.plot <- df.selAA %>% tidyr::pivot_longer(-c(Year,Block),
              names_to = "Age", 
              names_prefix = "Age_",
              names_ptypes = list(Age = character()),
              values_to = "Selectivity")
    df.plot$Age <- as.factor(as.integer(df.plot$Age))
    levels(df.plot$Age) = ages.lab
    print(dim(df.plot))
    print(block.names)
    print(unique(as.character(df.plot$Block)))
    df.plot$Block <- factor(as.character(df.plot$Block), levels=block.names)
    fn <- "SelAA_tile"
      print(ggplot2::ggplot(df.plot, ggplot2::aes(x=Year, y=Age, fill=Selectivity)) + 
        ggplot2::geom_tile() +
        ggplot2::scale_x_continuous(expand=c(0,0)) +
        ggplot2::scale_y_discrete(expand=c(0,0)) + #, breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +        
#        ggplot2::scale_y_continuous(expand=c(0,0), breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +        
        ggplot2::theme_bw() + 
        ggplot2::facet_wrap(~Block, dir="v") +
        ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(size = ggplot2::rel(1.5)), 
          axis.title = ggplot2::element_text(size = ggplot2::rel(2)), axis.text = ggplot2::element_text(size = ggplot2::rel(2)),
          legend.title = ggplot2::element_text(size = ggplot2::rel(1.5), hjust=0.5), legend.box.just = "center",
          legend.text = ggplot2::element_text(size = ggplot2::rel(1.5))) +#, strip.text.y = element_blank()) +
        ggplot2::scale_fill_viridis_c(begin = 0.2, end = 0.8, option = "turbo"))

}

get.brp.status.results <- function(mod, static = FALSE, alpha = 0.05){
  n_stocks <- mod$env$data$n_stocks
  n_fleets <- mod$env$data$n_fleets
  n_regions <- mod$env$data$n_regions
  percentSPR = mod$env$data$percentSPR
  n_yrs <- length(mod$years_full)
  std <- summary(mod$sdrep, "report")
  inds <- list()
  inds$ssb <- matrix(c(which(rownames(std) == "log_SSB"), which(rownames(std) == "log_SSB_all")), ncol = n_stocks+1)
  inds$full.f <- matrix(which(rownames(std) == "log_Fbar"), ncol = n_fleets + n_regions + 1)
  inds$F.t <- matrix(which(rownames(std) == "log_Fbar_XSPR"), ncol = n_fleets + n_regions + 1)
  inds$SSB.t <- matrix(which(rownames(std) == "log_SSB_FXSPR"), ncol = n_stocks+1)
  if(static){
    inds$F.t <- matrix(which(rownames(std) == "log_Fbar_XSPR_static"), nrow = length(mod$years_full), ncol = n_fleets + n_regions + 1, byrow = TRUE) #only 1 value
    inds$SSB.t <- matrix(which(rownames(std) == "log_SSB_FXSPR_static"), nrow = length(mod$years_full), ncol = n_stocks+1, byrow=TRUE)
  }
  brps.f.est <- matrix(std[c(inds$F.t),1], nrow = n_yrs, ncol = NCOL(inds$full.f))
  brps.f.cv <- matrix(std[c(inds$F.t),2], nrow = n_yrs, ncol = NCOL(inds$full.f))
  brps.ssb.est <- matrix(std[c(inds$SSB.t),1], nrow = n_yrs, ncol = NCOL(inds$ssb))
  brps.ssb.cv <- matrix(std[c(inds$SSB.t),2], nrow = n_yrs, ncol = NCOL(inds$ssb))
  log.rel.f.vals <- matrix(std[c(inds$full.f),1] - std[c(inds$F.t),1], nrow = n_yrs, ncol = NCOL(inds$full.f))
  log.rel.ssb.vals <- matrix(std[c(inds$ssb),1] - std[c(inds$SSB.t),1], nrow = n_yrs, ncol = NCOL(inds$ssb))
  cov <- mod$sdrep$cov
  log.rel.ssb.rel.F.cov <- lapply(1:n_yrs, function(x) {
    K <- cbind(c(1,-1,0,0),c(0,0,1,-1))
    f.index <- n_fleets + 1:(n_regions+1) #THIS IS NOT GENERAL, JUST FOR THIS BSB APPLICATION
    if(length(f.index) != NCOL(inds$SSB.t)) stop("in kobe.plot function: number of stocks isn't equal to the number of regions")
    status.cov <- lapply(1:length(f.index), \(i) {
      ind <- c(inds$ssb[x,i],inds$SSB.t[x,i],inds$full.f[x,f.index[i]],inds$F.t[x,f.index[i]])
      tcov <- cov[ind,ind]
      return(t(K) %*% tcov %*% K)
      })
    return(status.cov)
  })
  if(mod$env$data$n_years_proj>0) { #check whether projecting at F40/FMSY because the ratio to status those years will be 1 and variance 0, but numerical accuracy might be an issue.
    proj_F40 <- which(mod$env$data$proj_F_opt==3)
    if(length(proj_F40)){
      proj_F40 <- mod$env$data$n_years_model + proj_F40
      proj_F40 <- which((1:n_yrs) %in% proj_F40)
    }
    if(length(proj_F40)){
      log.rel.ssb.rel.F.cov[proj_F40] <- lapply(log.rel.ssb.rel.F.cov[proj_F40], function(x) {
        out <- lapply(1:(n_regions+1), \(i) {
          x[[i]][cbind(c(1,2,2),c(2,2,1))] <- 0
          return(x[[i]])
        })
        return(out)
      })
    }
  }
  region_names <- c(mod$input$region_names, "Total")
  for (i in 1:(n_regions+1)){
    if(i == 1) {
      df.brps <- cbind.data.frame(Year = mod$years_full, region = region_names[i], type = "SSB", val = brps.ssb.est[,i], se = brps.ssb.cv[,i]) 
      df.brps <- rbind(df.brps, cbind.data.frame(Year = mod$years_full, region = region_names[i], type = "Fbar", val = brps.f.est[,n_fleets + i], se = brps.f.cv[,n_fleets + i]))
      df.status <- cbind.data.frame(Year = mod$years_full, region = region_names[i], type = "SSB", val = log.rel.ssb.vals[,i], se = sqrt(sapply(log.rel.ssb.rel.F.cov, function(x) x[[i]][1,1])))
      df.status <- rbind(df.status, cbind.data.frame(Year = mod$years_full, region = region_names[i], type = "Fbar", val = log.rel.f.vals[,i], se = sqrt(sapply(log.rel.ssb.rel.F.cov, function(x) x[[i]][2,2]))))
    } else{
      df.brps <- rbind(df.brps, cbind.data.frame(Year = mod$years_full, region = region_names[i], type = "SSB", val = brps.ssb.est[,i], se = brps.ssb.cv[,i]))
      df.brps <- rbind(df.brps, cbind.data.frame(Year = mod$years_full, region = region_names[i], type = "Fbar", val = brps.f.est[,n_fleets + i], se = brps.f.cv[,n_fleets + i]))
      df.status <- rbind(df.status, cbind.data.frame(Year = mod$years_full, region = region_names[i], type = "SSB", val = log.rel.ssb.vals[,i], se = sqrt(sapply(log.rel.ssb.rel.F.cov, function(x) x[[i]][1,1]))))
      df.status <- rbind(df.status, cbind.data.frame(Year = mod$years_full, region = region_names[i], type = "Fbar", val = log.rel.f.vals[,i], se = sqrt(sapply(log.rel.ssb.rel.F.cov, function(x) x[[i]][2,2]))))
    }
  }
  for(x in 1:n_yrs) for (i in 1:(n_regions+1)){
    if(is.na(log.rel.f.vals[x,n_fleets +i]) | any(diag(log.rel.ssb.rel.F.cov[[x]][[i]])<0)) temp <- matrix(NA,100,2)
    else temp <- ellipse::ellipse(log.rel.ssb.rel.F.cov[[x]][[i]], centre = c(log.rel.ssb.vals[x,i],log.rel.f.vals[x,n_fleets +i]), level = 1-alpha)
    if(x == 1 & i == 1) {
      df.ellipse <- cbind.data.frame(Year = mod$years_full[x], region = region_names[i], ptype = "center", x = log.rel.ssb.vals[x,i], y = log.rel.f.vals[x,n_fleets +i]) 
    } else {
      df.ellipse <- rbind(df.ellipse, cbind.data.frame(Year = mod$years_full[x], region = region_names[i], ptype = "center", x = log.rel.ssb.vals[x,i], y = log.rel.f.vals[x,n_fleets +i]))
    }
    df.ellipse <- rbind(df.ellipse, cbind.data.frame(Year = mod$years_full[x],  region = region_names[i], ptype = "ellipse", x = temp[,"x"], y = temp[,"y"]))
  }
  return(list(df.brps = df.brps, df.status = df.status, df.ellipse = df.ellipse))
}

kobe.plot <- function(mod, status.years=NULL, static = FALSE, regions = 1:2,single.plot = FALSE, alpha = 0.05, max.x=NULL, max.y=NULL, stock = NULL){
  if(is.null(status.years)) {
    status.years <- length(mod$years)
    if(length(mod$years_full)> status.years) status.years <- c(status.years, length(mod$years_full))
  }
  n_stocks <- mod$env$data$n_stocks
  n_fleets <- mod$env$data$n_fleets
  n_regions <- mod$env$data$n_regions
  percentSPR = mod$env$data$percentSPR
  std <- summary(mod$sdrep, "report")
  inds <- list()
  inds$ssb <- matrix(c(which(rownames(std) == "log_SSB"), which(rownames(std) == "log_SSB_all")), ncol = n_stocks+1)
  inds$full.f <- matrix(which(rownames(std) == "log_Fbar"), ncol = n_fleets + n_regions + 1)
  inds$F.t <- matrix(which(rownames(std) == "log_Fbar_XSPR"), ncol = n_fleets + n_regions + 1)
  inds$SSB.t <- matrix(which(rownames(std) == "log_SSB_FXSPR"), ncol = n_stocks+1)
  if(static){
    inds$F.t <- matrix(which(rownames(std) == "log_Fbar_XSPR_static"), nrow = length(mod$years_full), ncol = n_fleets + n_regions + 1, byrow = TRUE) #only 1 value
    inds$SSB.t <- matrix(which(rownames(std) == "log_SSB_FXSPR_static"), nrow = length(mod$years_full), ncol = n_stocks+1, byrow=TRUE)
  }
  log.rel.f.vals <- matrix(std[inds$full.f[status.years,],1] - std[inds$F.t[status.years,],1], nrow = length(status.years), ncol = NCOL(inds$full.f))
  print(status.years)
  log.rel.ssb.vals <- matrix(std[inds$ssb[status.years,],1] - std[inds$SSB.t[status.years,],1], nrow = length(status.years), ncol = NCOL(inds$ssb))
  cov <- mod$sdrep$cov
  log.rel.ssb.rel.F.cov <- lapply(status.years, function(x) {
    K <- cbind(c(1,-1,0,0),c(0,0,1,-1))
    f.index <- n_fleets + 1:(n_regions+1) #THIS IS NOT GENERAL, JUST FOR THIS BSB APPLICATION
    if(length(f.index) != NCOL(inds$SSB.t)) stop("in kobe.plot function: number of stocks isn't equal to the number of regions")
    status.cov <- lapply(1:length(f.index), \(i) {
      ind <- c(inds$ssb[x,i],inds$SSB.t[x,i],inds$full.f[x,f.index[i]],inds$F.t[x,f.index[i]])
      tcov <- cov[ind,ind]
      return(t(K) %*% tcov %*% K)
      })
    return(status.cov)
  })
  if(mod$env$data$n_years_proj>0) { #check whether projecting at F40/FMSY because the ratio to status those years will be 1 and variance 0, but numerical accuracy might be an issue.
    proj_F40 <- which(mod$env$data$proj_F_opt==3)
    if(length(proj_F40)){
      proj_F40 <- mod$env$data$n_years_model + proj_F40
      proj_F40 <- which(status.years %in% proj_F40)
    }
    if(length(proj_F40)){
      log.rel.ssb.rel.F.cov[proj_F40] <- lapply(log.rel.ssb.rel.F.cov[proj_F40], function(x) {
        out <- lapply(1:(n_regions+1), \(i) {
          x[[i]][cbind(c(1,2,2),c(2,2,1))] <- 0
          return(x[[i]])
        })
        return(out)
      })
    }
  }
  region_names <- c(mod$input$region_names, "Total")
  for (i in 1:(n_regions+1)){
    if(i == 1) {
      df.status <- cbind.data.frame(Year = mod$years_full, region = region_names[i], type = "SSB", val = log.rel.ssb.vals[,i], cv = sapply(log.rel.ssb.rel.F.cov, function(x) x[[i]][1,1])) 
      df.status <- rbind(df.status, cbind.data.frame(Year = mod$years_full, region = region_names[i], type = "Fbar", val = log.rel.f.vals[,i], cv = sapply(log.rel.ssb.rel.F.cov, function(x) x[[i]][2,2])))
    } else{
      df.status <- rbind(df.status, cbind.data.frame(Year = mod$years_full, region = region_names[i], type = "SSB", val = log.rel.ssb.vals[,i], cv = sapply(log.rel.ssb.rel.F.cov, function(x) x[[i]][1,1])))
      df.status <- rbind(df.status, cbind.data.frame(Year = mod$years_full, region = region_names[i], type = "Fbar", val = log.rel.f.vals[,i], cv = sapply(log.rel.ssb.rel.F.cov, function(x) x[[i]][2,2])))
    }
  }
  for(x in 1:length(status.years)) for (i in 1:(n_regions+1)){
    if(is.na(log.rel.f.vals[x,n_fleets +i]) | any(diag(log.rel.ssb.rel.F.cov[[x]][[i]])<0)) temp <- matrix(NA,100,2)
    else temp <- ellipse::ellipse(log.rel.ssb.rel.F.cov[[x]][[i]], centre = c(log.rel.ssb.vals[x,i],log.rel.f.vals[x,n_fleets +i]), level = 1-alpha)
    if(x == 1 & i == 1) {
      df.ellipse <- cbind.data.frame(Year = status.years[x], region = region_names[i], ptype = "center", x = log.rel.ssb.vals[x,i], y = log.rel.f.vals[x,n_fleets +i]) 
    } else {
      df.ellipse <- rbind(df.ellipse, cbind.data.frame(Year = status.years[x], region = region_names[i], ptype = "center", x = log.rel.ssb.vals[x,i], y = log.rel.f.vals[x,n_fleets +i]))
    }
    df.ellipse <- rbind(df.ellipse, cbind.data.frame(Year = status.years[x],  region = region_names[i], ptype = "ellipse", x = temp[,"x"], y = temp[,"y"]))
  }

  print(paste0("log.rel.ssb.rel.F.cov: ", log.rel.ssb.rel.F.cov, collapse = ", "))
  do.kobe <- which(sapply(log.rel.ssb.rel.F.cov, \(x) !all(sapply(x, \(y) !is.finite(y))))) # only if some non-infinite values for at least some status years 
  if(length(do.kobe)<length(status.years)){
    no.kobe <- which(!status.years %in% do.kobe)
    print(paste0("status confidence region not available for years: ", mod$years_full[status.years[no.kobe]]))
    if(length(do.kobe)==0) return()
  }
  print(paste0("do.kobe: ", do.kobe, collapse = ", "))
  if(is.null(mod$input$region_names)) mod$input$region_names <- paste0("Region ", 1:n_regions)
  rel.ssb.rel.F.cr <- lapply(1:length(status.years), function(x){ 
    out <- lapply(1:(n_regions+1), \(i) {
      if(is.na(log.rel.f.vals[x,n_fleets +i]) | any(diag(log.rel.ssb.rel.F.cov[[x]][[i]])<0)) return(matrix(NA,100,2))
      else return(exp(ellipse::ellipse(log.rel.ssb.rel.F.cov[[x]][[i]], centre = c(log.rel.ssb.vals[x,i],log.rel.f.vals[x,n_fleets +i]), level = 1-alpha)))
    })
  })
  x <- lapply(1:3, \(x) list())
  p.status <- list(p.ssb.lo.f.lo = x, p.ssb.lo.f.hi = x, p.ssb.hi.f.lo = x, p.ssb.hi.f.hi = x)
  for(i in 1:length(status.years)){
    for(j in 1:(n_regions + 1)){
      check.bad.sd <- diag(log.rel.ssb.rel.F.cov[[i]][[j]])==0 | diag(log.rel.ssb.rel.F.cov[[i]][[j]]) < 0
      if(!any(is.na(check.bad.sd))) if(!any(check.bad.sd)){
        p.status$p.ssb.lo.f.lo[[i]][[j]] <- mnormt::sadmvn(lower = c(-Inf,-Inf), upper = c(-log(2), 0), mean = c(log.rel.ssb.vals[i,j],log.rel.f.vals[i,n_fleets +j]), 
          varcov = log.rel.ssb.rel.F.cov[[i]][[j]])
        p.status$p.ssb.lo.f.hi[[i]][[j]] <- mnormt::sadmvn(lower = c(-Inf,0), upper = c(-log(2), Inf), mean = c(log.rel.ssb.vals[i,j],log.rel.f.vals[i,n_fleets +j]), 
          varcov = log.rel.ssb.rel.F.cov[[i]][[j]])
        p.status$p.ssb.hi.f.lo[[i]][[j]] <- mnormt::sadmvn(lower = c(-log(2),-Inf), upper = c(Inf, 0), mean = c(log.rel.ssb.vals[i,j],log.rel.f.vals[i,n_fleets +j]), 
          varcov = log.rel.ssb.rel.F.cov[[i]][[j]])
        p.status$p.ssb.hi.f.hi[[i]][[j]] <- mnormt::sadmvn(lower = c(-log(2),0), upper = c(Inf, Inf), mean = c(log.rel.ssb.vals[i,j],log.rel.f.vals[i,n_fleets +j]), 
          varcov = log.rel.ssb.rel.F.cov[[i]][[j]])
      }
    }
  }
  if(is.null(max.x)) max.x <- max(sapply(rel.ssb.rel.F.cr, \(x) sapply(x, \(y) max(y[,1],na.rm = TRUE))),1.5)
  if(is.null(max.y)) max.y <- max(sapply(rel.ssb.rel.F.cr, \(x) sapply(x, \(y) max(y[,2],na.rm = TRUE))),1.5)

  # lims <- par("usr")
  lims <- c(-1e5,1e5,-1e5,1e5)
  print(lims)
  tcol <- col2rgb('red')
  tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
  pal <- c(topleft = tcol)
  borders <- cbind.data.frame(region = "topleft", x = c(lims[1],0.5,0.5,lims[1]), y = c(1,1,lims[4],lims[4]))
  tcol <- col2rgb('green')
  tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
  pal <- c(pal, bottomright = tcol)
  borders <- rbind(borders, cbind.data.frame(region = "bottomright", x = c(0.5,lims[2],lims[2],0.5), y = c(lims[3],lims[3],1,1)))
  tcol <- col2rgb('yellow')
  tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
  pal <- c(pal, topright = tcol, bottomleft = tcol)
  borders <- rbind(borders, cbind.data.frame(region = "bottomleft", x = c(lims[1],0.5,0.5,lims[1]), y = c(lims[3],lims[3],1,1)))
  borders <- rbind(borders, cbind.data.frame(region = "topright", x = c(0.5,lims[2],lims[2],0.5), y = c(1,1,lims[4],lims[4])))
  p.names <- c(topleft = "p.ssb.lo.f.hi", topright = "p.ssb.hi.f.hi", bottomleft = "p.ssb.lo.f.lo", bottomright = "p.ssb.hi.f.lo")
  print(borders)
  print(p.names)
  print(c(max.x, max.y))
  print(pal)
  if(length(do.kobe)){
    pcols <- n_regions + 1
    if(single.plot){
      prows <- 1
    } else{
      prows <- length(do.kobe)
    }
    par(mfrow = c(prows,pcols), mar = c(1,1,0,0), oma = c(5,5,3,1))
    n.plots <- 1
    if(!single.plot) n.plots <- length(do.kobe)
    for(k in 1:prows) for(j in 1:pcols) {
      st.yrs.plt <- status.years[do.kobe]
      print(st.yrs.plt)
      plot(log.rel.ssb.vals[do.kobe[k],j],log.rel.f.vals[do.kobe[k],n_fleets +j], ylim = c(0,max.y), xlim = c(0,max.x), xlab = "", ylab = "", type = 'n', axes = F)
      box(lwd = 2)
      if(k > (prows-1)*pcols) axis(1, lwd = 2)
      else axis(1, lwd = 2, labels = FALSE)
      if((j + pcols*(k-1)) %in% seq(1,prows*pcols, pcols)) axis(2,lwd = 2)
      else axis(2, lwd = 2, labels = FALSE)
      for(i in 1:4){ #quadrants of kobe plot
        df <- subset(borders, region == names(pal)[i])
        print(paste0("i: ", i))
        print(df)
        polygon(df$x,df$y, border = pal[i], col = pal[i])
        legend(names(pal)[i], legend = paste0("Year: ", mod$years_full[st.yrs.plt], ", Prob = ", 
          round(p.status[[p.names[names(pal)[i]]]][[do.kobe[k]]][[j]],2)), bty = "n", text.font=1)
      }
      # stop()
      mtext(side = 3, ifelse(j <= n_regions, mod$input$region_names[j], "Total"), line = 1)
      print(log.rel.ssb.vals)
      print(log.rel.f.vals)
      print(do.kobe[k])
      text(exp(log.rel.ssb.vals[do.kobe[k],j]),exp(log.rel.f.vals[do.kobe[k],n_fleets +j]), substr(mod$years_full[st.yrs.plt],3,4), font=1)
      if(prows == 1) for(i in do.kobe) {
        polygon(rel.ssb.rel.F.cr[[k]][[j]][,1], rel.ssb.rel.F.cr[[k]][[j]][,2], lwd=1)#, border = gray(0.7))
      } else {
        polygon(rel.ssb.rel.F.cr[[k]][[j]][,1], rel.ssb.rel.F.cr[[k]][[j]][,2], lwd=1)#, border = gray(0.7))
      }
    }
    mtext(side = 1, outer = TRUE, bquote(SSB*"/"*SSB[.(percentSPR)*"%"]), line = 2)
    mtext(side = 2, outer = TRUE, bquote(paste(italic(F),"/", italic(F)[paste(.(percentSPR),"%")])), line = 2)
  }

  return(list(p.status=p.status, rel.ssb.rel.F.cr=rel.ssb.rel.F.cr, log.rel.ssb.vals=log.rel.ssb.vals, log.rel.f.vals= log.rel.f.vals, 
    df.ellipse = df.ellipse, borders = borders, poly.colors=pal))
}

make.status.plots <- function(mod) {
  lims <- par("usr")
  print(lims)
  tcol <- col2rgb('red')
  tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
  pal <- c(topleft = tcol)
  borders <- cbind.data.frame(region = "topleft", x = c(lims[1],0.5,0.5,lims[1]), y = c(1,1,lims[4],lims[4]))
  tcol <- col2rgb('green')
  tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
  pal <- c(pal, bottomright = tcol)
  borders <- rbind(borders, cbind.data.frame(region = "bottomright", x = c(0.5,lims[2],lims[2],0.5), y = c(lims[3],lims[3],1,1)))
  tcol <- col2rgb('yellow')
  tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
  pal <- c(pal, topright = tcol, bottomleft = tcol)
  borders <- rbind(borders, cbind.data.frame(region = "bottomleft", x = c(lims[1],0.5,0.5,lims[1]), y = c(lims[3],lims[3],1,1)))
  borders <- rbind(borders, cbind.data.frame(region = "topright", x = c(0.5,lims[2],lims[2],0.5), y = c(1,1,lims[4],lims[4])))
  p.names <- c(topleft = "p.ssb.lo.f.hi", topright = "p.ssb.hi.f.hi", bottomleft = "p.ssb.lo.f.lo", bottomright = "p.ssb.hi.f.lo")
  print(borders)
  print(p.names)
  print(c(max.x, max.y))
  print(pal)
  if(length(do.kobe)){
    pcols <- n_regions + 1
    if(single.plot){
      prows <- 1
    } else{
      prows <- length(do.kobe)
    }
    par(mfrow = c(prows,pcols), mar = c(1,1,0,0), oma = c(5,5,3,1))
    n.plots <- 1
    if(!single.plot) n.plots <- length(do.kobe)
    for(k in 1:prows) for(j in 1:pcols) {
      st.yrs.plt <- status.years[do.kobe]
      print(st.yrs.plt)
      plot(log.rel.ssb.vals[do.kobe[k],j],log.rel.f.vals[do.kobe[k],n_fleets +j], ylim = c(0,max.y), xlim = c(0,max.x), xlab = "", ylab = "", type = 'n', axes = F)
      box(lwd = 2)
      if(k > (prows-1)*pcols) axis(1, lwd = 2)
      else axis(1, lwd = 2, labels = FALSE)
      if((j + pcols*(k-1)) %in% seq(1,prows*pcols, pcols)) axis(2,lwd = 2)
      else axis(2, lwd = 2, labels = FALSE)
      for(i in 1:4){ #quadrants of kobe plot
        df <- subset(borders, region == names(pal)[i])
        print(paste0("i: ", i))
        print(df)
        polygon(df$x,df$y, border = pal[i], col = pal[i])
        legend(names(pal)[i], legend = paste0("Year: ", mod$years_full[st.yrs.plt], ", Prob = ", 
          round(p.status[[p.names[names(pal)[i]]]][[do.kobe[k]]][[j]],2)), bty = "n", text.font=1)
      }
      mtext(side = 3, ifelse(j <= n_regions, mod$input$region_names[j], "Total"), line = 1)
      print(log.rel.ssb.vals)
      print(log.rel.f.vals)
      print(do.kobe[k])
      text(exp(log.rel.ssb.vals[do.kobe[k],j]),exp(log.rel.f.vals[do.kobe[k],n_fleets +j]), substr(mod$years_full[st.yrs.plt],3,4), font=1)
      if(prows == 1) for(i in do.kobe) {
        polygon(rel.ssb.rel.F.cr[[k]][[j]][,1], rel.ssb.rel.F.cr[[k]][[j]][,2], lwd=1)#, border = gray(0.7))
      } else {
        polygon(rel.ssb.rel.F.cr[[k]][[j]][,1], rel.ssb.rel.F.cr[[k]][[j]][,2], lwd=1)#, border = gray(0.7))
      }
    }
    mtext(side = 1, outer = TRUE, bquote(SSB*"/"*SSB[.(percentSPR)*"%"]), line = 2)
    mtext(side = 2, outer = TRUE, bquote(paste(italic(F),"/", italic(F)[paste(.(percentSPR),"%")])), line = 2)
  }

  return(list(p.status=p.status, rel.ssb.rel.F.cr=rel.ssb.rel.F.cr, log.rel.ssb.vals=log.rel.ssb.vals, log.rel.f.vals= log.rel.f.vals, 
    df.ellipse = df.ellipse, borders = borders, poly.colors=pal))
}
