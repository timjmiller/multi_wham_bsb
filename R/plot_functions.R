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

make_pred_R_cv <- function(proj, fjp, type = 1){

  x <- TMB:::as.list.sdreport(proj$sdrep, report = FALSE, what = "Est")
  cov_full <- solve(fjp$jointPrecision)
  if(type == 1){ #continue ecov process
    #theta = mean_rec_par, ecov_mu, Ecov_beta_R, Ecov_re
    theta <- c(x$mean_rec_pars[1,1], x$Ecov_process_pars[1,1], x$Ecov_beta_R[1,1,1], x$Ecov_re[,1])
    #log_pred_R = mean_rec_pars + Ecov_beta_R * ecov_mu + Ecov_beta_R * Ecov_re
    log_pred_R <- theta[1] + theta[3] * (theta[2] + x$Ecov_re[,1])
    jac <- cbind(1,theta[3], theta[2]+x$Ecov_re[,1], diag(theta[3], length(x$Ecov_re[,1])))
    ind <- match(theta, proj$env$last.par.best) #sketchy but seems to work
    print(paste("type: ", type))
    print(ind)
    cov <- cov_full[ind,ind]
    sd <- sqrt(diag(jac %*% cov %*% t(jac)))
    log_pred_R <- log_pred_R[which(proj$input$years_Ecov %in% proj$years_full)]
    sd <- sd[which(proj$input$years_Ecov %in% proj$years_full)]
  }
  if(type %in% 2:3){ # use average ecov in recent years
    #theta = mean_rec_par, ecov_mu, Ecov_beta_R, Ecov_re
    ind_m_yrs <- which(proj$input$years_Ecov %in% proj$years) #which Ecov years are
    print(ind_m_yrs)
    theta <- c(x$mean_rec_pars[1,1], x$Ecov_process_pars[1,1], x$Ecov_beta_R[1,1,1], x$Ecov_re[ind_m_yrs,1])
    ind <- match(theta, proj$env$last.par.best) #sketchy but seems to work
    print(paste("type: ", type))
    print(ind)
    #log_pred_R = mean_rec_pars + Ecov_beta_R * ecov_mu + Ecov_beta_R * Ecov_re
    log_pred_R <- theta[1] + theta[3] * (theta[2] + x$Ecov_re[ind_m_yrs,1])
    jac <- cbind(1,theta[3], theta[2]+x$Ecov_re[ind_m_yrs,1], diag(theta[3], length(ind_m_yrs)))
    print(dim(jac))
    cov <- cov_full[ind,ind]
    sd <- sqrt(diag(jac %*% cov %*% t(jac)))
    ind_p_yrs <- which(proj$input$years_full > max(proj$years)) #which projection years of Ecov_out_R to use
    if(type == 2){
      ind_avg_yrs <- which(proj$input$years_Ecov %in% tail(proj$years,5)) #which model years were used for average
      print(ind_p_yrs)
      theta <- c(x$mean_rec_pars[1,1], x$Ecov_process_pars[1,1], x$Ecov_beta_R[1,1,1], x$Ecov_re[ind_avg_yrs,1])
      ind <- match(theta, proj$env$last.par.best) #sketchy but seems to work
      print(ind)
      log_pred_R <- c(log_pred_R, theta[1] + theta[3] * (theta[2] + rep(mean(x$Ecov_re[ind_avg_yrs,1]),length(ind_p_yrs))))
      jac <- c(1,theta[3], theta[2]+mean(x$Ecov_re[ind_p_yrs,1]), rep(theta[3]/5, length(ind_avg_yrs)))
      jac <- matrix(jac, nrow = length(ind_p_yrs), ncol = length(jac), byrow = TRUE)
    }
    if(type == 3){
      ind_p_yrs <- which(proj$input$years_full > max(proj$years)) #which projection years of Ecov_out_R to use
      print(ind_p_yrs)
      theta <- c(x$mean_rec_pars[1,1], x$Ecov_beta_R[1,1,1])
      ind <- match(theta, proj$env$last.par.best) #sketchy but seems to work
      print(ind)
      log_pred_R <- c(log_pred_R, theta[1] + theta[2] * proj$rep$Ecov_out_R[1,ind_p_yrs,1])
      jac <- cbind(1,proj$rep$Ecov_out_R[1,ind_p_yrs,1])
    }
    print(dim(jac))
    print(length(ind))
    cov <- cov_full[ind,ind]
    sd <- c(sd, sqrt(diag(jac %*% cov %*% t(jac))))
    print(length(log_pred_R))
    print(length(sd))
  }
  return(cbind(pred_R = exp(log_pred_R), cv = sd))
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
    # block.names <- paste0("Block ",1:n_selblocks,": ", sel_mod,"\n(",sel_re," random effects)")
    # block.fleets.indices <- lapply(blocks, function(x){
    #   y <- dat$selblock_pointer_fleets
    #   z <- matrix(as.integer(y == x), NROW(y), NCOL(y))
    #   fleet_ind <- apply(z,2,any)
    #   out <- mod$input$fleet_names[which(fleet_ind)]
    #   y <- dat$selblock_pointer_indices
    #   z <- matrix(as.integer(y == x), NROW(y), NCOL(y))
    #   index_ind <- apply(z,2,any)
    #   out <- c(out, mod$input$index_names[which(index_ind)])
    # })
    # include.selblock <- sapply(block.fleets.indices, length) > 0
    # for(i in 1:n_selblocks) if(include.selblock[i]){
    #   block.names[i] <- paste0(block.names[i], "\n", paste(block.fleets.indices[[i]], collapse = ", "))
    # }
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
    # if(do.tex) cairo_pdf(file.path(od, paste0(fn, ".pdf")), family = fontfam, height = 10, width = 10)
    # if(do.png) png(filename = file.path(od, paste0(fn, ".png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
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
       viridis::scale_fill_viridis())

}