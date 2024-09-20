fits_M_re <- readRDS(here("results", paste0("fits_M_re_better.RDS")))

ny <- fits_M_re[[1]]$input$data$n_years_model
plot(x[,1],x[,2])
plot(x[,1],x[,2]/x[,1])

x <- cbind(SSB = fits_M_re[[1]]$rep$SSB[-ny,2],R = fits_M_re[[1]]$rep$NAA[2,2,-1,1])
plot(x[,1],x[,2])
plot(x[,1],x[,2]/x[,1])
plot(x[,1],log(x[,2]/x[,1]))

x <- cbind(SSB = fits_M_re[[1]]$rep$SSB[-ny,1],R = fits_M_re[[1]]$rep$NAA[1,1,-1,1])
fn <- function(par){
	mu <- par[1] + log(x[,1]) - log(1 + exp(par[2])*x[,1])
	return(-sum(dnorm(log(x[,2]), mu, sd = exp(par[3]), log = TRUE)))
}
phi0 <- mean(exp(fits_M_re[[1]]$rep$log_SPR0[,1]))
opt <- nlminb(c(0,-5,0), fn)

h <- exp(opt$par[1])*phi0/(4 + exp(opt$par[1])*phi0)

x <- cbind(SSB = fits_M_re[[1]]$rep$SSB[-ny,2],R = fits_M_re[[1]]$rep$NAA[2,2,-1,1])
phi0 <- c(phi0, mean(exp(fits_M_re[[1]]$rep$log_SPR0[,2])))

fn <- function(par){
	mu <- par[1] + log(x[,1]) - log(1 + exp(-6)*x[,1])
	return(-sum(dnorm(log(x[,2]), mu, sd = exp(par[2]), log = TRUE)))
}
opt <- nlminb(c(0,0), fn)

fn <- function(par){
	h <- 0.2 + 0.8/(1 + exp(-par[1]))
	R0 <- exp(par[2])
	mu <- log(4 * R0 * h * x[,1]/((1-h)*R0*phi0[2] + (5*h - 1)* x[,1]))
	return(-sum(dnorm(log(x[,2]), mu, sd = exp(par[3]), log = TRUE)))
}
opt <- nlminb(c(0,0,0), fn)

h <- c(h,exp(opt$par[1])*phi0/(4 + exp(opt$par[1])*phi0))


#############################################
#try BH
nms <- names(input_0$par)
nms <- nms[which(nms != "mean_rec_pars")]

NAA_re_1 <- NAA_re
NAA_re_1$recruit_model <- 3

input_bh_0 <- prepare_wham_input(asap, selectivity = sel, NAA_re = NAA_re_1, basic_info = basic_info, move = move, ecov = ecov_0, catch_info = catch_info,
	index_info = index_info, age_comp = age_comp)
input_bh_0$par[nms] <- fit_0$parList[nms]
input_bh_0$par$mean_rec_pars[,2] <- -7
fit_bh_0 <- fit_wham(input_bh_0, do.brps = FALSE, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE)
#doesn't work for either north or south.
#############################################
