asap <- read_asap3_dat(here("data",c("north.dat","south.dat")))
temp <- prepare_wham_input(asap)

#basic_info
basic_info <- list(region_names = c("North", "South"), stock_names = paste0("BSB_", c("North", "South"))) #, NAA_where = array(1, dim = c(2,2,6)))
#11 seasons each 1 month long except a 2 month interval in the model (June,July)
seasons = c(rep(1,5),2,rep(1,5))/12
basic_info$fracyr_seasons <- seasons
#each age other than 1 (recruitment) for north stock can be in either region on Jan 1 
basic_info$NAA_where <- array(1, dim = c(2,2,8))
basic_info$NAA_where[1,2,1] = 0 #stock 1, age 1 can't be in region 2 
basic_info$NAA_where[2,1,] = 0 #stock 2, any age can't be in region 1 (stock 2 doesn't move) 
#average recruitment over years 2000+ for SSB40 BRPs
basic_info$XSPR_R_avg_yrs <- which(temp$years>1999)
basic_info$XSPR_R_opt <- 2 #use average of recruitments (random effects), not expected/predicted given last time step
#############################################

#############################################
#ecov
north_bt <- read.csv(here("data","bsb_bt_temp_nmab_1959-2022.csv"))
south_bt <- read.csv(here("data","bsb_bt_temp_smab_1959-2022.csv"))
ecov <- list(label = c("North_BT","South_BT"))
ecov$mean <- cbind(north_bt[,'mean'], south_bt[,'mean'])
ecov$mean <- t(t(ecov$mean) - apply(ecov$mean,2,mean))
ecov$logsigma <- log(cbind(north_bt[,'se'], south_bt[,'se']))
ecov$year <- north_bt[,'year']
ecov$use_obs <- matrix(1, NROW(ecov$mean),NCOL(ecov$mean))
ecov$process_model <- "ar1"
ecov$process_mean_vals <- apply(ecov$mean, 2, mean)

ecov_0 <- ecov 
ecov_0$recruitment_how <- matrix(c("none","none","none","none"), 2,2)
ecov_1 <- ecov_2 <- ecov_0
ecov_1$recruitment_how[1,1] <- "controlling-lag-0-linear" #north
ecov_2$recruitment_how[2,2] <- "controlling-lag-0-linear" #south
ecov_3 <- ecov_1
ecov_3$recruitment_how[2,2] <- "controlling-lag-0-linear" #both

ecov_4 <- ecov_5 <- ecov
ecov_4$M_how <- ecov_5$M_how <- array("none", dim = c(2,2,8,2)) #n_ecov, n_stocks, n_ages, n_regions
ecov_4$M_how[1,1,1,1] <- "lag-0-linear" #north

ecov_5$M_how[2,2,1,2] <- "lag-0-linear" #south
ecov_6 <- ecov_4
ecov_6$M_how[2,2,1,2] <- "lag-0-linear" #both

#############################################

#############################################
#NAA_re
NAA_re = list(sigma = list("rec+1","rec+1"), cor = list("2dar1","2dar1"), N1_model = rep("equilibrium",2))

#Wasn't decoupled in RT
#NAA_re$decouple_recruitment <- FALSE

# set fixed values for NAA re for stock 1 in south on Jan 1
NAA_re$sigma_vals <- array(1, dim = c(2,2,temp$data$n_ages))
NAA_re$sigma_vals[1,2,2:temp$data$n_ages] <- 0.05 #2+ fixed to SCAA-ish for North fish in south on Jan 1.
NAA_re$sigma_map <- array(NA, dim = c(2,2,temp$data$n_ages))
NAA_re$sigma_map[1,1:2,1] <- 1
NAA_re$sigma_map[2,2,1] <- 3
NAA_re$sigma_map[1,1,2:temp$data$n_ages] <- 2
NAA_re$sigma_map[2,2,2:temp$data$n_ages] <- 4


#turn off estimation of AR1 cor parameters for north population in the south on Jan 1
x <- array(NA, dim = dim(temp$par$trans_NAA_rho))
x[1,1,1:2] <- 1:2
x[2,,1] <- 3
x[2,,2] <- 4
NAA_re$cor_map <- x

x <- array(NA, dim = dim(temp$par$trans_NAA_rho))
x[1,1,1:3] <- 1:3
x[2,2,1:3] <- 4:6
NAA_re$cor_map <- x
# input_1$map$trans_NAA_rho <- factor(x)
#############################################

#############################################
#move
move = list(stock_move = c(TRUE,FALSE), separable = TRUE) #north moves, south doesn't
move$must_move = array(0,dim = c(2,length(seasons),2))	
#if north stock in region 2 (south) must move back to region 1 (north) at the end of interval 5 right before spawning
move$must_move[1,5,2] <- 1 
move$can_move = array(0, dim = c(2,length(seasons),2,2))
move$can_move[1,c(1:4),2,1] <- 1 #only north stock can move and in seasons prior to spawning and after spawning
move$can_move[1,c(7:11),1,2] <- 1 #only north stock can move and in seasons prior to spawning and after spawning
move$can_move[1,5,2,] <- 1 #north stock can (and must) move in last season prior to spawning back to north 
mus <- array(0, dim = c(2,length(seasons),2,1))
mus[1,1:11,1,1] <- 0.02214863 #see here("2023.RT.Runs","transform_SS_move_rates.R") for how these numbers are derived.
mus[1,1:11,2,1] <- 0.3130358
move$mean_vals <- mus 
move$mean_model = matrix("stock_constant", 2,1)
#prior distribution on movement parameters 
move$use_prior <- array(0, dim = c(2,length(seasons),2,1))
move$use_prior[1,1,1,1] <- 1
move$use_prior[1,1,2,1] <- 1
move$prior_sigma <- array(0, dim = c(2,length(seasons),2,1))
move$prior_sigma[1,1,1,1] <- 0.2
move$prior_sigma[1,1,2,1] <- 0.2
#############################################

#############################################
#selectivity
sel <- list(n_selblocks = 8,
	model = rep(c("age-specific","logistic","age-specific"),
	c(2,2,4)))
sel$initial_pars <- list(
	rep(c(0.5,1),c(3,5)), #north comm
	rep(c(0.5,1),c(6,2)), #north rec
	c(5,1), #south comm
	c(5,1),	#south rec
	rep(c(0.5,1,1),c(1,1,6)), #north rec cpa
	rep(c(0.5,1),c(4,4)), #north vast
	rep(c(0.5,1,1),c(2,4,2)), #south rec cpa
	rep(c(0.5,1),c(1,7)) #south vast
)
sel$fix_pars <- list(
	4:8, #north comm
	7:8, #north rec
	NULL, #south comm
	NULL, #south rec
	2:8, #north rec cpa
	5:8, #north vast
	3:8, #south rec cpa
	2:8 #south vast
)
sel$re <- rep(c(
	"2dar1", "2dar1", "none", "ar1_y", "2dar1","none"), 
	c(1, 1, 2, 1, 1, 2))

#############################################

#############################################
#catch_info
catch_Neff <- temp$data$catch_Neff
catch_Neff[] <- 1000
catch_info <- list(catch_Neff = catch_Neff)
x <- temp$data$selblock_pointer_fleets
x[] <- rep(1:4, each = NROW(x))
catch_info$selblock_pointer_fleets = x
#############################################

#############################################
#index_info
index_Neff <- temp$data$index_Neff
index_Neff[] <- 1000
index_info <- list(index_Neff = index_Neff)
x <- temp$data$selblock_pointer_indices
x[] <- rep(5:8, each = NROW(x))
index_info$selblock_pointer_indices = x
#estimate obs sd scalar for Rec CPA in north and south
index_info$initial_index_sd_scale <- c(5,1,5,1)
index_info$map_index_sd_scale <- c(1,NA,2,NA)

#############################################

#############################################
#age_comp
age_comp = list(
	fleets = c("dir-mult","logistic-normal-miss0","logistic-normal-ar1-miss0","logistic-normal-ar1-miss0"), 
	indices = c("logistic-normal-miss0","dir-mult","logistic-normal-ar1-miss0","logistic-normal-ar1-miss0"))
#############################################

#############################################
#M random effects just on age 1
M_re_map <- array(NA, c(2,2,8))
M_re_map[1,1,1] <- 1
M_re_map[2,2,1] <- 2
M_sigma_vals <- diag(c(1,1))
M_sigma_map <- diag(1:2)
M_re_model <- matrix("none", 2,2)
M_re_model[1,1] <- "iid_y"
M_re_model[2,2] <- "iid_y"
M_list <- list(re_model = M_re_model, re_map = M_re_map, sigma_map = M_sigma_map, sigma_vals = M_sigma_vals)
#############################################

#############################################
#M random effects just on age 1 for North stock only (South stock variance goes to 0)
M_re_map <- array(NA, c(2,2,8))
M_re_map[1,1,1] <- 1
M_re_map[2,2,1] <- 2
M_sigma_vals <- diag(c(1,1))
M_sigma_map <- diag(c(1,NA))
M_re_model <- matrix("none", 2,2)
M_re_model[1,1] <- "iid_y"
# M_re_model[2,2] <- "iid_y"
M_list_2 <- list(re_model = M_re_model, re_map = M_re_map, sigma_map = M_sigma_map, sigma_vals = M_sigma_vals)
#############################################
