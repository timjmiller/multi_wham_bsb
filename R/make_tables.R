
mod_names <- paste0("$M_{", 0:13,"}$")
col_names <- c("Model", "Temperature effect (north)", "Temperature effect (south)", "$M$ at age 1 random effects")

mod_table <- matrix(nrow = length(mod_names), ncol = length(col_names), dimnames = list(mod_names,col_names))

mod_table[,1] <- mod_names
mod_table[,2:4] <- "--"
mod_table[,4] <- rep(c("none", "time-varying"), each = 7)
mod_table[c(2,4),2] <- "Recruitment"
mod_table[c(3,4),3] <- "Recruitment"
mod_table[c(5,7),2] <- "$M$ at age 1"
mod_table[c(6,7),3] <- "$M$ at age 1"
mod_table[8:14,2:3] <- mod_table[1:7,2:3]


fits <- readRDS(file.path("results","fits_no_M_re.RDS"))
fits_M_re <- readRDS(file.path("results","fits_M_re_better.RDS"))

# AIC

aic1 <- sapply(fits, function(x) {
	2*c(x$opt$obj + length(x$opt$par), sapply(x$peels, function(y) y$opt$obj + length(y$opt$par)))
})
aic2 <- sapply(fits_M_re, function(x) {
	2*c(x$opt$obj + length(x$opt$par), sapply(x$peels, function(y) y$opt$obj + length(y$opt$par)))
})

diff_aic <- apply(cbind(aic1,aic2),1, function(x) x - min(x))
aic_wts <- apply(aic,2, function(x) exp(-x/2)/sum(exp(-x/2)))
rownames(diff_aic) <- rownames(aic_wts) <- mod_table[,1]
colnames(diff_aic) <- colnames(aic_wts) <- paste0("Peel ", 0:7)
saveRDS(diff_aic, file.path("results","diff_aic.RDS"))
saveRDS(aic_wts, file.path("results","aic_wts.RDS"))


rho_ssb <- rbind(t(sapply(fits, function(x) mohns_rho(x)$SSB)),
 t(sapply(fits_M_re, function(x) mohns_rho(x)$SSB)))
rho_Fbar <- rbind(t(sapply(fits, function(x) mohns_rho(x)$Fbar)),
 t(sapply(fits_M_re, function(x) mohns_rho(x)$Fbar)))
rho_R <- rbind(t(sapply(fits, function(x) c(mohns_rho(x)$naa[1,1,1],mohns_rho(x)$naa[2,2,1]))),
 t(sapply(fits_M_re, function(x) c(mohns_rho(x)$naa[1,1,1],mohns_rho(x)$naa[2,2,1]))))

rho_table <- cbind.data.frame(round(rho_ssb,3), round(rho_Fbar,3), round(rho_R,3))
rho_table <- cbind(mod_names, rho_table)
saveRDS(rho_table, file.path("results", "rho_table.RDS"))


# mohns_rho <- function (model) {
#     npeels = length(model$peels)
#     data <- model$env$data
#     ny = data$n_years_model
#     na = data$n_ages
#     if (npeels) {
#         rho = list()
#         rho$SSB <- rho$Fbar <- numeric()
#         for (i in 1:data$n_stocks) rho$SSB[i] <- mean(sapply(1:npeels,
#             function(x) model$peels[[x]]$rep$SSB[ny - x, i]/model$rep$SSB[ny -
#                 x, i] - 1))
#         for (i in 1:data$n_regions) rho$Fbar[i] <- mean(sapply(1:npeels,
#             function(x) {
#                 mean(model$peels[[x]]$rep$Fbar[ny - x, i])/mean(model$rep$Fbar[ny -
#                   x, i]) - 1
#             }))
#         rho$naa <- array(NA, c(data$n_stocks, data$n_regions,
#             data$n_ages))
#         for (s in 1:data$n_stocks) for (r in 1:data$n_regions) for (a in 1:data$n_ages) if (data$NAA_where[s,
#             r, a]) {
#             rho$naa[s, r, a] = mean(sapply(1:npeels, function(x) model$peels[[x]]$rep$NAA[s,
#                 r, ny - x, a]/model$rep$NAA[s, r, ny - x, a] -
#                 1))
#         }
#         dimnames(rho$naa)[[3]] = c("R", paste0("N", model$ages.lab[2:na]))
#         rho$FAA <- matrix(NA,nrow = data$n_fleets, ncol = data$n_ages)
#         for(i in 1:data$n_fleets) for(j in 1:data$n_ages) rho$FAA[i,j] <- mean(sapply(1:npeels, function(x) model$peels[[x]]$rep$FAA[i,ny - x, j]/model$rep$FAA[i, ny - x, j] - 1))
#         dimnames(rho$FAA)[[1]] = model$input$fleet_names
#         dimnames(rho$FAA)[[2]] = c(paste0("F", model$ages.lab[1:na]))
#         rho$FAA_by_region <- matrix(NA,nrow = data$n_regions, ncol = data$n_ages)
#         for(i in 1:data$n_regions) for(j in 1:data$n_ages) rho$FAA_by_region[i,j] <- mean(sapply(1:npeels, function(x) {
#         	model$peels[[x]]$rep$FAA_by_region[i, ny - x, j]/model$rep$FAA_by_region[i, ny - x, j] - 1
#       	}))
#         dimnames(rho$FAA_by_region)[[1]] = model$input$region_names
#         dimnames(rho$FAA_by_region)[[2]] = c(paste0("F", model$ages.lab[1:na]))
#         return(rho)
#     }
#     else stop("There are no peels in this model")
# }
