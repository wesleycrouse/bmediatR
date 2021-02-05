#devtools::install_github("wesleycrouse/bmediatR")
source("mediation_bf.R")

#load data
y <- readRDS("data/Y_DO_example.RDS")
M <- readRDS("data/M_DO_example.RDS")
X <- readRDS("data/X_DO_example.RDS")
Z <- readRDS("data/covar_DO_example.RDS")

#call mediation_bf
results <- mediation_bf(y, M, X, Z)

#plot Bayes factors
plot(results$lnBF)

#significance by permutation
n_perms <- 1000

max_lnBF <- rep(NA, n_perms)

for (i in 1:n_perms){
  print(i)
  M_perm <- M[sample(1:nrow(M)),]
  results_perm <- mediation_bf(y, M_perm, X, Z, verbose=F)
  max_lnBF[i] <- max(results_perm$lnBF)
}

#save(max_lnBF, file="max_lnBF.RData")
#load("max_lnBF.RData")

evd_pars <- as.numeric(evir::gev(max_lnBF)$par.est)
threshold <- evir::qgev(0.95, xi=evd_pars[1], sigma=evd_pars[2], mu=evd_pars[3])
abline(h=threshold, col="red")
