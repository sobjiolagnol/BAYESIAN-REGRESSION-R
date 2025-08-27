# Installer le package si besoin
install.packages("rrBLUP")
# Charger la librairie
library(rrBLUP)

# Charger les donn�es
telecat <- read.csv("telecat.csv")
telecat <- telecat[ , -1]  # Suppression de la colonne Unnamed

# S�parer Y et les variables explicatives
Y <- telecat$Y
X <- telecat[, setdiff(names(telecat), "Y")]

# Normaliser X
X <- scale(X)

# Fixer la seed et d�couper en training (100) / test (50)
set.seed(1234)
index_train <- sample(1:nrow(X), 100)
index_test <- setdiff(1:nrow(X), index_train)

Y_train <- Y[index_train]
Y_test <- Y[index_test]
X_train <- X[index_train, ]
X_test <- X[index_test, ]

# Appliquer RR-BLUP
resBLUP <- mixed.solve(Y_train, Z = X_train, X = matrix(1, nrow = 100, ncol = 1), method = "REML")

# Estimations
muchap <- resBLUP$beta
betachap <- resBLUP$u
Ve <- resBLUP$Ve
Vu <- resBLUP$Vu

# Affichage
print(muchap)
print(Ve)
print(Vu)

# Pr�dictions
predBLUP <- as.vector(muchap) + X_test %*% betachap

# Corr�lation pr�dictions vs observations
cor_pred <- cor(predBLUP, Y_test)
print(paste("Corr�lation pr�dictions vs observations :", round(cor_pred, 4)))

# Visualisation
plot(predBLUP ~ Y_test, main="Pr�dictions vs Observations", xlab="Y test", ylab="Pr�dictions")
abline(0, 1, col="red")


# Visualisation
plot(predBLUP ~ Y_test, main="Pr�dictions vs Observations", xlab="Y test", ylab="Pr�dictions")
abline(0, 1, col="red")

# S�lection des variables importantes
boxplot(betachap, main="Boxplot des coefficients")
bb <- boxplot(betachap, plot = FALSE)
subsetSelected <- function(resAlgo, varexpl, mini, maxi) {
  numselected <- c(which(resAlgo < mini), which(resAlgo > maxi))
  selected <- names(varexpl)[numselected]
  valeurs <- resAlgo[numselected]
  data.frame(selected, valeurs)
}
varselectedBLUP <- subsetSelected(betachap, as.data.frame(X_train), bb$stats[1, ], bb$stats[5, ])
print(varselectedBLUP)

# Chargement des packages n�cessaires
install.packages("MCMCpack")
install.packages("LearnBayes")
library(MCMCpack)
library(LearnBayes)
install.packages("MatrixModels")
library(MatrixModels)
rinvgamma <- function(n, a, b) {
  return(1 / rgamma(n, shape = a, rate = b))
}
# Fonction BayesA vue en TD
BayesA <- function(y,X,a,b,c,d,muinit,nbiter,nburn) {
  p <- ncol(X)
  n <- nrow(X)
  resbeta <- matrix(0, p, nbiter-nburn)
  ressigma2beta <- matrix(0, p, nbiter-nburn)
  resmu <- rep(0, nbiter-nburn)
  ressigma2eps <- rep(0, nbiter-nburn)
  
  beta <- rep(0, p)
  mu <- muinit
  sigma2beta <- rinvgamma(p, a, b)
  sigma2eps <- rinvgamma(1, c, d)
  
  for (iter in 1:nbiter) {
    Sigmabeta <- solve(t(X)%*%X/sigma2eps + diag(1/sigma2beta))
    beta <- as.numeric(rmnorm(1, Sigmabeta %*% t(X) %*% (y - mu)/sigma2eps, Sigmabeta))
    mu <- rnorm(1, mean(y - X %*% beta), sqrt(sigma2eps / n))
    for (j in 1:p) {
      sigma2beta[j] <- rinvgamma(1, a + 0.5, b + 0.5 * beta[j]^2)
    }
    sigma2eps <- rinvgamma(1, c + n/2, d + 0.5 * sum((y - mu - X %*% beta)^2))
    
    if (iter > nburn) {
      resbeta[, iter - nburn] <- beta
      ressigma2beta[, iter - nburn] <- sigma2beta
      resmu[iter - nburn] <- mu
      ressigma2eps[iter - nburn] <- sigma2eps
    }
  }
  return(list(resbeta, ressigma2beta, resmu, ressigma2eps))
}
install.packages("mnormt")
library(mnormt)


# Lancer la cha�ne MCMC
priora <- 2 ; priorb <- 1
priorc <- 2 ; priord <- 1
resBAYESA <- BayesA(Y_train, X_train, priora, priorb, priorc, priord, mean(Y_train), 2000, 1500)

# Estimations a posteriori
moybeta <- apply(resBAYESA[[1]], 1, mean)
moymu <- mean(resBAYESA[[3]])
moysigma2eps <- mean(resBAYESA[[4]])

# Pr�dictions sur test
predBAYESA <- moymu + X_test %*% moybeta
cor(predBAYESA, Y_test)

# Traces de quelques param�tres
par(mfrow=c(2,2))
plot(resBAYESA[[1]][50, ], type='l', main="Trace beta_50")
plot(resBAYESA[[2]][50, ], type='l', main="Trace sigma2_beta_50")
plot(resBAYESA[[3]], type='l', main="Trace mu")
plot(resBAYESA[[4]], type='l', main="Trace sigma2_eps")

# S�lection des variables importantes via boxplot
boxplot(moybeta, main="Boxplot des coefficients post�rieurs")
bb <- boxplot(moybeta, plot=FALSE)
subsetSelected <- function(resAlgo, varexpl, mini, maxi) {
  numselected <- c(which(resAlgo < mini), which(resAlgo > maxi))
  selected <- names(varexpl)[numselected]
  valeurs <- resAlgo[numselected]
  data.frame(selected, valeurs)
}
varselectedBAYESA <- subsetSelected(moybeta, as.data.frame(X_train), bb$stats[1,], bb$stats[5,])
print(varselectedBAYESA)

library(BLR)

resLASSO <- BLR(y=Y_train, XL=as.matrix(X_train),
                prior=list(
                  varE=list(df=2, S=1),
                  lambda=list(shape=10, rate=0.1, type="random", value=2)
                ),
                nIter=4000, burnIn=3000, saveAt="BLR_"
)

# Estimations
mu_hat <- resLASSO$mu
beta_hat <- resLASSO$bL
varE_hat <- resLASSO$varE
lambda_hat <- resLASSO$lambda
trace_lambda <- scan("BLR_lambda.dat")
plot(density(trace_lambda), col="blue", main="lambda (posterior)")
curve(dgamma(x, 10, 0.1), add=TRUE, col="red", lty=2, lwd=2)
legend("topright", c("Posterior", "Prior"), col=c("blue", "red"), lty=c(1,2))
predLASSO <- mu_hat + X_test %*% beta_hat
cor(predLASSO, Y_test)

boxplot(beta_hat)
bb <- boxplot(beta_hat, plot=FALSE)
varselectedLASSO <- subsetSelected(beta_hat, as.data.frame(X_train), -0.3, 0.3)

# Fonction SSVS 
selection_SSVS <- function(y, X, tau0, tau1, pi, niter) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)
  sigma2 <- 1
  gamma <- rep(1, p)
  gammamat <- matrix(0, p, niter)
  
  for (iter in 1:niter) {
    for (j in 1:p) {
      tauj2 <- ifelse(gamma[j] == 1, tau1^2, tau0^2)
      Vj <- 1 / (t(X[, j]) %*% X[, j] / sigma2 + 1 / tauj2)
      mj <- Vj * (t(X[, j]) %*% (y - X[, -j] %*% beta[-j])) / sigma2
      beta[j] <- rnorm(1, mean = mj, sd = sqrt(Vj))
      
      # Update gamma
      dens1 <- dnorm(beta[j], 0, tau1)
      dens0 <- dnorm(beta[j], 0, tau0)
      prob <- (pi * dens1) / (pi * dens1 + (1 - pi) * dens0)
      gamma[j] <- rbinom(1, 1, prob)
    }
    gammamat[, iter] <- gamma
  }
  return(list(gamma = gammamat))
}

# Lancer SSVS
resSSVS <- selection_SSVS(Y_train, X_train, tau0 = 0.1, tau1 = 1, pi = 0.5, niter = 2000)

# Proba a posteriori d'inclusion
gamma_mean <- apply(resSSVS$gamma, 1, mean)

# S�lection des variables avec proba > 0.5
selected_vars <- which(gamma_mean > 0.5)
selected_names <- colnames(X_train)[selected_vars]

# Affichage
print(selected_names)
plot(gamma_mean, type = "h", main = "Proba. inclusion SSVS", ylab = "P(gamma = 1)")
abline(h = 0.5, col = "red", lty = 2)

# Sparse model (pi = 0.2)
res_sparse <- selection_SSVS(Y_train, X_train, tau0 = 0.05, tau1 = 1, pi = 0.2, niter = 2000)
gamma_sparse <- apply(res_sparse$gamma, 1, mean)
sum(gamma_sparse > 0.5)
#4.3 variation de quelques hyper param�tres 

compare_ssvs <- function(pi_values, tau0_values, tau1_values, y, X, niter = 2000, seuil = 0.5) {
  resultats <- list()
  k <- 1
  for (pi in pi_values) {
    for (tau0 in tau0_values) {
      for (tau1 in tau1_values) {
        cat("??? SSVS avec pi =", pi, ", tau0 =", tau0, ", tau1 =", tau1, "\n")
        res <- selection_SSVS(y, X, tau0, tau1, pi, niter)
        gamma_mean <- apply(res$gamma, 1, mean)
        nb_selected <- sum(gamma_mean > seuil)
        resultats[[k]] <- list(pi=pi, tau0=tau0, tau1=tau1, nb_selected=nb_selected)
        k <- k + 1
      }
    }
  }
  return(do.call(rbind, lapply(resultats, as.data.frame)))
}
# Valeurs d'hyperparam�tres � tester
pi_vals <- c(0.2, 0.5, 0.8)
tau0_vals <- c(0.01, 0.1)
tau1_vals <- c(1)

# Comparaison
res_hyper <- compare_ssvs(pi_vals, tau0_vals, tau1_vals, Y_train, X_train)

# Affichage
print(res_hyper)
# 5.6 - Intervalles de confiance a posteriori des coefficients importants
top_vars <- order(abs(moybeta), decreasing = TRUE)[1:5]
for (j in top_vars) {
  beta_samples <- resBAYESA[[1]][j, ]
  ci <- quantile(beta_samples, c(0.025, 0.975))
  cat(paste("V", j, " : [", round(ci[1], 3), ", ", round(ci[2], 3), "]\n", sep=""))
}

# 5.7 - R�gression p�nalis�e non bay�sienne (LASSO, Ridge.)
install.packages("glmnet")
library(glmnet)

# LASSO
lasso_model <- cv.glmnet(X_train, Y_train, alpha = 1) 
coef(lasso_model, s = "lambda.min")
#  Mod�le de r�gression standard (lin�aire)
colnames(X_train)[1:20] 
X_selected <- X_train[, c("Sport.10", "Music.13", "Sport.15", "Film.8", "Film.10")]

lm_model <- lm(Y_train ~ ., data = as.data.frame(X_selected))
summary(lm_model)

# 5.9 analyse non lineaire 
plot(lm_model$residuals, main = "R�sidus")
abline(h = 0, col = "red")
# 6 ajout de la variable sexe 
# Effet fixe : colonne sexe du training (matrice 100 x 1)
XF_sexe <- matrix(telecat$sexe[index_train], ncol = 1)
colnames(XF_sexe) <- "Sexe"

#Lancer BLR avec :
  #XL : tes 161 covariables normalis�es

# XF : l'effet fixe � Sexe �
library(BLR)

set.seed(1234)
mod_BLR_sexe <- BLR(
  y      = Y_train,
  XF     = XF_sexe,               # effet fixe (non p�nalis�)
  XL     = as.matrix(X_train),    # LASSO bay�sien
  prior  = list(
    varE  = list(df = 2, S = 1),
    lambda = list(shape = 10, rate = 0.1, type = "random", value = 2)
  ),
  nIter  = 4000,
  burnIn = 3000,
  saveAt = "BLR_sexe_",
  thin   = 1
)

# Effet a posteriori de Sexe

beta_sexe <- mod_BLR_sexe$betaF["Sexe"]       
tau_sexe  <- mod_BLR_sexe$se.betaF["Sexe"]   
CI_sexe   <- beta_sexe + qnorm(c(0.025,0.975)) * tau_sexe

cat(sprintf("Effet a posteriori de Sexe : %.3f (IC95 %% : [%.3f ; %.3f])\n",
            beta_sexe, CI_sexe[1], CI_sexe[2]))
# 1) Effet fixe "Sexe" : convertis-le en scalaire
beta_sexe_hat <- as.numeric(mod_BLR_sexe$betaF)   

# 2) Pr�dictions sur le test
pred_BLR_sexe <- as.numeric(mod_BLR_sexe$mu) +       
  as.numeric(XF_sexe_test) * beta_sexe_hat +  
  X_test %*% mod_BLR_sexe$bL                 

str(mod_BLR_sexe$mu)    
str(mod_BLR_sexe$betaF)  
str(mod_BLR_sexe$bL)  

mu_hat        <- as.numeric(mod_BLR_sexe$mu)          # scalaire
beta_sexe_hat <- as.numeric(mod_BLR_sexe$betaF)       # scalaire
bL_hat        <- as.numeric(mod_BLR_sexe$bL)          # vecteur 161
# (50) vecteur du sexe
sexe_test_vec   <- as.numeric(XF_sexe_test)       

# (50) contribution LASSO
lasso_test_vec  <- as.numeric(X_test %*% bL_hat)  

# R�sum� des longueurs
length(sexe_test_vec) 
length(lasso_test_vec)

pred_BLR_sexe <- mu_hat +
  beta_sexe_hat * sexe_test_vec +
  lasso_test_vec

length(pred_BLR_sexe)   

cor_pred <- cor(pred_BLR_sexe, as.numeric(Y_test))
cat(sprintf("Corr�lation (avec Sexe) : %.3f\n", cor_pred))
# m�me seuil que pr�c�demment
bb  <- boxplot(mod_BLR_sexe$bL, plot = FALSE)
sel <- which(mod_BLR_sexe$bL < bb$stats[1] | mod_BLR_sexe$bL > bb$stats[5])
important <- colnames(X_train)[sel]
print(important)

# 9.2 Application avec EBLglmnet
install.packages("glmnet")
library(glmnet)

# Ajustement Elastic Net avec alpha = 0.5 (m�lange Ridge/LASSO)
fit_glmnet <- cv.glmnet(X_train, Y_train, alpha = 0.5)  # alpha = 1: LASSO, 0: Ridge

# Coefficients
coef(fit_glmnet, s = "lambda.min")

# Pr�diction
pred_enet <- predict(fit_glmnet, newx = X_test, s = "lambda.min")
cor_pred <- cor(pred_enet, Y_test)
cat(sprintf("Corr�lation pr�dictions (glmnet ElasticNet) : %.3f\n", cor_pred))

# 12
# Rappel de l'algorithme ABC classique (rejet simple)
Y_obs <- rnorm(100, mean = 2, sd = 1)

# Fonction ABC classique avec seuil epsilon
ABC_mean_var <- function(N, epsilon) {
  res_mu <- c()
  res_sigma2 <- c()
  
  for (i in 1:N) {
    mu_sim <- runif(1, -5, 10)             
    sigma2_sim <- 1 / rgamma(1, 2, 1)      
    
    Y_sim <- rnorm(length(Y_obs), mu_sim, sqrt(sigma2_sim))
    
    # Crit�re d'acceptation classique : moyenne et variance proches
    dist <- sqrt((mean(Y_sim) - mean(Y_obs))^2 + (var(Y_sim) - var(Y_obs))^2)
    if (dist < epsilon) {
      res_mu <- c(res_mu, mu_sim)
      res_sigma2 <- c(res_sigma2, sigma2_sim)
    }
  }
  
  return(data.frame(mu = res_mu, sigma2 = res_sigma2))
}

# Transformation du crit�re : test d'ad�quation loi simul�e vs observ�e
ABC_KS <- function(N, alpha = 0.05) {
  res_mu <- c()
  res_sigma2 <- c()
  
  for (i in 1:N) {
    mu_sim <- runif(1, -5, 10)
    sigma2_sim <- 1 / rgamma(1, 2, 1)
    
    Y_sim <- rnorm(length(Y_obs), mu_sim, sqrt(sigma2_sim))
    
    ks <- suppressWarnings(ks.test(Y_obs, Y_sim))
    
    if (ks$p.value > alpha) {  
      res_mu <- c(res_mu, mu_sim)
      res_sigma2 <- c(res_sigma2, sigma2_sim)
    }
  }
  
  return(data.frame(mu = res_mu, sigma2 = res_sigma2))
}
res <- ABC_KS(10000, alpha = 0.05)

hist(res$mu, main = "Posterior de mu (ABC-KS)", xlab = "mu")
hist(res$sigma2, main = "Posterior de sigma� (ABC-KS)", xlab = "sigma�")

