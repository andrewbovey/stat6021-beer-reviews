setwd("~/Desktop/Fall 2018/STAT6021/Project")
getwd()

library(tidyverse)
library(car)
library(perturb)

##### User Functions #####
prin.comp.coefficients <- function(y.vect, X0.mat) {
  # Collect parameters
  n <- dim(X0.mat)[1]
  k <- dim(X0.mat)[2]
  p <- k + 1
  # Unit-length scaling
  y.bar <- mean(y.vect)
  y.cent <- y.vect - y.bar
  SS.T <- sum(y.cent^2)
  y.vect.scl <- y.vect / sqrt(SS.T)
  X.mat.scl <- matrix(data=NA, nrow=n, ncol=k)
  x.bar.list <- numeric(length=k)
  css.list <- numeric(length=k)
  for (j in 1:k) {
    x.bar.list[j] <- mean(X0.mat[,j])
    xj.cent <- X0.mat[,j] - x.bar.list[j]
    css.list[j] <- sum(xj.cent^2)
    X.mat.scl[,j] <- xj.cent / sqrt(css.list[j])
  }
  # Calculate principal components and coefficient estimates
  XpX.mat.scl <- t(X.mat.scl) %*% X.mat.scl
  eig.out <- eigen(XpX.mat.scl)
  Lambda.mat <- diag(eig.out$values)
  Lambda.inv <- diag(1/eig.out$values)
  T.mat <- eig.out$vectors
  Z.mat <- X.mat.scl %*% T.mat
  Zpy.mat <- t(Z.mat) %*% y.vect.scl
  a.hat.scl <- Lambda.inv %*% Zpy.mat[1:j,]
  a.hat.PC.scl.list <- matrix(data=0, nrow=k, ncol=k)
  b.hat.PC.scl.list <- matrix(data=0, nrow=k, ncol=k)
  for (j in 1:k) {
    a.hat.PC.scl.list[1:j,j] <- a.hat.scl[1:j,]
    b.hat.PC.scl.list[,j] <- T.mat %*% a.hat.PC.scl.list[,j]
  }
  # Convert to the original scale and calculate MS.Res and R2.
  X.mat <- as.matrix(cbind(rep(x=1, times=n), X0.mat))
  grid.size <- dim(b.hat.PC.scl.list)[2]
  b.hat.PC.list <- matrix(data=0, nrow=p, ncol=grid.size)
  SS.Res.list <- numeric(length=grid.size)
  R2.list <- numeric(length=grid.size)
  for (iVAL in 1:grid.size) {
    b.hat.PC.list[1,iVAL] <- y.bar
    for (j in 1:k) {
      b.hat.PC.list[j+1,iVAL] <- b.hat.PC.scl.list[j,iVAL] / sqrt(css.list[j] / SS.T)
      b.hat.PC.list[1,iVAL] <- b.hat.PC.list[1,iVAL] - b.hat.PC.list[j+1,iVAL]*x.bar.list[j]
    }
    SS.Res.list[iVAL] <- sum((y.vect - X.mat %*% b.hat.PC.list[,iVAL])^2)
    R2.list[iVAL] <- 1 - SS.Res.list[iVAL] / SS.T
  }
  MS.Res.list <- SS.Res.list / (n-p)
  out.list <- list(b.hat.PC.list=b.hat.PC.list, MS.Res.list=MS.Res.list, R2.list=R2.list)
  return(out.list)
}
ridge.reg.coefficients <- function(y.vect, X0.mat, plot=TRUE, grid.size=25, grid.st=0.001, grid.fn=0.5) {
  # Collect parameters
  n <- dim(X0.mat)[1]
  k <- dim(X0.mat)[2]
  p <- k + 1
  # Unit-length scaling
  y.bar <- mean(y.vect)
  y.cent <- y.vect - y.bar
  SS.T <- sum(y.cent^2)
  y.vect.scl <- y.vect / sqrt(SS.T)
  X.mat.scl <- matrix(data=NA, nrow=n, ncol=k)
  x.bar.list <- numeric(length=k)
  css.list <- numeric(length=k)
  for (j in 1:k) {
    x.bar.list[j] <- mean(X0.mat[,j])
    xj.cent <- X0.mat[,j] - x.bar.list[j]
    css.list[j] <- sum(xj.cent^2)
    X.mat.scl[,j] <- xj.cent / sqrt(css.list[j])
  }
  # Calculate ridge trace diagram
  ridge.k.grid <- exp(seq(from=log(grid.st), to=log(grid.fn), length.out=grid.size))
  b.hat.R.scl.list <- matrix(data=0, nrow=k, ncol=grid.size)
  y.vect.aug <- rbind(y.vect.scl, matrix(data=0, nrow=k, ncol=1))
  for (iVAL in 1:grid.size) {
    ridge.k.val <- ridge.k.grid[iVAL]
    X.mat.aug <- rbind(X.mat.scl, sqrt(ridge.k.val)*diag(k))
    XpX.mat.aug <- t(X.mat.aug) %*% X.mat.aug
    Xpy.mat.aug <- t(X.mat.aug) %*% y.vect.aug
    XpX.inv.aug <- solve(XpX.mat.aug)
    b.hat.R.scl.list[,iVAL] <- XpX.inv.aug %*% Xpy.mat.aug
  }
  if (plot) {
    plot(ridge.k.grid, rep(x=0, times=grid.size), pch=3, cex=1, ylim=c(min(b.hat.R.scl.list), max(b.hat.R.scl.list)), xlab="ridge constant, k", ylab="fitted ridge regression coefficient", main = "Ridge trace diagram")
    abline(h=0, lty=1, lwd=1)
    for (j in 1:k) {
      lines(ridge.k.grid, b.hat.R.scl.list[j,], type="l", lty=1, lwd=3)
    }
  }
  # Convert to the original scale and calculate MS.Res and R2.
  X.mat <- as.matrix(cbind(rep(x=1, times=n), X0.mat))
  b.hat.R.list <- matrix(data=0, nrow=p, ncol=grid.size)
  SS.Res.list <- numeric(length=grid.size)
  R2.list <- numeric(length=grid.size)
  for (iVAL in 1:grid.size) {
    b.hat.R.list[1,iVAL] <- y.bar
    for (j in 1:k) {
      b.hat.R.list[j+1,iVAL] <- b.hat.R.scl.list[j,iVAL] / sqrt(css.list[j] / SS.T)
      b.hat.R.list[1,iVAL] <- b.hat.R.list[1,iVAL] - b.hat.R.list[j+1,iVAL]*x.bar.list[j]
    }
    SS.Res.list[iVAL] <- sum((y.vect - X.mat %*% b.hat.R.list[,iVAL])^2)
    R2.list[iVAL] <- 1 - SS.Res.list[iVAL] / SS.T
  }
  MS.Res.list <- SS.Res.list / (n-p)
  out.list <- list(ridge.k.grid=ridge.k.grid, b.hat.R.list=b.hat.R.list, MS.Res.list=MS.Res.list, R2.list=R2.list)
  return(out.list)
}

##### Data Read #####
beer <- read_csv("beerReviews.csv")

names(beer)
str(beer)

##### Data Types #####
beer$beerID <- as.factor(beer$beerID)
beer$brewerID <- as.factor(beer$brewerID)
beer$name <- as.factor(beer$name)
head(beer)
beer <- beer[,c(1:2,4,6:11)]
head(beer)

##### Fix the Beer Style Feature #####
style <- read_csv("style.csv", col_names = F)
names(style) <- c("style","style2","category")
beer <- merge(beer, style, type="left", all.x = T)
head(beer)
beer <- beer[,c(2:11)]
beer$ABV <- as.numeric(beer$ABV)

a <- which(duplicated(c(beer$beerID, beer$ABV)))
beer3 <- beer[-a,]
abvstyle <- beer3 %>% group_by(style2) %>% summarise(length(which(is.na(ABV))), mean(ABV,na.rm=T))

abvstyle <- abvstyle[,c(1,3)]

beer <- merge(beer, abvstyle, all.x = T)
beer[is.na(beer$ABV),c("ABV")] <- beer[is.na(beer$ABV),c("mean(ABV, na.rm = T)")]



# Subset the quantitative columns of beer, after grouping by beerID
beerQuant <- beer %>% group_by(beerID) %>% summarise(ABV = mean(ABV), appearance = mean(appearance), aroma = mean(aroma), palate = mean(palate), taste = mean(taste), overall = mean(overall))

beerQuant <- beerQuant[,c(2:7)]

beerQuant <- as.data.frame(beerQuant)
head(beerQuant)

##### Multicollinearity #####

beerQuant.lm <- lm(overall ~ ., data = as.data.frame(beerQuant))
summary(beerQuant.lm)

vif(beerQuant.lm)

colldiag(beerQuant[,c(1:5)], scale=T, center=TRUE, add.intercept=FALSE)

cor(beerQuant[,c(1:5)])
cor(beerQuant2[,c(1:5)])

##### Unit Length Scaled #####
beerQuantu <- beerQuant2
beerQuantu <- as.data.frame(beerQuantu)
beerQuantu <- apply(beerQuantu,2,function(x) (x - mean(x)) / sd(x) )
beerQuantu.lm <- lm(overall ~ ., data = as.data.frame(beerQuantu))
summary(beerQuantu.lm)

vif(beerQuantu.lm)
colldiag(beerQuantu, scale=F, center=F, add.intercept=FALSE)
cor(beerQuantu[,c(1:5)])


##### Ridge #####
library(glmnet)
library(broom)

lambdas <- 10^seq(3,-3,by=-0.1)
lambdas

ridge <- glmnet(as.matrix(beerQuant[,1:5]), beerQuant[,6], alpha = 0, lambda =lambdas)
ridge.cv <- cv.glmnet(as.matrix(beerQuant[,1:5]), beerQuant[,6], alpha = 0, lambda =lambdas)

plot(ridge.cv)
summary(ridge.cv)

ridge.reg.coefficients(as.matrix(beerQuant[,6]),as.matrix(beerQuant[,c(1:5)]))

##### Principle Components #####

pcaCharts <- function(x) {
  x.var <- x$sdev ^ 2
  x.pvar <- x.var/sum(x.var)
  print("proportions of variance:")
  print(x.pvar)
  
  par(mfrow=c(2,2))
  plot(x.pvar,xlab="Principal component", ylab="Proportion of variance explained", ylim=c(0,1), type='b')
  plot(cumsum(x.pvar),xlab="Principal component", ylab="Cumulative Proportion of variance explained", ylim=c(0,1), type='b')
  screeplot(x)
  screeplot(x,type="l")
  par(mfrow=c(1,1))
}
pcaCharts(prcomp(beerQuant[,c(1:5)],scale. = T, center = T))

prcomp(beerQuant[,c(1:5)],scale. = T, center = T)

summary(prcomp(beerQuant[,c(1:5)],scale. = T, center = T))
       

##### Influence Diagnostics #####

influence <- influence.measures(beerQuant.lm)
influence.measures(beerQuantu.lm)

library(olsrr)

# Cooks D plot
ols_plot_cooksd_chart(beerQuant.lm)
# 40451, 35855, 35921, 52202, 45939, 54100, 47350, 52203, 35857, 52314, 44802, 37474, 99009

# DF Betas plot
ols_plot_dfbetas(beerQuant.lm)

# DF Fits plot
ols_plot_dffits(beerQuant.lm)

# Residuals pllot
ols_plot_resid_lev(beerQuant.lm)
# 103126 stands out


##### Influence Points #####
beer.inf <- covratio(beerQuant.lm)
beer.inf

# The associated cutoff value for identifying an         | 
# influential observation is                             |

n <- dim(beerQuant)[1]
p <- length(coefficients(beerQuant.lm))
cut.inf.lo <- 1 - 3*p/n
cut.inf.lo
cut.inf.hi <- 1 + 3*p/n
cut.inf.hi

# The results from comparing against the cutoff values   |
# are displayed as follows                               |

cbind(beer.inf, (beer.inf < cut.inf.lo) | (beer.inf > cut.inf.hi))

length(rownames(beerQuant[(beer.inf < cut.inf.lo) | (beer.inf > cut.inf.hi),]))

# Eliminate the points outside of the cutoff range for covratio
beerQuant2 <- beerQuant[(beer.inf >= cut.inf.lo) & (beer.inf <= cut.inf.hi),]
beerQuant.lm <- lm(overall ~ -aroma, data = as.data.frame(beerQuant2))

# Update Cook's D plot with influence points removed
ols_plot_cooksd_bar(beerQuant2.lm)
# c(34724, 76429, 55564, 83450, 47698, 48464, 27054, 14289)


beerQuant2[c(34724, 76429, 55564, 83450, 47698, 48464, 27054, 14289),]



##### Testing the Model #####
#Generating index numbers for training and test data creation
smp = sample(nrow(beerQuant2)*0.8)

#Creating test and training datasets
data_tr = beerQuant2[smp,]
data_tst = beerQuant2[-smp,]

overall = data_tst[,6]
data_tst  = data_tst[,-6]

#Fitting the model, and predicting for test dataset
model = lm(overall~ ABV+appearance+palate+taste, data = data_tr) # Includes aroma
model2 = lm(overall~ ABV+appearance+palate+taste+aroma, data = data_tr) # Excludes aroma

pred_val = predict(model, newdata = data_tst)
pred_val2 = predict(model2, newdata = data_tst)

data_tst$predicted_overall = pred_val
data_tst$predicted_overall2 = pred_val2
data_tst$actual_overall = overall

summary(model)
summary(model2)

# MSE for both models
mean((data_tst$actual_overall - data_tst$predicted_overall)^2)
mean((data_tst$actual_overall - data_tst$predicted_overall2)^2)

r2 = cor( data_tst$actual_overall , data_tst$predicted_overall )^2
r2  # R-squared for first model
r22 = cor( data_tst$actual_overall , data_tst$predicted_overall2 )^2
r22 # R-squared for second model (no aroma)

# Generate a citation for the package
citation(package = "olsrr")
