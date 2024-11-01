# Q1
# A

# No R code needed

# B
# Generate a 111 x 111 real symmetric matrix
p <- 111
M <- matrix(rnorm(p^2), nrow=p)
M <- (M + t(M)) / 2

# Verify that the inverse of the matrix of eigenvectors is its transpose
eigen_matrix <- eigen(M)$vectors
inverse_eigen_matrix <- solve(eigen_matrix)

all.equal(inverse_eigen_matrix, t(eigen_matrix))

# C

# No R code needed


############################################################################################

# Q2
# A

rnorm_multi <- function(n, mu, sigma) {
  p <- length(mu)
  x <- matrix(rnorm(n * p), ncol = p)
  return(t(apply(x, 1, function(z) mu + chol(sigma) %*% z)))
}

# Code to test function

n <- 100
p <- 3
mu <- c(1, 2, 3)
sigma <- matrix(c(2, 0.5, 0.3, 0.5, 1, 0.4, 0.3, 0.4, 3), nrow = p)
data <- rnorm_multi(n, mu, sigma)
data


# B

# No R code needed

# C

generate_random_mean_cov <- function(p) {
  # Generate p-dimensional mean vector
  mu <- rnorm(p)
  
  # Generate p × p-dimensional random square matrix A with N(0,1) entries
  A <- matrix(rnorm(p^2), nrow = p, ncol = p)
  
  # Calculate covariance matrix
  Sigma <- A %*% t(A)
  
  return(list(mu = mu, Sigma = Sigma))
}

# Generate random mean vector and covariance matrix
result <- generate_random_mean_cov(p = 50)
mu <- result$mu
Sigma1 <- result$Sigma

# Simulate data matrix X
X <- rnorm_multi(n = 1000, mu = mu, sigma = Sigma1)

p=50

# Define matrix C
C <- cbind(rep(1, p), c(1:p)^2)

# Evaluate data-matrix Y = XC'
Y <- X %*% C


# Generate random mean vector and covariance matrix
result <- generate_random_mean_cov(p = 50)
mu <- result$mu
Sigma <- result$Sigma


# Plot the histogram of Y[1] including an overlay of the data quantiles vs the theoretical quantiles

# Plot the Histogram of Y[1]
hist(Y[,1], main = paste("Histogram of Y[1]", sep=""), xlab="Y", ylab="Density", freq = FALSE, ylim = c(0, 0.009))

# Plot the curve of the theoretical density which is a multivariate normal curve
curve(dnorm(x, mean = mean(Y[,1]), sd = sd(Y[,1])), add = TRUE, col = "red")

# Compute the Quantiles of the data and the Quantiles of the Theoretical data
data_quantiles <- quantile(Y[,1], probs = seq(0, 1, length.out = 100))
theoretical_quantiles <- qnorm(seq(0, 1, length.out = 100), mean = mean(Y[,1]), sd = sd(Y[,1]))

# Use par(new=True) to overlay a plot of the data quantiles vs the theoretical quantiles
par(new = TRUE)
# Include axes = False so the data is readable
plot(data_quantiles, theoretical_quantiles, axes = FALSE, xlab = '', ylab = '')
abline(a=0, b=1)

# Introduce a legend to guide the reader through the plot
legend("topleft", c("Histogram", "Theoretical Density", "Quantiles vs\n Theoretical Densities"), fill = c("grey", "red", "black"), cex = 0.65)


# Plot a histogram of Y[2] in the same format as Y[1]
hist(Y[,2], main = paste("Histogram of Y[2]", sep=""), xlab="Y", ylab="Density", freq = FALSE, ylim = c(0, 0.000012))
curve(dnorm(x, mean = mean(Y[,2]), sd = sd(Y[,2])), add = TRUE, col = "red")

data_quantiles <- quantile(Y[,2], probs = seq(0, 1, length.out = 100))
theoretical_quantiles <- qnorm(seq(0, 1, length.out = 100), mean = mean(Y[,2]), sd = sd(Y[,2]))

par(new = TRUE)
plot(data_quantiles, theoretical_quantiles, axes = FALSE, xlab = '', ylab = '')
abline(a=0, b=1)

legend("topleft", c("Histogram", "Theoretical Density", "Quantiles vs\n Theoretical Densities"), fill = c("grey", "red", "black"), cex = 0.65)

# Separate plots of the quantiles
plot(data_quantiles, theoretical_quantiles, xlab = "Data Quantiles [2]", ylab = "Theoretical Quantiles [2]", main = "Data Quantiles vs Theoretical Quantiles")
abline(a=0, b=1)


##################################################################################################################################

# Q3

# A
load("C:\\Users\\Padraig\\Downloads\\HWKA1.Rdata")
library(MASS)

# Set parameters
N <- 900
P <- 60
mu <- data$mu
Sigma <- data$Sigma
A <- data$A
V <- data$V

# Compute the theoretical mean and covariance of Y1 = AX1
mu_y = A %*% mu # calculate the theoretical mean of Y1
sigma = A %*% Sigma %*% t(A) # calculate the theoretical covariance of Y1

# Simulate a random sample from the multivariate normal distribution
X = mvrnorm(N, mu, Sigma) # simulate a random sample of X values

# Compute the sample mean and covariance of Y1 = AX1
y_hat = A %*% apply(X,2,mean) # calculate the sample mean of Y1
cov_hat = A %*% cov(X) %*% t(A) # calculate the sample covariance of Y1

# Create plots comparing the theoretical and sample mean and covariance
plot(mu_y, y_hat, xlab='Theoretical Mean of Y', ylab="Sample Mean of Y", main="Theoretical Mean vs Sample Mean")
abline(a=0, b=1, col = 'blue') # add a diagonal line to the plot for comparison

plot(sigma, cov_hat, xlab="Theoretical Covariance of Y", ylab = "Sample Covariance of Y", main="Theoretical Covariance vs Sample Covariance")
abline(a=0, b=1, col='blue') # add a diagonal line to the plot for comparison


# B
# Calculate the Mahalanobis Distanec as per the formula
Z = numeric(N)
for (i in 1:N){
  Zi = t(as.matrix((X[i,]-mu), nrow = 1, ncol = 60)) %*% solve(Sigma) %*% as.matrix((X[i,]-mu), nrow = 60, ncol = 1)
  Z[i] = Zi
}

# Create a formula of the plot
hist(Z, freq = FALSE, xlab = "Mahalanobis Distance", ylim = c(0,0.038))

# Plot the theoretical density function
lines(dchisq(seq(1,max(100), by=1), df = P), col = 'red', lwd=2)

# C

# Calculate Mahalanobis distance of V
d2 <- t(as.matrix((V-mu), nrow = 1, ncol = P)) %*% solve(Sigma) %*% as.matrix((V-mu), nrow = P, ncol = 1)

# Calculate p-value
p_value <- 1 - pchisq(d2, df = P)

# Print results
cat("Mahalanobis distance:", d2, "\n")
cat("p-value:", p_value, "\n")
if (p_value <= 0.05) {
  cat("Reject null hypothesis at 0.05 level of significance\n")
} else {
  cat("Fail to reject null hypothesis at 0.05 level of significance\n")
}


#############################################################################################################################

# Q4
data(iris)
m1 = apply(iris[1:50, 1:4], 2, mean)
m2 = apply(iris[51:100, 1:4], 2, mean)
m3 = apply(iris[101:150, 1:4], 2, mean)
Sigma = (cov(iris[1:50, 1:4]) + cov(iris[51:100, 1:4]) + cov(iris[101:150, 1:4]))/3 #covariance of the measurements

# (A)

mean_deviation <- function(X, Sigma) { 
  grp = X[,5]
  m1 = apply(X[grp == 'setosa', 1:4], 2, mean)
  m2 = apply(X[grp == 'versicolor', 1:4], 2, mean)
  m3 = apply(X[grp == 'virginica', 1:4], 2, mean)
  m = apply(X[,1:4],2, mean)
  dev2 = sum((m1-m)*solve(Sigma, m1-m)) + sum((m2-m)*solve(Sigma, m2-m)) + sum((m3-m)*solve(Sigma, m3-m))
  dev2
  
}

ns = 1000 #number permutations
devs = rep(0, ns)#where deviations for each permutation will be stored of length ns
for (k in 1:ns) {
  devs[k] = mean_deviation(X = cbind(iris[,1:4], iris[sample.int(n = 150, size = 150, replace = FALSE),5]), Sigma = (cov(iris[1:50, 1:4]) + cov(iris[51:100, 1:4]) + cov(iris[101:150, 1:4]))/3)
}
# (B)
hist(devs, breaks = 50, main = "Permutation Distribution", xlab="Deviation from the mean")

# (C)
dev_value = mean_deviation(iris, Sigma)

pval = length(devs[devs>=dev_value])/ns
pval



####################################################################################################
# Q5
# (A)

# Load the USJudgeRatings dataset
data(USJudgeRatings)

# Select the numerical columns of the dataset
data_matrix <- USJudgeRatings

# Compute the correlation matrix
num_vars <- ncol(data_matrix)
R <- matrix(0, nrow = num_vars, ncol = num_vars)
for (i in 1:num_vars) {
  for (j in 1:num_vars) {
    if (i == j) {
      R[i,j] <- 1
    } else if (j > i) {
      x <- data_matrix[, i]
      y <- data_matrix[, j]
      n <- length(x)
      x_mean <- sum(x) / n
      y_mean <- sum(y) / n
      cov_xy <- sum((x - x_mean) * (y - y_mean)) / (n - 1)
      x_var <- sum((x - x_mean)^2) / (n - 1)
      y_var <- sum((y - y_mean)^2) / (n - 1)
      R[i,j] <- cov_xy / sqrt(x_var * y_var)
      R[j,i] <- R[i,j]
    }
  }
}

# (B)
# Spectral decomposition of the correlation matrix
eigen_result <- eigen(R)
G <- eigen_result$vectors
G_4 <- G[,1:4]

reduced_dimensional_df = as.matrix(R) %*% G_4

y = matrix(NA, nrow=nrow(R), ncol=4)
x = reduced_dimensional_df

mu_vector = unname(as.matrix(apply(G_4, 2, mean), nrow=1))
sd_vector = unname(as.matrix(apply(G_4, 2, sd), nrow=1))

for (i in 1:nrow(reduced_dimensional_df)){
  y[i,]  =as.numeric(((as.numeric(x[i,]) - mu_vector / sd_vector)))
}

lm.1 = lm(y~G_4-1)
summary(lm.1)


yGi <- t(G) %*% y
bi = lm.1$coefficients

plot(bi, yGi[1:4,], main = "Bi vs G'yi")
abline(a=0, b=1)

# C
# Compute the linear model residual sum of squares
RSS <- sum(resid(lm.1)^2)

# Compute the sum of the last 8 eigenvalues of the correlation matrix
sum_eigenvalues <- sum(eigen_result$values[5:12])

RSS
sum_eigenvalues



#######################################################################################################

# Q6
# (A)
library(dplyr)
set.seed(123)  # for reproducibility
n_pixels <- 10000

# Randomly sample pixel indices
pixels <- sample(length(ncd2), n_pixels)

if (max(pixels) > length(ncd2)) {
  stop("Invalid pixel index")
}

# Create a matrix X with rows corresponding to sampled pixel indices and columns corresponding to the sampled time-points
X <- matrix(0, n_pixels, 481)
for (i in 1:n_pixels) {
  row_index <- (pixels[i]) %% nrow(ncd2)+1
  col_index <- (pixels[i]) %% ncol(ncd2)+1
  # check that row_index and col_index are within valid range
  if (row_index <= nrow(ncd2) & col_index <= ncol(ncd2)) {
    X[i, ] <- ncd2[row_index, col_index, 1:481]
  } else {
    # handle out-of-bounds case
    cat(row_index, col_index)
  }
}

# (B)

# Remove time-point 327
X <- X[,-327]

# Compute correlation matrix
R <- cor(X)

# Make image of correlation matrix and overlay a plot of the temporal mean and standard deviation
image(R)
col_mean <- apply(X,2,mean)
col_sd <- apply(X, 2, sd)

par(new=TRUE)
plot(col_mean, col = 'blue', axes = FALSE)
par(new=TRUE)
plot(col_sd, col = 'red', axes = FALSE)

# Plot temporal mean and standard deviation separate
mean_X <- rowMeans(X)
sd_X <- apply(X, 1, sd)
plot(mean_X, type = "l", xlab = "Pixel index", ylab = "Temporal mean")
plot(sd_X, type = "l", xlab = "Pixel index", ylab = "Temporal standard deviation")


# Plot the temporal mean and standard deviation of X
plot(1:480, colMeans(R), type = "l", col = "blue", ylim = range(R-0.2), main = " temporal mean and standard deviation")
lines(1:480, colMeans(R) + sd(R), col = "red")
lines(1:480, colMeans(R) - sd(R), col = "red")
# legend("top", lty=1, col=c('blue', 'red'), legend=c("Mean", "Standard Deviation"))
legend("bottomright",legend = c("Mean","Standard Deviation"),col=c('blue', 'red'), lty=1)


# Compute eigen-decomposition of the correlation matrix
eigen_result <- eigen(R)

# Plot percent of variance explained by number of components
plot(cumsum(eigen_result$values)/sum(eigen_result$values)*100, xlab = "Number of components", ylab = "Percent variance explained")

# (C)


# Get the loading vectors
loading_vectors <- eigen_result$vectors

# Compute the absolute values of the loading vectors
abs_loading_vectors <- abs(loading_vectors)

# Adjust the sign of the loading vectors
loadings <- eigen_result$vectors
max_idx <- apply(abs(loadings), 2, which.max)
loadings <- loadings * sign(loadings[max_idx,])

# Center X by subtracting column means
X_centered <- scale(X, center = TRUE, scale = FALSE)

# Compute the scores for the first two principal components
scores <- X_centered %*% loading_vectors[, 1:2]

# Make a scatterplot of the scores
plot(scores, col = "blue", pch = 20, xlab = "PC1", ylab = "PC2", main = "PC Scores")


# (D)
# Standardize the data matrix X==
X_standardized <- scale(X, center = TRUE, scale = TRUE)

# Use the first 7 loading vectors to project the data
G_sub <- loading_vectors[, 1:7]
X_projected <- X_standardized %*% G_sub

par(mar=c(5,5,2,2))  # adjust the plot margins

# plot the first vector
plot(loading_vectors[,1], type="l", col="blue", ylim=c(-0.12, 0.12), xlab="Pixel Index", ylab="Loading Value", main = "Loading Vectors")
ncol(loading_vectors)

lines(loading_vectors[,2], col="red")
lines(loading_vectors[,3], col="green")
lines(loading_vectors[,4], col="#fc03e7")
lines(loading_vectors[,5], col="black")
lines(loading_vectors[,6], col="#00FFFF")
lines(loading_vectors[,7], col="#f4fc03")



# (E)

# Plot the images of the mean and standard deviation of the ncd2 datatset

mean_output = matrix(data = NA, nrow = 1080, ncol = 1280)
for(row in 1:1080){
  for(col in 1:1280){
    mean_output[row, col] = mean(ncd2[row, col,])
  }
}

sd_output = matrix(data = NA, nrow = 1080, ncol = 1280)
for(row in 1:1080){
  for(col in 1:1280){
    sd_output[row, col] = sd(ncd2[row, col,])
  }
}

par(mfrow=c(3,3))
image(mean_output, main = 'mean intensity', axes= FALSE)
image(sd_output, main = 'sd-intensity', axes= FALSE)


# Plot the images of layers 1 through 7 with an overlay of the corresponding loading vectors

layer<-ncd2[ , ,1]
image(layer, main = "1", axes= FALSE)
par(new=TRUE)
abline(h = 0.5, lty=2)
par(new=TRUE)
plot(loading_vectors[,1], type = "l", axes=FALSE, col="blue", lwd=2)

layer<-ncd2[ , ,2]
image(layer, main = "2", axes= FALSE)
par(new=TRUE)
abline(h = 0.5, lty=2)
par(new=TRUE)
plot(loading_vectors[,2], type = "l", axes=FALSE, col="red", lwd=2)

layer<-ncd2[ , ,3]
image(layer, main = "3", axes= FALSE)
par(new=TRUE)
abline(h = 0.5, lty=2)
par(new=TRUE)
plot(loading_vectors[,3], type = "l", axes=FALSE, col="green", lwd=2)

layer<-ncd2[ , ,4]
image(layer, main = "4", axes= FALSE)
par(new=TRUE)
abline(h = 0.5, lty=2)
par(new=TRUE)
plot(loading_vectors[,4], type = "l", axes=FALSE, col="#fc03e7", lwd=2)

layer<-ncd2[ , ,5]
image(layer, main = "5", axes= FALSE)
par(new=TRUE)
abline(h = 0.5, lty=2)
par(new=TRUE)
plot(loading_vectors[,5], type = "l", axes=FALSE, col="black", lwd=2)

layer<-ncd2[ , ,6]
image(layer, main = "6", axes= FALSE)
par(new=TRUE)
abline(h = 0.5, lty=2)
par(new=TRUE)
plot(loading_vectors[,6], type = "l", axes=FALSE, col="#00FFFF", lwd=2)

layer<-ncd2[ , ,7]
image(layer, main = '7', axes= FALSE)
par(new=TRUE)
abline(h = 0.5, lty=2)
par(new=TRUE)
plot(loading_vectors[,7], type = "l", axes=FALSE, col='#79787a', lwd=2)

# End of Assignment Script
#########################################################################################