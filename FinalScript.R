###################################################################################################################################
# Student ID: 119378796
# Name : Padraig O'Brien

###################################################################################################################################
# Q1
library(av)
library(magick)

# Specify the input movie file
mp4file <- "C:/Users/Padraig/Desktop/DSA4_SEM2/Multivariate/CA2/Materials/waveclip.mp4"
# Specify the destination directory for the output jpg files
junkd <- "C:/Users/Padraig/Desktop/DSA4_SEM2/Multivariate/CA2/Materials/waveclipimages/"

# Use av to create jpg files of individual movie frames
o <- av_video_images(mp4file, destdir = junkd)
setwd(junkd)

# Get the file names of the jpg files
jpgfilenames <- list.files(junkd, pattern = ".jpg$")

o <- image_read(paste0(junkd, "/", jpgfilenames[1]))

# initialize empty arrays for each color channel
red_channel <- array(dim = c(154, 140, 1790))
green_channel <- array(dim = c(154, 140, 1790))
blue_channel <- array(dim = c(154, 140, 1790))

# loop through each frame and extract the color channel data
for (k in 1:length(jpgfilenames)) {
  
  # read the k'th frame using magick
  o <- image_read(paste0(junkd, "/", jpgfilenames[k]))
  
  # extract the color channel data
  red_data <- as.numeric(image_data(o)[,,1])
  green_data <- as.numeric(image_data(o)[,,2])
  blue_data <- as.numeric(image_data(o)[,,3])
  
  # store the data in the corresponding arrays
  red_channel[,,k] <- array(red_data, dim = c(154, 140))
  green_channel[,,k] <- array(green_data, dim = c(154, 140))
  blue_channel[,,k] <- array(blue_data, dim = c(154, 140))
}

# save the color channel arrays as an .Rdata file
save(list = c("red_channel", "green_channel", "blue_channel"), file = "movie_data.Rdata")

# Reshape red_channel into a 2D array
red_matrix <- matrix(red_channel, nrow = 154 * 140, ncol = 1790)

# Compute the covariance matrix of red_matrix
red_cov <- cov(red_matrix)

eig <- eigen(red_cov)

# Extract eigenvalues and eigenvectors
eig_values <- eig$values ; eig_vectors <- eig$vectors

# Sort eigenvalues and eigenvectors in decreasing order of eigenvalue
order <- order(-eig_values)
eig_values <- eig_values[order]
eig_vectors <- eig_vectors[, order]; eig_diag <- diag(eig_values)

# Compute the spectral decomposition
spectral_decomp <- eig_vectors %*% eig_diag %*% t(eig_vectors)

J <- 10
plot(eig$values[1:J], type = "b", xlab = "Eigenvalue index", ylab = "Eigenvalue", 
     col = rainbow(J), pch = 19, main = "J eigenvalues of Σ_T")


###############################################################################################################################################
# Q2
# Compute 3 dim colour space
# Load the waves.Rdata file
load("C:/Users/Padraig/Desktop/DSA4_SEM2/Multivariate/CA2/Materials/Waves.Rdata")

red <- data$r ; green <- data$g ; blue <- data$b

col_space <- cbind(red, green, blue)

# var explained by each component
cov_mat <- cov(col_space); pca_res <- prcomp(col_space, scale. = TRUE)

var_expl <- pca_res$sdev^2/sum(pca_res$sdev^2)

print(var_expl)
# Normalise loading vector and proj
loading_vectors <- pca_res$rotation
first_loading_vector <- loading_vectors[,1]
normalized_loading_vector <- first_loading_vector / max(abs(first_loading_vector))
print(normalized_loading_vector)

proj <- (loading_vectors[,1][1] * red[,,T/2]) + (loading_vectors[,1][2] * green[,,T/2]) + (loading_vectors[,1][3] * blue[,,T/2])

# image plots
T = dim(red)[3]

par(mfrow=c(2,2))
image(red[,,T/2], main = "Red Colour Channel") ; image(blue[,,T/2], main = "Blue Colour Channel")
image(green[,,T/2], main = "Green Colour Channel") ; image(as.matrix(proj), main = "PCA1 Projection")
# Make Z
# Get the number of pixels and time points 
N <- 128*128
T <- 1790

# Project each pixel time-series onto the first loading vector
proj_data <- array(col_space %*% as.matrix(normalized_loading_vector), dim = c(16384,1790))

# Adjust each pixel time-series for a linear trend in time
trend <- seq(1, T) / T
trend_mat <- matrix(rep(trend, N), ncol=T, byrow=TRUE)

Z <- proj_data - trend_mat

###########################################
# B

# Temporal Covariance
temporal_cov <- cov(Z)
eigen_res <- eigen(temporal_cov)
# Relationship
TSS = sum(eigen_res$values)
var_explained = sum(eigen_res$values[1:5])
TSS - var_explained; sum(eigen_res$values[6:length(eigen_res$values)])

(TSS - var_explained) - (sum(eigen_res$values[6:length(eigen_res$values)]))

Z_fft <- apply(Z, 2, fft)
Z_fft <- Re(Z_fft) # My computer is acting up and outputting complex numbers where there is none so I had to run this
# compute variance explained by each Fourier component
variance_explained <- apply(Z_fft^2, 2, sum)/nrow(Z)

# sort components by the amount of variance explained
sorted_indices <- order(variance_explained, decreasing=TRUE)
variance_explained_sorted <- variance_explained[sorted_indices]

# PCA
temporal_eigenvalues <- eigen_res$values ; temporal_eigenvectors <- eigen_res$vectors

# Plot
# plot fraction of variance explained as a function of the number of components
cumulative_variance_explained_fourier <- cumsum(variance_explained_sorted)/sum(variance_explained_sorted)
par(mfrow = c(1,1))
plot(1:ncol(Z), cumulative_variance_explained_fourier, type='l', ylim = c(0,1)
     , xlab='Number of components', ylab='Fraction of variance explained', main='Fourier vs. PCA % Variance Explained', lwd = 3)
cumsum_eigvals <- cumsum(temporal_eigenvalues)
par(mfrow=c(1,1))
# Plot the cumulative sum of the eigenvalues and the fraction of variance explained
lines(cumsum_eigvals/sum(temporal_eigenvalues), col = 'blue', lwd = 3)

abline(h = 0.95, col = 'red', lty = 2)
# 95%
n_components_PCA <- which(cumsum_eigvals/sum(temporal_eigenvalues) >= 0.95)[1]
n_components_PCA
n_components_FT <- which(cumulative_variance_explained_fourier >= 0.95)[1]
n_components_FT

abline(v = n_components_PCA, lty = 3) ; abline(v = n_components_FT, lty = 3)
text(x = n_components_PCA-100, y = 0.01, labels = paste0("n = ", n_components_PCA))
text(x = n_components_FT-100, y = 0.01, labels = paste0("n = ", n_components_FT))

legend("bottom", legend = c("PCA", "Fourier", "95%"), fill = c("blue", "black", 'red'), horiz = TRUE)
###############################################
# C
uccstudentno <- 119378796
set.seed(uccstudentno+2000)

rand_seq = rbinom(n = 128*128*1790, p = .5, size = 1)
U = as.matrix(rand_seq, nrow = 128*128, ncol = 't')

ZNA = Z; ZNA[U==0] = NA

# Count the number of missing values per row
num_missing <- apply(ZNA, 1, function(x) sum(is.na(x)))

# Create a frequency table of missing cases per row
freq_table <- table(num_missing)

# Print the frequency table and summary statistics
print(freq_table)
cat("\n")
cat("Median: ", median(num_missing), "\n") ;cat("Mean: ", mean(num_missing), "\n") ; cat("Maximum: ", max(num_missing), "\n")

par(mar = c(3,3,3,3))

# Create a histogram of missing cases per row
hist(num_missing, xlab = "Number of missing cases per row", main = "Distribution of missing cases")

par(mfrow = c(1,1))

covZNA = cov(ZNA, use = 'pairwise.complete.obs')

covZNA <- scale(covZNA)

egZNA = eigen(covZNA)

percent_var_explained_ZNA = cumsum(abs(egZNA$values))/sum(abs(egZNA$values))

plot(percent_var_explained_ZNA, type='l', ylim = c(0,1),
     xlab = 'Number of Components', ylab = 'Cumulative sum of Varaince Explained',
     col = 'red', main = 'Plot of Variance Explained by Z and ZNA', lwd = 2)
lines(cumsum_eigvals/sum(temporal_eigenvalues), type = 'l', col = 'blue', lwd = 2)
legend("bottomright", legend = c("Z", "ZNA"), fill = c("blue", 'red'), cex = 1.5, horiz = TRUE)
###############################################################################################################################################
# Q3

uccstudentno <- 119378796 ; set.seed(uccstudentno+2000)  # set seed for reproducibility

# A
theta <- sample(seq(0, 2*pi, length=27), replace=FALSE) ; b <- 0.2
x <- cos(theta) ; y <- b*sin(theta)
A <- cbind(x, y)
A <- scale(A, scale=FALSE)
par(mforw = c(1,1))
plot(A, xlab="", ylab="", xlim=c(-1, 1)*1.1, ylim=c(-1, 1)*0.5, cex=2.5, main='Input Data')
for(i in 1:26) {
  text(A[i, 1], A[i, 2], as.character(i))
}

# B
E <- dist(A)

#C
B <- cmdscale(E, k = 2)
x <- -B[,1] ; y <- -B[,2]
plot(x, y, xlab="", ylab="", xlim=c(-1, 1)*1.1, ylim=c(-1, 1)*0.5, cex=2.5, main="2D Representation from MDS") # Flips top to bottom instead of left to right in example
for(i in 1:26) {
  text(x[i], y[i], as.character(i))
}

# D
library(bigutilsr)
out <- procrustes(A, B)
R_hat <- out$R
B_hat <- B %*% R_hat

# E
plot(B_hat, xlab="", ylab="", xlim=c(-1, 1)*1.1, ylim=c(-1, 1)*0.5, main = "Map after Procrustes Transformation", cex = 2.5)
for(i in 1:26) {
  text(B_hat[i, 1], B_hat[i, 2], as.character(i))
}

# F
par(mfrow = c(2,2), mar = c(2,2,2,2))
boxplot(A[,1] - B[,1], main = "Column 1: A-B", col = "#32a852")
boxplot(A[,2] - B[,2], main = "Column 2: A-B", col = "#32a89d")
boxplot(A[,1] - B_hat[,1], main = "Column 1: A-B_hat", col = "#6d32a8")
boxplot(A[,2] - B_hat[,2], main = "Column 2: A-B_hat", col = "#f5be40")


plot(A, xlab="", ylab="", xlim=c(-1, 1)*1.1, ylim=c(-1, 1)*0.5, cex=1.5, main='Input Data')
for(i in 1:26) {
  text(A[i, 1], A[i, 2], as.character(i), cex=0.5)
}
plot(x, y, xlab="", ylab="", xlim=c(-1, 1)*1.1, ylim=c(-1, 1)*0.5, cex=1.5, main="2D Representation from MDS") # Flips top to bottom instead of left to right in example
for(i in 1:26) {
  text(x[i], y[i], as.character(i), cex=0.5)
}
plot(B_hat, xlab="", ylab="", xlim=c(-1, 1)*1.1, ylim=c(-1, 1)*0.5, main = "Map after Procrustes Transformation", cex = 1.5)
for(i in 1:26) {
  text(B_hat[i, 1], B_hat[i, 2], as.character(i), cex=0.5)
}

boxplot(A[,1] - B[,1],A[,2] - B[,2], A[,1] - B_hat[,1], A[,2] - B_hat[,2], col="#e35d6f", main="Raw and Final Differences")


##################################################################################################################################################
# Q4

load("C:/Users/Padraig/Desktop/DSA4_SEM2/Multivariate/CA2/Materials/Digits.Rdata")

X=t(matrix(c(digits),ncol=10*1000))
y= c(matrix(rep(rep(0:9),1000),ncol=10,byrow=T))

# Compute the TSS matrix
TSS <- t(X) %*% X

# Compute the eigendecomposition of the TSS matrix
eigen_TSS <- eigen(TSS)

eig_values <- eigen_TSS$values ; eig_vectors <- eigen_TSS$vectors

# Sort eigenvalues and eigenvectors in decreasing order of eigenvalue
order <- order(-eig_values)
eig_values <- eig_values[order] ; eig_vectors <- eig_vectors[, order]

# Compute the diagonal matrix of eigenvalues
eig_diag <- diag(eig_values)

# Compute the spectral decomposition
spectral_decomp <- eig_vectors %*% eig_diag %*% t(eig_vectors)

# Extract the eigenvectors of the TSS matrix
v <- eigen_TSS$values
Gamma <- eigen_TSS$vectors


# Transform the data matrix X using the eigenvectors Gamma
colmeans_ <- colMeans(X)
Z <- scale(X, scale = F, center = colmeans_) %*% Gamma

digit_means=NULL ; for(k in 0:9) { digit_means=rbind(digit_means,apply(Z[y==k,],2,mean))}


par(mfrow = c(2,2))
# Plot Z1 vs Z2
plot(Z[,1], Z[,2], xlab = "Z1", ylab = "Z2", main = "Z1 vs Z2")
for (i in 0:9){
  points(digit_means[i,1], digit_means[i,2], col = rainbow(10)[i], pch = as.character(i), cex = 1.25)
}

plot(Z[,1], Z[,3], xlab = "Z1", ylab = "Z3", main = "Z1 vs Z3")
for (i in 0:9){
  points(digit_means[i,1], digit_means[i,3], col = rainbow(10)[i], pch = as.character(i), cex = 1.25)
}

plot(Z[,2], Z[,3], xlab = "Z2", ylab = "Z3", main = "Z2 vs Z3")
for (i in 0:9){
  points(digit_means[i,2], digit_means[i,3], col = rainbow(10)[i], pch = as.character(i), cex = 1.25)
}

cov_matrix <- cov(X)
eigen_decomp <- eigen(cov_matrix)

total_var <- sum(eigen_decomp$values)
cumulative_var <- cumsum(eigen_decomp$values) / total_var

# Calculate the fraction of variance explained by each eigenvector
variance_fraction <- eigen_TSS$values / sum(eigen_TSS$values)

# Calculate the cumulative amount of variance explained
cumulative_variance <- cumsum(variance_fraction)

# Create the plot
plot(variance_fraction, type="o", xlab="Eigenvectors", ylab="Fraction of variance explained",
     main="variance explained", ylim = c(0,1), col = "blue", log = "x"); par(new = T)
plot(cumulative_var, type="o",col = "red", axes = F, xlab = "", ylab = "", log = "x"); axis(4)

legend("topleft", legend = c("Fraction", "Cumulative"),
       fill = c("blue", "red"), cex = 0.8)

# B
cluster_means=NULL

for(d in 0:9) {
  
  clust=hclust(dist(X[y==d,])); l=cutree(clust,k=5)
  
  for(i in 1:5) {cluster_means=cbind(cluster_means,apply((X[y==d,])[l==i,],2,mean)) }
  
  if(d==2|d==9) {
    
    par(mfrow=c(2,3))
    
    plot(clust,main=paste0("Digit: ",as.character(d)),cex=.3)
    
    for(i in 1:5) { num_images=matrix(cluster_means[,d*5+i],ncol=28) ; image(num_images) }
  }}

####################################################################################
# Q5 
K = 50
Bp = Gamma[,1:K]
ids = rep(0:9, each = 1000)
oo.p = summary(lm(t(X)~Bp))
Bc = cluster_means
oo.c = summary(lm(t(X)~Bc))

rsq_p <- numeric(50) ; rsq_c <- numeric(50)
rsq_p_k <- numeric(50) ; rsq_c_k <- numeric(50)

for (i in 1:10000) {
  response_col <- paste0("Response Y", i)
  rsq_p[i] <- oo.p[[response_col]]$r.squared
  rsq_c[i] <- oo.c[[response_col]]$r.squared
}

for (i in 1:50) {
  response_col <- paste0("Response Y", i)
  rsq_p_k[i] <- oo.p[[response_col]]$r.squared
  rsq_c_k[i] <- oo.c[[response_col]]$r.squared
}


# print results
cat("Proportion of variance explained by Bp model:", mean(rsq_p), "\n") ; cat("Proportion of variance explained by Bc model:", mean(rsq_c), "\n")
cat("Proportion of variance explained by Bp K-dimensional model:", mean(rsq_p_k), "\n") ; cat("Proportion of variance explained by Bc K-dimensional model:", mean(rsq_c_k), "\n")

prop_var_expl <- data.frame(Bp = rsq_p,
                             Bc = rsq_c,
                             Bp_k = rsq_p_k,
                             Bc_k = rsq_c_k)

prop_var_expl
####################################################################################################
library(MASS)

Xc <- as.matrix(X %*% Bc) ; Xp <- as.matrix(X %*% Bp)

y <- factor(rep(0:9, each = 1000))
fit_lda_c <- lda(Xc, y, cv = TRUE) ; fit_lda_p <- lda(Xp, y, cv = TRUE)

fit_qda_c <- qda(Xc, y, cv = TRUE) ; fit_qda_p <- qda(Xp, y, cv = TRUE)

predict_lda_c <- predict(fit_lda_c, Xc) ; predict_lda_p <- predict(fit_lda_p, Xp)

predict_qda_c <- predict(fit_qda_c, Xc) ; predict_qda_p <- predict(fit_qda_p, Xp)

misclass_Calc <- function(table) {
  n <- ncol(table)
  rates <- c()
  for (i in 1:n) {
    col = table[,i]
    correct = col[i]
    sum = sum(col)
    misclass = sum - correct
    rates <- c(rates, misclass/sum)
  }
  return(rates)
}


lda_c_pred <- table(predict_lda_c$class, y)
lda_c <-misclass_Calc(lda_c_pred) ; lda_c_norm <- apply(lda_c_pred, MARGIN = 2, function(x) x/sum(x))
lda_c_norm

lda_p_pred <- table(predict_lda_p$class, y)
lda_p_misclass <- misclass_Calc(lda_p_pred) ; lda_p_norm <- apply(lda_p_pred, MARGIN = 2, function(x) x/sum(x))
lda_p_norm

qda_c_pred <- table(predict_qda_c$class, y)
qda_c_misclass <- misclass_Calc(qda_c_pred) ; qda_c_norm <- apply(qda_c_pred, MARGIN = 2, function(x) x/sum(x))
qda_c_norm

qda_p_pred <- table(predict_qda_p$class, y)
qda_p_misclass <- misclass_Calc(qda_p_pred) ; qda_p_norm <- apply(qda_p_pred, MARGIN = 2, function(x) x/sum(x))
qda_p_norm

misclass_table <- data.frame(Type = unique(y),
                             LDA_C = lda_c,
                             LDA_P = lda_p_misclass,
                             QDA_C = qda_c_misclass,
                             QDA_P = qda_p_misclass)
                             

rownames(misclass_table) <- NULL ; misclass_table
library(knitr)
kable(misclass_table, digits = 3, align = "c", caption = "Misclassification Rates by Type")
kable(lda_c_norm, digits = 3, align = "c", caption = "Misclassification Rates by Type")
kable(lda_p_norm, digits = 3, align = "c", caption = "Misclassification Rates by Type")
kable(qda_c_norm, digits = 3, align = "c", caption = "Misclassification Rates by Type")
kable(qda_p_norm, digits = 3, align = "c", caption = "Misclassification Rates by Type")

lda_acc_c <- mean(predict_lda_c$class == y) ; lda_acc_p <- mean(predict_lda_p$class == y)

cat("Linear Discriminant Analysis C Misclass:", round((1-lda_acc_c)*100, 2), "%\n")
cat("Linear Discriminant Analysis P Misclass:", round((1-lda_acc_p)*100, 2), "%\n")
qda_acc_c <- mean(predict_qda_c$class == y) ; qda_acc_p <- mean(predict_qda_p$class == y)

cat("Quadratic Discriminant Analysis C Misclass:", round((1-qda_acc_c)*100, 2), "%\n")
cat("Quadratic Discriminant Analysis P Misclass:", round((1-qda_acc_p)*100, 2), "%\n")

####################################################################################
# Q6 
# fix tables

data <- read.table("C:/Users/Padraig/Downloads/table2.dat")

colnames(data) <- c("ID", "Seq", "RAdeg", "DEdeg", 'Per', "R21", "phi21", "T0", "gmag", 'rmag', 'Per-g', 'Per-r', 'Ng', 'Nr', 'R21-g', 'R21-r', 'phi21-g',
                    'phi21-r', 'R2-g', 'R2-r', 'gAmp', 'rAmp', 'log(FAP-g)', 'log(FAP-r)', 'Type', 'Dmin-g','Dmin-r')

x <- data[3:23]
oo= lda(x=x, grouping=data[,"Type"])
# extract the loading vectors from the lda object
loadings <- oo$scaling

# get the indices of the best 4 discriminant variables
best_vars <- order(oo$svd, decreasing = TRUE)[1:4]

# subset the loading vectors for the best variables
loadings_best <- loadings[, best_vars]

par(mfrow = c(4,3), mar = c(2,2,2,2))

# plot the loading vectors
matplot(loadings_best, type = "l", xlab = "Predictor variables", ylab = "Loading vectors", main = "Loading Vectors")
#legend("bottom", legend = colnames(loadings_best), col = 1:4, lty = 1, horiz = T)


lda_predict <- predict(oo, x) ; discriminant_scores <- lda_predict$x

f_vals <- oo$svd^2
# print the F-test values
print(f_vals) ; plot(1:10, f_vals, main = "F-Tests", type = "b", pch = 19)

for (i in 1:10) {
  boxplot(discriminant_scores[,i] ~ data[,"Type"], main=paste("Discriminant Variable", i), ylab="Score", outline = FALSE , xlab = "", col = rainbow(11))
}

# B

y <- as.factor(data["Type"])

lda_model_resub= lda(x=x, grouping=data[,"Type"]) ; lda_model_cv = lda(x=x, grouping=data[,"Type"], cv = TRUE)

lda_predict_resub <- predict(lda_model_resub, x) ; lda_predict_cv <- predict(lda_model_cv, x)

qda_model_resub= qda(x=x, grouping=data[,"Type"]) ; qda_model_cv = qda(x=x, grouping=data[,"Type"], cv = TRUE)

qda_predict_resub <- predict(qda_model_resub, x) ; qda_predict_cv <- predict(qda_model_cv, x)

lda_resub_table <- table(lda_predict_resub$class, data[,"Type"]) ; lda_cv_table <- table(lda_predict_cv$class, data[,"Type"])

qda_resub_table <- table(qda_predict_resub$class, data[,"Type"]) ; qda_cv_table <- table(qda_predict_cv$class, data[,"Type"])

lda_resub <- apply(lda_resub_table, MARGIN = 2, function(x) x/sum(x)) ; lda_cv <- apply(lda_cv_table, MARGIN = 2, function(x) x/sum(x))
qda_resub <- apply(qda_resub_table, MARGIN = 2, function(x) x/sum(x)) ; qda_cv <- apply(qda_cv_table, MARGIN = 2, function(x) x/sum(x))

kable(lda_resub, digits = 3, align = "c", caption = "LDA Resub")
kable(lda_cv, digits = 3, align = "c", caption = "LDA CV")
kable(qda_resub, digits = 3, align = "c", caption = "QDA Resub")
kable(qda_cv, digits = 3, align = "c", caption = "QDA CV")

misclass_table <- data.frame(Type = unique(data[,"Type"]),
                             LDA_Resub = 1-diag(lda_resub),
                             LDA_CV = 1-diag(lda_cv),
                             QDA_Resub = 1-diag(qda_resub),
                             QDA_CV = 1-diag(qda_cv))
rownames(misclass_table) <- NULL ; misclass_table
library(knitr)
kable(misclass_table, digits = 3, align = "c", caption = "Misclassification Rates by Type")

cat("Linear Discriminant Analysis Resub Misclass:", round(mean(1-diag(lda_resub))*100, 2), "%\n")
cat("Linear Discriminant Analysis Cross Validation Misclass:", round(mean(1-diag(lda_cv))*100, 2), "%\n")
cat("Quadratic Discriminant Analysis Resub Misclass:", round(mean(1-diag(qda_resub))*100, 2), "%\n")
cat("Quadratic Discriminant Analysis Cross Validation Misclass:", round(mean(1-diag(qda_cv))*100, 2), "%\n")

# C

qda_scores <- qda_predict_resub$posterior

lda_qda_post <- lda(qda_scores, grouping = data[,"Type"])

loadings_lqda <- lda_qda_post$scaling[,1:4]


par(mfrow = c(4,3), mar = c(2,2,2,2))

# plot the loading vectors
matplot(loadings_lqda, type = "l", xlab = "Predictor variables", ylab = "Loading vectors", main = "Loading Vectors")

f_vals <- lda_qda_post$svd^2
lqda_discrim_scores <- predict(lda_qda_post, qda_scores)$x

lqda_discrim_scores <- as.matrix(qda_scores) %*% lda_qda_post$scaling

plot(1:10, f_vals, main = "F-Tests", type = 'b', pch = 19)
for (i in 1:10) {
  boxplot(lqda_discrim_scores[,i] ~ data[,"Type"], main=paste("Discriminant Variable: ", i), ylab="Class Probability", outline = FALSE, col = rainbow(11))
}
