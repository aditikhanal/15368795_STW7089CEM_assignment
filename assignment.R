
#List of Libraries used
install.packages('matlib')#library for Linear Algebra and Multivariate Statistics
library(matlib)

install.packages("ggplot2") #Library for creating time series plot
library(ggplot2)

install.packages("rsample")
library(rsample)

gene_data = (read.csv(file = "E:/gene_data.csv", header = T))

#converting the gene data into matrix
M = as.matrix(gene_data)
M

colnames(M) = c("T", "X1", "X2", "X3","X4","X5")

#extracting Time from matrix
T <- M[, "T"]
T

#extracting gene from matrix
Xs<- M[, c("X1", "X2", "X3", "X4", "X5")]
Xs


XSignals = ts(Xs)#ts function create time series object
XSignals


# Task 1: • Time series plots (of each gene against sampling time) 
plot_time_series <- function(matrix, time_col, gene_col) {
  plot(matrix[, time_col], matrix[, gene_col], type = "l", 
       xlab = "Time (minutes)", ylab = "Expression Level", 

       main = paste("Time Series of", colnames(matrix)[gene_col]))
}

# Plot time series for each gene
par(mfrow = c(1, 1))  # Arrange plots in a 3x2 grid
for (i in 2:ncol(M)) {
  plot_time_series(M, 1, i)
}

#All Input Signals
# Extract the gene expression columns
Xs <- M[, c("X1", "X2", "X3", "X4", "X5")]

# Set up plotting layout for histograms and density plots
par(mfrow = c(3,2)) # Arrange plots in a 3x2 grid


# Function to create histogram and density plot
plot_hist_density <- function(data, col_name) {
  hist(data, freq = FALSE, col = "blue", xlab = col_name, main = paste("Histogram and density curve of", col_name, "Signal"))
  lines(density(data), lwd = 2, col = "red")
}

# Plot histograms and density curves for each gene
plot_hist_density(Xs[, "X1"], "X1")
plot_hist_density(Xs[, "X2"], "X2")
plot_hist_density(Xs[, "X3"], "X3")
plot_hist_density(Xs[, "X4"], "X4")
plot_hist_density(Xs[, "X5"], "X5")



# Plotting correlations and scatter plots
par(mfrow = c(3,2))# Arrange plots in a 3x2 grid

# Define a function to create scatter plots with regression lines
plot_correlation <- function(X1, X2, name1, name2) {
  plot(X1, X2, pch = 1, col = "blue", 
       main = paste("Correlation between", name1, "and", name2),
       xlab = name1, ylab = name2)
  abline(lm(X2 ~ X1), col = "red", lwd = 1)
}

# Plotting pairs of genes
plot_correlation(Xs[, "X1"], Xs[, "X2"], "X1", "X2")
plot_correlation(Xs[, "X1"], Xs[, "X3"], "X1", "X3")
plot_correlation(Xs[, "X1"], Xs[, "X4"], "X1", "X4")
plot_correlation(Xs[, "X1"], Xs[, "X5"], "X1", "X5")
plot_correlation(Xs[, "X2"], Xs[, "X3"], "X2", "X3")
plot_correlation(Xs[, "X2"], Xs[, "X4"], "X2", "X4")
plot_correlation(Xs[, "X2"], Xs[, "X5"], "X2", "X5")
plot_correlation(Xs[, "X3"], Xs[, "X4"], "X3", "X4")
plot_correlation(Xs[, "X3"], Xs[, "X5"], "X3", "X5")
plot_correlation(Xs[, "X4"], Xs[, "X5"], "X4", "X5")



#--------------------------------------------------------------------------------
# Task 2: Regression – modeling the relationship between gene expression signals
#--------------------------------------------------------------------------------


# Creating a matrix of ones
onesMatrix = matrix(1, nrow(Xs), 1) # Creating matrix of ones
#Since X2 is the output gene
Y <- Xs[, "X2"]


# Creating Xata for model 1
XData_model1 = cbind(onesMatrix, Xs[, "X4"], Xs[, "X3"]^2)
# Creating Xdata for model 2
XData_model2 = cbind(onesMatrix, Xs[, "X4"], Xs[, "X3"]^2, Xs[, "X5"])
# Creating XData for model 3
XData_model3 = cbind(Xs[, "X3"], Xs[, "X4"], Xs[, "X5"]^3)
# Creating Data for model 4
XData_model4 = cbind(onesMatrix, Xs[, "X4"], Xs[, "X3"]^2, Xs[, "X5"]^3)
# Creating Data for model 5 from
XData_model5 = cbind(onesMatrix, Xs[, "X4"], Xs[, "X1"]^2, Xs[, "X3"]^2)

#--------------------------------------------------------------------------------
# Task 2.1: 
#--------------------------------------------------------------------------------

# Calculation of the least square (thetaHat) for each model

# Model 1
model1_Thetahat = solve(t(XData_model1) %*% XData_model1) %*% t(XData_model1) %*% Y
model1_Thetahat

# Model 2
model2_Thetahat = solve(t(XData_model2) %*% XData_model2) %*% t(XData_model2) %*% Y
model2_Thetahat

# Model 3
model3_Thetahat = solve(t(XData_model3) %*% XData_model3) %*% t(XData_model3) %*% Y
model3_Thetahat

# Model 4
model4_Thetahat = solve(t(XData_model4) %*% XData_model4) %*% t(XData_model4) %*% Y
model4_Thetahat

# Model 5
model5_Thetahat = solve(t(XData_model5) %*% XData_model5) %*% t(XData_model5) %*% Y
model5_Thetahat

#-------------------------------------------------------------------------------- 
#Task 2.2: Model Residual Error
#--------------------------------------------------------------------------------
# Calculating y-hat
# Model 1
model1_YHat = XData_model1 %*% model1_Thetahat
model1_YHat

# Model 2
model2_YHat = XData_model2 %*% model2_Thetahat
model2_YHat

# Model 3
model3_YHat = XData_model3 %*% model3_Thetahat
model3_YHat

# Model 4
model4_YHat = XData_model4 %*% model4_Thetahat
model4_YHat

# Model 5
model5_YHat = XData_model5 %*% model5_Thetahat
model5_YHat


# Calculating RSS
# Model 1
RSS_model1 = sum((Y-model1_YHat)^2)
RSS_model1

# Model 2
RSS_model2 = sum((Y-model2_YHat)^2)
RSS_model2

# Model 3
RSS_model3 = sum((Y-model3_YHat)^2)
RSS_model3

# Model 4
RSS_model4 = sum((Y-model4_YHat)^2)
RSS_model4

# Model 5
RSS_model5 = sum((Y-model5_YHat)^2)
RSS_model5



#-------------------------------------------------------------------------------- 
#Task 2.3: Calculating Likelihood and Variance for all models
#--------------------------------------------------------------------------------
n = length(Y) #Calculating length of Y
# Calculating Variance for
# Model 1
VAR_model1 = RSS_model1/(n-1)
VAR_model1

# Model 2
VAR_model2 = RSS_model2/(n-1)
VAR_model2

# Model 3
VAR_model3 = RSS_model3/(n-1)
VAR_model3

# Model 4
VAR_model4 = RSS_model4/(n-1)
VAR_model4

# Model 5
VAR_model5 = RSS_model5/(n-1)
VAR_model5

# Calculating likelihood for
# Model 1
Likelihood_model1 = -(n/2)*(log(2*pi))-(n/2)*(log(VAR_model1))-(1/(2*VAR_model1))*RSS_model1
Likelihood_model1

# Model 2
Likelihood_model2 = -(n/2)*(log(2*pi))-(n/2)*(log(VAR_model2))-(1/(2*VAR_model2))*RSS_model2
Likelihood_model2

# Model 3
Likelihood_model3 = -(n/2)*(log(2*pi))-(n/2)*(log(VAR_model3))-(1/(2*VAR_model3))*RSS_model3
Likelihood_model3

# Model 4
Likelihood_model4 = -(n/2)*(log(2*pi))-(n/2)*(log(VAR_model4))-(1/(2*VAR_model4))*RSS_model4
Likelihood_model4

# Model 5
Likelihood_model5 = -(n/2)*(log(2*pi))-(n/2)*(log(VAR_model5))-(1/(2*VAR_model5))*RSS_model5
Likelihood_model5


#-------------------------------------------------------------------------------- 
#Task 2.4: Calculating AIC and BIC for all models
#--------------------------------------------------------------------------------
# Calculating AIC
# Model 1
AIC_model1 = 2*(length(model1_Thetahat))-2*Likelihood_model1
AIC_model1

# Model 2
AIC_model2 = 2*(length(model2_Thetahat))-2*Likelihood_model2
AIC_model2

# Model 3
AIC_model3 = 2*(length(model3_Thetahat))-2*Likelihood_model3
AIC_model3

# Model 4
AIC_model4 = 2*(length(model4_Thetahat))-2*Likelihood_model4
AIC_model4

# Model 5
AIC_model5 = 2*(length(model1_Thetahat))-2*Likelihood_model5
AIC_model5


# Calculating BIC
# Model 1
BIC_model1 = length(model1_Thetahat)*log(n)-2*Likelihood_model1
BIC_model1

# Model 2
BIC_model2 = length(model2_Thetahat)*log(n)-2*Likelihood_model2
BIC_model2

# Model 3
BIC_model3 = length(model3_Thetahat)*log(n)-2*Likelihood_model3
BIC_model3

# Model 4
BIC_model4 = length(model4_Thetahat)*log(n)-2*Likelihood_model4
BIC_model4

# Model 5
BIC_model5 = length(model5_Thetahat)*log(n)-2*Likelihood_model5
BIC_model5

BIC_model1
BIC_model2
BIC_model3
BIC_model4
BIC_model5
#-------------------------------------------------------------------------------- 
#Task 2.5: Calculating Error for all models and Plotting Q-Q plot with Q-Q line for them
#--------------------------------------------------------------------------------

par(mfrow = c(3, 2))
# Model 1
Error_model1 = Y - model1_YHat #Error
qqnorm(Error_model1, col = "#336600", main = "Q-Q plot of Model 1") # Plots Graph
qqline(Error_model1, col = "#e60000", lwd = 1) # Adds Q-Q line on graph

# Model 2
Error_model2 = Y - model2_YHat #Error
qqnorm(Error_model2, col = "#336600", main = "Q-Q plot of Model 2") # Plots Graph
qqline(Error_model2, col = "#e60000", lwd = 1) # Adds Q-Q line on graph

# Model 3
Error_model3 = Y - model3_YHat #Error
qqnorm(Error_model3, col = "#336600", main = "Q-Q plot of Model 3") # Plots Graph
qqline(Error_model3, col = "#e60000", lwd = 1) # Adds Q-Q line on graph

# Model 4
Error_model4 = Y - model4_YHat #Error
qqnorm(Error_model4, col = "#336600", main = "Q-Q plot of Model 4") # Plots Graph
qqline(Error_model4, col = "#e60000", lwd = 1) # Adds Q-Q line on graph

# Model 5
Error_model5 = Y - model5_YHat #Error
qqnorm(Error_model5, col = "#336600", main = "Q-Q plot of Model 5") # Plots Graph
qqline(Error_model5, col = "#e60000", lwd = 1) # Adds Q-Q line on graph

#--------------------------------------------------------------------------------
#-------------------------------------------------------------------------------- 
#Task 2.6: Selecting a model
#--------------------------------------------------------------------------------
# Calculating Mean Squared Error (MSE)
# Model 1
MSE_model1 = sum(Error_model1^2)/length(Error_model1)
MSE_model1

# Model 2
MSE_model2 = sum(Error_model2^2)/length(Error_model2)
MSE_model2

# Model 3
MSE_model3 = sum(Error_model3^2)/length(Error_model3)
MSE_model3

# Model 4
MSE_model4 = sum(Error_model4^2)/length(Error_model4)
MSE_model4

# Model 5
MSE_model5 = sum(Error_model5^2)/length(Error_model5)
MSE_model5


