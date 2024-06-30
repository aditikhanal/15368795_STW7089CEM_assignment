
#Libraries used
install.packages('matlib')#library for Linear Algebra and Multivariate Statistics
library(matlib)

install.packages("ggplot2") #Library for creating time series plot
library(ggplot2)

install.packages("rsample")
library(rsample)


#Reading the csv file into gene_data
gene_data = (read.csv(file = "E:/gene_data.csv", header = TRUE))

#converting the gene data into matrix
M = as.matrix(gene_data)
M

#Column names for Time series and gene expression
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

# To extract the gene expression columns
Xs <- M[, c("X1", "X2", "X3", "X4", "X5")]

#layout for histograms and density plots
par(mfrow = c(3,2)) # Arrange plots in a 3x2 grid


# This function creates histogram and density plot
plot_hist_density <- function(data, col_name) {
  hist(data, freq = FALSE, col = "blue", xlab = col_name, 
       main = paste("Histogram and density curve of", col_name, "Signal"))
  lines(density(data), lwd = 2, col = "red")
}

# Plotting histograms and density plots for each gene
plot_hist_density(Xs[, "X1"], "X1")
plot_hist_density(Xs[, "X2"], "X2")
plot_hist_density(Xs[, "X3"], "X3")
plot_hist_density(Xs[, "X4"], "X4")
plot_hist_density(Xs[, "X5"], "X5")



# To plot correlation and scatter plots
par(mfrow = c(3,2))# Arrange plots in a 3x2 grid

# Defining function for creation of scatter plots with regression lines
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




# Task 2: Regression – modeling the relationship between gene expression signals


# Creating a matrix of ones
onesMatrix = matrix(1, nrow(Xs), 1) # Creating matrix of ones
#Since X2 is the output gene
Y <- M[, c("X2")]
Y
#input gene expressions
Xin<- M[, c("X1","X3", "X4", "X5")]
Xin


# To create Xdata for model 1
XData_model1 = cbind(onesMatrix, Xs[, "X4"], Xs[, "X3"]^2)
# To create Xdata for model 2
XData_model2 = cbind(onesMatrix, Xs[, "X4"], Xs[, "X3"]^2, Xs[, "X5"])
# To create XData for model 3
XData_model3 = cbind(Xs[, "X3"], Xs[, "X4"], Xs[, "X5"]^3)
# To create Data for model 4
XData_model4 = cbind(onesMatrix, Xs[, "X4"], Xs[, "X3"]^2, Xs[, "X5"]^3)
# To create Data for model 5 from
XData_model5 = cbind(onesMatrix, Xs[, "X4"], Xs[, "X1"]^2, Xs[, "X3"]^2)


# Task 2.1: 
# To calculate least square (thetaHat) for each model

# thetaHat for Model 1
model1_Thetahat = solve(t(XData_model1) %*% XData_model1) %*% t(XData_model1) %*% Y
model1_Thetahat

# thetaHat for Model 2
model2_Thetahat = solve(t(XData_model2) %*% XData_model2) %*% t(XData_model2) %*% Y
model2_Thetahat

# thetaHat for Model 3
model3_Thetahat = solve(t(XData_model3) %*% XData_model3) %*% t(XData_model3) %*% Y
model3_Thetahat

# thetaHat for Model 4
model4_Thetahat = solve(t(XData_model4) %*% XData_model4) %*% t(XData_model4) %*% Y
model4_Thetahat

# thetaHat for Model 5
model5_Thetahat = solve(t(XData_model5) %*% XData_model5) %*% t(XData_model5) %*% Y
model5_Thetahat


#Task 2.2: Model Residual Error

# Calculating y-hat for model 1
model1_YHat = XData_model1 %*% model1_Thetahat
model1_YHat

# Calculating y-hat for model 2
model2_YHat = XData_model2 %*% model2_Thetahat
model2_YHat

# Calculating y-hat for model 3
model3_YHat = XData_model3 %*% model3_Thetahat
model3_YHat

# Calculating y-hat for model 4
model4_YHat = XData_model4 %*% model4_Thetahat
model4_YHat

# Calculating y-hat for model 5
model5_YHat = XData_model5 %*% model5_Thetahat
model5_YHat


# Calculating RSS for model 1
RSS_model1 = sum((Y-model1_YHat)^2)
RSS_model1

# Calculating RSS for Model 2
RSS_model2 = sum((Y-model2_YHat)^2)
RSS_model2

# Calculating RSS for Model 3
RSS_model3 = sum((Y-model3_YHat)^2)
RSS_model3

# Calculating RSS for Model 4
RSS_model4 = sum((Y-model4_YHat)^2)
RSS_model4

#Calculating RSS for Model 5
RSS_model5 = sum((Y-model5_YHat)^2)
RSS_model5


#Task 2.3: Calculating Likelihood and Variance for all models
n = length(Y) #Calculating length of Y(output gene X2)

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



#Task 2.4: Calculating AIC and BIC for all models

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

#Task 2.5: Calculating Error for all models and Plotting Q-Q plot with Q-Q line for them


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

 
#Task 2.6: To select a model

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


#Task 2.7: To split the data sets into training and test set

# Splitting the data (Training Data set)
XSplit = initial_split(data = as.data.frame(Xin), prop = .7)
YSplit = initial_split(data = as.data.frame(Y), prop = .7)
# Training data For output gene Y=X2
Y_Training_Set = training(YSplit) # Y Training data set
Y_Training_Data = as.matrix(Y_Training_Set) # Y Training data
# Training data For input gene expressions
X_Training_Set = training(XSplit) # X Training data set
X_Training_Data = as.matrix(X_Training_Set) # X Training data

# Test Data for output gene Y=X2
Y_Testing_Set = testing(YSplit) # Y Testing data set
Y_Testing_Data = as.matrix(Y_Testing_Set) # Output Testing data
# Test data for input gene expressions
X_Testing_Set = testing(XSplit) # X Testing data set
X_Testing_Data = as.matrix(X_Testing_Set) # Input Testing data

# Selecting Model 5 and estimating model parameters
TrainingOneXMatrix = matrix(1, length(X_Training_Set$X1), 1) # ones matrix for training set
TrainingXModel = cbind(TrainingOneXMatrix, X_Training_Set[, "X4"],
                       (X_Training_Set[, "X1"])^2, (X_Training_Set[, "X3"])^2) #To Train Model5
TrainingThetaHat = solve(t(TrainingXModel) %*% TrainingXModel) %*% t(TrainingXModel) %*% Y_Training_Data
TrainingThetaHat

# Computing output/prediction of Model5 using testing data set
X_Testing_Data
TestingYHat = X_Testing_Data %*% TrainingThetaHat
RSStesting = sum((Y_Testing_Set - TestingYHat)^2)
RSStesting


t.test(Y_Training_Data, mu = 1, alternative = "two.sided", conf.level = 0.95)

C_I1 = 1.062184
C_I2 = 1.173627
meu = 1.117906

# With 95% of confidence interval, predicting the model and plotting them with testing data and error bars
par(mfrow = c(1, 1))
TrainingDensity = density(Y_Training_Data) # Density of training data of output signal
TrainingDensity
plot(TrainingDensity, col="#336600", lwd = 2, main="Distribution of Output Signal Training Data")
abline(v = C_I1, col = "#e60000", lty=2)
abline(v = C_I2, col = "#e60000", lty=2)
abline(v = meu, col = "#1a1a1a", lty=2)

residual = ((Y_Testing_Set - TestingYHat)) # Calculating Error
residual

# plotting Error Bars
# Calculating Standard Deviation (Sigma)
Sigma = sqrt(VAR_model5) #Variance of model5 from task 2.3
Sigma
XData_model5 #Data model 5 from task 2.1

dataFrame = data.frame(
  xAxis = XData_model5,
  yAxis = Y
)

dataFrame


par(mfrow = c(2,2))

ggplot(dataFrame) +
  geom_bar( aes(x=xAxis.1, y=Y), stat="identity", fill="#336600", alpha=0.7) +
  geom_errorbar( aes(x=xAxis.1, ymin=Y-Sigma, ymax=Y+Sigma), width=0.4, colour="#e60000", alpha=0.9, linewidth=1) +
  labs(title="Error Bar (Model 5 - X1)", x="Model 5 - X1", y = "Output Signal Data")

ggplot(dataFrame) +
  geom_bar( aes(x=xAxis.2, y=Y), stat="identity", fill="#336600", alpha=0.7) +
  geom_errorbar( aes(x=xAxis.2, ymin=Y-Sigma, ymax=Y+Sigma), width=0.4, colour="#e60000", alpha=0.9, linewidth=1) +
  labs(title="Error Bar (Model 5 - X2)", x="Model 5 - X2", y = "Output Signal Data")


ggplot(dataFrame) +
  geom_bar( aes(x=xAxis.3, y=Y), stat="identity", fill="#336600", alpha=0.7) +
  geom_errorbar( aes(x=xAxis.3, ymin=Y-Sigma, ymax=Y+Sigma), width=0.4, colour="#e60000", alpha=0.9, linewidth=1) +
  labs(title="Error Bar (Model 5 - X3)", x="Model 5 - X3", y = "Output Signal Data")

ggplot(dataFrame) +
  geom_bar( aes(x=xAxis.4, y=Y), stat="identity", fill="#336600", alpha=0.7) +
  geom_errorbar( aes(x=xAxis.4, ymin=Y-Sigma, ymax=Y+Sigma), width=0.4, colour="#e60000", alpha=0.9, linewidth=1) +
  labs(title="Error Bar (Model 5 - X4)", x="Model 5 - X4", y = "Output Signal Data")


#Task 3: Computing the Approximate Bayesian Computation (ABC)

array1 = 0
array2 = 0
f_value = 0
s_value = 0

# Model 5 thetahat values from
ThetaBias = 1.2951518 # Since this is the highest value we set it as chosen parameter
ThetaA = 0.8312983 # Since this is the highest value we set it as chosen parameter
ThetaB = 0.5385828 # The value for ThetaB is set constant
ThetaC = 0.1096679 # The value for ThetaC is set constant
Epsilon = RSS_model5 * 2 ## fixing value of epsilon, RSS_model5 from task 2.2
num = 100 # number of iteration
##Calculating Y-hat for performing rejection ABC
counter <- 0
for (i in 1:num) {
  range1 = runif(1, -1.295158, 1.2951518) # calculating the range
  range1
  range2 = runif(1, -0.8312983, 0.8312983)
  range2
  NewThetahat = matrix(c(range1, range2, ThetaB, ThetaC))
  NewYHat = XData_model5 %*% NewThetahat ## New Y hat and model5_Thetahat from task2.1
  NewRSS = sum((Y - NewYHat)^2)
  NewRSS
  if (NewRSS > Epsilon){ #Performing rejection ABC
    array1[i] = range1
    array2[i] = range2
    counter = counter + 1
    Fvalue = matrix(array1)
    Svalue = matrix(array2)
  }
}
# Plotting the graph for joint and marginal posterior distribution of model 5
plot(Fvalue, Svalue, col = c("#ff1a1a", "#3366ff"), main = "Joint and Marginal Posterior Distribution Model 5")
