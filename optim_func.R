rm(list = ls())


# Example Data
set.seed(123)  # For reproducibility
n <- 100
p <- 5
x <- matrix(rnorm(n * p), nrow = n, ncol = p)
y <- x %*% c(1, 2, -1, 0.5, 0) + rnorm(n)

# Elastic Net Objective Function (Modified for LOOCV)
elastic_net_objective_loocv <- function(beta, x, y, lambda1, lambda2, leave_out_index) {
  # Leave out the specified observation
  x_train <- x[-leave_out_index, ]
  y_train <- y[-leave_out_index]
  
  # Calculate predicted values on the training set
  y_pred <- x_train %*% beta
  
  # Calculate squared error loss on the training set
  loss <- sum((y_train - y_pred)^2)
  
  # Calculate L1 and L2 penalties
  l1_penalty <- lambda1 * sum(abs(beta))
  l2_penalty <- lambda2 * sum(beta^2)
  
  # Return the objective function value
  return(loss + l1_penalty + l2_penalty)
}

# Function to Perform LOOCV for a Given Lambda
loocv_elastic_net <- function(x, y, lambda1, lambda2) {
  n <- nrow(x)
  cv_errors <- numeric(n)  # Store cross-validation errors
  
  for (i in 1:n) {
    result <- optim(par = rep(0, ncol(x)),
                    fn = elastic_net_objective_loocv,
                    x = x, 
                    y = y, 
                    lambda1 = lambda1,
                    lambda2 = lambda2,
                    leave_out_index = i,
                    method = "BFGS")
    
    # Prediction on the left-out observation
    y_pred <- x[i, ] %*% result$par
    
    # Calculate squared error for this fold
    cv_errors[i] <- (y[i] - y_pred)^2
  }
  
  # Mean squared error across all folds
  return(mean(cv_errors))
}

prostate_data = read.table("data/prostate_data.txt")
head(prostate_data)
in_sample_obs = prostate_data%>%as_tibble() %>%
  filter(train==TRUE)
out_of_sample_obs = prostate_data%>%as_tibble() %>%
  filter(train==FALSE)

prostate_data = read.table("data/prostate_data.txt")
head(prostate_data)
in_sample_obs = prostate_data%>%as_tibble() %>%
  filter(train==TRUE)
out_of_sample_obs = prostate_data%>%as_tibble() %>%
  filter(train==FALSE)

in_sample_obs_std = scale(in_sample_obs[,-10], center=TRUE, scale=TRUE)%>%as_tibble()
Y_train = as.matrix(in_sample_obs_std$lpsa)
colnames(Y_train) = "lpsa"
X_train = as.matrix(in_sample_obs_std[, -9])

# LOOCV for Ridge (lambda1 = 0)
lambda_values <- seq(0, 5, length.out = 100)  # Range of lambda2 values
cv_errors_ridge <- sapply(lambda_values, function(lambda) loocv_elastic_net(X_train, Y_train, 0, lambda))

# LOOCV for Lasso (lambda2 = 0)
cv_errors_lasso <- sapply(lambda_values, function(lambda) loocv_elastic_net(X_train, Y_train, lambda, 0))

# Find the optimal lambda values
optimal_lambda_ridge <- lambda_values[which.min(cv_errors_ridge)]
optimal_lambda_lasso <- lambda_values[which.min(cv_errors_lasso)]

cat("Optimal Lambda (Ridge):", optimal_lambda_ridge, "\n")
cat("Optimal Lambda (Lasso):", optimal_lambda_lasso, "\n")

# Plotting the LOOCV Results
library(ggplot2)  # For plotting

# Create a data frame for plotting
df <- data.frame(
  lambda = rep(lambda_values, 2),
  mse = c(cv_errors_ridge, cv_errors_lasso),
  model = rep(c("Ridge", "Lasso"), each = length(lambda_values))
)

ggplot(df[df$model == "Ridge", ], aes(x = lambda, y = mse)) +
  geom_line() +
  geom_point() +
  labs(x = "Lambda", y = "Mean Squared Error (MSE)", title = "LOOCV Results: Ridge Regression") +
  theme_minimal()


ggplot(df[df$model == "Lasso", ], aes(x = lambda, y = mse)) +
  geom_line() +
  geom_point() +
  labs(x = "Lambda", y = "Mean Squared Error (MSE)", title = "LOOCV Results: Lasso Regression LOOCV") +
  theme_minimal()

