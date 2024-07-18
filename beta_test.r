R-code for the used $L^2$-type statistical test


## Data: 

x = # data to use

## Theoretical parameters of the stationary Beta distribution, Beta(a,b)

a <- # enter de value of s

b <- 1-a


## Function to calculate Tn according to Eq.4 in Ebner & Liebenberg (2021)

calculate_Tn <- function(x, a, b) {
n <- length(x)

# First term

sum1 <- 0
for (j in 1:n) {
for (k in 1:n) {
sum1 <- sum1 + ((a + b)^2 * x[j] * x[k] - a * (a + b) * (x[j] + x[k]) + a^2) * min(x[j], x[k])
}
}
term1 <- sum1 / n
cat("Term 1:", term1, "\n")

# Second term

B1 <- beta(a + 1, b + 1) / beta(a, b)
sum2 <- 0
for (j in 1:n) {
sum2 <- sum2 + ((a + b) * x[j] - a) * pbeta(x[j], a + 1, b + 1)
}
term2 <- 2 * B1 * sum2
cat("Term 2:", term2, "\n")

# Third term

B2 <- beta(2 * a + 1, 2 * b + 1) / beta(a, b)^2
term3 <- n * B2
cat("Term 3:", term3, "\n")

# Calculate Tn

Tn <- term1 - term2 + term3
return(Tn)
}


Tn_original <- calculate_Tn(x, a, b)
cat("Tn from the data:", Tn_original, "\n")


## Generate B bootstrap samples of length n from a Beta(a, b) distribution, and calculate Tn from each one

B <- 5000  
Tn_bootstrap <- numeric(B)
n <- 200

for (i in 1:B) {
sample_beta <- rbeta(n, a, b)
Tn_bootstrap[i] <- calculate_Tn(sample_beta, a, b)
}

## Obtain 1-alpha quantiles from the bootstrap Tn's, with alpha = 0.01, 0.05, 0.1

alpha_levels <- c(0.01, 0.05, 0.1)
quantiles_1_alpha <- sapply(alpha_levels, function(alpha) quantile(Tn_bootstrap, 1 - alpha))

cat("1-alpha quantiles from the distribution of Tn under the null hypothesis:\n")
print(quantiles_1_alpha) 
cat("Tn from the data:", Tn_original, "\n") # Hypothesis of goodness-of-fit is rejected at the level alpha if the observed Tn from the data is greater than the corresponding 1-alpha quantile under the null hypothesis 

## Generate graphical hypothesis testing

library(ggplot2)

df <- data.frame(Tn_bootstrap = Tn_bootstrap)

dens <- density(Tn_bootstrap)
y_max <- max(dens$y)


ggplot(df, aes(x = Tn_bootstrap)) + 
geom_density(fill = "blue", alpha = 0.3) + 
geom_vline(aes(xintercept = Tn_original), color = "red", linetype = "dashed", size = 1) +
geom_vline(aes(xintercept = quantiles_1_alpha[1]), color = "green", linetype = "dotted", size = 1) +
geom_vline(aes(xintercept = quantiles_1_alpha[2]), color = "purple", linetype = "dotted", size = 1) +
geom_vline(aes(xintercept = quantiles_1_alpha[3]), color = "orange", linetype = "dotted", size = 1) +
labs(x = "Tn", y = "Probability density") +


theme_minimal() +
theme(
plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
axis.title = element_text(size = 12),
axis.text = element_text(size = 10)
) +
annotate("text", x = Tn_original, y = y_max / 2, 
label = paste("Tn obs =", round(Tn_original, 5)), 
vjust = -1, hjust = -0.2, color = "red", size = 3.5, angle = 90) +
annotate("text", x = quantiles_1_alpha[1], y = y_max / 2, 
label = paste("quantile 99% =", round(quantiles_1_alpha[1], 5)), 
vjust = -1, hjust = -0.2, color = "green", size = 3.5, angle = 90) +
annotate("text", x = quantiles_1_alpha[2], y = y_max / 2, 
label = paste("quantile 95% =", round(quantiles_1_alpha[2], 5)), 
vjust = -1, hjust = -0.2, color = "purple", size = 3.5, angle = 90) +
annotate("text", x = quantiles_1_alpha[3], y = y_max / 2, 
label = paste("quantile 90% =", round(quantiles_1_alpha[3], 5)), 
vjust = -1, hjust = -0.2, color = "orange", size = 3.5, angle = 90)
