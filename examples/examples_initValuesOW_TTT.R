# Example 1
# Bathtuh hazard and its corresponding TTT plot
y1 <- rOW(n=50, mu=0.022, sigma=8, nu=0.01)
my_initial_guess <- initValuesOW_TTT(formula=y1~1)
summary(my_initial_guess)
plot(my_initial_guess, legend_options=list(pos=1.04),
     par_plot=list(mar=c(3.7,3.7,1,1.5), mgp=c(2.5,1,0)))

curve(hOW(x, mu = 0.022, sigma = 8, nu = 0.01), from = 0, 
      to = 80, ylim = c(0, 0.04), col = "red", 
      ylab = "Hazard function", las = 1)

# Example 2
# Parameters explained with a covariate
set.seed(789)
x <- c(rep(0,25), rep(1,25))
nu <- 0.1 + 1.8*x
y2 <- rOW(n=50, mu=0.05, sigma=2, nu=nu)
x <- as.factor(x)

# partitions <- list(method='quantile-based', folds=5)
# my_initial_val <- initValuesOW_TTT(y2 ~ x, partition_method=partitions)
my_initial_val <- initValuesOW_TTT(y2 ~ x)
summary(my_initial_val)
plot(my_initial_val, legend_options=list(pos=1.03),
     par_plot=list(mar=c(3.7,3.7,1,2.8), mgp=c(2.5,1,0)))
