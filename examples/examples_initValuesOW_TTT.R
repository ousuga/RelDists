# Example 1
# Bathtuh hazard and its corresponding TTT plot
y1 <- rOW(n = 1000, mu = 0.1, sigma = 7, nu = 0.08)
my_initial_guess1 <- initValuesOW_TTT(formula=y1~1)
summary(my_initial_guess1)
plot(my_initial_guess1, par_plot=list(mar=c(3.7,3.7,1,2.5),
                                     mgp=c(2.5,1,0)))

curve(hOW(x, mu = 0.022, sigma = 8, nu = 0.01), from = 0, 
      to = 80, ylim = c(0, 0.04), col = "red", 
      ylab = "Hazard function", las = 1)

# Example 2
# Bathtuh hazard and its corresponding TTT plot with right censored data
\dontrun{
y2 <- rOW(n = 1000, mu = 0.1, sigma = 7, nu = 0.08)
status <- c(rep(1, 980), rep(0, 20))
my_initial_guess2 <- initValuesOW_TTT(formula=Surv(y2, status)~1)
summary(my_initial_guess2)
plot(my_initial_guess2, par_plot=list(mar=c(3.7,3.7,1,2.5),
                                     mgp=c(2.5,1,0)))

curve(hOW(x, mu = 0.022, sigma = 8, nu = 0.01), from = 0, 
      to = 80, ylim = c(0, 0.04), col = "red", 
      ylab = "Hazard function", las = 1)
}