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