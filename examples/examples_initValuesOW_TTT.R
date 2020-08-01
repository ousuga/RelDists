# Example 1
# Bathtuh hazard and its corresponding TTT plot
y <- rOW(n=50, mu=0.022, sigma=8, nu=0.01)

my_initial_guess <- initValuesOW_TTT(formula=y~1)
curve(hOW(x, mu=0.022, sigma=8, nu=0.01), from=0, to=80, ylim=c(0, 0.04),
      col="red", ylab="Hazard function", las=1)
