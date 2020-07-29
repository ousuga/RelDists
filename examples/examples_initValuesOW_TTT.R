# Example 1
# Bathtuh hazard and its corresponding TTT plot
y <- c(5, 11, 21, 31, 46, 75, 98, 122, 145, 165, 195, 224, 245, 293, 
       321, 330, 350, 420)

my_initial_guess <- initValuesOW_TTT(formula=y~1)
curve(hOW(x, mu=0.022, sigma=8, nu=0.01), from=0, to=80, ylim=c(0, 0.04),
      col="red", ylab="Hazard function", las=1)
