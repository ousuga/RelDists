#Este es un mensaje de prueba (Jaime Mosquera)
.onAttach <- function(libname, pkgname){
  initmessage <- "<<<<<<<<<<<<<<<<<<<<<   RelDists Version 1.0.0  >>>>>>>>>>>>>>>>>>>>>
Feel free to report bugs in https://github.com/ousuga/RelDists/issues"
  packageStartupMessage(initmessage)
  invisible()
}