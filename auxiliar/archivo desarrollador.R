devtools::load_all()
devtools::document()

devtools::install()

require(RelDists)

# Para hacer un chequeo
devtools::check(remote=TRUE)

# Para crear la pagina web
pkgdown::build_site()
