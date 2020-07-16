devtools::load_all()
devtools::document()

devtools::install()

require(RelDists)

devtools::check(remote=TRUE)

pkgdown::build_site()
