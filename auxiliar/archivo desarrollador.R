devtools::load_all()
devtools::document()

devtools::install()

library(RelDists)

# Para hacer un chequeo
devtools::check(remote=TRUE)

# Para crear la pagina web
pkgdown::build_site()

# Badges
usethis::use_github_action_check_standard()

