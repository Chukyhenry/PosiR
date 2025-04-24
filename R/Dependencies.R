#install.packages(c("usethis", "devtools", "roxygen2", "testthat"))
library(usethis)
use_r("simultaneous_ci")
use_testthat()
use_roxygen_md()
use_github()
use_readme_rmd()
use_package("MASS")
use_package("stats")
use_github_action("check-standard")
usethis::use_package("pbapply", "Suggests")
renv::install("pbapply")
renv::snapshot()
devtools::document()
usethis::gitcreds_get()
usethis::use_readme_rmd()
# For parallel processing (base R, but good to declare import)
usethis::use_package("parallel")
# pbapply is already in Suggests, which is fine as it's optional

# For the new plot method
usethis::use_package("graphics", type = "Imports") # For basic plotting functions
# For parallel processing (base R, but good to declare import)
usethis::use_package("parallel")
# pbapply is already in Suggests, which is fine as it's optional

# For the new plot method
usethis::use_package("graphics", type = "Imports") # For basic plotting functions
# Creates vignettes/PosiR-introduction.Rmd
usethis::use_vignette("PosiR-introduction")
# Creates tests/testthat.R and tests/testthat/test-simultaneous_ci.R
usethis::use_test("simultaneous_ci")

usethis::use_git_remote("origin", url = NULL, overwrite = TRUE)

devtools::document()

usethis::use_git()
usethis::use_github()

