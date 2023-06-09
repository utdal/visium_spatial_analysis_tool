# List of required packages
required_packages <- c("Seurat", "ggplot2", "dplyr", "future")


check_install_packages <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

check_install_packages(required_packages)