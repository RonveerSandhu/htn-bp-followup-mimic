# R 4.6
# Installs the R package versions used in the paper (called from Python through rpy2).

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

remotes::install_version("lcmm",     version = "2.2.2",  upgrade = "never")
remotes::install_version("cmprsk",   version = "2.2-12", upgrade = "never")
remotes::install_version("survival", version = "3.8-3",  upgrade = "never")
