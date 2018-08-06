
#' mailto:j.s.hoyer@wustl.edu
#' Started 2016-06-19.
#' Adapted 2016-08-29.

#' ** Load packages; Define axis label strings
library(dplyr)

indx <- "Index of first true leaf with at least one abaxial trichome"

bw <- "Blade width (mm)"
pl <- "Petiole length (mm)"
bl <- "Blade length (mm)"
blw  <- "Blade length:Blade width ratio"
blpl <- "Blade length:Petiole length ratio"


#' ** Workup
workupTrichCaliper <- function(
                               dataframe,
                               ...)
{
    dataframe <- mutate(dataframe, lwratio = bladeL / bladeW)
    dataframe <- mutate(dataframe, bpratio = bladeL / petioleL)

    return(dataframe)
}
