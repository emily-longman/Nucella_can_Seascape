# ....

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Set path as main Github repo
install.packages(c('rprojroot'))
library(rprojroot)

# List all files and directories below the root
dir(find_root_file(criterion = has_file("README.md")))
root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
setwd(root_path)

# ================================================================================== #

# Load packages (https://diyabc.github.io/gui/)
install.packages("devtools")
library(devtools)
devtools::install_github(
 "diyabc/diyabcGUI",
 subdir = "R-pkg"
)
# Download required binary files (Note: only need to do once)
library(diyabcGUI)
diyabcGUI::dl_all_latest_bin()
# Launch the interfact
library(diyabcGUI)
diyabcGUI::diyabc()

# Note: At the moment, the standalone app is not available for Linux and MacOS users. 
#Nonetheless, Linux and MacOS users can install the diyabcGUI package, c.f. below, and run the DIYABC-RF GUI as a standard shiny app.

# ================================================================================== #
