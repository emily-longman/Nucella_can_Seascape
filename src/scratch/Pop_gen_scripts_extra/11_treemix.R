# Graph Treemix output

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Set path as main Github repo
install.packages(c('rprojroot'))
library(rprojroot)

# List all files and directories below the root
dir(find_root_file(criterion = has_file("README.md")))

# Define root path
root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
setwd(root_path)

# Set working directory as Treemix path from root
#treemix_path_from_root <- find_root_file("data", "processed", "pop_structure", "treemix", "treemix_output_test", criterion = has_file("README.md"))
#setwd(treemix_path_from_root)

# ================================================================================== #

# Load packages
install.packages(c('ggplot2', 'RColorBrewer'))
library(ggplot2)
library(RColorBrewer)

# Get Treemix R plotting script (in src folder of Treemix source code)
source("src/02_pop_structure/11_treemix_plotting_functions.R")

# ================================================================================== #

# Populations
pops <- c("STR", "SLR", "SH", "ARA", "BMR", "CBL", "FC", "FR", "HZD", "SBR", "PSG", "KH", "STC", "PL", "VD", "OCT", "PB", "PGP", "PSN")
write.table(pops, "output/figures/pop_structure/treemix/Treemix_pop_order.txt", quote = F)

# ================================================================================== #

# Set m = 0

# Graph tree
pdf("output/figures/pop_structure/treemix/Treemix_0.pdf", width = 8, height = 8)
plot_tree("data/processed/pop_structure/treemix/treemix_output_test/Treemix.Output.0")
dev.off()

# Residual visualization
pdf("output/figures/pop_structure/treemix/Treemix_0_resid.pdf", width = 8, height = 8)
plot_resid("data/processed/pop_structure/treemix/treemix_output_test/Treemix.Output.0", pop_order = "output/figures/pop_structure/treemix/Treemix_pop_order.txt")
dev.off()


# Set m = 1

# Graph tree
pdf("output/figures/pop_structure/treemix/Treemix_1.pdf", width = 8, height = 8)
plot_tree("data/processed/pop_structure/treemix/treemix_output_test/Treemix.Output.1")
dev.off()
# Residual visualization
pdf("output/figures/pop_structure/treemix/Treemix_1_resid.pdf", width = 8, height = 8)
plot_resid("data/processed/pop_structure/treemix/treemix_output_test/Treemix.Output.1", pop_order = "output/figures/pop_structure/treemix/Treemix_pop_order.txt")
dev.off()

# Set m = 2

# Graph tree
pdf("output/figures/pop_structure/treemix/Treemix_2.pdf", width = 8, height = 8)
plot_tree("data/processed/pop_structure/treemix/treemix_output_test/Treemix.Output.2")
dev.off()
# Residual visualization
pdf("output/figures/pop_structure/treemix/Treemix_2_resid.pdf", width = 8, height = 8)
plot_resid("data/processed/pop_structure/treemix/treemix_output_test/Treemix.Output.2", pop_order = "output/figures/pop_structure/treemix/Treemix_pop_order.txt")
dev.off()

# Set m = 3

# Graph tree
pdf("output/figures/pop_structure/treemix/Treemix_3.pdf", width = 8, height = 8)
plot_tree("data/processed/pop_structure/treemix/treemix_output_test/Treemix.Output.3")
dev.off()
# Residual visualization
pdf("output/figures/pop_structure/treemix/Treemix_3_resid.pdf", width = 8, height = 8)
plot_resid("data/processed/pop_structure/treemix/treemix_output_test/Treemix.Output.3", pop_order = "output/figures/pop_structure/treemix/Treemix_pop_order.txt")
dev.off()

# Set m = 4

# Graph tree
pdf("output/figures/pop_structure/treemix/Treemix_4.pdf", width = 8, height = 8)
plot_tree("data/processed/pop_structure/treemix/treemix_output_test/Treemix.Output.4")
dev.off()
# Residual visualization
pdf("output/figures/pop_structure/treemix/Treemix_4_resid.pdf", width = 8, height = 8)
plot_resid("data/processed/pop_structure/treemix/treemix_output_test/Treemix.Output.4", pop_order = "output/figures/pop_structure/treemix/Treemix_pop_order.txt")
dev.off()

# Set m = 5

# Graph tree
pdf("output/figures/pop_structure/treemix/Treemix_5.pdf", width = 8, height = 8)
plot_tree("data/processed/pop_structure/treemix/treemix_output_test/Treemix.Output.5")
dev.off()
# Residual visualization
pdf("output/figures/pop_structure/treemix/Treemix_5_resid.pdf", width = 8, height = 8)
plot_resid("data/processed/pop_structure/treemix/treemix_output_test/Treemix.Output.5", pop_order = "output/figures/pop_structure/treemix/Treemix_pop_order.txt")
dev.off()

# Set m = 6

# Graph tree
pdf("output/figures/pop_structure/treemix/Treemix_6.pdf", width = 8, height = 8)
plot_tree("data/processed/pop_structure/treemix/treemix_output_test/Treemix.Output.6")
dev.off()
# Residual visualization
pdf("output/figures/pop_structure/treemix/Treemix_6_resid.pdf", width = 8, height = 8)
plot_resid("data/processed/pop_structure/treemix/treemix_output_test/Treemix.Output.6", pop_order = "output/figures/pop_structure/treemix/Treemix_pop_order.txt")
dev.off()