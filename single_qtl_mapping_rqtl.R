# Estimate QTL using interval mapping with R/qtl package
# (Much of this code originates from rqtl.org tutorials.)              

# Change these file paths to reflect your system:
input_file_directory = file.path("./")
input_file_name_cross = file.path("./R07018xR07020_genetic_map_with_phenotypes.csv")

working_directory = getwd()
if (!is.null(working_directory)) setwd(input_file_directory)

# Generation interval of cross
generation_interval = 6
# Backcross interval of cross
backcross_interval = 0

# Available phenotypes:
phenotype_list = c("angle_leaf_3_avg_gh204A_2013_normalized",
                   "angle_leaf_4_avg_gh204A_2013_normalized",
                   "angle_leaf_3_avg_csfield_2014_rep1_normalized",
                   "angle_leaf_4_avg_csfield_2014_rep1_normalized",
                   "angle_leaf_5_avg_csfield_2014_rep1_normalized",
                   "angle_leaf_3_avg_csfield_2014_rep2_normalized",
                   "angle_leaf_4_avg_csfield_2014_rep2_normalized",
                   "angle_leaf_5_avg_csfield_2014_rep2_normalized")

# Permutations (R/qtl developers expect at least 1000) & alpha
number_permutations = 10
alpha_value = 0.05

library(qtl)

# Read in cross
cross_inputcross = read.cross(format="csvr", 
                               dir=input_file_directory,
                               file=input_file_name_cross, 
                               BC.gen=backcross_interval, 
                               F.gen=generation_interval,
                               genotypes=c("AA","AB","BB","notAA","notBB"))

# If the cross type is considered a RIL keep the next line, if not comment out with "#"
# cross_inputcross = convert2riself(cross_inputcross)

cross_inputcross = calc.genoprob(cross_inputcross,
                                 map.function = "haldane")
cross_inputcross = sim.geno(cross_inputcross, 
                            map.function = "haldane")

# For each phenotype:
for (phenotype in 1:(length(phenotype_list))) {
  ptm = proc.time()
  # Create directory
  directory = paste0(phenotype_list[phenotype], "_interval_mapping_folder")
  dir.create(directory)
  setwd(file.path(directory))
  
  # Interval Mapping
  print(phenotype_list[phenotype])
  cat("\t interval mapping \n")
  scanone_inputcross = scanone(cross = cross_inputcross,
                                method = c("em"),
                                pheno.col = phenotype_list[phenotype])
  
  # Permutations
  cat("\t permutations \n")
  scanone_inputcross_perm = scanone(cross = cross_inputcross,
                                     method=c("em"),
                                     pheno.col=phenotype_list[phenotype],
                                     n.perm=number_permutations)
  
  cat("\t", summary(scanone_inputcross_perm, alpha = alpha_value), "\n")
  print(proc.time() - ptm)
  
  cat("\t creating plots \n")
  scanone_plot_name = paste("scanone", phenotype_list[phenotype], sep = " - ")
  scanone_file_plot_name = paste0("scanone", "-", phenotype_list[phenotype], ".png")
  
  # Set y-limit on plot to be above the largest lod (either from permutations or single-QTL analysis)
  ylimit = max(summary(scanone_inputcross_perm, alpha = alpha_value), max(scanone_inputcross$lod))
  
  png(scanone_file_plot_name, width=1600, height=800, res=200)
  plot(scanone_inputcross, col = "blue", 
       lty=1, 
       ylim=c(0, ylimit), 
       main=scanone_plot_name)
  add.threshold(scanone_inputcross, 
                perms=scanone_inputcross_perm, 
                alpha=alpha_value, 
                col="red", 
                lty=2)
  dev.off()
  
  # Additive and dominance (if not RIL) effect plot
  effect_plot_name = paste("effect", phenotype_list[phenotype], sep = " - ")
  effect_file_plot_name = paste0("effect","-", phenotype_list[phenotype], ".png")

  png(effect_file_plot_name, width=1600, height=800, res=200)
  effectscan(cross_inputcross, 
             pheno.col=phenotype_list[phenotype], 
             get.se=FALSE, 
             draw=TRUE, 
             main=effect_plot_name)
  dev.off() 
  
  # Construct list of significant markers (alpha=0.05)
  significant_markers_file_name = paste0("signficant_markers_", phenotype_list[phenotype],".csv")
  sink(file=significant_markers_file_name, append = FALSE)
  cat("Significant markers for \n",
      "\t phenotype:", phenotype_list[phenotype], "\n",
      "\t cross object:", input_file_name_cross, "\n",
      "\t alpha-value:", alpha_value, "\n",
      "\t permutations:", number_permutations, "\n", 
      "\t lod threshold:", summary(scanone_inputcross_perm, alpha = alpha_value), "\n\n")
  cat("chr_position(bp)_variant", ",", "position(cM)",	",", "lod","\n")
  sink()
  for (marker in 1:(length(scanone_inputcross$lod))) {
    if (scanone_inputcross$lod[marker] >= summary(scanone_inputcross_perm, alpha = alpha_value)) {
      sink(file = significant_markers_file_name, append = TRUE)
      cat(rownames(scanone_inputcross)[marker], ",", scanone_inputcross$pos[marker], ",", scanone_inputcross$lod[marker],"\n")
      sink()
    }
  }
  
  # Construct folder containing phenotype x genotype plots of significant markers (alpha=0.05)
  # Create p x g directory
  pxg_directory = "pxg_plots"
  dir.create(pxg_directory)
  setwd(file.path(pxg_directory))
  for (marker in 1:(length(scanone_inputcross$lod))) {
    if (scanone_inputcross$lod[marker] >= summary(scanone_inputcross_perm, alpha = alpha_value)) {
      PXG_name = paste0("pxg_", rownames(scanone_inputcross)[marker], ".png")
      png(PXG_name, width=1600, height=800, res=200)
      plotPXG(cross_inputcross, 
              rownames(scanone_inputcross)[marker],
              pheno.col=phenotype_list[phenotype],
              jitter=1,
              infer=TRUE,
              main=rownames(scanone_inputcross)[marker])
      dev.off()
    }
  }
  setwd("../../")
}
