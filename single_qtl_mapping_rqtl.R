# Estimate QTL using interval mapping with R/qtl package
# (Much of this code originates from rqtl.org tutorials.)              

# Load cross object
load("./R07018_R07020_2017-01-30")

# Cross name
input_file_name_cross = "R07018_R07020_2017-01-30"

# Available phenotypes:
phenotype_list = c("angle_leaf_3_avg_gh204A_2013_normalized",
                   "angle_leaf_4_avg_gh204A_2013_normalized")

alpha_value = 0.05
number_permutations = 25000
main_permutation_threshold = 3.205457

library(qtl)

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
  
  cat("\t creating plots \n")
  scanone_plot_name = paste("scanone", phenotype_list[phenotype], sep = " - ")
  scanone_file_plot_name = paste0("scanone", "-", phenotype_list[phenotype], ".png")
  
  # Set y-limit on plot to be above the largest lod (either from permutations or single-QTL analysis)
  ylimit = max(main_permutation_threshold, max(scanone_inputcross$lod))
  
  png(scanone_file_plot_name, width=1600, height=800, res=200)
  plot(scanone_inputcross, col = "green", 
       lty=1, 
       ylim=c(0, ylimit), 
       main=scanone_plot_name)
  dev.off()
  
  # Additive and dominance (if not RIL) effect plot
  effect_plot_name = paste("effect", phenotype_list[phenotype], sep = " - ")
  effect_file_plot_name = paste0("effect","-", phenotype_list[phenotype], ".png")
  
  png(effect_file_plot_name, width=1600, height=800, res=200)
  effectscan_values = effectscan(cross_inputcross, 
                                 pheno.col=phenotype_list[phenotype], 
                                 get.se=TRUE, 
                                 draw=TRUE, 
                                 main=effect_plot_name)
  dev.off()
  
  # Construct list of significant markers (alpha=0.05)
  significant_markers_file_name = paste0("significant_markers_", phenotype_list[phenotype],".csv")
  sink(file=significant_markers_file_name, append = FALSE)
  cat("Significant markers for \n",
      "\t phenotype:", phenotype_list[phenotype], "\n",
      "\t cross object:", input_file_name_cross, "\n",
      "\t alpha-value:", alpha_value, "\n",
      "\t permutations:", number_permutations, "\n", 
      "\t lod threshold:", main_permutation_threshold, "\n\n")
  cat("chr_position(bp)_variant", ",", "position(cM)",	",", "lod", ",", "effect(a)", ",", "SE","\n")
  sink()
  for (marker in 1:(length(scanone_inputcross$lod))) {
    if (scanone_inputcross$lod[marker] >= main_permutation_threshold) {
      sink(file = significant_markers_file_name, append = TRUE)
      cat(rownames(scanone_inputcross)[marker], ",",
          scanone_inputcross$pos[marker], ",",
          scanone_inputcross$lod[marker], ",",
          effectscan_values$a[marker], ",",
          effectscan_values$se.a[marker],"\n")
      sink()
    }
  }
  # Save the lod intervals. Note that this assumes 10 chromosomes.
  sink(file=significant_markers_file_name, append = TRUE)
  cat("\nLOD-2 intervals for all chromosomes:\n")
  for (k in 1:10) {
    print(lodint(scanone_inputcross, drop=2, chr=k, expandtomarkers=TRUE))
  }
  sink()
  setwd("../")
}
