# Estimate QTL using multiple mapping with R/qtl package
# (Much of this code originates from rqtl.org tutorials.)    

# Load cross object
load("./R07018xR07020_2017-01-30")

# Phenotype: 
phenotype = c("angle_leaf_3_avg_gh204A_2013_normalized")

library(qtl)

################################################################################
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
# Initialize the multiple-QTL-model traversal by adding in the markers' positions (chromosome and cM)
# Here we used the results from single QTL analysis from single_qtl_mapping_rqtl.R
init_qtl = makeqtl(cross_inputcross, 
                    chr=c(1,6),
                    pos=c(117.6491,103.1357))

# Put in the formula that you would like to initialize with the QTL used in init_qtl
# Penalities calcluated from scantwo (permutations = 24,000, alpha = 0.05, Truong et al (2015))
stepout2 = stepwiseqtl(cross_inputcross, 
                        pheno.col=phenotype, 
                        penalties=c(3.679856, 5.156491, 2.713436),  
                        qtl=init_qtl,          
                        formula=y~Q1+Q2,
                        max.qtl=3,
                        scan.pairs=FALSE,
                        refine.locations=FALSE,
                        keeptrace=TRUE,
                        verbose=TRUE)
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
################################################################################
phenotype_directory = paste0(phenotype, "_multiple_qtl_mapping_folder")
dir.create(file.path(phenotype_directory))
setwd(file.path(phenotype_directory))

# all models traversed
thetrace = attr(stepout2, "trace")

# Make file for model selection of best fit model
all_pLOD_model_file_name = paste0("pLOD_models_", phenotype, ".txt")
sink(file = all_pLOD_model_file_name)
cat("QTL models traversed with associated pLOD for phenotype: \n")
cat(phenotype)
cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
print(thetrace)
cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
sink()

qtl_formula = attr(stepout2, "formula")

# write .txt lod2interval file
lod2int_file_name = paste0("lod2interval_for_model_",
                           phenotype,
                           ".txt")
sink(file = lod2int_file_name)
cat("Refined QTL from chosen model (within model selection) with pLOD for phenotype: \n")
cat(phenotype)
cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Fit model statistics\n\n")
print(summary(fitqtl(cross=cross_inputcross, pheno.col=phenotype, qtl=stepout2, formula = qtl_formula)))
cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
cat("LOD 2 intervals of the QTL model\n\n")
for(qtl in 1:length(stepout2$name)){ 
  cat("Q")
  cat(qtl)
  cat("\n")
  print(lodint(stepout2, qtl.index=qtl , drop=2, expandtomarkers=TRUE))
  cat("\n")
}
sink()

# write lod score of entire model
lod2int_file_name = paste0("lods_for_multiple-QTL_model_",
                           phenotype,
                           ".txt")
sink(file = lod2int_file_name)
print(attr(stepout2, "lodprofile"))
sink()

# plot the qtl model's lod scores
qtl_model_lod_file_name = paste0("lod_plot_for_multiple-QTL_model_",
                                 phenotype,
                                 ".png")
png(file = qtl_model_lod_file_name,  
    width = 20, 
    height = 10, 
    units ="in", 
    res = 300,
    g = "white")
par(mfrow=c(1,1))
plotLodProfile(stepout2, showallchr=TRUE)
dev.off()

setwd("../")
