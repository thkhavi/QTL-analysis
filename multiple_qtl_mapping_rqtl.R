# Estimate QTL using multiple mapping with R/qtl package
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
# Phenotype: 
phenotype = c("angle_leaf_3_avg_gh204A_2013_normalized")

library(qtl)

# Read in cross
cross_inputcross = read.cross(format="csvr", 
                              dir=input_file_directory,
                              file=input_file_name_cross, 
                              BC.gen=backcross_interval, 
                              F.gen=generation_interval,
                              genotypes=c("AA","AB","BB","notAA","notBB"))

# If the cross type is considered a RIL keep the next line, if not comment out with "#"
cross_inputcross = convert2riself(cross_inputcross)

################################################################################
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
# multiple-QTL-model traversal as currently implemented (10/2014) is unable to handle > 1400 markers
# If your genetic map is made of more than 1400 markers, you may need to thin out markers:
# Choose the distance (in cM) to thin out
# marker_distance = 0.5
# If you don't want to thin markers, then comment out with "#"
# cross_inputcross_map = pull.map(cross_inputcross) 
# markers2keep = lapply(cross_inputcross_map, pickMarkerSubset, min.distance=marker_distance) 
# cross_sub = pull.markers(cross_inputcross, unlist(markers2keep))
# cross_inputcross = cross_sub
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
################################################################################

cross_inputcross = calc.genoprob(cross_inputcross,
                                 map.function = "haldane")
cross_inputcross = sim.geno(cross_inputcross, 
                            map.function = "haldane")

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

# model selection traversal images
model_plot_directory = "./model_plots"
dir.create(file.path(model_plot_directory))
setwd(file.path(model_plot_directory))

i_models_in_row = 2
i_models_in_columns = 1

model_selection_i_file_name = paste0("models_",
                                     1,
                                     "-",
                                     (1+(i_models_in_row * i_models_in_columns)-1),
                                     ".png")


png(file = model_selection_i_file_name,  
    width = 6.5,
    height = 9, 
    units ="in", 
    res = 300, 
    bg = "white")

par(mfrow=c(i_models_in_row,i_models_in_columns))

i_models_in_image = 0
for(i in seq(along=thetrace)){
  if ((i_models_in_row * i_models_in_columns) == i_models_in_image){
    dev.off()
    i_models_in_image = 0
    model_selection_i_file_name = paste0("models_",
                                         i,
                                         "-",
                                         (i+(i_models_in_row * i_models_in_columns)-1),
                                         ".png")
    png(file = model_selection_i_file_name,  
        width = 6.5, 
        height = 9, 
        units ="in",
        res = 300, 
        bg = "white")
    par(mfrow=c(i_models_in_row,i_models_in_columns))
    
  }
  plotModel(thetrace[[i]], 
            chronly=FALSE, 
            circrad.rel=0.5, 
            col = "thistle1", 
            main = paste("model", i, ": pLOD =", round(attr(thetrace[[i]], "pLOD"), 2)))
  i_models_in_image = i_models_in_image + 1
  
}
dev.off()
setwd("../")

# Make file for model selection of best fit model
all_pLOD_model_file_name = paste0("pLOD_models_", phenotype, ".txt")
sink(file = all_pLOD_model_file_name)
cat("QTL models traversed with associated pLOD for phenotype: \n")
cat(phenotype)
cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
print(thetrace)
cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
sink()

################################################################################
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
# Choose a model to examine from "pLOD_models_file"
# Refine QTL positions
# Here I chose a model 2 with pLOD = -1.451
# even though it's not the largest pLOD 
# to illustrate analyses of the interaction in later steps
# Q#  chr pos
# chr       pos
# Q1   1 117.64908
# Q2   6  89.02274
# Q3   2 121.63643
# Q4   8  98.63018
# Formula: y ~ Q1 + Q2 + Q3 + Q4 + Q3:Q4 
# pLOD:  -1.451 
make_qtl = makeqtl(cross=cross_inputcross, 
                   chr=c(1,6,2,8),
                   pos=c(117.64908,89.02274,121.63643,98.63018))

qtl_formula = y ~ Q1 + Q2 + Q3 + Q4 + Q3:Q4

model_chosen = 2
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
################################################################################
model_stat_directory = paste0("./model_", model_chosen)
dir.create(file.path(model_stat_directory))
setwd(file.path(model_stat_directory))

# refine positions of qtl
refine_qtl = refineqtl(cross = cross_inputcross,
                        pheno.col=phenotype,
                        qtl = make_qtl,
                        formula = qtl_formula,
                        incl.markers=TRUE,
                        keeplodprofile=TRUE)

# main qtl phenotype by genotype effect images
main_effect_directory = file.path("./main_qtl_effect_plot")
dir.create(file.path(main_effect_directory))
setwd(file.path(main_effect_directory))

i_QTL_in_row = 4
i_QTL_in_columns = 1

main_QTL_i_name = paste0(1,
                         "-",
                         ((i_QTL_in_row * i_QTL_in_columns)),
                         ".png")
png(file = main_QTL_i_name, 
    width = 6.5, 
    height = 9, 
    units ="in", 
    res = 300, 
    bg = "white")
par(mfrow=c(i_QTL_in_row,i_QTL_in_columns))

i_QTL_in_image = 0
for(i in 1:length(refine_qtl$name)){
  if ((i_QTL_in_row * i_QTL_in_columns) == i_QTL_in_image){
    dev.off()
    i_QTL_in_image = 0
    main_QTL_i_name = paste(i,
                            "-",
                            (i+(i_QTL_in_row * i_QTL_in_columns)-1),
                            ".png")
    png(file = main_QTL_i_name,  
        width = 6.5, 
        height = 9, 
        units ="in", 
        res = 300, 
        bg = "white")
    par(mfrow=c(i_QTL_in_row,i_QTL_in_columns))
  }
  marker_name = find.marker(cross_inputcross, refine_qtl$chr[i], refine_qtl$pos[i])
  plotPXG(cross_inputcross, 
          marker_name, 
          pheno.col= phenotype)
  i_QTL_in_image = i_QTL_in_image + 1
}
dev.off()
setwd("../")

# write .txt lod2interval file
lod2int_file_name = paste0("lod2interval_for_model_",
                           model_chosen,
                           ".txt")
sink(file = lod2int_file_name)
cat("Refined QTL from chosen model (within model selection) with pLOD for phenotype: \n")
cat(phenotype)
cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Fit model statistics\n\n")
print(summary(fitqtl(cross=cross_inputcross, pheno.col=phenotype, qtl=refine_qtl, formula = qtl_formula)))
cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
cat("LOD 2 intervals of the QTL model\n\n")
for(qtl in 1:length(refine_qtl$name)){ 
  cat("Q")
  cat(qtl)
  cat("\n")
  print(lodint(refine_qtl, qtl.index=qtl , drop=2, expandtomarkers=TRUE))
  cat("\n")
}
cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
cat("All lod scores for markers on chromosomes of the QTL model\n\n")
print(attr(refine_qtl, "lodprofile"))
cat("\n")
cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
sink()

# plot the qtl model's lod scores
qtl_model_lod_file_name = paste0("lod_plot_for_model_",
                                 model_chosen,
                                 ".png")
png(file = qtl_model_lod_file_name,  
    width = 20, 
    height = 10, 
    units ="in", 
    res = 300,
    g = "white")
par(mfrow=c(1,1))
plotLodProfile(refine_qtl)
dev.off()

################################################################################
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
################################################################################
# If there were interactions from the model selected; output respective plots: #
################################################################################
# In this case there was an interaction between Q3 and Q7
# Formula: y ~ Q1 + Q2 + Q3 + Q4 + Q3:Q4 
# Change to reflect QTL that have interactions:
QTL_a = 3
QTL_b = 4
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
################################################################################

interaction_directory = paste0("./interaction_", refine_qtl$name[QTL_a], "-", refine_qtl$name[QTL_b])
dir.create(file.path(interaction_directory))
setwd(file.path(interaction_directory))

# effect/interaction plot
i_effects_in_row = 1
i_effects_in_columns = 2

effect_plot_file_name = paste0("effect_",
                                    refine_qtl$name[QTL_a],
                                    "-", 
                                    refine_qtl$name[QTL_b], 
                                    ".png")
png(file = effect_plot_file_name,  
    width = 20, 
    height = 10, 
    units ="in", 
    res = 300,
    g = "white")
par(mfrow=c(i_effects_in_row,i_effects_in_columns))
effectplot(cross_inputcross,pheno.col=phenotype,mname1=refine_qtl$name[QTL_a],mname2=refine_qtl$name[QTL_b])
effectplot(cross_inputcross,pheno.col=phenotype,mname1=refine_qtl$name[QTL_b],mname2=refine_qtl$name[QTL_a])
dev.off()

# phenotype by the two genotype plots
pxg_plot_file_name = paste0("pxg_",
                            refine_qtl$name[QTL_a],
                            "-",
                            refine_qtl$name[QTL_b],
                            ".png")
png(file = pxg_plot_file_name,  
    width = 20, 
    height = 10, 
    units ="in", 
    res = 300, 
    bg = "white")
plotPXG(x = cross_inputcross,
        pheno.col = phenotype,
        marker = c(find.marker(cross_inputcross,refine_qtl$chr[QTL_a],refine_qtl$pos[QTL_a]), 
                   find.marker(cross_inputcross,refine_qtl$chr[QTL_b],refine_qtl$pos[QTL_b])))
dev.off()
setwd("../../")
