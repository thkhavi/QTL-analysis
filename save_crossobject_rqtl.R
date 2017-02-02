# Write the crossobject after sim geno for downstream multiple-QTL mapping with R/qtl package
# (Much of this code originates from rqtl.org tutorials.)    

# Change these file paths to reflect your system:
input_file_directory = file.path("./")
input_file_name_cross = file.path("./R07018xR07020_genetic_map_with_phenotypes.csv")

# Generation interval of cross
generation_interval = 6
# Backcross interval of cross
backcross_interval = 0

library(qtl)

# Read in cross
cross_inputcross = read.cross(format="csvr", 
                              dir=input_file_directory,
                              file=input_file_name_cross, 
                              BC.gen=backcross_interval, 
                              F.gen=generation_interval,
                              genotypes=c("AA","AB","BB","notAA","notBB"))

# Choose the distance (in cM) to thin out
marker_distance = 1.0
# If you don't want to thin markers, then comment out with "#"
cross_inputcross_map = pull.map(cross_inputcross) 
markers2keep = lapply(cross_inputcross_map, pickMarkerSubset, min.distance=marker_distance) 
cross_sub = pull.markers(cross_inputcross, unlist(markers2keep))
cross_inputcross = cross_sub

# If the cross type is considered a RIL keep the next line, if not comment out with "#"
cross_inputcross = convert2riself(cross_inputcross)
cross_inputcross = jittermap(cross_inputcross)
cross_inputcross = fill.geno(cross_inputcross, 
                             map.function = "haldane")
cross_inputcross = calc.genoprob(cross_inputcross,
                                 map.function = "haldane")
cross_inputcross = sim.geno(cross_inputcross, 
                            map.function = "haldane")

save(cross_inputcross, file="R07018_R07020_2017-01-30")
