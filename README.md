## QTL-analysis
Script(s) for standard QTL analyses.
### Input data:
The single and multiple QTL mapping R scripts (single\_qtl\_mapping\_rqtl.R, multiple\_qtl\_mapping\_rqtl.R) both use as input an R/qtl cross object. The input R/qtl cross object is constructed and saved with save\_crossobject\_rqtl.R. The [specific cross](https://github.com/MulletLab/leafangle_supplement/blob/master/h2_and_qtl/R07018_x_R07020/R07018xR07020_genetic_map_with_phenotypes.csv) used in the script can be found at the [MulletLab GitHub](https://github.com/MulletLab/leafangle_supplement) which is part of the Supplemental Information from [Truong, McCormick, Rooney & Mullet (Genetics, 2015)](http://www.genetics.org/content/201/3/1229):
* [R07018xR07020_genetic_map_with_phenotypes.csv](https://github.com/MulletLab/leafangle_supplement/blob/master/h2_and_qtl/R07018_x_R07020/R07018xR07020_genetic_map_with_phenotypes.csv)

The scripts used to calculate penalties for multiple QTL model traversal [(as described by Manichaikul _et al._ (Genetics, 2009))](http://www.genetics.org/content/181/3/1077) is also found on the [MulletLab GitHub](https://github.com/MulletLab/leafangle_supplement).

### Requirements:
These scripts (save\_crossobject\_rqtl.R, single\_qtl\_mapping\_rqtl.R, multiple\_qtl\_mapping\_rqtl.R) both require the [R/qtl library](http://rqtl.org/) developed by [Professor Karl Broman](http://kbroman.org/).

## Permutations for QTL-analysis
Scripts to calculate permutations on the Texas A&M Institute for Genome Science & Society High Performance Cluster (tigss-hpc).
This workflow is specific to the tigss-hpc, however can be modified to be of use elsewhere.
### Input data:
This takes as input the cross object that is constructed by save\_crossobject\_rqtl.R.


