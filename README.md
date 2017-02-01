## QTL-analysis
Script(s) for standard QTL analyses. QTL analyses use penalties calculated from one-,two- dimensional scans/permutations for single-, multiple- QTL mapping via R/qtl. This includes bash scripts to calculate permutations on the Texas A&M Institute for Genome Science & Society High Performance Cluster (tigss-hpc). 

### Input data:
The single and multiple QTL mapping R scripts (single\_qtl\_mapping\_rqtl.R, multiple\_qtl\_mapping\_rqtl.R) both use as input an R/qtl cross object. The input R/qtl cross object is constructed and saved with save\_crossobject\_rqtl.R. The [specific cross](https://github.com/MulletLab/leafangle_supplement/blob/master/h2_and_qtl/R07018_x_R07020/R07018xR07020_genetic_map_with_phenotypes.csv) used in the script can be found at the [MulletLab GitHub](https://github.com/MulletLab/leafangle_supplement) which is part of the Supplemental Information from [Truong, McCormick, Rooney & Mullet (Genetics, 2015)](http://www.genetics.org/content/201/3/1229):
* [R07018xR07020_genetic_map_with_phenotypes.csv](https://github.com/MulletLab/leafangle_supplement/blob/master/h2_and_qtl/R07018_x_R07020/R07018xR07020_genetic_map_with_phenotypes.csv)

The scripts used to calculate penalties for multiple QTL model traversal [(as described by Manichaikul _et al._ (Genetics, 2009))](http://www.genetics.org/content/181/3/1077) is also found on the [MulletLab GitHub](https://github.com/MulletLab/leafangle_supplement).

### Requirements:
These scripts (save\_crossobject\_rqtl.R, single\_qtl\_mapping\_rqtl.R, multiple\_qtl\_mapping\_rqtl.R) require the [R/qtl library](http://rqtl.org/) developed by [Professor Karl Broman](http://kbroman.org/), and are not tied to any architecture.

Permutation bash scripts found in the tigss-hpc\_rqtl\_penalty\_folders.zip are specific to the tigss-hpc architecture, however the workflow can be modified to be of use elsewhere.

### Typical workflow:
1. Construct cross object with single\_qtl\_mapping\_rqtl.R
2. Calculate penalites for QTL mapping in compute cluster (tigss-hpc, see configuration\_tigss-hpc\_steps.txt to set up)
 1. Modify parameters and filenames in folders:
    * /data/{username}/rqtl\_crosses/
    * /data/{username}/rqtl\_mqm\_scripts/
 2. Submit jobs (`$ ./data/rqtl\_mqm\_scripts/run\_two-dimensional\_scans\_job\_array\_{crossname}.sh`)
 3. Combine penalties (`$ ./data/rqtl\_mqm\_scripts/run\_combine\_scans\_{crossname}.sh`)
 4. Get penalty scores from (`$ vim /data/rqtl\_mqm\_output/{crossname}//mqm\_scantwo\_penalties/0p05\_penalities.txt`)
3. Run single-QTL mapping with single\_qtl\_mapping\_rqtl.R with 
 * crossobject from step 1 
 * main penalties from step 2
4. Run multiple-QTL mapping with multiple\_qtl\_mapping\_rqtl.R with
 * crossobject from step 1 
 * main penalties from step 2
 * initial model with main QTL identified from step 3
 
