# rhizoscale
Rcodes for a study on grapevine rhizosphere scales


This README and repository is currently under development.


A publication that these codes were developed for is under review. If you use any of these codes for a publication, please check back prior to publishing for an updated citation. 

All RCodes in this project are based on Phyloseq objects with ASVs for bacteria (16s) and fungi (ITS). There are 5 separate codes:

1PhyloseqObjectCreation
- 1)creates MEMs based on the coordinates used in this study for later use
- 2) converts categorical variables to dummy variables,
- 3) takes initial 16s and ITS phyloseq objects creates 11 separate Phyloseq objects that vary in either filtration parameters or taxonomic level


2VariationPartitioningAnalysis
- Based on Borcard et al., 2018 (https://doi.org/10.1007/978-1-4419-7976-6)
- Takes each of the 11 Phyloseq objects and completes variation partitioning on them.

3ShannonDiversity
- Computes unrarefied Shannon diversity on the phyloseq objects and tests differences with mixed models

4DifferentialAbundance
- uses the CornCob package in R to complete differential abundance on phyloseq objects at the ASV taxonomic level

5BCdispersion
- calculates bray curtis dissimilarities and tests differences in community dispersion between old and young roots

