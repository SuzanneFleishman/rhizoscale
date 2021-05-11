# rhizoscale
Rcodes for a study on grapevine rhizosphere scales


This README and repository is currently under development.


A publication that these codes were developed for is under review. If you use any of these codes for a publication, please check back prior to publishing for an updated citation. 

This project involved several analyses based on phyloseq objects containing ASVs. 

1PhyloseqObjectCreation
--This code 1)creates MEMs based on the coordinates used in this study for later use, 2) converts categorical variables to dummy variables,
3) takes initial 16s and ITS phyloseq objects creates 11 separate Phyloseq objects that vary in either filtration parameters or taxonomic level


2VariationPartition
--This code is based on Borcard et al., 2018.
--This code takes each of the 11 Phyloseq objects and completes variation partitioning on them.

