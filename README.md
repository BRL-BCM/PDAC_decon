# R script for extended EDec methodology 

Make sure the original EDec package (https://github.com/BRL-BCM/EDec) is installed in your R library prior to running.

This script is intended to group samples with similar cell-type specific expression of a gene of interest for stage 2 of EDec and returns estimated average cell type specific gene expression profiles for each group.

**Arguments** <br>
goi = A string indicating the gene of interest for which you would like to group. Note: This must be exact spelling and correspond to the row in the bulk expression matrix. <br>
Tum = A matrix containing the bulk tissue gene expression values for all genes (rows) across all samples (columns). Note this is the same input as the original EDec stage 2.  <br>
propTum = A matrix containing all the EDec stage 1 predicted proportions of all cell types (columns) for all samples (rows). Note this is the same input as the original EDec stage 2. <br>
k = An integer indicating the number of groups you would like to split into (Default = 2).<br>

**Value** <br>
A list containing the original output of EDec stage 2 for each group separately. 

**Example usage** <br>
split_groups <- EDec_stage2_aroundGOI(goi = "MSLN",Tum = bulk_gene_expression_matrix ,propTum = cell_type_props_from_stage1)
