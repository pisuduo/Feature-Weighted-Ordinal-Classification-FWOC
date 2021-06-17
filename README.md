# Feature-Weighted-Ordinal-Classification-FWOC-
The source R code for FWOC. 

## Required R packages:

- MASS
- caret
- class
- Biobase
- GEOquery
- ordinalgmifs
- BhGLM
- penalizedLDA

### Helper Functions. R
This file contains the helper functions that are necessary for the implementation of FWOC.

### Main Function. R
This file contains the main function of FWOC, need to run Helper Functions. R first.

### Cross Validation.R
This file contains the function to run cross validation for FWOC, need to run Helper Functions. R and Main Function.R first.

### Projection_GSE9782.R (Projection_GSE68871.R)
This file contains the process of testing FWOC on the GSE9782 (GSE68871) dataset, including the data retrival, preprocessing and modeling and achieving the final projections shown in the paper.

### Implementaions.R
This file contains the scripts for implmenting all the four methods based on given training and test datasets, need to run Helper Functions. R first.

### opt_parameters.xlsx
This file contains the chosen parameters for each training for the four methods. The final results are also included.


