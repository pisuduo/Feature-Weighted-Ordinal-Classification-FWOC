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

### Application_GSE9782.R (Application_GSE68871.R)
This file contains the process of testing FWOC on the GSE9782 (GSE68871) dataset, including the data retrival, preprocessing and modeling. 

