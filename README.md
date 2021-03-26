# brgeEnrich

This package contains functions to use online tools to perform enrichments. *brgeEnrich* allows you to work with online enrichment results in R.
The implemented online resources are : 

* **ConsensusPathDB-human (CPDB)**,  that integrates interaction networks in Homo sapiens including binary and complex protein-protein, genetic, metabolic, signaling, gene regulatory and drug-target interactions, as well as biochemical pathways

## Install : 

Install in R as : 

```r
# Install devtools
install.packages("devtools")

# Install required packages 
devtools::source_url("https://raw.githubusercontent.com/isglobal-brge/brgeEnrich/HEAD/installer.R")

# Install brgeEnrich package
devtools::install_github("isglobal-brge/brgeEnrich@HEAD")
```
