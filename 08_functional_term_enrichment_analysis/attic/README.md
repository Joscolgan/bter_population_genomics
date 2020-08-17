### Gene Ontology analyses:  
This subdirectory contains scripts for identifying gene ontology terms enriched within a subset of genes of interest.  

For gene ontology enrichment analysis, two approaches were taken:  
1. [TopGO](http://bioconductor.org/packages/release/bioc/html/topGO.html) analysis:  
This software tool performs Gene Ontology enrichment analysis using Fisher and KS-based tests.  
The software also comes with built-in algorithms, which reduce redundancy in identification of enriched terms.   

2. [Gowinda](https://sourceforge.net/p/gowinda/wiki/Main/) analysis:  
This software tool performs Gene Ontology enrichment analysis using a Fisher's exact test.  
The software tool coorects for length-based biases, whereby longer genes may be annotated with greater numbers of gene ontology terms.  
The paper describing the software tool is [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3400962/).
