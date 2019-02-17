## TJ Colgan 2019  

### Gene Ontology analyses:  
This subdirectory contains scripts for identifying gene ontology terms enriched within a subset of genes of interest.  

For gene ontology enrichment analysis, two approaches were taken:  
1. TopGO analysis:  
This software tool performs Gene Ontology enrichment analysis using Fisher and KS-based tests.  
The software also comes with built-in algorithms, which reduce redundancy in identification of enriched terms.   

2. Gowinda analysis:  
This software tool performs Gene Ontology enrichment analysis using a Fisher's exact test.  
The software tool coorects for length-based biases, whereby longer genes may be annotated with greater numbers of gene ontology terms.  
