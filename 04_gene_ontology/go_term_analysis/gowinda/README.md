#### TJ Colgan 2019

For the running of Gowinda, a geneset file is required, which contains the following format:
GO term name, GO term description, genes annotated with that GO term.  
For example:
```
GO:0000002      mitochondrial genome maintenance        CG11077 CG33650 CG4337 CG5924
GO:0000003      reproduction    CG10112 CG10128 CG1262 CG13873 CG14034 CG15117 CG15616 CG1656
```
To obtain this script, the following steps can be taken:  

Using the input file used for topGO (``` ```) that contains two columns:
Column 1: Gene ID
Column 2: List of comma separated assigned gene ontology terms.  

The script ```generate_gowinda_genelist.sh``` creates a unique list of gene ontology terms and each gene assigne to it.

The second script ```go_term_converter.R``` takes the output file of ```generate_gowinda_genelist.sh``` and converts into an output which switches the column information.  
The output will be a tab-delimited text file containing two columns:  
Column 1: GO term ID  
Column 2: A list of comma separated genes assigned that gene ontology term within the input file.  

