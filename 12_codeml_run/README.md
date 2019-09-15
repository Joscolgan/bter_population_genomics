# MSc-research-project-2017

The majority of genomic studies of social insects focused on certain lineages such as honeybees and 
ants. For the bumblebee, an important commercial and ecological pollinator, our understanding of their 
genomics is restricted to two bumblebee genome studies as well as a comparative study that examined signatures
of selection across 10 bee species (incorporating members of the Apidae and Megachilidae) However, there are 
approximately more than 250 species worldwide, 25 of which, are present within the British environment majority 
of which lack genomic resources. In addition, recent wild bumblebee population declines have also been identified 
with habitat fragmentation, novel and reemerging pathogens, as well as pesticide exposure, highlighted as potential 
contributing factors. Therefore, advancing our understanding of local adaptation at the genomic level through examination
of sites undergoing selection is warranted.

In this project we aim to used the branch-site model of CodeMl as part of as part of the phylogenetic analysis of
maximum likelihood (PAML) toolkit to assessed the presence of signatures of selection acting on the protein coding 
sites of four Bombus species. This includes the buff-tailed bumblebee (B. terrestris); the white-tailed bumblebee 
(Bombus lucorum); the cryptic bumblebee (Bombus cryptarum); and the eastern bumblebee (B. impatiens). 

First script: find_recipocal.py:
```
python find_recipocal.py tblastn.output.csv blastx.output.csv output.txt
```

Second script: PRANK can complain
```
runPrank.qsub
```
Third script: Check outside of PRANK - justification of why PRANK did not align samples:
```
python seq_length.py alignment_PAML.fasta
```
Fourth step: Filtering of output alignments:  
- Main issues where have an abundance of gaps  
- Gblocks was used for filtering <update with command>  

Post-filtering: 
```
./launch.sh
```
Alters control file per directory.  


For renaming of trees within each directory containing an aligned orthologue, use:  
```
renameTree.sh
```
