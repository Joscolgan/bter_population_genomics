## For detecting deletions within the British <i> Bombus terrestris audax </i> population.  

## Stage One: Running snakefile to generate intersected BAM files for each deleted loci:  
If you do not check for the completion of output, the script runs fine and will generate empty bams for regions with no read coverage.  
The script will stop and complain about empty files when trying to calculate read depth due to the files containing zero read depth.  
To get around this issue, the file ```empty_BAM_filler.py```, which will populate the empty depth.txt with a zero. This will also need to be performed for the third step for the generation of median count files.  

