#!/bin/sh
# the code loops in each directory to format the codeml control
# file using the input alignment 
 
for line in */
	do
          	[ -d $line ] && cd "$line"
                for f in *.fas-gb
                        do
                        b="$f".out
                        cp ../../codeml_M0.ctl codeml.ctl
                        perl -pi -e "s/infile/$f/g" codeml.ctl
                        perl -pi -e "s/output/$b/g" codeml.ctl
                        done
                cd ../
        done < fasta_files
