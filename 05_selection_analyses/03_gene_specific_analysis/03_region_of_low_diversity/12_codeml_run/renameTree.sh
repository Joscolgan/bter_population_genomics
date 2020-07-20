#!/bin/sh
#the code formats the codeml tree base on branch or clade selected

for line in */
do
 [ -d $line ] && cd "$line" && cp ../../speciesTree.txt species.tree
 cd ..
done
