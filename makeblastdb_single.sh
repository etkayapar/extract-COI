#!/bin/bash

fullpath=$1
sp=$(basename $fullpath)
sp=${sp%%.fna}

echo "processing $sp"
makeblastdb -in $fullpath -dbtype nucl -out databases/${sp}
