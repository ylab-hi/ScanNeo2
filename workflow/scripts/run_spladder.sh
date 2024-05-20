#!/bin/bash

input=$1
threads=$2
anno=$3
output=$4
confidence=$5
iterations=$6
edgelimit=$7
log=$8

# run spladder
spladder build \
    -b $input \
    --parallel $threads \
    -a $anno \
    -o $output \
    --filter-overlap-exons \
    --no-primary-only \
    --quantify-graph \
    --confidence $confidence \
    --iterations $iterations \
    --ase-edge-limit $edgelimit \
    --qmode all > $log 2>&1

# check if spladder generated output (if not create folder)
if [ ! -d $output ]; then
    echo "Spladder did not generate output. Creating empty folder."
    mkdir -p $output
fi

