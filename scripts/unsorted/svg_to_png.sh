#!/bin/bash

inDir=${1:-"."}
outDir=${2:-"$inDir"}
dpi=${3:300}

for file in "${inDir}/"*".svg"; do
    base=$(basename $file .svg)
    inkscape $file --export-dpi 300 --export-png "${outDir}/${base}.png"
done
