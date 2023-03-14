# Prints min and max region sizes of a BED file
cat [IN] | awk '{print $3 - $2}' | awk 'NR == 1 { max=$1; min=$1 } { if ($1>max) max=$1; if ($1<min) min=$1; } END {printf "Min: %d\tMax: %d\n", min, max}'

# Size filter BED regions
awk '{ if( ($3 - $2) >= 22 && ($3 - $2) <= 24 ) print $0}' < [IN] > [OUT]


