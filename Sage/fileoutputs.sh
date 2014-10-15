#!/bin/sh
for f in ./outputs/enumerations*.csv
do
  echo "sort $f"
  sort -t\, -k 1,1n -k 2,2n -k 3,3n -k 4,4n $f -o $f
done
